# =============================================================================
# 12  Diagnosing-site accuracy on the analysis cohort (assessment only)
# -----------------------------------------------------------------------------
# Brings the site-concordance work together on the ACTUAL analysis cohort for
# the endoscopy-to-decision-to-treat positive deviance study, rather than the
# whole registry cohort. Everything is restricted to:
#   stage 1-3, curative intent, a valid endoscopy-to-DTT interval, and a
#   5-character endoscopy site.
# On that group it reports, for each diagnosing-site measure:
#   - concordance with the registry diagnosing site (diag_hosp), exact and trust
#   - agreement with the registry diagnosing trust (diag_trust)
#   - coverage and full-site share
# and sets the two CWT/HES measures head to head, so the accuracy figures quoted
# in the analysis are the ones that apply to the patients actually analysed.
#
# The two measures compared:
#   sitetret_endoscopy_apc  HES site of the diagnostic endoscopy (5-char, APC)
#   cwt_diag_site_v2        CWT-derived diagnosing site (earliest event, trust-
#                           refereed; built and saved by script 09)
#
# Reads : the CWT-merged cohort, and cwt_diag_site_v2.rds (from 09)
# Writes: diag_site_accuracy_analysis_cohort.csv
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

library(dplyr)
library(stringr)

first3   <- function(x) str_sub(x, 1, 3)
is_site5 <- function(x) !is.na(x) & nchar(x) == 5

out_dir <- if (exists("dir_out")) dir_out else
  if (exists("dir_icon")) dir_icon else "."

endo_dtt_max_days <- 365   # upper bound on a sensible endoscopy-to-DTT gap

# -----------------------------------------------------------------------------
# 1. Read the cohort and attach the CWT-derived diagnosing site.
# -----------------------------------------------------------------------------
cohort_path <- if (exists("f_cohort_cwt") && file.exists(f_cohort_cwt)) f_cohort_cwt else f_cohort
dat <- readRDS(cohort_path) %>%
  mutate(across(any_of(c("diag_hosp", "diag_trust", "sitetret_endoscopy_apc",
                         "cwt_diag_site_v2")),
                ~ na_if(str_trim(as.character(.x)), "")),
         pseudo_patientid = as.character(pseudo_patientid))

cwt_site_file <- file.path(out_dir, "cwt_diag_site_v2.rds")
if (!("cwt_diag_site_v2" %in% names(dat)) && file.exists(cwt_site_file)) {
  cwt_v2 <- readRDS(cwt_site_file) %>%
    mutate(pseudo_patientid = as.character(pseudo_patientid),
           cwt_diag_site_v2 = na_if(str_trim(as.character(cwt_diag_site_v2)), ""))
  dat <- left_join(dat, cwt_v2, by = "pseudo_patientid")
}
have_cwt <- "cwt_diag_site_v2" %in% names(dat)

# -----------------------------------------------------------------------------
# 2. Define the analysis cohort. Stage 1-3, curative intent, a valid endoscopy-
#    to-DTT interval, and a 5-character endoscopy site (the grouping variable).
# -----------------------------------------------------------------------------
analysis <- dat %>%
  filter(
    stage_clean %in% c("1", "2", "3"),
    received_curative_tx %in% TRUE,
    !is.na(endoscopy_date), !is.na(cwt_dtt_date),
    !is.na(wt_endo_to_dtt), wt_endo_to_dtt >= 0, wt_endo_to_dtt <= endo_dtt_max_days,
    is_site5(sitetret_endoscopy_apc))

n_full     <- nrow(dat)
n_analysis <- nrow(analysis)
cat("Whole cohort            :", n_full, "patients\n")
cat("Analysis cohort         :", n_analysis,
    sprintf("(%.1f%% of the whole cohort)\n", 100 * n_analysis / n_full))
cat("  stage 1-3, curative intent, valid endoscopy-to-DTT, 5-char endoscopy site\n\n")

# -----------------------------------------------------------------------------
# 3. Concordance helper (exact + trust + 5-char breakdown), as in 09/10.
# -----------------------------------------------------------------------------
concordance <- function(data, site_col, ref_col, label) {
  d <- data %>%
    transmute(site = .data[[site_col]], ref = .data[[ref_col]]) %>%
    filter(!is.na(site), !is.na(ref))
  n_pair <- nrow(d)
  if (n_pair == 0) {
    cat(sprintf("%-30s  no overlapping records\n", label))
    return(tibble(measure = label, n_compared = 0L, pct_exact = NA_real_,
                  pct_trust = NA_real_, exact_within_5char = NA_real_))
  }
  exact <- mean(d$site == d$ref)
  trust <- mean(first3(d$site) == first3(d$ref))
  both5 <- is_site5(d$site) & is_site5(d$ref)
  exact_in5 <- if (any(both5)) mean(d$site[both5] == d$ref[both5]) else NA_real_
  cat(sprintf("%-30s  n=%6d   exact %5.1f%%   trust %5.1f%%   site-within-5char %5.1f%%\n",
              label, n_pair, 100*exact, 100*trust, 100*exact_in5))
  tibble(measure = label, n_compared = n_pair,
         pct_exact = round(100*exact,1), pct_trust = round(100*trust,1),
         exact_within_5char = round(100*exact_in5,1))
}

cat("Concordance with the registry diagnosing site (diag_hosp), analysis cohort\n")
cat(strrep("-", 82), "\n")
res_endo <- concordance(analysis, "sitetret_endoscopy_apc", "diag_hosp",
                        "HES endoscopy site")
res_cwt  <- if (have_cwt)
  concordance(analysis, "cwt_diag_site_v2", "diag_hosp", "CWT derived site") else NULL
cat("\n")

# -----------------------------------------------------------------------------
# 4. Coverage + trust accuracy within the analysis cohort, side by side.
# -----------------------------------------------------------------------------
perf <- function(data, site_col, label) {
  covered <- sum(!is.na(data[[site_col]]))
  d <- data %>% transmute(site = .data[[site_col]], diag_hosp, diag_trust) %>%
    filter(!is.na(site))
  tibble(measure = label,
         coverage_pct      = round(100 * covered / nrow(data), 1),
         n_covered         = covered,
         correct_trust_pct = round(100 * mean(first3(d$site) == d$diag_trust, na.rm = TRUE), 1),
         exact_site_pct    = round(100 * mean(d$site == d$diag_hosp, na.rm = TRUE), 1),
         full_site_pct     = round(100 * mean(is_site5(d$site)), 1))
}

perf_tbl <- bind_rows(
  perf(analysis, "sitetret_endoscopy_apc", "HES endoscopy site"),
  if (have_cwt) perf(analysis, "cwt_diag_site_v2", "CWT derived site"))

cat("Coverage and accuracy within the analysis cohort\n")
cat(strrep("-", 82), "\n")
print(as.data.frame(perf_tbl))
cat("\n")

# -----------------------------------------------------------------------------
# 5. Head to head: where both measures exist in the analysis cohort, how often
#    do they agree with each other? This is the robustness check for the ranking.
# -----------------------------------------------------------------------------
if (have_cwt) {
  both <- analysis %>%
    filter(!is.na(sitetret_endoscopy_apc), !is.na(cwt_diag_site_v2))
  if (nrow(both) > 0) {
    agree_site  <- mean(both$sitetret_endoscopy_apc == both$cwt_diag_site_v2)
    agree_trust <- mean(first3(both$sitetret_endoscopy_apc) == first3(both$cwt_diag_site_v2))
    cat("HES endoscopy site vs CWT derived site, head to head (analysis cohort)\n")
    cat(strrep("-", 82), "\n")
    cat(sprintf("both present for %d patients: agree on site %5.1f%%, on trust %5.1f%%\n\n",
                nrow(both), 100*agree_site, 100*agree_trust))
  }
}

# -----------------------------------------------------------------------------
# 6. Save the accuracy table for the analysis cohort.
# -----------------------------------------------------------------------------
accuracy_tbl <- perf_tbl %>%
  left_join(bind_rows(res_endo, res_cwt) %>%
              select(measure, pct_exact, pct_trust, exact_within_5char),
            by = "measure")

write.csv(accuracy_tbl, file.path(out_dir, "diag_site_accuracy_analysis_cohort.csv"),
          row.names = FALSE)
cat("Saved:", file.path(out_dir, "diag_site_accuracy_analysis_cohort.csv"), "\n")
cat("12 assessment complete. No cohort files were changed.\n")