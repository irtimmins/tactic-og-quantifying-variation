# =============================================================================
# 10  Endoscopy-site concordance against the registry (assessment only)
# -----------------------------------------------------------------------------
# The endoscopy-to-decision-to-treat positive deviance analysis groups patients
# by the hospital that did the diagnostic endoscopy. Script 04 now carries that
# site from HES APC (sitetret_endoscopy_apc, a 5-character site) through to the
# cohort. This script checks how well it agrees with the registry's own
# diagnosing site/trust, and sets it directly against the CWT-derived diagnosing
# site (cwt_diag_site_v2) so we can see which is the better grouping variable.
#
# It writes a small summary and prints the comparison. It does NOT change the
# cohort - which site variable the deviance analysis ends up using is a later,
# separate decision.
#
# Why compare against the registry at all: the registry diag_hosp / diag_trust
# is the reference for where diagnosis happened. The HES endoscopy site is the
# site of the diagnostic procedure itself. For a straightforward pathway these
# are the same hospital; where they differ, it is usually a referred patient
# (scoped locally, decision taken at a centre) or a coding-granularity gap.
#
# Code formats, which drive how agreement is measured (same as script 09):
#   sitetret_endoscopy_apc  5-character site (APC); NA for the OP-sourced minority
#   diag_hosp               MIXED - mostly 5-character site, some 3-character trust
#   diag_trust              3-character trust, essentially complete
#   cwt_diag_site_v2        5-character site (derived in 09)
# We report agreement at two levels, exact (full string) and trust (first three
# characters), because a plain match understates agreement across mixed formats.
#
# Reads : Data/ICON/og_cohort_cwt_2015_2022.rds (the cohort), and
#         cwt_diag_site_v2.rds (the derived CWT diagnosing site, saved by 09).
#         Run 09 first so that side table exists; if it is missing, the CWT
#         comparison is skipped and only the HES endoscopy site is scored.
# Writes: cwt_vs_endoscopy_site_concordance.csv
#         endoscopy_site_performance.csv
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

library(dplyr)
library(stringr)

first3   <- function(x) str_sub(x, 1, 3)
is_site5 <- function(x) !is.na(x) & nchar(x) == 5

out_dir <- if (exists("dir_out")) dir_out else
  if (exists("dir_icon")) dir_icon else "."

# -----------------------------------------------------------------------------
# 1. Read the cohort. Prefer the CWT-merged cohort so cwt_diag_site_v2 can be
#    compared alongside; fall back to the plain cohort if the merged one is not
#    present, in which case the CWT row is simply skipped.
# -----------------------------------------------------------------------------
cohort_path <- if (exists("f_cohort_cwt") && file.exists(f_cohort_cwt)) f_cohort_cwt else f_cohort
dat <- readRDS(cohort_path) %>%
  mutate(across(any_of(c("diag_hosp", "diag_trust", "sitetret_endoscopy_apc",
                         "procode3_endo_apc", "cwt_diag_site_v2")),
                ~ na_if(str_trim(as.character(.x)), "")))

# the CWT-derived diagnosing site is built and saved by script 09 as a side
# table, not carried in the cohort. Join it on if it is there; if 09 has not
# been run, the comparison against it is simply skipped.
cwt_site_file <- file.path(out_dir, "cwt_diag_site_v2.rds")
if (!("cwt_diag_site_v2" %in% names(dat)) && file.exists(cwt_site_file)) {
  cwt_v2 <- readRDS(cwt_site_file) %>%
    mutate(pseudo_patientid = as.character(pseudo_patientid),
           cwt_diag_site_v2 = na_if(str_trim(as.character(cwt_diag_site_v2)), ""))
  dat <- dat %>%
    mutate(pseudo_patientid = as.character(pseudo_patientid)) %>%
    left_join(cwt_v2, by = "pseudo_patientid")
}

have_cwt <- "cwt_diag_site_v2" %in% names(dat)
cat("Cohort read from:", cohort_path, "\n")
if (have_cwt) cat("cwt_diag_site_v2 joined from:", cwt_site_file, "\n")
cat("rows:", nrow(dat), " cwt_diag_site_v2 present:", have_cwt, "\n\n")

# -----------------------------------------------------------------------------
# 2. Agreement helper: exact and trust level, plus the 5-character breakdown.
#    (Identical in spirit to the concordance() in script 09.)
# -----------------------------------------------------------------------------
concordance <- function(data, cwt_col, ref_col, label) {
  d <- data %>%
    transmute(cwt = .data[[cwt_col]], ref = .data[[ref_col]]) %>%
    filter(!is.na(cwt), !is.na(ref))
  
  n_pair <- nrow(d)
  if (n_pair == 0) {
    cat(sprintf("%-46s  no overlapping records\n", label))
    return(tibble(comparison = label, n_compared = 0L,
                  pct_exact = NA_real_, pct_trust = NA_real_,
                  pct_both_5char = NA_real_, exact_within_5char = NA_real_))
  }
  
  exact <- mean(d$cwt == d$ref)
  trust <- mean(first3(d$cwt) == first3(d$ref))
  both5 <- is_site5(d$cwt) & is_site5(d$ref)
  exact_in5 <- if (any(both5)) mean(d$cwt[both5] == d$ref[both5]) else NA_real_
  
  cat(sprintf("%-46s  n=%6d\n", label, n_pair))
  cat(sprintf("    all pairs   : exact %5.1f%%   trust %5.1f%%\n", 100*exact, 100*trust))
  cat(sprintf("    both 5-char : %5.1f%% of pairs, site agreement within those %5.1f%%\n",
              100*mean(both5), 100*exact_in5))
  
  tibble(comparison = label, n_compared = n_pair,
         pct_exact = round(100*exact, 1), pct_trust = round(100*trust, 1),
         pct_both_5char = round(100*mean(both5), 1),
         exact_within_5char = round(100*exact_in5, 1))
}

# -----------------------------------------------------------------------------
# 3. Endoscopy site vs the registry, and the CWT-derived site for contrast.
# -----------------------------------------------------------------------------
cat("Diagnosing/endoscopy site vs registry diag_hosp\n")
cat(strrep("-", 78), "\n")
res_endo <- concordance(dat, "sitetret_endoscopy_apc", "diag_hosp",
                        "HES endoscopy site (APC)")
res_cwt  <- if (have_cwt)
  concordance(dat, "cwt_diag_site_v2", "diag_hosp",
              "CWT derived diagnosing site") else NULL
cat("\n")

# -----------------------------------------------------------------------------
# 4. Trust-level decomposition: of the endoscopy sites, how many sit in the
#    registry diagnosing trust, and how many are a genuinely different trust.
# -----------------------------------------------------------------------------
trust_report <- function(data, cwt_col, label) {
  d <- data %>%
    transmute(cwt = .data[[cwt_col]], diag_trust) %>%
    filter(!is.na(cwt), !is.na(diag_trust))
  if (nrow(d) == 0) { cat(sprintf("%-30s no records\n", label)); return(invisible()) }
  in_trust <- mean(first3(d$cwt) == d$diag_trust)
  cat(sprintf("%-30s  in registry trust %5.1f%%   different trust %5.1f%%  (n=%d)\n",
              label, 100*in_trust, 100*(1-in_trust), nrow(d)))
  invisible()
}

cat("Trust-level agreement with the registry diagnosing trust\n")
cat(strrep("-", 78), "\n")
trust_report(dat, "sitetret_endoscopy_apc", "HES endoscopy site")
if (have_cwt) trust_report(dat, "cwt_diag_site_v2", "CWT derived site")
cat("\n")

# -----------------------------------------------------------------------------
# 5. Performance and coverage of the two candidate grouping variables, side by
#    side. This is the table to read when choosing which to build the deviance
#    analysis on: coverage over the whole cohort, then agreement among covered.
# -----------------------------------------------------------------------------
n_cohort <- nrow(dat)

perf <- function(data, site_col, label) {
  covered <- sum(!is.na(data[[site_col]]))
  d <- data %>%
    transmute(site = .data[[site_col]], diag_hosp, diag_trust) %>%
    filter(!is.na(site))
  tibble(variable = label,
         coverage_pct      = round(100 * covered / n_cohort, 1),
         n_covered         = covered,
         correct_trust_pct = round(100 * mean(first3(d$site) == d$diag_trust, na.rm = TRUE), 1),
         exact_site_pct    = round(100 * mean(d$site == d$diag_hosp, na.rm = TRUE), 1),
         full_site_pct     = round(100 * mean(is_site5(d$site)), 1))
}

perf_tbl <- bind_rows(
  perf(dat, "sitetret_endoscopy_apc", "HES endoscopy site (APC)"),
  if (have_cwt) perf(dat, "cwt_diag_site_v2", "CWT derived diagnosing site"))

cat("Candidate grouping variables, side by side\n")
cat(strrep("-", 78), "\n")
cat("coverage_pct is of the whole cohort; agreement columns are among covered\n")
cat("patients that also have a registry code.\n\n")
print(as.data.frame(perf_tbl))
cat("\n")

# -----------------------------------------------------------------------------
# 6. Where the two derived sites are both present, do THEY agree with each other?
#    High mutual agreement means the deviance ranking is robust to the choice.
# -----------------------------------------------------------------------------
if (have_cwt) {
  both <- dat %>%
    filter(!is.na(sitetret_endoscopy_apc), !is.na(cwt_diag_site_v2))
  if (nrow(both) > 0) {
    agree_site  <- mean(both$sitetret_endoscopy_apc == both$cwt_diag_site_v2)
    agree_trust <- mean(first3(both$sitetret_endoscopy_apc) == first3(both$cwt_diag_site_v2))
    cat("HES endoscopy site vs CWT derived site, head to head\n")
    cat(strrep("-", 78), "\n")
    cat(sprintf("both present for %d patients: agree on site %5.1f%%, on trust %5.1f%%\n\n",
                nrow(both), 100*agree_site, 100*agree_trust))
  }
}

# -----------------------------------------------------------------------------
# 7. Save the summaries (counts and percentages only, no patient data).
# -----------------------------------------------------------------------------
concordance_tbl <- bind_rows(res_endo, res_cwt)
write.csv(concordance_tbl, file.path(out_dir, "cwt_vs_endoscopy_site_concordance.csv"),
          row.names = FALSE)
write.csv(perf_tbl,        file.path(out_dir, "endoscopy_site_performance.csv"),
          row.names = FALSE)

cat("Saved:\n")
cat("  ", file.path(out_dir, "cwt_vs_endoscopy_site_concordance.csv"), "\n")
cat("  ", file.path(out_dir, "endoscopy_site_performance.csv"), "\n")
cat("10 assessment complete. No cohort files were changed.\n")