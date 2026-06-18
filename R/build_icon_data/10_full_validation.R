# =============================================================================
# 10  Full validation
# -----------------------------------------------------------------------------
# Concrete pass/fail tests on the real built cohort, for absolute reassurance
# that a rebuild has not drifted. Three tiers:
#   structural  - row counts, keys, required columns, valid factor levels
#   logical     - invariants that must hold by construction (flag consistency,
#                 date ordering, the NA patterns that define each pathway)
#   benchmark   - distributions within tolerance of the known audit figures
#
# Every check is collected; the script prints a tier-by-tier result and then
# STOPS with a non-zero count if any check fails. Structural and logical checks
# are hard (a failure is a bug); benchmark checks use deliberately wide bands so
# they flag real drift, not noise. Tune the bands to the published audit if you
# want them tighter.
#
# Reads: Data/ICON/og_cohort_cwt_2015_2022.rds  (and the chemo anchor, if present)
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

og <- readRDS(f_cohort_cwt)
n  <- nrow(og)

# --- harness: record, classify by tier, print, then stop on any failure ------
.checks <- list()
chk <- function(tier, label, cond, detail = "") {
  ok <- isTRUE(cond)
  .checks[[length(.checks) + 1]] <<-
    list(tier = tier, label = label, ok = ok, detail = detail)
}
# benchmark helper: value within [lo, hi]
within <- function(x, lo, hi) is.finite(x) && x >= lo && x <= hi

stage13 <- og %>% filter(stage_clean %in% c("1","2","3"))
curative_pathways <- c("EMR/ESD only","EMR/ESD then surgery",
                       "Surgery + neoadjuvant chemoRT","Surgery + neoadjuvant chemo",
                       "Surgery + neoadjuvant RT","Surgery + adjuvant chemo","Surgery only",
                       "Surgery + other","Definitive chemoRT","Curative RT only")
noncurative_pathways <- c("Palliative chemo + RT","SACT only",
                          "Palliative RT only","No treatment recorded")

# =============================================================================
# Structural
# =============================================================================
chk("structural", "cohort size in plausible range (35k-45k)", within(n, 35000, 45000),
    paste(n, "patients"))
chk("structural", "no duplicate patient ids",
    sum(duplicated(og$pseudo_patientid)) == 0)
req_cols <- c("pseudo_patientid","diagmdy","stage_clean","cancer_subtype",
              "tumour_site_grp","tx_pathway","first_tx_date","tx_trust","cwt_dtt_date",
              "cwt_treat_date","cwt_modality","dtt_valid","tx_modality_audit",
              "tx_intent_audit","received_any_tx","received_curative_tx_audit",
              "wt_dx_to_dtt","wt_dtt_to_tx","wt_dx_to_tx","cci_group")
miss <- setdiff(req_cols, names(og))
chk("structural", "all required columns present", length(miss) == 0,
    if (length(miss)) paste("missing:", paste(miss, collapse = ", ")) else "")
chk("structural", "tx_pathway only uses defined levels",
    length(setdiff(unique(og$tx_pathway), tx_pathway_levels)) == 0)
chk("structural", "stage_clean is 1-3 (the analysis cohort)",
    all(og$stage_clean %in% c("1","2","3")))
chk("structural", "cci_group only uses 0/1/2/3+",
    all(as.character(og$cci_group) %in% c("0","1","2","3+")))
chk("structural", "diagnosis dates within the study window 2015-2022",
    all(og$diagmdy >= as.Date("2015-01-01") & og$diagmdy <= as.Date("2022-12-31"), na.rm = TRUE))

# =============================================================================
# Logical invariants
# =============================================================================
# first_tx_date is missing for exactly the non-curative pathways (Surgery+other
# is curative and anchored on surgery, so it is not missing)
chk("logical", "first_tx_date present for every curative pathway",
    all(!is.na(og$first_tx_date[og$tx_pathway %in% curative_pathways])))
chk("logical", "first_tx_date missing for every non-curative pathway",
    all(is.na(og$first_tx_date[og$tx_pathway %in% noncurative_pathways])))

# tx_trust present exactly when the pathway has a curative clock-stop with a
# provider (all curative pathways here), missing otherwise
chk("logical", "tx_trust present for curative pathways",
    all(!is.na(og$tx_trust[og$tx_pathway %in% curative_pathways])))
chk("logical", "tx_trust missing for non-curative pathways",
    all(is.na(og$tx_trust[og$tx_pathway %in% noncurative_pathways])))

# curative implies any-treatment (never the reverse failure)
chk("logical", "received_curative implies received_any (no row curative-not-any)",
    sum(og$received_curative_tx_audit & !og$received_any_tx, na.rm = TRUE) == 0)
chk("logical", "received_any and received_curative are not identical",
    !identical(og$received_any_tx, og$received_curative_tx_audit))

# date ordering where the relevant dates are present
chk("logical", "no decision-to-treat before diagnosis",
    sum(og$wt_dx_to_dtt < 0, na.rm = TRUE) == 0)
chk("logical", "no first treatment before diagnosis",
    sum(og$first_tx_date < og$diagmdy, na.rm = TRUE) == 0)
chk("logical", "treatment not before DTT beyond tolerance",
    sum(og$wt_dtt_to_tx < -treat_tol_days, na.rm = TRUE) == 0)

# dtt_valid is NA exactly for the EMR/ESD pathways (by design)
emr_pw <- c("EMR/ESD only","EMR/ESD then surgery")
chk("logical", "dtt_valid is NA for EMR/ESD pathways",
    all(is.na(og$dtt_valid[og$tx_pathway %in% emr_pw])))
chk("logical", "dtt_valid is non-NA for non-EMR pathways with a DTT",
    all(!is.na(og$dtt_valid[!og$tx_pathway %in% emr_pw & !is.na(og$cwt_dtt_date)])))

# audit category and intent are internally consistent
chk("logical", "every treated patient has an audit modality",
    all(!is.na(og$tx_modality_audit[og$tx_pathway != "No treatment recorded"])))
chk("logical", "curative audit flag matches Curative intent",
    all(og$received_curative_tx_audit[og$tx_intent_audit != "Curative"] == FALSE |
          is.na(og$tx_intent_audit)))

# pathway proportions form a complete partition
chk("logical", "pathway mix sums to 100%",
    abs(sum(prop.table(table(og$tx_pathway))) - 1) < 1e-9)

# =============================================================================
# Benchmark distributions (wide tolerance bands -> flag drift, not noise)
# =============================================================================
pc <- function(x) 100 * mean(x, na.rm = TRUE)
cur_rate <- pc(stage13$received_curative_tx_audit)
any_rate <- pc(stage13$received_any_tx)
chk("benchmark", "stage 1-3 curative rate in 40-60%", within(cur_rate, 40, 60),
    sprintf("%.0f%%", cur_rate))
chk("benchmark", "stage 1-3 any-treatment rate in 65-85%", within(any_rate, 65, 85),
    sprintf("%.0f%%", any_rate))

scc <- stage13 %>% filter(cancer_subtype == "Oes SCC")
gas <- stage13 %>% filter(cancer_subtype == "Gast")
scc_dcrt <- pc(scc$tx_modality_audit == "Definitive chemoRT")
gas_surg <- pc(gas$tx_modality_audit == "Surgery only")
chk("benchmark", "SCC definitive chemoRT in 12-30%", within(scc_dcrt, 12, 30),
    sprintf("%.0f%%", scc_dcrt))
chk("benchmark", "gastric surgery-only in 8-20%", within(gas_surg, 8, 20),
    sprintf("%.0f%%", gas_surg))

cwt_cov <- pc(!is.na(og$cwt_dtt_date))
chk("benchmark", "CWT decision-to-treat coverage in 80-95%", within(cwt_cov, 80, 95),
    sprintf("%.0f%%", cwt_cov))

# no-treatment leakage: active-modality CWT in window among "No treatment"
leak <- og %>% filter(tx_pathway == "No treatment recorded", !is.na(cwt_treat_date),
                      cwt_modality %in% c("01","23","24","02","04","05")) %>%
  mutate(d = as.integer(cwt_treat_date - diagmdy)) %>%
  filter(d >= 0, d <= cwt_window_days)
leak_pct <- 100 * nrow(leak) / sum(og$tx_pathway == "No treatment recorded")
chk("benchmark", "no-treatment leakage below 3%", leak_pct < 3,
    sprintf("%.1f%%", leak_pct))

# neoadjuvant anchoring: CWT lands on chemo/RT, not surgery
neo <- og %>% filter(tx_pathway %in% c("Surgery + neoadjuvant chemo",
                                       "Surgery + neoadjuvant chemoRT","Surgery + neoadjuvant RT"), !is.na(cwt_modality))
neo_ok <- pc(!neo$cwt_modality %in% c("01","23","24"))
chk("benchmark", "neoadjuvant anchored on chemo/RT >= 90%", neo_ok >= 90,
    sprintf("%.0f%%", neo_ok))

# dtt_valid health on the non-EMR curative pathways
dtt_ok <- pc(og$dtt_valid[!og$tx_pathway %in% emr_pw])
chk("benchmark", "dtt_valid TRUE share (non-EMR) >= 80%", dtt_ok >= 80,
    sprintf("%.0f%%", dtt_ok))

# =============================================================================
# Report and stop on any failure
# =============================================================================
res <- tibble(tier = map_chr(.checks, "tier"), label = map_chr(.checks, "label"),
              ok = map_lgl(.checks, "ok"), detail = map_chr(.checks, "detail"))
for (t in c("structural","logical","benchmark")) {
  cat("\n==========  ", t, "  ==========\n")
  res %>% filter(tier == t) %>%
    transmute(result = if_else(ok, "pass", "FAIL"), label,
              detail = if_else(detail == "", "", paste0("(", detail, ")"))) %>%
    pwalk(function(result, label, detail)
      cat(sprintf("  %-5s %s %s\n", result, label, detail)))
}
cat(sprintf("\n==========  %d passed, %d failed, %d total  ==========\n",
            sum(res$ok), sum(!res$ok), nrow(res)))
if (any(!res$ok))
  stop(sum(!res$ok), " validation check(s) failed - the build has drifted; ",
       "see the FAIL lines above. Diagnostics: 11_check_treatment_classification.R",
       call. = FALSE)
cat("full validation passed - the build is sound.\n")