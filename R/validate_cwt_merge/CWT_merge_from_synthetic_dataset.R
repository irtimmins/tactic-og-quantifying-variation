# =============================================================================
# OG Waiting Times - end-to-end CWT merge test on synthetic data
# -----------------------------------------------------------------------------
# This is the real test: read both SYNTH files cold from disk and run the
# complete merge exactly as a partner would on their real data, with all the
# validation and summary tables.
#
# Inputs:  og_cohort_precwt_SYNTH.rds   (Table A)
#          cwt_records_SYNTH.rds         (Table B, raw dd/mm/yyyy dates)
#          og_pipeline_spec.rds          (globals)
# Output:  og_cohort_cwt_SYNTH.rds       (merged, analysis-ready)
#
# The only difference from the real merge script is that the synthetic CWT
# comes from readRDS() rather than open_dataset() - everything else is
# identical, so this is a faithful proxy for what the partner will run.
# =============================================================================

library(tidyverse)
library(lubridate)

base_dir <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/"
spec_obj       <- readRDS(paste0(base_dir, "og_pipeline_spec.rds"))
og_icd10       <- spec_obj$og_icd10
tx_window_days <- spec_obj$tx_window_days

og_cohort      <- readRDS(paste0(base_dir, "og_cohort_precwt_SYNTH.rds"))
ncras_og_ids   <- og_cohort$pseudo_patientid

cat("Pre-CWT cohort loaded:", nrow(og_cohort), "patients\n")
cat("Columns:", ncol(og_cohort), "\n\n")

# =============================================================================
# 1. Read + parse CWT records
# =============================================================================
cwt_og <- readRDS(paste0(base_dir, "cwt_records_SYNTH.rds")) %>%
  filter(site_icd10 %in% og_icd10) %>%
  mutate(
    pseudo_patientid  = as.character(pseudo_patientid),
    cwt_dtt_date      = as.Date(treat_period_start, format = "%d/%m/%Y"),
    cwt_treat_date    = as.Date(treat_start,        format = "%d/%m/%Y"),
    cwt_referral_date = as.Date(crtp_date,          format = "%d/%m/%Y"),
    cwt_first_seen    = as.Date(date_first_seen,    format = "%d/%m/%Y"),
    cwt_mdt_date      = as.Date(mdt_date,           format = "%d/%m/%Y")
  ) %>%
  filter(pseudo_patientid %in% ncras_og_ids)

cat("CWT rows (synthetic):", nrow(cwt_og),
    "| unique patients:", n_distinct(cwt_og$pseudo_patientid), "\n")

# Records per patient
cat("\nRecords per patient:\n")
cwt_og %>%
  count(pseudo_patientid, name = "k") %>%
  count(k, name = "n_patients") %>%
  mutate(pct = round(100 * n_patients / sum(n_patients), 1)) %>%
  print()

# Modality distribution
cat("\nModality distribution:\n")
cwt_og %>% count(modality, sort = TRUE) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>% print()

# =============================================================================
# 2. CWT anchor: earliest valid DTT per patient within treatment window
# =============================================================================
cwt_anchor <- og_cohort %>%
  select(pseudo_patientid, diagmdy) %>%
  left_join(
    cwt_og %>%
      filter(
        !(modality %in% c("23", "24") & cwt_treat_date < as.Date("2020-06-01")),
        !modality %in% c("97", "98", "99")
      ) %>%
      select(pseudo_patientid, cwt_dtt_date, cwt_treat_date,
             cwt_mdt_date, modality),
    by = "pseudo_patientid"
  ) %>%
  mutate(days_dx_to_dtt = as.integer(cwt_dtt_date - diagmdy)) %>%
  filter(
    !is.na(days_dx_to_dtt),
    days_dx_to_dtt >= -30,
    days_dx_to_dtt <= tx_window_days
  ) %>%
  arrange(pseudo_patientid, cwt_dtt_date) %>%
  distinct(pseudo_patientid, .keep_all = TRUE) %>%
  select(pseudo_patientid, cwt_dtt_date, cwt_treat_date,
         cwt_mdt_date, modality, days_dx_to_dtt)

cat("\nCWT anchor patients (real target ~36,197):", nrow(cwt_anchor), "\n")
cat("DTT completeness:  ", round(100 * mean(!is.na(cwt_anchor$cwt_dtt_date)), 1), "%\n")
cat("MDT completeness:  ", round(100 * mean(!is.na(cwt_anchor$cwt_mdt_date)), 1),
    "% (real ~40%)\n")
cat("\ndays_dx_to_dtt (real: median 39, IQR 24-56):\n")
print(summary(cwt_anchor$days_dx_to_dtt))

# =============================================================================
# 3. Validation against pipeline treatment dates
# =============================================================================
cwt_validation <- cwt_anchor %>%
  left_join(
    og_cohort %>%
      select(pseudo_patientid, diagmdy, first_tx_date,
             surgery_date, sact_date, rt_date, tx_pathway),
    by = "pseudo_patientid"
  ) %>%
  mutate(
    dtt_to_cwt_treat = as.integer(cwt_treat_date - cwt_dtt_date),
    dtt_to_tx        = as.integer(first_tx_date  - cwt_dtt_date),
    cwt_vs_first_tx  = as.integer(cwt_treat_date - first_tx_date)
  )

cat("\n--- DTT to CWT treat date (internal CWT consistency) ---\n")
cat("Real: median 11, IQR 3-18, no negatives\n")
cwt_validation %>%
  filter(!is.na(dtt_to_cwt_treat)) %>%
  summarise(n = n(),
            n_negative = sum(dtt_to_cwt_treat < 0),
            pct_neg    = round(100 * mean(dtt_to_cwt_treat < 0), 1),
            median = median(dtt_to_cwt_treat),
            p25    = quantile(dtt_to_cwt_treat, 0.25),
            p75    = quantile(dtt_to_cwt_treat, 0.75),
            max    = max(dtt_to_cwt_treat)) %>% print()

cat("\n--- DTT to first_tx_date ---\n")
cat("Real: median 14, IQR 7-22, 5.3% negative\n")
cwt_validation %>%
  filter(!is.na(dtt_to_tx)) %>%
  summarise(n = n(),
            n_negative = sum(dtt_to_tx < 0),
            pct_neg    = round(100 * mean(dtt_to_tx < 0), 1),
            median = median(dtt_to_tx),
            p25    = quantile(dtt_to_tx, 0.25),
            p75    = quantile(dtt_to_tx, 0.75),
            max    = max(dtt_to_tx)) %>% print()

cat("\n--- CWT treat date vs first_tx_date ---\n")
cat("Real: 71.1% exact, 85.6% within 14d\n")
v <- cwt_validation %>% filter(!is.na(cwt_vs_first_tx))
v %>%
  summarise(n = n(),
            pct_exact     = round(100 * mean(cwt_vs_first_tx == 0),      1),
            pct_within_5  = round(100 * mean(abs(cwt_vs_first_tx) <= 5), 1),
            pct_within_14 = round(100 * mean(abs(cwt_vs_first_tx) <= 14),1),
            median_diff   = median(cwt_vs_first_tx),
            p25 = quantile(cwt_vs_first_tx, 0.25),
            p75 = quantile(cwt_vs_first_tx, 0.75)) %>% print()

cat("\n--- Negative dtt_to_tx by pathway ---\n")
cat("Real: EMR/ESD then surgery ~50%, EMR/ESD only ~14.5%\n")
cwt_validation %>%
  filter(!is.na(dtt_to_tx)) %>%
  group_by(tx_pathway) %>%
  summarise(n = n(),
            n_neg   = sum(dtt_to_tx < 0),
            pct_neg = round(100 * mean(dtt_to_tx < 0), 1),
            .groups = "drop") %>%
  arrange(desc(pct_neg)) %>% print()

# =============================================================================
# 4. Merge DTT anchor onto cohort + derive intervals and validity flag
# =============================================================================
og_cohort <- og_cohort %>%
  left_join(
    cwt_anchor %>%
      select(pseudo_patientid, cwt_dtt_date, cwt_mdt_date, cwt_treat_date),
    by = "pseudo_patientid"
  ) %>%
  mutate(
    wt_endo_to_dtt = as.integer(cwt_dtt_date - endoscopy_date),
    wt_dtt_to_tx   = as.integer(first_tx_date - cwt_dtt_date),
    wt_dx_to_dtt   = as.integer(cwt_dtt_date - diagmdy),
    dtt_valid = !is.na(cwt_dtt_date) &
      wt_dx_to_dtt >= 0 &
      wt_dtt_to_tx >= -14,
    dtt_valid = if_else(
      tx_pathway %in% c("EMR/ESD only", "EMR/ESD then surgery"),
      NA, dtt_valid
    )
  )

saveRDS(og_cohort, paste0(base_dir, "og_cohort_cwt_SYNTH.rds"))
cat("\nSaved og_cohort_cwt_SYNTH.rds\n")

# =============================================================================
# 5. Post-merge summaries
# =============================================================================
cat("\n--- dtt_valid by pathway ---\n")
cat("Real: curative RT ~99.6%, neoadj chemo ~95.7%, surgery only ~96.3%\n")
og_cohort %>%
  filter(!is.na(cwt_dtt_date)) %>%
  count(dtt_valid, tx_pathway) %>%
  group_by(tx_pathway) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  filter(!is.na(dtt_valid)) %>%
  arrange(tx_pathway) %>%
  print(n = 30)

cat("\n--- Staging (endo->DTT) and scheduling (DTT->tx) intervals ---\n")
cat("Real: endo->DTT median 44 (IQR 34-59), DTT->tx median 15 (IQR 9-23)\n")
og_cohort %>%
  filter(dtt_valid == TRUE,
         !is.na(wt_endo_to_dtt), wt_endo_to_dtt >= 0,
         !is.na(wt_dtt_to_tx),   wt_dtt_to_tx   >= 0) %>%
  summarise(
    n               = n(),
    median_endo_dtt = median(wt_endo_to_dtt),
    p25_endo_dtt    = quantile(wt_endo_to_dtt, 0.25),
    p75_endo_dtt    = quantile(wt_endo_to_dtt, 0.75),
    median_dtt_tx   = median(wt_dtt_to_tx),
    p25_dtt_tx      = quantile(wt_dtt_to_tx, 0.25),
    p75_dtt_tx      = quantile(wt_dtt_to_tx, 0.75)
  ) %>% print()

cat("\n--- Intervals by pathway ---\n")
og_cohort %>%
  filter(dtt_valid == TRUE,
         !is.na(wt_endo_to_dtt), wt_endo_to_dtt >= 0,
         !is.na(wt_dtt_to_tx),   wt_dtt_to_tx   >= 0) %>%
  group_by(tx_pathway) %>%
  summarise(
    n               = n(),
    median_endo_dtt = median(wt_endo_to_dtt),
    p25_endo_dtt    = quantile(wt_endo_to_dtt, 0.25),
    p75_endo_dtt    = quantile(wt_endo_to_dtt, 0.75),
    median_dtt_tx   = median(wt_dtt_to_tx),
    p25_dtt_tx      = quantile(wt_dtt_to_tx, 0.25),
    p75_dtt_tx      = quantile(wt_dtt_to_tx, 0.75),
    .groups         = "drop"
  ) %>%
  arrange(desc(n)) %>% print()

cat("\n--- Deprivation gradient (staging vs scheduling component) ---\n")
cat("Real: endo->tx medians 62/61/60/61/59 across IMD 1-> 5 (very flat)\n")
og_cohort %>%
  filter(dtt_valid == TRUE,
         !is.na(wt_endo_to_dtt), wt_endo_to_dtt >= 0,
         !is.na(wt_dtt_to_tx),   wt_dtt_to_tx   >= 0,
         !is.na(NHSE_reversed_imd_quintile_lsoas)) %>%
  group_by(NHSE_reversed_imd_quintile_lsoas) %>%
  summarise(
    n               = n(),
    median_endo_tx  = median(wt_endo_to_tx,  na.rm = TRUE),
    median_endo_dtt = median(wt_endo_to_dtt, na.rm = TRUE),
    median_dtt_tx   = median(wt_dtt_to_tx,   na.rm = TRUE),
    .groups         = "drop"
  ) %>% print()

cat("\n--- Missing reason summary ---\n")
og_cohort %>%
  mutate(missing_reason = case_when(
    is.na(endoscopy_date) & is.na(first_tx_date)  ~ "no endoscopy, no treatment",
    is.na(endoscopy_date) & !is.na(first_tx_date) ~ "no endoscopy, has treatment",
    !is.na(endoscopy_date) & is.na(first_tx_date) ~ "has endoscopy, no treatment",
    wt_endo_to_tx < 0                             ~ "negative wait",
    TRUE                                          ~ "complete"
  )) %>%
  count(missing_reason) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>% print()


