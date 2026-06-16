#############

# =============================================================================
# OG Cancer Waiting Times - CWT merge (script 4, clean)
# Joins CWT decision-to-treat (DTT) anchor onto the canonical pre-CWT cohort.
# Requires globals: og_icd10, tx_window_days (load from og_pipeline_spec.rds)
# =============================================================================

library(tidyverse)
library(arrow)
library(lubridate)

base_dir <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/"

# --- Globals + cohort --------------------------------------------------------
spec_obj       <- readRDS(paste0(base_dir, "og_pipeline_spec.rds"))
og_icd10       <- spec_obj$og_icd10
tx_window_days <- spec_obj$tx_window_days

ncras_og     <- readRDS(paste0(base_dir, "ncras_og_2015_2022.rds"))
ncras_og_ids <- ncras_og %>% distinct(pseudo_patientid) %>% pull()

og_cohort <- readRDS(paste0(base_dir, "og_cohort_precwt_spec_2015_2022.rds"))

# =============================================================================
# 1. Read + parse CWT for OG patients
# =============================================================================
cwt_og <- open_dataset(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/CWT/11_CWT_data_partitioned"
) %>%
  filter(site_icd10 %in% og_icd10) %>%
  collect() %>%
  mutate(
    pseudo_patientid  = as.character(pseudo_patientid),
    cwt_dtt_date      = as.Date(treat_period_start, format = "%d/%m/%Y"),
    cwt_treat_date    = as.Date(treat_start,        format = "%d/%m/%Y"),
    cwt_referral_date = as.Date(crtp_date,          format = "%d/%m/%Y"),
    cwt_first_seen    = as.Date(date_first_seen,    format = "%d/%m/%Y"),
    cwt_mdt_date      = as.Date(mdt_date,           format = "%d/%m/%Y")
  ) %>%
  filter(pseudo_patientid %in% ncras_og_ids)

cat("CWT OG rows:", nrow(cwt_og),
    "| Patients:", n_distinct(cwt_og$pseudo_patientid), "\n")

# =============================================================================
# 2. CWT anchor: earliest valid DTT per patient within treatment window
# =============================================================================
cwt_anchor <- ncras_og %>%
  select(pseudo_patientid, diagmdy) %>%
  left_join(
    cwt_og %>%
      filter(
        # endoscopic codes only valid from mid-2020
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

cat("CWT anchor patients:", nrow(cwt_anchor), "\n")
cat("DTT completeness:   ", round(100 * mean(!is.na(cwt_anchor$cwt_dtt_date)), 1), "%\n")
cat("MDT completeness:   ", round(100 * mean(!is.na(cwt_anchor$cwt_mdt_date)), 1), "%\n")
summary(cwt_anchor$days_dx_to_dtt)

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
    dtt_to_cwt_treat = as.integer(cwt_treat_date - cwt_dtt_date),  # CWT internal
    dtt_to_tx        = as.integer(first_tx_date  - cwt_dtt_date),  # DTT -> our tx
    cwt_vs_first_tx  = as.integer(cwt_treat_date - first_tx_date)  # CWT tx vs ours
  )

# --- DTT to CWT treatment date (internal CWT consistency) --------------------
cat("\nDTT to CWT treat date:\n")
cwt_validation %>%
  filter(!is.na(dtt_to_cwt_treat)) %>%
  summarise(
    n = n(),
    n_negative = sum(dtt_to_cwt_treat < 0),
    pct_neg    = round(100 * mean(dtt_to_cwt_treat < 0), 1),
    n_zero     = sum(dtt_to_cwt_treat == 0),
    median     = median(dtt_to_cwt_treat),
    p25        = quantile(dtt_to_cwt_treat, 0.25),
    p75        = quantile(dtt_to_cwt_treat, 0.75),
    max        = max(dtt_to_cwt_treat)
  ) %>% print()

# --- DTT to our first_tx_date ------------------------------------------------
cat("\nDTT to first_tx_date:\n")
cwt_validation %>%
  filter(!is.na(dtt_to_tx)) %>%
  summarise(
    n = n(),
    n_negative = sum(dtt_to_tx < 0),
    pct_neg    = round(100 * mean(dtt_to_tx < 0), 1),
    n_zero     = sum(dtt_to_tx == 0),
    median     = median(dtt_to_tx),
    p25        = quantile(dtt_to_tx, 0.25),
    p75        = quantile(dtt_to_tx, 0.75),
    max        = max(dtt_to_tx)
  ) %>% print()

# --- CWT treat date vs our first_tx_date -------------------------------------
cat("\nCWT treat date vs first_tx_date (days difference):\n")
cwt_validation %>%
  filter(!is.na(cwt_vs_first_tx)) %>%
  summarise(
    n = n(),
    n_exact_match = sum(cwt_vs_first_tx == 0),
    pct_exact     = round(100 * mean(cwt_vs_first_tx == 0), 1),
    n_within_5    = sum(abs(cwt_vs_first_tx) <= 5),
    pct_within_5  = round(100 * mean(abs(cwt_vs_first_tx) <= 5), 1),
    n_within_14   = sum(abs(cwt_vs_first_tx) <= 14),
    pct_within_14 = round(100 * mean(abs(cwt_vs_first_tx) <= 14), 1),
    median_diff   = median(cwt_vs_first_tx),
    p25           = quantile(cwt_vs_first_tx, 0.25),
    p75           = quantile(cwt_vs_first_tx, 0.75)
  ) %>% print()

# --- Negative DTT-to-treatment by pathway ------------------------------------
cat("\nNegative dtt_to_tx by pathway:\n")
cwt_validation %>%
  filter(!is.na(dtt_to_tx)) %>%
  mutate(dtt_issue = dtt_to_tx < 0) %>%
  group_by(tx_pathway) %>%
  summarise(n = n(), n_neg = sum(dtt_issue),
            pct_neg = round(100 * mean(dtt_issue), 1), .groups = "drop") %>%
  arrange(desc(pct_neg)) %>%
  print()

# =============================================================================
# 4. Merge DTT anchor onto cohort + derive DTT intervals and validity
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
    
    # DTT after diagnosis and at/before treatment (14-day tolerance)
    dtt_valid = !is.na(cwt_dtt_date) &
      wt_dx_to_dtt >= 0 &
      wt_dtt_to_tx >= -14,
    
    # EMR/ESD pathways: DTT less meaningful -> flag as NA
    dtt_valid = if_else(
      tx_pathway %in% c("EMR/ESD only", "EMR/ESD then surgery"),
      NA, dtt_valid
    )
  )

saveRDS(og_cohort, paste0(base_dir, "og_cohort_cwt_2015_2022.rds"))

# =============================================================================
# 5. Key interval summaries (where DTT valid)
# =============================================================================
cat("\nValidity by pathway:\n")
og_cohort %>%
  filter(!is.na(cwt_dtt_date)) %>%
  count(dtt_valid, tx_pathway) %>%
  group_by(tx_pathway) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  filter(!is.na(dtt_valid)) %>%
  arrange(tx_pathway) %>%
  print(n = 30)

cat("\nStaging (endo->DTT) and scheduling (DTT->tx) intervals:\n")
og_cohort %>%
  filter(dtt_valid == TRUE,
         !is.na(wt_endo_to_dtt), wt_endo_to_dtt >= 0,
         !is.na(wt_dtt_to_tx),   wt_dtt_to_tx >= 0) %>%
  summarise(
    n = n(),
    median_endo_dtt = median(wt_endo_to_dtt),
    p25_endo_dtt    = quantile(wt_endo_to_dtt, 0.25),
    p75_endo_dtt    = quantile(wt_endo_to_dtt, 0.75),
    median_dtt_tx   = median(wt_dtt_to_tx),
    p25_dtt_tx      = quantile(wt_dtt_to_tx, 0.25),
    p75_dtt_tx      = quantile(wt_dtt_to_tx, 0.75)
  ) %>% print()

# By deprivation: which component carries the gradient?
og_cohort %>%
  filter(dtt_valid == TRUE,
         !is.na(wt_endo_to_dtt), wt_endo_to_dtt >= 0,
         !is.na(wt_dtt_to_tx),   wt_dtt_to_tx >= 0,
         !is.na(NHSE_reversed_imd_quintile_lsoas)) %>%
  group_by(NHSE_reversed_imd_quintile_lsoas) %>%
  summarise(
    n = n(),
    median_endo_tx  = median(wt_endo_to_tx,  na.rm = TRUE),
    median_endo_dtt = median(wt_endo_to_dtt, na.rm = TRUE),
    median_dtt_tx   = median(wt_dtt_to_tx,   na.rm = TRUE),
    .groups = "drop"
  ) %>% print()

