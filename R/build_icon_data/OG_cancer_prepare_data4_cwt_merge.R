
library(tidyverse)
library(arrow)
library(haven)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(purrr)
library(lubridate)

ncras_og <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/ncras_og_2015_2022.rds")
ncras_og_ids <- ncras_og %>% distinct(pseudo_patientid) %>% pull()

og_cohort <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/og_cohort_2015_2022.rds")

# attach the RCS Charlson comorbidity lookup built in script 1d
og_cci <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/og_cci_2015_2022.rds")
og_cohort <- og_cohort %>%
  left_join(og_cci, by = "pseudo_patientid")

summary(as.factor(og_cohort$tx_pathway))
names(og_cohort)
# What CWT files do you have?
list.files(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/CWT/",
  full.names = FALSE
)

# Check structure - what ICD-10 field is used for site filtering?
# And what's the treatment modality coding?
cwt_test <- open_dataset(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/CWT/11_CWT_data_partitioned"
) %>%
  filter(site_icd10 %in% og_icd10) %>%
  collect() %>%
  head(1000)

names(cwt_test)

#cwt_test %>% View()

cwt_test %>%
  count(modality, sort = TRUE)

# Date fields
cwt_test %>%
  select(contains("date"), contains("Date")) %>%
  names()


# What CWT files do you have?
list.files(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/CWT/",
  full.names = FALSE
)

# Check structure - what ICD-10 field is used for site filtering?
# And what's the treatment modality coding?
cwt_test <- open_dataset(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/CWT/11_CWT_data_partitioned"
) %>%
  filter(site_icd10 %in% og_icd10) %>%
  collect() %>%
  head(1000)

names(cwt_test)

cwt_test %>%
  count(modality, sort = TRUE)

# Date fields
cwt_test %>%
  select(contains("date"), contains("Date")) %>%
  names()


# Standard CWT modality codes:
# 01 = Surgery
# 02 = Anti-cancer drug (chemo)
# 03 = Radiotherapy
# 04 = Concurrent chemoRT
# 05 = Other (includes active monitoring)
# 06 = Brachytherapy
# 07 = Surgery + drug
# 08 = Surgery + RT
# 09 = Surgery + chemoRT
# 23 = Endoscopic (EMR/ESD) -- post-2020
# 24 = Endoscopic + other

cwt_test %>%
  select(modality, treat_start, treat_period_start, 
         mdt_date, crtp_date, date_first_seen) %>%
  head(20)

# Date formats
cwt_test %>%
  select(crtp_date, treat_period_start, treat_start, 
         mdt_date, date_first_seen) %>%
  summary()

# How complete are the key date fields?
cwt_test %>%
  summarise(
    n                     = n(),
    pct_crtp              = round(100 * mean(!is.na(crtp_date) & 
                                               crtp_date != ""), 1),
    pct_first_seen        = round(100 * mean(!is.na(date_first_seen) & 
                                               date_first_seen != ""), 1),
    pct_dtt               = round(100 * mean(!is.na(treat_period_start) & 
                                               treat_period_start != ""), 1),
    pct_treat_start       = round(100 * mean(!is.na(treat_start) & 
                                               treat_start != ""), 1),
    pct_mdt               = round(100 * mean(!is.na(mdt_date) & 
                                               mdt_date != ""), 1)
  )


###################################################
# Read CWT for OG patients.
cwt_og <- open_dataset(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/CWT/11_CWT_data_partitioned"
) %>%
  filter(site_icd10 %in% og_icd10) %>%
  collect()  %>%
  mutate(
    pseudo_patientid   = as.character(pseudo_patientid),
    cwt_dtt_date       = as.Date(treat_period_start, format = "%d/%m/%Y"),
    cwt_treat_date     = as.Date(treat_start,        format = "%d/%m/%Y"),
    cwt_referral_date  = as.Date(crtp_date,          format = "%d/%m/%Y"),
    cwt_first_seen     = as.Date(date_first_seen,    format = "%d/%m/%Y"),
    cwt_mdt_date       = as.Date(mdt_date,           format = "%d/%m/%Y")
  ) %>%
  filter(pseudo_patientid %in% ncras_og_ids)

cat("CWT OG rows:", nrow(cwt_og), "\n")
cat("Patients:   ", n_distinct(cwt_og$pseudo_patientid), "\n")

# Records per patient
cwt_og %>%
  count(pseudo_patientid, sort = TRUE) %>%
  count(n, name = "n_patients") %>%
  print(n = 10)

cwt_og %>%
  group_by(modality) %>%
  summarise(n = n())


# Modality distribution for OG
cwt_og %>%
  count(modality, sort = TRUE)
names(cwt_og)
# Year coverage
cwt_og %>%
  mutate(year = year(cwt_treat_date)) %>%
  count(year, sort = FALSE) %>%
  print(n = 20)

cwt_anchor <- ncras_og %>%
  select(pseudo_patientid, diagmdy) %>%
  left_join(
    cwt_og %>%
      filter(
        !(modality %in% c("23","24") &
            cwt_treat_date < as.Date("2020-06-01")),
        !modality %in% c("97","98","99")
      ) %>%
      select(pseudo_patientid, cwt_dtt_date,
             cwt_treat_date, cwt_mdt_date, modality),
    by = "pseudo_patientid"
  ) %>%
  # Restrict to records within treatment window of diagnosis
  mutate(days_dx_to_dtt = as.integer(cwt_dtt_date - diagmdy)) %>%
  filter(
    !is.na(days_dx_to_dtt),
    days_dx_to_dtt >= -30,
    days_dx_to_dtt <= tx_window_days
  ) %>%
  # Keep earliest DTT within window per patient
  arrange(pseudo_patientid, cwt_dtt_date) %>%
  distinct(pseudo_patientid, .keep_all = TRUE) %>%
  select(pseudo_patientid, cwt_dtt_date, cwt_treat_date,
         cwt_mdt_date, modality, days_dx_to_dtt)

cat("CWT anchor patients:", nrow(cwt_anchor), "\n")
cat("DTT completeness:   ",
    round(100 * mean(!is.na(cwt_anchor$cwt_dtt_date)), 1), "%\n")
cat("MDT completeness:   ",
    round(100 * mean(!is.na(cwt_anchor$cwt_mdt_date)), 1), "%\n")

# DTT timing distribution
summary(cwt_anchor$days_dx_to_dtt)


################################################
# Validation checks.
################################################

# Join to og_cohort treatment dates for validation
cwt_validation <- cwt_anchor %>%
  left_join(
    og_cohort %>%
      select(pseudo_patientid, diagmdy, first_tx_date,
             surgery_date, sact_date, rt_date, tx_pathway),
    by = "pseudo_patientid"
  ) %>%
  mutate(
    # CWT internal consistency: DTT before treatment
    dtt_to_cwt_treat  = as.integer(cwt_treat_date - cwt_dtt_date),
    
    # CWT vs HES/SACT/RTDS treatment date comparison
    dtt_to_first_tx   = as.integer(first_tx_date - cwt_dtt_date),
    
    # CWT vs HES/SACT/RTDS treatment date comparison
    dtt_to_tx   = as.integer(first_tx_date - cwt_dtt_date),
    
    # CWT treat date vs our first_tx_date
    cwt_vs_first_tx     = as.integer(cwt_treat_date - first_tx_date)
  )

# --- DTT to CWT treatment date (internal CWT consistency) ------------------
cat("DTT to CWT treat date:\n")
cwt_validation %>%
  filter(!is.na(dtt_to_cwt_treat)) %>%
  summarise(
    n          = n(),
    n_negative = sum(dtt_to_cwt_treat < 0),
    pct_neg    = round(100 * mean(dtt_to_cwt_treat < 0), 1),
    n_zero     = sum(dtt_to_cwt_treat == 0),
    median     = median(dtt_to_cwt_treat),
    p25        = quantile(dtt_to_cwt_treat, 0.25),
    p75        = quantile(dtt_to_cwt_treat, 0.75),
    max        = max(dtt_to_cwt_treat)
  ) %>% print()

# --- DTT to our first treatment date ---------------------------------------
cat("\nDTT to our first_tx_date:\n")
cwt_validation %>%
  filter(!is.na(dtt_to_first_tx)) %>%
  summarise(
    n          = n(),
    n_negative = sum(dtt_to_first_tx < 0),
    pct_neg    = round(100 * mean(dtt_to_first_tx < 0), 1),
    n_zero     = sum(dtt_to_first_tx == 0),
    median     = median(dtt_to_first_tx),
    p25        = quantile(dtt_to_first_tx, 0.25),
    p75        = quantile(dtt_to_first_tx, 0.75),
    max        = max(dtt_to_first_tx)
  ) %>% print()

# --- CWT treat date vs our first_tx_date -----------------------------------
cat("\nCWT treat date vs our first_tx_date (days difference):\n")
cwt_validation %>%
  filter(!is.na(cwt_vs_first_tx)) %>%
  summarise(
    n              = n(),
    n_exact_match  = sum(cwt_vs_first_tx == 0),
    pct_exact      = round(100 * mean(cwt_vs_first_tx == 0), 1),
    n_within_5     = sum(abs(cwt_vs_first_tx) <= 5),
    pct_within_5   = round(100 * mean(abs(cwt_vs_first_tx) <= 5), 1),
    n_within_14    = sum(abs(cwt_vs_first_tx) <= 14),
    pct_within_14  = round(100 * mean(abs(cwt_vs_first_tx) <= 14), 1),
    median_diff    = median(cwt_vs_first_tx),
    p25            = quantile(cwt_vs_first_tx, 0.25),
    p75            = quantile(cwt_vs_first_tx, 0.75)
  ) %>% print()

# --- Negative DTT-to-treatment by pathway ----------------------------------
cat("\nNegative dtt_to_first_tx by pathway:\n")
cwt_validation %>%
  filter(!is.na(dtt_to_first_tx)) %>%
  mutate(dtt_issue = dtt_to_first_tx < 0) %>%
  group_by(tx_pathway) %>%
  summarise(
    n         = n(),
    n_neg     = sum(dtt_issue),
    pct_neg   = round(100 * mean(dtt_issue), 1),
    .groups   = "drop"
  ) %>%
  arrange(desc(pct_neg)) %>%
  print()

#################################################

# Add validity flag to cwt_anchor join
og_cohort <- og_cohort %>%
  left_join(
    cwt_anchor %>%
      select(pseudo_patientid, cwt_dtt_date, cwt_mdt_date, cwt_treat_date),
    by = "pseudo_patientid"
  ) %>%
  mutate(
    # DTT-based waiting time intervals
    wt_endo_to_dtt = as.integer(cwt_dtt_date - endoscopy_date),
    wt_dtt_to_tx   = as.integer(first_tx_date - cwt_dtt_date),
    wt_dx_to_dtt   = as.integer(cwt_dtt_date - diagmdy),
    
    # Validity flag: DTT must be after diagnosis and before/on treatment
    # Allow 14-day tolerance for minor date discrepancies
    dtt_valid = !is.na(cwt_dtt_date) &
      wt_dx_to_dtt >= 0 &
      wt_dtt_to_tx >= -14,
    
    # For EMR/ESD pathways DTT is less meaningful - flag separately
    dtt_valid = if_else(
      tx_pathway %in% c("EMR/ESD only", "EMR/ESD then surgery"),
      NA,
      dtt_valid
    )
  )


saveRDS(og_cohort,
        "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/og_cohort_cwt_2015_2022.rds")


names(og_cohort)

# Check validity
og_cohort %>%
  filter(!is.na(cwt_dtt_date)) %>%
  count(dtt_valid, tx_pathway) %>%
  group_by(tx_pathway) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  filter(!is.na(dtt_valid)) %>%
  arrange(tx_pathway) %>%
  print(n = 30)

#################################################

# Check the key intervals where dtt_valid
og_cohort %>%
  filter(dtt_valid == TRUE,
         !is.na(wt_endo_to_dtt), wt_endo_to_dtt >= 0,
         !is.na(wt_dtt_to_tx),   wt_dtt_to_tx >= 0) %>%
  summarise(
    n               = n(),
    # Endoscopy to DTT (staging interval)
    median_endo_dtt = median(wt_endo_to_dtt),
    p25_endo_dtt    = quantile(wt_endo_to_dtt, 0.25),
    p75_endo_dtt    = quantile(wt_endo_to_dtt, 0.75),
    # DTT to treatment (scheduling interval)
    median_dtt_tx   = median(wt_dtt_to_tx),
    p25_dtt_tx      = quantile(wt_dtt_to_tx, 0.25),
    p75_dtt_tx      = quantile(wt_dtt_to_tx, 0.75)
  )

# By pathway
og_cohort %>%
  filter(dtt_valid == TRUE,
         !is.na(wt_endo_to_dtt), wt_endo_to_dtt >= 0,
         !is.na(wt_dtt_to_tx),   wt_dtt_to_tx >= 0) %>%
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
  arrange(desc(n))


#####################################
# Does deprivation affect staging interval, scheduling interval, or both?
og_cohort %>%
  filter(dtt_valid == TRUE,
         !is.na(wt_endo_to_dtt), wt_endo_to_dtt >= 0,
         !is.na(wt_dtt_to_tx),   wt_dtt_to_tx >= 0,
         !is.na(NHSE_reversed_imd_quintile_lsoas)) %>%
  group_by(NHSE_reversed_imd_quintile_lsoas) %>%
  summarise(
    n               = n(),
    # Full pathway
    median_endo_tx  = median(wt_endo_to_tx,  na.rm = TRUE),
    # Staging component
    median_endo_dtt = median(wt_endo_to_dtt, na.rm = TRUE),
    # Scheduling component
    median_dtt_tx   = median(wt_dtt_to_tx,   na.rm = TRUE),
    .groups         = "drop"
  )



####################################################################

# Does staging gradient persist within TWW?

og_cohort %>%
  filter(dtt_valid == TRUE,
         !is.na(wt_endo_to_dtt), wt_endo_to_dtt > 0,
         !is.na(wt_dtt_to_tx),   wt_dtt_to_tx > 0,
         !is.na(NHSE_reversed_imd_quintile_lsoas),
         as.character(route_combined) == "TWW") %>%
  group_by(NHSE_reversed_imd_quintile_lsoas) %>%
  summarise(
    n               = n(),
    median_endo_tx  = median(wt_endo_to_tx,  na.rm = TRUE),
    median_endo_dtt = median(wt_endo_to_dtt, na.rm = TRUE),
    median_dtt_tx   = median(wt_dtt_to_tx,   na.rm = TRUE),
    .groups         = "drop"
  )

# And non-TWW routes
og_cohort %>%
  filter(dtt_valid == TRUE,
         !is.na(wt_endo_to_dtt), wt_endo_to_dtt > 0,
         !is.na(wt_dtt_to_tx),   wt_dtt_to_tx > 0,
         !is.na(NHSE_reversed_imd_quintile_lsoas),
         as.character(route_combined) == "GP referral") %>%
  group_by(NHSE_reversed_imd_quintile_lsoas) %>%
  summarise(
    n               = n(),
    median_endo_tx  = median(wt_endo_to_tx,  na.rm = TRUE),
    median_endo_dtt = median(wt_endo_to_dtt, na.rm = TRUE),
    median_dtt_tx   = median(wt_dtt_to_tx,   na.rm = TRUE),
    .groups         = "drop"
  )





