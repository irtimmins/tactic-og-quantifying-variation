

ncras_og <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/ncras_og_2015_2022.rds")
ncras_og_ids <- ncras_og %>% distinct(pseudo_patientid) %>% pull()
#names(ncras_og)
# Read in data from previous script.

og_cohort <- readRDS(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/og_cohort_ncras_hes_2015_2022.rds"
)


surgery_anchor <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/og_surgery_anchor_2015_2022.rds")



# =============================================================================
# 8. SACT (SYSTEMIC ANTI-CANCER THERAPY) - PLACEHOLDER
#    Source: SACT dataset (NDRS) linked via pseudo_patientid / STUDY_ID
#    NOGCA Appendices 9-10; benchmark_group classifies palliative regimens
#
#    Key SACT variables:
#      start_date_of_cycle  - cycle start date (used for tx date; PI9/PI10)
#      benchmark_group      - regimen classification (Appendix 10)
#      primary_diagnosis    - confirm C15 or C16
#      morphology_clean     - fallback epithelial filter
#
#    Palliative regimens (Appendix 10):
#      1. Immunotherapy: pembrolizumab or nivolumab (alone or in combination)
#      2. Trastuzumab-containing regimens
#      3. Triplet regimens: cisplatin/capecitabine/epirubicin combinations
#      4. Doublet platinum+5FU: oxaliplatin/cisplatin + fluorouracil/
#         capecitabine/tegafur; carboplatin+paclitaxel
#
#    GIST exclusions: imatinib, sunitinib, regorafenib
#    NE exclusions:   carboplatin/cisplatin + etoposide/sunitinib/everolimus
# =============================================================================



library(readr)
test <- read_csv("E:/Data_PHE/Raw data files received from PHE READ ONLY/SACT/6_SACT_data_treat_yr_2024.csv",n_max =  100000)
names(test)
#test
View(test)
# How many files and what years?
list.files(
  "E:/Data_PHE/Raw data files received from PHE READ ONLY/SACT/",
  pattern = "*.csv"
)

# Check the parsing issues
problems(test)

# Check all column names before deciding what to keep
names(test)

# How many OG records in this one file?
test %>%
  filter(str_starts(PRIMARY_DIAGNOSIS, "C15") |
           str_starts(PRIMARY_DIAGNOSIS, "C16")) %>%
  nrow()

# Check the fields you need
test %>%
  select(
    PSEUDO_PATIENTID,
    PRIMARY_DIAGNOSIS,
    MORPHOLOGY_CLEAN,
    BENCHMARK_GROUP,
    ANALYSIS_GROUP,
    INTENT_OF_TREATMENT,
    INTENT_OF_TREATMENT_V3,
    START_DATE_OF_REGIMEN,
    DATE_DECISION_TO_TREAT,
    ORGANISATION_CODE_OF_PROVIDER,
    STAGE_AT_START
  ) %>%
  glimpse()

# What OG diagnoses look like
test %>%
  filter(str_starts(PRIMARY_DIAGNOSIS, "C15") |
           str_starts(PRIMARY_DIAGNOSIS, "C16")) %>%
  count(PRIMARY_DIAGNOSIS, BENCHMARK_GROUP, sort = TRUE) %>%
  print(n = 20)

# Intent of treatment codes
test %>%
  count(INTENT_OF_TREATMENT, INTENT_OF_TREATMENT_V3, sort = TRUE)

# Date format
test %>%
  select(START_DATE_OF_REGIMEN, DATE_DECISION_TO_TREAT) %>%
  head(10)


sact_cols <- c(
  "PSEUDO_PATIENTID",
  "PRIMARY_DIAGNOSIS",
  "MORPHOLOGY_CLEAN",
  "BENCHMARK_GROUP",
  "ANALYSIS_GROUP",
  "INTENT_OF_TREATMENT_V3",
  "START_DATE_OF_REGIMEN",
  "START_DATE_OF_CYCLE",
  "DATE_DECISION_TO_TREAT",
  "DATE_OF_FINAL_TREATMENT",
  "CYCLE_NUMBER",
  "NUMBER_OF_CYCLES_PLANNED",
  "ORGANISATION_CODE_OF_PROVIDER",
  "PERF_STAT_START_OF_REG_ADULT",
  "STAGE_AT_START",
  "REGIMEN_MOD_STOPPED_EARLY",
  "REGIMEN_OUTCOME_SUMMARY",
  "CHEMO_RADIATION"
)

sact_file_list <- list.files(
  "E:/Data_PHE/Raw data files received from PHE READ ONLY/SACT/",
  pattern    = "*.csv",
  full.names = TRUE
)

sact_og <- map_dfr(
  sact_file_list,
  ~{
    read_csv(
      .x,
      col_select  = all_of(sact_cols),
      col_types   = cols(.default = col_character()),
      show_col_types = FALSE
    ) %>%
      filter(
        str_starts(PRIMARY_DIAGNOSIS, "C15") |
          str_starts(PRIMARY_DIAGNOSIS, "C16")
      ) %>%
      mutate(
        pseudo_patientid    = as.character(as.integer(PSEUDO_PATIENTID)),
        sact_regimen_date   = as.Date(START_DATE_OF_REGIMEN, format = "%d/%m/%Y"),
        sact_cycle_date     = as.Date(START_DATE_OF_CYCLE,   format = "%d/%m/%Y"),
        date_decision_treat = as.Date(DATE_DECISION_TO_TREAT, format = "%d/%m/%Y"),
        date_final_tx       = as.Date(DATE_OF_FINAL_TREATMENT, format = "%d/%m/%Y"),
        cycle_number        = as.integer(CYCLE_NUMBER),
        benchmark_group_lwr = tolower(trimws(BENCHMARK_GROUP))
      )
  },
  .progress = TRUE
) %>%
  filter(pseudo_patientid %in% ncras_og_ids)

cat("SACT OG rows (cohort patients):", nrow(sact_og), "\n")
cat("Patients with any SACT:        ", n_distinct(sact_og$pseudo_patientid), "\n")

# Quick look at regimens
cat("\nTop BENCHMARK_GROUP values:\n")
sact_og %>%
  count(BENCHMARK_GROUP, sort = TRUE) %>%
  print(n = 20)

# Intent distribution
cat("\nINTENT_OF_TREATMENT_V3:\n")
count(sact_og, INTENT_OF_TREATMENT_V3, sort = TRUE)

# Date range
cat("\nSACT regimen date range:\n")
summary(sact_og$sact_regimen_date)

saveRDS(
  sact_og,
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/sact_og_2012_2024.rds"
)

sact_og <- readRDS( "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/sact_og_2012_2024.rds")
#View(sact_og[1:100,])

# --- GIST and NE exclusions (NOGCA ?3) --------------------------------------
gist_regimens <- c("imatinib", "sunitinib", "regorafenib")

ne_regimens <- c(
  "carboplatin + etoposide", "cisplatin + etoposide",
  "carboplatin + sunitinib", "cisplatin + sunitinib",
  "carboplatin + everolimus", "cisplatin + everolimus"
)

# --- First SACT date per patient within treatment window --------------------
sact_anchor <- ncras_og %>%
  select(pseudo_patientid, diagmdy) %>%
  left_join(
    sact_og %>%
      # Exclude GIST and NE regimens
      filter(
        !benchmark_group_lwr %in% gist_regimens,
        !benchmark_group_lwr %in% ne_regimens,
        !BENCHMARK_GROUP %in% c("NOT CHEMO", "TRIAL UNSPECIFIED")
      ) %>%
      select(pseudo_patientid, sact_regimen_date, sact_cycle_date,
             BENCHMARK_GROUP, benchmark_group_lwr,
             INTENT_OF_TREATMENT_V3, CYCLE_NUMBER, cycle_number,
             ORGANISATION_CODE_OF_PROVIDER, CHEMO_RADIATION),
    by = "pseudo_patientid"
  ) %>%
  mutate(days_dx_to_sact = as.integer(sact_regimen_date - diagmdy)) %>%
  filter(
    !is.na(days_dx_to_sact),
    days_dx_to_sact >= -30,        # small pre-diagnosis window
    days_dx_to_sact <= tx_window_days
  ) %>%
  arrange(pseudo_patientid, sact_regimen_date) %>%
  distinct(pseudo_patientid, .keep_all = TRUE) %>%
  rename(sact_date = sact_regimen_date) %>%
  select(pseudo_patientid, sact_date, days_dx_to_sact,
         BENCHMARK_GROUP, benchmark_group_lwr,
         INTENT_OF_TREATMENT_V3, CHEMO_RADIATION,
         ORGANISATION_CODE_OF_PROVIDER)

saveRDS(sact_anchor, "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/og_sact_anchor_2015_2022.rds")


cat("Patients with SACT in treatment window:", 
    n_distinct(sact_anchor$pseudo_patientid), "\n")

# Date distribution
summary(sact_anchor$sact_date)
summary(sact_anchor$days_dx_to_sact)

# How many patients have both SACT and surgery?
both <- sum(sact_anchor$pseudo_patientid %in% surgery_anchor$pseudo_patientid)
sact_only <- sum(!sact_anchor$pseudo_patientid %in% surgery_anchor$pseudo_patientid)
surg_only <- sum(!surgery_anchor$pseudo_patientid %in% sact_anchor$pseudo_patientid)

cat("Both SACT and surgery:  ", both, "\n")
cat("SACT only:              ", sact_only, "\n")
cat("Surgery only:           ", surg_only, "\n")

# Of those with both - how many had SACT before surgery (neoadjuvant)?
ncras_og %>%
  select(pseudo_patientid) %>%
  left_join(sact_anchor %>% select(pseudo_patientid, sact_date), 
            by = "pseudo_patientid") %>%
  left_join(surgery_anchor %>% select(pseudo_patientid, surgery_date),
            by = "pseudo_patientid") %>%
  filter(!is.na(sact_date), !is.na(surgery_date)) %>%
  mutate(
    sact_before_surgery = sact_date < surgery_date,
    days_sact_to_surg   = as.integer(surgery_date - sact_date)
  ) %>%
  summarise(
    n                   = n(),
    n_neoadjuvant       = sum(sact_before_surgery),
    n_adjuvant          = sum(!sact_before_surgery),
    median_sact_to_surg = median(days_sact_to_surg[sact_before_surgery], na.rm = TRUE),
    p25                 = quantile(days_sact_to_surg[sact_before_surgery], 0.25, na.rm = TRUE),
    p75                 = quantile(days_sact_to_surg[sact_before_surgery], 0.75, na.rm = TRUE)
  ) %>%
  print()

# =============================================================================
# 9. RADIOTHERAPY (RTDS) - PLACEHOLDER
#    Source: RTDS (Radiotherapy Dataset, NDRS)
#    Filter: radiotherapydiagnosisicd %in% og_icd10
#    Tx date: apptdate of first RT prescription within treatment window
#
#    Curative RT fractionation schedules (NOGCA Table 4.3 definition):
#      50 Gy / 25 fractions
#      50.4 Gy / 28 fractions
#      60 Gy / 30 fractions
#      50 Gy / 15 or 16 fractions
#      50-55 Gy / 20 fractions
#      45-52.5 Gy / 15 or 16 fractions
#    All other RT = palliative
# =============================================================================


# What RTDS files are available?
list.files(
  "E:/Data_PHE/Raw data files received from PHE READ ONLY/RTDS/",
  full.names = FALSE
)

# Read a sample
rtds_test <- read_csv(
  # or read_dta if Stata format - adjust path as needed
  "E:/Data_PHE/Raw data files received from PHE READ ONLY/RTDS/7_RTDS_data_treat_yr_2023.csv",
  n_max = 10000
)

names(rtds_test)
#View(rtds_test)
# Key fields to check
rtds_test %>%
  select(any_of(c(
    "PSEUDO_PATIENTID", "STUDY_ID",
    "RADIOTHERAPYDIAGNOSISICD",
    "APPTDATE", "PRESCRIBEDDOSE", "PRESCRIBEDFRACTIONS",
    "STARTDATEOFRADIOTHERAPY", "MANDATE_ID"
  ))) %>%
  glimpse()

# Intent values
rtds_test %>%
  count(RADIOTHERAPYINTENT, sort = TRUE)

# Dose/fractions distribution for OG
rtds_test %>%
  filter(str_starts(RADIOTHERAPYDIAGNOSISICD, "C15") |
           str_starts(RADIOTHERAPYDIAGNOSISICD, "C16")) %>%
  count(RTPRESCRIBEDDOSE, PRESCRIBEDFRACTIONS, sort = TRUE) %>%
  print(n = 20)

# Date format
rtds_test %>%
  select(APPTDATE, TREATMENTSTARTDATE, DECISIONTOTREATDATE) %>%
  head(10)

# How many OG rows in this file?
rtds_test %>%
  filter(str_starts(RADIOTHERAPYDIAGNOSISICD, "C15") |
           str_starts(RADIOTHERAPYDIAGNOSISICD, "C16")) %>%
  nrow()


rtds_cols <- c(
  "PSEUDO_PATIENTID",
  "RADIOTHERAPYDIAGNOSISICD",
  "TREATMENTSTARTDATE",
  "DECISIONTOTREATDATE",
  "RADIOTHERAPYINTENT",
  "RTPRESCRIBEDDOSE",
  "PRESCRIBEDFRACTIONS",
  "RTTREATMENTREGION",
  "RTTREATMENTANATOMICALSITE",
  "ORGCODEPROVIDER"
)

rtds_file_list <- list.files(
  "E:/Data_PHE/Raw data files received from PHE READ ONLY/RTDS/",
  pattern    = "*.csv",
  full.names = TRUE
)

rtds_og <- map_dfr(
  rtds_file_list,
  ~{
    read_csv(
      .x,
      col_select     = all_of(rtds_cols),
      col_types      = cols(.default = col_character()),
      show_col_types = FALSE
    ) %>%
      filter(
        str_starts(RADIOTHERAPYDIAGNOSISICD, "C15") |
          str_starts(RADIOTHERAPYDIAGNOSISICD, "C16")
      ) %>%
      mutate(
        pseudo_patientid = as.character(as.integer(as.numeric(PSEUDO_PATIENTID))),
        rt_start_date    = as.Date(TREATMENTSTARTDATE,   format = "%d/%m/%Y"),
        rt_decision_date = as.Date(DECISIONTOTREATDATE,  format = "%d/%m/%Y"),
        rt_dose          = as.numeric(RTPRESCRIBEDDOSE),
        rt_fractions     = as.integer(PRESCRIBEDFRACTIONS)
      )
  },
  .progress = TRUE
) %>%
  filter(pseudo_patientid %in% ncras_og_ids)

saveRDS(rtds_og,
        "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/rtds_og_2009_2024.rds")
rtds_og <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/rtds_og_2009_2024.rds")

cat("RTDS OG rows (cohort patients):", nrow(rtds_og), "\n")
cat("Patients with any RT:          ", n_distinct(rtds_og$pseudo_patientid), "\n")

# Intent distribution
cat("\nRADIOTHERAPYINTENT:\n")
count(rtds_og, RADIOTHERAPYINTENT, sort = TRUE)

# Dose/fractions for curative intent
cat("\nDose/fractions (curative intent = 01):\n")
rtds_og %>%
  filter(RADIOTHERAPYINTENT == "01") %>%
  count(rt_dose, rt_fractions, sort = TRUE) %>%
  print(n = 20)


# Build RT anchor using dose/fractionation for curative classification
rt_anchor <- ncras_og %>%
  select(pseudo_patientid, diagmdy) %>%
  left_join(
    rtds_og %>%
      arrange(pseudo_patientid, rt_start_date) %>%
      distinct(pseudo_patientid, rt_start_date, .keep_all = TRUE) %>%
      mutate(
        rt_curative = case_when(
          rt_dose == 50    & rt_fractions == 25           ~ TRUE,
          rt_dose == 50.4  & rt_fractions == 28           ~ TRUE,
          rt_dose == 60    & rt_fractions == 30           ~ TRUE,
          rt_dose == 50    & rt_fractions %in% c(15, 16) ~ TRUE,
          between(rt_dose, 50, 55) & rt_fractions == 20  ~ TRUE,
          between(rt_dose, 45, 52.5) &
            rt_fractions %in% c(15, 16)                  ~ TRUE,
          rt_dose == 41.4  & rt_fractions == 23           ~ TRUE,
          rt_dose == 45    & rt_fractions == 25           ~ TRUE,
          TRUE                                            ~ FALSE
        )
      ) %>%
      select(pseudo_patientid, rt_start_date, rt_curative,
             rt_dose, rt_fractions, ORGCODEPROVIDER),
    by = "pseudo_patientid"
  ) %>%
  mutate(days_dx_to_rt = as.integer(rt_start_date - diagmdy)) %>%
  filter(
    !is.na(days_dx_to_rt),
    days_dx_to_rt >= -30,
    days_dx_to_rt <= tx_window_days
  ) %>%
  # Prioritise curative RT course, then earliest within each group
  arrange(pseudo_patientid, desc(rt_curative), rt_start_date) %>%
  distinct(pseudo_patientid, .keep_all = TRUE) %>%
  rename(rt_date = rt_start_date) %>%
  select(pseudo_patientid, rt_date, rt_curative,
         rt_dose, rt_fractions, days_dx_to_rt, ORGCODEPROVIDER)

cat("Patients with RT in treatment window:", 
    n_distinct(rt_anchor$pseudo_patientid), "\n")
count(rt_anchor, rt_curative)

cat("Patients with RT in treatment window:", 
    n_distinct(rt_anchor$pseudo_patientid), "\n")

cat("\nCurative vs palliative RT:\n")
count(rt_anchor, rt_curative)

cat("\nCurative RT dose/fractions:\n")
rt_anchor %>%
  filter(rt_curative) %>%
  count(rt_dose, rt_fractions, sort = TRUE)

# Do any patients have both curative and palliative RT?
rtds_og %>%
  filter(pseudo_patientid %in% ncras_og_ids) %>%
  mutate(
    rt_curative = case_when(
      as.numeric(RTPRESCRIBEDDOSE) == 50   & as.integer(PRESCRIBEDFRACTIONS) == 25  ~ TRUE,
      as.numeric(RTPRESCRIBEDDOSE) == 41.4 & as.integer(PRESCRIBEDFRACTIONS) == 23  ~ TRUE,
      as.numeric(RTPRESCRIBEDDOSE) == 45   & as.integer(PRESCRIBEDFRACTIONS) == 25  ~ TRUE,
      between(as.numeric(RTPRESCRIBEDDOSE), 50, 55) & 
        as.integer(PRESCRIBEDFRACTIONS) == 20                                        ~ TRUE,
      as.numeric(RTPRESCRIBEDDOSE) == 50.4 & as.integer(PRESCRIBEDFRACTIONS) == 28  ~ TRUE,
      as.numeric(RTPRESCRIBEDDOSE) == 60   & as.integer(PRESCRIBEDFRACTIONS) == 30  ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  distinct(pseudo_patientid, rt_start_date, rt_curative) %>%
  group_by(pseudo_patientid) %>%
  summarise(
    has_curative  = any(rt_curative),
    has_palliative = any(!rt_curative)
  ) %>%
  filter(has_curative & has_palliative) %>%
  nrow()



saveRDS(rt_anchor,
        "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/rt_anchor_og.rds")

rtds_og <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/rtds_og_2009_2024.rds")
rt_anchor <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/rt_anchor_og.rds")
names(rtds_og)
sort(summary(as.factor(rtds_og$ORGCODEPROVIDER)))


# =============================================================================
# 10. BUILD OG_COHORT
# Join all treatment anchors to NCRAS anchor
# Derive treatment pathway classification and waiting time components
# =============================================================================

endoscopy_anchor_combined <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/OG_endoscopy_anchor_combined.rds")
emresd_anchor <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/OG_emresd_anchor.rds")

# =============================================================================
# 9b. HES chemotherapy anchor (supplements SACT)
# SACT misses chemo delivery for a material group (coverage gaps, trial
# regimens, earlier years); HES APC records the delivery episode via OPCS
# X70-X74 and ICD-10 Z51.1. HES is a supplement, not a replacement: SACT stays
# primary (it carries regimen, intent, benchmark group), HES fills where SACT
# is absent, and the chemo date is SACT-preferred. chemo_source (sact/hes/both)
# keeps provenance auditable; HES-only chemo carries no benchmark/intent.
# =============================================================================
hes_apc <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/hes_apc_og_2014_2022.rds")

opcs_chemo_delivery <- c("X701","X702","X703","X704","X709",
                         "X711","X712","X713","X714","X719",
                         "X721","X722","X723","X724","X725","X726","X727",
                         "X728","X729",
                         "X731","X732","X738","X739",
                         "X741","X742","X743","X744","X748","X749")
icd_chemo_attendance <- c("Z511")

hes_sub <- hes_apc %>%
  mutate(STUDY_ID = as.character(STUDY_ID)) %>%
  filter(STUDY_ID %in% ncras_og_ids)

hes_opcs <- hes_sub %>%
  select(STUDY_ID, EPISTART, starts_with("OPERTN_")) %>%
  pivot_longer(starts_with("OPERTN_"), names_to = "pos",
               values_to = "opcs", names_prefix = "OPERTN_") %>%
  left_join(
    hes_sub %>%
      select(STUDY_ID, EPISTART, starts_with("OPDATE_")) %>%
      pivot_longer(starts_with("OPDATE_"), names_to = "pos",
                   values_to = "opdate", names_prefix = "OPDATE_"),
    by = c("STUDY_ID", "EPISTART", "pos"), relationship = "many-to-many"
  ) %>%
  mutate(opcs4 = str_to_upper(str_remove_all(str_trim(opcs), "\\."))) %>%
  filter(opcs4 %in% opcs_chemo_delivery) %>%
  transmute(STUDY_ID, chemo_date = coalesce(as.Date(opdate), as.Date(EPISTART)))

hes_icd <- hes_sub %>%
  select(STUDY_ID, EPISTART, starts_with("DIAG_4_")) %>%
  pivot_longer(starts_with("DIAG_4_"), names_to = "pos",
               values_to = "icd", names_prefix = "DIAG_4_") %>%
  mutate(icd4 = str_to_upper(str_sub(str_remove_all(str_trim(icd), "\\."), 1, 4))) %>%
  filter(icd4 %in% icd_chemo_attendance) %>%
  transmute(STUDY_ID, chemo_date = as.Date(EPISTART))

hes_chemo_anchor <- bind_rows(hes_opcs, hes_icd) %>%
  filter(!is.na(chemo_date)) %>%
  rename(pseudo_patientid = STUDY_ID) %>%
  inner_join(ncras_og %>% select(pseudo_patientid, diagmdy),
             by = "pseudo_patientid") %>%
  mutate(days_dx_to_hes_chemo = as.integer(chemo_date - diagmdy)) %>%
  filter(days_dx_to_hes_chemo >= -30, days_dx_to_hes_chemo <= tx_window_days) %>%
  arrange(pseudo_patientid, chemo_date) %>%
  distinct(pseudo_patientid, .keep_all = TRUE) %>%
  transmute(pseudo_patientid, hes_chemo_date = chemo_date, days_dx_to_hes_chemo)

cat("Patients with in-window HES chemo:", nrow(hes_chemo_anchor), "\n")

# combine SACT and HES: SACT-preferred date, provenance kept
chemo_anchor <- ncras_og %>%
  select(pseudo_patientid, diagmdy) %>%
  left_join(sact_anchor %>% select(pseudo_patientid, sact_date, days_dx_to_sact,
                                   BENCHMARK_GROUP, benchmark_group_lwr,
                                   INTENT_OF_TREATMENT_V3, CHEMO_RADIATION,
                                   ORGANISATION_CODE_OF_PROVIDER),
            by = "pseudo_patientid") %>%
  left_join(hes_chemo_anchor, by = "pseudo_patientid") %>%
  filter(!is.na(sact_date) | !is.na(hes_chemo_date)) %>%
  mutate(
    chemo_source     = case_when(
      !is.na(sact_date) &  !is.na(hes_chemo_date) ~ "both",
      !is.na(sact_date) &   is.na(hes_chemo_date) ~ "sact",
      TRUE                                        ~ "hes"
    ),
    chemo_date       = coalesce(sact_date, hes_chemo_date),
    days_dx_to_chemo = as.integer(chemo_date - diagmdy)
  )

cat("Chemo anchor source split:\n")
chemo_anchor %>% count(chemo_source) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>% print()
saveRDS(chemo_anchor,
        "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/og_chemo_anchor_2015_2022.rds")


og_cohort <- ncras_og %>%
  
  # --- Diagnostic endoscopy (PI3 clock start) --------------------------------
left_join(
  endoscopy_anchor_combined %>%
    select(pseudo_patientid, endoscopy_date, days_endo_to_dx),
  by = "pseudo_patientid"
) %>%
  
  # --- EMR/ESD ---------------------------------------------------------------
left_join(
  emresd_anchor %>%
    select(pseudo_patientid, emresd_date, days_dx_to_emresd),
  by = "pseudo_patientid"
) %>%
  
  # --- Surgery ---------------------------------------------------------------
left_join(
  surgery_anchor %>%
    select(pseudo_patientid, surgery_date, surgery_type, surgery_class,
           opcs_primary, PROCODE3, SITETRET, days_dx_to_surg, curative_surgery),
  by = "pseudo_patientid"
) %>%
  
  # --- SACT (+ HES chemo supplement) -----------------------------------------
# chemo_anchor combines the SACT anchor with HES delivery codes (OPCS X70-74,
# ICD-10 Z51.1); see the HES chemo anchor block above. sact_date here is the
# SACT-preferred combined chemo date, had_sact below becomes SACT-or-HES, and
# chemo_source records provenance. hes_chemo_date is carried so the curative-RT
# concurrency guard can require HES-only chemo to sit near the RT before it
# reclassifies a patient to definitive chemoRT.
left_join(
  chemo_anchor %>%
    select(pseudo_patientid,
           sact_date = chemo_date, days_dx_to_sact = days_dx_to_chemo,
           chemo_source, hes_chemo_date,
           BENCHMARK_GROUP, benchmark_group_lwr,
           INTENT_OF_TREATMENT_V3, CHEMO_RADIATION,
           ORGANISATION_CODE_OF_PROVIDER),
  by = "pseudo_patientid"
) %>%
  
  # --- RT --------------------------------------------------------------------
left_join(
  rt_anchor %>%
    select(pseudo_patientid, rt_date, rt_curative,
           rt_dose, rt_fractions, days_dx_to_rt, ORGCODEPROVIDER),
  by = "pseudo_patientid"
) %>%
  # --- Treatment presence flags ----------------------------------------------
mutate(
  had_emresd           = !is.na(emresd_date),
  had_surgery          = !is.na(surgery_date),
  had_curative_surgery = !is.na(surgery_date) & curative_surgery == TRUE,
  had_sact             = !is.na(sact_date),
  had_rt               = !is.na(rt_date),
  had_curative_rt      = !is.na(rt_date) & rt_curative == TRUE,
  had_palliative_rt    = !is.na(rt_date) & rt_curative == FALSE,
  
  # chemo provenance ("sact"/"hes"/"both") arrives from chemo_anchor via the
  # join above; HES-only chemo carries no benchmark/intent, so regimen-specific
  # analysis should exclude chemo_source == "hes".
  # chemo eligible to define a non-surgical definitive-chemoRT pathway. SACT
  # chemo always counts; HES-only chemo counts only when it sits within 28 days
  # of the curative RT, so a temporally separate HES chemo episode does not
  # manufacture definitive chemoRT out of a curative-RT-only patient. The
  # surgery and no-treatment branches use had_sact directly and are unaffected.
  had_chemo_for_chemort = had_sact &
    ( coalesce(chemo_source, "sact") != "hes" |
        ( !is.na(hes_chemo_date) & !is.na(rt_date) &
            abs(as.integer(hes_chemo_date - rt_date)) <= 28 ) ),
  
  # --- Treatment sequencing flags ------------------------------------------
  sact_before_surgery = had_sact & had_surgery & sact_date < surgery_date,
  sact_after_surgery  = had_sact & had_surgery & sact_date > surgery_date,
  rt_before_surgery   = had_rt   & had_surgery & rt_date   < surgery_date,
  rt_after_surgery    = had_rt   & had_surgery & rt_date   > surgery_date,
  
  # Chemo-RT flag: SACT and curative RT within 14 days of each other
  concurrent_chemo_rt = had_sact & had_curative_rt &
    abs(as.integer(sact_date - rt_date)) <= 14,
  
  # --- Curative treatment received -----------------------------------------
  received_curative_tx = had_emresd | had_curative_surgery | had_curative_rt,
  
  # --- Treatment pathway classification ------------------------------------
  tx_pathway = case_when(
    # EMR/ESD pathways
    had_emresd & !had_surgery & !had_sact & !concurrent_chemo_rt
    ~ "EMR/ESD only",
    had_emresd & had_surgery                      ~ "EMR/ESD then surgery",
    
    # Surgical pathways
    had_surgery & sact_before_surgery & rt_before_surgery
    ~ "Surgery + neoadjuvant chemoRT",
    had_surgery & sact_before_surgery & !rt_before_surgery
    ~ "Surgery + neoadjuvant chemo",
    had_surgery & rt_before_surgery & !sact_before_surgery
    ~ "Surgery + neoadjuvant RT",
    had_surgery & sact_after_surgery & !sact_before_surgery
    ~ "Surgery + adjuvant chemo",
    had_surgery & !had_sact & !concurrent_chemo_rt ~ "Surgery only",
    had_surgery                                   ~ "Surgery + other",
    
    # Non-surgical curative. Definitive chemoRT requires chemo that is part of
    # the RT course: SACT chemo, or HES-only chemo within 28d of the RT (the
    # had_chemo_for_chemort guard). A curative-RT patient whose only chemo is a
    # temporally separate HES episode is "Curative RT only" here (the chemo is
    # kept as a had_sact/chemo_source flag but does not define the pathway),
    # rather than mislabelled definitive chemoRT.
    !had_surgery & had_curative_rt & had_chemo_for_chemort  ~ "Definitive chemoRT",
    !had_surgery & had_curative_rt & !had_chemo_for_chemort ~ "Curative RT only",
    
    # Palliative / non-curative
    !had_surgery & had_palliative_rt & had_sact   ~ "Palliative chemo + RT",
    !had_surgery & had_sact & !had_curative_rt    ~ "SACT only",
    !had_surgery & had_palliative_rt & !had_sact  ~ "Palliative RT only",
    
    TRUE                                          ~ "No treatment recorded"
  ),
  
  # --- First disease-targeted treatment date -------------------------------
  # first_tx_date = pmin(emresd_date, surgery_date, sact_date, rt_date,
  #                      na.rm = TRUE),
  
  # --- First curative treatment date ---------------------------------------
  # Curative modalities are: endoscopic resection, major surgery,
  # curative RT, and definitive chemoRT.
  # Neoadjuvant RT/chemoRT prior to surgery sets the date (curative intent
  # starts with RT); neoadjuvant chemo alone does not (surgery sets the date).
  first_tx_date = case_when(
    
    # Endoscopic resection
    tx_pathway %in% c("EMR/ESD only", "EMR/ESD then surgery")
    ~ emresd_date,
    
    # Neoadjuvant chemoRT then surgery: curative RT/chemoRT sets the date
    tx_pathway == "Surgery + neoadjuvant chemoRT" ~ pmin(sact_date, rt_date,
                                                         na.rm = TRUE),
    
    # Neoadjuvant RT then surgery: curative RT sets the date
    tx_pathway == "Surgery + neoadjuvant RT"      ~ rt_date,
    
    # Neoadjuvant chemo then surgery: SACT starts the curative sequence
    tx_pathway == "Surgery + neoadjuvant chemo"   ~ sact_date,
    
    # All other surgical pathways: surgery is the curative treatment
    tx_pathway %in% c("Surgery + adjuvant chemo",
                      "Surgery only",
                      "Surgery + other")          ~ surgery_date,
    
    # Definitive chemoRT: whichever of SACT/RT starts first
    tx_pathway == "Definitive chemoRT"            ~ pmin(sact_date, rt_date,
                                                         na.rm = TRUE),
    
    # Curative RT only
    tx_pathway == "Curative RT only"              ~ rt_date,
    
    # Palliative and no treatment
    TRUE                                          ~ as.Date(NA)
  ),
  tx_trust = case_when(
    
    # Endoscopic resection and all surgical pathways: trust from HES
    tx_pathway %in% c("EMR/ESD only",
                      "EMR/ESD then surgery",
                      "Surgery + neoadjuvant chemo",
                      "Surgery + adjuvant chemo",
                      "Surgery only",
                      "Surgery + other")          ~ substr(PROCODE3, 1, 3),
    
    # Neoadjuvant chemoRT/RT then surgery: curative act is the RT, trust from RT
    tx_pathway %in% c("Surgery + neoadjuvant chemoRT",
                      "Surgery + neoadjuvant RT")  ~ substr(ORGCODEPROVIDER, 1, 3),
    
    # Non-surgical curative RT pathways: trust from RT
    tx_pathway %in% c("Definitive chemoRT",
                      "Curative RT only")          ~ substr(ORGCODEPROVIDER, 1, 3),
    
    # Palliative and no treatment
    TRUE                                           ~ NA_character_
  ),
  # --- Waiting time components ---------------------------------------------
  wt_dx_to_tx       = as.integer(first_tx_date - diagmdy),
  wt_endo_to_tx     = as.integer(first_tx_date - endoscopy_date),
  wt_dx_to_tx       = as.integer(first_tx_date - diagmdy),
  wt_endo_to_tx     = as.integer(first_tx_date - endoscopy_date),
  wt_dx_to_surg     = as.integer(surgery_date - diagmdy),
  wt_endo_to_surg   = as.integer(surgery_date - endoscopy_date),
  wt_dx_to_sact     = as.integer(sact_date - diagmdy),
  wt_endo_to_sact   = as.integer(sact_date - endoscopy_date),
  wt_sact_to_surg   = as.integer(surgery_date - sact_date),
  wt_surg_to_sact   = as.integer(sact_date - surgery_date),
  wt_dx_to_rt       = as.integer(rt_date - diagmdy),
  wt_endo_to_rt     = as.integer(rt_date - endoscopy_date),
  wt_rt_to_surg     = as.integer(surgery_date - rt_date),
  
  # --- Survival from surgery (PI7/PI8) -------------------------------------
  surv_from_surg_days = as.integer(finmdy - surgery_date),
  alive_90d_post_surg = had_surgery &
    !is.na(surv_from_surg_days) &
    (surv_from_surg_days > 90  | died == 0L),
  alive_1yr_post_surg = had_surgery &
    !is.na(surv_from_surg_days) &
    (surv_from_surg_days > 365 | died == 0L)
)  


saveRDS(og_cohort,
        "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/og_cohort_2015_2022.rds")

# =============================================================================
# COHORT SUMMARY
# =============================================================================
#


og_cohort %>%
  mutate(
    missing_reason = case_when(
      is.na(endoscopy_date) & is.na(first_tx_date)  ~ "no endoscopy, no treatment",
      is.na(endoscopy_date) & !is.na(first_tx_date) ~ "no endoscopy, has treatment",
      !is.na(endoscopy_date) & is.na(first_tx_date) ~ "has endoscopy, no treatment",
      wt_endo_to_tx < 0                             ~ "negative wait",
      TRUE                                           ~ "complete"
    )
  ) %>%
  count(missing_reason) %>%
  mutate(pct = round(100 * n / sum(n), 1))

og_cohort %>%
  summarise(
    pct_curative_surgery = round(100 * mean(had_curative_surgery), 1),
    pct_curative_tx      = round(100 * mean(received_curative_tx), 1)
  )

cat("Total OG cohort (stage 1-3):", nrow(og_cohort), "\n")

og_cohort %>%
  filter(tx_pathway == "SACT only") %>%
  count(had_rt, rt_curative)

og_cohort %>%
  count(tx_pathway, sort = TRUE) %>%
  mutate(pct = round(100 * n / sum(n), 1))

cat("\nTreatment presence:\n")
og_cohort %>%
  summarise(
    n                    = n(),
    pct_emresd           = round(100 * mean(had_emresd),           1),
    pct_surgery          = round(100 * mean(had_surgery),          1),
    pct_curative_surgery = round(100 * mean(had_curative_surgery), 1),
    pct_sact             = round(100 * mean(had_sact),             1),
    pct_rt               = round(100 * mean(had_rt),               1),
    pct_curative_rt      = round(100 * mean(had_curative_rt),      1),
    pct_curative_tx      = round(100 * mean(received_curative_tx), 1)
  ) %>%
  print()

cat("\nTreatment pathway:\n")
og_cohort %>%
  count(tx_pathway, sort = TRUE) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

cat("\nBy cancer subtype and pathway:\n")
og_cohort %>%
  count(cancer_subtype, tx_pathway) %>%
  group_by(cancer_subtype) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(cancer_subtype, desc(n)) %>%
  print(n = 40)

cat("\nWaiting time summary (surgical patients):\n")
og_cohort %>%
  filter(had_surgery) %>%
  summarise(
    n                  = n(),
    median_dx_surg     = median(wt_dx_to_surg,   na.rm = TRUE),
    p25_dx_surg        = quantile(wt_dx_to_surg,  0.25, na.rm = TRUE),
    p75_dx_surg        = quantile(wt_dx_to_surg,  0.75, na.rm = TRUE),
    median_endo_surg   = median(wt_endo_to_surg,  na.rm = TRUE),
    p25_endo_surg      = quantile(wt_endo_to_surg, 0.25, na.rm = TRUE),
    p75_endo_surg      = quantile(wt_endo_to_surg, 0.75, na.rm = TRUE),
    median_sact_surg   = median(wt_sact_to_surg[sact_before_surgery],
                                na.rm = TRUE),
    p25_sact_surg      = quantile(wt_sact_to_surg[sact_before_surgery],
                                  0.25, na.rm = TRUE),
    p75_sact_surg      = quantile(wt_sact_to_surg[sact_before_surgery],
                                  0.75, na.rm = TRUE)
  ) %>%
  print()

cat("\nWaiting times by pathway:\n")
og_cohort %>%
  filter(had_surgery) %>%
  group_by(tx_pathway) %>%
  summarise(
    n              = n(),
    median_dx_surg = median(wt_dx_to_surg, na.rm = TRUE),
    p25            = quantile(wt_dx_to_surg, 0.25, na.rm = TRUE),
    p75            = quantile(wt_dx_to_surg, 0.75, na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  arrange(desc(n)) %>%
  print()





##################################################################################



# Overall - compare with NOGCA Table 3 median 64 days (IQR 49-84) England
og_cohort %>%
  filter(!is.na(wt_endo_to_tx), wt_endo_to_tx >= 0) %>%
  summarise(
    n      = n(),
    median = median(wt_endo_to_tx),
    p25    = quantile(wt_endo_to_tx, 0.25),
    p75    = quantile(wt_endo_to_tx, 0.75)
  )

# By treatment intent (Table 3 columns: curative vs non-curative)
og_cohort %>%
  filter(!is.na(wt_endo_to_tx), wt_endo_to_tx >= 0) %>%
  mutate(intent = if_else(received_curative_tx, "curative", "non-curative")) %>%
  group_by(intent) %>%
  summarise(
    n      = n(),
    median = median(wt_endo_to_tx),
    p25    = quantile(wt_endo_to_tx, 0.25),
    p75    = quantile(wt_endo_to_tx, 0.75),
    .groups = "drop"
  )

# By primary treatment modality (Table 3 columns)
og_cohort %>%
  filter(!is.na(wt_endo_to_tx), wt_endo_to_tx >= 0) %>%
  group_by(tx_pathway) %>%
  summarise(
    n      = n(),
    median = median(wt_endo_to_tx),
    p25    = quantile(wt_endo_to_tx, 0.25),
    p75    = quantile(wt_endo_to_tx, 0.75),
    .groups = "drop"
  ) %>%
  arrange(desc(n))


og_cohort %>%
  mutate(
    missing_reason = case_when(
      is.na(endoscopy_date) & is.na(first_tx_date) ~ "no endoscopy, no treatment",
      is.na(endoscopy_date) & !is.na(first_tx_date) ~ "no endoscopy, has treatment",
      !is.na(endoscopy_date) & is.na(first_tx_date) ~ "has endoscopy, no treatment",
      wt_endo_to_tx < 0                             ~ "negative wait",
      TRUE                                           ~ "complete"
    )
  ) %>%
  count(missing_reason)


og_cohort %>%
  filter(ydiag %in% 2022:2023) %>%
  mutate(
    missing_reason = case_when(
      is.na(endoscopy_date) & is.na(first_tx_date) ~ "no endoscopy, no treatment",
      is.na(endoscopy_date) & !is.na(first_tx_date) ~ "no endoscopy, has treatment",
      !is.na(endoscopy_date) & is.na(first_tx_date) ~ "has endoscopy, no treatment",
      wt_endo_to_tx < 0                             ~ "negative wait",
      TRUE                                           ~ "complete"
    )
  ) %>%
  count(missing_reason) %>%
  mutate(pct = round(100 * n / sum(n), 1))


list.files(
  "E:/Data_PHE/Raw data files received from PHE READ ONLY/HES",
  full.names = FALSE
)

# What files are in OP folder?
list.files(
  "E:/Data_PHE/Raw data files received from PHE READ ONLY/HES/OP/",
  full.names = FALSE
) %>% head(20)


#op_test


op_test <- read_delim(
  "E:/Data_PHE/Raw data files received from PHE READ ONLY/HES/OP/FILE0224580_NIC656757_HES_OP_202299.txt",
  delim          = "|",
  n_max          = 1000,
  col_types      = cols(.default = col_character()),
  show_col_types = FALSE
)

names(op_test)
glimpse(op_test)

# Check operation fields
op_test %>%
  select(STUDY_ID, APPTDATE, starts_with("OPERTN_")) %>%
  head(10)

# Check how many OPERTN fields have real data vs NA
op_test %>%
  select(starts_with("OPERTN_")) %>%
  summarise(across(everything(), ~sum(!is.na(.) & . != ""))) %>%
  pivot_longer(everything(), names_to = "field", values_to = "n_populated") %>%
  filter(n_populated > 0)





