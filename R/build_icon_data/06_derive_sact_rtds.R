# =============================================================================
# 06  Derive SACT and RTDS
# -----------------------------------------------------------------------------
# Builds the systemic-therapy and radiotherapy anchors from the SACT and RTDS
# extracts. SACT and RTDS are kept together here because they are the two
# non-surgical treatment modalities and the pathway logic (script 07) reasons
# about them jointly (concurrent chemoradiotherapy, definitive chemoRT).
#
#   chemo anchor : first eligible SACT regimen in the window (GIST / NE regimens
#                  excluded), supplemented by HES chemotherapy delivery (OPCS
#                  X70-74, ICD-10 Z51.1) so chemo coded only in HES is not missed.
#                  Provenance (sact / hes / both) and the HES chemo date are kept.
#   RT anchor    : the RT course in the window, classified curative vs palliative
#                  by dose and fractionation, preferring the curative course.
#
# Reads : ncras cohort, sact / rtds / hes_apc extracts
# Writes: og_chemo_anchor_2015_2022.rds, og_rt_anchor.rds
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

ncras_og     <- readRDS(f_ncras_cohort)
ncras_og_ids <- ncras_og %>% distinct(pseudo_patientid) %>% pull()
sact_og      <- readRDS(f_sact_extract)
rtds_og      <- readRDS(f_rtds_extract)
hes_apc      <- readRDS(f_hes_apc_extract)

# the extracts must carry the derived columns this script needs (check_extract is
# defined in 01); a stale extract is reported clearly rather than failing later.
check_extract(sact_og, c("pseudo_patientid","sact_regimen_date","benchmark_group_lwr",
                         "BENCHMARK_GROUP","INTENT_OF_TREATMENT_V3",
                         "ORGANISATION_CODE_OF_PROVIDER","CHEMO_RADIATION"),
              "SACT", f_sact_extract)
check_extract(rtds_og, c("pseudo_patientid","rt_start_date","rt_dose",
                         "rt_fractions","ORGCODEPROVIDER"),
              "RTDS", f_rtds_extract)

# =============================================================================
# Chemotherapy: SACT anchor, HES supplement, combined
# =============================================================================

# regimens excluded as not OG chemotherapy (GIST agents; neuroendocrine doublets)
gist_regimens <- c("imatinib", "sunitinib", "regorafenib")
ne_regimens   <- c("carboplatin + etoposide", "cisplatin + etoposide",
                   "carboplatin + sunitinib", "cisplatin + sunitinib",
                   "carboplatin + everolimus", "cisplatin + everolimus")

sact_anchor <- ncras_og %>%
  select(pseudo_patientid, diagmdy) %>%
  left_join(
    sact_og %>%
      filter(!benchmark_group_lwr %in% gist_regimens,
             !benchmark_group_lwr %in% ne_regimens,
             !BENCHMARK_GROUP %in% c("NOT CHEMO", "TRIAL UNSPECIFIED")) %>%
      select(pseudo_patientid, sact_regimen_date, sact_cycle_date,
             BENCHMARK_GROUP, benchmark_group_lwr, INTENT_OF_TREATMENT_V3,
             CYCLE_NUMBER, cycle_number, ORGANISATION_CODE_OF_PROVIDER, CHEMO_RADIATION),
    by = "pseudo_patientid") %>%
  mutate(days_dx_to_sact = as.integer(sact_regimen_date - diagmdy)) %>%
  filter(!is.na(days_dx_to_sact),
         days_dx_to_sact >= -30, days_dx_to_sact <= tx_window_days) %>%
  arrange(pseudo_patientid, sact_regimen_date) %>%
  distinct(pseudo_patientid, .keep_all = TRUE) %>%
  rename(sact_date = sact_regimen_date) %>%
  select(pseudo_patientid, sact_date, days_dx_to_sact, BENCHMARK_GROUP,
         benchmark_group_lwr, INTENT_OF_TREATMENT_V3, CHEMO_RADIATION,
         ORGANISATION_CODE_OF_PROVIDER)

# HES chemotherapy supplement: delivery via OPCS X70-74 or ICD-10 Z51.1.
# Carries no benchmark/intent (provenance is flagged "hes" so regimen-specific
# analysis can exclude it).
hes_sub <- hes_apc %>%
  mutate(STUDY_ID = as.character(STUDY_ID)) %>%
  filter(STUDY_ID %in% ncras_og_ids)

hes_opcs <- hes_sub %>%
  select(STUDY_ID, EPISTART, starts_with("OPERTN_")) %>%
  pivot_longer(starts_with("OPERTN_"), names_to = "pos", values_to = "opcs",
               names_prefix = "OPERTN_") %>%
  left_join(hes_sub %>% select(STUDY_ID, EPISTART, starts_with("OPDATE_")) %>%
              pivot_longer(starts_with("OPDATE_"), names_to = "pos",
                           values_to = "opdate", names_prefix = "OPDATE_"),
            by = c("STUDY_ID", "EPISTART", "pos"), relationship = "many-to-many") %>%
  mutate(opcs4 = str_to_upper(str_remove_all(str_trim(opcs), "\\."))) %>%
  filter(opcs4 %in% opcs_chemo_delivery) %>%
  transmute(STUDY_ID, chemo_date = coalesce(as.Date(opdate), as.Date(EPISTART)))

hes_icd <- hes_sub %>%
  select(STUDY_ID, EPISTART, starts_with("DIAG_4_")) %>%
  pivot_longer(starts_with("DIAG_4_"), names_to = "pos", values_to = "icd",
               names_prefix = "DIAG_4_") %>%
  mutate(icd4 = str_to_upper(str_sub(str_remove_all(str_trim(icd), "\\."), 1, 4))) %>%
  filter(icd4 %in% icd_chemo_attendance) %>%
  transmute(STUDY_ID, chemo_date = as.Date(EPISTART))

hes_chemo_anchor <- bind_rows(hes_opcs, hes_icd) %>%
  filter(!is.na(chemo_date)) %>%
  rename(pseudo_patientid = STUDY_ID) %>%
  inner_join(ncras_og %>% select(pseudo_patientid, diagmdy), by = "pseudo_patientid") %>%
  mutate(days_dx_to_hes_chemo = as.integer(chemo_date - diagmdy)) %>%
  filter(days_dx_to_hes_chemo >= -30, days_dx_to_hes_chemo <= tx_window_days) %>%
  arrange(pseudo_patientid, chemo_date) %>%
  distinct(pseudo_patientid, .keep_all = TRUE) %>%
  transmute(pseudo_patientid, hes_chemo_date = chemo_date, days_dx_to_hes_chemo)

# combine: SACT-preferred date, provenance and HES date kept for the guard in 07
chemo_anchor <- ncras_og %>%
  select(pseudo_patientid, diagmdy) %>%
  left_join(sact_anchor, by = "pseudo_patientid") %>%
  left_join(hes_chemo_anchor, by = "pseudo_patientid") %>%
  filter(!is.na(sact_date) | !is.na(hes_chemo_date)) %>%
  mutate(
    chemo_source = case_when(
      !is.na(sact_date) &  !is.na(hes_chemo_date) ~ "both",
      !is.na(sact_date) &   is.na(hes_chemo_date) ~ "sact",
      TRUE                                        ~ "hes"),
    chemo_date       = coalesce(sact_date, hes_chemo_date),
    days_dx_to_chemo = as.integer(chemo_date - diagmdy))

saveRDS(chemo_anchor, f_chemo_anchor)
cat("Chemo anchor:", nrow(chemo_anchor), "patients (",
    paste(names(table(chemo_anchor$chemo_source)),
          table(chemo_anchor$chemo_source), collapse = " "), ")\n")

# =============================================================================
# Radiotherapy: RT anchor, curative vs palliative by dose / fractionation
# =============================================================================
# Curative-intent dose-fractionation schedules (NOGCA radiotherapy appendix):
# everything else is treated as palliative.
rt_anchor <- ncras_og %>%
  select(pseudo_patientid, diagmdy) %>%
  left_join(
    rtds_og %>%
      arrange(pseudo_patientid, rt_start_date) %>%
      distinct(pseudo_patientid, rt_start_date, .keep_all = TRUE) %>%
      mutate(rt_curative = case_when(
        rt_dose == 50    & rt_fractions == 25            ~ TRUE,
        rt_dose == 50.4  & rt_fractions == 28            ~ TRUE,
        rt_dose == 60    & rt_fractions == 30            ~ TRUE,
        rt_dose == 50    & rt_fractions %in% c(15, 16)   ~ TRUE,
        between(rt_dose, 50, 55)   & rt_fractions == 20  ~ TRUE,
        between(rt_dose, 45, 52.5) & rt_fractions %in% c(15, 16) ~ TRUE,
        rt_dose == 41.4  & rt_fractions == 23            ~ TRUE,
        rt_dose == 45    & rt_fractions == 25            ~ TRUE,
        TRUE                                             ~ FALSE)) %>%
      select(pseudo_patientid, rt_start_date, rt_curative,
             rt_dose, rt_fractions, ORGCODEPROVIDER),
    by = "pseudo_patientid") %>%
  mutate(days_dx_to_rt = as.integer(rt_start_date - diagmdy)) %>%
  filter(!is.na(days_dx_to_rt),
         days_dx_to_rt >= -30, days_dx_to_rt <= tx_window_days) %>%
  # prefer the curative course, then the earliest within each group
  arrange(pseudo_patientid, desc(rt_curative), rt_start_date) %>%
  distinct(pseudo_patientid, .keep_all = TRUE) %>%
  rename(rt_date = rt_start_date) %>%
  select(pseudo_patientid, rt_date, rt_curative,
         rt_dose, rt_fractions, days_dx_to_rt, ORGCODEPROVIDER)

saveRDS(rt_anchor, f_rt_anchor)
cat("RT anchor:", n_distinct(rt_anchor$pseudo_patientid), "patients (curative",
    sum(rt_anchor$rt_curative, na.rm = TRUE), ").",
    "Next: 07_build_pathways.R\n")

# ---- optional checks (uncomment to inspect) ---------------------------------
# count(rt_anchor, rt_curative) %>% print()
# rt_anchor %>% filter(rt_curative) %>% count(rt_dose, rt_fractions, sort = TRUE) %>% print()