# =============================================================================
# 03  Extract the raw sources  (the optional, slow step)
# -----------------------------------------------------------------------------
# Reads the patient-restricted ultra-raw sources - HES APC, HES OP, SACT and
# RTDS - filtered to the cohort, and saves one extract per source under
# Data/ICON. These reads are the slow part of the whole build (each loops over
# many annual files and takes minutes), and they are pure data ingestion: they
# have nothing to do with the derivation logic that follows.
#
# Each source is gated independently: it is read from raw only when its extract
# is missing or refresh_raw = TRUE; otherwise this script does nothing for that
# source and the derivation scripts read the existing extract. So a normal
# rebuild skips this script's slow work entirely.
#
# Reads : HES APC parquet, HES OP txt, SACT csv, RTDS csv  (all gated)
#         Data/ICON/ncras_og_2015_2022.rds  (for the cohort id filter)
# Writes: hes_apc / hes_op / sact / rtds extracts under Data/ICON
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

ncras_og_ids <- readRDS(f_ncras_cohort) %>% distinct(pseudo_patientid) %>% pull()

# explain why a gated source is being read: either the extract is missing (first
# run) or a refresh was forced. Makes the slow reads self-explanatory on screen.
why_reading <- function(label, extract) {
  reason <- if (!file.exists(extract)) "extract not found - first run"
  else "refresh_raw = TRUE - rebuilding"
  message(sprintf("Reading %s from raw (%s); this is slow...", label, reason))
}

# -----------------------------------------------------------------------------
# HES APC  (operations, dates and secondary diagnoses; 2014-2024)
# -----------------------------------------------------------------------------
if (refresh_raw || !file.exists(f_hes_apc_extract)) {
  why_reading("HES APC", f_hes_apc_extract)
  hes_apc_files <- list.files(path_hes_apc_dir, pattern = "FILE*", full.names = TRUE) %>%
    keep(~{
      yr <- str_extract(.x, "(?<=HES_APC_)\\d{4}") %>% as.integer()
      !is.na(yr) && yr %in% 2014:2024
    })
  stopifnot(length(hes_apc_files) > 0)
  hes_cols_select <- c("STUDY_ID","ADMIDATE","ADMIMETH","PROCODE3","SITETRET",
                       "EPISTART","EPIORDER","EPITYPE", op_cols, opdate_cols, diag_cols)
  hes_apc <- map_dfr(hes_apc_files, ~{
    read_parquet(.x, col_select = all_of(hes_cols_select)) %>%
      filter(STUDY_ID %in% ncras_og_ids) %>%
      mutate(STUDY_ID = as.character(STUDY_ID), ADMIMETH = as.character(ADMIMETH),
             EPISTART = as.Date(EPISTART), ADMIDATE = as.Date(ADMIDATE),
             across(all_of(op_cols), as.character),
             across(all_of(opdate_cols), as.Date),
             across(any_of(diag_cols), as.character))
  }, .progress = TRUE)
  saveRDS(hes_apc, f_hes_apc_extract)
  cat("HES APC extract:", nrow(hes_apc), "rows,",
      n_distinct(hes_apc$STUDY_ID), "patients\n")
} else {
  cat("HES APC extract present - skipping read.\n")
}

# -----------------------------------------------------------------------------
# HES OP  (outpatient endoscopy; attended appointments only; 2014+)
# -----------------------------------------------------------------------------
if (refresh_raw || !file.exists(f_hes_op_extract)) {
  why_reading("HES OP", f_hes_op_extract)
  hes_op_files <- list.files(path_hes_op_dir, pattern = "*.txt", full.names = TRUE) %>%
    keep(~{
      yr <- str_extract(.x, "(?<=HES_OP_)\\d{4}") %>% as.integer()
      !is.na(yr) && yr >= 2014
    })
  op_cols_select <- c("STUDY_ID","APPTDATE","ATTENDED",
                      paste0("OPERTN_0", 1:9), paste0("OPERTN_", 10:24),
                      "PROCODET","TRETSPEF","MAINSPEF")
  hes_op <- map_dfr(hes_op_files, ~{
    read_delim(.x, delim = "|", col_select = any_of(op_cols_select),
               col_types = cols(.default = col_character()), show_col_types = FALSE) %>%
      filter(STUDY_ID %in% ncras_og_ids, ATTENDED %in% c("5", "6")) %>%
      mutate(STUDY_ID = as.character(STUDY_ID), appt_date = as.Date(APPTDATE))
  }, .progress = TRUE)
  saveRDS(hes_op, f_hes_op_extract)
  cat("HES OP extract:", nrow(hes_op), "rows,",
      n_distinct(hes_op$STUDY_ID), "patients\n")
} else {
  cat("HES OP extract present - skipping read.\n")
}

# -----------------------------------------------------------------------------
# SACT  (systemic anti-cancer therapy)
# -----------------------------------------------------------------------------
if (refresh_raw || !file.exists(f_sact_extract)) {
  why_reading("SACT", f_sact_extract)
  sact_cols <- c("PSEUDO_PATIENTID","PRIMARY_DIAGNOSIS","MORPHOLOGY_CLEAN",
                 "BENCHMARK_GROUP","ANALYSIS_GROUP","INTENT_OF_TREATMENT_V3",
                 "START_DATE_OF_REGIMEN","START_DATE_OF_CYCLE",
                 "DATE_DECISION_TO_TREAT","DATE_OF_FINAL_TREATMENT",
                 "CYCLE_NUMBER","NUMBER_OF_CYCLES_PLANNED",
                 "ORGANISATION_CODE_OF_PROVIDER","PERF_STAT_START_OF_REG_ADULT",
                 "STAGE_AT_START","REGIMEN_MOD_STOPPED_EARLY",
                 "REGIMEN_OUTCOME_SUMMARY","CHEMO_RADIATION")
  sact_files <- list.files(path_sact_dir, pattern = "*.csv", full.names = TRUE)
  sact_og <- map_dfr(sact_files, ~{
    read_csv(.x, col_select = all_of(sact_cols),
             col_types = cols(.default = col_character()), show_col_types = FALSE) %>%
      filter(str_starts(PRIMARY_DIAGNOSIS, "C15") | str_starts(PRIMARY_DIAGNOSIS, "C16")) %>%
      mutate(pseudo_patientid    = as.character(as.integer(PSEUDO_PATIENTID)),
             sact_regimen_date   = as.Date(START_DATE_OF_REGIMEN,   "%d/%m/%Y"),
             sact_cycle_date     = as.Date(START_DATE_OF_CYCLE,     "%d/%m/%Y"),
             date_decision_treat = as.Date(DATE_DECISION_TO_TREAT,  "%d/%m/%Y"),
             date_final_tx       = as.Date(DATE_OF_FINAL_TREATMENT, "%d/%m/%Y"),
             cycle_number        = as.integer(CYCLE_NUMBER),
             benchmark_group_lwr = tolower(trimws(BENCHMARK_GROUP))) %>%
      filter(pseudo_patientid %in% ncras_og_ids)
  }, .progress = TRUE)
  saveRDS(sact_og, f_sact_extract)
  cat("SACT extract:", nrow(sact_og), "rows,",
      n_distinct(sact_og$pseudo_patientid), "patients\n")
} else {
  cat("SACT extract present - skipping read.\n")
}

# -----------------------------------------------------------------------------
# RTDS  (radiotherapy)
# -----------------------------------------------------------------------------
if (refresh_raw || !file.exists(f_rtds_extract)) {
  why_reading("RTDS", f_rtds_extract)
  rtds_cols <- c("PSEUDO_PATIENTID","RADIOTHERAPYDIAGNOSISICD","TREATMENTSTARTDATE",
                 "DECISIONTOTREATDATE","RADIOTHERAPYINTENT","RTPRESCRIBEDDOSE",
                 "PRESCRIBEDFRACTIONS","RTTREATMENTREGION",
                 "RTTREATMENTANATOMICALSITE","ORGCODEPROVIDER")
  rtds_files <- list.files(path_rtds_dir, pattern = "*.csv", full.names = TRUE)
  rtds_og <- map_dfr(rtds_files, ~{
    read_csv(.x, col_select = all_of(rtds_cols),
             col_types = cols(.default = col_character()), show_col_types = FALSE) %>%
      filter(str_starts(RADIOTHERAPYDIAGNOSISICD, "C15") |
               str_starts(RADIOTHERAPYDIAGNOSISICD, "C16")) %>%
      mutate(pseudo_patientid = as.character(as.integer(as.numeric(PSEUDO_PATIENTID))),
             rt_start_date    = as.Date(TREATMENTSTARTDATE,  "%d/%m/%Y"),
             rt_decision_date = as.Date(DECISIONTOTREATDATE, "%d/%m/%Y"),
             rt_dose          = as.numeric(RTPRESCRIBEDDOSE),
             rt_fractions     = as.integer(PRESCRIBEDFRACTIONS)) %>%
      filter(pseudo_patientid %in% ncras_og_ids)
  }, .progress = TRUE)
  saveRDS(rtds_og, f_rtds_extract)
  cat("RTDS extract:", nrow(rtds_og), "rows,",
      n_distinct(rtds_og$pseudo_patientid), "patients\n")
} else {
  cat("RTDS extract present - skipping read.\n")
}

cat("03 complete. Next: 04_derive_comorbidities.R\n")