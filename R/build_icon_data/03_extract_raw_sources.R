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
  hes_cols_select <- c("STUDY_ID","ADMIDATE","ADMIMETH","PROCODE3","SITETRET",
                       "EPISTART","EPIORDER","EPITYPE", op_cols, opdate_cols, diag_cols)
  # partitioned Parquet store: the year filter prunes partitions, the id filter
  # and column projection are pushed into the scan, so only the needed slice is read
  hes_apc <- open_dataset(path_hes_apc_dir) %>%
    filter(year >= 2014, year <= 2024, STUDY_ID %in% ncras_og_ids) %>%
    select(any_of(hes_cols_select)) %>%
    collect() %>%
    mutate(STUDY_ID = as.character(STUDY_ID), ADMIMETH = as.character(ADMIMETH),
           EPISTART = as.Date(EPISTART), ADMIDATE = as.Date(ADMIDATE),
           across(all_of(op_cols), as.character),
           across(all_of(opdate_cols), as.Date),
           across(any_of(diag_cols), as.character))
  if (nrow(hes_apc) == 0)
    stop("HES APC read returned 0 rows - likely an id-format mismatch between ",
         "STUDY_ID and the NCRAS ids. Not saving an empty extract.", call. = FALSE)
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
  op_cols_select <- c("STUDY_ID","APPTDATE","appt_date","ATTENDED",
                      "PROCODET","TRETSPEF","MAINSPEF",
                      paste0("OPERTN_", str_pad(1:24, 2, pad = "0")))
  # same partitioned-Parquet read; ATTENDED 5/6 (attended) pushed into the scan
  hes_op <- open_dataset(path_hes_op_dir) %>%
    filter(year >= 2014, STUDY_ID %in% ncras_og_ids, ATTENDED %in% c("5", "6")) %>%
    select(any_of(op_cols_select)) %>%
    collect() %>%
    mutate(STUDY_ID = as.character(STUDY_ID),
           appt_date = if ("appt_date" %in% names(.)) as.Date(appt_date)
           else as.Date(APPTDATE))
  if (nrow(hes_op) == 0)
    warning("HES OP read returned 0 rows. OP endoscopy can be sparse, but check ",
            "this is not an id-format mismatch (STUDY_ID vs NCRAS ids).", call. = FALSE)
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
  if (nrow(sact_og) == 0)
    stop("SACT read returned 0 rows. Likely an id-format mismatch between SACT's ",
         "PSEUDO_PATIENTID and the NCRAS ids, or no C15/C16 rows in the source. ",
         "Not saving an empty extract.", call. = FALSE)
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
  if (nrow(rtds_og) == 0)
    stop("RTDS read returned 0 rows. Likely an id-format mismatch or no C15/C16 ",
         "rows in the source. Not saving an empty extract.", call. = FALSE)
  saveRDS(rtds_og, f_rtds_extract)
  cat("RTDS extract:", nrow(rtds_og), "rows,",
      n_distinct(rtds_og$pseudo_patientid), "patients\n")
} else {
  cat("RTDS extract present - skipping read.\n")
}

cat("03 complete. Next: 04_derive_hes_treatments.R\n")