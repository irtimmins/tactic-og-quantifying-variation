# =============================================================================
# 01  Define parameters
# -----------------------------------------------------------------------------
# Defines everything the build needs before any data is touched: the file paths,
# the refresh-raw switch, the interval and window constants (each explained), the
# clinical code lists (ICD-10, morphology, OPCS), and the shared helpers. Sourced
# by every later script, and by 00_master; it reads and writes nothing itself.
#
# Keeping all of these in one place means the analytic decisions - which windows,
# which codes - are visible and auditable up front, rather than scattered through
# the build.
# =============================================================================

suppressPackageStartupMessages({
  library(arrow)
  library(haven)
  library(tidyverse)
  library(lubridate)
})

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
# dir_raw  : the ultra-raw source extracts (parquet / txt / dta) supplied by PHE.
#            Only script 02 reads these, and only when a refresh is requested.
# dir_icon : the project data folder. All derived extracts, anchors and the final
#            cohort are written here, and every downstream script reads from here.
# Both can be pre-set before sourcing this file (a test harness or 00_master can
# point them elsewhere); only the defaults are filled in when they are unset.
if (!exists("dir_raw"))  dir_raw  <- "E:/Data_PHE"
if (!exists("dir_icon")) dir_icon <- "Data/ICON"
dir.create(dir_icon, recursive = TRUE, showWarnings = FALSE)

# specific raw source locations (only used by 02_extract_raw_sources.R)
path_ncras_parquet <- file.path(dir_raw,
                                "Extracts/#2045_ICON_TACTIC/NCRAS/NCRAS_clean_1995_2022_route_sitestr")
path_cosd_dta      <- file.path(dir_raw,
                                "Raw data files received from PHE READ ONLY/NCRAS/Stata files/18_COSD_data.dta")
path_hes_apc_dir   <- file.path(dir_raw, "Extracts/#2045_ICON_TACTIC/HES/APC")
path_hes_op_dir    <- file.path(dir_raw,
                                "Raw data files received from PHE READ ONLY/HES/OP")
path_sact_dir      <- file.path(dir_raw, "Extracts/#2045_ICON_TACTIC/SACT")
path_rtds_dir      <- file.path(dir_raw, "Extracts/#2045_ICON_TACTIC/RTDS")
if (!exists("path_cwt_partition"))
  path_cwt_partition <- file.path(dir_raw,
                                  "Extracts/#2045_ICON_TACTIC/CWT/11_CWT_data_partitioned")

# derived extract / anchor / cohort files (all under dir_icon)
f_hes_apc_extract  <- file.path(dir_icon, "hes_apc_og_2014_2022.rds")
f_hes_op_extract   <- file.path(dir_icon, "hes_op_og_2014_2024.rds")
f_sact_extract     <- file.path(dir_icon, "sact_og_2012_2024.rds")
f_rtds_extract     <- file.path(dir_icon, "rtds_og_2009_2024.rds")
f_ncras_cohort     <- file.path(dir_icon, "ncras_og_2015_2022.rds")
f_cci              <- file.path(dir_icon, "og_cci_2015_2022.rds")
f_endoscopy_anchor <- file.path(dir_icon, "og_endoscopy_anchor.rds")
f_emresd_anchor    <- file.path(dir_icon, "og_emresd_anchor.rds")
f_surgery_anchor   <- file.path(dir_icon, "og_surgery_anchor.rds")
f_chemo_anchor     <- file.path(dir_icon, "og_chemo_anchor_2015_2022.rds")
f_rt_anchor        <- file.path(dir_icon, "og_rt_anchor.rds")
f_cohort           <- file.path(dir_icon, "og_cohort_2015_2022.rds")
f_cohort_cwt       <- file.path(dir_icon, "og_cohort_cwt_2015_2022.rds")

# -----------------------------------------------------------------------------
# Refresh switch
# -----------------------------------------------------------------------------
# The ultra-raw reads in script 02 (NCRAS, HES APC, HES OP, SACT, RTDS) take
# several minutes each. They only need to run when the source data changes. With
# refresh_raw = FALSE (the default) script 02 reuses the saved extracts under
# dir_icon and skips the slow reads; set TRUE to re-pull from the raw sources.
# 00_master.R can override this for a whole run.
if (!exists("refresh_raw")) refresh_raw <- FALSE

# -----------------------------------------------------------------------------
# Interval and window constants  (the analytic time decisions, all in days)
# -----------------------------------------------------------------------------
# tx_window_days   : the treatment-capture window. A treatment counts towards the
#                    cohort if it falls within 9 months (~275 days) of diagnosis,
#                    the NOGCA standard for treatment ascertainment.
tx_window_days   <- 275L

# endoscopy_pre_dx_days : a diagnostic endoscopy is taken as the PI3 clock start
#                    if it falls on or up to 30 days before the diagnosis date.
endoscopy_pre_dx_days <- 30L

# dtt_min_offset   : in the CWT merge, the decision-to-treat may sit at most this
#                    many days before diagnosis and still be valid (allows for the
#                    odd pre-diagnosis MDT decision).
dtt_min_offset   <- -30L

# cwt_window_days  : the CWT decision-to-treat must sit within this many days of
#                    diagnosis to be an eligible anchor (matches tx ascertainment;
#                    270 in the merge, slightly tighter than the 275 capture).
cwt_window_days  <- 270L

# treat_tol_days   : a CWT treatment date may precede the decision-to-treat by up
#                    to this many days and the record still count as valid; beyond
#                    it the timing is treated as inconsistent.
treat_tol_days   <- 14L

# chemo_rt_concurrent_days : SACT and curative RT within this many days are taken
#                    as a concurrent chemoradiotherapy course.
chemo_rt_concurrent_days <- 14L

# hes_chemo_near_rt_days : HES-only chemotherapy (no SACT record) counts towards a
#                    definitive-chemoRT classification only when it sits within
#                    this many days of the curative RT, so a separate HES chemo
#                    episode does not manufacture chemoRT.
hes_chemo_near_rt_days <- 28L

# surg_switch_date / surg_01_rule : CWT surgery code 01 was retired on
#                    2020-10-01 and replaced by 23/24. date_split counts 01 as
#                    surgery only before the switch and 23/24 only on/after it.
surg_switch_date <- as.Date("2020-10-01")
surg_01_rule     <- "date_split"   # one of: date_split | always | never

# treatment pathway levels, in reporting order
tx_pathway_levels <- c(
  "EMR/ESD only", "EMR/ESD then surgery",
  "Surgery + neoadjuvant chemoRT", "Surgery + neoadjuvant chemo",
  "Surgery + neoadjuvant RT", "Surgery + adjuvant chemo",
  "Surgery only", "Surgery + other",
  "Definitive chemoRT", "Curative RT only",
  "Palliative chemo + RT", "SACT only", "Palliative RT only",
  "No treatment recorded")

# =============================================================================
# Clinical code lists
# =============================================================================

# --- NCRAS columns to keep ---------------------------------------------------
ncras_cols <- c(
  "pseudo_patientid", "pseudo_tumourid", "diagmdy", "ydiag",
  "cancer", "sitestr", "typestr", "basisofdiagnosis", "grade", "behav",
  "morphology",
  "stage_best", "stage_best_system", "t_best", "n_best", "m_best",
  "t_path", "n_path", "m_path",
  "sex", "agediag", "birthmdy", "ethnicity_group_broad",
  "lsoa11_code", "NHSE_reversed_imd_quintile_lsoas",
  "canalliance_2024_code", "canalliance_2024_name",
  "diag_trust", "diag_trust_name",
  "first_trust", "first_trust_name", "first_hosp_date", "diag_hosp",
  "route_bjc", "final_route", "route_code",
  "tumour_performancestatus", "clinicalnursespecialist", "firstmdtmeetingdate",
  "sg_flag", "rt_flag", "ct_flag",
  "dead", "finmdy", "dco")

# --- ICD-10 site (oesophageal C15x, gastric C16x) ----------------------------
og_icd10 <- c("C15","C150","C151","C152","C153","C154","C155","C158","C159",
              "C16","C160","C161","C162","C163","C164","C165","C166","C168","C169")

# --- Morphology: epithelial inclusion (Appendix 4), neuroendocrine exclusion --
morph_epithelial <- c(
  8005, 8010, 8020, 8021, 8032, 8033, 8050, 8051, 8052,
  8070, 8071, 8072, 8073, 8074, 8075, 8076, 8077, 8078, 8083, 8084,
  8140, 8141, 8142, 8143, 8144, 8145,
  8190, 8210, 8211, 8213, 8214, 8231,
  8255, 8260, 8261, 8262, 8263, 8310, 8323, 8430, 8440,
  8480, 8481, 8490, 8510, 8512, 8560, 8562,
  8570, 8571, 8572, 8573, 8574, 8576, 8982)
morph_neuroendocrine <- c(
  8013, 8041, 8042, 8043, 8044, 8045,
  8150, 8151, 8152, 8153, 8154, 8155, 8156, 8157, 8158,
  8240, 8241, 8242, 8243, 8244, 8245, 8246, 8247, 8249, 9091)

# --- Subtype morphology (oesophageal SCC vs ACA; Appendix 3) ------------------
morph_oes_scc <- c(8033, 8051, 8052, 8070, 8071, 8072, 8073, 8074, 8075,
                   8076, 8077, 8078, 8083, 8084)
morph_oes_aca <- c(8005, 8140, 8141, 8142, 8143, 8144, 8145,
                   8190, 8210, 8211, 8213, 8214, 8255, 8260, 8261, 8262, 8263,
                   8310, 8323, 8440, 8480, 8481,
                   8570, 8571, 8572, 8573, 8574, 8576)

# --- OPCS-4: diagnostic endoscopy (Appendix 6) -------------------------------
opcs_diagnostic_endoscopy <- c(
  "G142","G143","G145","G147",
  "G152","G153","G154","G156","G157","G158","G159",
  "G161","G162","G168","G169",
  "G172","G173","G188","G189",
  "G191","G198","G199",
  "G201","G202","G208","G209",
  "G214","G215","G218","G219",
  "G422","G432","G433","G435",
  "G441","G443","G445","G446","G448","G449",
  "G451","G452","G454","G458","G459",
  "G462","G463","G468","G469")

# --- OPCS-4: EMR / ESD endotherapy (Appendix 7) ------------------------------
opcs_emresd <- c(
  "G121","G128","G129", "G141","G146","G148","G149", "G171","G178","G179",
  "G421","G423","G428","G429", "G431","G438","G439",
  "G143","G145","G433","G435")   # ablation codes for HGD

# --- OPCS-4: major OG resection (Appendix 8) ---------------------------------
opcs_og_surgery <- list(
  oesophagectomy = c("G011","G018","G019",
                     "G021","G022","G023","G024","G025","G028","G029",
                     "G031","G032","G033","G034","G035","G036","G038","G039"),
  oesophagogastrectomy_jejunum = c("G012","G013"),
  total_gastrectomy   = c("G271","G272","G273","G274","G275","G278","G279"),
  partial_gastrectomy = c("G281","G282","G283","G288","G289"))
opcs_og_surgery_all    <- unique(unlist(opcs_og_surgery))
opcs_oesophagectomy    <- unique(c(opcs_og_surgery$oesophagectomy,
                                   opcs_og_surgery$oesophagogastrectomy_jejunum))
opcs_gastrectomy_total <- opcs_og_surgery$total_gastrectomy
opcs_gastrectomy_part  <- opcs_og_surgery$partial_gastrectomy
opcs_og_surgery_lookup <- c(
  setNames(rep("oesophagectomy",     length(opcs_oesophagectomy)),    opcs_oesophagectomy),
  setNames(rep("total_gastrectomy",  length(opcs_gastrectomy_total)), opcs_gastrectomy_total),
  setNames(rep("partial_gastrectomy",length(opcs_gastrectomy_part)),  opcs_gastrectomy_part))

# --- OPCS-4: HES chemotherapy delivery (supplements SACT) --------------------
opcs_chemo_delivery <- c(
  "X701","X702","X703","X704","X709", "X711","X712","X713","X714","X719",
  "X721","X722","X723","X724","X725","X726","X727","X728","X729",
  "X731","X732","X738","X739", "X741","X742","X743","X744","X748","X749")
icd_chemo_attendance <- c("Z511")   # encounter for antineoplastic chemotherapy

# --- HES emergency admission methods -----------------------------------------
admimeth_emerg <- c("21","22","23","24","25","28","2A","2B","2C","2D")

# =============================================================================
# Shared helpers
# =============================================================================

# normalise an OPCS code for matching (upper case, dots removed)
normalise_opcs <- function(x) str_replace_all(str_to_upper(as.character(x)), "\\.", "")

# guard a saved extract against staleness: the file-existence gate in 03 reuses
# whatever is on disk, so an extract written before a change to 03 may lack a
# derived column. This stops with a clear, actionable message naming the missing
# columns rather than letting a downstream step fail with "object not found".
check_extract <- function(df, needed, label, extract_path) {
  miss <- setdiff(needed, names(df))
  if (length(miss))
    stop(label, " extract is missing columns: ", paste(miss, collapse = ", "),
         "\n  - it is stale. Delete ", basename(extract_path),
         " and re-run 03_extract_raw_sources.R (or set refresh_raw <- TRUE).",
         call. = FALSE)
}

# pivot HES OPERTN_/OPDATE_ fields long and keep episodes whose OPCS is in a list.
# Used by the endoscopy, EMR/ESD and surgery anchors so the pivot lives once.
match_opcs_episodes <- function(hes_data, opcs_list, op_cols, opdate_cols) {
  ops_long <- hes_data %>%
    filter(!is.na(OPERTN_01), OPERTN_01 != "-") %>%
    pivot_longer(all_of(op_cols), names_to = "op_position", values_to = "opcs_code") %>%
    filter(!is.na(opcs_code), opcs_code != "-") %>%
    mutate(opcs4 = normalise_opcs(opcs_code),
           op_position_n = as.integer(str_extract(op_position, "[0-9]+"))) %>%
    filter(opcs4 %in% opcs_list)
  if (nrow(ops_long) == 0) return(tibble())
  dates_long <- hes_data %>%
    pivot_longer(all_of(opdate_cols), names_to = "opdate_position", values_to = "op_date") %>%
    mutate(op_position_n = as.integer(str_extract(opdate_position, "[0-9]+")),
           op_date = as.Date(op_date)) %>%
    select(STUDY_ID, EPISTART, EPIORDER, op_position_n, op_date)
  ops_long %>%
    left_join(dates_long,
              by = c("STUDY_ID","EPISTART","EPIORDER","op_position_n"),
              relationship = "many-to-many") %>%
    rename(pseudo_patientid = STUDY_ID)
}

# HES wide-field name vectors (24 operation slots, 20 diagnosis slots)
op_cols     <- paste0("OPERTN_", str_pad(1:24, 2, pad = "0"))
opdate_cols <- paste0("OPDATE_", str_pad(1:24, 2, pad = "0"))
diag_cols   <- paste0("DIAG_4_", str_pad(1:20, 2, pad = "0"))

cat("01 parameters defined: paths under", dir_icon,
    "| refresh_raw =", refresh_raw, "| tx window", tx_window_days, "days\n")