# =============================================================================
# 04  Derive HES treatments
# -----------------------------------------------------------------------------
# From the HES APC and HES OP extracts, build the three HES-coded treatment
# anchors, one row per patient:
#   - diagnostic endoscopy : first qualifying endoscopy on/up to 30 days before
#                            diagnosis (the PI3 clock start), from APC then OP
#   - EMR/ESD              : first endotherapy in the treatment window
#   - major surgery        : first elective OG resection in the window, with the
#                            surgery type, provider and curative-intent flag
#
# Reads : Data/ICON/ncras_og_2015_2022.rds, hes_apc / hes_op extracts
# Writes: og_endoscopy_anchor.rds, og_emresd_anchor.rds, og_surgery_anchor.rds
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

ncras_og <- readRDS(f_ncras_cohort)
hes_apc  <- readRDS(f_hes_apc_extract)
hes_op   <- readRDS(f_hes_op_extract)

# the extracts must carry the fields the anchors need (check_extract is in 01);
# a stale extract is reported clearly rather than failing later.
check_extract(hes_apc, c("STUDY_ID","EPISTART","EPIORDER","ADMIDATE","ADMIMETH",
                         "PROCODE3","SITETRET","EPITYPE", op_cols, opdate_cols),
              "HES APC", f_hes_apc_extract)
check_extract(hes_op, c("STUDY_ID","appt_date","ATTENDED"),
              "HES OP", f_hes_op_extract)

# -----------------------------------------------------------------------------
# Diagnostic endoscopy: APC anchor, then supplement from OP, then combine
# -----------------------------------------------------------------------------
hes_endoscopy_apc <- match_opcs_episodes(hes_apc, opcs_diagnostic_endoscopy,
                                         op_cols, opdate_cols)

endoscopy_anchor_apc <- ncras_og %>%
  select(pseudo_patientid, diagmdy) %>%
  left_join(hes_endoscopy_apc %>%
              select(pseudo_patientid, EPISTART, EPIORDER, opcs4, op_date),
            by = "pseudo_patientid") %>%
  mutate(days_endo_to_dx = as.integer(diagmdy - op_date)) %>%
  filter(!is.na(days_endo_to_dx),
         days_endo_to_dx >= 0,
         days_endo_to_dx <= endoscopy_pre_dx_days) %>%
  arrange(pseudo_patientid, op_date) %>%
  distinct(pseudo_patientid, .keep_all = TRUE) %>%
  rename(endoscopy_date = op_date) %>%
  select(pseudo_patientid, endoscopy_date, days_endo_to_dx) %>%
  mutate(endo_source = "APC")

# OP endoscopy: attended outpatient appointments carrying an endoscopy OPCS code
op_op_cols <- names(hes_op)[str_starts(names(hes_op), "OPERTN_")]
hes_op_endoscopy <- hes_op %>%
  pivot_longer(all_of(op_op_cols), names_to = "op_position", values_to = "opcs_code") %>%
  filter(!is.na(opcs_code), opcs_code != "-") %>%
  mutate(opcs4 = normalise_opcs(opcs_code)) %>%
  filter(opcs4 %in% opcs_diagnostic_endoscopy) %>%
  select(STUDY_ID, appt_date, opcs4)

endoscopy_anchor_op <- ncras_og %>%
  select(pseudo_patientid, diagmdy) %>%
  left_join(hes_op_endoscopy %>%
              rename(pseudo_patientid = STUDY_ID, endoscopy_date = appt_date),
            by = "pseudo_patientid") %>%
  mutate(days_endo_to_dx = as.integer(diagmdy - endoscopy_date)) %>%
  filter(!is.na(days_endo_to_dx),
         days_endo_to_dx >= 0,
         days_endo_to_dx <= endoscopy_pre_dx_days) %>%
  arrange(pseudo_patientid, endoscopy_date) %>%
  distinct(pseudo_patientid, .keep_all = TRUE) %>%
  select(pseudo_patientid, endoscopy_date, days_endo_to_dx) %>%
  mutate(endo_source = "OP")

# combine: APC preferred, OP fills patients APC missed
endoscopy_anchor <- endoscopy_anchor_apc %>%
  bind_rows(endoscopy_anchor_op %>%
              filter(!pseudo_patientid %in% endoscopy_anchor_apc$pseudo_patientid)) %>%
  arrange(pseudo_patientid, endoscopy_date) %>%
  distinct(pseudo_patientid, .keep_all = TRUE)

saveRDS(endoscopy_anchor, f_endoscopy_anchor)
cat("Endoscopy anchor:", n_distinct(endoscopy_anchor$pseudo_patientid),
    "patients (APC", sum(endoscopy_anchor$endo_source == "APC"),
    "/ OP", sum(endoscopy_anchor$endo_source == "OP"), ")\n")

# -----------------------------------------------------------------------------
# EMR / ESD endotherapy: first in the -30d .. +9-month window
# -----------------------------------------------------------------------------
hes_emresd <- match_opcs_episodes(hes_apc, opcs_emresd, op_cols, opdate_cols)

emresd_anchor <- ncras_og %>%
  select(pseudo_patientid, diagmdy) %>%
  left_join(hes_emresd %>% select(pseudo_patientid, EPISTART, EPIORDER, opcs4, op_date),
            by = "pseudo_patientid") %>%
  mutate(days_dx_to_emresd = as.integer(op_date - diagmdy)) %>%
  filter(!is.na(days_dx_to_emresd),
         days_dx_to_emresd >= -30,
         days_dx_to_emresd <= tx_window_days) %>%
  arrange(pseudo_patientid, op_date) %>%
  distinct(pseudo_patientid, .keep_all = TRUE) %>%
  rename(emresd_date = op_date) %>%
  select(pseudo_patientid, emresd_date, days_dx_to_emresd)

saveRDS(emresd_anchor, f_emresd_anchor)
cat("EMR/ESD anchor:", n_distinct(emresd_anchor$pseudo_patientid), "patients\n")

# -----------------------------------------------------------------------------
# Major OG resection: earliest elective episode in the window, with intent
# -----------------------------------------------------------------------------
hes_og_surgery <- match_opcs_episodes(hes_apc, opcs_og_surgery_all, op_cols, opdate_cols)

proc_priority_og <- c("oesophagectomy", "total_gastrectomy", "partial_gastrectomy")
first_or_na <- function(x) if (length(x) == 0) NA_character_ else dplyr::first(x)

# collapse to one row per surgical episode; type hierarchy oesophagectomy >
# total gastrectomy > partial gastrectomy
hes_og_surgery_episodes <- hes_og_surgery %>%
  mutate(surgery_type = unname(opcs_og_surgery_lookup[opcs4])) %>%
  arrange(pseudo_patientid, EPISTART, op_position_n) %>%
  group_by(pseudo_patientid, ADMIDATE, EPISTART, EPIORDER,
           EPITYPE, PROCODE3, SITETRET, ADMIMETH) %>%
  summarise(
    surgery_type = first_or_na(intersect(proc_priority_og, unique(surgery_type))),
    opcs_primary = first_or_na(opcs4),
    all_og_opcs  = paste(unique(opcs4), collapse = "; "),
    surgery_date = as.Date(first_or_na(as.character(op_date[!is.na(op_date)]))),
    .groups      = "drop"
  ) %>%
  mutate(emergency = ADMIMETH %in% admimeth_emerg)

surgery_anchor <- ncras_og %>%
  select(pseudo_patientid, diagmdy, stage_clean) %>%
  left_join(hes_og_surgery_episodes, by = "pseudo_patientid") %>%
  mutate(days_dx_to_surg = as.integer(surgery_date - diagmdy)) %>%
  filter(!is.na(days_dx_to_surg),
         days_dx_to_surg >= -30,
         days_dx_to_surg <= tx_window_days) %>%
  mutate(
    # curative intent: stage 4 + partial gastrectomy is the non-curative case
    curative_surgery = !(stage_clean == "4" & surgery_type == "partial_gastrectomy"),
    surgery_class = case_when(
      surgery_type == "oesophagectomy"                            ~ "oesophagectomy",
      surgery_type %in% c("total_gastrectomy","partial_gastrectomy") ~ "gastrectomy",
      TRUE ~ NA_character_)
  ) %>%
  arrange(pseudo_patientid, surgery_date, EPIORDER) %>%
  distinct(pseudo_patientid, .keep_all = TRUE) %>%
  select(pseudo_patientid, surgery_date, surgery_type, surgery_class,
         opcs_primary, all_og_opcs, PROCODE3, SITETRET,
         days_dx_to_surg, curative_surgery, ADMIMETH, emergency)

saveRDS(surgery_anchor, f_surgery_anchor)
cat("Surgery anchor:", n_distinct(surgery_anchor$pseudo_patientid),
    "patients (curative", sum(surgery_anchor$curative_surgery, na.rm = TRUE), ")\n")
cat("04 complete. Next: 05_derive_comorbidities.R\n")

# ---- optional checks (uncomment to inspect) ---------------------------------
# count(surgery_anchor, surgery_type, surgery_class, curative_surgery) %>% print()
# summary(emresd_anchor$days_dx_to_emresd)