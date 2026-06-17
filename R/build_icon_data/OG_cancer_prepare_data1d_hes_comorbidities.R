# =============================================================================
# OG cancer - Royal College of Surgeons Charlson comorbidity index (CCI)
# -----------------------------------------------------------------------------
# Builds a patient-level comorbidity lookup for the OG cohort from HES APC
# secondary diagnoses, following the OG audit appendix definition.
#
# Method follows the audit's cited reference:
#   Armitage JN, van der Meulen JH. Identifying co-morbidity in surgical
#   patients using administrative data with the Royal College of Surgeons
#   Charlson Score. Br J Surg 2010;97:772-81. doi:10.1002/bjs.6930
#
# Appendix 2 conditions (all included here): myocardial infarction, congestive
# cardiac failure, peripheral vascular disease, cerebrovascular disease,
# dementia, chronic pulmonary disease, rheumatological disease, liver disease,
# diabetes mellitus, hemiplegia/paraplegia, renal disease, any malignancy,
# metastatic solid tumour, AIDS/HIV.
#
# Comorbidity window
#   The RCS score is built around the surgical admission. In this cohort, as in
#   most surgical cancer cohorts, a large share of pre-existing comorbidity is
#   first coded at the treatment (usually operative) admission rather than in
#   prior admissions. The window therefore has two parts:
#     - before diagnosis  [diagmdy - 365, diagmdy - 1]  : all codes count
#     - diagnosis to treatment [diagmdy, treatment_date] : non-acute codes only
#   Acute codes (e.g. acute MI, acute kidney injury) only ever count from the
#   pre-diagnosis part, so a peri-operative event is treated as a complication,
#   not a pre-existing comorbidity.
#
# OG-specific points relative to the bowel (NDRS) version:
#   * Malignancy and metastatic tumour ARE scored (the appendix lists them).
#     The cohort's own site codes (C15*, C16*) are dropped from the malignancy
#     match so the index OG cancer does not count as its own comorbidity.
#   * AIDS/HIV is in the definition but, per the appendix note, HIV/AIDS
#     diagnoses are legally removed from HES APC, so it will essentially never
#     be detected here. The match is kept for completeness.
#
# Run order: after scripts 1-2 (needs hes_apc_og with diagnosis fields, plus the
# treatment anchors for the treatment date). Output is joined onto the cohort in
# script 4.
# =============================================================================

library(tidyverse)

base_dir <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/"

# -----------------------------------------------------------------------------
# 0. Choices (kept explicit so the method is auditable)
# -----------------------------------------------------------------------------
lookback_days        <- 365L   # length of the pre-diagnosis window
include_treatment_adm<- TRUE   # also score the treatment admission?
tx_pre_days          <- 3L     # days before the treatment date the operative
                               # spell may start (admission a day or two early)
tx_post_days         <- 7L     # days after, to catch the rest of the spell's
                               # coding without reaching later readmissions
n_diag_fields        <- 7L    # diagnosis positions to scan (RCS paper used 1-7)
score_malignancy     <- TRUE   # OG appendix lists "any malignancy"
score_metastatic     <- TRUE   # OG appendix lists "metastatic solid tumour"
drop_own_site_codes  <- c("C15", "C16")  # excluded from the malignancy match

# -----------------------------------------------------------------------------
# 1. RCS Charlson code list (Armitage & van der Meulen 2010)
#    `acute` codes only count from a pre-diagnosis admission, never the
#    diagnosis-to-treatment part. Codes are dot-free 3- or 4-character prefixes
#    to match the HES DIAG_4_ format.
# -----------------------------------------------------------------------------
charlson_icd10 <- list(
  myocardial_infarction    = list(codes = c("I252"),                                  acute = FALSE),
  myocardial_infarction_ac = list(codes = c("I21", "I22", "I23"),                     acute = TRUE),
  congestive_heart_failure = list(codes = c("I11", "I13", "I255", "I42", "I43",
                                            "I50", "I517"),                           acute = FALSE),
  peripheral_vascular      = list(codes = c("I70", "I71", "I72", "I73", "I770",
                                            "I771", "K551", "K558", "K559", "R02",
                                            "Z958", "Z959"),                          acute = FALSE),
  cerebrovascular          = list(codes = c("G45", "G46", "I60", "I61", "I62", "I63",
                                            "I64", "I65", "I66", "I67", "I68", "I69"),acute = FALSE),
  dementia                 = list(codes = c("A810", "F00", "F01", "F02", "F03",
                                            "F051", "G30", "G31", "R54"),             acute = FALSE),
  chronic_pulmonary        = list(codes = c("I26", "I27", "J40", "J41", "J42", "J43",
                                            "J44", "J45", "J47", "J60", "J61", "J62",
                                            "J63", "J64", "J65", "J66", "J67", "J684",
                                            "J701", "J703"),                          acute = FALSE),
  chronic_pulmonary_ac     = list(codes = c("J46"),                                   acute = TRUE),
  rheumatological          = list(codes = c("M05", "M06", "M09", "M120", "M315",
                                            "M32", "M33", "M34", "M35", "M36"),        acute = FALSE),
  liver                    = list(codes = c("B18", "I85", "I864", "I982", "K70",
                                            "K71", "K721", "K729", "K76", "R162",
                                            "Z944"),                                  acute = FALSE),
  diabetes                 = list(codes = c("E10", "E11", "E12", "E13", "E14"),       acute = FALSE),
  hemiplegia_paraplegia    = list(codes = c("G114", "G81", "G82", "G83"),             acute = FALSE),
  renal                    = list(codes = c("I12", "I13", "N01", "N03", "N05", "N07",
                                            "N08", "N18", "N25", "Z49", "Z940",
                                            "Z992"),                                  acute = FALSE),
  renal_ac                 = list(codes = c("N171", "N172", "N19"),                   acute = TRUE),
  aids                     = list(codes = c("B20", "B21", "B22", "B23", "B24"),       acute = FALSE)
)

# malignancy (any) and metastatic solid tumour, built as 3-char C-code prefixes
make_c <- function(nums) sprintf("C%02d", nums)
malignancy_codes <- make_c(c(0:26, 30:34, 37:41, 43, 45:58, 60:76, 80:85, 88, 90:97))
malignancy_codes <- setdiff(malignancy_codes, drop_own_site_codes)  # OG: drop own site
metastatic_codes <- c("C77", "C78", "C79")

if (score_malignancy) charlson_icd10$malignancy <-
  list(codes = malignancy_codes, acute = FALSE)
if (score_metastatic) charlson_icd10$metastatic_solid_tumour <-
  list(codes = metastatic_codes, acute = FALSE)

# flat lookup, with a flag for whether each prefix is acute
charlson_lookup <- imap_dfr(charlson_icd10, function(x, name) {
  base_cond <- str_remove(name, "_ac$")   # acute sub-lists share the parent condition
  tibble(prefix = x$codes, condition = base_cond, acute = x$acute)
})

# -----------------------------------------------------------------------------
# 2. Inputs
# -----------------------------------------------------------------------------
ncras_og    <- readRDS(paste0(base_dir, "ncras_og_2015_2022.rds"))
hes_apc_raw <- readRDS(paste0(base_dir, "hes_apc_og_2014_2022.rds"))

diag_cols <- paste0("DIAG_4_", str_pad(1:n_diag_fields, 2, pad = "0"))
if (!all(diag_cols %in% names(hes_apc_raw))) {
  stop("hes_apc_og is missing DIAG_4_ fields. Re-run script 1 with the ",
       "diagnosis columns added to hes_cols_select.")
}

# treatment date per patient, from the anchors built in scripts 1-2.
# Priority is the admission that actually generates APC coding: surgery first,
# then EMR/ESD, then the systemic/RT date as a fallback window end.
read_anchor <- function(file, datecol) {
  readRDS(paste0(base_dir, file)) %>%
    transmute(STUDY_ID = as.character(pseudo_patientid),
              !!datecol := as.Date(.data[[datecol]]))
}
surgery_anchor <- read_anchor("og_surgery_anchor_2015_2022.rds", "surgery_date")
emresd_anchor  <- read_anchor("OG_emresd_anchor.rds",            "emresd_date")
sact_anchor    <- read_anchor("og_sact_anchor_2015_2022.rds",    "sact_date")
rt_anchor      <- read_anchor("rt_anchor_og.rds",                "rt_date")

tx_dates <- ncras_og %>%
  transmute(STUDY_ID = as.character(pseudo_patientid)) %>%
  left_join(surgery_anchor, by = "STUDY_ID") %>%
  left_join(emresd_anchor,  by = "STUDY_ID") %>%
  left_join(sact_anchor,    by = "STUDY_ID") %>%
  left_join(rt_anchor,      by = "STUDY_ID") %>%
  # only surgery and EMR/ESD are APC admissions that code comorbidity; the
  # systemic/RT dates are dropped so non-surgical patients keep the strict
  # pre-diagnosis window.
  mutate(tx_date = coalesce(surgery_date, emresd_date)) %>%
  select(STUDY_ID, tx_date)

# -----------------------------------------------------------------------------
# 3. Comorbidity windows
#    prev:      12 months up to the day before diagnosis (all codes)
#    treatment: a tight band around the surgery/EMR admission (non-acute only)
# -----------------------------------------------------------------------------
diag_dates <- ncras_og %>%
  transmute(STUDY_ID = as.character(pseudo_patientid),
            diagmdy  = as.Date(diagmdy)) %>%
  left_join(tx_dates, by = "STUDY_ID") %>%
  mutate(
    lookback_start = diagmdy - lookback_days,
    prev_end       = diagmdy - 1L,
    # tight band around the operative admission; NA for patients with no
    # surgery/EMR, who then contribute only the pre-diagnosis window
    tx_band_start  = tx_date - tx_pre_days,
    tx_band_end    = tx_date + tx_post_days
  )

prev_episodes <- hes_apc_raw %>%
  select(STUDY_ID, EPISTART, any_of(diag_cols)) %>%
  mutate(STUDY_ID = as.character(STUDY_ID), EPISTART = as.Date(EPISTART)) %>%
  inner_join(diag_dates, by = "STUDY_ID") %>%
  filter(EPISTART >= lookback_start, EPISTART <= prev_end) %>%
  select(STUDY_ID, any_of(diag_cols))

if (include_treatment_adm) {
  tx_episodes <- hes_apc_raw %>%
    select(STUDY_ID, EPISTART, any_of(diag_cols)) %>%
    mutate(STUDY_ID = as.character(STUDY_ID), EPISTART = as.Date(EPISTART)) %>%
    inner_join(diag_dates, by = "STUDY_ID") %>%
    filter(!is.na(tx_date),
           EPISTART >= tx_band_start, EPISTART <= tx_band_end,
           EPISTART >= diagmdy) %>%   # never before diagnosis (avoids overlap with prev)
    select(STUDY_ID, any_of(diag_cols))
}

# coverage checks
cat("Patients with >=1 HES admission before diagnosis:",
    n_distinct(prev_episodes$STUDY_ID), "\n")
if (include_treatment_adm) {
  cat("Patients with >=1 HES admission in the operative band:",
      n_distinct(tx_episodes$STUDY_ID), "\n")
}

# -----------------------------------------------------------------------------
# 4. Pivot diagnoses long and match to the Charlson conditions
# -----------------------------------------------------------------------------
pivot_and_match <- function(episodes, source_label) {
  episodes %>%
    pivot_longer(any_of(diag_cols), names_to = "diag_position", values_to = "icd_code") %>%
    filter(!is.na(icd_code), icd_code != "-", icd_code != "") %>%
    mutate(
      icd_code = str_remove_all(str_trim(icd_code), "\\."),
      source   = source_label,
      p3 = str_sub(icd_code, 1, 3),
      p4 = str_sub(icd_code, 1, 4)
    ) %>%
    left_join(charlson_lookup %>% filter(nchar(prefix) == 3) %>%
                rename(cond3 = condition, acute3 = acute),
              by = c("p3" = "prefix"), relationship = "many-to-many") %>%
    left_join(charlson_lookup %>% filter(nchar(prefix) == 4) %>%
                rename(cond4 = condition, acute4 = acute),
              by = c("p4" = "prefix"), relationship = "many-to-many") %>%
    mutate(
      charlson_condition = coalesce(cond4, cond3),
      is_acute           = coalesce(acute4, acute3)
    ) %>%
    filter(!is.na(charlson_condition)) %>%
    select(STUDY_ID, charlson_condition, is_acute, source)
}

prev_matched <- pivot_and_match(prev_episodes, "previous")
all_matched  <- prev_matched

if (include_treatment_adm) {
  # acute conditions are excluded from the diagnosis-to-treatment part
  tx_matched  <- pivot_and_match(tx_episodes, "treatment") %>% filter(!is_acute)
  all_matched <- bind_rows(prev_matched, tx_matched)
}

cat("Charlson-relevant rows - previous:", nrow(prev_matched))
if (include_treatment_adm) cat(" | treatment:", nrow(tx_matched))
cat("\n")

# -----------------------------------------------------------------------------
# 5. Collapse to patient level and score (RCS: each condition weight 1)
# -----------------------------------------------------------------------------
cci_hits <- all_matched %>%
  distinct(STUDY_ID, charlson_condition) %>%
  group_by(STUDY_ID) %>%
  summarise(
    cci_n_conditions = n_distinct(charlson_condition),
    cci_conditions   = paste(sort(unique(charlson_condition)), collapse = "; "),
    .groups = "drop"
  )

og_cci <- tibble(STUDY_ID = as.character(unique(ncras_og$pseudo_patientid))) %>%
  left_join(cci_hits, by = "STUDY_ID") %>%
  mutate(
    cci_n_conditions = replace_na(cci_n_conditions, 0L),
    cci_conditions   = replace_na(cci_conditions, "none"),
    rcs_ch_score     = pmin(cci_n_conditions, 3L),    # RCS score capped at 3+
    cci_any          = as.integer(cci_n_conditions >= 1),
    cci_group        = factor(
      case_when(
        cci_n_conditions == 0 ~ "0",
        cci_n_conditions == 1 ~ "1",
        cci_n_conditions == 2 ~ "2",
        TRUE                  ~ "3+"
      ),
      levels = c("0", "1", "2", "3+")
    )
  ) %>%
  rename(pseudo_patientid = STUDY_ID)

cat("\nRCS Charlson group distribution:\n")
print(table(og_cci$cci_group))
cat("\nProportions:\n")
print(round(prop.table(table(og_cci$cci_group)), 3))

cat("\nMost frequent conditions:\n")
all_matched %>% distinct(STUDY_ID, charlson_condition) %>%
  count(charlson_condition, sort = TRUE) %>% print(n = 20)

# -----------------------------------------------------------------------------
# 6. Save the patient-level lookup (joined onto the cohort in script 4)
# -----------------------------------------------------------------------------
saveRDS(og_cci, paste0(base_dir, "og_cci_2015_2022.rds"))
cat("\nSaved og_cci_2015_2022.rds (", nrow(og_cci), "patients )\n")

