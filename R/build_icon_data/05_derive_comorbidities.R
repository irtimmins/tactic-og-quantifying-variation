# =============================================================================
# 05  Derive comorbidities
# -----------------------------------------------------------------------------
# Builds the RCS Charlson comorbidity index (Armitage & van der Meulen 2010) for
# each patient from the HES APC secondary-diagnosis fields, over two windows:
#   - the 12 months before diagnosis (all Charlson codes), and
#   - a tight band around the surgery / EMR admission (chronic codes only; acute
#     codes are excluded from the diagnosis-to-treatment period).
# Each condition scores 1; the score is capped at 3+ (the RCS grouping).
#
# Reads : Data/ICON/ncras_og_2015_2022.rds, hes_apc extract,
#         og_surgery_anchor.rds, og_emresd_anchor.rds
# Writes: Data/ICON/og_cci_2015_2022.rds
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

# -----------------------------------------------------------------------------
# Comorbidity window choices (kept explicit so the method is auditable)
# -----------------------------------------------------------------------------
lookback_days         <- 365L  # length of the pre-diagnosis window
include_treatment_adm <- TRUE  # also score the operative admission?
tx_pre_days           <- 3L    # operative spell may start a few days early
tx_post_days          <- 7L    # ... and extends a week, not into readmissions
n_diag_fields         <- 7L    # diagnosis positions to scan (RCS used 1-7)
score_malignancy      <- TRUE  # OG appendix scores "any malignancy"
score_metastatic      <- TRUE  # ... and "metastatic solid tumour"
drop_own_site_codes   <- c("C15", "C16")   # the OG site itself is not a comorbidity

# -----------------------------------------------------------------------------
# RCS Charlson code list. `acute` codes only count from a pre-diagnosis
# admission. Prefixes are dot-free to match the HES DIAG_4_ format.
# -----------------------------------------------------------------------------
charlson_icd10 <- list(
  myocardial_infarction    = list(codes = c("I252"), acute = FALSE),
  myocardial_infarction_ac = list(codes = c("I21","I22","I23"), acute = TRUE),
  congestive_heart_failure = list(codes = c("I11","I13","I255","I42","I43","I50","I517"), acute = FALSE),
  peripheral_vascular      = list(codes = c("I70","I71","I72","I73","I770","I771",
                                            "K551","K558","K559","R02","Z958","Z959"), acute = FALSE),
  cerebrovascular          = list(codes = c("G45","G46","I60","I61","I62","I63","I64",
                                            "I65","I66","I67","I68","I69"), acute = FALSE),
  dementia                 = list(codes = c("A810","F00","F01","F02","F03","F051","G30","G31","R54"), acute = FALSE),
  chronic_pulmonary        = list(codes = c("I26","I27","J40","J41","J42","J43","J44","J45",
                                            "J47","J60","J61","J62","J63","J64","J65","J66",
                                            "J67","J684","J701","J703"), acute = FALSE),
  chronic_pulmonary_ac     = list(codes = c("J46"), acute = TRUE),
  rheumatological          = list(codes = c("M05","M06","M09","M120","M315","M32","M33","M34","M35","M36"), acute = FALSE),
  liver                    = list(codes = c("B18","I85","I864","I982","K70","K71","K721",
                                            "K729","K76","R162","Z944"), acute = FALSE),
  diabetes                 = list(codes = c("E10","E11","E12","E13","E14"), acute = FALSE),
  hemiplegia_paraplegia    = list(codes = c("G114","G81","G82","G83"), acute = FALSE),
  renal                    = list(codes = c("I12","I13","N01","N03","N05","N07","N08",
                                            "N18","N25","Z49","Z940","Z992"), acute = FALSE),
  renal_ac                 = list(codes = c("N171","N172","N19"), acute = TRUE),
  aids                     = list(codes = c("B20","B21","B22","B23","B24"), acute = FALSE)
)

# any-malignancy (own site dropped) and metastatic solid tumour, as C-code prefixes
make_c <- function(nums) sprintf("C%02d", nums)
malignancy_codes <- setdiff(
  make_c(c(0:26, 30:34, 37:41, 43, 45:58, 60:76, 80:85, 88, 90:97)),
  drop_own_site_codes)
metastatic_codes <- c("C77", "C78", "C79")
if (score_malignancy)
  charlson_icd10$malignancy <- list(codes = malignancy_codes, acute = FALSE)
if (score_metastatic)
  charlson_icd10$metastatic_solid_tumour <- list(codes = metastatic_codes, acute = FALSE)

# flat lookup; acute sub-lists ("_ac") share the parent condition name
charlson_lookup <- imap_dfr(charlson_icd10, function(x, name) {
  tibble(prefix = x$codes, condition = str_remove(name, "_ac$"), acute = x$acute)
})

# -----------------------------------------------------------------------------
# Inputs and the two comorbidity windows
# -----------------------------------------------------------------------------
ncras_og <- readRDS(f_ncras_cohort)
hes_apc  <- readRDS(f_hes_apc_extract)
cci_diag_cols <- paste0("DIAG_4_", str_pad(1:n_diag_fields, 2, pad = "0"))
# the APC extract must carry STUDY_ID, EPISTART and the diagnosis fields scanned
# for comorbidities (check_extract is in 01).
check_extract(hes_apc, c("STUDY_ID", "EPISTART", cci_diag_cols),
              "HES APC", f_hes_apc_extract)

# treatment date = the admission that generates APC coding: surgery, else EMR.
# Systemic/RT dates are not used, so non-surgical patients keep the strict
# pre-diagnosis window only.
surgery_anchor <- readRDS(f_surgery_anchor) %>%
  transmute(STUDY_ID = as.character(pseudo_patientid), surgery_date = as.Date(surgery_date))
emresd_anchor  <- readRDS(f_emresd_anchor) %>%
  transmute(STUDY_ID = as.character(pseudo_patientid), emresd_date = as.Date(emresd_date))

diag_dates <- ncras_og %>%
  transmute(STUDY_ID = as.character(pseudo_patientid), diagmdy = as.Date(diagmdy)) %>%
  left_join(surgery_anchor, by = "STUDY_ID") %>%
  left_join(emresd_anchor,  by = "STUDY_ID") %>%
  mutate(tx_date        = coalesce(surgery_date, emresd_date),
         lookback_start = diagmdy - lookback_days,
         prev_end       = diagmdy - 1L,
         tx_band_start  = tx_date - tx_pre_days,
         tx_band_end    = tx_date + tx_post_days)

prev_episodes <- hes_apc %>%
  select(STUDY_ID, EPISTART, any_of(cci_diag_cols)) %>%
  mutate(STUDY_ID = as.character(STUDY_ID), EPISTART = as.Date(EPISTART)) %>%
  inner_join(diag_dates, by = "STUDY_ID") %>%
  filter(EPISTART >= lookback_start, EPISTART <= prev_end) %>%
  select(STUDY_ID, any_of(cci_diag_cols))

if (include_treatment_adm) {
  tx_episodes <- hes_apc %>%
    select(STUDY_ID, EPISTART, any_of(cci_diag_cols)) %>%
    mutate(STUDY_ID = as.character(STUDY_ID), EPISTART = as.Date(EPISTART)) %>%
    inner_join(diag_dates, by = "STUDY_ID") %>%
    filter(!is.na(tx_date),
           EPISTART >= tx_band_start, EPISTART <= tx_band_end,
           EPISTART >= diagmdy) %>%   # never before diagnosis (avoids double counting)
    select(STUDY_ID, any_of(cci_diag_cols))
}

# -----------------------------------------------------------------------------
# Pivot diagnoses long and match to the Charlson conditions
# -----------------------------------------------------------------------------
pivot_and_match <- function(episodes, source_label) {
  episodes %>%
    pivot_longer(any_of(cci_diag_cols), names_to = "diag_position", values_to = "icd_code") %>%
    filter(!is.na(icd_code), icd_code != "-", icd_code != "") %>%
    mutate(icd_code = str_remove_all(str_trim(icd_code), "\\."),
           source = source_label,
           p3 = str_sub(icd_code, 1, 3), p4 = str_sub(icd_code, 1, 4)) %>%
    left_join(charlson_lookup %>% filter(nchar(prefix) == 3) %>%
                rename(cond3 = condition, acute3 = acute),
              by = c("p3" = "prefix"), relationship = "many-to-many") %>%
    left_join(charlson_lookup %>% filter(nchar(prefix) == 4) %>%
                rename(cond4 = condition, acute4 = acute),
              by = c("p4" = "prefix"), relationship = "many-to-many") %>%
    mutate(charlson_condition = coalesce(cond4, cond3),
           is_acute           = coalesce(acute4, acute3)) %>%
    filter(!is.na(charlson_condition)) %>%
    select(STUDY_ID, charlson_condition, is_acute, source)
}

prev_matched <- pivot_and_match(prev_episodes, "previous")
all_matched  <- prev_matched
if (include_treatment_adm) {
  tx_matched  <- pivot_and_match(tx_episodes, "treatment") %>% filter(!is_acute)
  all_matched <- bind_rows(prev_matched, tx_matched)
}

# -----------------------------------------------------------------------------
# Collapse to patient level and score (each condition weight 1, capped at 3+)
# -----------------------------------------------------------------------------
cci_hits <- all_matched %>%
  distinct(STUDY_ID, charlson_condition) %>%
  group_by(STUDY_ID) %>%
  summarise(cci_n_conditions = n_distinct(charlson_condition),
            cci_conditions   = paste(sort(unique(charlson_condition)), collapse = "; "),
            .groups = "drop")

og_cci <- tibble(STUDY_ID = as.character(unique(ncras_og$pseudo_patientid))) %>%
  left_join(cci_hits, by = "STUDY_ID") %>%
  mutate(cci_n_conditions = replace_na(cci_n_conditions, 0L),
         cci_conditions   = replace_na(cci_conditions, "none"),
         rcs_ch_score     = pmin(cci_n_conditions, 3L),
         cci_any          = as.integer(cci_n_conditions >= 1),
         cci_group        = factor(case_when(
           cci_n_conditions == 0 ~ "0", cci_n_conditions == 1 ~ "1",
           cci_n_conditions == 2 ~ "2", TRUE ~ "3+"),
           levels = c("0","1","2","3+"))) %>%
  rename(pseudo_patientid = STUDY_ID)

saveRDS(og_cci, f_cci)
cat("Saved", f_cci, "(", nrow(og_cci), "patients ).",
    "Next: 06_derive_sact_rtds.R\n")

# ---- optional checks (uncomment to inspect) ---------------------------------
# print(round(prop.table(table(og_cci$cci_group)), 3))
# all_matched %>% distinct(STUDY_ID, charlson_condition) %>%
#   count(charlson_condition, sort = TRUE) %>% print(n = 20)