# =============================================================================
# OG Cancer Waiting Times - PRE-CWT cohort builder + format spec
# -----------------------------------------------------------------------------
# Purpose
#   1. Define the column spec for the pre-CWT cohort (Table A)
#      and the raw CWT records (Table B) that the merging script expects.
#   2. Rebuild the pre-CWT cohort from existing ICON-derived objects in
#      exactly that format, so the CWT merge runs against a known-conformant
#      frame (and so the same spec can drive synthetic-data generation later).
#
# Clock-stop convention: a SINGLE curative-intent treatment date, `first_tx_date`
#   (NOT a separate first_curative_tx_date). Palliative / no-treatment pathways
#   get NA, which is intended for the curative waiting-times analysis.
#
# Run order: scripts 1-2 (ncras_og + anchors) must have been run first.
# =============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(tibble)

# -----------------------------------------------------------------------------
# 0. Globals (must match the rest of the pipeline)
# -----------------------------------------------------------------------------
tx_window_days <- 270L   # post-diagnosis treatment window (9 months)

og_icd10 <- c(
  "C15","C150","C151","C152","C153","C154","C155","C158","C159",
  "C16","C160","C161","C162","C163","C164","C165","C166","C168","C169"
)

# Canonical tx_pathway levels (the only permitted values)
tx_pathway_levels <- c(
  "EMR/ESD only", "EMR/ESD then surgery",
  "Surgery + neoadjuvant chemoRT", "Surgery + neoadjuvant chemo",
  "Surgery + neoadjuvant RT", "Surgery + adjuvant chemo",
  "Surgery only", "Surgery + other",
  "Definitive chemoRT", "Curative RT only",
  "Palliative chemo + RT", "SACT only", "Palliative RT only",
  "No treatment recorded"
)

# =============================================================================
# 1. SPEC MANIFESTS  (data dictionaries)
# =============================================================================
# tier:      "required" = needed for the CWT merge and/or core script-4 analysis
#            "core"     = needed to regenerate pathway/trust/PIs & realism
# missing_ok: TRUE if the column may legitimately contain NA for some patients

pre_cwt_spec <- tribble(
  ~name,                               ~type,        ~tier,      ~missing_ok, ~notes,
  # --- linkage & cohort ------------------------------------------------------
  "pseudo_patientid",                  "character",  "required", FALSE, "Join key to CWT; unique per row",
  "pseudo_tumourid",                   "character",  "core",     FALSE, "Tumour-level key (COSD join)",
  "diagmdy",                           "Date",       "required", FALSE, "Diagnosis date; anchors all wt_dx_*",
  "ydiag",                             "integer",    "core",     FALSE, "Diagnosis year (>=2015)",
  # --- tumour / clinical -----------------------------------------------------
  "tumour_site_grp",                   "character",  "core",     FALSE, "'oesophageal' | 'gastric'",
  "cancer_subtype",                    "character",  "core",     TRUE,  "'Oes SCC' | 'Oes ACA' | 'Gast'",
  "stage_clean",                       "character",  "core",     FALSE, "'1' | '2' | '3' (stage 1-3 cohort)",
  "morphology_num",                    "integer",    "core",     TRUE,  "ICD-O morphology (drives subtype)",
  # --- demographics / geography ---------------------------------------------
  "sex",                               "integer",    "core",     TRUE,  "1 male, 2 female",
  "agediag",                           "numeric",    "core",     TRUE,  "Age at diagnosis (years)",
  "ethnicity_group_broad",             "character",  "core",     TRUE,  "Broad ethnicity group",
  "NHSE_reversed_imd_quintile_lsoas",  "character",  "required", TRUE,  "IMD quintile 1-5 (deprivation analysis)",
  "canalliance_2024_code",             "character",  "core",     TRUE,  "Cancer Alliance code",
  "canalliance_2024_name",             "character",  "core",     TRUE,  "Cancer Alliance name",
  "lsoa11_code",                       "character",  "core",     TRUE,  "LSOA 2011 code",
  # --- route -----------------------------------------------------------------
  "route_combined",                    "factor",     "required", FALSE, "Diagnosis route; incl 'TWW','GP referral','Emergency presentation'",
  "emergency_admission",               "integer",    "core",     TRUE,  "0/1 emergency presentation flag",
  # --- COSD ------------------------------------------------------------------
  "ps_num",                            "integer",    "core",     TRUE,  "Performance status 0-4",
  "cnsinvolved",                       "integer",    "core",     TRUE,  "CNS involvement 0/1",
  # --- diagnosing/first organisation ----------------------------------------
  "diag_hosp",                         "character",  "core",     TRUE,  "Diagnosing hospital code",
  "diag_trust",                        "character",  "core",     TRUE,  "Diagnosing trust code",
  "first_trust",                       "character",  "core",     TRUE,  "First-seen trust code",
  # --- treatment anchor dates ------------------------------------------------
  "endoscopy_date",                    "Date",       "required", TRUE,  "Diagnostic endoscopy; PI3 clock start",
  "emresd_date",                       "Date",       "core",     TRUE,  "First EMR/ESD",
  "surgery_date",                      "Date",       "required", TRUE,  "First elective resection",
  "sact_date",                         "Date",       "required", TRUE,  "First in-window SACT regimen",
  "rt_date",                           "Date",       "required", TRUE,  "First in-window RT course",
  "first_tx_date",                     "Date",       "required", TRUE,  "Curative clock-stop (pathway-based); NA if palliative/none",
  # --- surgery attributes ----------------------------------------------------
  "surgery_type",                      "character",  "core",     TRUE,  "oesophagectomy | total_/partial_gastrectomy",
  "surgery_class",                     "character",  "core",     TRUE,  "'oesophagectomy' | 'gastrectomy'",
  "curative_surgery",                  "logical",    "core",     TRUE,  "FALSE only for stage4 + partial gastrectomy",
  "opcs_primary",                      "character",  "core",     TRUE,  "Primary OPCS-4 of resection",
  "PROCODE3",                          "character",  "core",     TRUE,  "HES provider (surgery trust source)",
  "SITETRET",                          "character",  "core",     TRUE,  "HES site of treatment",
  # --- RT attributes ---------------------------------------------------------
  "rt_curative",                       "logical",    "core",     TRUE,  "Curative by dose/fractionation schedule",
  "rt_dose",                           "numeric",    "core",     TRUE,  "Prescribed dose (Gy)",
  "rt_fractions",                      "integer",    "core",     TRUE,  "Prescribed fractions",
  "ORGCODEPROVIDER",                   "character",  "core",     TRUE,  "RTDS provider (RT trust source)",
  # --- SACT attributes -------------------------------------------------------
  "BENCHMARK_GROUP",                   "character",  "core",     TRUE,  "SACT regimen classification",
  "benchmark_group_lwr",               "character",  "core",     TRUE,  "Lower-cased benchmark group",
  "INTENT_OF_TREATMENT_V3",            "character",  "core",     TRUE,  "SACT intent",
  "CHEMO_RADIATION",                   "character",  "core",     TRUE,  "ChemoRT flag from SACT",
  "ORGANISATION_CODE_OF_PROVIDER",     "character",  "core",     TRUE,  "SACT provider",
  # --- interval days  -----------------------------------------
  "days_endo_to_dx",                   "integer",    "core",     TRUE,  "Endoscopy->diagnosis offset",
  "days_dx_to_emresd",                 "integer",    "core",     TRUE,  "Diagnosis->EMR/ESD",
  "days_dx_to_surg",                   "integer",    "core",     TRUE,  "Diagnosis->surgery",
  "days_dx_to_sact",                   "integer",    "core",     TRUE,  "Diagnosis->SACT",
  "days_dx_to_rt",                     "integer",    "core",     TRUE,  "Diagnosis->RT",
  # --- presence / sequencing flags ------------------------------------------
  "had_emresd",                        "logical",    "core",     FALSE, "",
  "had_surgery",                       "logical",    "core",     FALSE, "",
  "had_curative_surgery",              "logical",    "core",     FALSE, "",
  "had_sact",                          "logical",    "core",     FALSE, "",
  "had_rt",                            "logical",    "core",     FALSE, "",
  "had_curative_rt",                   "logical",    "core",     FALSE, "",
  "had_palliative_rt",                 "logical",    "core",     FALSE, "",
  "sact_before_surgery",              "logical",    "core",     FALSE, "",
  "sact_after_surgery",               "logical",    "core",     FALSE, "",
  "rt_before_surgery",                "logical",    "core",     FALSE, "",
  "rt_after_surgery",                 "logical",    "core",     FALSE, "",
  "concurrent_chemo_rt",              "logical",    "core",     FALSE, "SACT & curative RT within 14d",
  "received_curative_tx",             "logical",    "core",     FALSE, "EMR/ESD | curative surgery | curative RT",
  # --- pathway & trust -------------------------------------------------------
  "tx_pathway",     "character", "required", FALSE, "One of tx_pathway_levels",
  "tx_trust",       "character", "core",     TRUE,  "3-char trust of curative tx (PROCODE3 or provider)",
  "change_trust",   "logical",   "core",     TRUE,  "TRUE = diagnosis trust differs from treatment trust; NA if no curative tx",
  # --- waiting times (pre-CWT) ----------------------------------------------
  "wt_dx_to_tx",                       "integer",    "required", TRUE,  "Diagnosis->first_tx_date",
  "wt_endo_to_tx",                     "integer",    "required", TRUE,  "Endoscopy->first_tx_date",
  "wt_dx_to_surg",                     "integer",    "core",     TRUE,  "",
  "wt_endo_to_surg",                   "integer",    "core",     TRUE,  "",
  "wt_dx_to_sact",                     "integer",    "core",     TRUE,  "",
  "wt_endo_to_sact",                   "integer",    "core",     TRUE,  "",
  "wt_sact_to_surg",                   "integer",    "core",     TRUE,  "Neoadjuvant interval",
  "wt_surg_to_sact",                   "integer",    "core",     TRUE,  "Adjuvant interval",
  "wt_dx_to_rt",                       "integer",    "core",     TRUE,  "",
  "wt_endo_to_rt",                     "integer",    "core",     TRUE,  "",
  "wt_rt_to_surg",                     "integer",    "core",     TRUE,  "",
  # --- survival --------------------------------------------------------------
  "finmdy",                            "Date",       "core",     TRUE,  "Date of death or last follow-up (censor date)",
  "died",                              "integer",    "core",     FALSE, "0/1 death indicator"
)

# Raw CWT records spec (Table B) - pre-parsing, as in the partitioned dataset.
# Dates are CHARACTER "dd/mm/yyyy" so the merge script's as.Date() parsing runs.
cwt_spec <- tribble(
  ~name,                 ~type,       ~tier,      ~missing_ok, ~notes,
  "pseudo_patientid",    "character", "required", FALSE, "Join key; subset must intersect cohort IDs",
  "site_icd10",          "character", "required", FALSE, "C15x/C16x; filtered against og_icd10",
  "modality",            "character", "required", FALSE, "CWT modality code (see modality_codes)",
  "crtp_date",           "character", "required", TRUE,  "dd/mm/yyyy; referral / clock start",
  "date_first_seen",     "character", "core",     TRUE,  "dd/mm/yyyy; first seen",
  "mdt_date",            "character", "required", TRUE,  "dd/mm/yyyy; -> cwt_mdt_date",
  "treat_period_start",  "character", "required", TRUE,  "dd/mm/yyyy; -> cwt_dtt_date (decision to treat)",
  "treat_start",         "character", "required", TRUE,  "dd/mm/yyyy; -> cwt_treat_date (first treatment)"
)

# CWT modality code reference (documentation for generation/filtering)
modality_codes <- tribble(
  ~code, ~meaning,
  "01", "Surgery", "02", "Anti-cancer drug (chemo)", "03", "Radiotherapy",
  "04", "Concurrent chemoRT", "05", "Other (incl active monitoring)",
  "06", "Brachytherapy", "07", "Surgery + drug", "08", "Surgery + RT",
  "09", "Surgery + chemoRT", "23", "Endoscopic EMR/ESD (post-2020 only)",
  "24", "Endoscopic + other", "97", "Excluded", "98", "Excluded", "99", "Excluded"
)

# =============================================================================
# 2. Load existing ICON-derived objects
# =============================================================================
library(here)
icon_dir <- here("Data", "ICON")        # real ICON inputs (read-only source)
syn_dir  <- here("Data", "synthetic")   # synthetic-build outputs
dir.create(syn_dir, recursive = TRUE, showWarnings = FALSE)

# The spec object is built purely from the manifests above - no data needed -
# so save it now. This lets the rest of the synthetic pipeline run even when the
# full ICON anchor set is not present in Data/ICON.
saveRDS(list(pre_cwt_spec = pre_cwt_spec, cwt_spec = cwt_spec,
             modality_codes = modality_codes,
             tx_pathway_levels = tx_pathway_levels,
             og_icd10 = og_icd10, tx_window_days = tx_window_days),
        file.path(syn_dir, "og_pipeline_spec.rds"))
cat("Saved og_pipeline_spec.rds to Data/synthetic\n")

# The rest of this script rebuilds the REAL pre-CWT cohort from the ICON anchor
# objects, purely as an on-server conformance check. It is not consumed by the
# synthetic chain, so it only runs when every anchor file is present.
anchor_files <- c("ncras_og_2015_2022.rds", "OG_endoscopy_anchor_combined.rds",
                  "OG_emresd_anchor.rds", "og_surgery_anchor_2015_2022.rds",
                  "og_sact_anchor_2015_2022.rds", "rt_anchor_og.rds")
have_anchors <- all(file.exists(file.path(icon_dir, anchor_files)))

if (!have_anchors) {
  message("ICON anchor objects not found in Data/ICON - skipping the real ",
          "cohort rebuild. Spec is saved; synthetic generation can proceed.")
} else {

ncras_og                  <- readRDS(file.path(icon_dir, "ncras_og_2015_2022.rds"))
endoscopy_anchor_combined <- readRDS(file.path(icon_dir, "OG_endoscopy_anchor_combined.rds"))
emresd_anchor             <- readRDS(file.path(icon_dir, "OG_emresd_anchor.rds"))
surgery_anchor            <- readRDS(file.path(icon_dir, "og_surgery_anchor_2015_2022.rds"))
sact_anchor               <- readRDS(file.path(icon_dir, "og_sact_anchor_2015_2022.rds"))
rt_anchor                 <- readRDS(file.path(icon_dir, "rt_anchor_og.rds"))

# =============================================================================
# 3. Build the pre-CWT cohort  (joins + derivations)
# =============================================================================
og_cohort <- ncras_og %>%
  left_join(endoscopy_anchor_combined %>%
              select(pseudo_patientid, endoscopy_date, days_endo_to_dx),
            by = "pseudo_patientid") %>%
  left_join(emresd_anchor %>%
              select(pseudo_patientid, emresd_date, days_dx_to_emresd),
            by = "pseudo_patientid") %>%
  left_join(surgery_anchor %>%
              select(pseudo_patientid, surgery_date, surgery_type, surgery_class,
                     opcs_primary, PROCODE3, SITETRET, days_dx_to_surg,
                     curative_surgery),
            by = "pseudo_patientid") %>%
  left_join(sact_anchor %>%
              select(pseudo_patientid, sact_date, days_dx_to_sact,
                     BENCHMARK_GROUP, benchmark_group_lwr,
                     INTENT_OF_TREATMENT_V3, CHEMO_RADIATION,
                     ORGANISATION_CODE_OF_PROVIDER),
            by = "pseudo_patientid") %>%
  left_join(rt_anchor %>%
              select(pseudo_patientid, rt_date, rt_curative,
                     rt_dose, rt_fractions, days_dx_to_rt, ORGCODEPROVIDER),
            by = "pseudo_patientid") %>%
  mutate(
    # --- presence flags ------------------------------------------------------
    had_emresd           = !is.na(emresd_date),
    had_surgery          = !is.na(surgery_date),
    had_curative_surgery = !is.na(surgery_date) & curative_surgery == TRUE,
    had_sact             = !is.na(sact_date),
    had_rt               = !is.na(rt_date),
    had_curative_rt      = !is.na(rt_date) & rt_curative == TRUE,
    had_palliative_rt    = !is.na(rt_date) & rt_curative == FALSE,
    # --- sequencing ----------------------------------------------------------
    sact_before_surgery  = had_sact & had_surgery & sact_date < surgery_date,
    sact_after_surgery   = had_sact & had_surgery & sact_date > surgery_date,
    rt_before_surgery    = had_rt   & had_surgery & rt_date   < surgery_date,
    rt_after_surgery     = had_rt   & had_surgery & rt_date   > surgery_date,
    concurrent_chemo_rt  = had_sact & had_curative_rt &
      abs(as.integer(sact_date - rt_date)) <= 14,
    received_curative_tx = had_emresd | had_curative_surgery | had_curative_rt,
    # --- pathway -------------------------------------------------------------
    tx_pathway = case_when(
      had_emresd & !had_surgery & !had_sact & !concurrent_chemo_rt ~ "EMR/ESD only",
      had_emresd & had_surgery                                     ~ "EMR/ESD then surgery",
      had_surgery & sact_before_surgery & rt_before_surgery        ~ "Surgery + neoadjuvant chemoRT",
      had_surgery & sact_before_surgery & !rt_before_surgery       ~ "Surgery + neoadjuvant chemo",
      had_surgery & rt_before_surgery & !sact_before_surgery       ~ "Surgery + neoadjuvant RT",
      had_surgery & sact_after_surgery & !sact_before_surgery      ~ "Surgery + adjuvant chemo",
      had_surgery & !had_sact & !concurrent_chemo_rt               ~ "Surgery only",
      had_surgery                                                  ~ "Surgery + other",
      !had_surgery & had_curative_rt & had_sact                    ~ "Definitive chemoRT",
      !had_surgery & had_curative_rt & !had_sact                   ~ "Curative RT only",
      !had_surgery & had_palliative_rt & had_sact                  ~ "Palliative chemo + RT",
      !had_surgery & had_sact & !had_curative_rt                   ~ "SACT only",
      !had_surgery & had_palliative_rt & !had_sact                 ~ "Palliative RT only",
      TRUE                                                         ~ "No treatment recorded"
    ),
    # --- single curative clock-stop: first_tx_date ---------------------------
    # Surgery/EMR-ESD -> surgery/emresd date; neoadjuvant chemo -> sact_date;
    # neoadjuvant RT/chemoRT -> earliest of RT/SACT; definitive chemoRT ->
    # earliest of RT/SACT; curative RT only -> rt_date; else NA.
    first_tx_date = case_when(
      tx_pathway %in% c("EMR/ESD only", "EMR/ESD then surgery") ~ emresd_date,
      tx_pathway == "Surgery + neoadjuvant chemoRT"            ~ pmin(sact_date, rt_date, na.rm = TRUE),
      tx_pathway == "Surgery + neoadjuvant RT"                 ~ rt_date,
      tx_pathway == "Surgery + neoadjuvant chemo"              ~ sact_date,
      tx_pathway %in% c("Surgery + adjuvant chemo",
                        "Surgery only", "Surgery + other")     ~ surgery_date,
      tx_pathway == "Definitive chemoRT"                       ~ pmin(sact_date, rt_date, na.rm = TRUE),
      tx_pathway == "Curative RT only"                         ~ rt_date,
      TRUE                                                     ~ as.Date(NA)
    ),
    # --- trust of curative treatment (3-char) --------------------------------
    tx_trust = case_when(
      tx_pathway %in% c("EMR/ESD only", "EMR/ESD then surgery",
                        "Surgery + neoadjuvant chemo", "Surgery + adjuvant chemo",
                        "Surgery only", "Surgery + other")     ~ substr(PROCODE3, 1, 3),
      tx_pathway %in% c("Surgery + neoadjuvant chemoRT",
                        "Surgery + neoadjuvant RT")            ~ substr(ORGCODEPROVIDER, 1, 3),
      tx_pathway %in% c("Definitive chemoRT", "Curative RT only") ~ substr(ORGCODEPROVIDER, 1, 3),
      TRUE                                                     ~ NA_character_
    ),
    # Did the patient move trust between diagnosis and curative treatment?
    # diag_trust is trimmed to its first 3 chars to align with tx_trust (already
    # 3-char). The comparison returns NA automatically when tx_trust is NA
    # (palliative / no curative treatment) or diag_trust is missing - no guard
    # needed. TRUE = trusts differ, FALSE = same trust.
    change_trust = substr(diag_trust, 1, 3) != tx_trust,
    # --- waiting times (deduplicated) ----------------------------------------
    wt_dx_to_tx     = as.integer(first_tx_date - diagmdy),
    wt_endo_to_tx   = as.integer(first_tx_date - endoscopy_date),
    wt_dx_to_surg   = as.integer(surgery_date - diagmdy),
    wt_endo_to_surg = as.integer(surgery_date - endoscopy_date),
    wt_dx_to_sact   = as.integer(sact_date - diagmdy),
    wt_endo_to_sact = as.integer(sact_date - endoscopy_date),
    wt_sact_to_surg = as.integer(surgery_date - sact_date),
    wt_surg_to_sact = as.integer(sact_date - surgery_date),
    wt_dx_to_rt     = as.integer(rt_date - diagmdy),
    wt_endo_to_rt   = as.integer(rt_date - endoscopy_date),
    wt_rt_to_surg   = as.integer(surgery_date - rt_date),
    # --- survival ------------------------------------------------------------
    finmdy = as.Date(finmdy)
  )

# =============================================================================
# 4. Enforce the spec: light type coercion + column order
# =============================================================================
# Coerce only the columns where downstream code is type-sensitive.
og_cohort <- og_cohort %>%
  mutate(
    pseudo_patientid = as.character(pseudo_patientid),
    pseudo_tumourid  = as.character(pseudo_tumourid),
    diagmdy          = as.Date(diagmdy),
    finmdy           = as.Date(finmdy),
    ydiag            = as.integer(ydiag),
    tx_pathway       = factor(tx_pathway, levels = tx_pathway_levels) %>% as.character()
  )

# Keep spec columns that exist, in spec order (any_of tolerates absent optionals)
og_cohort_spec <- og_cohort %>% select(any_of(pre_cwt_spec$name))

# =============================================================================
# 5. Conformance check
# =============================================================================
check_conformance <- function(df, spec) {
  present  <- spec$name %in% names(df)
  missing_required <- spec$name[!present & spec$tier == "required"]
  missing_core     <- spec$name[!present & spec$tier == "core"]
  
  type_of <- function(x) {
    if (inherits(x, "Date"))      "Date"      else
      if (is.factor(x))             "factor"    else
        if (is.logical(x))            "logical"   else
          if (is.integer(x))            "integer"   else
            if (is.numeric(x))            "numeric"   else
              if (is.character(x))          "character" else class(x)[1]
  }
  type_rows <- spec %>% filter(name %in% names(df))
  type_obs  <- vapply(type_rows$name, function(n) type_of(df[[n]]), character(1))
  # integer/numeric treated as compatible
  compat <- function(exp, obs) (exp == obs) ||
    (exp %in% c("integer","numeric") && obs %in% c("integer","numeric"))
  mism <- type_rows %>%
    mutate(observed = type_obs) %>%
    filter(!mapply(compat, type, observed)) %>%
    select(name, expected = type, observed)
  
  bad_pathway <- setdiff(unique(na.omit(df$tx_pathway)), tx_pathway_levels)
  bad_stage   <- if ("stage_clean" %in% names(df))
    setdiff(unique(na.omit(df$stage_clean)), c("1","2","3")) else character(0)
  
  cat("== Conformance report ==\n")
  cat("Rows:", nrow(df), " Cols:", ncol(df), "\n")
  cat("Missing REQUIRED cols:", if (length(missing_required)) paste(missing_required, collapse=", ") else "none", "\n")
  cat("Missing core cols:    ", if (length(missing_core))     paste(missing_core, collapse=", ")     else "none", "\n")
  if (nrow(mism)) { cat("Type mismatches:\n"); print(mism) } else cat("Type mismatches:     none\n")
  cat("Unexpected tx_pathway:", if (length(bad_pathway)) paste(bad_pathway, collapse=", ") else "none", "\n")
  cat("Unexpected stage_clean:", if (length(bad_stage)) paste(bad_stage, collapse=", ") else "none", "\n")
  cat("Duplicate patient IDs:", sum(duplicated(df$pseudo_patientid)), "\n")
  invisible(list(missing_required = missing_required, type_mismatches = mism))
}

check_conformance(og_cohort_spec, pre_cwt_spec)

# =============================================================================
# 6. Save canonical pre-CWT cohort + spec objects
# =============================================================================
saveRDS(og_cohort_spec, file.path(syn_dir, "og_cohort_precwt_spec_2015_2022.rds"))

cat("\nSaved og_cohort_precwt_spec_2015_2022.rds (", nrow(og_cohort_spec),
    "patients, ", ncol(og_cohort_spec), "cols )\n")

}  # end anchor-guarded rebuild

