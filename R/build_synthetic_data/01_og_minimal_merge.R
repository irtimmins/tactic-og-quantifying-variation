# =============================================================================
# OG cancer - pathway derivation + CWT merge (shared engine)
# -----------------------------------------------------------------------------
# The two transparent coding stages the full pipeline performs, written once and
# shared between synthetic data and the real condensed cohort:
#
#   Stage 1  og_derive_pathway()  takes a RAW cohort - treatment dates, curative
#            descriptors, chemo provenance, and per-modality provider codes - and
#            derives the treatment flags, the sequencing flags, tx_pathway,
#            first_tx_date, and tx_trust. The pathway is a function of the flags
#            and dates alone; nothing is pre-supplied.
#
#   Stage 2  og_cwt_merge()  takes the derived cohort (Table A) plus the raw CWT
#            records (Table B) and attaches the decision-to-treat (DTT) node,
#            the waiting-time family, dtt_valid, and the audit categories.
#
# Table B holds one row per recorded CWT treatment event, dates as "dd/mm/yyyy"
# character exactly as the partitioned dataset stores them.
#
# Splitting the work this way means every coding stage is inspectable: the raw
# cohort shows the inputs, og_derive_pathway() shows how the pathway and trust
# are built from them, and og_cwt_merge() shows the linkage. The same two
# functions run on synthetic data and on the real cohort condensed to the raw
# column set (see condense_icon_to_minimal / condense_icon_to_raw below).
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# Constants - kept identical to the full pipeline so behaviour matches
# -----------------------------------------------------------------------------
og_merge_const <- list(
  tx_window_days   = 270L,
  dtt_min_offset   = -30L,
  treat_tol_days   = 14L,
  surg_switch_date = as.Date("2020-10-01"),
  surg_01_rule     = "date_split"   # date_split | always | never
)

tx_pathway_levels <- c(
  "EMR/ESD only", "EMR/ESD then surgery",
  "Surgery + neoadjuvant chemoRT", "Surgery + neoadjuvant chemo",
  "Surgery + neoadjuvant RT", "Surgery + adjuvant chemo",
  "Surgery only", "Surgery + other",
  "Definitive chemoRT", "Curative RT only",
  "Palliative chemo + RT", "SACT only", "Palliative RT only",
  "No treatment recorded"
)

# modality code -> broad group (CWT data dictionary; matches the full merge)
og_modality_group <- tribble(
  ~modality, ~mod_group,
  "01", "surgery", "23", "surgery", "24", "surgery",
  "02", "chemo",   "14", "chemo",   "15", "chemo",
  "03", "hormone",
  "04", "chemort",
  "05", "radiotherapy", "06", "radiotherapy", "13", "radiotherapy",
  "07", "palliative",   "08", "palliative",   "09", "palliative",
  "97", "other", "98", "declined"
)

# treatment groups that are plausible clock-stops for each pathway
og_pathway_groups <- list(
  "EMR/ESD only"                  = c("surgery", "other"),
  "EMR/ESD then surgery"          = c("surgery"),
  "Surgery + neoadjuvant chemoRT" = c("surgery", "chemort", "chemo", "radiotherapy"),
  "Surgery + neoadjuvant chemo"   = c("surgery", "chemo"),
  "Surgery + neoadjuvant RT"      = c("surgery", "radiotherapy", "chemort"),
  "Surgery + adjuvant chemo"      = c("surgery", "chemo"),
  "Surgery only"                  = c("surgery"),
  "Surgery + other"               = c("surgery", "other"),
  "Definitive chemoRT"            = c("chemort", "chemo", "radiotherapy"),
  "Curative RT only"              = c("radiotherapy", "chemort"),
  "Palliative chemo + RT"         = c("chemo", "radiotherapy", "chemort", "palliative"),
  "SACT only"                     = c("chemo", "hormone", "palliative"),
  "Palliative RT only"            = c("radiotherapy", "palliative"),
  "No treatment recorded"         = c("palliative", "other")
)

# the single defining clock-stop modality per pathway (the tie-break primary)
og_pathway_primary <- c(
  "EMR/ESD only" = "surgery", "EMR/ESD then surgery" = "surgery",
  "Surgery + neoadjuvant chemoRT" = "chemort",
  "Surgery + neoadjuvant chemo" = "chemo",
  "Surgery + neoadjuvant RT" = "radiotherapy",
  "Surgery + adjuvant chemo" = "surgery",
  "Surgery only" = "surgery", "Surgery + other" = "surgery",
  "Definitive chemoRT" = "chemort", "Curative RT only" = "radiotherapy",
  "Palliative chemo + RT" = "chemo", "SACT only" = "chemo",
  "Palliative RT only" = "radiotherapy", "No treatment recorded" = "palliative"
)

# =============================================================================
# Stage 1: derive the treatment pathway and the treatment trust from the raw
#          treatment flags, dates and provider codes
# -----------------------------------------------------------------------------
# This is the transparent classification step. It takes a "raw" cohort - one
# carrying the individual treatment dates, the curative descriptors, the chemo
# provenance, and the per-modality provider codes - and builds, in order:
#   - the treatment-presence flags (had_surgery, had_sact, ...)
#   - the sequencing flags (sact_before_surgery, concurrent_chemo_rt, ...)
#   - tx_pathway, via the case_when ladder
#   - first_tx_date, the clock-stop date for the pathway
#   - tx_trust, the provider of the clock-stop treatment
# Mirrors OG_cancer_prepare_data2_sact_rtds.R so the same derivation runs on the
# real condensed cohort as on synthetic data. Nothing here is pre-supplied: the
# pathway is a function of the flags and dates alone.
#
# The raw cohort must carry (see og_raw_cols): the registry descriptors, the
# treatment dates (endoscopy/emresd/surgery/sact/rt), the curative descriptors
# (curative_surgery, rt_curative), chemo_source (+ hes_chemo_date if present),
# and the provider codes surgery_provider (PROCODE3-style) and rt_provider
# (ORGCODEPROVIDER-style).
# -----------------------------------------------------------------------------
og_derive_pathway <- function(raw, const = og_merge_const) {
  
  # hes_chemo_date is optional; if absent treat it as missing so the chemo-RT
  # concurrency guard simply falls back to "SACT chemo always counts"
  if (!"hes_chemo_date" %in% names(raw)) raw$hes_chemo_date <- as.Date(NA)
  if (!"chemo_source"   %in% names(raw)) raw$chemo_source   <- NA_character_
  
  raw %>%
    mutate(
      # --- treatment-presence flags ----------------------------------------
      had_emresd           = !is.na(emresd_date),
      had_surgery          = !is.na(surgery_date),
      had_curative_surgery = !is.na(surgery_date) & curative_surgery == TRUE,
      had_sact             = !is.na(sact_date),
      had_rt               = !is.na(rt_date),
      had_curative_rt      = !is.na(rt_date) & rt_curative == TRUE,
      had_palliative_rt    = !is.na(rt_date) & rt_curative == FALSE,
      
      # chemo eligible to define a non-surgical definitive-chemoRT pathway:
      # SACT chemo always counts; HES-only chemo only when within 28 days of RT
      had_chemo_for_chemort = had_sact &
        ( coalesce(chemo_source, "sact") != "hes" |
            ( !is.na(hes_chemo_date) & !is.na(rt_date) &
                abs(as.integer(hes_chemo_date - rt_date)) <= 28 ) ),
      
      # --- sequencing flags ------------------------------------------------
      sact_before_surgery = had_sact & had_surgery & sact_date < surgery_date,
      sact_after_surgery  = had_sact & had_surgery & sact_date > surgery_date,
      rt_before_surgery   = had_rt   & had_surgery & rt_date   < surgery_date,
      rt_after_surgery    = had_rt   & had_surgery & rt_date   > surgery_date,
      concurrent_chemo_rt = had_sact & had_curative_rt &
        abs(as.integer(sact_date - rt_date)) <= 14,
      
      received_curative_tx = had_emresd | had_curative_surgery | had_curative_rt,
      
      # --- tx_pathway, from the flags only ---------------------------------
      tx_pathway = case_when(
        had_emresd & !had_surgery & !had_sact & !concurrent_chemo_rt ~ "EMR/ESD only",
        had_emresd & had_surgery                                     ~ "EMR/ESD then surgery",
        had_surgery & sact_before_surgery & rt_before_surgery        ~ "Surgery + neoadjuvant chemoRT",
        had_surgery & sact_before_surgery & !rt_before_surgery       ~ "Surgery + neoadjuvant chemo",
        had_surgery & rt_before_surgery & !sact_before_surgery       ~ "Surgery + neoadjuvant RT",
        had_surgery & sact_after_surgery & !sact_before_surgery      ~ "Surgery + adjuvant chemo",
        had_surgery & !had_sact & !concurrent_chemo_rt               ~ "Surgery only",
        had_surgery                                                  ~ "Surgery + other",
        !had_surgery & had_curative_rt & had_chemo_for_chemort       ~ "Definitive chemoRT",
        !had_surgery & had_curative_rt & !had_chemo_for_chemort      ~ "Curative RT only",
        !had_surgery & had_palliative_rt & had_sact                  ~ "Palliative chemo + RT",
        !had_surgery & had_sact & !had_curative_rt                   ~ "SACT only",
        !had_surgery & had_palliative_rt & !had_sact                 ~ "Palliative RT only",
        TRUE                                                         ~ "No treatment recorded"
      ),
      
      # --- first_tx_date: the clock-stop date for the pathway --------------
      first_tx_date = case_when(
        tx_pathway %in% c("EMR/ESD only", "EMR/ESD then surgery") ~ emresd_date,
        tx_pathway == "Surgery + neoadjuvant chemoRT" ~ pmin(sact_date, rt_date, na.rm = TRUE),
        tx_pathway == "Surgery + neoadjuvant RT"      ~ rt_date,
        tx_pathway == "Surgery + neoadjuvant chemo"   ~ sact_date,
        tx_pathway %in% c("Surgery + adjuvant chemo",
                          "Surgery only", "Surgery + other") ~ surgery_date,
        tx_pathway == "Definitive chemoRT" ~ pmin(sact_date, rt_date, na.rm = TRUE),
        tx_pathway == "Curative RT only"   ~ rt_date,
        TRUE                               ~ as.Date(NA)
      ),
      
      # --- tx_trust: provider of the clock-stop treatment ------------------
      # surgical/EMR pathways take the trust from the HES surgery provider; the
      # RT-anchored pathways take it from the RT provider. SACT's own provider is
      # never the trust source - even neoadjuvant chemo takes surgery's, because
      # the curative act is the surgery.
      tx_trust = case_when(
        tx_pathway %in% c("EMR/ESD only", "EMR/ESD then surgery",
                          "Surgery + neoadjuvant chemo",
                          "Surgery + adjuvant chemo",
                          "Surgery only", "Surgery + other") ~ substr(surgery_provider, 1, 3),
        tx_pathway %in% c("Surgery + neoadjuvant chemoRT",
                          "Surgery + neoadjuvant RT",
                          "Definitive chemoRT",
                          "Curative RT only")                ~ substr(rt_provider, 1, 3),
        TRUE                                                 ~ NA_character_
      )
    )
}

# the raw cohort contract: what og_derive_pathway() consumes
og_raw_cols <- c(
  "pseudo_patientid", "diagmdy", "ydiag",
  "sex", "agediag", "ethnicity_group_broad",
  "NHSE_reversed_imd_quintile_lsoas",
  "tumour_site_grp", "cancer_subtype", "stage_clean",
  "route_combined", "cci_group",
  "diag_trust", "diag_hosp",
  "endoscopy_date", "emresd_date", "surgery_date", "sact_date", "rt_date",
  "curative_surgery", "rt_curative", "chemo_source", "hes_chemo_date",
  "surgery_provider", "rt_provider", "sact_provider",
  "finmdy", "died"
)

# -----------------------------------------------------------------------------
# Stage 2: og_cwt_merge(A, cwt, const)
#   A    : Table A - the derived cohort from og_derive_pathway() (or any cohort
#          carrying tx_pathway + first_tx_date + the treatment/anchor dates)
#   cwt  : raw CWT records (Table B), character dd/mm/yyyy dates
#   const: merge constants (defaults to og_merge_const)
# Returns A with cwt_dtt_date, cwt_mdt_date, cwt_treat_date, cwt_modality, the
# waiting-time family, dtt_valid, and the audit categories attached.
# -----------------------------------------------------------------------------
og_cwt_merge <- function(A, cwt, const = og_merge_const) {
  
  # 1. parse CWT dates and assign each record a modality group, applying the
  #    surgery-01 date-split rule
  cwt_grouped <- cwt %>%
    mutate(
      pseudo_patientid = as.character(pseudo_patientid),
      cwt_dtt_date     = as.Date(treat_period_start, "%d/%m/%Y"),
      cwt_treat_date   = as.Date(treat_start,        "%d/%m/%Y"),
      cwt_mdt_date     = as.Date(mdt_date,           "%d/%m/%Y")
    ) %>%
    left_join(og_modality_group, by = "modality") %>%
    mutate(mod_group = case_when(
      modality == "01" & const$surg_01_rule == "never"      ~ NA_character_,
      modality == "01" & const$surg_01_rule == "date_split" &
        cwt_treat_date >= const$surg_switch_date            ~ NA_character_,
      modality %in% c("23","24") & const$surg_01_rule == "date_split" &
        cwt_treat_date <  const$surg_switch_date            ~ NA_character_,
      TRUE                                                  ~ mod_group
    )) %>%
    filter(!is.na(mod_group), mod_group != "declined", !is.na(cwt_dtt_date))
  
  pathway_group_long <- enframe(og_pathway_groups, name = "tx_pathway",
                                value = "ok_group") %>%
    unnest_longer(ok_group)
  
  pw <- A %>% select(pseudo_patientid, diagmdy, tx_pathway, first_tx_date)
  
  # 2. candidate rows: in-window, and (where any exist) pathway-consistent
  candidates <- cwt_grouped %>%
    inner_join(pw, by = "pseudo_patientid") %>%
    mutate(days_dx_to_dtt = as.integer(cwt_dtt_date - diagmdy)) %>%
    filter(days_dx_to_dtt >= const$dtt_min_offset,
           days_dx_to_dtt <= const$tx_window_days) %>%
    left_join(pathway_group_long %>% mutate(group_ok = TRUE),
              by = c("tx_pathway", "mod_group" = "ok_group")) %>%
    mutate(group_ok   = coalesce(group_ok, FALSE),
           is_primary = mod_group == unname(og_pathway_primary[tx_pathway])) %>%
    group_by(pseudo_patientid) %>%
    mutate(any_match = any(group_ok)) %>%
    filter(if (first(any_match)) group_ok else TRUE) %>%
    ungroup()
  
  # 3. anchor: primary modality wins, then earliest DTT
  anchor <- candidates %>%
    group_by(pseudo_patientid) %>%
    arrange(desc(is_primary), cwt_dtt_date, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(pseudo_patientid, cwt_dtt_date, cwt_mdt_date, cwt_treat_date,
              cwt_modality = modality)
  
  # 4. attach to A and derive the waiting-time family + validity
  out <- A %>%
    left_join(anchor, by = "pseudo_patientid") %>%
    mutate(
      # the six core intervals: each of diagnosis / endoscopy to each of the
      # decision-to-treat (DTT) and first treatment, plus the DTT->treatment link
      wt_endo_to_dx  = as.integer(diagmdy      - endoscopy_date),
      wt_dx_to_dtt   = as.integer(cwt_dtt_date  - diagmdy),
      wt_endo_to_dtt = as.integer(cwt_dtt_date  - endoscopy_date),
      wt_dx_to_tx    = as.integer(first_tx_date - diagmdy),
      wt_endo_to_tx  = as.integer(first_tx_date - endoscopy_date),
      wt_dtt_to_tx   = as.integer(first_tx_date - cwt_dtt_date),
      # per-modality component intervals (which treatment arm the clock stopped
      # on), matching the full pipeline; NA where that arm did not occur
      wt_dx_to_surg   = as.integer(surgery_date - diagmdy),
      wt_endo_to_surg = as.integer(surgery_date - endoscopy_date),
      wt_dx_to_sact   = as.integer(sact_date    - diagmdy),
      wt_endo_to_sact = as.integer(sact_date    - endoscopy_date),
      wt_dx_to_rt     = as.integer(rt_date      - diagmdy),
      wt_endo_to_rt   = as.integer(rt_date      - endoscopy_date),
      # treatment-sequencing gaps (neoadjuvant / adjuvant)
      wt_sact_to_surg = as.integer(surgery_date - sact_date),
      wt_surg_to_sact = as.integer(sact_date    - surgery_date),
      wt_rt_to_surg   = as.integer(surgery_date - rt_date),
      dtt_valid = !is.na(cwt_dtt_date) & wt_dx_to_dtt >= 0 &
        wt_dtt_to_tx >= -const$treat_tol_days,
      dtt_valid = if_else(tx_pathway %in% c("EMR/ESD only", "EMR/ESD then surgery"),
                          NA, dtt_valid)
    )
  
  # 5. audit categories (NOGCA Tables 3 & 4 groupings)
  out %>%
    mutate(
      tx_modality_audit = case_when(
        tx_pathway == "Surgery only"                              ~ "Surgery only",
        tx_pathway %in% c("Surgery + neoadjuvant chemo",
                          "Surgery + neoadjuvant chemoRT",
                          "Surgery + neoadjuvant RT",
                          "Surgery + adjuvant chemo",
                          "Surgery + other")                      ~ "Surgery plus SACT/RT",
        tx_pathway %in% c("EMR/ESD only", "EMR/ESD then surgery") ~ "EMR/ESD",
        tx_pathway == "Definitive chemoRT"                        ~ "Definitive chemoRT",
        tx_pathway == "Curative RT only"                          ~ "Curative RT only",
        tx_pathway %in% c("Palliative chemo + RT", "SACT only",
                          "Palliative RT only")                   ~ "Chemo/RT only (non-curative)",
        tx_pathway == "No treatment recorded"                     ~ "No treatment recorded",
        TRUE                                                      ~ NA_character_),
      tx_modality_audit = factor(tx_modality_audit, levels = c(
        "Surgery only", "Surgery plus SACT/RT", "EMR/ESD",
        "Definitive chemoRT", "Curative RT only",
        "Chemo/RT only (non-curative)", "No treatment recorded")),
      tx_intent_audit = case_when(
        tx_modality_audit %in% c("Surgery only", "Surgery plus SACT/RT",
                                 "EMR/ESD", "Definitive chemoRT",
                                 "Curative RT only")              ~ "Curative",
        tx_modality_audit == "Chemo/RT only (non-curative)"       ~ "Non-curative",
        tx_modality_audit == "No treatment recorded"              ~ "No treatment",
        TRUE                                                      ~ NA_character_),
      received_any_tx = tx_pathway != "No treatment recorded" &
        ( (!is.na(first_tx_date) & wt_dx_to_tx <= const$tx_window_days) |
            (!is.na(cwt_treat_date) &
               as.integer(cwt_treat_date - diagmdy) >= 0 &
               as.integer(cwt_treat_date - diagmdy) <= const$tx_window_days) ),
      received_curative_tx_audit = tx_intent_audit == "Curative" &
        !is.na(first_tx_date) & wt_dx_to_tx <= const$tx_window_days
    )
}

# -----------------------------------------------------------------------------
# The minimal Table A contract: the columns og_cwt_merge needs, plus the basic
# patient descriptors that make the cohort usable. Anything else is optional.
# -----------------------------------------------------------------------------
og_minimal_cols <- c(
  # identity & cohort
  "pseudo_patientid", "diagmdy", "ydiag",
  # basic patient / tumour descriptors
  "sex", "agediag", "ethnicity_group_broad",
  "NHSE_reversed_imd_quintile_lsoas",
  "tumour_site_grp", "cancer_subtype", "stage_clean",
  "route_combined", "cci_group",
  # diagnosing organisation
  "diag_trust", "diag_hosp",
  # treatment anchor dates (the merge needs endoscopy_date + first_tx_date;
  # the others let the pathway and intervals be re-derived/inspected)
  "endoscopy_date", "emresd_date", "surgery_date", "sact_date", "rt_date",
  "first_tx_date",
  # pathway, provenance, trust
  "tx_pathway", "chemo_source", "tx_trust",
  # survival
  "finmdy", "died"
)

# -----------------------------------------------------------------------------
# condense_icon_to_minimal(full)
#   Reduce the full ICON og_cohort_cwt (128 columns) to the minimal Table A so
#   the same minimal merge can be re-run on the real data. Keeps only the
#   minimal columns that are present; derives cci_group if absent.
# -----------------------------------------------------------------------------
condense_icon_to_minimal <- function(full) {
  out <- full
  if (!"cci_group" %in% names(out) && "rcs_ch_score" %in% names(out)) {
    out <- out %>% mutate(cci_group = cut(
      rcs_ch_score, breaks = c(-Inf, 0, 1, 2, Inf),
      labels = c("0", "1", "2", "3+"), right = TRUE) %>% as.character())
  }
  if (!"tumour_site_grp" %in% names(out) && "cancer_subtype" %in% names(out)) {
    out <- out %>% mutate(tumour_site_grp = if_else(
      grepl("^Gast", coalesce(cancer_subtype, "")), "gastric", "oesophageal"))
  }
  out %>% select(any_of(og_minimal_cols))
}

# -----------------------------------------------------------------------------
# condense_icon_to_raw(full)
#   Reduce the full ICON og_cohort to the RAW derivation input - the treatment
#   dates, curative descriptors, chemo provenance and provider codes - so
#   og_derive_pathway() can be re-run on the real data and checked against the
#   pipeline's own tx_pathway. Renames the real provider fields (PROCODE3,
#   ORGCODEPROVIDER, ORGANISATION_CODE_OF_PROVIDER) to the raw-contract names.
# -----------------------------------------------------------------------------
condense_icon_to_raw <- function(full) {
  out <- full
  if (!"cci_group" %in% names(out) && "rcs_ch_score" %in% names(out)) {
    out <- out %>% mutate(cci_group = cut(
      rcs_ch_score, breaks = c(-Inf, 0, 1, 2, Inf),
      labels = c("0", "1", "2", "3+"), right = TRUE) %>% as.character())
  }
  if (!"tumour_site_grp" %in% names(out) && "cancer_subtype" %in% names(out)) {
    out <- out %>% mutate(tumour_site_grp = if_else(
      grepl("^Gast", coalesce(cancer_subtype, "")), "gastric", "oesophageal"))
  }
  rename_if_present <- function(d, new, old)
    if (old %in% names(d) && !new %in% names(d)) rename(d, !!new := !!sym(old)) else d
  out <- out %>%
    rename_if_present("surgery_provider", "PROCODE3") %>%
    rename_if_present("rt_provider",      "ORGCODEPROVIDER") %>%
    rename_if_present("sact_provider",    "ORGANISATION_CODE_OF_PROVIDER")
  out %>% select(any_of(og_raw_cols))
}

cat("Loaded og_derive_pathway() [stage 1], og_cwt_merge() [stage 2],",
    "condense_icon_to_raw()/condense_icon_to_minimal(), and the OG constants.\n",
    "Raw cohort needs", length(og_raw_cols), "columns; derivation builds",
    "tx_pathway, first_tx_date and tx_trust from the flags and dates alone.\n")
