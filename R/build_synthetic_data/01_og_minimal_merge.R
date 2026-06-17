# =============================================================================
# OG cancer - minimal CWT merge (shared engine)
# -----------------------------------------------------------------------------
# A self-contained, minimal version of the CWT merge that the full pipeline
# (OG_cancer_prepare_data4_cwt_merge.R) performs. It operates on:
#
#   Table A  a minimal registry + treatment cohort (one row per patient), with
#            the treatment anchor dates, the derived tx_pathway, and first_tx_date
#   Table B  raw CWT records (one row per recorded treatment event, dates as
#            "dd/mm/yyyy" character, exactly as the partitioned dataset stores them)
#
# and returns Table A with the CWT decision-to-treat (DTT) node attached and the
# waiting-time + audit fields derived. The same function runs on:
#   - synthetic Table A / Table B from the generator, and
#   - the real ICON cohort once condensed to the minimal column set (see
#     condense_icon_to_minimal() below),
# so the merge logic is written once and shared.
#
# This is deliberately the *minimal* merge: modality-to-pathway matching, the
# neoadjuvant primary-modality tie-break, the in-window DTT filter, and the
# audit categorisation. It does not reproduce the descriptive/diagnostic blocks
# of the full script.
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# Constants - kept identical to the full merge so behaviour matches
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

# -----------------------------------------------------------------------------
# og_cwt_merge(A, cwt, const)
#   A    : minimal Table A (see og_minimal_cols below for the contract)
#   cwt  : raw CWT records (Table B), character dd/mm/yyyy dates
#   const: merge constants (defaults to og_merge_const)
# Returns A with cwt_dtt_date, cwt_mdt_date, cwt_treat_date, cwt_modality, the
# wt_*_dtt intervals, dtt_valid, and the audit categories attached.
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

cat("Loaded og_cwt_merge(), condense_icon_to_minimal(), and the OG merge",
    "constants/lookups.\n",
    "Minimal Table A needs:", length(og_minimal_cols), "columns;",
    "the merge itself uses diagmdy, endoscopy_date, first_tx_date, tx_pathway.\n")