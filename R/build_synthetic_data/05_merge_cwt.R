# =============================================================================
# 08  Merge CWT
# -----------------------------------------------------------------------------
# The CWT extract holds many rows per patient (one per recorded treatment event).
# This step picks the single record matching the patient's treatment pathway,
# takes the earliest valid Decision To Treat (DTT) among the matching records,
# and derives the waiting-time family and the audit categories. The result is the
# final analysis cohort.
#
# Surgery modality note: CWT code 01 (Surgery) was retired on 2020-10-01 and
# replaced by 23/24. surg_01_rule (in 01_define_parameters.R) controls how 01 is
# treated across the switch; date_split is the data-dictionary-faithful default.
#
# Reads : Data/ICON/og_cohort_2015_2022.rds, the partitioned CWT dataset
# Writes: Data/ICON/og_cohort_cwt_2015_2022.rds   (the final cohort)
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

# -----------------------------------------------------------------------------
# CWT modality lookups (specific to this merge)
# -----------------------------------------------------------------------------
# each CWT modality maps to a broad treatment group
modality_group <- tribble(
  ~modality, ~mod_group,
  "01", "surgery", "23", "surgery", "24", "surgery",
  "02", "chemo", "14", "chemo", "15", "chemo", "03", "hormone",
  "04", "chemort",
  "05", "radiotherapy", "06", "radiotherapy", "13", "radiotherapy",
  "07", "palliative", "08", "palliative", "09", "palliative",
  "97", "other", "98", "declined")

# treatment groups that are plausible clock-stops for each pathway
pathway_groups <- list(
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
  "No treatment recorded"         = c("palliative", "other"))

# the single defining clock-stop modality per pathway: it wins over a merely
# pathway-consistent row even with a later DTT, so a neoadjuvant patient anchors
# on their chemo/RT, not a surgical row carrying an earlier decision-to-treat.
pathway_primary <- c(
  "EMR/ESD only" = "surgery", "EMR/ESD then surgery" = "surgery",
  "Surgery + neoadjuvant chemoRT" = "chemort", "Surgery + neoadjuvant chemo" = "chemo",
  "Surgery + neoadjuvant RT" = "radiotherapy", "Surgery + adjuvant chemo" = "surgery",
  "Surgery only" = "surgery", "Surgery + other" = "surgery",
  "Definitive chemoRT" = "chemort", "Curative RT only" = "radiotherapy",
  "Palliative chemo + RT" = "chemo", "SACT only" = "chemo",
  "Palliative RT only" = "radiotherapy", "No treatment recorded" = "palliative")

# -----------------------------------------------------------------------------
# Inputs and CWT read
# -----------------------------------------------------------------------------
og_cohort    <- readRDS(f_cohort)
ncras_og_ids <- unique(as.character(og_cohort$pseudo_patientid))

# read_cwt() returns the raw CWT rows for the OG sites. A test harness can define
# read_cwt before sourcing this script to supply a fixture; by default it reads
# the partitioned dataset.
if (!exists("read_cwt"))
  read_cwt <- function() open_dataset(path_cwt_partition) %>%
  filter(site_icd10 %in% og_icd10) %>% collect()

cwt_og <- read_cwt() %>%
  mutate(pseudo_patientid  = as.character(pseudo_patientid),
         cwt_dtt_date      = as.Date(treat_period_start, "%d/%m/%Y"),
         cwt_treat_date    = as.Date(treat_start,        "%d/%m/%Y"),
         cwt_referral_date = as.Date(crtp_date,          "%d/%m/%Y"),
         cwt_first_seen    = as.Date(date_first_seen,    "%d/%m/%Y"),
         cwt_mdt_date      = as.Date(mdt_date,           "%d/%m/%Y")) %>%
  filter(pseudo_patientid %in% ncras_og_ids)

# -----------------------------------------------------------------------------
# Assign each CWT row a treatment group, applying the surgery (01/23/24) rule.
# transition_window: inside [start, end] all three surgery codes count; before
# the window only 01 counts; after it only 23/24. This rescues the 2020 overlap
# year, where a single hard switch date would drop early 23/24 and late 01.
# date_split is the older single-cliff behaviour, kept as a fallback.
# -----------------------------------------------------------------------------
cwt_grouped <- cwt_og %>%
  left_join(modality_group, by = "modality") %>%
  mutate(mod_group = case_when(
    # transition window
    modality == "01" & surg_01_rule == "transition_window" &
      cwt_treat_date > surg_transition_end                                 ~ NA_character_,
    modality %in% c("23","24") & surg_01_rule == "transition_window" &
      cwt_treat_date < surg_transition_start                               ~ NA_character_,
    # single-cliff fallback
    modality == "01" & surg_01_rule == "date_split" &
      cwt_treat_date >= surg_switch_date                                   ~ NA_character_,
    modality %in% c("23","24") & surg_01_rule == "date_split" &
      cwt_treat_date < surg_switch_date                                    ~ NA_character_,
    # always drop 01 if requested
    modality == "01" & surg_01_rule == "never"                            ~ NA_character_,
    TRUE                                                                  ~ mod_group)) %>%
  filter(!is.na(mod_group), mod_group != "declined", !is.na(cwt_dtt_date))

# -----------------------------------------------------------------------------
# Candidate rows (in-window, pathway-consistent) and the anchor
# -----------------------------------------------------------------------------
pathway_group_long <- enframe(pathway_groups, name = "tx_pathway", value = "ok_group") %>%
  unnest_longer(ok_group)
pw <- og_cohort %>% select(pseudo_patientid, diagmdy, tx_pathway, first_tx_date)

cwt_candidates <- cwt_grouped %>%
  inner_join(pw, by = "pseudo_patientid") %>%
  mutate(days_dx_to_dtt = as.integer(cwt_dtt_date - diagmdy)) %>%
  filter(days_dx_to_dtt >= dtt_min_offset, days_dx_to_dtt <= cwt_window_days) %>%
  left_join(pathway_group_long %>% mutate(group_ok = TRUE),
            by = c("tx_pathway", "mod_group" = "ok_group")) %>%
  mutate(group_ok   = coalesce(group_ok, FALSE),
         is_primary = mod_group == unname(pathway_primary[tx_pathway])) %>%
  group_by(pseudo_patientid) %>%
  mutate(any_match = any(group_ok)) %>%
  filter(if (first(any_match)) group_ok else TRUE) %>%
  ungroup()

# anchor: primary modality wins, then earliest DTT
cwt_anchor <- cwt_candidates %>%
  group_by(pseudo_patientid) %>%
  arrange(desc(is_primary), cwt_dtt_date, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(pseudo_patientid, cwt_dtt_date, cwt_treat_date, cwt_mdt_date,
            modality, mod_group, matched_pathway = any_match, days_dx_to_dtt)

cat("CWT anchored:", nrow(cwt_anchor), "patients (",
    round(100 * mean(cwt_anchor$matched_pathway), 1), "% pathway-consistent )\n")

# -----------------------------------------------------------------------------
# Attach the DTT node, derive the waiting-time family and validity
# -----------------------------------------------------------------------------
og_cohort <- og_cohort %>%
  left_join(cwt_anchor %>% select(pseudo_patientid, cwt_dtt_date, cwt_mdt_date,
                                  cwt_treat_date, cwt_modality = modality),
            by = "pseudo_patientid") %>%
  mutate(
    wt_endo_to_dtt = as.integer(cwt_dtt_date  - endoscopy_date),
    wt_dtt_to_tx   = as.integer(first_tx_date - cwt_dtt_date),
    wt_dx_to_dtt   = as.integer(cwt_dtt_date  - diagmdy),
    # valid if DTT is on/after diagnosis and treatment does not precede it by
    # more than the tolerance; NA for EMR/ESD where the DTT is less meaningful
    dtt_valid = !is.na(cwt_dtt_date) & wt_dx_to_dtt >= 0 &
      wt_dtt_to_tx >= -treat_tol_days,
    dtt_valid = if_else(tx_pathway %in% c("EMR/ESD only", "EMR/ESD then surgery"),
                        NA, dtt_valid))

# -----------------------------------------------------------------------------
# Audit categories (NOGCA Tables 3 and 4), derived from tx_pathway
# -----------------------------------------------------------------------------
og_cohort <- og_cohort %>%
  mutate(
    tx_modality_audit = case_when(
      tx_pathway == "Surgery only"                                ~ "Surgery only",
      tx_pathway %in% c("Surgery + neoadjuvant chemo",
                        "Surgery + neoadjuvant chemoRT",
                        "Surgery + neoadjuvant RT",
                        "Surgery + adjuvant chemo",
                        "Surgery + other")                        ~ "Surgery plus SACT/RT",
      tx_pathway %in% c("EMR/ESD only", "EMR/ESD then surgery")   ~ "EMR/ESD",
      tx_pathway == "Definitive chemoRT"                          ~ "Definitive chemoRT",
      tx_pathway == "Curative RT only"                            ~ "Curative RT only",
      tx_pathway %in% c("Palliative chemo + RT", "SACT only",
                        "Palliative RT only")                     ~ "Chemo/RT only (non-curative)",
      tx_pathway == "No treatment recorded"                       ~ "No treatment recorded",
      TRUE                                                        ~ NA_character_),
    tx_modality_audit = factor(tx_modality_audit, levels = c(
      "Surgery only", "Surgery plus SACT/RT", "EMR/ESD",
      "Definitive chemoRT", "Curative RT only",
      "Chemo/RT only (non-curative)", "No treatment recorded")),
    
    tx_intent_audit = case_when(
      tx_modality_audit %in% c("Surgery only", "Surgery plus SACT/RT", "EMR/ESD",
                               "Definitive chemoRT", "Curative RT only") ~ "Curative",
      tx_modality_audit == "Chemo/RT only (non-curative)"               ~ "Non-curative",
      tx_modality_audit == "No treatment recorded"                      ~ "No treatment",
      TRUE                                                              ~ NA_character_),
    
    # received any treatment within nine months. Counts the CWT treatment date as
    # well as first_tx_date, since first_tx_date is NA for palliative patients who
    # nonetheless received an in-window treatment.
    received_any_tx = tx_pathway != "No treatment recorded" &
      ( (!is.na(first_tx_date) & wt_dx_to_tx <= tx_window_days) |
          (!is.na(cwt_treat_date) &
             as.integer(cwt_treat_date - diagmdy) >= 0 &
             as.integer(cwt_treat_date - diagmdy) <= tx_window_days) ),
    received_curative_tx_audit = tx_intent_audit == "Curative" &
      !is.na(first_tx_date) & wt_dx_to_tx <= tx_window_days)

saveRDS(og_cohort, f_cohort_cwt)
cat("Saved final cohort:", f_cohort_cwt, "(", nrow(og_cohort), "patients ).\n")
cat("08 complete. Next: 09_validate_build.R\n")