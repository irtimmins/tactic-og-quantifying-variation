
# =============================================================================
# OG cancer - merge Cancer Waiting Times (CWT) and derive the DTT intervals
# -----------------------------------------------------------------------------
# The CWT extract holds many rows per patient (one per recorded treatment event,
# first and subsequent). This script picks the single record that matches the
# patient's derived treatment pathway, takes the earliest valid Decision To Treat
# (DTT) among the matching records, and derives the waiting-time components:
#   diagnosis/endoscopy -> DTT   (staging / work-up)
#   DTT -> treatment             (scheduling)
#
# Modality handling (NHS CWT data dictionary, National Cancer WT MDS):
#   The modality of each CWT row is matched to the pathway so that, for example,
#   a surgical patient is anchored on the surgical CWT record rather than an
#   earlier palliative or chemo row. Key point for a 2015-2022 cohort:
#     - 01 Surgery was RETIRED on 01 October 2020 and replaced by
#       23 Surgery (excluding enabling treatment) and 24 Surgery (enabling).
#     - So surgery appears as 01 before the switch and 23/24 after it; both are
#       surgery, not endoscopy. The switch date is a parameter below so the
#       effect of including / excluding 01 at different points can be explored.
#   Other groups: chemo 02/14/15 (+03 hormone), radiotherapy 05/06/13,
#   chemoradiotherapy 04, palliative/active monitoring 07/08/09, 98 declined,
#   97 other.
# =============================================================================

library(tidyverse)
library(arrow)
library(haven)
library(lubridate)

base_dir <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/"
cwt_path <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/CWT/11_CWT_data_partitioned"

# -----------------------------------------------------------------------------
# 0. Parameters (kept explicit so the matching rules are auditable)
# -----------------------------------------------------------------------------
tx_window_days   <- 270L              # DTT must sit within this many days of dx
dtt_min_offset   <- -30L              # earliest DTT relative to diagnosis (days)
treat_tol_days   <- 14L               # tolerance: treatment may precede DTT by up to this

# surgery modality handling -------------------------------------------------
# 01 (old surgery) retired 01 Oct 2020; 23/24 are the replacements.
surg_switch_date <- as.Date("2020-10-01")
# how to treat code 01: "date_split" counts 01 as surgery only before the
# switch (and 23/24 only on/after it); "always" counts 01 as surgery throughout;
# "never" ignores 01. date_split is the data-dictionary-faithful default.
surg_01_rule     <- "date_split"      # one of: date_split | always | never

og_icd10 <- c("C150","C151","C152","C153","C154","C155","C158","C159","C15",
              "C160","C161","C162","C163","C164","C165","C166","C168","C169","C16")

# -----------------------------------------------------------------------------
# 1. Modality -> pathway-group lookup
# -----------------------------------------------------------------------------
# Each CWT modality maps to a broad treatment group. The patient's tx_pathway
# (from the HES/SACT/RTDS build) maps to the same groups, and we keep CWT rows
# whose group is consistent with the pathway.
modality_group <- tribble(
  ~modality, ~mod_group,
  "01", "surgery",      # retired 2020-10; handled by surg_01_rule
  "23", "surgery",
  "24", "surgery",
  "02", "chemo",
  "14", "chemo",
  "15", "chemo",
  "03", "hormone",
  "04", "chemort",
  "05", "radiotherapy",
  "06", "radiotherapy",
  "13", "radiotherapy",
  "07", "palliative",
  "08", "palliative",
  "09", "palliative",
  "97", "other",
  "98", "declined"
)

# the treatment groups that are plausible clock-stops for each pathway
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
  "No treatment recorded"         = c("palliative", "other")
)

# the single defining clock-stop modality for each pathway. When more than one
# pathway-consistent CWT row is available, the row in this group wins over a
# merely-consistent row, even if the other has an earlier DTT. This anchors a
# neoadjuvant patient on their chemo/RT rather than a surgical row that happens
# to carry an earlier decision-to-treat date (the first definitive treatment for
# a neoadjuvant pathway is the neoadjuvant treatment, per CWT guidance 3.9.1).
pathway_primary <- c(
  "EMR/ESD only"                  = "surgery",
  "EMR/ESD then surgery"          = "surgery",
  "Surgery + neoadjuvant chemoRT" = "chemort",
  "Surgery + neoadjuvant chemo"   = "chemo",
  "Surgery + neoadjuvant RT"      = "radiotherapy",
  "Surgery + adjuvant chemo"      = "surgery",
  "Surgery only"                  = "surgery",
  "Surgery + other"               = "surgery",
  "Definitive chemoRT"            = "chemort",
  "Curative RT only"              = "radiotherapy",
  "Palliative chemo + RT"         = "chemo",
  "SACT only"                     = "chemo",
  "Palliative RT only"            = "radiotherapy",
  "No treatment recorded"         = "palliative"
)

# -----------------------------------------------------------------------------
# 2. Inputs
# -----------------------------------------------------------------------------
ncras_og     <- readRDS(paste0(base_dir, "ncras_og_2015_2022.rds"))
ncras_og_ids <- unique(as.character(ncras_og$pseudo_patientid))

og_cohort <- readRDS(paste0(base_dir, "og_cohort_2015_2022.rds"))

# RCS Charlson comorbidity lookup from script 1d
og_cci    <- readRDS(paste0(base_dir, "og_cci_2015_2022.rds"))
og_cohort <- og_cohort %>% left_join(og_cci, by = "pseudo_patientid")

# -----------------------------------------------------------------------------
# 3. Read the OG CWT records and parse dates
# -----------------------------------------------------------------------------
cwt_og <- open_dataset(cwt_path) %>%
  filter(site_icd10 %in% og_icd10) %>%
  collect() %>%
  mutate(
    pseudo_patientid  = as.character(pseudo_patientid),
    cwt_dtt_date      = as.Date(treat_period_start, format = "%d/%m/%Y"),
    cwt_treat_date    = as.Date(treat_start,        format = "%d/%m/%Y"),
    cwt_referral_date = as.Date(crtp_date,          format = "%d/%m/%Y"),
    cwt_first_seen    = as.Date(date_first_seen,    format = "%d/%m/%Y"),
    cwt_mdt_date      = as.Date(mdt_date,           format = "%d/%m/%Y")
  ) %>%
  filter(pseudo_patientid %in% ncras_og_ids)

cat("CWT OG rows:", nrow(cwt_og),
    "| patients:", n_distinct(cwt_og$pseudo_patientid), "\n")
cat("Modality distribution:\n")
cwt_og %>% count(modality, sort = TRUE) %>% print(n = 30)

# -----------------------------------------------------------------------------
# 3b. Describe the raw CWT extract: records per patient and field completeness
# -----------------------------------------------------------------------------
cat("\nCWT records per patient (raw extract, all modalities):\n")
recs_per <- cwt_og %>% count(pseudo_patientid, name = "n_records")
recs_per %>% count(n_records, name = "n_patients") %>%
  mutate(pct = round(100 * n_patients / sum(n_patients), 1)) %>%
  arrange(n_records) %>% print(n = 30)
cat("  mean records/patient:", round(mean(recs_per$n_records), 2),
    "| median:", median(recs_per$n_records),
    "| max:", max(recs_per$n_records), "\n")
cat("  patients with >1 record:", sum(recs_per$n_records > 1),
    paste0("(", round(100 * mean(recs_per$n_records > 1), 1), "%)"), "\n")

cat("\nDistinct DTT dates per patient (raw, any modality, parseable DTT):\n")
dtt_per <- cwt_og %>% filter(!is.na(cwt_dtt_date)) %>%
  group_by(pseudo_patientid) %>%
  summarise(n_dtt = n_distinct(cwt_dtt_date), .groups = "drop")
dtt_per %>% count(n_dtt, name = "n_patients") %>%
  mutate(pct = round(100 * n_patients / sum(n_patients), 1)) %>%
  arrange(n_dtt) %>% print(n = 30)
cat("  patients with >1 distinct DTT:", sum(dtt_per$n_dtt > 1),
    paste0("(", round(100 * mean(dtt_per$n_dtt > 1), 1), "%)"), "\n")

# -----------------------------------------------------------------------------
# 4. Assign each CWT row a treatment group, applying the surgery (01) rule
# -----------------------------------------------------------------------------
cwt_grouped <- cwt_og %>%
  left_join(modality_group, by = "modality") %>%
  mutate(
    mod_group = case_when(
      # code 01 only counts as surgery according to the chosen rule
      modality == "01" & surg_01_rule == "never"                          ~ NA_character_,
      modality == "01" & surg_01_rule == "date_split" &
        cwt_treat_date >= surg_switch_date                                ~ NA_character_,
      # 23/24 are post-switch surgery; under date_split ignore them before
      modality %in% c("23","24") & surg_01_rule == "date_split" &
        cwt_treat_date < surg_switch_date                                 ~ NA_character_,
      TRUE                                                                ~ mod_group
    )
  ) %>%
  # drop declined / unusable rows and rows with no usable date
  filter(!is.na(mod_group), mod_group != "declined", !is.na(cwt_dtt_date))

# -----------------------------------------------------------------------------
# 5. Match modality group to the patient's pathway, keep the earliest valid DTT
# -----------------------------------------------------------------------------
pathway_group_long <- enframe(pathway_groups, name = "tx_pathway",
                              value = "ok_group") %>%
  unnest_longer(ok_group)

pw <- og_cohort %>% select(pseudo_patientid, diagmdy, tx_pathway, first_tx_date)

# candidate rows: in-window, and (where any exist) pathway-consistent. These are
# the legitimate decision-to-treat records the anchor is chosen from.
cwt_candidates <- cwt_grouped %>%
  inner_join(pw, by = "pseudo_patientid") %>%
  mutate(days_dx_to_dtt = as.integer(cwt_dtt_date - diagmdy)) %>%
  filter(days_dx_to_dtt >= dtt_min_offset, days_dx_to_dtt <= tx_window_days) %>%
  left_join(pathway_group_long %>% mutate(group_ok = TRUE),
            by = c("tx_pathway", "mod_group" = "ok_group")) %>%
  mutate(group_ok   = coalesce(group_ok, FALSE),
         is_primary = mod_group == unname(pathway_primary[tx_pathway])) %>%
  group_by(pseudo_patientid) %>%
  mutate(any_match = any(group_ok)) %>%
  filter(if (first(any_match)) group_ok else TRUE) %>%
  ungroup()

# the anchor: primary modality wins, then earliest DTT
cwt_anchor <- cwt_candidates %>%
  group_by(pseudo_patientid) %>%
  arrange(desc(is_primary), cwt_dtt_date, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(pseudo_patientid, cwt_dtt_date, cwt_treat_date, cwt_mdt_date,
            modality, mod_group, matched_pathway = any_match, days_dx_to_dtt)

cat("\nCWT anchor patients:", nrow(cwt_anchor), "\n")
cat("  pathway-consistent modality:",
    sum(cwt_anchor$matched_pathway),
    paste0("(", round(100 * mean(cwt_anchor$matched_pathway), 1), "%)"), "\n")
cat("  DTT completeness:", round(100 * mean(!is.na(cwt_anchor$cwt_dtt_date)), 1), "%\n")
cat("  MDT completeness:", round(100 * mean(!is.na(cwt_anchor$cwt_mdt_date)), 1), "%\n")
cat("days_dx_to_dtt summary:\n"); print(summary(cwt_anchor$days_dx_to_dtt))

# -----------------------------------------------------------------------------
# 5b. Legitimate decision-to-treat multiplicity
# How many patients have more than one legitimate DTT (in-window AND pathway-
# consistent AND on a distinct date) recorded in their OG journey. This is the
# multiplicity the anchor selection collapses over - distinct from the raw
# records-per-patient in 3b, which includes subsequent and inconsistent rows.
# -----------------------------------------------------------------------------
legit_dtt <- cwt_candidates %>%
  group_by(pseudo_patientid) %>%
  summarise(n_legit_dtt    = n_distinct(cwt_dtt_date),
            n_legit_groups = n_distinct(mod_group),
            .groups = "drop")

cat("\nLegitimate (in-window, pathway-consistent) distinct DTT dates per patient:\n")
legit_dtt %>% count(n_legit_dtt, name = "n_patients") %>%
  mutate(pct = round(100 * n_patients / sum(n_patients), 1)) %>%
  arrange(n_legit_dtt) %>% print(n = 30)
cat("  patients with >1 legitimate DTT date:", sum(legit_dtt$n_legit_dtt > 1),
    paste0("(", round(100 * mean(legit_dtt$n_legit_dtt > 1), 1), "% of anchored)"), "\n")
cat("  of those, span >1 modality group:",
    sum(legit_dtt$n_legit_dtt > 1 & legit_dtt$n_legit_groups > 1), "\n")

cat("\nMultiple legitimate DTTs by pathway (share with >1 distinct DTT):\n")
cwt_candidates %>%
  group_by(tx_pathway, pseudo_patientid) %>%
  summarise(n_dtt = n_distinct(cwt_dtt_date), .groups = "drop") %>%
  group_by(tx_pathway) %>%
  summarise(n = n(),
            pct_multi = round(100 * mean(n_dtt > 1), 1),
            .groups = "drop") %>%
  arrange(desc(pct_multi)) %>% print(n = 20)

# -----------------------------------------------------------------------------
# 6. Validation against the HES/SACT/RTDS treatment dates
# -----------------------------------------------------------------------------
cwt_validation <- cwt_anchor %>%
  left_join(og_cohort %>% select(pseudo_patientid, first_tx_date, tx_pathway),
            by = "pseudo_patientid") %>%
  mutate(
    dtt_to_cwt_treat = as.integer(cwt_treat_date - cwt_dtt_date),
    dtt_to_tx        = as.integer(first_tx_date  - cwt_dtt_date),
    cwt_vs_first_tx  = as.integer(cwt_treat_date - first_tx_date)
  )

cat("\nDTT to CWT treat date (internal consistency):\n")
cwt_validation %>% filter(!is.na(dtt_to_cwt_treat)) %>%
  summarise(n = n(), pct_neg = round(100 * mean(dtt_to_cwt_treat < 0), 1),
            median = median(dtt_to_cwt_treat),
            p25 = quantile(dtt_to_cwt_treat, .25),
            p75 = quantile(dtt_to_cwt_treat, .75)) %>% print()

cat("\nDTT to first_tx_date:\n")
cwt_validation %>% filter(!is.na(dtt_to_tx)) %>%
  summarise(n = n(), pct_neg = round(100 * mean(dtt_to_tx < 0), 1),
            median = median(dtt_to_tx),
            p25 = quantile(dtt_to_tx, .25),
            p75 = quantile(dtt_to_tx, .75)) %>% print()

cat("\nCWT treat date vs first_tx_date:\n")
cwt_validation %>% filter(!is.na(cwt_vs_first_tx)) %>%
  summarise(n = n(), pct_exact = round(100 * mean(cwt_vs_first_tx == 0), 1),
            pct_within_14 = round(100 * mean(abs(cwt_vs_first_tx) <= 14), 1),
            median = median(cwt_vs_first_tx)) %>% print()

cat("\nNegative dtt_to_tx by pathway:\n")
cwt_validation %>% filter(!is.na(dtt_to_tx)) %>%
  group_by(tx_pathway) %>%
  summarise(n = n(), pct_neg = round(100 * mean(dtt_to_tx < 0), 1), .groups = "drop") %>%
  arrange(desc(pct_neg)) %>% print(n = 20)

# -----------------------------------------------------------------------------
# 7. Attach the DTT node to the cohort and derive intervals + validity
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
    # more than the tolerance; NA for EMR/ESD where DTT is less meaningful
    dtt_valid = !is.na(cwt_dtt_date) & wt_dx_to_dtt >= 0 &
      wt_dtt_to_tx >= -treat_tol_days,
    dtt_valid = if_else(tx_pathway %in% c("EMR/ESD only", "EMR/ESD then surgery"),
                        NA, dtt_valid)
  )

# -----------------------------------------------------------------------------
# 7b. Audit reporting categories, derived from tx_pathway
# Collapses the 13-level pathway onto the OG audit groupings (NOGCA Tables 3
# and 4) so the curative breakdown can be reproduced without losing the detail.
# Curative RT only is kept as its own bucket: radical RT alone can be curative
# for early oesophageal SCC, so it does not sit cleanly in either audit column.
# -----------------------------------------------------------------------------
og_cohort <- og_cohort %>%
  mutate(
    tx_modality_audit = case_when(
      tx_pathway == "Surgery only"                                          ~ "Surgery only",
      tx_pathway %in% c("Surgery + neoadjuvant chemo",
                        "Surgery + neoadjuvant chemoRT",
                        "Surgery + neoadjuvant RT",
                        "Surgery + adjuvant chemo",
                        "Surgery + other")                                  ~ "Surgery plus SACT/RT",
      tx_pathway %in% c("EMR/ESD only", "EMR/ESD then surgery")             ~ "EMR/ESD",
      tx_pathway == "Definitive chemoRT"                                    ~ "Definitive chemoRT",
      tx_pathway == "Curative RT only"                                      ~ "Curative RT only",
      tx_pathway %in% c("Palliative chemo + RT", "SACT only",
                        "Palliative RT only")                               ~ "Chemo/RT only (non-curative)",
      tx_pathway == "No treatment recorded"                                 ~ "No treatment recorded",
      TRUE                                                                  ~ NA_character_
    ),
    tx_modality_audit = factor(tx_modality_audit, levels = c(
      "Surgery only", "Surgery plus SACT/RT", "EMR/ESD",
      "Definitive chemoRT", "Curative RT only",
      "Chemo/RT only (non-curative)", "No treatment recorded")),
    
    tx_intent_audit = case_when(
      tx_modality_audit %in% c("Surgery only", "Surgery plus SACT/RT",
                               "EMR/ESD", "Definitive chemoRT",
                               "Curative RT only")                          ~ "Curative",
      tx_modality_audit == "Chemo/RT only (non-curative)"                   ~ "Non-curative",
      tx_modality_audit == "No treatment recorded"                         ~ "No treatment",
      TRUE                                                                  ~ NA_character_
    ),
    
    # received treatment within nine months of diagnosis, as the audit defines
    # it. Counts the CWT treatment date as well as first_tx_date, because
    # first_tx_date is the curative clock-stop (NA for palliative patients); a
    # palliative patient with an in-window CWT treatment still received treatment.
    received_any_tx            = tx_pathway != "No treatment recorded" &
      ( (!is.na(first_tx_date) & wt_dx_to_tx <= tx_window_days) |
          (!is.na(cwt_treat_date) &
             as.integer(cwt_treat_date - diagmdy) >= 0 &
             as.integer(cwt_treat_date - diagmdy) <= tx_window_days) ),
    received_curative_tx_audit = tx_intent_audit == "Curative" &
      !is.na(first_tx_date) & wt_dx_to_tx <= tx_window_days
  )

saveRDS(og_cohort, paste0(base_dir, "og_cohort_cwt_2015_2022.rds"))
cat("\nSaved og_cohort_cwt_2015_2022.rds (", nrow(og_cohort), "patients )\n")

# -----------------------------------------------------------------------------
# 8. Post-merge summaries
# -----------------------------------------------------------------------------
cat("\ndtt_valid by pathway:\n")
og_cohort %>% filter(!is.na(cwt_dtt_date)) %>%
  count(tx_pathway, dtt_valid) %>% group_by(tx_pathway) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  filter(!is.na(dtt_valid)) %>% arrange(tx_pathway) %>% print(n = 30)

cat("\nIntervals where dtt_valid:\n")
og_cohort %>%
  filter(dtt_valid == TRUE,
         !is.na(wt_endo_to_dtt), wt_endo_to_dtt >= 0,
         !is.na(wt_dtt_to_tx),   wt_dtt_to_tx   >= 0) %>%
  summarise(n = n(),
            median_endo_dtt = median(wt_endo_to_dtt),
            median_dtt_tx   = median(wt_dtt_to_tx)) %>% print()

# -----------------------------------------------------------------------------
# 9. Audit reproductions
# These mirror the NOGCA national audit tables so the cohort can be eyeballed
# against published figures. The audit covers diagnoses in 2022-2023; set
# audit_window below to restrict to those years, or leave NULL for the full span.
# -----------------------------------------------------------------------------
audit_window <- NULL   # e.g. c(2022L, 2023L) to match the published period

aud <- og_cohort
if (!is.null(audit_window)) {
  aud <- aud %>% filter(ydiag >= audit_window[1], ydiag <= audit_window[2])
  cat("\nAudit window restricted to", audit_window[1], "-", audit_window[2],
      "(", nrow(aud), "patients )\n")
}

# --- Table 3: endoscopy -> first treatment, by treatment modality ------------
# Median (IQR) days from diagnostic endoscopy to start of treatment, among
# treated patients with a usable endoscopy date and a non-negative interval.
cat("\n--- Table 3: endoscopy to treatment, days, by modality ---\n")
cat("(audit England all-curative ~ surgery only 69, surgery+chemo/RT 60, ",
    "EMR/ESD 78, definitive chemoRT 67; chemo/RT only 62)\n")
aud %>%
  filter(!is.na(tx_modality_audit),
         tx_modality_audit != "No treatment recorded",
         !is.na(wt_endo_to_tx), wt_endo_to_tx >= 0) %>%
  group_by(tx_modality_audit) %>%
  summarise(n      = n(),
            median = median(wt_endo_to_tx),
            p25    = quantile(wt_endo_to_tx, .25),
            p75    = quantile(wt_endo_to_tx, .75),
            .groups = "drop") %>%
  arrange(tx_modality_audit) %>% print()

# --- Table 4: % treated within nine months, by stage and subtype -------------
# The audit denominator is stage 1-3; cohort is already stage 1-3, but the
# stage filter is kept explicit. Percentages are of all patients in the
# stage/subtype cell (treated and untreated).
cat("\n--- Table 4: % receiving treatment within 9 months, by subtype ---\n")
cat("(audit England stage 1-3: overall curative ~53%, overall any ~76%)\n")

table4 <- aud %>%
  filter(stage_clean %in% c("1", "2", "3")) %>%
  mutate(subtype = coalesce(cancer_subtype, "Unknown")) %>%
  group_by(subtype) %>%
  summarise(
    n_people              = n(),
    pct_surgery_only      = round(100 * mean(tx_modality_audit == "Surgery only", na.rm = TRUE), 0),
    pct_surgery_plus      = round(100 * mean(tx_modality_audit == "Surgery plus SACT/RT", na.rm = TRUE), 0),
    pct_definitive_chemRT = round(100 * mean(tx_modality_audit == "Definitive chemoRT", na.rm = TRUE), 0),
    pct_curative_rt_only  = round(100 * mean(tx_modality_audit == "Curative RT only", na.rm = TRUE), 0),
    pct_emresd            = round(100 * mean(tx_modality_audit == "EMR/ESD", na.rm = TRUE), 0),
    pct_curative_overall  = round(100 * mean(received_curative_tx_audit, na.rm = TRUE), 0),
    pct_any_treatment     = round(100 * mean(received_any_tx, na.rm = TRUE), 0),
    .groups = "drop"
  )

# add an All row across subtypes
table4_all <- aud %>%
  filter(stage_clean %in% c("1", "2", "3")) %>%
  summarise(
    subtype               = "All",
    n_people              = n(),
    pct_surgery_only      = round(100 * mean(tx_modality_audit == "Surgery only", na.rm = TRUE), 0),
    pct_surgery_plus      = round(100 * mean(tx_modality_audit == "Surgery plus SACT/RT", na.rm = TRUE), 0),
    pct_definitive_chemRT = round(100 * mean(tx_modality_audit == "Definitive chemoRT", na.rm = TRUE), 0),
    pct_curative_rt_only  = round(100 * mean(tx_modality_audit == "Curative RT only", na.rm = TRUE), 0),
    pct_emresd            = round(100 * mean(tx_modality_audit == "EMR/ESD", na.rm = TRUE), 0),
    pct_curative_overall  = round(100 * mean(received_curative_tx_audit, na.rm = TRUE), 0),
    pct_any_treatment     = round(100 * mean(received_any_tx, na.rm = TRUE), 0)
  )

bind_rows(table4_all, table4) %>% print(width = Inf)
