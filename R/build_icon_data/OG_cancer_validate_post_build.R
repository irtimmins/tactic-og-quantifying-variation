# =============================================================================
# OG cancer - post-build validation
# -----------------------------------------------------------------------------
# Runs after the dataset is built (og_cohort_cwt_2015_2022.rds). Consolidates the
# leakage investigation and the audit reproduction into one print-only report, so
# a rebuild can be regression-checked in a single pass. Sections:
#
#   A. dataset shape and the merge validity (dtt_valid) by pathway
#   B. audit Table 3 (endoscopy -> treatment) and Table 4 (% treated)
#   C. leakage stage 1 - localise CWT-treated but pathway-untreated patients
#   D. leakage stage 2 - are the CWT-"surgical" leaks actually resections?
#   E. leakage stage 3 - characterise the endoscopy-only group (stent/enabling)
#   F. stage-1 stent sanity check - is the early-stage-untreated group real?
#
# Conclusions reached during the investigation, confirmed by the numbers below:
#   - the surgical "leak" is mostly palliative/enabling endoscopy (stents) that
#     CWT codes as surgery; these are correctly non-curative.
#   - the only recoverable curative cases were ~123 emergency resections, now
#     folded in by dropping the !emergency filter in the surgery anchor.
#   - received_any_tx must count cwt_treat_date (not just first_tx_date), or it
#     collapses onto the curative flag and drops every palliative patient.
#   - chemo/RT leakage is a SACT/RTDS coverage gap, not a pathway-logic fault.
# =============================================================================

library(tidyverse)

base_dir       <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/"
near_days      <- 60L     # proximity of a HES op to the CWT treatment date
tx_window_days <- 270L    # CWT-merge window (script 2 surgery window is 275)
surg_window    <- 275L

og      <- readRDS(paste0(base_dir, "og_cohort_cwt_2015_2022.rds"))
hes_apc <- readRDS(paste0(base_dir, "hes_apc_og_2014_2022.rds"))

section <- function(x) cat("\n\n==================  ", x, "  ==================\n")

# code lists, verbatim from the build (script 2)
opcs_resection <- c("G011","G012","G013","G018","G019",
                    "G021","G022","G023","G024","G025","G028","G029",
                    "G031","G032","G033","G034","G035","G036","G038","G039",
                    "G271","G272","G273","G274","G275","G278","G279",
                    "G281","G282","G283","G288","G289")
opcs_emresd <- c("G121","G128","G129","G141","G146","G148","G149",
                 "G171","G178","G179","G421","G423","G428","G429",
                 "G431","G438","G439","G143","G145","G433","G435")
opcs_diagnostic_endoscopy <- c(
  "G142","G143","G145","G147","G152","G153","G154","G156","G157","G158","G159",
  "G161","G162","G168","G169","G172","G173","G188","G189","G191","G198","G199",
  "G201","G202","G208","G209","G214","G215","G218","G219",
  "G422","G432","G433","G435","G441","G443","G445","G446","G448","G449",
  "G451","G452","G454","G458","G459","G462","G463","G468","G469")
opcs_stent_enabling <- c("G141","G142","G143","G144","G145","G146","G147","G148",
                         "G441","G442","G443","G444","G445","G446","G447","G448")
admimeth_emerg <- c("21","22","23","24","25","28","2A","2B","2C","2D")

# -----------------------------------------------------------------------------
# Shared helper: HES operations near a set of patients' CWT dates, long form
# Pivots OPERTN_01..24 with paired OPDATE_01..24, carries ADMIMETH, one row per
# (patient, operation). Used by every leakage section so the pivot lives once.
# -----------------------------------------------------------------------------
hes_ops_long <- function(ids) {
  op <- hes_apc %>%
    mutate(STUDY_ID = as.character(STUDY_ID),
           ADMIMETH = as.character(ADMIMETH)) %>%
    filter(STUDY_ID %in% ids) %>%
    select(STUDY_ID, EPISTART, ADMIMETH, starts_with("OPERTN_")) %>%
    pivot_longer(starts_with("OPERTN_"), names_to = "pos",
                 values_to = "opcs", names_prefix = "OPERTN_")
  dt <- hes_apc %>%
    mutate(STUDY_ID = as.character(STUDY_ID)) %>%
    filter(STUDY_ID %in% ids) %>%
    select(STUDY_ID, EPISTART, starts_with("OPDATE_")) %>%
    pivot_longer(starts_with("OPDATE_"), names_to = "pos",
                 values_to = "opdate", names_prefix = "OPDATE_")
  op %>%
    left_join(dt, by = c("STUDY_ID", "EPISTART", "pos"),
              relationship = "many-to-many") %>%
    filter(!is.na(opcs), opcs != "", opcs != "-", opcs != "&") %>%
    mutate(opcs4     = str_to_upper(str_remove_all(str_trim(opcs), "\\.")),
           op_date   = coalesce(as.Date(opdate), as.Date(EPISTART)),
           emergency = ADMIMETH %in% admimeth_emerg) %>%
    select(STUDY_ID, op_date, opcs4, ADMIMETH, emergency)
}

# =============================================================================
# A. Dataset shape and merge validity
# =============================================================================
section("A. dataset shape and dtt_valid")

cat("patients:", nrow(og),
    "| with CWT DTT:", sum(!is.na(og$cwt_dtt_date)),
    "| with CWT treat date:", sum(!is.na(og$cwt_treat_date)), "\n")

cat("\nstage_clean (audit denominator should be stage 1-3 only):\n")
og %>% count(stage_clean) %>% print()

cat("\ndtt_valid by pathway (TRUE share; EMR/ESD excluded as NA by design):\n")
og %>% filter(!is.na(cwt_dtt_date)) %>%
  count(tx_pathway, dtt_valid) %>% group_by(tx_pathway) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  filter(dtt_valid == TRUE) %>% arrange(desc(pct)) %>% print(n = 30)

cat("\nany-treatment vs curative flag - must NOT be identical:\n")
chk_flags <- identical(og$received_any_tx, og$received_curative_tx_audit)
cat("identical(received_any_tx, received_curative_tx_audit):", chk_flags,
    if (chk_flags) " <- REGRESSION: the cwt_treat_date fix has been lost\n"
    else " (ok)\n")

cat("\nneoadjuvant anchor: CWT modality among neoadjuvant pathways\n")
cat("(after the tie-break fix, surgery-modality (01/23/24) anchors should be\n",
    " rare - only patients with no in-window chemo/RT CWT row):\n")
og %>%
  filter(tx_pathway %in% c("Surgery + neoadjuvant chemo",
                           "Surgery + neoadjuvant chemoRT",
                           "Surgery + neoadjuvant RT"),
         !is.na(cwt_modality)) %>%
  mutate(anchor = if_else(cwt_modality %in% c("01","23","24"),
                          "surgery (residual)", "chemo/RT (expected)")) %>%
  count(tx_pathway, anchor) %>%
  group_by(tx_pathway) %>% mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(tx_pathway, desc(anchor)) %>% print(n = 20)

# Definitive chemoRT carries the lowest dtt_valid share. The tie-break anchors
# these on the chemo/chemoRT DTT, which can sit more than treat_tol_days before
# first_tx_date (= pmin(sact, rt)) because the two arms of a combined course
# start on different days. This characterises the invalid cases: clustered just
# past -treat_tol_days -> the tolerance is slightly tight for this pathway;
# scattered far negative -> genuine timing oddities the rule rightly excludes.
treat_tol_days <- 14L
cat("\nDefinitive chemoRT: wt_dtt_to_tx where dtt_valid is FALSE\n")
cat("(validity boundary is -", treat_tol_days, "d; values just past it suggest\n",
    " the tolerance is tight rather than the timing being wrong):\n", sep = "")
dchemort_bad <- og %>%
  filter(tx_pathway == "Definitive chemoRT", dtt_valid == FALSE,
         !is.na(wt_dtt_to_tx))
cat("n invalid with a wt_dtt_to_tx:", nrow(dchemort_bad), "\n")
print(summary(dchemort_bad$wt_dtt_to_tx))
cat("\nbanding of wt_dtt_to_tx among the invalid:\n")
dchemort_bad %>%
  mutate(band = cut(wt_dtt_to_tx,
                    breaks = c(-Inf, -60, -30, -15, -1, Inf),
                    labels = c("<= -60", "-59..-30", "-29..-15",
                               "-14..-1 (within wider tol)", ">= 0"))) %>%
  count(band) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% print()
cat("\nFor comparison, the same banding under a 30-day tolerance would reclass\n")
cat("  newly-valid (-29..-15):",
    sum(dchemort_bad$wt_dtt_to_tx >= -30 & dchemort_bad$wt_dtt_to_tx < -treat_tol_days),
    "of", nrow(dchemort_bad), "currently-invalid definitive chemoRT cases\n")

# =============================================================================
# B. Audit Table 3 and Table 4
# =============================================================================
section("B. audit reproduction (Tables 3 and 4)")
audit_window <- NULL   # set c(2022L, 2023L) to match the published period
aud <- og
if (!is.null(audit_window))
  aud <- aud %>% filter(ydiag >= audit_window[1], ydiag <= audit_window[2])

cat("Table 3 - endoscopy to treatment, days (audit England: surgery only 69,\n",
    "  surgery+chemo/RT 60, EMR/ESD 78, definitive chemoRT 67):\n")
aud %>%
  filter(!is.na(tx_modality_audit),
         tx_modality_audit != "No treatment recorded",
         !is.na(wt_endo_to_tx), wt_endo_to_tx >= 0) %>%
  group_by(tx_modality_audit) %>%
  summarise(n = n(), median = median(wt_endo_to_tx),
            p25 = quantile(wt_endo_to_tx, .25),
            p75 = quantile(wt_endo_to_tx, .75), .groups = "drop") %>%
  arrange(tx_modality_audit) %>% print()

cat("\nTable 4 - % treated within 9 months, stage 1-3 (audit England:\n",
    "  overall curative ~53%, overall any ~76%):\n")
t4_one <- function(d) d %>% summarise(
  n_people             = n(),
  pct_surgery_only     = round(100 * mean(tx_modality_audit == "Surgery only", na.rm = TRUE)),
  pct_surgery_plus     = round(100 * mean(tx_modality_audit == "Surgery plus SACT/RT", na.rm = TRUE)),
  pct_definitive_chemRT= round(100 * mean(tx_modality_audit == "Definitive chemoRT", na.rm = TRUE)),
  pct_curative_rt_only = round(100 * mean(tx_modality_audit == "Curative RT only", na.rm = TRUE)),
  pct_emresd           = round(100 * mean(tx_modality_audit == "EMR/ESD", na.rm = TRUE)),
  pct_curative_overall = round(100 * mean(received_curative_tx_audit, na.rm = TRUE)),
  pct_any_treatment    = round(100 * mean(received_any_tx, na.rm = TRUE))
)
aud13 <- aud %>% filter(stage_clean %in% c("1", "2", "3"))
bind_rows(
  aud13 %>% t4_one() %>% mutate(subtype = "All"),
  aud13 %>% mutate(subtype = coalesce(cancer_subtype, "Unknown")) %>%
    group_by(subtype) %>% t4_one() %>% ungroup()
) %>% relocate(subtype) %>% print(width = Inf)

# =============================================================================
# C. Leakage stage 1 - localise CWT-treated but pathway-untreated patients
# =============================================================================
section("C. leakage stage 1 - where treated-in-CWT but untreated-in-pathway")

leak <- og %>%
  filter(tx_pathway == "No treatment recorded",
         !is.na(cwt_treat_date),
         cwt_modality %in% c("01","23","24","02","04","05")) %>%
  mutate(cwt_group = case_when(
    cwt_modality %in% c("01","23","24") ~ "surgery",
    cwt_modality == "02"                ~ "chemo",
    cwt_modality %in% c("04","05")      ~ "radiotherapy")) %>%
  select(pseudo_patientid, diagmdy, cwt_treat_date, cwt_modality, cwt_group)

cat("leaked patients by CWT treatment group:\n")
leak %>% count(cwt_group, sort = TRUE) %>% print()

# build the HES op table once for the surgical-leak set, reused in C/D/E.
# C's leak_surg and D's leak_surg_full apply the same filter (No treatment +
# cwt_treat_date + modality 01/23/24), so this id set covers both sections.
leak_surg <- leak %>% filter(cwt_group == "surgery")
hes_long  <- hes_ops_long(leak_surg$pseudo_patientid)

surg_dx <- leak_surg %>% select(pseudo_patientid, diagmdy, cwt_treat_date)

surg_evidence <- surg_dx %>%
  left_join(hes_long, by = c("pseudo_patientid" = "STUDY_ID")) %>%
  mutate(near      = !is.na(op_date) &
           abs(as.integer(op_date - cwt_treat_date)) <= near_days,
         is_resn   = opcs4 %in% opcs_resection,
         dx_to_op  = as.integer(op_date - diagmdy),
         in_window = !is.na(dx_to_op) & dx_to_op >= -30 & dx_to_op <= surg_window) %>%
  group_by(pseudo_patientid) %>%
  summarise(any_opcs_near        = any(near, na.rm = TRUE),
            any_resection_near   = any(near & is_resn, na.rm = TRUE),
            resn_window_emerg    = any(is_resn & in_window &  emergency, na.rm = TRUE),
            resn_elective_oow    = any(is_resn & !in_window & !emergency, na.rm = TRUE),
            resn_anchor_eligible = any(is_resn & in_window & !emergency, na.rm = TRUE),
            .groups = "drop")

cat("\nsurgical leak attribution (n =", nrow(leak_surg), "):\n")
surg_evidence %>% summarise(
  has_any_opcs_near         = sum(any_opcs_near),
  has_resection_near        = sum(any_resection_near),
  resn_window_emergency     = sum(resn_window_emerg),
  resn_elective_outofwindow = sum(resn_elective_oow),
  resn_anchor_eligible      = sum(resn_anchor_eligible),
  no_hes_activity_near      = sum(!any_opcs_near)
) %>% print(width = Inf)
cat("  resn_window_emergency  -> recoverable (emergency resection; now folded in\n",
    "                            by dropping the !emergency filter in script 2).\n",
    "  resn_anchor_eligible   -> should be ~0; non-zero means an ID/link issue.\n")

# chemo / RT leakage: coverage vs filter, against the real anchor objects
chk_source <- function(grp, raw_path, anchor_path, date_col) {
  lk <- leak %>% filter(cwt_group == grp)
  cat("\n", grp, "leak (n =", nrow(lk), "):\n")
  if (!file.exists(raw_path) || nrow(lk) == 0) {
    cat("  raw object not found:", raw_path, "\n"); return(invisible())
  }
  raw <- readRDS(raw_path)
  anc <- if (file.exists(anchor_path)) readRDS(anchor_path)$pseudo_patientid else character(0)
  in_raw <- lk$pseudo_patientid %in% raw$pseudo_patientid
  lkw <- raw %>%
    filter(pseudo_patientid %in% lk$pseudo_patientid) %>%
    left_join(lk %>% select(pseudo_patientid, diagmdy), by = "pseudo_patientid") %>%
    mutate(d = as.integer(.data[[date_col]] - diagmdy),
           inw = !is.na(d) & d >= -30 & d <= surg_window)
  cat("  in raw extract:", sum(in_raw), "of", nrow(lk),
      "| in anchor:", sum(lk$pseudo_patientid %in% anc),
      "| in raw & in-window:", n_distinct(lkw$pseudo_patientid[lkw$inw]), "\n")
  cat("  (low 'in raw extract' -> coverage gap, not a pathway fault.)\n")
}
chk_source("chemo",
           paste0(base_dir, "sact_og_2012_2024.rds"),
           paste0(base_dir, "og_sact_anchor_2015_2022.rds"),
           "sact_regimen_date")
chk_source("radiotherapy",
           paste0(base_dir, "rtds_og_2009_2024.rds"),
           paste0(base_dir, "rt_anchor_og.rds"),
           "rt_start_date")

# =============================================================================
# D. Leakage stage 2 - are the CWT-"surgical" leaks actually resections?
# =============================================================================
section("D. leakage stage 2 - surgical leak composition")

leak_surg_full <- og %>%
  filter(tx_pathway == "No treatment recorded",
         !is.na(cwt_treat_date),
         cwt_modality %in% c("01","23","24")) %>%
  select(pseudo_patientid, cwt_treat_date, emresd_date,
         endoscopy_date, surgery_date)

cat("anchors already present (any > 0 = a pathway-logic miss, not data):\n")
leak_surg_full %>% summarise(
  has_emresd_date    = sum(!is.na(emresd_date)),
  has_endoscopy_date = sum(!is.na(endoscopy_date)),
  has_surgery_date   = sum(!is.na(surgery_date))) %>% print(width = Inf)

near_surg <- leak_surg_full %>%
  select(pseudo_patientid, cwt_treat_date) %>%
  left_join(hes_long, by = c("pseudo_patientid" = "STUDY_ID")) %>%
  filter(!is.na(op_date),
         abs(as.integer(op_date - cwt_treat_date)) <= near_days)

per_patient <- near_surg %>%
  group_by(pseudo_patientid) %>%
  summarise(any_resection = any(opcs4 %in% opcs_resection),
            any_emresd    = any(opcs4 %in% opcs_emresd),
            any_diag_endo = any(opcs4 %in% opcs_diagnostic_endoscopy),
            .groups = "drop")

cat("\nlayered attribution over all", nrow(leak_surg_full),
    "leaked CWT-surgical patients:\n")
leak_surg_full %>%
  left_join(per_patient, by = "pseudo_patientid") %>%
  mutate(category = case_when(
    !is.na(surgery_date)           ~ "1 has surgery_date (pathway-logic miss)",
    !is.na(emresd_date)            ~ "2 has emresd_date (pathway-logic miss)",
    coalesce(any_resection, FALSE) ~ "3 resection in HES, not anchored (emerg/window)",
    coalesce(any_emresd, FALSE)    ~ "4 EMR/ESD in HES, not anchored",
    coalesce(any_diag_endo, FALSE) ~ "5 only diagnostic endoscopy near CWT date",
    TRUE                           ~ "6 no resection/EMR/endoscopy trace")) %>%
  count(category) %>% mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(category) %>% print(width = Inf)
cat("  expectation: category 5 dominates (~90%); categories 1-2 should be ~0.\n")

# =============================================================================
# E. Leakage stage 3 - characterise the endoscopy-only group
# =============================================================================
section("E. leakage stage 3 - endoscopy-only group profile")

endo_only_ids <- per_patient %>%
  filter(!any_resection, !any_emresd) %>% pull(pseudo_patientid)

endo_only <- og %>%
  filter(pseudo_patientid %in% endo_only_ids) %>%
  select(pseudo_patientid, cwt_modality, stage_clean, cancer_subtype, agediag)

cat("endoscopy-only patients:", nrow(endo_only), "\n")
cat("\nstage:\n");    endo_only %>% count(stage_clean) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>% print()
cat("\nCWT modality (01 = old surgery, 23/24 incl. enabling):\n")
endo_only %>% count(cwt_modality) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>% print()
cat("\nsubtype:\n");  endo_only %>% count(cancer_subtype) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>% print()
cat("\nage at diagnosis:\n"); print(summary(endo_only$agediag))

cat("\nOPCS at the CWT date (top 15; stent/enabling flagged):\n")
near_surg %>%
  filter(pseudo_patientid %in% endo_only_ids) %>%
  count(opcs4, sort = TRUE) %>%
  mutate(stent_enabling = opcs4 %in% opcs_stent_enabling) %>%
  print(n = 15)
cat("  expectation: stage 3 + high age + stent codes (G441 oesophageal stent)\n",
    "  -> palliative/enabling endoscopy, correctly non-curative.\n")

# =============================================================================
# F. Stage-1 stent sanity check
# =============================================================================
section("F. stage-1 endoscopy-only sanity check")
# A stage-1 patient receiving only a stent is clinically incongruous: either
# genuinely early-stage but undertreated (worth a footnote), or mis-staged
# advanced disease (a data artefact). basisofdiagnosis and short-term mortality
# tell the two apart.

stage1_endo <- og %>%
  filter(pseudo_patientid %in% endo_only_ids, stage_clean == "1") %>%
  mutate(surv_days = as.integer(finmdy - diagmdy),
         died_90d  = !is.na(surv_days) & died == 1L & surv_days <= 90)

cat("stage-1 endoscopy-only patients:", nrow(stage1_endo), "\n")
if (nrow(stage1_endo) > 0) {
  cat("\nbasis of diagnosis:\n")
  stage1_endo %>% count(basisofdiagnosis, sort = TRUE) %>% print(n = 20)
  cat("\nsubtype:\n")
  stage1_endo %>% count(cancer_subtype) %>% print()
  cat("\ndied within 90 days of diagnosis:\n")
  stage1_endo %>% count(died_90d) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>% print()
  cat("\nsurvival from diagnosis (days):\n")
  print(summary(stage1_endo$surv_days))
  cat("\n  high 90-day mortality or non-histological basis -> likely mis-staged\n",
      "  advanced disease (artefact). Long survival + histological basis ->\n",
      "  genuinely early-stage with only palliative/enabling care (a real,\n",
      "  reportable finding worth a methods footnote).\n")
}

# =============================================================================
# G. CWT modality vs pathway-derived treatment - concordance
# -----------------------------------------------------------------------------
# Direct answer to: "is there perfect alignment between CWT modality and what we
# flagged as treatment from HES/SACT/RTDS?" No - and this quantifies the gap.
# For every cohort patient we compare the CWT anchor's modality group against the
# treatment group implied by tx_pathway, and report agreement, the direction of
# disagreement, and what the disagreements are.
#
# Two distinct sources have a CWT modality view and a pathway view:
#   - CWT side : the modality of the anchored CWT record (cwt_modality), mapped
#                to the same broad groups used in the merge.
#   - source side : tx_pathway, mapped to the modality group of its clock-stop.
# =============================================================================
section("G. CWT modality vs HES/SACT/RTDS treatment concordance")

# ANSWER to the colleague's question (is CWT modality perfectly aligned with the
# HES/SACT/RTDS treatment flags?): no, and the gap is structural, not error.
# Among patients with both a CWT record and a pathway treatment, modality groups
# agree ~87%. The disagreement has two benign causes: CWT counts palliative /
# enabling events the curative flags exclude (the ~18% CWT-only group), and CWT
# anchors only first-definitive-treatment so misses some subsequent / unmatched
# treatments the flags catch (the ~4% pathway-only group). Neoadjuvant sequencing
# (CWT clock-stops on surgery, pathway on the neoadjuvant chemo/RT) accounts for
# most of the remaining off-diagonal and is correct behaviour.

cwt_group_of <- function(m) case_when(
  m %in% c("01","23","24")      ~ "surgery",
  m %in% c("02","14","15")      ~ "chemo",
  m == "03"                     ~ "hormone",
  m == "04"                     ~ "chemoRT",
  m %in% c("05","06","13")      ~ "radiotherapy",
  m %in% c("07","08","09")      ~ "palliative/AM",
  m == "97"                     ~ "other",
  m == "98"                     ~ "declined",
  TRUE                          ~ NA_character_)

# the dominant modality the pathway implies (its curative/defining clock-stop)
pathway_group_of <- function(p) case_when(
  p %in% c("EMR/ESD only","EMR/ESD then surgery","Surgery only",
           "Surgery + other","Surgery + adjuvant chemo")        ~ "surgery",
  p %in% c("Surgery + neoadjuvant chemo")                       ~ "chemo",
  p %in% c("Surgery + neoadjuvant chemoRT")                     ~ "chemoRT",
  p %in% c("Surgery + neoadjuvant RT")                          ~ "radiotherapy",
  p == "Definitive chemoRT"                                     ~ "chemoRT",
  p %in% c("Curative RT only","Palliative RT only")             ~ "radiotherapy",
  p == "SACT only"                                              ~ "chemo",
  p == "Palliative chemo + RT"                                  ~ "chemo",
  p == "No treatment recorded"                                  ~ "none",
  TRUE                                                          ~ NA_character_)

conc <- og %>%
  mutate(cwt_grp  = cwt_group_of(cwt_modality),
         path_grp = pathway_group_of(tx_pathway),
         has_cwt  = !is.na(cwt_modality),
         has_path = tx_pathway != "No treatment recorded")

# 1. coverage overlap: who has a CWT record vs a pathway treatment, 2x2
cat("Coverage overlap (any CWT anchor vs any pathway treatment):\n")
conc %>% count(has_cwt, has_path) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>% print()
cat("  has_cwt & !has_path -> CWT says treated, pathway says none (the leak).\n",
    " !has_cwt &  has_path -> pathway found treatment with no CWT anchor.\n")

# 2. among patients with BOTH, do the modality groups agree?
both <- conc %>% filter(has_cwt, has_path, !is.na(cwt_grp), !is.na(path_grp))
cat("\nPatients with both a CWT modality and a pathway treatment:", nrow(both), "\n")
cat("modality-group agreement:\n")
both %>%
  mutate(agree = cwt_grp == path_grp |
           # neoadjuvant chemo/RT before surgery legitimately differs from
           # a CWT surgical clock-stop; count these as concordant
           (cwt_grp == "surgery" & path_grp %in% c("chemo","chemoRT","radiotherapy"))) %>%
  count(agree) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% print()

# 3. the full crosstab: CWT modality group (rows) x pathway group (cols)
cat("\nCWT modality group (row) x pathway-implied group (col), patients with both:\n")
both %>% count(cwt_grp, path_grp) %>%
  pivot_wider(names_from = path_grp, values_from = n, values_fill = 0) %>%
  arrange(cwt_grp) %>% print(width = Inf)

# 3b. same crosstab as row percentages: each CWT modality group sums to 100%,
#     so the diagonal cell is that modality's agreement rate with the pathway.
cat("\nRow %% (each CWT modality group sums to 100; diagonal = agreement rate):\n")
both %>% count(cwt_grp, path_grp) %>%
  group_by(cwt_grp) %>% mutate(pct = round(100 * n / sum(n))) %>% ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = path_grp, values_from = pct, values_fill = 0) %>%
  arrange(cwt_grp) %>% print(width = Inf)

# 4. the headline numbers to quote back
cat("\nHeadline figures:\n")
hl <- conc %>% summarise(
  n                    = n(),
  cwt_present          = sum(has_cwt),
  pathway_treated      = sum(has_path),
  cwt_only             = sum(has_cwt & !has_path),
  pathway_only         = sum(!has_cwt & has_path),
  both_present         = sum(has_cwt & has_path)
)
print(hl, width = Inf)
cat(sprintf(
  "  - %d patients (%.1f%%) have a CWT treatment record but no pathway treatment;\n    these are dominated by palliative/enabling endoscopy CWT scores as surgery.\n",
  hl$cwt_only, 100 * hl$cwt_only / hl$n))
cat(sprintf(
  "  - %d patients (%.1f%%) have a pathway treatment with no anchored CWT record;\n    CWT records first-definitive-treatment only, so subsequent-only or\n    unmatched-modality rows do not anchor.\n",
  hl$pathway_only, 100 * hl$pathway_only / hl$n))
cat("  - the two systems use different definitions (CWT: first definitive\n",
    "   treatment incl. palliative/enabling; pathway: curative HES/SACT/RTDS\n",
    "   anchors), so partial disagreement is expected, not an error.\n")

# 5. of the CWT-only group, the modality breakdown (what CWT calls treatment
#    that the pathway does not)
cat("\nCWT modality among CWT-treated / pathway-untreated patients:\n")
conc %>% filter(has_cwt, !has_path) %>%
  mutate(grp = cwt_group_of(cwt_modality)) %>%
  count(cwt_modality, grp, sort = TRUE) %>% print(n = 25)

section("end of validation")
