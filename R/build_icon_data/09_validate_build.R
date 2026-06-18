# =============================================================================
# 09  Validate the build  (print-only)
# -----------------------------------------------------------------------------
# Regression checks to run after a rebuild. Print-only - it reads the final
# cohort and reports, but writes nothing. Three sections:
#   A. shape and merge validity (dtt_valid by pathway; the anti-regression
#      guards for received_any_tx and the neoadjuvant anchor)
#   B. audit Tables 3 and 4 against the published NOGCA figures
#   C. CWT modality vs pathway-derived treatment concordance
#
# The one-off leakage investigation that established the surgery/chemo/RT
# coverage gaps is not re-run here; it lives in the archive. These three sections
# are what a routine rebuild needs.
#
# Reads: Data/ICON/og_cohort_cwt_2015_2022.rds
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

og <- readRDS(f_cohort_cwt)
section <- function(x) cat("\n\n==========  ", x, "  ==========\n")

# =============================================================================
# A. Shape and merge validity
# =============================================================================
section("A. shape and dtt_valid")

cat("patients:", nrow(og),
    "| with CWT DTT:", sum(!is.na(og$cwt_dtt_date)),
    "| with CWT treat date:", sum(!is.na(og$cwt_treat_date)), "\n")

cat("\nstage_clean (audit denominator is stage 1-3):\n")
og %>% count(stage_clean) %>% print()

cat("\ndtt_valid TRUE share by pathway (EMR/ESD is NA by design):\n")
og %>% filter(!is.na(cwt_dtt_date)) %>%
  count(tx_pathway, dtt_valid) %>% group_by(tx_pathway) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  filter(dtt_valid == TRUE) %>% arrange(desc(pct)) %>% print(n = 30)

# anti-regression: the any-treatment and curative flags must differ. If they are
# identical, the cwt_treat_date term in received_any_tx has been lost and every
# palliative patient has dropped out.
chk_flags <- identical(og$received_any_tx, og$received_curative_tx_audit)
cat("\nidentical(received_any_tx, received_curative_tx_audit):", chk_flags,
    if (chk_flags) " <- REGRESSION\n" else " (ok - they differ as intended)\n")

# anti-regression: neoadjuvant patients should anchor on a chemo/RT CWT record,
# not surgery (the primary-modality tie-break). Surgery anchors should be rare.
cat("\nneoadjuvant CWT anchor (surgery anchors should be rare):\n")
og %>%
  filter(tx_pathway %in% c("Surgery + neoadjuvant chemo",
                           "Surgery + neoadjuvant chemoRT",
                           "Surgery + neoadjuvant RT"),
         !is.na(cwt_modality)) %>%
  mutate(anchor = if_else(cwt_modality %in% c("01","23","24"),
                          "surgery (residual)", "chemo/RT (expected)")) %>%
  count(tx_pathway, anchor) %>% group_by(tx_pathway) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(tx_pathway, desc(anchor)) %>% print(n = 20)

# =============================================================================
# B. Audit Tables 3 and 4
# =============================================================================
section("B. audit reproduction (Tables 3 and 4)")

cat("Table 3 - endoscopy to treatment, days",
    "(England: surgery only 69, surgery+chemo/RT 60, EMR/ESD 78, def chemoRT 67):\n")
og %>%
  filter(!is.na(tx_modality_audit), tx_modality_audit != "No treatment recorded",
         !is.na(wt_endo_to_tx), wt_endo_to_tx >= 0) %>%
  group_by(tx_modality_audit) %>%
  summarise(n = n(), median = median(wt_endo_to_tx),
            p25 = quantile(wt_endo_to_tx, .25),
            p75 = quantile(wt_endo_to_tx, .75), .groups = "drop") %>%
  arrange(tx_modality_audit) %>% print()

cat("\nTable 4 - % treated within 9 months, stage 1-3",
    "(England: curative ~53%, any ~76%):\n")
t4_one <- function(d) d %>% summarise(
  n_people              = n(),
  pct_surgery_only      = round(100 * mean(tx_modality_audit == "Surgery only", na.rm = TRUE)),
  pct_surgery_plus      = round(100 * mean(tx_modality_audit == "Surgery plus SACT/RT", na.rm = TRUE)),
  pct_definitive_chemRT = round(100 * mean(tx_modality_audit == "Definitive chemoRT", na.rm = TRUE)),
  pct_curative_rt_only  = round(100 * mean(tx_modality_audit == "Curative RT only", na.rm = TRUE)),
  pct_emresd            = round(100 * mean(tx_modality_audit == "EMR/ESD", na.rm = TRUE)),
  pct_curative_overall  = round(100 * mean(received_curative_tx_audit, na.rm = TRUE)),
  pct_any_treatment     = round(100 * mean(received_any_tx, na.rm = TRUE)))
aud13 <- og %>% filter(stage_clean %in% c("1", "2", "3"))
bind_rows(
  aud13 %>% t4_one() %>% mutate(subtype = "All"),
  aud13 %>% mutate(subtype = coalesce(cancer_subtype, "Unknown")) %>%
    group_by(subtype) %>% t4_one() %>% ungroup()
) %>% relocate(subtype) %>% print(width = Inf)

# =============================================================================
# C. CWT modality vs pathway-derived treatment concordance
# =============================================================================
section("C. CWT modality vs HES/SACT/RTDS concordance")

cwt_group_of <- function(m) case_when(
  m %in% c("01","23","24") ~ "surgery", m %in% c("02","14","15") ~ "chemo",
  m == "03" ~ "hormone", m == "04" ~ "chemoRT",
  m %in% c("05","06","13") ~ "radiotherapy", m %in% c("07","08","09") ~ "palliative/AM",
  m == "97" ~ "other", m == "98" ~ "declined", TRUE ~ NA_character_)
pathway_group_of <- function(p) case_when(
  p %in% c("EMR/ESD only","EMR/ESD then surgery","Surgery only",
           "Surgery + other","Surgery + adjuvant chemo")        ~ "surgery",
  p == "Surgery + neoadjuvant chemo"                            ~ "chemo",
  p == "Surgery + neoadjuvant chemoRT"                          ~ "chemoRT",
  p == "Surgery + neoadjuvant RT"                               ~ "radiotherapy",
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

cat("Coverage overlap (any CWT anchor vs any pathway treatment):\n")
conc %>% count(has_cwt, has_path) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>% print()

both <- conc %>% filter(has_cwt, has_path, !is.na(cwt_grp), !is.na(path_grp))
cat("\nPatients with both:", nrow(both), "- modality-group agreement:\n")
both %>%
  mutate(agree = cwt_grp == path_grp |
           # neoadjuvant chemo/RT before a CWT surgical clock-stop is concordant
           (cwt_grp == "surgery" & path_grp %in% c("chemo","chemoRT","radiotherapy"))) %>%
  count(agree) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% print()

cat("\n09 validation complete.\n")
