# =============================================================================
# 11  Check treatment classification  (optional, print-only)
# -----------------------------------------------------------------------------
# Complements 09. Where 09 checks the shape and the audit totals, this checks the
# classification decisions hold up under scrutiny - the things worth confirming
# before trusting the pathway assignment and the CWT anchoring:
#   A. is "No treatment recorded" genuinely untreated, or leaking CWT-treated cases?
#   B. CWT anchor and treatment-date coverage by pathway
#   C. consistency invariant: any-treatment is never below curative
#   D. neoadjuvant clock-stop: does the CWT date sit on the earlier chemo/RT event
#      (per CWT guidance 3.9.1), not the later surgery?
#   E. HES-only chemo: do the reclassified patients sit in clinically plausible
#      positions (neoadjuvant before surgery, concurrent with RT)?
#
# Reads the final cohort and the chemo anchor; writes nothing. Left out of
# 00_master by default - source it directly when you want the deeper assurance.
# The one-off investigations that established these (the residual-case
# characterisation, the original received_any_tx bug confirmation) are not
# reproduced; they are closed.
#
# Reads: Data/ICON/og_cohort_cwt_2015_2022.rds, og_chemo_anchor_2015_2022.rds
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

og  <- readRDS(f_cohort_cwt)
sec <- function(x) cat("\n\n==========  ", x, "  ==========\n")

# =============================================================================
# A. No-treatment leakage
# =============================================================================
sec("A. 'No treatment recorded' vs presence of a CWT treatment")

og %>%
  mutate(no_tx = tx_pathway == "No treatment recorded",
         has_cwt_treat = !is.na(cwt_treat_date)) %>%
  count(no_tx, has_cwt_treat) %>%
  group_by(no_tx) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% print()

cat("\nCWT modality among 'No treatment recorded' patients who have a CWT record\n",
    "(palliative 07/08/09 and other 97 are expected; surgery 01/23/24 or\n",
    " chemo/RT 02/04/05 would indicate genuine leakage from script 07):\n")
og %>%
  filter(tx_pathway == "No treatment recorded", !is.na(cwt_treat_date)) %>%
  count(cwt_modality, sort = TRUE) %>% print(n = 30)

# of any curative/active modalities that do leak, are they even in-window?
leak <- og %>%
  filter(tx_pathway == "No treatment recorded", !is.na(cwt_treat_date),
         cwt_modality %in% c("01","23","24","02","04","05")) %>%
  mutate(dx_to_cwttx = as.integer(cwt_treat_date - diagmdy))
cat("\nactive-modality leakage cases:", nrow(leak),
    "| of which within", cwt_window_days, "days of dx:",
    sum(leak$dx_to_cwttx >= 0 & leak$dx_to_cwttx <= cwt_window_days, na.rm = TRUE), "\n")

# =============================================================================
# B. CWT anchor and treatment-date coverage by pathway
# =============================================================================
sec("B. coverage by pathway (DTT, CWT treat date, first_tx_date)")

og %>%
  group_by(tx_pathway) %>%
  summarise(n = n(),
            pct_with_dtt      = round(100 * mean(!is.na(cwt_dtt_date)), 1),
            pct_with_cwt_tx   = round(100 * mean(!is.na(cwt_treat_date)), 1),
            pct_with_first_tx = round(100 * mean(!is.na(first_tx_date)), 1),
            .groups = "drop") %>%
  arrange(desc(n)) %>% print(n = 20)

# =============================================================================
# C. Consistency invariant: any-treatment >= curative in every subtype
# =============================================================================
sec("C. invariant - any-treatment never below curative")

inv <- og %>% filter(stage_clean %in% c("1","2","3")) %>%
  mutate(subtype = coalesce(cancer_subtype, "Unknown")) %>%
  group_by(subtype) %>%
  summarise(pct_curative = round(100 * mean(received_curative_tx_audit, na.rm = TRUE)),
            pct_any      = round(100 * mean(received_any_tx, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(ok = pct_any >= pct_curative)
inv_all <- og %>% filter(stage_clean %in% c("1","2","3")) %>%
  summarise(subtype = "All",
            pct_curative = round(100 * mean(received_curative_tx_audit, na.rm = TRUE)),
            pct_any      = round(100 * mean(received_any_tx, na.rm = TRUE))) %>%
  mutate(ok = pct_any >= pct_curative)
bind_rows(inv_all, inv) %>% print()
cat("(every ok should be TRUE; FALSE means the any/curative flags are inconsistent)\n")

# =============================================================================
# D. Neoadjuvant clock-stop: CWT date on the earlier chemo/RT, not surgery
# =============================================================================
sec("D. neoadjuvant clock-stop lands on the earlier treatment (guidance 3.9.1)")

neo <- og %>%
  filter(tx_pathway %in% c("Surgery + neoadjuvant chemo",
                           "Surgery + neoadjuvant chemoRT",
                           "Surgery + neoadjuvant RT"),
         !is.na(cwt_treat_date)) %>%
  mutate(neoadj_date = pmin(sact_date, rt_date, na.rm = TRUE),
         d_to_neoadj = as.integer(cwt_treat_date - neoadj_date),
         d_to_surg   = as.integer(cwt_treat_date - surgery_date),
         closest = case_when(
           is.na(neoadj_date) & is.na(surgery_date) ~ "neither",
           is.na(surgery_date)                      ~ "neoadjuvant",
           is.na(neoadj_date)                       ~ "surgery",
           abs(d_to_neoadj) <= abs(d_to_surg)       ~ "neoadjuvant",
           TRUE                                     ~ "surgery"))

cat("neoadjuvant patients with a CWT treat date:", nrow(neo), "\n\n")
cat("which event is cwt_treat_date closest to (neoadjuvant expected):\n")
neo %>% count(tx_pathway, closest) %>% group_by(tx_pathway) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(tx_pathway, desc(n)) %>% print(n = 30)

cat("\ndays from CWT treat date to the neoadjuvant event (0 = on it):\n")
neo %>% filter(!is.na(d_to_neoadj)) %>% group_by(tx_pathway) %>%
  summarise(n = n(), median = median(d_to_neoadj),
            pct_within_14 = round(100 * mean(abs(d_to_neoadj) <= 14), 1),
            .groups = "drop") %>% print()

# =============================================================================
# E. HES-only chemo: are the reclassified patients clinically plausible?
# =============================================================================
sec("E. HES-only chemo reclassification timing")

chemo_path <- f_chemo_anchor
if (file.exists(chemo_path)) {
  hes_only <- readRDS(chemo_path) %>%
    filter(chemo_source == "hes") %>%
    select(pseudo_patientid, hes_chemo_date, days_dx_to_chemo) %>%
    left_join(og %>% select(pseudo_patientid, tx_pathway, surgery_date,
                            rt_date, rt_curative, stage_clean),
              by = "pseudo_patientid")
  cat("HES-only chemo patients:", nrow(hes_only), "\n")
  
  # surgery patients: chemo before surgery is genuine neoadjuvant
  surg <- hes_only %>% filter(!is.na(surgery_date)) %>%
    mutate(chemo_to_surg = as.integer(surgery_date - hes_chemo_date),
           position = case_when(
             chemo_to_surg >  14  ~ "before surgery (neoadjuvant-like)",
             chemo_to_surg >= -14 ~ "around surgery (+/-14d)",
             TRUE                 ~ "after surgery (adjuvant/later)"))
  cat("\nin a surgery pathway (n =", nrow(surg), "):\n")
  surg %>% count(position) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% print()
  
  # curative-RT patients: chemo concurrent with RT is genuine definitive chemoRT
  rt <- hes_only %>% filter(!is.na(rt_date), coalesce(rt_curative, FALSE)) %>%
    mutate(chemo_to_rt = as.integer(rt_date - hes_chemo_date),
           concurrent  = abs(chemo_to_rt) <= hes_chemo_near_rt_days)
  cat("\nin a curative-RT pathway (n =", nrow(rt), "):\n")
  if (nrow(rt) > 0)
    rt %>% count(concurrent) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% print()
  
  # the clean rescue group: chemo where nothing else was recorded
  notx <- hes_only %>% filter(tx_pathway == "No treatment recorded")
  cat("\n'No treatment recorded' rescued by HES chemo (n =", nrow(notx),
      "): reclassify to SACT only / palliative chemo\n")
  cat("(reading: surgery cases mostly 'before surgery' and RT cases mostly\n",
      " 'concurrent' confirms the HES-chemo reclassification is sound)\n")
} else {
  cat("chemo anchor not found - skipping (run 06_derive_sact_rtds.R first)\n")
}

cat("\n11 classification checks complete.\n")