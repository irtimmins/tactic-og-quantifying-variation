# =============================================================================
# OG cancer - neoadjuvant clock-stop check
# -----------------------------------------------------------------------------
# CWT records the FIRST definitive treatment, and for a neoadjuvant-chemo-then-
# surgery patient that is the chemotherapy (guidance 3.9.1), i.e. the EARLIER
# event. The section-G crosstab showed off-diagonal cells where the CWT anchor's
# modality is "surgery" but the pathway is a neoadjuvant chemo/RT pathway. This
# script asks, for those patients, whether cwt_treat_date actually sits on the
# earlier neoadjuvant treatment (CWT behaving per guidance, merge captured it) or
# on the later surgery (CWT/merge anchored on the wrong event).
#
# It compares cwt_treat_date against sact_date, rt_date and surgery_date for the
# neoadjuvant pathways, and reports which event the CWT date is closest to.
# =============================================================================

library(tidyverse)

base_dir <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/"
og <- readRDS(paste0(base_dir, "og_cohort_cwt_2015_2022.rds"))

neoadj_pathways <- c("Surgery + neoadjuvant chemo",
                     "Surgery + neoadjuvant chemoRT",
                     "Surgery + neoadjuvant RT")

neo <- og %>%
  filter(tx_pathway %in% neoadj_pathways, !is.na(cwt_treat_date)) %>%
  mutate(
    # earliest neoadjuvant systemic/RT event the pathway is built on
    neoadj_date = pmin(sact_date, rt_date, na.rm = TRUE),
    d_to_neoadj = as.integer(cwt_treat_date - neoadj_date),
    d_to_surg   = as.integer(cwt_treat_date - surgery_date),
    # which event is the CWT treatment date closest to?
    closest = case_when(
      is.na(neoadj_date) & is.na(surgery_date) ~ "neither date",
      is.na(surgery_date)                      ~ "neoadjuvant",
      is.na(neoadj_date)                       ~ "surgery",
      abs(d_to_neoadj) <= abs(d_to_surg)       ~ "neoadjuvant",
      TRUE                                     ~ "surgery"
    )
  )

cat("Neoadjuvant-pathway patients with a CWT treat date:", nrow(neo), "\n\n")

cat("Which event is cwt_treat_date closest to, by pathway:\n")
neo %>% count(tx_pathway, closest) %>%
  group_by(tx_pathway) %>% mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(tx_pathway, desc(n)) %>% print(n = 30)

cat("\nDays from CWT treat date to the neoadjuvant event (negative = CWT earlier):\n")
neo %>% filter(!is.na(d_to_neoadj)) %>%
  group_by(tx_pathway) %>%
  summarise(n = n(),
            median = median(d_to_neoadj),
            p25 = quantile(d_to_neoadj, .25),
            p75 = quantile(d_to_neoadj, .75),
            pct_within_14 = round(100 * mean(abs(d_to_neoadj) <= 14), 1),
            .groups = "drop") %>% print()

cat("\nDays from CWT treat date to surgery (negative = CWT earlier than surgery):\n")
neo %>% filter(!is.na(d_to_surg)) %>%
  group_by(tx_pathway) %>%
  summarise(n = n(),
            median = median(d_to_surg),
            p25 = quantile(d_to_surg, .25),
            p75 = quantile(d_to_surg, .75),
            .groups = "drop") %>% print()

cat("\nFor reference, CWT modality among these neoadjuvant patients:\n")
neo %>% count(tx_pathway, cwt_modality) %>%
  group_by(tx_pathway) %>% mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(tx_pathway, desc(n)) %>% print(n = 40)

cat("\nReading:\n",
    " 'closest = neoadjuvant' dominant + CWT date near the SACT/RT date + CWT\n",
    "   modality 02/04/05 -> CWT anchors on the earlier neoadjuvant treatment as\n",
    "   the guidance intends, and the merge captured it; the section-G surgery\n",
    "   off-diagonal is then a small minority, not the rule.\n",
    " 'closest = surgery' dominant + CWT modality 01/23/24 -> for these patients\n",
    "   CWT (or the merge selection) anchored on the later surgery; worth checking\n",
    "   whether a chemo CWT row existed but fell outside the merge's DTT window.\n")