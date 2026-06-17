# =============================================================================
# OG cancer - HES-only chemo: timing check before reclassification
# -----------------------------------------------------------------------------
# The HES chemo supplement adds 905 net-new chemo patients, but only ~272 are
# currently "No treatment recorded". The rest already carry a surgery or
# curative-RT anchor, so adding a chemo date RECLASSIFIES their pathway. Whether
# that is correct depends on the chemo sitting in a plausible position:
#   - surgery patients: chemo before surgery -> genuine neoadjuvant (reclassify
#     to neoadjuvant chemo); chemo well after -> adjuvant or later-line, less
#     clearly an index treatment.
#   - curative-RT patients: chemo concurrent with RT -> genuine definitive
#     chemoRT; chemo far from RT -> unrelated.
# This profiles those gaps so the reclassification can be trusted or guarded.
# =============================================================================

library(tidyverse)

base_dir <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/"

og <- readRDS(paste0(base_dir, "og_cohort_2015_2022.rds"))
# chemo_anchor must be in scope from the anchor-build block; if running stand-
# alone, source that block first so chemo_anchor (with chemo_source, chemo_date)
# exists. Here we assume it is available.

hes_only <- chemo_anchor %>%
  filter(chemo_source == "hes") %>%
  select(pseudo_patientid, hes_chemo_date, days_dx_to_chemo) %>%
  left_join(og %>% select(pseudo_patientid, tx_pathway, surgery_date,
                          rt_date, rt_curative, stage_clean),
            by = "pseudo_patientid")

# --- surgery patients: where does the HES chemo sit relative to surgery? -----
cat("HES-only chemo patients currently in a surgery pathway:\n")
surg <- hes_only %>%
  filter(!is.na(surgery_date)) %>%
  mutate(chemo_to_surg = as.integer(surgery_date - hes_chemo_date),
         position = case_when(
           chemo_to_surg >  14            ~ "before surgery (neoadjuvant-like)",
           chemo_to_surg >= -14           ~ "around surgery (+/-14d)",
           TRUE                           ~ "after surgery (adjuvant/later)"
         ))
cat("  n =", nrow(surg), "\n")
surg %>% count(position) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% print()
cat("\n  days from HES chemo to surgery (positive = chemo before surgery):\n")
print(summary(surg$chemo_to_surg))

# --- curative-RT patients: is the HES chemo concurrent with the RT? ----------
cat("\nHES-only chemo patients currently in a curative-RT pathway:\n")
rt <- hes_only %>%
  filter(!is.na(rt_date), coalesce(rt_curative, FALSE)) %>%
  mutate(chemo_to_rt = as.integer(rt_date - hes_chemo_date),
         concurrent  = abs(chemo_to_rt) <= 28)
cat("  n =", nrow(rt), "\n")
if (nrow(rt) > 0) {
  rt %>% count(concurrent) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% print()
  cat("\n  days between HES chemo and RT start (concurrent if within 28d):\n")
  print(summary(rt$chemo_to_rt))
}

# --- no-treatment patients: the clean rescue group ---------------------------
cat("\nHES-only chemo patients currently 'No treatment recorded':\n")
notx <- hes_only %>% filter(tx_pathway == "No treatment recorded")
cat("  n =", nrow(notx), " (these reclassify to SACT only / palliative chemo)\n")
cat("  stage split:\n")
notx %>% count(stage_clean) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% print()
cat("  days dx to HES chemo:\n")
print(summary(notx$days_dx_to_chemo))

cat("\nReading:\n",
    " surgery patients mostly 'before surgery' -> genuine neoadjuvant SACT\n",
    "   missed; reclassification to neoadjuvant chemo is sound.\n",
    " surgery patients mostly 'after surgery' -> adjuvant/later-line; consider\n",
    "   whether these should set first_tx_date or only the treatment flag.\n",
    " curative-RT patients mostly 'concurrent' -> genuine definitive chemoRT.\n",
    " no-treatment group is the clean rescue: chemo where none was recorded.\n")
