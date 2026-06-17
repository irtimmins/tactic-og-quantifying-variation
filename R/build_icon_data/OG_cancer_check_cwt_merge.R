# =============================================================================
# OG cancer - checks on the CWT merge and audit reproduction
# -----------------------------------------------------------------------------
# Run after script 4. Reads og_cohort_cwt_2015_2022.rds and works through the
# things worth confirming before trusting the audit numbers:
#   1. what stage_clean actually contains (is the Table 4 denominator real?)
#   2. the received_any_tx bug (does it collapse onto the curative flag?)
#   3. whether the "No treatment recorded" group is genuinely untreated
#   4. modality-to-pathway matching quality
#   5. the corrected any-treatment flag and a re-run of Table 4
#
# Nothing is saved; this only prints. Apply the fixes in script 4 once the
# numbers below confirm them.
# =============================================================================

library(tidyverse)

base_dir <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/"
tx_window_days <- 270L

og <- readRDS(paste0(base_dir, "og_cohort_cwt_2015_2022.rds"))

line <- function(x) cat("\n========== ", x, " ==========\n")

# -----------------------------------------------------------------------------
# 1. What does stage_clean hold? Is the stage 1-3 filter doing anything?
# -----------------------------------------------------------------------------
line("1. stage_clean contents")
og %>% count(stage_clean) %>% arrange(desc(n)) %>% print(n = 30)

cat("\nrows matching the literal c('1','2','3') filter:",
    sum(og$stage_clean %in% c("1", "2", "3")), "of", nrow(og), "\n")
cat("If these differ, the Table 4 stage filter is a no-op and needs to match\n",
    "the real coding (e.g. 'Stage 1' or numeric).\n")

# -----------------------------------------------------------------------------
# 2. The received_any_tx bug: does it collapse onto the curative flag?
# -----------------------------------------------------------------------------
line("2. received_any_tx vs received_curative_tx_audit")

# rebuild the flags exactly as script 4 currently does (the suspect version)
chk <- og %>%
  mutate(
    any_tx_OLD = tx_pathway != "No treatment recorded" &
      !is.na(first_tx_date) & wt_dx_to_tx <= tx_window_days,
    cur_tx     = tx_intent_audit == "Curative" &
      !is.na(first_tx_date) & wt_dx_to_tx <= tx_window_days
  )

cat("identical(any_tx_OLD, cur_tx):",
    identical(chk$any_tx_OLD, chk$cur_tx), "\n")
cat("(TRUE confirms the bug: 'any treatment' is really 'curative treatment')\n\n")

# the mechanism: first_tx_date is the curative clock-stop, NA for palliative
cat("first_tx_date missingness by treatment intent:\n")
og %>% group_by(tx_intent_audit) %>%
  summarise(n = n(),
            n_first_tx_na = sum(is.na(first_tx_date)),
            pct_na = round(100 * mean(is.na(first_tx_date)), 1),
            .groups = "drop") %>% print()
cat("\nNon-curative rows have first_tx_date NA, so the nine-month gate drops\n",
    "every palliative patient from 'any treatment'.\n")

# -----------------------------------------------------------------------------
# 3. Is "No treatment recorded" genuinely untreated, or leaking treated cases?
# -----------------------------------------------------------------------------
line("3. No treatment recorded vs presence of a CWT anchor")

og %>%
  mutate(no_tx_pathway = tx_pathway == "No treatment recorded",
         has_cwt_treat = !is.na(cwt_treat_date)) %>%
  count(no_tx_pathway, has_cwt_treat) %>%
  group_by(no_tx_pathway) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

cat("\nOf the 'No treatment recorded' group, the share with a CWT treatment\n",
    "record is the leakage estimate: if non-trivial, treated patients are\n",
    "landing in the no-treatment bucket (a script-2 pathway issue, not here).\n")

# what modalities do those supposedly-untreated-but-CWT-present patients have?
cat("\nCWT modality among 'No treatment recorded' patients who DO have a CWT record:\n")
og %>%
  filter(tx_pathway == "No treatment recorded", !is.na(cwt_treat_date)) %>%
  count(cwt_modality, sort = TRUE) %>% print(n = 30)

# are these mostly palliative codes (expected) or curative ones (a problem)?
cat("\n('07/08/09' palliative and '97' other are expected here; a lot of\n",
    "'23/24/01' surgery or '02/04/05' would indicate genuine leakage.)\n")

# -----------------------------------------------------------------------------
# 4. Matching quality: how often did the modality match the pathway?
# -----------------------------------------------------------------------------
line("4. CWT anchor coverage and treatment-date completeness by pathway")

og %>%
  group_by(tx_pathway) %>%
  summarise(
    n               = n(),
    pct_with_dtt    = round(100 * mean(!is.na(cwt_dtt_date)), 1),
    pct_with_cwttx  = round(100 * mean(!is.na(cwt_treat_date)), 1),
    pct_with_firsttx= round(100 * mean(!is.na(first_tx_date)), 1),
    .groups = "drop"
  ) %>% arrange(desc(n)) %>% print(n = 20)

# -----------------------------------------------------------------------------
# 5. Corrected any-treatment flag, and Table 4 re-run with it
# -----------------------------------------------------------------------------
line("5. Corrected received_any_tx and Table 4")

og_fix <- og %>%
  mutate(
    days_dx_to_cwttx = as.integer(cwt_treat_date - diagmdy),
    received_any_tx_FIX = tx_pathway != "No treatment recorded" &
      ( (!is.na(first_tx_date)  & wt_dx_to_tx <= tx_window_days) |
          (!is.na(cwt_treat_date) & days_dx_to_cwttx <= tx_window_days &
             days_dx_to_cwttx >= 0) )
  )

cat("any-treatment rate, OLD vs FIXED (whole cohort):\n")
og_fix %>% summarise(
  old = round(100 * mean(tx_pathway != "No treatment recorded" &
                           !is.na(first_tx_date) & wt_dx_to_tx <= tx_window_days), 1),
  fixed = round(100 * mean(received_any_tx_FIX), 1)
) %>% print()
cat("(audit England any treatment ~76%; expect the fixed figure to rise toward it)\n")

# Table 4 re-run, stage 1-3, using a stage filter that adapts to the coding
stage_levels <- og_fix %>% count(stage_clean) %>% pull(stage_clean)
use_stage <- intersect(c("1", "2", "3"), stage_levels)
if (length(use_stage) == 0) {
  cat("\nstage_clean is not '1'/'2'/'3' - using all rows (cohort is stage 1-3).\n")
  aud <- og_fix
} else {
  aud <- og_fix %>% filter(stage_clean %in% use_stage)
}

cat("\n--- Table 4 (corrected any-treatment), by subtype ---\n")
t4 <- aud %>%
  mutate(subtype = coalesce(cancer_subtype, "Unknown")) %>%
  group_by(subtype) %>%
  summarise(
    n_people             = n(),
    pct_surgery_only     = round(100 * mean(tx_modality_audit == "Surgery only", na.rm = TRUE), 0),
    pct_surgery_plus     = round(100 * mean(tx_modality_audit == "Surgery plus SACT/RT", na.rm = TRUE), 0),
    pct_definitive_chemRT= round(100 * mean(tx_modality_audit == "Definitive chemoRT", na.rm = TRUE), 0),
    pct_curative_rt_only = round(100 * mean(tx_modality_audit == "Curative RT only", na.rm = TRUE), 0),
    pct_emresd           = round(100 * mean(tx_modality_audit == "EMR/ESD", na.rm = TRUE), 0),
    pct_curative_overall = round(100 * mean(received_curative_tx_audit, na.rm = TRUE), 0),
    pct_any_treatment    = round(100 * mean(received_any_tx_FIX, na.rm = TRUE), 0),
    .groups = "drop"
  )
t4_all <- aud %>%
  summarise(
    subtype = "All", n_people = n(),
    pct_surgery_only     = round(100 * mean(tx_modality_audit == "Surgery only", na.rm = TRUE), 0),
    pct_surgery_plus     = round(100 * mean(tx_modality_audit == "Surgery plus SACT/RT", na.rm = TRUE), 0),
    pct_definitive_chemRT= round(100 * mean(tx_modality_audit == "Definitive chemoRT", na.rm = TRUE), 0),
    pct_curative_rt_only = round(100 * mean(tx_modality_audit == "Curative RT only", na.rm = TRUE), 0),
    pct_emresd           = round(100 * mean(tx_modality_audit == "EMR/ESD", na.rm = TRUE), 0),
    pct_curative_overall = round(100 * mean(received_curative_tx_audit, na.rm = TRUE), 0),
    pct_any_treatment    = round(100 * mean(received_any_tx_FIX, na.rm = TRUE), 0)
  )
bind_rows(t4_all, t4) %>% print(width = Inf)

# -----------------------------------------------------------------------------
# 6. Sense check: any-treatment should never be below curative
# -----------------------------------------------------------------------------
line("6. consistency: any-treatment >= curative in every subtype")
bind_rows(t4_all, t4) %>%
  mutate(ok = pct_any_treatment >= pct_curative_overall) %>%
  select(subtype, pct_curative_overall, pct_any_treatment, ok) %>%
  print()
cat("(all ok should be TRUE; any FALSE means the flags are still inconsistent)\n")


leak <- og %>%
  filter(tx_pathway == "No treatment recorded",
         !is.na(cwt_treat_date),
         cwt_modality %in% c("01","23","24","02","04","05")) %>%
  mutate(dx_to_cwttx = as.integer(cwt_treat_date - diagmdy))
summary(leak$dx_to_cwttx)        # are these within 270 days of diagnosis?
leak %>% count(cwt_modality, sort = TRUE)

