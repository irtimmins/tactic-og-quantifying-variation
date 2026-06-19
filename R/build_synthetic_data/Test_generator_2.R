# =============================================================================
# Run the synthetic two-step workflow end to end  (convenience driver)
# -----------------------------------------------------------------------------
# Reads the raw synthetic cohort and CWT records that 03 produced, then runs the
# two coding stages in one go:
#   step 1  og_derive_pathway()  raw -> derived (flags, tx_pathway, first_tx_date, tx_trust)
#   step 2  og_cwt_merge()       derived + CWT -> merged analysis cohort
#
# This does the same work as 04_derive_pathway.R + 05_merge_cwt.R, but in a
# single script for quick end-to-end runs and inspection. The functions live in
# 01_define_functions.R. Run 03_generate_raw_synthetic_datasets.R first.
# =============================================================================

library(tidyverse)

dir_fns <- "R/build_synthetic_data"
dir_syn <- "Data/synthetic"

source(file.path(dir_fns, "01_define_functions.R"))   # og_derive_pathway(), og_cwt_merge()

# --- read the inputs --------------------------------------------------------
# The .rds files preserve R Date types, so the date arithmetic in the two
# functions works directly. (The .dta exports store dates as Stata daily values;
# use those for Stata, not here.)
raw <- readRDS(file.path(dir_syn, "og_ncras_treatment_synthetic.rds"))
cwt <- readRDS(file.path(dir_syn, "og_cwt_records_synthetic.rds"))
cat("raw cohort:", nrow(raw), "patients,", ncol(raw), "columns\n")
cat("CWT records:", nrow(cwt), "rows\n")

# --- step 1: derive pathway, first treatment date and treatment trust -------
derived <- og_derive_pathway(raw)
cat("\nderived pathway mix:\n")
derived %>% count(tx_pathway) %>% mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(desc(n)) %>% print(n = 20)

# --- step 2: merge the CWT records onto the derived cohort ------------------
merged <- og_cwt_merge(derived, cwt)
cat("\nmerged cohort:", nrow(merged), "patients,", ncol(merged), "columns\n")
cat("CWT coverage (with a DTT):",
    round(100 * mean(!is.na(merged$cwt_dtt_date)), 1), "%\n")

cat("\naudit Table 4 (stage 1-3): % curative and % any treatment\n")
merged %>% filter(stage_clean %in% c("1", "2", "3")) %>%
  summarise(pct_curative = round(100 * mean(received_curative_tx_audit, na.rm = TRUE)),
            pct_any_tx   = round(100 * mean(received_any_tx, na.rm = TRUE))) %>%
  print()

# --- optional: save the two products under distinct names -------------------
# Named *_from_workflow so they do not clash with the per-step outputs that
# 04 and 05 write (og_derived_synthetic.rds / og_cohort_synthetic.rds).
saveRDS(derived, file.path(dir_syn, "og_derived_from_workflow.rds"))
saveRDS(merged,  file.path(dir_syn, "og_cohort_from_workflow.rds"))
cat("\nSaved og_derived_from_workflow.rds and og_cohort_from_workflow.rds\n")

