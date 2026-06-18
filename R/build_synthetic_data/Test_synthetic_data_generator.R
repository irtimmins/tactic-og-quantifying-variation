# =============================================================================
# OG cancer - run the two-step workflow in R
# -----------------------------------------------------------------------------
# Reads the raw registry+treatment cohort and the raw CWT records, then runs the
# two coding stages as functions:
#   step 1  og_derive_pathway()  raw  -> derived (flags, tx_pathway, first_tx_date, tx_trust)
#   step 2  og_cwt_merge()       derived + CWT -> merged analysis cohort
# Both functions live in 01_og_minimal_merge.R. This script just wires them
# together so the pipeline can be run and inspected end to end.
# =============================================================================

library(tidyverse)

source("R/build_synthetic_data/01_og_minimal_merge.R")   # og_derive_pathway(), og_cwt_merge()

base_dir <- "Data/synthetic/"

# --- read the inputs --------------------------------------------------------
# The .rds files preserve R Date types, so the date arithmetic in the two
# functions works directly. (The .dta exports store dates as Stata daily values;
# use those for Stata, not here.)
raw <- readRDS(paste0(base_dir, "og_ncras_treatment_synthetic.rds"))
cwt <- readRDS(paste0(base_dir, "og_cwt_records_synthetic.rds"))

cat("raw cohort:", nrow(raw), "patients,", ncol(raw), "columns\n")
cat("CWT records:", nrow(cwt), "rows\n")

# --- step 1: derive the pathway, first treatment date and treatment trust ---
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

# --- optional: save the two products ----------------------------------------
saveRDS(derived, paste0(base_dir, "og_derived_from_workflow.rds"))
saveRDS(merged,  paste0(base_dir, "og_cohort_from_workflow.rds"))
cat("\nSaved og_derived_from_workflow.rds and og_cohort_from_workflow.rds\n")
