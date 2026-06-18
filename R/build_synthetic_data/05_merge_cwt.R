# =============================================================================
# 05  Merge the CWT records
# -----------------------------------------------------------------------------
# Step 2 of the pipeline applied to the synthetic data. Reads the derived cohort
# and the raw CWT records, runs og_cwt_merge() to attach the decision-to-treat
# node, the waiting-time family, dtt_valid and the audit categories, and saves
# the merged analysis cohort.
#
# Reads : Data/synthetic/og_derived_synthetic.rds        (derived cohort, from 04)
#         Data/synthetic/og_cwt_records_synthetic.rds     (CWT records, from 03)
# Writes: Data/synthetic/og_cohort_synthetic.rds          (merged cohort)
# =============================================================================

library(tidyverse)
library(haven)

# paths, relative to the project root (the .Rproj working directory)
dir_fns <- "R/build_synthetic_data"
dir_syn <- "Data/synthetic"

source(file.path(dir_fns, "01_define_functions.R"))   # og_cwt_merge()

derived <- readRDS(file.path(dir_syn, "og_derived_synthetic.rds"))
cwt     <- readRDS(file.path(dir_syn, "og_cwt_records_synthetic.rds"))

# --- merge -----------------------------------------------------------------
merged <- og_cwt_merge(derived, cwt)

cat("Merged cohort:", nrow(merged), "patients,", ncol(merged), "columns\n")
cat("CWT coverage (with a DTT):",
    round(100 * mean(!is.na(merged$cwt_dtt_date)), 1), "%\n")
cat("dtt_valid TRUE share (non-EMR pathways):",
    round(100 * mean(merged$dtt_valid, na.rm = TRUE), 1), "%\n")

cat("\nAudit Table 4 (stage 1-3): % curative and % any treatment\n")
merged %>% filter(stage_clean %in% c("1","2","3")) %>%
  summarise(pct_curative = round(100*mean(received_curative_tx_audit, na.rm=TRUE)),
            pct_any_tx   = round(100*mean(received_any_tx, na.rm=TRUE))) %>%
  print()

# --- save ------------------------------------------------------------------
saveRDS(merged, file.path(dir_syn, "og_cohort_synthetic.rds"))
write_dta(
  merged %>% mutate(across(where(is.factor), as.character),
                    across(where(is.logical), as.integer)),
  file.path(dir_syn, "og_cohort_synthetic.dta"))

cat("\nSaved merged cohort (", nrow(merged), "patients ).",
    "Next: 06_validate_outputs.R\n")