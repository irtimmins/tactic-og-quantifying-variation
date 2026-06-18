# =============================================================================
# 04  Derive the treatment pathway
# -----------------------------------------------------------------------------
# Step 1 of the pipeline applied to the synthetic data. Reads the raw cohort,
# runs og_derive_pathway() to build the treatment flags, the sequencing flags,
# tx_pathway, first_tx_date and tx_trust, and saves the derived cohort.
#
# Reads : Data/synthetic/og_ncras_treatment_synthetic.rds   (raw cohort, from 03)
# Writes: Data/synthetic/og_derived_synthetic.rds            (derived cohort)
# =============================================================================

library(tidyverse)
library(haven)

# paths, relative to the project root (the .Rproj working directory)
dir_fns <- "R/build_synthetic_data"
dir_syn <- "Data/synthetic"

source(file.path(dir_fns, "01_define_functions.R"))   # og_derive_pathway()

raw <- readRDS(file.path(dir_syn, "og_ncras_treatment_synthetic.rds"))

# --- derive ----------------------------------------------------------------
derived <- og_derive_pathway(raw)

cat("Derived pathway mix:\n")
derived %>% count(tx_pathway) %>% mutate(pct = round(100*n/sum(n),1)) %>%
  arrange(desc(n)) %>% print(n = 20)

# --- QC: does the derivation recover the generator's intended pathway? ------
# (only available for synthetic data, where 03 saved the intended labels)
qc_path <- file.path(dir_syn, "og_intended_pathway_qc.rds")
if (file.exists(qc_path)) {
  qc <- readRDS(qc_path)
  chk <- derived %>% select(pseudo_patientid, tx_pathway) %>%
    inner_join(qc, by = "pseudo_patientid")
  cat("\nDerived vs intended pathway agreement:",
      round(100 * mean(chk$tx_pathway == chk$tx_pathway_intended), 1), "%\n")
}

# --- save ------------------------------------------------------------------
saveRDS(derived, file.path(dir_syn, "og_derived_synthetic.rds"))
write_dta(
  derived %>% mutate(across(where(is.factor), as.character),
                     across(where(is.logical), as.integer)),
  file.path(dir_syn, "og_derived_synthetic.dta"))

cat("\nSaved derived cohort (", nrow(derived), "patients ).",
    "Next: 05_merge_cwt.R\n")