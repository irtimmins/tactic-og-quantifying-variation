# =============================================================================
# OG cancer - validate the minimal synthetic dataset
# -----------------------------------------------------------------------------
# Print-only checks that the synthetic cohort (a) conforms to the minimal column
# contract, (b) is internally consistent (the audit categories re-derive from
# tx_pathway, dates are ordered), and (c) reproduces the real cohort's pathway
# mix and audit targets within tolerance. Run after 03_og_generate.R.
# =============================================================================

library(tidyverse)
source("01_og_minimal_merge.R")   # og_minimal_cols, tx_pathway_levels

base_dir <- "Data/synthetic/"
A   <- readRDS(paste0(base_dir, "og_ncras_treatment_synthetic.rds"))
coh <- readRDS(paste0(base_dir, "og_cohort_synthetic.rds"))
cwt <- readRDS(paste0(base_dir, "og_cwt_records_synthetic.rds"))
profile <- if (file.exists(paste0(base_dir, "og_profile_for_synthetic.rds")))
  readRDS(paste0(base_dir, "og_profile_for_synthetic.rds")) else NULL

sec <- function(x) cat("\n== ", x, " ==\n")

# -----------------------------------------------------------------------------
sec("1. Conformance")
cat("Table A rows:", nrow(A), " cols:", ncol(A), "\n")
cat("Merged rows: ", nrow(coh), " CWT records:", nrow(cwt), "\n")
missing_req <- setdiff(c("pseudo_patientid","diagmdy","tx_pathway",
                         "first_tx_date","endoscopy_date"), names(A))
cat("Missing merge-critical cols:",
    if (length(missing_req)) paste(missing_req, collapse=", ") else "none", "\n")
cat("Duplicate patient IDs:", sum(duplicated(A$pseudo_patientid)), "\n")
bad_pw <- setdiff(unique(na.omit(A$tx_pathway)), tx_pathway_levels)
cat("Unexpected tx_pathway:",
    if (length(bad_pw)) paste(bad_pw, collapse=", ") else "none", "\n")
bad_stage <- setdiff(unique(na.omit(A$stage_clean)), c("1","2","3"))
cat("Unexpected stage_clean:",
    if (length(bad_stage)) paste(bad_stage, collapse=", ") else "none", "\n")

# -----------------------------------------------------------------------------
sec("2. Internal consistency")
# first_tx_date should be NA exactly for the non-curative pathways
noncurative <- c("Palliative chemo + RT","SACT only","Palliative RT only",
                 "No treatment recorded")
ftx_na_ok <- all(is.na(A$first_tx_date[A$tx_pathway %in% noncurative])) &&
             all(!is.na(A$first_tx_date[!A$tx_pathway %in% noncurative &
                                        A$tx_pathway != "Surgery + other"]))
cat("first_tx_date NA pattern matches pathway intent:", ftx_na_ok, "\n")

# date ordering: endoscopy <= diagnosis-ish, first_tx >= diagnosis
cat("first_tx_date before diagnosis (should be 0):",
    sum(!is.na(coh$first_tx_date) & coh$first_tx_date < coh$diagmdy), "\n")
cat("negative wt_dx_to_dtt (should be ~0):",
    sum(coh$wt_dx_to_dtt < 0, na.rm = TRUE), "\n")

# audit flags must not be identical (the regression we guarded against upstream)
cat("received_any_tx identical to curative (should be FALSE):",
    identical(coh$received_any_tx, coh$received_curative_tx_audit), "\n")

# neoadjuvant patients should anchor on a chemo/RT CWT modality, not surgery
neo <- coh %>%
  filter(tx_pathway %in% c("Surgery + neoadjuvant chemo",
                           "Surgery + neoadjuvant chemoRT",
                           "Surgery + neoadjuvant RT"),
         !is.na(cwt_modality))
if (nrow(neo))
  cat("neoadjuvant patients anchored on chemo/RT modality (02/04/05):",
      round(100 * mean(neo$cwt_modality %in% c("02","04","05")), 1), "%\n")

# -----------------------------------------------------------------------------
sec("3. Pathway mix: synthetic vs profile")
syn_mix <- A %>% count(tx_pathway) %>% mutate(syn = round(100*n/sum(n),1)) %>%
  select(tx_pathway, syn)
if (!is.null(profile$pathway_overall)) {
  real_mix <- profile$pathway_overall %>%
    transmute(tx_pathway = level, real = round(100*prop,1))
  full_join(syn_mix, real_mix, by = "tx_pathway") %>%
    mutate(diff = syn - real) %>% arrange(desc(real)) %>% print(n = 20)
} else {
  print(syn_mix %>% arrange(desc(syn)), n = 20)
  cat("(no profile loaded - showing synthetic mix only)\n")
}

# -----------------------------------------------------------------------------
sec("4. Audit targets: synthetic vs real")
syn_aud <- coh %>% filter(stage_clean %in% c("1","2","3")) %>%
  summarise(curative = round(100*mean(received_curative_tx_audit, na.rm=TRUE)),
            any_tx   = round(100*mean(received_any_tx, na.rm=TRUE)))
cat("synthetic  curative:", syn_aud$curative, "  any-tx:", syn_aud$any_tx, "\n")
if (!is.null(profile$audit_targets))
  cat("real       curative:", round(100*profile$audit_targets$pct_curative),
      "  any-tx:", round(100*profile$audit_targets$pct_any_tx), "\n")

sec("5. Audit Table 4 by subtype (synthetic)")
coh %>% filter(stage_clean %in% c("1","2","3")) %>%
  mutate(subtype = coalesce(cancer_subtype,"Unknown")) %>%
  group_by(subtype) %>%
  summarise(n = n(),
            surgery_only = round(100*mean(tx_modality_audit=="Surgery only", na.rm=TRUE)),
            definitive_chemoRT = round(100*mean(tx_modality_audit=="Definitive chemoRT", na.rm=TRUE)),
            emresd = round(100*mean(tx_modality_audit=="EMR/ESD", na.rm=TRUE)),
            curative = round(100*mean(received_curative_tx_audit, na.rm=TRUE)),
            .groups = "drop") %>% print()

cat("\nValidation complete.\n")
