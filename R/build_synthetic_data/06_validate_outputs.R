# =============================================================================
# 06  Validate the outputs
# -----------------------------------------------------------------------------
# Print-only checks that the synthetic data (a) conforms to the raw and derived
# contracts, (b) is internally consistent - the pathway and trust re-derive from
# the raw flags/dates, the audit categories follow tx_pathway, dates are ordered
# - and (c) reproduces the real cohort's pathway mix and audit targets. Run last,
# after 03 (generate), 04 (derive) and 05 (merge).
#
# Reads: the raw, derived and merged cohorts + the CWT records from Data/synthetic
# =============================================================================

library(tidyverse)

# paths, relative to the project root (the .Rproj working directory)
dir_fns <- "R/build_synthetic_data"
dir_syn <- "Data/synthetic"

source(file.path(dir_fns, "01_define_functions.R"))   # og_derive_pathway, contracts

raw <- readRDS(file.path(dir_syn, "og_ncras_treatment_synthetic.rds"))  # raw (03)
der <- readRDS(file.path(dir_syn, "og_derived_synthetic.rds"))          # derived (04)
coh <- readRDS(file.path(dir_syn, "og_cohort_synthetic.rds"))           # merged (05)
cwt <- readRDS(file.path(dir_syn, "og_cwt_records_synthetic.rds"))
profile <- if (file.exists(file.path(dir_syn, "og_profile_for_synthetic.rds")))
  readRDS(file.path(dir_syn, "og_profile_for_synthetic.rds")) else NULL

sec <- function(x) cat("\n== ", x, " ==\n")

# -----------------------------------------------------------------------------
sec("1. Conformance")
cat("Raw rows:", nrow(raw), " cols:", ncol(raw), "\n")
cat("Derived rows:", nrow(der), " | Merged rows:", nrow(coh),
    " | CWT records:", nrow(cwt), "\n")
missing_raw <- setdiff(c("pseudo_patientid","diagmdy","endoscopy_date",
                         "surgery_date","sact_date","rt_date","emresd_date",
                         "curative_surgery","rt_curative",
                         "surgery_provider","rt_provider"), names(raw))
cat("Missing raw-contract cols:",
    if (length(missing_raw)) paste(missing_raw, collapse=", ") else "none", "\n")
cat("Raw carries a pre-built pathway (should be FALSE):",
    "tx_pathway" %in% names(raw), "\n")
cat("Duplicate patient IDs:", sum(duplicated(raw$pseudo_patientid)), "\n")
bad_pw <- setdiff(unique(na.omit(der$tx_pathway)), tx_pathway_levels)
cat("Unexpected tx_pathway:",
    if (length(bad_pw)) paste(bad_pw, collapse=", ") else "none", "\n")

# -----------------------------------------------------------------------------
sec("2. Derivation reproduces (re-run stage 1 on the raw inputs)")
# the headline transparency check: deriving the pathway from the raw flags/dates
# again must reproduce the saved derived cohort exactly
re <- og_derive_pathway(raw %>% select(any_of(og_raw_cols)))
cat("re-derived tx_pathway matches saved:",
    round(100 * mean(re$tx_pathway == der$tx_pathway), 2), "%\n")
cat("re-derived first_tx_date matches saved:",
    round(100 * mean(re$first_tx_date == der$first_tx_date |
                       (is.na(re$first_tx_date) & is.na(der$first_tx_date))), 2), "%\n")
cat("re-derived tx_trust matches saved:",
    round(100 * mean(re$tx_trust == der$tx_trust |
                       (is.na(re$tx_trust) & is.na(der$tx_trust))), 2), "%\n")

# tx_trust comes from the right provider: surgical pathways from surgery,
# RT-anchored from RT, NA for non-curative
trust_ok_surg <- der %>%
  filter(tx_pathway %in% c("Surgery only","Surgery + neoadjuvant chemo",
                           "Surgery + adjuvant chemo","EMR/ESD only",
                           "EMR/ESD then surgery")) %>%
  summarise(ok = mean(tx_trust == substr(surgery_provider,1,3), na.rm = TRUE)) %>%
  pull(ok)
trust_ok_rt <- der %>%
  filter(tx_pathway %in% c("Definitive chemoRT","Curative RT only",
                           "Surgery + neoadjuvant RT","Surgery + neoadjuvant chemoRT")) %>%
  summarise(ok = mean(tx_trust == substr(rt_provider,1,3), na.rm = TRUE)) %>%
  pull(ok)
cat("tx_trust from surgery provider on surgical pathways:",
    round(100*trust_ok_surg,1), "%\n")
cat("tx_trust from RT provider on RT-anchored pathways:",
    round(100*trust_ok_rt,1), "%\n")
cat("tx_trust missing on no-treatment (should be 100):",
    round(100*mean(is.na(der$tx_trust[der$tx_pathway == "No treatment recorded"])),1), "%\n")
bad_pw <- setdiff(unique(na.omit(A$tx_pathway)), tx_pathway_levels)
cat("Unexpected tx_pathway:",
    if (length(bad_pw)) paste(bad_pw, collapse=", ") else "none", "\n")
bad_stage <- setdiff(unique(na.omit(A$stage_clean)), c("1","2","3"))
cat("Unexpected stage_clean:",
    if (length(bad_stage)) paste(bad_stage, collapse=", ") else "none", "\n")

# -----------------------------------------------------------------------------
sec("3. Internal consistency")
# first_tx_date should be NA exactly for the non-curative pathways
noncurative <- c("Palliative chemo + RT","SACT only","Palliative RT only",
                 "No treatment recorded")
ftx_na_ok <- all(is.na(der$first_tx_date[der$tx_pathway %in% noncurative])) &&
  all(!is.na(der$first_tx_date[!der$tx_pathway %in% noncurative &
                                 der$tx_pathway != "Surgery + other"]))
cat("first_tx_date NA pattern matches pathway intent:", ftx_na_ok, "\n")

# date ordering: first_tx >= diagnosis, no negative dx -> dtt
cat("first_tx_date before diagnosis (should be 0):",
    sum(!is.na(coh$first_tx_date) & coh$first_tx_date < coh$diagmdy), "\n")
cat("negative wt_dx_to_dtt (should be ~0):",
    sum(coh$wt_dx_to_dtt < 0, na.rm = TRUE), "\n")

# audit flags must not be identical (the regression guarded against upstream)
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
sec("4. Pathway mix: synthetic vs profile")
syn_mix <- der %>% count(tx_pathway) %>% mutate(syn = round(100*n/sum(n),1)) %>%
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
sec("5. Audit targets: synthetic vs real")
syn_aud <- coh %>% filter(stage_clean %in% c("1","2","3")) %>%
  summarise(curative = round(100*mean(received_curative_tx_audit, na.rm=TRUE)),
            any_tx   = round(100*mean(received_any_tx, na.rm=TRUE)))
cat("synthetic  curative:", syn_aud$curative, "  any-tx:", syn_aud$any_tx, "\n")
if (!is.null(profile$audit_targets))
  cat("real       curative:", round(100*profile$audit_targets$pct_curative),
      "  any-tx:", round(100*profile$audit_targets$pct_any_tx), "\n")

sec("6. Audit Table 4 by subtype (synthetic)")
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