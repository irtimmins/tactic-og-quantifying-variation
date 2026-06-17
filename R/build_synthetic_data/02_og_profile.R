# =============================================================================
# OG cancer - profiling for minimal synthetic data
# -----------------------------------------------------------------------------
# Runs on the secure server against the real post-merge cohort. Extracts the
# aggregate, disclosure-controlled distributions the generator needs to build a
# synthetic registry+treatment cohort (Table A) and a synthetic CWT records
# table (Table B) that reproduce the OG audit pathways and the CWT merge.
#
# Pathway-aware by design: OG treatment is multi-modality, so the profile is
# keyed on tx_pathway (and stage x subtype), not a single surgery date. That is
# what lets the generator reproduce the audit categorisation realistically.
#
# Disclosure control: counts rounded to nearest 5, cells < 10 suppressed,
# intervals reported as quantiles only. Output is aggregate; still send through
# output checking before release.
#
# Produces: og_profile_for_synthetic.rds   (distributions)
#           og_minimal_spec.rds            (column manifest + merge constants)
# Input:    og_cohort_cwt_2015_2022.rds    (post-merge analysis cohort)
#           the partitioned CWT dataset    (optional, for per-record stats)
# =============================================================================

library(tidyverse)
library(arrow)

base_dir <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/"
save_dir <- "Data/synthetic/"
cwt_path <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/CWT/11_CWT_data_partitioned"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
SDC_MIN  <- 10L
og_icd   <- c("C150","C151","C152","C153","C154","C155","C158","C159","C15",
              "C160","C161","C162","C163","C164","C165","C166","C168","C169","C16")

cohort <- readRDS(paste0(base_dir, "og_cohort_cwt_2015_2022.rds"))

# -----------------------------------------------------------------------------
# Disclosure-safe helpers
# -----------------------------------------------------------------------------
r5  <- function(x) round(x / 5) * 5
sup <- function(n) ifelse(n < SDC_MIN, NA_integer_, as.integer(r5(n)))

cat_marg <- function(df, var) {
  if (!var %in% names(df)) return(NULL)
  df %>% count(.data[[var]], name = "n") %>%
    filter(!is.na(.data[[var]])) %>%
    mutate(level = as.character(.data[[var]]),
           n_safe = sup(n), prop = round(n / sum(n), 4)) %>%
    select(level, n_safe, prop)
}

q_sum <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NULL)
  tibble(n = length(x), mean = round(mean(x), 1), sd = round(sd(x), 1),
         p05 = quantile(x, .05), p25 = quantile(x, .25), p50 = quantile(x, .50),
         p75 = quantile(x, .75), p95 = quantile(x, .95))
}

profile <- list()

# =============================================================================
# A. Scale and temporal structure
# =============================================================================
profile$n_patients <- r5(nrow(cohort))
profile$by_year    <- cat_marg(cohort, "ydiag")

# =============================================================================
# B. Patient / tumour marginals (drive covariate sampling)
# =============================================================================
cohort2 <- cohort %>%
  mutate(age_grp = cut(agediag,
                       breaks = c(-Inf,50,55,60,65,70,75,80,85,Inf),
                       labels = c("<50","50-54","55-59","60-64","65-69",
                                  "70-74","75-79","80-84","85+"),
                       right = FALSE))

profile$marginals <- list(
  sex            = cat_marg(cohort,  "sex"),
  age_grp        = cat_marg(cohort2, "age_grp"),
  ethnicity      = cat_marg(cohort,  "ethnicity_group_broad"),
  imd_quintile   = cat_marg(cohort,  "NHSE_reversed_imd_quintile_lsoas"),
  site_grp       = cat_marg(cohort,  "tumour_site_grp"),
  subtype        = cat_marg(cohort,  "cancer_subtype"),
  stage          = cat_marg(cohort,  "stage_clean"),
  route_combined = cat_marg(cohort,  "route_combined"),
  cci_group      = cat_marg(cohort,  "cci_group")
)

# =============================================================================
# C. Pathway structure - the core of the OG model
#    tx_pathway conditional on stage x subtype, so the generator reproduces the
#    audit modality mix (surgery / neoadjuvant / definitive chemoRT / palliative
#    / none) with the right clinical dependence.
# =============================================================================
profile$pathway_overall <- cat_marg(cohort, "tx_pathway")

profile$pathway_by_stage_subtype <- cohort %>%
  mutate(subtype = coalesce(cancer_subtype, "Unknown")) %>%
  count(stage_clean, subtype, tx_pathway, name = "n") %>%
  group_by(stage_clean, subtype) %>%
  mutate(prop = round(n / sum(n), 4), n_safe = sup(n)) %>%
  ungroup() %>%
  select(stage_clean, subtype, tx_pathway, n_safe, prop)

# chemo provenance split (SACT vs HES-supplemented), among chemo-bearing pathways
if ("chemo_source" %in% names(cohort)) {
  profile$chemo_source <- cohort %>%
    filter(!is.na(chemo_source)) %>%
    cat_marg("chemo_source")
}

# audit summary the generator should reproduce (targets for the validator)
profile$audit_targets <- cohort %>%
  filter(stage_clean %in% c("1","2","3")) %>%
  summarise(
    pct_curative   = round(mean(received_curative_tx_audit, na.rm = TRUE), 4),
    pct_any_tx     = round(mean(received_any_tx, na.rm = TRUE), 4)
  )

# =============================================================================
# D. Treatment-timing intervals, per pathway where it matters
#    The generator places anchor dates off diagnosis using these.
# =============================================================================
profile$intervals_overall <- list(
  days_dx_to_endo   = q_sum(-as.integer(cohort$endoscopy_date - cohort$diagmdy)),
  wt_dx_to_dtt      = q_sum(cohort$wt_dx_to_dtt[cohort$wt_dx_to_dtt >= 0]),
  wt_dtt_to_tx      = q_sum(cohort$wt_dtt_to_tx[cohort$wt_dtt_to_tx >= 0]),
  wt_dx_to_tx       = q_sum(cohort$wt_dx_to_tx[cohort$wt_dx_to_tx >= 0]),
  surv_from_dx_days = q_sum(as.integer(cohort$finmdy - cohort$diagmdy))
)

# dx -> first_tx_date by pathway (curative pathways only have a first_tx_date)
profile$dx_to_tx_by_pathway <- cohort %>%
  filter(!is.na(first_tx_date)) %>%
  mutate(d = as.integer(first_tx_date - diagmdy)) %>%
  filter(d >= 0) %>%
  group_by(tx_pathway) %>%
  filter(n() >= SDC_MIN) %>%
  summarise(n = n(), mean = round(mean(d),1), sd = round(sd(d),1),
            p25 = quantile(d,.25), p50 = quantile(d,.50), p75 = quantile(d,.75),
            .groups = "drop")

# neoadjuvant chemo/RT -> surgery gap (defines the two-stage sequence)
profile$neoadj_to_surg <- cohort %>%
  filter(tx_pathway %in% c("Surgery + neoadjuvant chemo",
                           "Surgery + neoadjuvant chemoRT",
                           "Surgery + neoadjuvant RT"),
         !is.na(surgery_date), !is.na(first_tx_date)) %>%
  mutate(d = as.integer(surgery_date - first_tx_date)) %>%
  pull(d) %>% q_sum()

# =============================================================================
# E. Trust / hospital volume structure
# =============================================================================
vol <- function(v) {
  v <- v[!is.na(v) & v != ""]
  if (!length(v)) return(NULL)
  sizes <- as.integer(table(v))
  tibble(n_distinct = length(sizes),
         vol_p25 = quantile(sizes, .25), vol_p50 = quantile(sizes, .50),
         vol_p75 = quantile(sizes, .75), vol_p90 = quantile(sizes, .90))
}
profile$volume <- list(
  diag_trust = vol(cohort$diag_trust),
  diag_hosp  = vol(cohort$diag_hosp[grepl("^R[A-Z0-9]{4}$", cohort$diag_hosp)])
)

# between-hospital SD of the mean DTT wait (random-intercept signal)
hosp_sd <- function(num, hosp) {
  d <- tibble(num = num, hosp = hosp) %>%
    filter(!is.na(num), num >= 0, grepl("^R[A-Z0-9]{4}$", hosp)) %>%
    group_by(hosp) %>% filter(n() >= SDC_MIN) %>%
    summarise(m = mean(num), .groups = "drop")
  if (nrow(d) < 2) return(NA_real_)
  round(sd(d$m), 2)
}
profile$between_hosp_sd <- list(
  wt_dx_to_dtt = hosp_sd(cohort$wt_dx_to_dtt, cohort$diag_hosp),
  wt_dtt_to_tx = hosp_sd(cohort$wt_dtt_to_tx, cohort$diag_hosp)
)

# =============================================================================
# F. CWT per-record structure and merge-glue distributions
# =============================================================================
if (dir.exists(cwt_path)) {
  ids <- cohort$pseudo_patientid
  cwt <- open_dataset(cwt_path) %>%
    filter(site_icd10 %in% og_icd) %>% collect() %>%
    mutate(pseudo_patientid = as.character(pseudo_patientid)) %>%
    filter(pseudo_patientid %in% ids)

  profile$cwt_records_per_patient <- cwt %>%
    count(pseudo_patientid, name = "k") %>% count(k, name = "n_pat") %>%
    mutate(n_safe = sup(n_pat), prop = round(n_pat / sum(n_pat), 4)) %>%
    select(records = k, n_safe, prop)
  profile$cwt_modality_overall <- cat_marg(cwt, "modality")
  profile$cwt_coverage <- tibble(pct_any_cwt = round(mean(ids %in% cwt$pseudo_patientid), 4))
} else {
  profile$cwt_records_per_patient <- NULL
  profile$cwt_modality_overall    <- NULL
  profile$cwt_coverage            <- NULL
}

# modality of the ANCHORED CWT record, by pathway - lets the generator emit a
# CWT modality consistent with each synthetic patient's pathway
profile$cwt_modality_by_pathway <- cohort %>%
  filter(!is.na(cwt_modality)) %>%
  count(tx_pathway, cwt_modality, name = "n") %>%
  group_by(tx_pathway) %>%
  mutate(prop = round(n / sum(n), 4), n_safe = sup(n)) %>%
  ungroup() %>%
  select(tx_pathway, cwt_modality, n_safe, prop)

# DTT completeness and agreement glue (signed offsets reproduce linkage)
profile$cwt_completeness <- tibble(
  pct_dtt = round(mean(!is.na(cohort$cwt_dtt_date)), 4),
  pct_mdt = round(mean(!is.na(cohort$cwt_mdt_date)), 4)
)
glue <- cohort %>%
  filter(!is.na(cwt_dtt_date)) %>%
  mutate(dtt_to_cwt_treat = as.integer(cwt_treat_date - cwt_dtt_date),
         cwt_vs_first_tx  = as.integer(cwt_treat_date - first_tx_date))
profile$cwt_glue <- list(
  dtt_to_cwt_treat = q_sum(glue$dtt_to_cwt_treat),
  cwt_vs_first_tx  = q_sum(glue$cwt_vs_first_tx)
)
profile$cwt_agreement <- glue %>%
  filter(!is.na(cwt_vs_first_tx)) %>%
  summarise(pct_exact     = round(mean(cwt_vs_first_tx == 0), 4),
            pct_within_14 = round(mean(abs(cwt_vs_first_tx) <= 14), 4))
profile$mdt_to_dtt <- cohort %>%
  filter(!is.na(cwt_mdt_date), !is.na(cwt_dtt_date)) %>%
  mutate(d = as.integer(cwt_dtt_date - cwt_mdt_date)) %>% pull(d) %>% q_sum()

# =============================================================================
# Minimal spec: column manifest + merge constants for the generator/merge
# =============================================================================
minimal_spec <- list(
  table_a_cols = c(
    "pseudo_patientid","diagmdy","ydiag",
    "sex","agediag","ethnicity_group_broad","NHSE_reversed_imd_quintile_lsoas",
    "tumour_site_grp","cancer_subtype","stage_clean","route_combined","cci_group",
    "diag_trust","diag_hosp",
    "endoscopy_date","emresd_date","surgery_date","sact_date","rt_date",
    "first_tx_date","tx_pathway","chemo_source","tx_trust","finmdy","died"),
  cwt_cols = c("pseudo_patientid","site_icd10","modality","crtp_date",
               "date_first_seen","mdt_date","treat_period_start","treat_start"),
  tx_pathway_levels = tx_pathway_levels,
  og_icd        = og_icd,
  stage_levels  = c("1","2","3"),
  merge_const   = list(
    tx_window_days   = 270L, dtt_min_offset = -30L, treat_tol_days = 14L,
    surg_switch_date = "2020-10-01", surg_01_rule = "date_split")
)

saveRDS(profile,      paste0(save_dir, "og_profile_for_synthetic.rds"))
saveRDS(minimal_spec, paste0(save_dir, "og_minimal_spec.rds"))
cat("Saved og_profile_for_synthetic.rds (", length(profile), "sections) and",
    "og_minimal_spec.rds\n")
