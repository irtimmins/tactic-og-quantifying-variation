########

# =============================================================================
# OG Waiting Times - profiling for synthetic data generation
# -----------------------------------------------------------------------------
# Extracts the distributional parameters needed to generate a realistic
# synthetic pre-CWT cohort + CWT dataset that merge successfully.
#
# DISCLOSURE CONTROL: counts rounded to nearest 5, cells < 10 suppressed (NA),
# interval summaries reported as quantiles (no raw min/max, no individual rows).
# Output is aggregate only. STILL run the saved object past your output-checking
# / DAR process before it leaves the environment.
#
# Produces: og_profile (list) -> og_profile_for_synthetic.rds
# Inputs:   og_cohort_cwt_2015_2022.rds  (post-merge cohort)
#           the partitioned CWT dataset   (for raw per-record stats)
# =============================================================================

library(tidyverse)
library(arrow)
library(lubridate)

library(here)
icon_dir <- here("Data", "ICON")        # real ICON inputs
syn_dir  <- here("Data", "synthetic")   # synthetic-build outputs
dir.create(syn_dir, recursive = TRUE, showWarnings = FALSE)

SDC_MIN  <- 10L   # suppress counts below this
spec_obj <- readRDS(file.path(syn_dir, "og_pipeline_spec.rds"))
og_icd10 <- spec_obj$og_icd10

og <- readRDS(file.path(icon_dir, "og_cohort_cwt_2015_2022.rds"))

# -----------------------------------------------------------------------------
# Disclosure-safe helpers
# -----------------------------------------------------------------------------
r5  <- function(x) round(x / 5) * 5                      # round counts to 5
sup <- function(n) ifelse(n < SDC_MIN, NA_real_, r5(n)) # suppress small cells

# Categorical marginal (safe counts + proportions)
cat_marg <- function(df, var) {
  df %>%
    count(.data[[var]], name = "n") %>%
    mutate(n_safe = sup(n),
           prop   = round(n / sum(n), 4)) %>%       # prop from true n, then drop n
    select(level = 1, n_safe, prop)
}

# Quantile summary of a numeric/interval vector
q_sum <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < SDC_MIN) return(NULL)
  tibble(
    n_nonmiss = r5(length(x)),
    p05 = quantile(x, .05), p10 = quantile(x, .10), p25 = quantile(x, .25),
    p50 = quantile(x, .50), p75 = quantile(x, .75), p90 = quantile(x, .90),
    p95 = quantile(x, .95),
    mean = round(mean(x), 1), sd = round(sd(x), 1)
  )
}

# Quantile summary of an interval, grouped by pathway
q_by_pathway <- function(df, var) {
  df %>%
    group_by(tx_pathway) %>%
    filter(sum(!is.na(.data[[var]])) >= SDC_MIN) %>%
    summarise(
      n   = r5(sum(!is.na(.data[[var]]))),
      p10 = quantile(.data[[var]], .10, na.rm = TRUE),
      p25 = quantile(.data[[var]], .25, na.rm = TRUE),
      p50 = quantile(.data[[var]], .50, na.rm = TRUE),
      p75 = quantile(.data[[var]], .75, na.rm = TRUE),
      p90 = quantile(.data[[var]], .90, na.rm = TRUE),
      mean = round(mean(.data[[var]], na.rm = TRUE), 1),
      .groups = "drop"
    )
}

og_profile <- list()

# =============================================================================
# A. Scale & temporal structure
# =============================================================================
og_profile$n_patients <- r5(nrow(og))
og_profile$by_year    <- cat_marg(og, "ydiag")

# =============================================================================
# B. Univariate registry marginals (drives covariate sampling)
# =============================================================================
og2 <- og %>%
  mutate(age_grp = cut(agediag,
                       breaks = c(-Inf,50,55,60,65,70,75,80,85,Inf),
                       labels = c("<50","50-54","55-59","60-64","65-69",
                                  "70-74","75-79","80-84","85+"),
                       right = FALSE))

og_profile$marginals <- list(
  tumour_site_grp = cat_marg(og,  "tumour_site_grp"),
  cancer_subtype  = cat_marg(og,  "cancer_subtype"),
  stage_clean     = cat_marg(og,  "stage_clean"),
  sex             = cat_marg(og,  "sex"),
  age_grp         = cat_marg(og2, "age_grp"),
  ethnicity       = cat_marg(og,  "ethnicity_group_broad"),
  imd_quintile    = cat_marg(og,  "NHSE_reversed_imd_quintile_lsoas"),
  route_combined  = cat_marg(og,  "route_combined"),
  ps_num          = cat_marg(og,  "ps_num"),
  cnsinvolved     = cat_marg(og,  "cnsinvolved"),
  died            = cat_marg(og,  "died")
)

# =============================================================================
# C. Pathway structure - marginal + conditional (the central driver)
# =============================================================================
og_profile$tx_pathway        <- cat_marg(og, "tx_pathway")

safe_xtab <- function(df, rowvar, colvar = "tx_pathway") {
  df %>%
    count(.data[[rowvar]], .data[[colvar]], name = "n") %>%
    group_by(.data[[rowvar]]) %>%
    mutate(prop = round(n / sum(n), 4), n_safe = sup(n)) %>%
    ungroup() %>%
    select(row = 1, pathway = 2, n_safe, prop)
}
og_profile$pathway_by_subtype <- safe_xtab(og, "cancer_subtype")
og_profile$pathway_by_stage   <- safe_xtab(og, "stage_clean")
og_profile$pathway_by_site    <- safe_xtab(og, "tumour_site_grp")
og_profile$pathway_by_year    <- safe_xtab(og, "ydiag")  # captures EMR/ESD post-2020 rise

# Presence-flag rates (cross-check; derivable from pathway)
og_profile$presence_rates <- og %>%
  summarise(across(c(had_emresd, had_surgery, had_curative_surgery, had_sact,
                     had_rt, had_curative_rt, had_palliative_rt,
                     received_curative_tx),
                   ~round(mean(.x, na.rm = TRUE), 4)))

# =============================================================================
# D. Anchor-date timing (relative to diagnosis) + key intervals, by pathway
# =============================================================================
og_profile$intervals_overall <- list(
  days_endo_to_dx     = q_sum(og$days_endo_to_dx),
  wt_dx_to_tx         = q_sum(og$wt_dx_to_tx),
  wt_endo_to_tx       = q_sum(og$wt_endo_to_tx),
  wt_dx_to_surg       = q_sum(og$wt_dx_to_surg),
  wt_dx_to_sact       = q_sum(og$wt_dx_to_sact),
  wt_dx_to_rt         = q_sum(og$wt_dx_to_rt),
  wt_sact_to_surg     = q_sum(og$wt_sact_to_surg),
  wt_rt_to_surg       = q_sum(og$wt_rt_to_surg),
  surv_from_dx_days = q_sum(as.integer(og$finmdy - og$diagmdy))
)
og_profile$intervals_by_pathway <- list(
  wt_dx_to_tx     = q_by_pathway(og, "wt_dx_to_tx"),
  wt_endo_to_tx   = q_by_pathway(og, "wt_endo_to_tx"),
  wt_dx_to_surg   = q_by_pathway(og, "wt_dx_to_surg"),
  wt_dx_to_sact   = q_by_pathway(og, "wt_dx_to_sact"),
  wt_sact_to_surg = q_by_pathway(og, "wt_sact_to_surg")
)
# Interval x IMD (to reproduce the mild deprivation gradient you saw)
og_profile$wt_endo_to_tx_by_imd <- og %>%
  filter(!is.na(NHSE_reversed_imd_quintile_lsoas)) %>%
  group_by(NHSE_reversed_imd_quintile_lsoas) %>%
  summarise(n = r5(n()),
            p50 = median(wt_endo_to_tx, na.rm = TRUE),
            p25 = quantile(wt_endo_to_tx, .25, na.rm = TRUE),
            p75 = quantile(wt_endo_to_tx, .75, na.rm = TRUE),
            .groups = "drop")

# =============================================================================
# E. Treatment-attribute realism (RT schedules, SACT regimens)
# =============================================================================
og_profile$rt_curative_rate <- round(mean(og$rt_curative, na.rm = TRUE), 4)
og_profile$rt_dose_fractions <- og %>%
  filter(!is.na(rt_dose), !is.na(rt_fractions)) %>%
  count(rt_dose, rt_fractions, rt_curative, name = "n") %>%
  mutate(n_safe = sup(n), prop = round(n / sum(n), 4)) %>%
  filter(!is.na(n_safe)) %>% select(-n) %>% arrange(desc(prop))
og_profile$benchmark_group <- og %>%
  filter(had_sact) %>% cat_marg("BENCHMARK_GROUP") %>% filter(!is.na(n_safe))
og_profile$chemo_radiation <- og %>% filter(had_sact) %>% cat_marg("CHEMO_RADIATION")
og_profile$sact_intent     <- og %>% filter(had_sact) %>% cat_marg("INTENT_OF_TREATMENT_V3")
og_profile$surgery_type    <- og %>% filter(had_surgery) %>% cat_marg("surgery_type")
og_profile$curative_surgery_rate <- round(mean(og$curative_surgery, na.rm = TRUE), 4)

# =============================================================================
# F. Missingness - anchor-date presence by pathway (drives NA structure)
# =============================================================================
og_profile$date_presence_by_pathway <- og %>%
  group_by(tx_pathway) %>%
  summarise(n = r5(n()),
            pct_endoscopy = round(mean(!is.na(endoscopy_date)), 4),
            pct_emresd    = round(mean(!is.na(emresd_date)),    4),
            pct_surgery   = round(mean(!is.na(surgery_date)),   4),
            pct_sact      = round(mean(!is.na(sact_date)),      4),
            pct_rt        = round(mean(!is.na(rt_date)),        4),
            pct_first_tx  = round(mean(!is.na(first_tx_date)),  4),
            .groups = "drop")

# =============================================================================
# G. Trust / provider volume structure (AGGREGATE ONLY - no codes)
# =============================================================================
trust_vol <- function(df, var) {
  v <- df[[var]]; v <- v[!is.na(v)]
  if (length(v) == 0) return(NULL)
  sizes <- as.integer(table(v))
  tibble(n_distinct = length(sizes),
         vol_p10 = quantile(sizes, .10), vol_p25 = quantile(sizes, .25),
         vol_p50 = quantile(sizes, .50), vol_p75 = quantile(sizes, .75),
         vol_p90 = quantile(sizes, .90), vol_max = max(sizes))
}
og_profile$trust_volume <- list(
  tx_trust        = trust_vol(og, "tx_trust"),
  PROCODE3        = trust_vol(og, "PROCODE3"),
  ORGCODEPROVIDER = trust_vol(og, "ORGCODEPROVIDER"),
  diag_hosp       = trust_vol(og, "diag_hosp")    # ADD THIS
)
# Change-of-trust rate (diagnosis vs treatment), on 3-char codes
og_profile$change_trust_rate <- og %>%
  filter(!is.na(diag_trust), !is.na(tx_trust)) %>%
  summarise(rate = round(mean(substr(diag_trust,1,3) != tx_trust), 4)) %>% pull(rate)
# Per-pathway breakdown
og_profile$change_trust_by_pathway <- og %>%
  filter(!is.na(diag_trust), !is.na(tx_trust)) %>%
  group_by(tx_pathway) %>%
  filter(n() >= SDC_MIN) %>%
  summarise(
    n          = r5(n()),
    pct_change = round(mean(substr(diag_trust, 1, 3) != tx_trust), 4),
    .groups    = "drop"
  )
# =============================================================================
# H. CWT structure - per-record stats from the raw partitioned dataset
# =============================================================================
ncras_og_ids <- og$pseudo_patientid
cwt_og <- open_dataset(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/CWT/11_CWT_data_partitioned"
) %>%
  filter(site_icd10 %in% og_icd10) %>%
  collect() %>%
  mutate(
    pseudo_patientid = as.character(pseudo_patientid),
    cwt_dtt_date     = as.Date(treat_period_start, format = "%d/%m/%Y"),
    cwt_treat_date   = as.Date(treat_start,        format = "%d/%m/%Y"),
    cwt_mdt_date     = as.Date(mdt_date,           format = "%d/%m/%Y")
  ) %>%
  filter(pseudo_patientid %in% ncras_og_ids)

# Records per patient distribution
og_profile$cwt_records_per_patient <- cwt_og %>%
  count(pseudo_patientid, name = "k") %>%
  count(k, name = "n_patients") %>%
  mutate(n_safe = sup(n_patients), prop = round(n_patients / sum(n_patients), 4)) %>%
  select(records = k, n_safe, prop)

# Modality marginal, and modality x pathway (conditional)
og_profile$cwt_modality <- cat_marg(cwt_og, "modality") %>% filter(!is.na(n_safe))
og_profile$cwt_modality_by_pathway <- cwt_og %>%
  left_join(og %>% select(pseudo_patientid, tx_pathway), by = "pseudo_patientid") %>%
  count(tx_pathway, modality, name = "n") %>%
  group_by(tx_pathway) %>%
  mutate(prop = round(n / sum(n), 4), n_safe = sup(n)) %>% ungroup() %>%
  filter(!is.na(n_safe)) %>% select(tx_pathway, modality, n_safe, prop)
og_profile$cwt_site_icd10 <- cat_marg(cwt_og, "site_icd10") %>% filter(!is.na(n_safe))

# Coverage: share of cohort with any CWT record / any in-window DTT
og_profile$cwt_coverage <- og %>%
  summarise(pct_any_cwt = round(mean(pseudo_patientid %in% cwt_og$pseudo_patientid), 4),
            pct_dtt_anchor = round(mean(!is.na(cwt_dtt_date)), 4))

# Field completeness on the anchored record
og_profile$cwt_completeness <- og %>%
  filter(!is.na(cwt_dtt_date)) %>%
  summarise(pct_mdt = round(mean(!is.na(cwt_mdt_date)), 4))

# MDT timing relative to DTT (for the ~40% with an MDT date)
og_profile$mdt_to_dtt <- og %>%
  filter(!is.na(cwt_mdt_date), !is.na(cwt_dtt_date)) %>%
  mutate(d = as.integer(cwt_dtt_date - cwt_mdt_date)) %>% pull(d) %>% q_sum()

# =============================================================================
# I. The merge-glue distributions (signed) - overall and by pathway
#    These reproduce the validation behaviour you saw.
# =============================================================================
val <- og %>%
  filter(!is.na(cwt_dtt_date)) %>%
  mutate(
    days_dx_to_dtt   = as.integer(cwt_dtt_date - diagmdy),
    dtt_to_cwt_treat = as.integer(cwt_treat_date - cwt_dtt_date),
    dtt_to_tx        = as.integer(first_tx_date - cwt_dtt_date),
    cwt_vs_first_tx  = as.integer(cwt_treat_date - first_tx_date)
  )

og_profile$cwt_glue_overall <- list(
  days_dx_to_dtt   = q_sum(val$days_dx_to_dtt),
  dtt_to_cwt_treat = q_sum(val$dtt_to_cwt_treat),
  cwt_vs_first_tx  = q_sum(val$cwt_vs_first_tx)
)
og_profile$cwt_glue_by_pathway <- list(
  days_dx_to_dtt   = q_by_pathway(val, "days_dx_to_dtt"),
  dtt_to_cwt_treat = q_by_pathway(val, "dtt_to_cwt_treat"),
  cwt_vs_first_tx  = q_by_pathway(val, "cwt_vs_first_tx")
)
# Agreement rates by pathway (exact / within 5 / within 14)
og_profile$cwt_agreement_by_pathway <- val %>%
  filter(!is.na(cwt_vs_first_tx)) %>%
  group_by(tx_pathway) %>%
  filter(n() >= SDC_MIN) %>%
  summarise(n = r5(n()),
            pct_exact     = round(mean(cwt_vs_first_tx == 0), 4),
            pct_within_5  = round(mean(abs(cwt_vs_first_tx) <= 5), 4),
            pct_within_14 = round(mean(abs(cwt_vs_first_tx) <= 14), 4),
            .groups = "drop")
# dtt_valid rate by pathway
og_profile$dtt_valid_by_pathway <- og %>%
  filter(!is.na(cwt_dtt_date), !is.na(dtt_valid)) %>%
  group_by(tx_pathway) %>%
  filter(n() >= SDC_MIN) %>%
  summarise(n = r5(n()), pct_valid = round(mean(dtt_valid), 4), .groups = "drop")

# =============================================================================
# Save + quick console view
# =============================================================================
saveRDS(og_profile, file.path(syn_dir, "og_profile_for_synthetic.rds"))
cat("Saved og_profile_for_synthetic.rds with", length(og_profile), "sections.\n")
str(og_profile, max.level = 1)



