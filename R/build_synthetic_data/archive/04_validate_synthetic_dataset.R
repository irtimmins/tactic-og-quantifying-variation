# =============================================================================
# OG Waiting Times - synthetic data validation
# -----------------------------------------------------------------------------
# Two checks on the generated pre-CWT cohort:
#   1. CONFORMANCE - does og_cohort_precwt_SYNTH.rds match the spec (columns,
#      types, pathway levels, unique IDs)? Confirms it is a drop-in for the
#      real og_cohort_precwt_spec_*.rds.
#   2. INTERNAL CONSISTENCY - re-derive tx_pathway and first_tx_date from the
#      stored anchor dates and confirm they match the stored columns. Confirms
#      the saved artefact stands on its own (independent of the generator).
#
# Run OFF-server on the synthetic files. Needs og_pipeline_spec.rds (the spec
# manifest is aggregate/non-disclosive and travels with the bundle).
# =============================================================================

library(dplyr)
library(tidyr)
library(stringr)

library(here)
syn_dir <- here("Data", "synthetic")
spec_obj          <- readRDS(file.path(syn_dir, "og_pipeline_spec.rds"))
pre_cwt_spec      <- spec_obj$pre_cwt_spec
tx_pathway_levels <- spec_obj$tx_pathway_levels
syn               <- readRDS(file.path(syn_dir, "og_cohort_precwt_SYNTH.rds"))

# -----------------------------------------------------------------------------
# 1. Conformance
# -----------------------------------------------------------------------------
check_conformance <- function(df, spec) {
  present          <- spec$name %in% names(df)
  missing_required <- spec$name[!present & spec$tier == "required"]
  missing_core     <- spec$name[!present & spec$tier == "core"]
  
  type_of <- function(x) {
    if (inherits(x, "Date")) "Date" else if (is.factor(x)) "factor" else
      if (is.logical(x)) "logical" else if (is.integer(x)) "integer" else
        if (is.numeric(x)) "numeric" else if (is.character(x)) "character" else class(x)[1]
  }
  rows <- spec %>% filter(name %in% names(df))
  obs  <- vapply(rows$name, function(n) type_of(df[[n]]), character(1))
  compat <- function(e, o) (e == o) ||
    (e %in% c("integer","numeric") && o %in% c("integer","numeric"))
  mism <- rows %>% mutate(observed = obs) %>%
    filter(!mapply(compat, type, observed)) %>% select(name, expected = type, observed)
  
  bad_pathway <- setdiff(unique(na.omit(df$tx_pathway)), tx_pathway_levels)
  bad_stage   <- setdiff(unique(na.omit(df$stage_clean)), c("1","2","3"))
  
  cat("== 1. Conformance ==\n")
  cat("Rows:", nrow(df), " Cols:", ncol(df), "\n")
  cat("Missing REQUIRED cols:", if (length(missing_required)) paste(missing_required, collapse=", ") else "none", "\n")
  cat("Missing core cols:    ", if (length(missing_core)) paste(missing_core, collapse=", ") else "none", "\n")
  if (nrow(mism)) { cat("Type mismatches:\n"); print(mism) } else cat("Type mismatches:     none\n")
  cat("Unexpected tx_pathway:", if (length(bad_pathway)) paste(bad_pathway, collapse=", ") else "none", "\n")
  cat("Unexpected stage_clean:", if (length(bad_stage)) paste(bad_stage, collapse=", ") else "none", "\n")
  cat("Duplicate patient IDs:", sum(duplicated(df$pseudo_patientid)), "\n\n")
}
check_conformance(syn, pre_cwt_spec)

# -----------------------------------------------------------------------------
# 2. Internal consistency: re-derive pathway + first_tx_date from stored dates
# -----------------------------------------------------------------------------
re <- syn %>% transmute(
  pseudo_patientid,
  had_emresd          = !is.na(emresd_date),
  had_surgery         = !is.na(surgery_date),
  had_curative_surg   = !is.na(surgery_date) & curative_surgery == TRUE,
  had_sact            = !is.na(sact_date),
  had_curative_rt     = !is.na(rt_date) & rt_curative == TRUE,
  had_palliative_rt   = !is.na(rt_date) & rt_curative == FALSE,
  sact_before_surgery = had_sact & had_surgery & sact_date < surgery_date,
  sact_after_surgery  = had_sact & had_surgery & sact_date > surgery_date,
  rt_before_surgery   = (!is.na(rt_date)) & had_surgery & rt_date < surgery_date,
  concurrent_chemo_rt = had_sact & had_curative_rt &
    abs(as.integer(sact_date - rt_date)) <= 14,
  re_pathway = case_when(
    had_emresd & !had_surgery & !had_sact & !concurrent_chemo_rt ~ "EMR/ESD only",
    had_emresd & had_surgery                                     ~ "EMR/ESD then surgery",
    had_surgery & sact_before_surgery & rt_before_surgery        ~ "Surgery + neoadjuvant chemoRT",
    had_surgery & sact_before_surgery & !rt_before_surgery       ~ "Surgery + neoadjuvant chemo",
    had_surgery & rt_before_surgery & !sact_before_surgery       ~ "Surgery + neoadjuvant RT",
    had_surgery & sact_after_surgery & !sact_before_surgery      ~ "Surgery + adjuvant chemo",
    had_surgery & !had_sact & !concurrent_chemo_rt               ~ "Surgery only",
    had_surgery                                                  ~ "Surgery + other",
    !had_surgery & had_curative_rt & had_sact                    ~ "Definitive chemoRT",
    !had_surgery & had_curative_rt & !had_sact                   ~ "Curative RT only",
    !had_surgery & had_palliative_rt & had_sact                  ~ "Palliative chemo + RT",
    !had_surgery & had_sact & !had_curative_rt                   ~ "SACT only",
    !had_surgery & had_palliative_rt & !had_sact                 ~ "Palliative RT only",
    TRUE                                                         ~ "No treatment recorded"
  )
) %>%
  left_join(syn %>% select(pseudo_patientid, emresd_date, surgery_date,
                           sact_date, rt_date, tx_pathway, first_tx_date),
            by = "pseudo_patientid") %>%
  mutate(re_first_tx = case_when(
    re_pathway %in% c("EMR/ESD only","EMR/ESD then surgery") ~ emresd_date,
    re_pathway == "Surgery + neoadjuvant chemoRT"           ~ pmin(sact_date, rt_date, na.rm = TRUE),
    re_pathway == "Surgery + neoadjuvant RT"                ~ rt_date,
    re_pathway == "Surgery + neoadjuvant chemo"             ~ sact_date,
    re_pathway %in% c("Surgery + adjuvant chemo","Surgery only","Surgery + other") ~ surgery_date,
    re_pathway == "Definitive chemoRT"                      ~ pmin(sact_date, rt_date, na.rm = TRUE),
    re_pathway == "Curative RT only"                        ~ rt_date,
    TRUE                                                    ~ as.Date(NA)
  ))

cat("== 2. Internal consistency ==\n")
cat("tx_pathway re-derivation match: ",
    round(100 * mean(re$re_pathway == re$tx_pathway), 2), "%\n")
cat("first_tx_date match:            ",
    round(100 * mean(re$re_first_tx == re$first_tx_date |
                       (is.na(re$re_first_tx) & is.na(re$first_tx_date))), 2), "%\n")
mismatch <- re %>% filter(re_pathway != tx_pathway)
if (nrow(mismatch)) {
  cat("Pathway mismatches by type:\n")
  mismatch %>% count(tx_pathway, re_pathway, sort = TRUE) %>% print(n = 20)
} else cat("No pathway mismatches.\n")

# -----------------------------------------------------------------------------
# 3. (Optional) synthetic vs profile pathway mix - quick eyeball
# -----------------------------------------------------------------------------

if (file.exists(file.path(syn_dir, "og_profile_for_synthetic.rds"))) {
  prof <- readRDS((file.path(syn_dir, "og_profile_for_synthetic.rds")))
  cmp <- syn %>% count(tx_pathway, name = "n_syn") %>%
    mutate(prop_syn = round(n_syn / sum(n_syn), 4)) %>%
    left_join(prof$tx_pathway %>% select(tx_pathway = level, prop_real = prop),
              by = "tx_pathway") %>%
    arrange(desc(prop_syn))
  cat("\n== 3. Pathway mix: synthetic vs profile ==\n")
  print(cmp, n = 20)
}
