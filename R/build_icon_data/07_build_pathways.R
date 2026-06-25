# =============================================================================
# 07  Build pathways
# -----------------------------------------------------------------------------
# Assembles the cohort with every treatment anchor and the comorbidity index,
# then derives the treatment classification: the presence flags, the sequencing
# flags, tx_pathway, the first curative-treatment date, and the treating trust.
# This is the step where the individual treatment dates become a pathway.
#
# The pathway is a function of the flags and dates alone. tx_trust is the trust
# of the clock-stop treatment (surgery's provider for surgical pathways, RT's for
# RT-anchored pathways; SACT's provider is never the trust source).
#
# Reads : ncras cohort, all anchors (endoscopy / emresd / surgery / chemo / rt),
#         og_cci
# Writes: Data/ICON/og_cohort_2015_2022.rds
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

ncras_og         <- readRDS(f_ncras_cohort)
endoscopy_anchor <- readRDS(f_endoscopy_anchor)
emresd_anchor    <- readRDS(f_emresd_anchor)
surgery_anchor   <- readRDS(f_surgery_anchor)
chemo_anchor     <- readRDS(f_chemo_anchor)
rt_anchor        <- readRDS(f_rt_anchor)
og_cci           <- readRDS(f_cci)

# -----------------------------------------------------------------------------
# Assemble: cohort + all anchors + comorbidity index
# -----------------------------------------------------------------------------
og_cohort <- ncras_og %>%
  left_join(endoscopy_anchor %>%
              select(pseudo_patientid, endoscopy_date, days_endo_to_dx),
            by = "pseudo_patientid") %>%
  left_join(emresd_anchor %>%
              select(pseudo_patientid, emresd_date, emresd_provider, days_dx_to_emresd),
            by = "pseudo_patientid") %>%
  left_join(surgery_anchor %>%
              select(pseudo_patientid, surgery_date, surgery_type, surgery_class,
                     opcs_primary, PROCODE3, SITETRET, days_dx_to_surg,
                     curative_surgery),
            by = "pseudo_patientid") %>%
  # chemo_anchor supplies the SACT-preferred combined chemo date as sact_date,
  # plus the HES chemo date and provenance for the chemoRT guard below
  left_join(chemo_anchor %>%
              select(pseudo_patientid,
                     sact_date = chemo_date, days_dx_to_sact = days_dx_to_chemo,
                     chemo_source, hes_chemo_date,
                     BENCHMARK_GROUP, benchmark_group_lwr,
                     INTENT_OF_TREATMENT_V3, CHEMO_RADIATION,
                     ORGANISATION_CODE_OF_PROVIDER),
            by = "pseudo_patientid") %>%
  left_join(rt_anchor %>%
              select(pseudo_patientid, rt_date, rt_curative,
                     rt_dose, rt_fractions, days_dx_to_rt, ORGCODEPROVIDER),
            by = "pseudo_patientid") %>%
  left_join(og_cci %>%
              select(pseudo_patientid, rcs_ch_score, cci_any, cci_group,
                     cci_n_conditions, cci_conditions),
            by = "pseudo_patientid")

# -----------------------------------------------------------------------------
# Treatment flags, sequencing, pathway, dates and trust
# -----------------------------------------------------------------------------
og_cohort <- og_cohort %>%
  mutate(
    # presence flags
    had_emresd           = !is.na(emresd_date),
    had_surgery          = !is.na(surgery_date),
    had_curative_surgery = !is.na(surgery_date) & curative_surgery == TRUE,
    had_sact             = !is.na(sact_date),
    had_rt               = !is.na(rt_date),
    had_curative_rt      = !is.na(rt_date) & rt_curative == TRUE,
    had_palliative_rt    = !is.na(rt_date) & rt_curative == FALSE,
    
    # chemo eligible to define non-surgical definitive chemoRT: SACT chemo always
    # counts; HES-only chemo only when within hes_chemo_near_rt_days of the RT, so
    # a separate HES chemo episode cannot manufacture definitive chemoRT.
    had_chemo_for_chemort = had_sact &
      ( coalesce(chemo_source, "sact") != "hes" |
          ( !is.na(hes_chemo_date) & !is.na(rt_date) &
              abs(as.integer(hes_chemo_date - rt_date)) <= hes_chemo_near_rt_days ) ),
    
    # sequencing flags
    sact_before_surgery = had_sact & had_surgery & sact_date < surgery_date,
    sact_after_surgery  = had_sact & had_surgery & sact_date > surgery_date,
    rt_before_surgery   = had_rt   & had_surgery & rt_date   < surgery_date,
    rt_after_surgery    = had_rt   & had_surgery & rt_date   > surgery_date,
    concurrent_chemo_rt = had_sact & had_curative_rt &
      abs(as.integer(sact_date - rt_date)) <= chemo_rt_concurrent_days,
    
    received_curative_tx = had_emresd | had_curative_surgery | had_curative_rt,
    
    # pathway classification (first matching branch wins)
    tx_pathway = case_when(
      had_emresd & !had_surgery & !had_sact & !concurrent_chemo_rt ~ "EMR/ESD only",
      had_emresd & had_surgery                                     ~ "EMR/ESD then surgery",
      had_surgery & sact_before_surgery & rt_before_surgery        ~ "Surgery + neoadjuvant chemoRT",
      had_surgery & sact_before_surgery & !rt_before_surgery       ~ "Surgery + neoadjuvant chemo",
      had_surgery & rt_before_surgery & !sact_before_surgery       ~ "Surgery + neoadjuvant RT",
      had_surgery & sact_after_surgery & !sact_before_surgery      ~ "Surgery + adjuvant chemo",
      had_surgery & !had_sact & !concurrent_chemo_rt               ~ "Surgery only",
      had_surgery                                                  ~ "Surgery + other",
      !had_surgery & had_curative_rt & had_chemo_for_chemort       ~ "Definitive chemoRT",
      !had_surgery & had_curative_rt & !had_chemo_for_chemort      ~ "Curative RT only",
      !had_surgery & had_palliative_rt & had_sact                  ~ "Palliative chemo + RT",
      !had_surgery & had_sact & !had_curative_rt                   ~ "SACT only",
      !had_surgery & had_palliative_rt & !had_sact                 ~ "Palliative RT only",
      TRUE                                                         ~ "No treatment recorded"),
    
    # first curative-treatment date (the clock-stop). Neoadjuvant RT/chemoRT sets
    # the date; neoadjuvant chemo alone does not (surgery sets it).
    first_tx_date = case_when(
      tx_pathway %in% c("EMR/ESD only", "EMR/ESD then surgery") ~ emresd_date,
      tx_pathway == "Surgery + neoadjuvant chemoRT" ~ pmin(sact_date, rt_date, na.rm = TRUE),
      tx_pathway == "Surgery + neoadjuvant RT"      ~ rt_date,
      tx_pathway == "Surgery + neoadjuvant chemo"   ~ sact_date,
      tx_pathway %in% c("Surgery + adjuvant chemo",
                        "Surgery only", "Surgery + other") ~ surgery_date,
      tx_pathway == "Definitive chemoRT" ~ pmin(sact_date, rt_date, na.rm = TRUE),
      tx_pathway == "Curative RT only"   ~ rt_date,
      TRUE                               ~ as.Date(NA)),
    
    # treating trust = provider of the clock-stop treatment. EMR-only patients
    # have no surgery record, so their trust comes from the EMR provider; EMR-then-
    # surgery and the surgical pathways take the surgery provider (PROCODE3).
    tx_trust = case_when(
      tx_pathway == "EMR/ESD only"                         ~ substr(emresd_provider, 1, 3),
      tx_pathway %in% c("EMR/ESD then surgery",
                        "Surgery + neoadjuvant chemo", "Surgery + adjuvant chemo",
                        "Surgery only", "Surgery + other") ~ substr(PROCODE3, 1, 3),
      tx_pathway %in% c("Surgery + neoadjuvant chemoRT", "Surgery + neoadjuvant RT",
                        "Definitive chemoRT", "Curative RT only") ~ substr(ORGCODEPROVIDER, 1, 3),
      TRUE                                                 ~ NA_character_),
    
    # waiting-time components (the CWT-anchored family is added in 08)
    wt_dx_to_tx     = as.integer(first_tx_date - diagmdy),
    wt_endo_to_tx   = as.integer(first_tx_date - endoscopy_date),
    wt_dx_to_surg   = as.integer(surgery_date - diagmdy),
    wt_endo_to_surg = as.integer(surgery_date - endoscopy_date),
    wt_dx_to_sact   = as.integer(sact_date - diagmdy),
    wt_endo_to_sact = as.integer(sact_date - endoscopy_date),
    wt_sact_to_surg = as.integer(surgery_date - sact_date),
    wt_surg_to_sact = as.integer(sact_date - surgery_date),
    wt_dx_to_rt     = as.integer(rt_date - diagmdy),
    wt_endo_to_rt   = as.integer(rt_date - endoscopy_date),
    wt_rt_to_surg   = as.integer(surgery_date - rt_date),
    
    # survival from surgery (PI7/PI8)
    surv_from_surg_days = as.integer(finmdy - surgery_date),
    alive_90d_post_surg = had_surgery & !is.na(surv_from_surg_days) &
      (surv_from_surg_days > 90  | died == 0L),
    alive_1yr_post_surg = had_surgery & !is.na(surv_from_surg_days) &
      (surv_from_surg_days > 365 | died == 0L)
  )

saveRDS(og_cohort, f_cohort)
cat("Saved", f_cohort, "(", nrow(og_cohort), "patients ).",
    "Next: 08_merge_cwt.R\n")

# ---- optional checks (uncomment to inspect) ---------------------------------
# og_cohort %>% count(tx_pathway) %>% mutate(pct = round(100*n/sum(n),1)) %>%
#   arrange(desc(n)) %>% print(n = 20)
# og_cohort %>% filter(stage_clean %in% c("1","2","3")) %>%
#   summarise(curative = round(100*mean(received_curative_tx)),
#             any_tx = round(100*mean(tx_pathway != "No treatment recorded"))) %>% print()