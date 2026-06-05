# =============================================================================
# OG Waiting Times - synthetic data generators (pre-CWT cohort + CWT)
# -----------------------------------------------------------------------------
# Driven entirely by og_profile_for_synthetic.rds (aggregate, disclosure-checked).
# Runs OFF the ICON server - output is fully synthetic and shareable.
#
# Produces:
#   syn_cohort  -> og_cohort_precwt_SYNTH.rds   (Table A, spec-conformant)
#   syn_cwt     -> cwt_records_SYNTH.rds         (Table B, raw dd/mm/yyyy dates)
# and runs the merge on the synthetic data to compare against the real targets.
#
# Deps: tidyverse, lubridate.  Input: og_profile_for_synthetic.rds
# =============================================================================

library(tidyverse)
library(lubridate)
base_dir <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/"
profile <- readRDS(paste0(base_dir, "og_profile_for_synthetic.rds"))
set.seed(20260601)

N_TOTAL        <- profile$n_patients      # match real scale; override if desired
tx_window_days <- 270L
tx_pathway_levels <- profile$tx_pathway$level

# =============================================================================
# Samplers
# =============================================================================
# Categorical draw from a {level, prop} table (keeps NA levels)
sample_cat <- function(tbl, n) {
  p <- tbl$prop; p[is.na(p)] <- 0; p <- p / sum(p)
  sample(tbl$level, n, replace = TRUE, prob = p)
}

# Quantile-function sampler from 5 inner quantiles (p10..p90)
qv5 <- function(v5, n) {
  pr <- c(.10, .25, .50, .75, .90)
  u  <- runif(n, .10, .90)
  round(approx(pr, v5, u, rule = 2)$y)
}
# from overall row (p05..p95)
qov <- function(row, n) {
  pr <- c(.05, .10, .25, .50, .75, .90, .95)
  v  <- c(row$p05, row$p10, row$p25, row$p50, row$p75, row$p90, row$p95)
  u  <- runif(n, .05, .95)
  round(approx(pr, v, u, rule = 2)$y)
}
# by-pathway interval row -> 5 quantiles, with overall fallback
iv5 <- function(by_tbl, ov_row, pw) {
  r <- by_tbl %>% filter(tx_pathway == pw)
  if (nrow(r) == 1) c(r$p10, r$p25, r$p50, r$p75, r$p90)
  else              c(ov_row$p10, ov_row$p25, ov_row$p50, ov_row$p75, ov_row$p90)
}

# Conditional sampler: P(level | pathway) from P(pathway|level)*N(level)
build_cond <- function(tbl, Nlev) {
  tbl %>%
    mutate(N = Nlev[as.character(row)], joint = prop * N) %>%
    filter(!is.na(joint), joint > 0) %>%
    group_by(pathway) %>%
    summarise(levels = list(row), w = list(joint / sum(joint)), .groups = "drop")
}
sample_cond <- function(cond, pathways) {
  out <- vector("list", length(pathways))
  for (p in unique(pathways)) {
    idx <- which(pathways == p)
    r   <- cond %>% filter(pathway == p)
    if (nrow(r) == 0) { out[idx] <- NA; next }
    out[idx] <- as.list(sample(r$levels[[1]], length(idx), TRUE, r$w[[1]]))
  }
  unlist(out)
}

# Level totals from a marginal table
Nlev_of <- function(marg) setNames(round(marg$prop * N_TOTAL), as.character(marg$level))

# RT schedule sampler (curative TRUE/FALSE) -> dose & fractions
sample_rt <- function(n, curative) {
  s <- profile$rt_dose_fractions %>% filter(rt_curative == curative)
  i <- sample(seq_len(nrow(s)), n, TRUE, s$prop / sum(s$prop))
  list(dose = s$rt_dose[i], fr = s$rt_fractions[i])
}

# =============================================================================
# Trust pool (fake 3-char codes with a skewed volume distribution)
# =============================================================================
n_trust <- profile$trust_volume$tx_trust$n_distinct
trust_codes <- sprintf("T%02d", seq_len(n_trust))
# skewed weights (a few high-volume trusts) ~ lognormal
trust_w <- sort(rlnorm(n_trust, meanlog = 1, sdlog = 1), decreasing = TRUE)
trust_w <- trust_w / sum(trust_w)
draw_trust <- function(n) sample(trust_codes, n, TRUE, trust_w)

# ---- Diagnosing-hospital pool, NESTED within trusts -------------------------
# Real structure: a Trust operates one or more hospital sites; each hospital
# belongs to exactly one Trust. We allocate the H distinct hospitals across the
# T trusts -- every trust gets >=1, and the remainder go preferentially to
# larger trusts (so big trusts have more sites). A patient's diag_hosp is then
# drawn ONLY from hospitals inside their diag_trust, which guarantees nesting.
n_hosp <- tryCatch(max(profile$trust_volume$diag_hosp$n_distinct, n_trust),
                   error = function(e) n_trust * 2L)

hosp_per_trust <- rep(1L, n_trust)                       # every trust gets >=1
extra <- n_hosp - n_trust
if (extra > 0) {
  # hand the leftover hospitals to trusts in proportion to trust size
  add <- tabulate(sample(seq_len(n_trust), extra, replace = TRUE, prob = trust_w),
                  nbins = n_trust)
  hosp_per_trust <- hosp_per_trust + add
}
# NHS ODS convention: trust = "R" + 2 alphanumeric chars (e.g. "RXX", "RA9")
#                   hospital = trust code + 2-char site suffix (e.g. "RXXA1", "RXX01")
# First 3 chars of any hospital code always equal its parent trust code.
chars     <- c(LETTERS, 0:9)                          # A-Z + 0-9 = 36 chars
all_2char <- as.vector(outer(chars, chars, paste0))   # 1296 unique 2-char combos

# REPLACE: trust_codes <- sprintf("T%02d", seq_len(n_trust))
trust_codes <- paste0("R", sample(all_2char, n_trust))   # e.g. "RXX", "RA9", "RJ1"

# REPLACE: hosp_codes <- sprintf("H%04d", seq_len(sum(hosp_per_trust)))
#          hosp_trust <- rep(trust_codes, times = hosp_per_trust)
hosp_codes <- unlist(mapply(
  function(t, nh) paste0(t, sample(all_2char, nh)),   # trust prefix + unique site suffix
  trust_codes, hosp_per_trust, SIMPLIFY = FALSE
))
hosp_trust <- rep(trust_codes, times = hosp_per_trust)

# within-trust site weights: one main site usually handles most activity
hosp_w <- rlnorm(length(hosp_codes), meanlog = 0, sdlog = 0.8)
hosp_w <- ave(hosp_w, hosp_trust, FUN = function(x) x / sum(x))

# draw a hospital for each patient, nested inside their (diagnosing) trust
draw_hosp <- function(trusts) {
  out <- character(length(trusts))
  for (t in unique(trusts)) {
    idx <- which(trusts == t)
    h   <- hosp_codes[hosp_trust == t]
    w   <- hosp_w[hosp_trust == t]
    out[idx] <- sample(h, length(idx), replace = TRUE, prob = w)
  }
  out
}

# =============================================================================
# 1. Patient covariates
# =============================================================================
cohort <- tibble(
  pseudo_patientid = sprintf("S%07d", seq_len(N_TOTAL)),
  pseudo_tumourid  = sprintf("T%07d", seq_len(N_TOTAL))
)

# pathway ~ marginal
cohort$tx_pathway_asg <- sample_cat(profile$tx_pathway, N_TOTAL)

# stage | pathway ; subtype | pathway ; year | pathway
cond_stage <- build_cond(profile$pathway_by_stage,
                         Nlev_of(profile$marginals$stage_clean))
cond_sub   <- build_cond(profile$pathway_by_subtype,
                         Nlev_of(profile$marginals$cancer_subtype))
cond_year  <- build_cond(profile$pathway_by_year,
                         Nlev_of(profile$by_year))

cohort$stage_clean    <- sample_cond(cond_stage, cohort$tx_pathway_asg)
cohort$cancer_subtype <- sample_cond(cond_sub,   cohort$tx_pathway_asg)
cohort$ydiag          <- as.integer(sample_cond(cond_year, cohort$tx_pathway_asg))

# site from subtype; morphology proxy
cohort <- cohort %>%
  mutate(
    tumour_site_grp = if_else(cancer_subtype == "Gast" & !is.na(cancer_subtype),
                              "gastric", "oesophageal"),
    morphology_num = case_when(
      cancer_subtype == "Oes SCC" ~ 8070L,
      cancer_subtype == "Oes ACA" ~ 8140L,
      cancer_subtype == "Gast"    ~ 8140L,
      TRUE                        ~ 8140L
    )
  )

# other covariates ~ marginals (independent)
cohort$sex                              <- as.integer(sample_cat(profile$marginals$sex, N_TOTAL))
cohort$ethnicity_group_broad            <- sample_cat(profile$marginals$ethnicity, N_TOTAL)
cohort$NHSE_reversed_imd_quintile_lsoas <- sample_cat(profile$marginals$imd_quintile, N_TOTAL)
cohort$route_combined                   <- sample_cat(profile$marginals$route_combined, N_TOTAL)
cohort$ps_num                           <- as.integer(sample_cat(profile$marginals$ps_num, N_TOTAL))
cohort$cnsinvolved                      <- as.integer(sample_cat(profile$marginals$cnsinvolved, N_TOTAL))
cohort$died                             <- as.integer(sample_cat(profile$marginals$died, N_TOTAL))

# age: sample band then a value within it
age_band <- sample_cat(profile$marginals$age_grp, N_TOTAL)
band_lo <- c("<50"=40,"50-54"=50,"55-59"=55,"60-64"=60,"65-69"=65,
             "70-74"=70,"75-79"=75,"80-84"=80,"85+"=85)
band_hi <- c("<50"=49,"50-54"=54,"55-59"=59,"60-64"=64,"65-69"=69,
             "70-74"=74,"75-79"=79,"80-84"=84,"85+"=95)
cohort$agediag <- round(runif(N_TOTAL, band_lo[age_band], band_hi[age_band]))

# diagnosis date: random day in sampled year
cohort$diagmdy <- as.Date(paste0(cohort$ydiag, "-01-01")) +
  as.integer(runif(N_TOTAL, 0, 364))

# per-pathway probability that the diagnosis trust differs from the treatment trust
ctp <- setNames(profile$change_trust_by_pathway$pct_change,
                profile$change_trust_by_pathway$tx_pathway)
p_change <- ctp[cohort$tx_pathway_asg]
p_change[is.na(p_change)] <- profile$change_trust_rate    # fallback = overall rate

cohort$treat_trust <- draw_trust(N_TOTAL)                 # the treating trust
diff_trust         <- runif(N_TOTAL) < p_change           # did they move trust?
cohort$diag_trust  <- if_else(diff_trust, draw_trust(N_TOTAL), cohort$treat_trust)
cohort$first_trust <- cohort$diag_trust
cohort$diag_hosp   <- draw_hosp(cohort$diag_trust)        # NESTED within diag_trust
cohort$canalliance_2024_code <- paste0("CA", sample(1:21, N_TOTAL, TRUE))
cohort$canalliance_2024_name <- cohort$canalliance_2024_code

# LSOA 2011 codes: realistic format E01xxxxxx (England), with IMD quintile
# baked in via the suffix range so deprivation geography is internally consistent.
# Each quintile gets its own block of LSOA numbers:
#   Q1 (most deprived) -> E01000001-E01010000
#   Q2                 -> E01010001-E01020000
#   Q3                 -> E01020001-E01030000
#   Q4                 -> E01030001-E01040000
#   Q5 (least)         -> E01040001-E01050000
imd_lsoa_base <- c(
  "1 - most deprived" = 0L,
  "2"                 = 10000L,
  "3"                 = 20000L,
  "4"                 = 30000L,
  "5 - least deprived"= 40000L
)
lsoa_base  <- imd_lsoa_base[cohort$NHSE_reversed_imd_quintile_lsoas]
lsoa_base[is.na(lsoa_base)] <- 0L    # fallback for missing IMD
cohort$lsoa11_code <- sprintf("E01%06d",
                              lsoa_base + sample(1:10000, N_TOTAL, replace = TRUE))

# survival
surv <- qov(profile$intervals_overall$surv_from_dx_days, N_TOTAL)
cohort$finmdy <- cohort$diagmdy + pmax(surv, 1L)

# =============================================================================
# 2. Treatment dates per pathway  (constructed so the build re-derives pathway)
# =============================================================================
endo <- emresd <- surgery <- sact <- rt <- as.Date(rep(NA, N_TOTAL))
rt_dose <- rt_fr <- rep(NA_real_, N_TOTAL); rt_cur <- rep(NA, N_TOTAL)

pe_tbl <- profile$date_presence_by_pathway
ov_tx   <- profile$intervals_overall$wt_dx_to_tx
ov_surg <- profile$intervals_overall$wt_dx_to_surg
ov_sact <- profile$intervals_overall$wt_dx_to_sact
ov_rt   <- profile$intervals_overall$wt_dx_to_rt
by_tx   <- profile$intervals_by_pathway$wt_dx_to_tx
by_surg <- profile$intervals_by_pathway$wt_dx_to_surg
by_sact <- profile$intervals_by_pathway$wt_dx_to_sact
by_ss   <- profile$intervals_by_pathway$wt_sact_to_surg

asg <- cohort$tx_pathway_asg
dx  <- cohort$diagmdy

for (pw in unique(asg)) {
  idx <- which(asg == pw); m <- length(idx); if (m == 0) next
  d <- dx[idx]
  
  # endoscopy present at observed per-pathway rate
  pe <- pe_tbl$pct_endoscopy[pe_tbl$tx_pathway == pw]
  if (length(pe) == 1) {
    has_e <- runif(m) < pe
    e <- d - qov(profile$intervals_overall$days_endo_to_dx, m)
    e[!has_e] <- NA; endo[idx] <- e
  }
  
  if (pw == "EMR/ESD only") {
    emresd[idx] <- d + qv5(iv5(by_tx, ov_tx, pw), m)
    
  } else if (pw == "EMR/ESD then surgery") {
    em <- d + qv5(iv5(by_tx, ov_tx, pw), m)
    sg <- d + qv5(iv5(by_surg, ov_surg, pw), m)
    sg <- pmax(sg, em + 1L)
    emresd[idx] <- em; surgery[idx] <- sg
    
  } else if (pw == "Surgery + neoadjuvant chemo") {
    sa <- d + qv5(iv5(by_sact, ov_sact, pw), m)
    sg <- sa + pmax(qv5(iv5(by_ss, ov_surg, pw), m), 1L)
    sact[idx] <- sa; surgery[idx] <- sg
    
  } else if (pw == "Surgery + neoadjuvant chemoRT") {
    sa <- d + qv5(iv5(by_sact, ov_sact, pw), m)
    rr <- sa + sample(-7:7, m, TRUE)                 # concurrent (<=14d)
    sg <- sa + pmax(qv5(iv5(by_ss, ov_surg, pw), m), 1L)
    sg <- pmax(sg, pmax(sa, rr) + 1L)
    s  <- sample_rt(m, TRUE)
    sact[idx] <- sa; rt[idx] <- rr; surgery[idx] <- sg
    rt_dose[idx] <- s$dose; rt_fr[idx] <- s$fr; rt_cur[idx] <- TRUE
    
  } else if (pw == "Surgery + neoadjuvant RT") {
    rr <- d + qv5(iv5(by_tx, ov_tx, pw), m)
    sg <- d + qv5(iv5(by_surg, ov_surg, pw), m)
    sg <- pmax(sg, rr + 1L)
    s  <- sample_rt(m, TRUE)
    rt[idx] <- rr; surgery[idx] <- sg
    rt_dose[idx] <- s$dose; rt_fr[idx] <- s$fr; rt_cur[idx] <- TRUE
    
  } else if (pw == "Surgery + adjuvant chemo") {
    sg <- d + qv5(iv5(by_surg, ov_surg, pw), m)
    gap <- pmax(-qv5(iv5(by_ss, ov_sact, pw), m), 1L)   # wt_sact_to_surg is negative
    surgery[idx] <- sg; sact[idx] <- sg + gap
    
  } else if (pw %in% c("Surgery only", "Surgery + other")) {
    surgery[idx] <- d + qv5(iv5(by_surg, ov_surg, pw), m)
    
  } else if (pw == "Definitive chemoRT") {
    sa <- d + qv5(iv5(by_sact, ov_sact, pw), m)
    rr <- sa + sample(-7:7, m, TRUE)
    s  <- sample_rt(m, TRUE)
    sact[idx] <- sa; rt[idx] <- rr
    rt_dose[idx] <- s$dose; rt_fr[idx] <- s$fr; rt_cur[idx] <- TRUE
    
  } else if (pw == "Curative RT only") {
    rr <- d + qv5(iv5(by_tx, ov_tx, pw), m)
    s  <- sample_rt(m, TRUE)
    rt[idx] <- rr; rt_dose[idx] <- s$dose; rt_fr[idx] <- s$fr; rt_cur[idx] <- TRUE
    
  } else if (pw == "Palliative chemo + RT") {
    sa <- d + qv5(iv5(by_sact, ov_sact, pw), m)
    rr <- sa + sample(-7:7, m, TRUE)
    s  <- sample_rt(m, FALSE)
    sact[idx] <- sa; rt[idx] <- rr
    rt_dose[idx] <- s$dose; rt_fr[idx] <- s$fr; rt_cur[idx] <- FALSE
    
  } else if (pw == "SACT only") {
    sact[idx] <- d + qv5(iv5(by_sact, ov_sact, pw), m)
    
  } else if (pw == "Palliative RT only") {
    rr <- d + qov(ov_rt, m)
    s  <- sample_rt(m, FALSE)
    rt[idx] <- rr; rt_dose[idx] <- s$dose; rt_fr[idx] <- s$fr; rt_cur[idx] <- FALSE
  }
  # "No treatment recorded": leave all dates NA
}

cohort <- cohort %>%
  mutate(
    endoscopy_date = endo, emresd_date = emresd, surgery_date = surgery,
    sact_date = sact, rt_date = rt,
    rt_dose = rt_dose, rt_fractions = as.integer(rt_fr), rt_curative = rt_cur,
    # surgery attributes
    surgery_type = if_else(!is.na(surgery_date),
                           sample_cat(profile$surgery_type, n()), NA_character_),
    surgery_class = case_when(
      surgery_type == "oesophagectomy" ~ "oesophagectomy",
      surgery_type %in% c("total_gastrectomy","partial_gastrectomy") ~ "gastrectomy",
      TRUE ~ NA_character_),
    curative_surgery = if_else(!is.na(surgery_date), TRUE, NA),
    opcs_primary = if_else(!is.na(surgery_date), "G021", NA_character_),
    PROCODE3 = if_else(!is.na(surgery_date) | !is.na(emresd_date), treat_trust, NA_character_),
    SITETRET = if_else(!is.na(PROCODE3), draw_hosp(treat_trust), NA_character_),
    ORGCODEPROVIDER = if_else(!is.na(rt_date), treat_trust, NA_character_),
    # SACT attributes
    BENCHMARK_GROUP = if_else(!is.na(sact_date),
                              sample_cat(profile$benchmark_group, n()), NA_character_),
    benchmark_group_lwr = tolower(BENCHMARK_GROUP),
    INTENT_OF_TREATMENT_V3 = if_else(!is.na(sact_date),
                                     sample_cat(profile$sact_intent, n()), NA_character_),
    CHEMO_RADIATION = if_else(!is.na(sact_date),
                              sample_cat(profile$chemo_radiation, n()), NA_character_),
    ORGANISATION_CODE_OF_PROVIDER = if_else(!is.na(sact_date), treat_trust, NA_character_),
    emergency_admission = as.integer(as.character(route_combined) == "Emergency presentation")
  )

# =============================================================================
# 3. Re-derive flags / pathway / first_tx_date / tx_trust / waits (REAL logic)
# =============================================================================
build_pre_cwt <- function(df) {
  df %>% mutate(
    had_emresd           = !is.na(emresd_date),
    had_surgery          = !is.na(surgery_date),
    had_curative_surgery = !is.na(surgery_date) & curative_surgery == TRUE,
    had_sact             = !is.na(sact_date),
    had_rt               = !is.na(rt_date),
    had_curative_rt      = !is.na(rt_date) & rt_curative == TRUE,
    had_palliative_rt    = !is.na(rt_date) & rt_curative == FALSE,
    sact_before_surgery  = had_sact & had_surgery & sact_date < surgery_date,
    sact_after_surgery   = had_sact & had_surgery & sact_date > surgery_date,
    rt_before_surgery    = had_rt   & had_surgery & rt_date   < surgery_date,
    rt_after_surgery     = had_rt   & had_surgery & rt_date   > surgery_date,
    concurrent_chemo_rt  = had_sact & had_curative_rt &
      abs(as.integer(sact_date - rt_date)) <= 14,
    received_curative_tx = had_emresd | had_curative_surgery | had_curative_rt,
    tx_pathway = case_when(
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
    ),
    first_tx_date = case_when(
      tx_pathway %in% c("EMR/ESD only", "EMR/ESD then surgery") ~ emresd_date,
      tx_pathway == "Surgery + neoadjuvant chemoRT"            ~ pmin(sact_date, rt_date, na.rm = TRUE),
      tx_pathway == "Surgery + neoadjuvant RT"                 ~ rt_date,
      tx_pathway == "Surgery + neoadjuvant chemo"              ~ sact_date,
      tx_pathway %in% c("Surgery + adjuvant chemo","Surgery only","Surgery + other") ~ surgery_date,
      tx_pathway == "Definitive chemoRT"                       ~ pmin(sact_date, rt_date, na.rm = TRUE),
      tx_pathway == "Curative RT only"                         ~ rt_date,
      TRUE                                                     ~ as.Date(NA)
    ),
    tx_trust = case_when(
      tx_pathway %in% c("EMR/ESD only","EMR/ESD then surgery","Surgery + neoadjuvant chemo",
                        "Surgery + adjuvant chemo","Surgery only","Surgery + other") ~ substr(PROCODE3,1,3),
      tx_pathway %in% c("Surgery + neoadjuvant chemoRT","Surgery + neoadjuvant RT") ~ substr(ORGCODEPROVIDER,1,3),
      tx_pathway %in% c("Definitive chemoRT","Curative RT only")                    ~ substr(ORGCODEPROVIDER,1,3),
      TRUE                                                                          ~ NA_character_
    ),
    # trust change between diagnosis and curative treatment
    change_trust = substr(diag_trust, 1, 3) != tx_trust,
    days_endo_to_dx   = as.integer(diagmdy - endoscopy_date),
    days_dx_to_emresd = as.integer(emresd_date - diagmdy),
    days_dx_to_surg   = as.integer(surgery_date - diagmdy),
    days_dx_to_sact   = as.integer(sact_date - diagmdy),
    days_dx_to_rt     = as.integer(rt_date - diagmdy),
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
    wt_rt_to_surg   = as.integer(surgery_date - rt_date)
    )
}

syn_cohort <- build_pre_cwt(cohort)

# QC: did the build re-derive the assigned pathway?
cat("Pathway re-derivation match:",
    round(100 * mean(syn_cohort$tx_pathway == syn_cohort$tx_pathway_asg), 1), "%\n")
syn_cohort %>% filter(tx_pathway != tx_pathway_asg) %>%
  count(tx_pathway_asg, tx_pathway, sort = TRUE) %>% print(n = 20)

spec_cols <- c(
  "pseudo_patientid","pseudo_tumourid","diagmdy","ydiag","tumour_site_grp",
  "cancer_subtype","stage_clean","morphology_num","sex","agediag",
  "ethnicity_group_broad","NHSE_reversed_imd_quintile_lsoas",
  "canalliance_2024_code","canalliance_2024_name","lsoa11_code",
  "route_combined",  "emergency_admission","ps_num","cnsinvolved",
  "diag_hosp","diag_trust","first_trust",
  "endoscopy_date","emresd_date","surgery_date","sact_date","rt_date","first_tx_date",
  "surgery_type","surgery_class","curative_surgery","opcs_primary","PROCODE3","SITETRET",
  "rt_curative","rt_dose","rt_fractions","ORGCODEPROVIDER",
  "BENCHMARK_GROUP","benchmark_group_lwr","INTENT_OF_TREATMENT_V3","CHEMO_RADIATION",
  "ORGANISATION_CODE_OF_PROVIDER","days_endo_to_dx","days_dx_to_emresd","days_dx_to_surg",
  "days_dx_to_sact","days_dx_to_rt","had_emresd","had_surgery","had_curative_surgery",
  "had_sact","had_rt","had_curative_rt","had_palliative_rt","sact_before_surgery",
  "sact_after_surgery","rt_before_surgery","rt_after_surgery","concurrent_chemo_rt",
  "received_curative_tx","tx_pathway","tx_trust","change_trust",
  "wt_dx_to_tx","wt_endo_to_tx","wt_dx_to_surg","wt_endo_to_surg","wt_dx_to_sact",
  "wt_endo_to_sact","wt_sact_to_surg","wt_surg_to_sact","wt_dx_to_rt","wt_endo_to_rt",
  "wt_rt_to_surg","finmdy","died"
)
syn_cohort_out <- syn_cohort %>% select(any_of(spec_cols))
saveRDS(syn_cohort_out, paste0(base_dir, "og_cohort_precwt_SYNTH.rds") )

# =============================================================================
# 4. CWT generator (anchored to the cohort so the merge reproduces validation)
# =============================================================================
site_c15 <- profile$cwt_site_icd10 %>% filter(str_starts(level, "C15"))
site_c16 <- profile$cwt_site_icd10 %>% filter(str_starts(level, "C16"))
mod_by_pw <- profile$cwt_modality_by_pathway
g_vs   <- profile$cwt_glue_by_pathway$cwt_vs_first_tx
g_off  <- profile$cwt_glue_by_pathway$dtt_to_cwt_treat
g_ddt  <- profile$cwt_glue_by_pathway$days_dx_to_dtt
g_exact <- profile$cwt_agreement_by_pathway   # tx_pathway, pct_exact, ...
ov_off <- profile$cwt_glue_overall$dtt_to_cwt_treat
ov_ddt <- profile$cwt_glue_overall$days_dx_to_dtt
recs   <- profile$cwt_records_per_patient %>% filter(!is.na(prop), prop > 0)


fmt <- function(d) ifelse(is.na(d), NA_character_, format(d, "%d/%m/%Y"))

gen_cwt_one <- function(pid, pw, site_grp, dx, ftx) {
  if (runif(1) >= profile$cwt_coverage$pct_any_cwt) return(NULL)  # no CWT record
  k <- sample(recs$records, 1, prob = recs$prop / sum(recs$prop))
  
  # anchored (earliest) record
  if (!is.na(ftx)) {
    pe <- g_exact$pct_exact[g_exact$tx_pathway == pw]
    pe <- if (length(pe) == 1) pe else 0
    if (runif(1) < pe) {
      delta <- 0L                                   # exact-agreement point mass
    } else {
      repeat {                                       # non-zero spread only
        delta <- qv5(iv5(g_vs, profile$cwt_glue_overall$cwt_vs_first_tx, pw), 1)
        if (delta != 0L) break
      }
    }
    treat <- ftx + delta
    off   <- pmax(qv5(iv5(g_off, ov_off, pw), 1), 0L)
    dtt   <- treat - off
  } else {
    ddt   <- qv5(iv5(g_ddt, ov_ddt, pw), 1)
    dtt   <- dx + ddt
    off   <- pmax(qv5(iv5(g_off, ov_off, pw), 1), 0L)
    treat <- dtt + off
  }
  # keep anchored DTT inside the window so it stays the merge anchor
  dtt   <- as.Date(pmin(pmax(as.integer(dtt - dx), -30L), tx_window_days), origin = dx)
  treat <- pmax(treat, dtt)
  
  mb <- mod_by_pw %>% filter(tx_pathway == pw)
  mod_draw <- function(nn) if (nrow(mb)) sample(mb$modality, nn, TRUE, mb$prop/sum(mb$prop))
  else sample(profile$cwt_modality$level, nn, TRUE,
              profile$cwt_modality$prop)
  st_tbl <- if (site_grp == "gastric") site_c16 else site_c15
  site   <- sample(st_tbl$level, 1, TRUE, st_tbl$prop / sum(st_tbl$prop))
  
  has_mdt <- runif(1) < profile$cwt_completeness$pct_mdt
  mdt <- if (has_mdt) dtt - qov(profile$mdt_to_dtt, 1) else as.Date(NA)
  crtp <- dtt - sample(20:60, 1)
  fseen <- crtp + sample(0:14, 1)
  
  rows <- tibble(
    pseudo_patientid = pid, site_icd10 = site, modality = mod_draw(1),
    crtp_date = fmt(crtp), date_first_seen = fmt(fseen), mdt_date = fmt(mdt),
    treat_period_start = fmt(dtt), treat_start = fmt(treat)
  )
  # extra (non-anchor) records with later DTTs
  if (k > 1) {
    gaps <- cumsum(sample(14:120, k - 1, TRUE))
    dtt2 <- dtt + gaps
    rows <- bind_rows(rows, tibble(
      pseudo_patientid = pid, site_icd10 = site, modality = mod_draw(k - 1),
      crtp_date = fmt(dtt2 - sample(20:60, k-1, TRUE)),
      date_first_seen = NA_character_, mdt_date = NA_character_,
      treat_period_start = fmt(dtt2),
      treat_start = fmt(dtt2 + pmax(qv5(iv5(g_off, ov_off, pw), k-1), 0L))
    ))
  }
  rows
}

syn_cwt <- pmap_dfr(
  list(syn_cohort$pseudo_patientid, syn_cohort$tx_pathway,
       syn_cohort$tumour_site_grp, syn_cohort$diagmdy, syn_cohort$first_tx_date),
  gen_cwt_one
)

saveRDS(syn_cwt, paste0(base_dir, "cwt_records_SYNTH.rds"))
cat("\nSynthetic CWT rows:", nrow(syn_cwt),
    "| patients with >=1 record:", n_distinct(syn_cwt$pseudo_patientid), "\n")

# =============================================================================
# 5. MERGE QC - run the real merge on synthetic data, compare to targets
# =============================================================================
cwt_og <- syn_cwt %>%
  mutate(cwt_dtt_date   = as.Date(treat_period_start, "%d/%m/%Y"),
         cwt_treat_date = as.Date(treat_start,        "%d/%m/%Y"),
         cwt_mdt_date   = as.Date(mdt_date,           "%d/%m/%Y"))

cwt_anchor <- syn_cohort %>% select(pseudo_patientid, diagmdy) %>%
  left_join(cwt_og %>%
              filter(!(modality %in% c("23","24") & cwt_treat_date < as.Date("2020-06-01")),
                     !modality %in% c("97","98","99")) %>%
              select(pseudo_patientid, cwt_dtt_date, cwt_treat_date, cwt_mdt_date, modality),
            by = "pseudo_patientid") %>%
  mutate(days_dx_to_dtt = as.integer(cwt_dtt_date - diagmdy)) %>%
  filter(!is.na(days_dx_to_dtt), days_dx_to_dtt >= -30, days_dx_to_dtt <= tx_window_days) %>%
  arrange(pseudo_patientid, cwt_dtt_date) %>%
  distinct(pseudo_patientid, .keep_all = TRUE)

val <- cwt_anchor %>%
  left_join(syn_cohort %>% select(pseudo_patientid, first_tx_date, tx_pathway),
            by = "pseudo_patientid") %>%
  mutate(dtt_to_cwt_treat = as.integer(cwt_treat_date - cwt_dtt_date),
         dtt_to_tx        = as.integer(first_tx_date - cwt_dtt_date),
         cwt_vs_first_tx  = as.integer(cwt_treat_date - first_tx_date))

cat("\n--- MERGE QC (synthetic vs real targets) ---\n")
cat("CWT anchor patients (real ~36,197):", nrow(cwt_anchor), "\n")
cat("days_dx_to_dtt median (real 39):",
    median(cwt_anchor$days_dx_to_dtt), "\n")
cat("dtt_to_cwt_treat median (real 11):",
    median(val$dtt_to_cwt_treat, na.rm = TRUE), "\n")
v <- val %>% filter(!is.na(cwt_vs_first_tx))
cat("cwt_vs_first_tx exact-match (real 71.1%):",
    round(100 * mean(v$cwt_vs_first_tx == 0), 1), "%\n")
cat("within 14d (real 85.6%):",
    round(100 * mean(abs(v$cwt_vs_first_tx) <= 14), 1), "%\n")
cat("negative dtt_to_tx (real 5.3%):",
    round(100 * mean(val$dtt_to_tx < 0, na.rm = TRUE), 1), "%\n\n")
cat("Negative dtt_to_tx by pathway (real: EMR/ESD then surgery ~50%):\n")
val %>% filter(!is.na(dtt_to_tx)) %>% group_by(tx_pathway) %>%
  summarise(n = n(), pct_neg = round(100*mean(dtt_to_tx < 0),1), .groups="drop") %>%
  arrange(desc(pct_neg)) %>% print()

# Save files for stata.
syn_cohort_out_stata <- syn_cohort_out %>% mutate(across(where(is.factor), as.character))
write_dta(syn_cohort_out_stata, paste0(base_dir, "og_cohort_precwt_SYNTH.dta"))
write_dta(syn_cwt, paste0(base_dir,  "cwt_records_SYNTH.dta"))