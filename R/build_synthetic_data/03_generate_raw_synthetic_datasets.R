# =============================================================================
# 03  Generate the raw synthetic datasets
# -----------------------------------------------------------------------------
# Generates the two raw synthetic inputs and stops - it does NOT derive the
# pathway or merge (those are scripts 04 and 05). Runs anywhere; no real data.
#
#   raw cohort   the synthetic registry + treatment cohort: identity, patient and
#                tumour descriptors, the treatment DATES, the curative descriptors
#                (curative_surgery, rt_curative), chemo provenance, and the
#                per-modality provider codes. No pathway - that is derived in 04.
#   CWT records  the synthetic raw CWT records, one row per treatment event,
#                dates as dd/mm/yyyy character.
#
# How it works: a latent "intended" pathway is drawn per patient (from the
# profile, or built-in defaults) purely to place a realistic set of treatment
# dates and to choose a consistent CWT modality. The latent label is then
# dropped - script 04 re-derives the pathway from the dates alone. The latent
# pathway is also saved on its own as a QC file so script 04 can confirm the
# derivation recovers what was intended.
#
# Reads (optional): Data/synthetic/og_profile_for_synthetic.rds  (from 02)
# Writes: Data/synthetic/og_ncras_treatment_synthetic.rds/.dta   (raw cohort)
#         Data/synthetic/og_cwt_records_synthetic.rds/.dta        (CWT records)
#         Data/synthetic/og_intended_pathway_qc.rds               (QC only)
# =============================================================================

library(tidyverse)
library(haven)

# paths, relative to the project root (the .Rproj working directory)
dir_fns <- "R/build_synthetic_data"   # location of 01_define_functions.R
dir_syn <- "Data/synthetic"           # synthetic inputs and outputs
dir.create(dir_syn, recursive = TRUE, showWarnings = FALSE)

# load the constants and lookups (tx_pathway_levels, og_raw_cols, og_merge_const)
source(file.path(dir_fns, "01_define_functions.R"))

set.seed(2045)
N_TOTAL  <- 40000L

prof_path <- file.path(dir_syn, "og_profile_for_synthetic.rds")
spec_path <- file.path(dir_syn, "og_minimal_spec.rds")
profile <- if (file.exists(prof_path)) readRDS(prof_path) else NULL
spec    <- if (file.exists(spec_path)) readRDS(spec_path) else NULL
mc      <- og_merge_const   # from 01

# -----------------------------------------------------------------------------
# Helpers - sampling with profile-or-default fallback
# -----------------------------------------------------------------------------
def_marg <- function(levels, props) tibble(level = as.character(levels),
                                           prop = props / sum(props))
sample_cat <- function(tbl, n) {
  p <- tbl$prop; p[is.na(p)] <- 0
  sample(tbl$level, n, replace = TRUE, prob = p / sum(p))
}
num_or <- function(x, d) if (is.null(x) || length(x) == 0 || is.na(x)) d else x
marg <- function(key, default) {
  m <- tryCatch(profile$marginals[[key]], error = function(e) NULL)
  if (is.null(m) || !nrow(m)) default else m %>% transmute(level, prop)
}
interval_ms <- function(key, m_def, s_def) {
  q <- tryCatch(profile$intervals_overall[[key]], error = function(e) NULL)
  if (is.null(q)) c(m_def, s_def) else c(num_or(q$mean, m_def), num_or(q$sd, s_def))
}
rgamma_ms <- function(n, m, s) {       # gamma with target mean/sd, >= 0 integer
  m <- max(m, 1); s <- max(s, 1)
  shape <- (m / s)^2; rate <- m / s^2
  pmax(0L, as.integer(round(rgamma(n, shape = shape, rate = rate))))
}
fmt <- function(d) ifelse(is.na(d), NA_character_, format(d, "%d/%m/%Y"))

# -----------------------------------------------------------------------------
# Defaults (used if the profile is absent or a section is missing)
# -----------------------------------------------------------------------------
defaults <- list(
  by_year   = def_marg(2015:2022, rep(1, 8)),
  sex       = def_marg(c(1,2), c(.74,.26)),
  age_grp   = def_marg(c("<50","50-54","55-59","60-64","65-69","70-74",
                         "75-79","80-84","85+"),
                       c(.03,.04,.07,.12,.17,.20,.18,.12,.07)),
  ethnicity = def_marg(c("White","Asian","Black","Other","Unknown"),
                       c(.88,.04,.02,.02,.04)),
  imd       = def_marg(1:5, c(.18,.19,.20,.21,.22)),
  site_grp  = def_marg(c("oesophageal","gastric"), c(.65,.35)),
  subtype   = def_marg(c("Oes ACA","Oes SCC","Gast"), c(.42,.21,.37)),
  stage     = def_marg(c("1","2","3"), c(.18,.28,.54)),
  route     = def_marg(c("TWW","GP referral","Emergency presentation",
                         "Other outpatient","Inpatient elective","Screening",
                         "Unknown"),
                       c(.39,.23,.15,.11,.09,.005,.025)),
  cci       = def_marg(c("0","1","2","3+"), c(.55,.26,.12,.07)),
  modality  = def_marg(c("01","02","04","05","23","07","08","09","24"),
                       c(.31,.31,.08,.15,.08,.02,.02,.005,.005)),
  records_per_patient = def_marg(1:4, c(.48,.32,.14,.06)),
  cwt_coverage = .895, pct_exact = .73, pct_within_14 = .87,
  pct_dtt = 1.0, pct_mdt = .39,
  n_trust = 140L, n_hosp = 220L, sd_hosp_dtt = 8, sd_hosp_tx = 6
)

# =============================================================================
# Trust / hospital geography (shared by registry and CWT org codes)
# =============================================================================
n_trust <- as.integer(num_or(profile$volume$diag_trust$n_distinct, defaults$n_trust))
n_hosp  <- max(n_trust, as.integer(num_or(profile$volume$diag_hosp$n_distinct, defaults$n_hosp)))
all_2char   <- as.vector(outer(LETTERS, c(0:9, LETTERS), paste0))
trust_codes <- paste0("R", sample(all_2char, n_trust))
trust_w     <- sort(rlnorm(n_trust, 1, 1), decreasing = TRUE); trust_w <- trust_w / sum(trust_w)
hosp_per_trust <- rep(1L, n_trust)
extra <- n_hosp - n_trust
if (extra > 0)
  hosp_per_trust <- hosp_per_trust +
  tabulate(sample(seq_len(n_trust), extra, TRUE, prob = trust_w), nbins = n_trust)
hosp_codes <- unlist(mapply(function(t, k) paste0(t, sample(all_2char, k)),
                            trust_codes, hosp_per_trust, SIMPLIFY = FALSE))
hosp_trust <- rep(trust_codes, times = hosp_per_trust)
draw_trust <- function(n) sample(trust_codes, n, TRUE, trust_w)
draw_hosp  <- function(trusts) {
  out <- character(length(trusts))
  for (tr in unique(trusts)) {
    idx <- which(trusts == tr); h <- hosp_codes[hosp_trust == tr]
    out[idx] <- sample(h, length(idx), TRUE)
  }
  out
}
re_dtt <- setNames(rnorm(length(hosp_codes), 0,
                         num_or(profile$between_hosp_sd$wt_dx_to_dtt, defaults$sd_hosp_dtt)),
                   hosp_codes)

# =============================================================================
# Table A - registry covariates
# =============================================================================
n <- N_TOTAL
A <- tibble(pseudo_patientid = sprintf("SYN%07d", seq_len(n)))

yr_tbl <- if (!is.null(profile$by_year))
  profile$by_year %>% transmute(level, prop = prop / sum(prop)) else defaults$by_year
A$ydiag <- as.integer(as.character(sample_cat(yr_tbl, n)))
doy <- as.integer(runif(n, 0, 364))
A$diagmdy <- as.Date(paste0(A$ydiag, "-01-01")) + doy

age_grp <- sample_cat(marg("age_grp", defaults$age_grp), n)
lo <- c("<50"=40,"50-54"=50,"55-59"=55,"60-64"=60,"65-69"=65,"70-74"=70,"75-79"=75,"80-84"=80,"85+"=85)
hi <- c("<50"=49,"50-54"=54,"55-59"=59,"60-64"=64,"65-69"=69,"70-74"=74,"75-79"=79,"80-84"=84,"85+"=95)
A$agediag <- round(runif(n, lo[age_grp], hi[age_grp]) + runif(n), 2)

A$sex                              <- as.integer(sample_cat(marg("sex", defaults$sex), n))
A$ethnicity_group_broad            <- sample_cat(marg("ethnicity", defaults$ethnicity), n)
A$NHSE_reversed_imd_quintile_lsoas <- as.character(sample_cat(marg("imd_quintile", defaults$imd), n))
A$tumour_site_grp                  <- sample_cat(marg("site_grp", defaults$site_grp), n)
A$stage_clean                      <- as.character(sample_cat(marg("stage", defaults$stage), n))
A$route_combined                   <- sample_cat(marg("route_combined", defaults$route), n)
A$cci_group                        <- as.character(sample_cat(marg("cci_group", defaults$cci), n))

# subtype consistent with site: gastric site -> "Gast"; oesophageal -> SCC/ACA
sub_oes <- def_marg(c("Oes ACA","Oes SCC"), c(.67,.33))
A$cancer_subtype <- if_else(
  A$tumour_site_grp == "gastric", "Gast",
  sample_cat(sub_oes, n))

# diagnosing trust / hospital
A$diag_trust <- draw_trust(n)
A$diag_hosp  <- draw_hosp(A$diag_trust)

# =============================================================================
# Pathway assignment - conditional on stage x subtype (the audit core)
# Uses profile$pathway_by_stage_subtype when present, else a clinically-shaped
# default that puts neoadjuvant+surgery in earlier stage, palliative/none in
# stage 3 and the elderly.
# =============================================================================
pbss <- profile$pathway_by_stage_subtype

default_pathway_props <- function(stage, subtype) {
  # baseline mix; shifted by stage. props need not sum to 1 (normalised later)
  base <- c(
    "Surgery + neoadjuvant chemo"   = 19,
    "No treatment recorded"         = 25,
    "SACT only"                     = 11,
    "Palliative RT only"            = 9,
    "Definitive chemoRT"            = 9,
    "Surgery only"                  = 8,
    "EMR/ESD only"                  = 6,
    "Palliative chemo + RT"         = 6,
    "Curative RT only"              = 2,
    "Surgery + neoadjuvant chemoRT" = 3,
    "Surgery + adjuvant chemo"      = 1,
    "EMR/ESD then surgery"          = 1,
    "Surgery + neoadjuvant RT"      = 0.2,
    "Surgery + other"               = 0.1)
  if (stage == "1") {           # more curative, more EMR/ESD, less palliative
    base["EMR/ESD only"] <- 22; base["Surgery only"] <- 14
    base["No treatment recorded"] <- 16; base["SACT only"] <- 5
    base["Palliative RT only"] <- 4; base["Palliative chemo + RT"] <- 2
  } else if (stage == "3") {    # less surgery, more palliative / none
    base["Surgery + neoadjuvant chemo"] <- 16; base["No treatment recorded"] <- 28
    base["SACT only"] <- 13; base["Palliative RT only"] <- 11
    base["EMR/ESD only"] <- 2
  }
  if (identical(subtype, "Oes SCC")) {   # SCC -> more definitive chemoRT, less surgery
    base["Definitive chemoRT"] <- base["Definitive chemoRT"] * 2.5
    base["Surgery + neoadjuvant chemo"] <- base["Surgery + neoadjuvant chemo"] * 0.6
  }
  if (identical(subtype, "Gast")) {      # gastric -> more surgery, less chemoRT
    base["Surgery only"] <- base["Surgery only"] * 1.8
    base["Definitive chemoRT"] <- base["Definitive chemoRT"] * 0.2
  }
  base
}

assign_pathway <- function(stage, subtype) {
  if (!is.null(pbss)) {
    sub <- pbss %>% filter(stage_clean == stage, subtype == !!subtype, !is.na(prop))
    if (nrow(sub) >= 1 && sum(sub$prop) > 0)
      return(sample(sub$tx_pathway, 1, prob = sub$prop / sum(sub$prop)))
  }
  p <- default_pathway_props(stage, subtype)
  sample(names(p), 1, prob = p / sum(p))
}

# vectorise over the (stage, subtype) cells for speed. This is the LATENT
# intended pathway - it drives date placement but is dropped before output; the
# real tx_pathway is re-derived from the dates/flags by og_derive_pathway().
A$tx_pathway_latent <- {
  key <- paste(A$stage_clean, coalesce(A$cancer_subtype, "Unknown"))
  out <- character(n)
  for (k in unique(key)) {
    idx <- which(key == k)
    parts <- strsplit(k, " ", fixed = TRUE)[[1]]
    st <- parts[1]; sb <- paste(parts[-1], collapse = " ")
    if (!is.null(pbss)) {
      sub <- pbss %>% filter(stage_clean == st, subtype == sb, !is.na(prop))
      if (nrow(sub) >= 1 && sum(sub$prop) > 0) {
        out[idx] <- sample(sub$tx_pathway, length(idx), TRUE, sub$prop / sum(sub$prop))
        next
      }
    }
    p <- default_pathway_props(st, sb)
    out[idx] <- sample(names(p), length(idx), TRUE, p / sum(p))
  }
  factor(out, levels = tx_pathway_levels) %>% as.character()
}

# =============================================================================
# Treatment dates and descriptors, placed off diagnosis per the latent pathway
# -----------------------------------------------------------------------------
# The intended pathway (p) is a latent draw used only to place a realistic set
# of treatment dates, curative descriptors and provider codes. It is NOT stored
# on the cohort: og_derive_pathway() re-derives tx_pathway from the dates and
# flags alone, so the generator and the derivation are kept honest against each
# other. The realised raw fields are what the rest of the pipeline sees.
# =============================================================================
ms_endo <- interval_ms("days_dx_to_endo", 2, 5)     # endoscopy ~ at/just before dx
ms_dtt  <- interval_ms("wt_dx_to_dtt", 42, 30)
neoadj_gap_ms <- {
  q <- profile$neoadj_to_surg
  if (is.null(q)) c(100, 30) else c(num_or(q$mean,100), num_or(q$sd,30))
}

# endoscopy date: a few days before/around diagnosis
A$endoscopy_date <- A$diagmdy - rgamma_ms(n, ms_endo[1], ms_endo[2])

# per-pathway dx -> primary-treatment offset (clock-stop timing)
dx_to_primary <- rgamma_ms(n, ms_dtt[1] + 14, ms_dtt[2]) +
  round(re_dtt[A$diag_hosp])
dx_to_primary <- pmax(0L, pmin(dx_to_primary, mc$tx_window_days))

# raw treatment fields: dates, curative descriptors, chemo provenance, providers
A$emresd_date  <- as.Date(NA); A$surgery_date <- as.Date(NA)
A$sact_date    <- as.Date(NA); A$rt_date      <- as.Date(NA)
A$curative_surgery <- NA; A$rt_curative <- NA; A$chemo_source <- NA_character_
A$hes_chemo_date   <- as.Date(NA)
A$surgery_provider <- NA_character_; A$rt_provider <- NA_character_
A$sact_provider    <- NA_character_

p <- A$tx_pathway_latent   # the latent intended pathway (dropped before output)
prim_date  <- A$diagmdy + dx_to_primary
neoadj_gap <- rgamma_ms(n, neoadj_gap_ms[1], neoadj_gap_ms[2])

# EMR/ESD pathways: emresd_date is the clock-stop
i <- p %in% c("EMR/ESD only","EMR/ESD then surgery")
A$emresd_date[i] <- prim_date[i]
j <- i & p == "EMR/ESD then surgery"
A$surgery_date[j] <- prim_date[j] + neoadj_gap[j]
A$curative_surgery[j] <- TRUE

# neoadjuvant chemo: SACT before surgery
i <- p == "Surgery + neoadjuvant chemo"
A$sact_date[i] <- prim_date[i]; A$surgery_date[i] <- prim_date[i] + neoadj_gap[i]
A$curative_surgery[i] <- TRUE; A$chemo_source[i] <- "sact"

# neoadjuvant chemoRT: sact + rt before surgery (rt close to sact)
i <- p == "Surgery + neoadjuvant chemoRT"
A$sact_date[i] <- prim_date[i]; A$rt_date[i] <- prim_date[i] + sample(0:10, sum(i), TRUE)
A$surgery_date[i] <- prim_date[i] + neoadj_gap[i]
A$curative_surgery[i] <- TRUE; A$rt_curative[i] <- TRUE; A$chemo_source[i] <- "sact"

# neoadjuvant RT: rt before surgery, no chemo
i <- p == "Surgery + neoadjuvant RT"
A$rt_date[i] <- prim_date[i]; A$surgery_date[i] <- prim_date[i] + neoadj_gap[i]
A$curative_surgery[i] <- TRUE; A$rt_curative[i] <- TRUE

# adjuvant chemo: surgery first, chemo after
i <- p == "Surgery + adjuvant chemo"
A$surgery_date[i] <- prim_date[i]; A$curative_surgery[i] <- TRUE
A$sact_date[i] <- prim_date[i] + sample(40:80, sum(i), TRUE); A$chemo_source[i] <- "sact"

# surgery only / surgery + other: surgery is the only curative act
i <- p %in% c("Surgery only","Surgery + other")
A$surgery_date[i] <- prim_date[i]; A$curative_surgery[i] <- TRUE

# definitive chemoRT: concurrent sact + curative rt, no surgery
i <- p == "Definitive chemoRT"
A$sact_date[i] <- prim_date[i]; A$rt_date[i] <- prim_date[i] + sample(0:14, sum(i), TRUE)
A$rt_curative[i] <- TRUE; A$chemo_source[i] <- "sact"

# curative RT only: curative rt, no chemo, no surgery
i <- p == "Curative RT only"
A$rt_date[i] <- prim_date[i]; A$rt_curative[i] <- TRUE

# palliative chemo + RT: chemo + palliative rt, no surgery, no curative clock-stop
i <- p == "Palliative chemo + RT"
A$sact_date[i] <- prim_date[i]; A$rt_date[i] <- prim_date[i] + sample(0:30, sum(i), TRUE)
A$rt_curative[i] <- FALSE; A$chemo_source[i] <- "sact"

# SACT only: chemo, no rt, no surgery
i <- p == "SACT only"
A$sact_date[i] <- prim_date[i]; A$chemo_source[i] <- "sact"

# palliative RT only: palliative rt, no chemo, no surgery
i <- p == "Palliative RT only"
A$rt_date[i] <- prim_date[i]; A$rt_curative[i] <- FALSE

# "No treatment recorded": no treatment dates at all

# provider codes for whichever arms occurred (3-char trust embedded in code 1)
A$surgery_provider[!is.na(A$surgery_date)] <- draw_trust(sum(!is.na(A$surgery_date)))
A$surgery_provider[!is.na(A$emresd_date) & is.na(A$surgery_provider)] <-
  draw_trust(sum(!is.na(A$emresd_date) & is.na(A$surgery_provider)))
A$rt_provider[!is.na(A$rt_date)]     <- draw_trust(sum(!is.na(A$rt_date)))
A$sact_provider[!is.na(A$sact_date)] <- draw_trust(sum(!is.na(A$sact_date)))

# =============================================================================
# Survival + death
# =============================================================================
ms_surv <- interval_ms("surv_from_dx_days", 850, 600)
surv_days <- rgamma_ms(n, ms_surv[1], ms_surv[2])
# worse survival for no-treatment / palliative, better for curative
surv_mult <- case_when(
  p == "No treatment recorded"                 ~ 0.45,
  p %in% c("SACT only","Palliative RT only",
           "Palliative chemo + RT")            ~ 0.65,
  TRUE                                          ~ 1.15)
surv_days <- pmax(1L, as.integer(surv_days * surv_mult))
A$died   <- as.integer(runif(n) < pmin(0.95, 0.55 + 0.0002 * surv_days * (surv_mult < 1)))
A$finmdy <- A$diagmdy + surv_days

# =============================================================================
# CWT records (Table B) - built from the latent pathway and the raw dates
# -----------------------------------------------------------------------------
# No derivation is needed here: the latent intended pathway (a generation
# parameter) chooses a consistent CWT modality, and the earliest actual
# treatment date anchors the record. Script 05 will re-match modality to the
# derived pathway during the merge.
# =============================================================================
cov     <- num_or(profile$cwt_coverage$pct_any_cwt, defaults$cwt_coverage)
p_exact <- num_or(profile$cwt_agreement$pct_exact,  defaults$pct_exact)
p_w14   <- num_or(profile$cwt_agreement$pct_within_14, defaults$pct_within_14)
p_mdt   <- num_or(profile$cwt_completeness$pct_mdt, defaults$pct_mdt)
recs    <- if (!is.null(profile$cwt_records_per_patient))
  profile$cwt_records_per_patient %>% transmute(level = records, prop) else
    defaults$records_per_patient

mod_by_pathway <- profile$cwt_modality_by_pathway
default_mod_for_pathway <- function(pw) {
  switch(pw,
         "EMR/ESD only" = "23", "EMR/ESD then surgery" = "01",
         "Surgery + neoadjuvant chemo" = "02",
         "Surgery + neoadjuvant chemoRT" = "04",
         "Surgery + neoadjuvant RT" = "05",
         "Surgery + adjuvant chemo" = "01", "Surgery only" = "01",
         "Surgery + other" = "01",
         "Definitive chemoRT" = "04", "Curative RT only" = "05",
         "Palliative chemo + RT" = "02", "SACT only" = "02",
         "Palliative RT only" = "05", "No treatment recorded" = "07", "01")
}
draw_modality <- function(pw_vec) {
  out <- character(length(pw_vec))
  for (pw in unique(pw_vec)) {
    ix <- which(pw_vec == pw)
    tbl <- if (!is.null(mod_by_pathway))
      mod_by_pathway %>% filter(tx_pathway == pw, !is.na(prop)) else NULL
    if (!is.null(tbl) && nrow(tbl) >= 1 && sum(tbl$prop) > 0)
      out[ix] <- sample(tbl$cwt_modality, length(ix), TRUE, tbl$prop / sum(tbl$prop))
    else
      out[ix] <- default_mod_for_pathway(pw)
  }
  out
}

# the CWT event date is the earliest actual treatment the patient received,
# computed straight from the raw dates (no derivation)
cwt_event_date <- pmin(A$emresd_date, A$surgery_date, A$sact_date, A$rt_date,
                       na.rm = TRUE)

has_cwt <- runif(n) < cov & !is.na(cwt_event_date)
idx     <- which(has_cwt)
m       <- length(idx)

# DTT precedes the treatment date by a short interval
dtt_lead   <- rgamma_ms(m, 14, 10)
dtt_anchor <- cwt_event_date[idx] - dtt_lead

# offset between cwt_treat_date and the event date, reproducing the linkage
# agreement seen in the real data (exact / within-14 / further out)
u <- runif(m); delta <- integer(m)
mid <- u >= p_exact & u < p_w14; far <- u >= p_w14
delta[mid] <- sample(c(-14:-1, 1:14), sum(mid), TRUE)
delta[far] <- sample(c(-60:-15, 15:60), sum(far), TRUE)
treat_anchor <- cwt_event_date[idx] + delta

mdt_have   <- runif(m) < p_mdt
mdt_anchor <- as.Date(rep(NA, m), origin = "1970-01-01")
mdt_anchor[mdt_have] <- dtt_anchor[mdt_have] - sample(0:21, sum(mdt_have), TRUE)

crtp  <- dtt_anchor - sample(20:60, m, TRUE)
fseen <- crtp + sample(0:14, m, TRUE)
site  <- ifelse(A$tumour_site_grp[idx] == "gastric",
                sample(c("C160","C161","C162","C163","C164","C165","C166","C169"), m, TRUE),
                sample(c("C150","C151","C152","C153","C154","C155","C159"), m, TRUE))
modal <- draw_modality(A$tx_pathway_latent[idx])

anchor <- tibble(
  pseudo_patientid   = A$pseudo_patientid[idx],
  site_icd10         = site,
  modality           = modal,
  crtp_date          = fmt(crtp),
  date_first_seen    = fmt(fseen),
  mdt_date           = fmt(mdt_anchor),
  treat_period_start = fmt(dtt_anchor),
  treat_start        = fmt(treat_anchor)
)

# a minority get extra, later (non-anchor) records - subsequent treatments and
# noise the merge's in-window + pathway-consistency rules should screen out
k <- as.integer(sample_cat(recs, m))
extra_idx <- which(k > 1)
extra <- map_dfr(extra_idx, function(j) {
  mm <- k[j] - 1L
  base_dtt <- cwt_event_date[idx[j]] + cumsum(sample(30:120, mm, TRUE))
  tibble(
    pseudo_patientid   = A$pseudo_patientid[idx[j]],
    site_icd10         = anchor$site_icd10[j],
    modality           = sample(c("02","05","07","08","09"), mm, TRUE),
    crtp_date          = fmt(base_dtt - sample(20:60, mm, TRUE)),
    date_first_seen    = NA_character_,
    mdt_date           = NA_character_,
    treat_period_start = fmt(base_dtt),
    treat_start        = fmt(base_dtt + sample(10:40, mm, TRUE)))
})
syn_cwt <- bind_rows(anchor, extra) %>% arrange(pseudo_patientid)

# =============================================================================
# Save the two raw inputs (+ the intended-pathway QC file)
# =============================================================================
A_raw <- A %>% select(any_of(og_raw_cols))
qc_intended <- A %>% transmute(pseudo_patientid,
                               tx_pathway_intended = tx_pathway_latent)

cat("\nIntended (latent) pathway mix:\n")
qc_intended %>% count(tx_pathway_intended) %>%
  mutate(pct = round(100*n/sum(n),1)) %>% arrange(desc(n)) %>% print(n = 20)
cat("\nRaw cohort:", nrow(A_raw), "patients,", ncol(A_raw), "columns",
    "| CWT records:", nrow(syn_cwt), "\n")

saveRDS(A_raw,       file.path(dir_syn, "og_ncras_treatment_synthetic.rds"))
saveRDS(syn_cwt,     file.path(dir_syn, "og_cwt_records_synthetic.rds"))
saveRDS(qc_intended, file.path(dir_syn, "og_intended_pathway_qc.rds"))

to_stata <- function(df) df %>%
  mutate(across(where(is.factor), as.character),
         across(where(is.logical), as.integer))
write_dta(to_stata(A_raw),   file.path(dir_syn, "og_ncras_treatment_synthetic.dta"))
write_dta(to_stata(syn_cwt), file.path(dir_syn, "og_cwt_records_synthetic.dta"))

cat("Saved raw cohort + CWT records (and the intended-pathway QC file).\n",
    "Next: 04_derive_pathway.R\n")