# =============================================================================
# OG cancer - generate minimal synthetic data (Table A, Table B, merged cohort)
# -----------------------------------------------------------------------------
# Builds a synthetic OG cohort that reproduces the audit treatment pathways and
# the CWT merge, off the disclosure-safe profile. Runs anywhere - no real data.
#
#   Table A  synthetic registry + treatment cohort (one row per patient), with
#            the treatment anchor dates, derived tx_pathway and first_tx_date.
#            This is the LEFT side of the merge and is realistic on its own.
#   Table B  synthetic raw CWT records, modality consistent with each patient's
#            pathway, dates as dd/mm/yyyy character (the RIGHT side).
#   merged   Table A passed through the shared og_cwt_merge() - identical engine
#            to the one the real condensed cohort uses.
#
# Source 01_og_minimal_merge.R first (provides og_cwt_merge + constants/lookups).
# =============================================================================

library(tidyverse)
library(haven)

# the shared merge engine and OG lookups
source("R/build_synthetic_data/01_og_minimal_merge.R")

base_dir <- "Data/synthetic/"
set.seed(2045)
N_TOTAL  <- 40000L

prof_path <- paste0(base_dir, "og_profile_for_synthetic.rds")
spec_path <- paste0(base_dir, "og_minimal_spec.rds")
profile <- if (file.exists(prof_path)) readRDS(prof_path) else NULL
spec    <- if (file.exists(spec_path)) readRDS(spec_path) else NULL
mc      <- og_merge_const   # from the sourced merge engine

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
# Stage 1 derivation: build tx_pathway, first_tx_date and tx_trust from the raw
# treatment fields. The latent intended pathway is dropped here - everything
# downstream uses the DERIVED pathway, exactly as the real pipeline does.
# =============================================================================
A_raw <- A %>% select(any_of(og_raw_cols))
A_der <- og_derive_pathway(A_raw)

# how often does the derived pathway match the latent intended one? (a check
# that the generator's date placement and the derivation are consistent)
cat("\nDerived vs latent pathway agreement:",
    round(100 * mean(A_der$tx_pathway == A$tx_pathway_latent), 1), "%\n")

# =============================================================================
# Table B - synthetic raw CWT records, modality consistent with the pathway
# =============================================================================
cov     <- num_or(profile$cwt_coverage$pct_any_cwt, defaults$cwt_coverage)
p_exact <- num_or(profile$cwt_agreement$pct_exact,  defaults$pct_exact)
p_w14   <- num_or(profile$cwt_agreement$pct_within_14, defaults$pct_within_14)
p_mdt   <- num_or(profile$cwt_completeness$pct_mdt, defaults$pct_mdt)
recs    <- if (!is.null(profile$cwt_records_per_patient))
  profile$cwt_records_per_patient %>% transmute(level = records, prop) else
    defaults$records_per_patient

# modality of the anchored CWT record, by pathway (profile or sensible default)
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
    idx <- which(pw_vec == pw)
    tbl <- if (!is.null(mod_by_pathway))
      mod_by_pathway %>% filter(tx_pathway == pw, !is.na(prop)) else NULL
    if (!is.null(tbl) && nrow(tbl) >= 1 && sum(tbl$prop) > 0)
      out[idx] <- sample(tbl$cwt_modality, length(idx), TRUE, tbl$prop / sum(tbl$prop))
    else
      out[idx] <- default_mod_for_pathway(pw)
  }
  out
}

# which CWT clock-stop date does the anchor sit on? for curative pathways it is
# first_tx_date; for palliative/none it is the earliest recorded treatment date.
# All fields come from the DERIVED cohort.
A_der$cwt_event_date <- A_der$first_tx_date
i <- is.na(A_der$cwt_event_date)
A_der$cwt_event_date[i] <- pmin(A_der$sact_date[i], A_der$rt_date[i],
                                A_der$surgery_date[i], A_der$emresd_date[i],
                                na.rm = TRUE)

has_cwt <- runif(n) < cov & !is.na(A_der$cwt_event_date)
idx     <- which(has_cwt)
m       <- length(idx)

# DTT precedes the treatment date by the dx->dtt vs dx->tx gap; approximate by
# placing DTT a short interval before the event date
dtt_lead <- rgamma_ms(m, 14, 10)
dtt_anchor <- A_der$cwt_event_date[idx] - dtt_lead

# agreement offset between cwt_treat_date and first_tx_date
u <- runif(m); delta <- integer(m)
mid <- u >= p_exact & u < p_w14; far <- u >= p_w14
delta[mid] <- sample(c(-14:-1, 1:14), sum(mid), TRUE)
delta[far] <- sample(c(-60:-15, 15:60), sum(far), TRUE)
treat_anchor <- A_der$cwt_event_date[idx] + delta

mdt_have <- runif(m) < p_mdt
mdt_anchor <- as.Date(rep(NA, m), origin = "1970-01-01")
mdt_anchor[mdt_have] <- dtt_anchor[mdt_have] - sample(0:21, sum(mdt_have), TRUE)

crtp  <- dtt_anchor - sample(20:60, m, TRUE)
fseen <- crtp + sample(0:14, m, TRUE)
site  <- ifelse(A_der$tumour_site_grp[idx] == "gastric",
                sample(c("C160","C161","C162","C163","C164","C165","C166","C169"), m, TRUE),
                sample(c("C150","C151","C152","C153","C154","C155","C159"), m, TRUE))
modal <- draw_modality(A_der$tx_pathway[idx])

anchor <- tibble(
  pseudo_patientid   = A_der$pseudo_patientid[idx],
  site_icd10         = site,
  modality           = modal,
  crtp_date          = fmt(crtp),
  date_first_seen    = fmt(fseen),
  mdt_date           = fmt(mdt_anchor),
  treat_period_start = fmt(dtt_anchor),
  treat_start        = fmt(treat_anchor)
)

# extra non-anchor records for a minority (subsequent treatments / noise),
# placed later so the merge's in-window + pathway-consistency picks the anchor
k <- as.integer(sample_cat(recs, m))
extra_idx <- which(k > 1)
extra <- map_dfr(extra_idx, function(j) {
  mm <- k[j] - 1L
  base_dtt <- A_der$cwt_event_date[idx[j]] + cumsum(sample(30:120, mm, TRUE))
  tibble(
    pseudo_patientid   = A_der$pseudo_patientid[idx[j]],
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
# Stage 2: run the shared merge on the derived cohort, then report
# =============================================================================
syn_cohort <- og_cwt_merge(A_der, syn_cwt)

cat("\nSynthetic pathway mix (derived):\n")
A_der %>% count(tx_pathway) %>% mutate(pct = round(100*n/sum(n),1)) %>%
  arrange(desc(n)) %>% print(n = 20)

cat("\nAudit Table 4 (synthetic, stage 1-3):\n")
syn_cohort %>%
  filter(stage_clean %in% c("1","2","3")) %>%
  summarise(pct_curative  = round(100*mean(received_curative_tx_audit, na.rm = TRUE)),
            pct_any_tx     = round(100*mean(received_any_tx, na.rm = TRUE))) %>%
  print()
if (!is.null(profile$audit_targets))
  cat("  (real targets: curative",
      round(100*profile$audit_targets$pct_curative), "any",
      round(100*profile$audit_targets$pct_any_tx), ")\n")

cat("\nCWT coverage among synthetic cohort:",
    round(100*mean(!is.na(syn_cohort$cwt_dtt_date)),1), "%\n")
cat("dtt_valid TRUE share (non-EMR pathways):",
    round(100*mean(syn_cohort$dtt_valid, na.rm = TRUE),1), "%\n")

# =============================================================================
# Save the three cohort stages + the CWT records
#   _raw     the raw registry+treatment inputs (dates, descriptors, providers)
#   _derived the cohort after og_derive_pathway() (flags, tx_pathway, tx_trust)
#   _cohort  the merged analysis cohort after og_cwt_merge()
# =============================================================================
saveRDS(A_raw,      paste0(base_dir, "og_ncras_treatment_synthetic.rds"))
saveRDS(A_der,      paste0(base_dir, "og_derived_synthetic.rds"))
saveRDS(syn_cwt,    paste0(base_dir, "og_cwt_records_synthetic.rds"))
saveRDS(syn_cohort, paste0(base_dir, "og_cohort_synthetic.rds"))

to_stata <- function(df) df %>%
  mutate(across(where(is.factor), as.character),
         across(where(is.logical), as.integer))
write_dta(to_stata(A_raw),      paste0(base_dir, "og_ncras_treatment_synthetic.dta"))
write_dta(to_stata(A_der),      paste0(base_dir, "og_derived_synthetic.dta"))
write_dta(to_stata(syn_cwt),    paste0(base_dir, "og_cwt_records_synthetic.dta"))
write_dta(to_stata(syn_cohort), paste0(base_dir, "og_cohort_synthetic.dta"))

cat("\nSaved raw Table A (", nrow(A_raw), "rows), derived cohort (",
    nrow(A_der), "rows), Table B CWT (", nrow(syn_cwt),
    "rows), merged cohort (", nrow(syn_cohort), "rows).\n")
