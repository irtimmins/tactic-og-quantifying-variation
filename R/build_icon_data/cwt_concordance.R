# =============================================================================
# 09  CWT site concordance (assessment only)
# -----------------------------------------------------------------------------
# Derives two site variables from the CWT record and checks how well they agree
# with the site/provider codes already held in the cohort from other sources.
# This is an assessment script: it prints agreement rates and writes a small
# summary, but it does NOT change the cohort. Merging these variables into the
# main build (script 08) is a later, separate decision.
#
# The derived variables:
#   cwt_treat_site    the organisation recorded against the start of treatment
#                     (org_treat_start) on the anchor row - the treating site.
#   cwt_diag_site_v2  the DIAGNOSING site, built from the earliest plausible CWT
#                     event across all of a patient's rows (not just the anchor),
#                     with the registry diagnosing trust (diag_trust) as referee.
#                     This is the recommended diagnosing-site variable; the raw
#                     decision-to-treat and first-seen orgs are also reported as
#                     the single-source starting points it improves on.
#
# The treating-site variable is taken from the same single CWT row that script 08
# anchors on, so it lines up with the treatment date already in use. The
# diagnosing-site variable deliberately does NOT use the anchor row, because the
# anchor is the treatment clock-stop and its organisation has usually drifted to
# the specialist centre by then; the diagnosing hospital is better recovered from
# the earliest event in the pathway.
#
# What we compare them against, and why the comparison is split:
#   cwt_diag_site  vs  diag_hosp   (registry diagnosing site)
#   cwt_treat_site vs  SITETRET                        for STRAIGHT-TO-SURGERY patients
#                                                       only (tx_pathway "Surgery only")
#                      ORGANISATION_CODE_OF_PROVIDER   for NEOADJUVANT SACT only (chemo
#                                                       given before surgery, marked by
#                                                       sact_before_surgery)
#                      ORGCODEPROVIDER                 for radiotherapy      (RTDS)
# Each treatment source only covers the patients who had that treatment, so each
# comparison is restricted to the relevant group.
#
# The surgery comparison is restricted to "Surgery only" patients, not everyone
# who had surgery. For a neoadjuvant-then-surgery patient the anchored CWT row
# is the treatment-PERIOD start, which is the chemo or radiotherapy episode, so
# cwt_treat_site would be the SACT/RT site rather than the surgery site, and
# comparing it to SITETRET (the surgery site) would compare two different places
# and understate agreement. Straight-to-surgery patients have no upstream
# episode to displace the anchor, so their cwt_treat_site is the surgery site.
#
# The SACT comparison is deliberately restricted to sact_before_surgery rather
# than all SACT patients. SACT given after surgery, or on its own, is often
# delivered locally rather than at the specialist centre that did the surgery,
# so including it would understate agreement for reasons that have nothing to
# do with data quality - a mismatch there can be a genuine difference in care
# setting, not an error. Neoadjuvant SACT is given as part of the same curative
# pathway as the anchored surgery, so it is the fair comparison.
#
# A note on code formats, which drives how agreement is measured:
#   CWT org codes are 5-character site codes (e.g. R0A02).
#   diag_hosp, SITETRET and the SACT provider code are MIXED - mostly 5-character
#   sites, but some 3-character trust codes.
#   The RTDS code (ORGCODEPROVIDER) is 3-character trust only.
# So a plain string match understates real agreement: a patient can be at the
# same trust but recorded once at site level and once at trust level. We report
# agreement at two levels:
#   exact   full-string match (only fair where both sides are full site codes)
#   trust   first three characters match (the common denominator, and the only
#           level at which RTDS can be compared at all)
# The trust-level rate is the one to lean on for a like-for-like read.
#
# Reads : Data/ICON/og_cohort_2015_2022.rds (cohort after 07), the CWT dataset
# Writes: cwt_site_concordance_summary.csv    (agreement rates per comparison)
#         cwt_site_derived_performance.csv    (coverage + accuracy of the final
#                                              derived variables)
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

library(dplyr)
library(tidyr)
library(stringr)

# how close to the registry diagnosis date a CWT event must sit to be a
# plausible diagnosing-hospital event. Wide enough to catch the real diagnosing
# visit, tight enough to drop unrelated later episodes and date errors.
diag_window_lo <- -60    # up to ~2 months before the registry diagnosis date
diag_window_hi <- 365    # up to a year after

# -----------------------------------------------------------------------------
# 1. Rebuild the CWT anchor exactly as script 08 does, then keep the two org
#    codes off that same anchored row.
# -----------------------------------------------------------------------------
# This mirrors 08's logic so the site codes correspond to the chosen DTT/treat
# record. If 08's anchoring changes, this block should be kept in step with it.

modality_group <- tribble(
  ~modality, ~mod_group,
  "01", "surgery", "23", "surgery", "24", "surgery",
  "02", "chemo", "14", "chemo", "15", "chemo", "03", "hormone",
  "04", "chemort",
  "05", "radiotherapy", "06", "radiotherapy", "13", "radiotherapy",
  "07", "palliative", "08", "palliative", "09", "palliative",
  "97", "other", "98", "declined")

pathway_groups <- list(
  "EMR/ESD only"                  = c("surgery", "other"),
  "EMR/ESD then surgery"          = c("surgery"),
  "Surgery + neoadjuvant chemoRT" = c("surgery", "chemort", "chemo", "radiotherapy"),
  "Surgery + neoadjuvant chemo"   = c("surgery", "chemo"),
  "Surgery + neoadjuvant RT"      = c("surgery", "radiotherapy", "chemort"),
  "Surgery + adjuvant chemo"      = c("surgery", "chemo"),
  "Surgery only"                  = c("surgery"),
  "Surgery + other"               = c("surgery", "other"),
  "Definitive chemoRT"            = c("chemort", "chemo", "radiotherapy"),
  "Curative RT only"              = c("radiotherapy", "chemort"),
  "Palliative chemo + RT"         = c("chemo", "radiotherapy", "chemort", "palliative"),
  "SACT only"                     = c("chemo", "hormone", "palliative"),
  "Palliative RT only"            = c("radiotherapy", "palliative"),
  "No treatment recorded"         = c("palliative", "other"))

pathway_primary <- c(
  "EMR/ESD only" = "surgery", "EMR/ESD then surgery" = "surgery",
  "Surgery + neoadjuvant chemoRT" = "chemort", "Surgery + neoadjuvant chemo" = "chemo",
  "Surgery + neoadjuvant RT" = "radiotherapy", "Surgery + adjuvant chemo" = "surgery",
  "Surgery only" = "surgery", "Surgery + other" = "surgery",
  "Definitive chemoRT" = "chemort", "Curative RT only" = "radiotherapy",
  "Palliative chemo + RT" = "chemo", "SACT only" = "chemo",
  "Palliative RT only" = "radiotherapy", "No treatment recorded" = "palliative")

og_cohort    <- readRDS(f_cohort)
ncras_og_ids <- unique(as.character(og_cohort$pseudo_patientid))

if (!exists("read_cwt"))
  read_cwt <- function() open_dataset(path_cwt_partition) %>%
  filter(site_icd10 %in% og_icd10) %>% collect()

# read CWT, keeping the two org columns we need alongside the anchoring fields
cwt_og <- read_cwt() %>%
  mutate(pseudo_patientid = as.character(pseudo_patientid),
         cwt_dtt_date     = as.Date(treat_period_start, "%d/%m/%Y"),
         cwt_treat_date   = as.Date(treat_start,        "%d/%m/%Y")) %>%
  filter(pseudo_patientid %in% ncras_og_ids)

cwt_grouped <- cwt_og %>%
  left_join(modality_group, by = "modality") %>%
  mutate(mod_group = case_when(
    modality == "01" & surg_01_rule == "transition_window" &
      cwt_treat_date > surg_transition_end                                 ~ NA_character_,
    modality %in% c("23","24") & surg_01_rule == "transition_window" &
      cwt_treat_date < surg_transition_start                               ~ NA_character_,
    modality == "01" & surg_01_rule == "date_split" &
      cwt_treat_date >= surg_switch_date                                   ~ NA_character_,
    modality %in% c("23","24") & surg_01_rule == "date_split" &
      cwt_treat_date < surg_switch_date                                    ~ NA_character_,
    modality == "01" & surg_01_rule == "never"                            ~ NA_character_,
    TRUE                                                                  ~ mod_group)) %>%
  filter(!is.na(mod_group), mod_group != "declined", !is.na(cwt_dtt_date))

pathway_group_long <- enframe(pathway_groups, name = "tx_pathway", value = "ok_group") %>%
  unnest_longer(ok_group)
pw <- og_cohort %>% select(pseudo_patientid, diagmdy, tx_pathway, first_tx_date)

cwt_candidates <- cwt_grouped %>%
  inner_join(pw, by = "pseudo_patientid") %>%
  mutate(days_dx_to_dtt = as.integer(cwt_dtt_date - diagmdy)) %>%
  filter(days_dx_to_dtt >= dtt_min_offset, days_dx_to_dtt <= cwt_window_days) %>%
  left_join(pathway_group_long %>% mutate(group_ok = TRUE),
            by = c("tx_pathway", "mod_group" = "ok_group")) %>%
  mutate(group_ok   = coalesce(group_ok, FALSE),
         is_primary = !is.na(mod_group) &
           mod_group == coalesce(unname(pathway_primary[tx_pathway]), "")) %>%
  group_by(pseudo_patientid) %>%
  mutate(any_match = coalesce(any(group_ok), FALSE)) %>%
  filter(if (isTRUE(first(any_match))) group_ok else TRUE) %>%
  ungroup()

# the anchor row, carrying every org code we might use. org_dec_to_treat is the
# treatment-period-start org (always present); org_first_seen is where the
# patient was first seen, which sits earlier in the pathway and may line up
# better with the registry's DIAGNOSING hospital; org_ppi is the pathway-owning
# provider (trust level).
cwt_anchor <- cwt_candidates %>%
  group_by(pseudo_patientid) %>%
  arrange(desc(is_primary), cwt_dtt_date, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(pseudo_patientid,
            cwt_diag_site      = na_if(str_trim(org_dec_to_treat), ""),
            cwt_treat_site     = na_if(str_trim(org_treat_start),  ""),
            cwt_first_seen_org = na_if(str_trim(org_first_seen),   ""),
            cwt_ppi_org        = na_if(str_trim(org_ppi),          ""))

cat("CWT anchor rows with a decision-to-treat org:", sum(!is.na(cwt_anchor$cwt_diag_site)),
    "of", nrow(cwt_anchor), "\n")
cat("CWT anchor rows with a treating-site org:    ", sum(!is.na(cwt_anchor$cwt_treat_site)),
    "of", nrow(cwt_anchor), "\n")
cat("CWT anchor rows with a first-seen org:       ", sum(!is.na(cwt_anchor$cwt_first_seen_org)),
    "of", nrow(cwt_anchor), "\n\n")

first3   <- function(x) str_sub(x, 1, 3)
is_site5 <- function(x) !is.na(x) & nchar(x) == 5

# -----------------------------------------------------------------------------
# 1b. Build the diagnosing-site variable from the EARLIEST plausible CWT event.
#
# The anchor row is chosen for the treatment clock-stop, so its organisation has
# usually drifted to the specialist centre by the time of the decision to treat.
# The diagnosing hospital is better recovered from the earliest event in the
# pathway. We therefore look across ALL of a patient's CWT rows (not just the
# anchor), treat each dated event (first-seen / decision / treatment) as its own
# record, keep those plausibly near the registry diagnosis date, and take the
# organisation from the earliest one. The registry diagnosing trust (diag_trust)
# then acts as a referee: if the earliest org is not in that trust but a later
# plausible event is, we take the in-trust one instead.
#
# This block reads every CWT row afresh (the anchor above kept only one row per
# patient), so it is self-contained.
# -----------------------------------------------------------------------------
reg_dx <- og_cohort %>%
  transmute(pseudo_patientid = as.character(pseudo_patientid),
            diagmdy,
            diag_trust = na_if(str_trim(as.character(diag_trust)), ""))

cwt_events <- read_cwt() %>%
  mutate(pseudo_patientid = as.character(pseudo_patientid)) %>%
  filter(pseudo_patientid %in% ncras_og_ids) %>%
  transmute(
    pseudo_patientid,
    d_first_seen = as.Date(date_first_seen,    "%d/%m/%Y"),
    d_decision   = as.Date(treat_period_start, "%d/%m/%Y"),
    d_treatment  = as.Date(treat_start,        "%d/%m/%Y"),
    o_first_seen = na_if(str_trim(org_first_seen),   ""),
    o_decision   = na_if(str_trim(org_dec_to_treat), ""),
    o_treatment  = na_if(str_trim(org_treat_start),  ""))

# one row per (patient, dated event) that carries an organisation
events_long <- bind_rows(
  cwt_events %>% transmute(pseudo_patientid, edate = d_first_seen, eorg = o_first_seen),
  cwt_events %>% transmute(pseudo_patientid, edate = d_decision,   eorg = o_decision),
  cwt_events %>% transmute(pseudo_patientid, edate = d_treatment,  eorg = o_treatment)) %>%
  filter(!is.na(edate), !is.na(eorg)) %>%
  left_join(reg_dx, by = "pseudo_patientid") %>%
  mutate(days_from_dx = as.integer(edate - diagmdy)) %>%
  filter(days_from_dx >= diag_window_lo, days_from_dx <= diag_window_hi)  # plausible only

# earliest plausible org-bearing event per patient
earliest_org <- events_long %>%
  group_by(pseudo_patientid) %>%
  arrange(edate, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(pseudo_patientid, earliest_org = eorg)

# earliest event whose org sits in the registry diagnosing trust (the referee)
intrust_org <- events_long %>%
  filter(first3(eorg) == diag_trust) %>%
  group_by(pseudo_patientid) %>%
  arrange(edate, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(pseudo_patientid, intrust_org = eorg)

# the derived diagnosing site: earliest org if it is already in the right trust,
# else a later in-trust org if one exists, else the earliest org (out of trust).
cwt_diag_site_derived <- reg_dx %>%
  left_join(earliest_org, by = "pseudo_patientid") %>%
  left_join(intrust_org,  by = "pseudo_patientid") %>%
  mutate(cwt_diag_site_v2 = case_when(
    !is.na(earliest_org) & first3(earliest_org) == diag_trust ~ earliest_org,
    !is.na(intrust_org)                                       ~ intrust_org,
    !is.na(earliest_org)                                      ~ earliest_org,
    TRUE                                                      ~ NA_character_)) %>%
  select(pseudo_patientid, cwt_diag_site_v2)

# persist the derived diagnosing site (one row per patient) so other assessment
# scripts - notably 10, the HES-vs-CWT comparison - can read it without having
# to rebuild it. This is a derived-variable side table, not the cohort itself.
out_dir_09 <- if (exists("dir_out")) dir_out else
  if (exists("dir_icon")) dir_icon else "."
saveRDS(cwt_diag_site_derived, file.path(out_dir_09, "cwt_diag_site_v2.rds"))
cat("Saved derived diagnosing site:",
    file.path(out_dir_09, "cwt_diag_site_v2.rds"), "\n\n")

# -----------------------------------------------------------------------------
# 2. Attach the derived codes to the cohort (in memory only - not saved).
#    cwt_diag_site / cwt_treat_site come from the anchor row; cwt_diag_site_v2
#    is the earliest-event, trust-refereed diagnosing site built in 1b.
# -----------------------------------------------------------------------------
dat <- og_cohort %>%
  left_join(cwt_anchor,             by = "pseudo_patientid") %>%
  left_join(cwt_diag_site_derived,  by = "pseudo_patientid") %>%
  mutate(across(c(diag_hosp, SITETRET, ORGANISATION_CODE_OF_PROVIDER, ORGCODEPROVIDER),
                ~ na_if(str_trim(as.character(.x)), "")))

# -----------------------------------------------------------------------------
# 3. A helper to score agreement between two code columns.
#
#    Reports several things, because the mixed code lengths mean a single number
#    hides what is going on:
#      n_pair       rows where both sides have a value
#      pct_exact    full-string match across all pairs
#      pct_trust    first-three-character (trust) match across all pairs
#      pct_both5    share of pairs where BOTH sides are full 5-character site
#                   codes (as opposed to a 3-character trust code on one side)
#      exact_in5    site-level agreement WITHIN those full-5-character pairs -
#                   this is the cleanest read of true site concordance, because
#                   it removes the cases where one side only recorded a trust
#    A code is treated as "full site" if it is 5 characters long.
#    (first3 and is_site5 are defined earlier, in the derivation block.)
# -----------------------------------------------------------------------------
concordance <- function(data, cwt_col, ref_col, label) {
  d <- data %>%
    transmute(cwt = .data[[cwt_col]], ref = .data[[ref_col]]) %>%
    filter(!is.na(cwt), !is.na(ref))
  
  n_pair <- nrow(d)
  if (n_pair == 0) {
    cat(sprintf("%-46s  no overlapping records to compare\n", label))
    return(tibble(comparison = label, n_compared = 0L,
                  pct_exact = NA_real_, pct_trust = NA_real_,
                  pct_both_5char = NA_real_, exact_within_5char = NA_real_))
  }
  
  exact <- mean(d$cwt == d$ref)
  trust <- mean(first3(d$cwt) == first3(d$ref))
  
  # subset where both sides carry a full 5-character site code
  both5 <- is_site5(d$cwt) & is_site5(d$ref)
  pct_both5 <- mean(both5)
  exact_in5 <- if (any(both5)) mean(d$cwt[both5] == d$ref[both5]) else NA_real_
  
  cat(sprintf("%-46s  n=%6d\n", label, n_pair))
  cat(sprintf("    all pairs      : exact %5.1f%%   trust %5.1f%%\n",
              100 * exact, 100 * trust))
  cat(sprintf("    both 5-char    : %5.1f%% of pairs, and within those, site agreement %5.1f%%\n",
              100 * pct_both5, 100 * exact_in5))
  
  tibble(comparison = label, n_compared = n_pair,
         pct_exact = round(100 * exact, 1), pct_trust = round(100 * trust, 1),
         pct_both_5char = round(100 * pct_both5, 1),
         exact_within_5char = round(100 * exact_in5, 1))
}

# -----------------------------------------------------------------------------
# 4. Diagnosing-site concordance.
#
# The registry diag_hosp is the DIAGNOSING hospital. In many pathways diagnosis
# (endoscopy) happens at a local hospital and the decision-to-treat is recorded
# later at the specialist centre the patient is referred to - so comparing
# diag_hosp against the decision-to-treat org can show a genuine pathway
# difference rather than a data error. CWT records an earlier point too, the
# org where the patient was first seen (org_first_seen), which should sit closer
# to diagnosis. We compare diag_hosp against both, to see whether the earlier
# point aligns better and would be the better basis for a diagnosing-site
# variable.
# -----------------------------------------------------------------------------
cat("Diagnosing site vs registry diag_hosp\n")
cat(strrep("-", 78), "\n")
res_diag_dtt <- concordance(dat, "cwt_diag_site", "diag_hosp",
                            "raw: decision-to-treat org")
res_diag_fs  <- concordance(dat, "cwt_first_seen_org", "diag_hosp",
                            "raw: first-seen org")
res_diag_v2  <- concordance(dat, "cwt_diag_site_v2", "diag_hosp",
                            "DERIVED: earliest event, trust-refereed")
cat("\nThe derived variable (cwt_diag_site_v2) is the one to judge - it takes the\n")
cat("earliest plausible CWT organisation and uses the registry trust as referee.\n")
cat("The two raw rows above show the single-source starting points it improves on.\n\n")

# -----------------------------------------------------------------------------
# 4b. Stabilising the CWT diagnosing site against the registry trust.
#
# We already hold diag_trust: the registry's own diagnosing TRUST (a reliable
# 3-character code, essentially complete). We can use it as an anchor to
# strengthen the CWT site code.
#
# The idea: when CWT and diag_hosp disagree, that disagreement is one of two
# very different things -
#   wrong trust        the CWT org is a different organisation entirely - a real
#                      mismatch, and the CWT site should not be trusted.
#   right trust,       the CWT org is in the correct trust but carries a
#   different site      different site code (often just which building within
#                      the trust was recorded) - a minor, usually harmless
#                      labelling difference.
# By checking each CWT org against diag_trust we can separate these, and we can
# build a stabilised site variable that only accepts a CWT site when it sits in
# the correct trust, and otherwise falls back to the trust code. This trades a
# little site-level detail for a variable that is never in the wrong trust.
#
# The breakdown below is computed for whichever CWT org agreed best above - by
# default org_first_seen, with org_dec_to_treat shown alongside for comparison.
# -----------------------------------------------------------------------------
trust_anchor_report <- function(data, cwt_col, label) {
  d <- data %>%
    transmute(cwt = .data[[cwt_col]], diag_trust, diag_hosp) %>%
    filter(!is.na(cwt), !is.na(diag_trust))
  
  n <- nrow(d)
  if (n == 0) { cat(sprintf("%-34s  no records\n", label)); return(invisible(NULL)) }
  
  in_trust      <- first3(d$cwt) == d$diag_trust          # CWT org sits in the registry trust
  d_hosp        <- !is.na(d$diag_hosp)
  site_match    <- d_hosp & d$cwt == d$diag_hosp          # full site agreement (where diag_hosp present)
  right_t_wrong_s <- in_trust & d_hosp & !site_match       # right trust, different site
  
  cat(sprintf("%s (n=%d with a registry trust)\n", label, n))
  cat(sprintf("    CWT org sits in the registry trust    : %5.1f%%\n", 100*mean(in_trust)))
  cat(sprintf("    of the rest, wrong trust entirely     : %5.1f%%\n", 100*mean(!in_trust)))
  cat(sprintf("    right trust but a different site code : %5.1f%% (of all)\n", 100*mean(right_t_wrong_s)))
  
  # a stabilised site: accept the CWT site only when it is in the correct trust,
  # otherwise fall back to the (reliable) trust code.
  stabilised <- if_else(in_trust, d$cwt, d$diag_trust)
  covered_site <- mean(is_site5(stabilised))
  cat(sprintf("    stabilised variable is a full site for: %5.1f%% (rest fall back to trust)\n\n",
              100*covered_site))
  invisible(NULL)
}

cat("Stabilising the CWT diagnosing site against the registry trust (diag_trust)\n")
cat(strrep("-", 78), "\n")
trust_anchor_report(dat, "cwt_first_seen_org", "first-seen org vs diag_trust")
trust_anchor_report(dat, "cwt_diag_site",      "decision-to-treat org vs diag_trust")
cat("A high 'sits in the registry trust' figure means most apparent site\n")
cat("disagreement is within the right trust, so the trust-anchored variable is\n")
cat("safe to build. A high 'wrong trust' figure means CWT is genuinely pointing\n")
cat("elsewhere and the registry site should be preferred.\n\n")

# -----------------------------------------------------------------------------
# 5. Treating-site concordance, split by treatment type so each row is compared
#    against the source that actually covers it.
#
#    The surgery comparison uses tx_pathway == "Surgery only" (straight to
#    surgery), and the SACT comparison uses sact_before_surgery, not had_sact -
#    see the note at the top of the file for why each is restricted this way.
# -----------------------------------------------------------------------------
cat("Treating site (cwt_treat_site vs the source covering each treatment)\n")
cat(strrep("-", 78), "\n")

res_surg <- dat %>% filter(tx_pathway == "Surgery only") %>%
  concordance("cwt_treat_site", "SITETRET",
              "straight-to-surgery patients vs HES SITETRET")

res_sact <- dat %>% filter(sact_before_surgery %in% TRUE) %>%
  concordance("cwt_treat_site", "ORGANISATION_CODE_OF_PROVIDER",
              "neoadjuvant SACT patients vs SACT provider code")

res_rt <- dat %>% filter(had_rt %in% TRUE) %>%
  concordance("cwt_treat_site", "ORGCODEPROVIDER",
              "radiotherapy patients vs RTDS provider (trust-level only)")

cat("\nNote: the RTDS code is trust-level (3 characters) only, so its exact\n")
cat("column is not a like-for-like site comparison - read its trust column.\n\n")

# -----------------------------------------------------------------------------
# 6. Performance and coverage of the two derived variables.
#
# This pulls the whole assessment together into one place: for each derived
# variable, how much of the cohort it covers, and - among covered patients - how
# often it lands in the right trust and on the exact registry site. This is the
# table to read when deciding whether a variable is good enough to merge.
#
#   cwt_diag_site_v2  the derived diagnosing site (earliest event, trust-refereed)
#                     judged against diag_hosp / diag_trust across all patients.
#   cwt_treat_site    the treating site from the anchor row, judged against the
#                     source that covers each treatment group, pooled here into a
#                     single treating-site performance line per group.
# -----------------------------------------------------------------------------
# a small scorer: coverage over the whole cohort, then trust/site agreement among
# covered patients that also have a reference code to check against.
perf <- function(data, cwt_col, ref_col, ref_trust_col, label, denom_n) {
  covered <- sum(!is.na(data[[cwt_col]]))
  d <- data %>%
    transmute(cwt = .data[[cwt_col]], ref = .data[[ref_col]],
              ref_trust = .data[[ref_trust_col]]) %>%
    filter(!is.na(cwt))
  in_trust <- mean(first3(d$cwt) == d$ref_trust, na.rm = TRUE)
  exact    <- mean(d$cwt == d$ref, na.rm = TRUE)
  full5    <- mean(is_site5(d$cwt))
  tibble(variable = label,
         coverage_pct    = round(100 * covered / denom_n, 1),
         n_covered       = covered,
         correct_trust_pct = round(100 * in_trust, 1),
         exact_site_pct    = round(100 * exact, 1),
         full_site_pct     = round(100 * full5, 1))
}

n_cohort <- nrow(dat)

# diagnosing site: derived variable, whole cohort, vs registry diag_hosp/diag_trust
perf_diag <- perf(dat, "cwt_diag_site_v2", "diag_hosp", "diag_trust",
                  "cwt_diag_site_v2 (diagnosing)", n_cohort)

# treating site: the anchor org, scored within each treatment group against that
# group's source. diag_trust here is the group's own reference trust (first 3 of
# the reference code), so correct_trust is like-for-like.
perf_treat <- bind_rows(
  dat %>% filter(tx_pathway == "Surgery only") %>%
    mutate(ref_trust = first3(SITETRET)) %>%
    perf("cwt_treat_site", "SITETRET", "ref_trust",
         "cwt_treat_site: straight-to-surgery", n_cohort),
  dat %>% filter(sact_before_surgery %in% TRUE) %>%
    mutate(ref_trust = first3(ORGANISATION_CODE_OF_PROVIDER)) %>%
    perf("cwt_treat_site", "ORGANISATION_CODE_OF_PROVIDER", "ref_trust",
         "cwt_treat_site: neoadjuvant SACT", n_cohort),
  dat %>% filter(had_rt %in% TRUE) %>%
    mutate(ref_trust = first3(ORGCODEPROVIDER)) %>%
    perf("cwt_treat_site", "ORGCODEPROVIDER", "ref_trust",
         "cwt_treat_site: radiotherapy", n_cohort))

perf_tbl <- bind_rows(perf_diag, perf_treat)

cat("Derived variable performance and coverage\n")
cat(strrep("-", 78), "\n")
cat("coverage_pct is of the whole cohort; the agreement columns are among\n")
cat("covered patients who also have a reference code.\n")
cat("(for radiotherapy, exact_site is not meaningful - RTDS is trust-level only)\n\n")
print(as.data.frame(perf_tbl))
cat("\n")

# -----------------------------------------------------------------------------
# 7. Save the aggregate summaries (counts and percentages only, no patient data).
# -----------------------------------------------------------------------------
# outputs folder: 01 may not define dir_out, so fall back rather than error.
out_dir <- if (exists("dir_out")) dir_out else
  if (exists("dir_icon")) dir_icon else "."

concordance_tbl <- bind_rows(
  res_diag_dtt %>% mutate(group = "diagnosing site"),
  res_diag_fs  %>% mutate(group = "diagnosing site"),
  res_diag_v2  %>% mutate(group = "diagnosing site"),
  res_surg     %>% mutate(group = "treating site"),
  res_sact     %>% mutate(group = "treating site"),
  res_rt       %>% mutate(group = "treating site")) %>%
  select(group, comparison, n_compared, pct_exact, pct_trust,
         pct_both_5char, exact_within_5char)

write.csv(concordance_tbl, file.path(out_dir, "cwt_site_concordance_summary.csv"),
          row.names = FALSE)
write.csv(perf_tbl,        file.path(out_dir, "cwt_site_derived_performance.csv"),
          row.names = FALSE)

cat("Saved:\n")
cat("  ", file.path(out_dir, "cwt_site_concordance_summary.csv"), "\n")
cat("  ", file.path(out_dir, "cwt_site_derived_performance.csv"), "\n")
cat("09 assessment complete. No cohort files were changed.\n")