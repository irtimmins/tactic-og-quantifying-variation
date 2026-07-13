# =============================================================================
# 09b  Finding the diagnosing site in CWT (exploration only)
# -----------------------------------------------------------------------------
# Question this answers: across ALL of a patient's CWT rows - not just the one
# treatment row we anchor on - is there an earlier event that sits at the
# DIAGNOSING hospital, before the referral to a specialist centre moves the
# recorded organisation elsewhere?
#
# Each CWT row carries several dated events, each with its own organisation:
#   first-seen   date_first_seen   / org_first_seen    (first outpatient visit)
#   referral     crtp_date         / (no org)          (referral into the period)
#   decision     treat_period_start/ org_dec_to_treat  (decision to treat)
#   treatment    treat_start       / org_treat_start   (treatment starts)
# The pathway owner org_ppi (trust level) is also on the row.
#
# We reshape CWT so there is one row per (patient, event), work out how far each
# event is from the registry diagnosis date, drop implausible ones, and then ask
# for each event type: how close to diagnosis does it sit, how complete is it,
# and how often does its organisation fall in the registry's diagnosing trust
# (diag_trust). The event that is both close to diagnosis and usually in the
# right trust is the best source for a diagnosing-site variable.
#
# This is exploration: it prints tables and writes nothing to the cohort.
#
# Reads : Data/ICON/og_cohort_2015_2022.rds (cohort after 07), the CWT dataset
# Writes: outputs/cwt_diag_site_event_summary.csv  (small aggregate table)
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

library(dplyr)
library(tidyr)
library(stringr)

first3   <- function(x) str_sub(x, 1, 3)
is_site5 <- function(x) !is.na(x) & nchar(x) == 5

# where to write the summary. 01 may call the outputs folder something other
# than dir_out, so fall back gracefully rather than error at the end.
out_dir <- if (exists("dir_out")) dir_out else
  if (exists("dir_icon")) dir_icon else "."

# -----------------------------------------------------------------------------
# 1. Inputs. We keep the registry diagnosis date and diagnosing trust to anchor
#    proximity and trust checks.
# -----------------------------------------------------------------------------
og_cohort <- readRDS(f_cohort)
reg <- og_cohort %>%
  transmute(pseudo_patientid = as.character(pseudo_patientid),
            diagmdy,
            diag_trust = na_if(str_trim(as.character(diag_trust)), ""),
            diag_hosp  = na_if(str_trim(as.character(diag_hosp)),  ""))
ncras_og_ids <- unique(reg$pseudo_patientid)

if (!exists("read_cwt"))
  read_cwt <- function() open_dataset(path_cwt_partition) %>%
  filter(site_icd10 %in% og_icd10) %>% collect()

# read every CWT row (not just the anchor), parse the event dates and keep the
# organisation attached to each event.
cwt_raw <- read_cwt() %>%
  mutate(pseudo_patientid = as.character(pseudo_patientid)) %>%
  filter(pseudo_patientid %in% ncras_og_ids) %>%
  transmute(
    pseudo_patientid,
    d_first_seen = as.Date(date_first_seen,     "%d/%m/%Y"),
    d_referral   = as.Date(crtp_date,           "%d/%m/%Y"),
    d_decision   = as.Date(treat_period_start,  "%d/%m/%Y"),
    d_treatment  = as.Date(treat_start,         "%d/%m/%Y"),
    o_first_seen = na_if(str_trim(org_first_seen),   ""),
    o_decision   = na_if(str_trim(org_dec_to_treat), ""),
    o_treatment  = na_if(str_trim(org_treat_start),  ""),
    o_ppi        = na_if(str_trim(org_ppi),          ""))

cat("CWT rows read:", nrow(cwt_raw), "over",
    n_distinct(cwt_raw$pseudo_patientid), "patients\n")
cat("(so about", round(nrow(cwt_raw) / n_distinct(cwt_raw$pseudo_patientid), 1),
    "CWT rows per patient on average)\n\n")

# -----------------------------------------------------------------------------
# 2. Reshape to one row per (patient, event). Each event has a date and, where
#    the event carries one, an organisation. Referral has a date but no org, so
#    it contributes to the timing picture only.
# -----------------------------------------------------------------------------
events <- bind_rows(
  cwt_raw %>% transmute(pseudo_patientid, event = "first_seen",
                        edate = d_first_seen, eorg = o_first_seen),
  cwt_raw %>% transmute(pseudo_patientid, event = "referral",
                        edate = d_referral,  eorg = NA_character_),
  cwt_raw %>% transmute(pseudo_patientid, event = "decision",
                        edate = d_decision,  eorg = o_decision),
  cwt_raw %>% transmute(pseudo_patientid, event = "treatment",
                        edate = d_treatment, eorg = o_treatment)) %>%
  filter(!is.na(edate)) %>%
  left_join(reg, by = "pseudo_patientid") %>%
  mutate(days_from_dx = as.integer(edate - diagmdy))

# -----------------------------------------------------------------------------
# 3. Drop implausible dates. A diagnosing-hospital event should sit near the
#    diagnosis, not years away and not far before it. Keep a generous window
#    around diagnosis; the tails are data errors or unrelated episodes.
#    (Window is deliberately wide so we can see the distribution before tightening.)
# -----------------------------------------------------------------------------
lo <- -60      # up to ~2 months before the registry diagnosis date
hi <- 365      # up to a year after
plausible <- events %>% filter(days_from_dx >= lo, days_from_dx <= hi)

cat("events before plausibility filter:", nrow(events), "\n")
cat("events kept (", lo, "to", hi, "days from diagnosis):", nrow(plausible), "\n\n")

# -----------------------------------------------------------------------------
# 4. For each event type: how close to diagnosis does it sit, how complete is
#    its organisation, and how often is that organisation in the registry trust?
#    This is the table that tells us which event to use.
# -----------------------------------------------------------------------------
event_summary <- plausible %>%
  group_by(event) %>%
  summarise(
    n_events        = n(),
    median_days_dx  = median(days_from_dx),
    p25_days_dx     = quantile(days_from_dx, .25),
    p75_days_dx     = quantile(days_from_dx, .75),
    has_org         = mean(!is.na(eorg)),
    org_is_site5    = mean(is_site5(eorg), na.rm = TRUE),
    in_diag_trust   = mean(first3(eorg) == diag_trust, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(median_days_dx)

cat("Event timing and agreement with the registry diagnosing trust\n")
cat(strrep("-", 78), "\n")
print(as.data.frame(event_summary %>%
                      mutate(across(c(has_org, org_is_site5, in_diag_trust), ~ round(100 * .x, 1)),
                             across(c(median_days_dx, p25_days_dx, p75_days_dx), as.integer))))
cat("\nRead: the event with a small median distance from diagnosis AND a high\n")
cat("in_diag_trust is the best diagnosing-site source. has_org shows how many\n")
cat("patients it would actually cover.\n\n")

# -----------------------------------------------------------------------------
# 5. The refinement you asked for: instead of one event, take the EARLIEST
#    plausible organisation-bearing CWT event per patient, then apply the
#    diagnosing trust to accept or repair it. Compare that against just using
#    the decision (treatment-anchor) org.
#
#    earliest_org  : org from the earliest plausible dated event that has an org.
#    Coalescing rule: prefer the earliest event's org, but if that org is not in
#    the registry trust and a LATER plausible event's org IS, take the in-trust
#    one - the aim is the code that best matches where the patient was diagnosed,
#    and the registry trust is the referee.
# -----------------------------------------------------------------------------
with_org <- plausible %>% filter(!is.na(eorg))

# earliest org-bearing event per patient
earliest <- with_org %>%
  group_by(pseudo_patientid) %>%
  arrange(edate, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(pseudo_patientid, earliest_org = eorg,
            earliest_days = days_from_dx)

# best in-trust org among a patient's plausible events (any event), earliest first.
# with_org already carries diag_trust from the reg join above, so no re-join.
in_trust_pick <- with_org %>%
  filter(first3(eorg) == diag_trust) %>%
  group_by(pseudo_patientid) %>%
  arrange(edate, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(pseudo_patientid, intrust_org = eorg)

# coalesced diagnosing site: earliest org, but swap in an in-trust org where the
# earliest one is out of trust and an in-trust one exists.
combined <- reg %>%
  left_join(earliest,      by = "pseudo_patientid") %>%
  left_join(in_trust_pick, by = "pseudo_patientid") %>%
  mutate(
    earliest_in_trust = !is.na(earliest_org) & first3(earliest_org) == diag_trust,
    cwt_diag_site_v2  = case_when(
      earliest_in_trust                    ~ earliest_org,   # earliest is already right trust
      !is.na(intrust_org)                  ~ intrust_org,    # a later in-trust event rescues it
      !is.na(earliest_org)                 ~ earliest_org,   # keep earliest (out of trust, flagged below)
      TRUE                                 ~ NA_character_))

# how good is the resulting variable?
cat("Refined diagnosing site: earliest plausible CWT org, trust-refereed\n")
cat(strrep("-", 78), "\n")
have <- combined %>% filter(!is.na(cwt_diag_site_v2))
cat(sprintf("patients with a derived diagnosing site : %d of %d (%.1f%%)\n",
            nrow(have), nrow(combined), 100 * nrow(have) / nrow(combined)))
cat(sprintf("of those, organisation in registry trust: %.1f%%\n",
            100 * mean(first3(have$cwt_diag_site_v2) == have$diag_trust, na.rm = TRUE)))
cat(sprintf("of those, exact site match to diag_hosp  : %.1f%% (where diag_hosp present)\n",
            100 * mean(have$cwt_diag_site_v2 == have$diag_hosp, na.rm = TRUE)))
cat(sprintf("resulting variable is a full 5-char site : %.1f%%\n",
            100 * mean(is_site5(have$cwt_diag_site_v2))))

# for contrast, the same three numbers for the plain decision (anchor) org
dec_only <- reg %>%
  left_join(with_org %>% filter(event == "decision") %>%
              group_by(pseudo_patientid) %>% arrange(edate) %>% slice(1) %>% ungroup() %>%
              transmute(pseudo_patientid, dec_org = eorg),
            by = "pseudo_patientid") %>%
  filter(!is.na(dec_org))
cat("\nfor contrast, the plain decision-to-treat org:\n")
cat(sprintf("  in registry trust: %.1f%%   exact to diag_hosp: %.1f%%\n",
            100 * mean(first3(dec_only$dec_org) == dec_only$diag_trust, na.rm = TRUE),
            100 * mean(dec_only$dec_org == dec_only$diag_hosp, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# 6. Save the small event summary (aggregate counts and rates only).
# -----------------------------------------------------------------------------
out_csv <- file.path(out_dir, "cwt_diag_site_event_summary.csv")
write.csv(event_summary, out_csv, row.names = FALSE)
cat("\nSaved event summary:", out_csv, "\n")
cat("09b exploration complete. No cohort files were changed.\n")