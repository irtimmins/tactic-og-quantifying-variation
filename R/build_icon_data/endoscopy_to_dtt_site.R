# =============================================================================
# 11  Endoscopy-to-DTT analysis cohort: how many have a 5-digit endoscopy site
# -----------------------------------------------------------------------------
# The positive deviance analysis needs three things per patient: an endoscopy
# date, a valid decision-to-treat, and the hospital that did the endoscopy at
# 5-character site resolution. This reports how many patients survive each
# condition and, crucially, the intersection - the actual analysis n.
#
# "Valid endoscopy-to-DTT" is reported under a few definitions, because the
# number depends on how strict you are and it is better to see them side by side
# than to bake one in:
#   has both dates   : an endoscopy date and a CWT DTT date both present
#   interval sensible: the above, and wt_endo_to_dtt is non-negative and not
#                      implausibly long
#   dtt_valid flag   : the cohort's own dtt_valid, which also asks that curative
#                      treatment followed the decision sensibly
#
# Reads : the CWT-merged cohort (needs endoscopy_date, cwt_dtt_date,
#         wt_endo_to_dtt, dtt_valid, sitetret_endoscopy_apc)
# Writes: nothing - prints only
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

library(dplyr)
library(stringr)

is_site5 <- function(x) !is.na(x) & nchar(x) == 5

# upper bound on a sensible endoscopy-to-DTT gap. Endoscopy is the clock start,
# DTT the decision; a gap beyond this is almost certainly a linkage artefact
# rather than a real wait. Tunable.
endo_dtt_max_days <- 365

cohort_path <- if (exists("f_cohort_cwt") && file.exists(f_cohort_cwt)) f_cohort_cwt else f_cohort
dat <- readRDS(cohort_path) %>%
  mutate(sitetret_endoscopy_apc = na_if(str_trim(as.character(sitetret_endoscopy_apc)), ""))

n_cohort <- nrow(dat)
has_site  <- is_site5(dat$sitetret_endoscopy_apc)

# the three validity definitions, from loosest to strictest
has_both  <- !is.na(dat$endoscopy_date) & !is.na(dat$cwt_dtt_date)
sensible  <- has_both & !is.na(dat$wt_endo_to_dtt) &
  dat$wt_endo_to_dtt >= 0 & dat$wt_endo_to_dtt <= endo_dtt_max_days
flag_valid <- if ("dtt_valid" %in% names(dat)) dat$dtt_valid %in% TRUE else rep(NA, n_cohort)

pct <- function(x) round(100 * sum(x, na.rm = TRUE) / n_cohort, 1)

cat("Whole cohort:", n_cohort, "patients\n\n")

cat("Each condition on its own (share of whole cohort)\n")
cat(strrep("-", 68), "\n")
cat(sprintf("  has a 5-digit endoscopy site        : %6d  (%.1f%%)\n", sum(has_site), pct(has_site)))
cat(sprintf("  has both endoscopy and DTT dates    : %6d  (%.1f%%)\n", sum(has_both), pct(has_both)))
cat(sprintf("  ... and interval sensible (0-%dd)   : %6d  (%.1f%%)\n", endo_dtt_max_days, sum(sensible), pct(sensible)))
if ("dtt_valid" %in% names(dat))
  cat(sprintf("  dtt_valid flag TRUE                 : %6d  (%.1f%%)\n",
              sum(flag_valid, na.rm = TRUE), pct(flag_valid)))
cat("\n")

# the intersections that actually matter: valid endoscopy-to-DTT AND a 5-digit site
report_intersection <- function(valid, valid_label) {
  both_conditions <- valid & has_site
  cat(sprintf("%-34s valid: %6d   valid AND 5-digit site: %6d  (%.1f%% of cohort, %.1f%% of the valid group)\n",
              valid_label,
              sum(valid, na.rm = TRUE),
              sum(both_conditions, na.rm = TRUE),
              round(100 * sum(both_conditions, na.rm = TRUE) / n_cohort, 1),
              round(100 * sum(both_conditions, na.rm = TRUE) / sum(valid, na.rm = TRUE), 1)))
}

cat("Analysis n: valid endoscopy-to-DTT AND a 5-digit endoscopy site\n")
cat(strrep("-", 68), "\n")
report_intersection(has_both,  "has both dates")
report_intersection(sensible,  "interval sensible (0-365d)")
if ("dtt_valid" %in% names(dat))
  report_intersection(flag_valid, "dtt_valid flag")
cat("\n")

# of the valid group that LACKS a 5-digit site, where do they go? (OP-sourced,
# or no endoscopy site at all) - tells you what you would lose or need a fallback for
if ("endo_source" %in% names(dat)) {
  cat("Of the sensible-interval group, endoscopy site availability by source\n")
  cat(strrep("-", 68), "\n")
  dat %>%
    filter(sensible) %>%
    mutate(site5 = is_site5(sitetret_endoscopy_apc)) %>%
    count(endo_source, has_5digit_site = site5) %>%
    as.data.frame() %>% print()
}