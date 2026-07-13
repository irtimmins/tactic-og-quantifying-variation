# =============================================================================
# 04b  APC vs OP split of the endoscopy anchor
# -----------------------------------------------------------------------------
# The endoscopy anchor built in 04 takes APC first and fills in OP only for the
# patients APC missed. This matters for site resolution: APC carries a full
# 5-character SITETRET, OP only a 3-character provider code, so the site variable
# ends up mixed resolution. Before deciding whether that mix is workable, this
# just reports the split plainly - how many patients land in each source, and
# what share of the whole cohort has no endoscopy anchor at all.
#
# Reads : Data/ICON/og_endoscopy_anchor.rds, Data/ICON/ncras_og_2015_2022.rds
# Writes: nothing - prints only
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

library(dplyr)

endoscopy_anchor <- readRDS(f_endoscopy_anchor)
ncras_og         <- readRDS(f_ncras_cohort)

n_cohort <- n_distinct(ncras_og$pseudo_patientid)
n_anchor <- n_distinct(endoscopy_anchor$pseudo_patientid)

by_source <- endoscopy_anchor %>%
  count(endo_source, name = "n_patients") %>%
  mutate(pct_of_anchor = round(100 * n_patients / n_anchor, 1),
         pct_of_cohort = round(100 * n_patients / n_cohort, 1))

cat("Whole cohort              :", n_cohort, "patients\n")
cat("Have an endoscopy anchor  :", n_anchor, "patients (",
    round(100 * n_anchor / n_cohort, 1), "% of the cohort)\n")
cat("No endoscopy anchor found :", n_cohort - n_anchor, "patients (",
    round(100 * (n_cohort - n_anchor) / n_cohort, 1), "% of the cohort)\n\n")

cat("Split of the anchor by source\n")
cat(strrep("-", 60), "\n")
print(as.data.frame(by_source))
cat("\npct_of_anchor is share of patients WITH an anchor; pct_of_cohort is\n")
cat("share of the whole cohort, including those with no anchor at all.\n\n")

# the number that actually matters for site resolution: what share of the
# WHOLE cohort would get a 5-character endoscopy site (APC) versus a
# 3-character trust-level one (OP) versus none at all.
site_resolution <- tibble(
  resolution = c("5-character site (APC)", "3-character trust only (OP)", "no endoscopy anchor"),
  n_patients = c(sum(endoscopy_anchor$endo_source == "APC"),
                 sum(endoscopy_anchor$endo_source == "OP"),
                 n_cohort - n_anchor)) %>%
  mutate(pct_of_cohort = round(100 * n_patients / n_cohort, 1))

cat("Site resolution this gives you, as a share of the WHOLE cohort\n")
cat(strrep("-", 60), "\n")
print(as.data.frame(site_resolution))

cat("\nIf the OP row is small, most of the cohort gets a full site code and the\n")
cat("trust-level tail is a minor caveat. If it is large, a substantial share of\n")
cat("any site-level analysis would be resting on trust-level-only records.\n")