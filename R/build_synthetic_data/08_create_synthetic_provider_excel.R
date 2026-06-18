# =============================================================================
# Synthetic provider-characteristics Excel for the site-level waiting-times merge
# -----------------------------------------------------------------------------
# Builds a small Excel that mimics NHSHospitals_services_*.xlsx but contains ONLY
# the columns the Stata provider script actually reads, one row per site code in
# NHSHospitals_valid_sites_SYNTH.dta.
#
# Columns the Stata script uses (and nothing else is needed):
#   Trust_Name, Trust_Name_colour, Hospital_site_code, Bowel_ca_surgery,
#   Comprehensive_centre, Teaching_hospitals, Latest_Rating,
#   Staff_engagement, Moral, mean (= bed occupancy rate)
#
# Trust-level fields (Trust_Name, colour, Staff_engagement, Moral, mean) are kept
# constant within a trust (= first 3 chars of the site code), as in the real file.
# Site-level fields (surgery flags, centre type, CQC rating) vary by site.
# =============================================================================

library(tidyverse)
library(haven)
library(writexl)   # install.packages("writexl") if needed; no Java required
#install.packages("writexl")

library(here)
syn_dir  <- here("Data", "synthetic")
prov_dir <- file.path(syn_dir, "provider_level")
dir.create(prov_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260601)

# ---- Site codes from the synthetic valid-sites file -------------------------
sites <- read_dta(file.path(prov_dir, "NHSHospitals_valid_sites_SYNTH.dta")) %>%
  distinct(diag_hosp) %>%
  transmute(Hospital_site_code = diag_hosp,
            trust = substr(diag_hosp, 1, 3))     # NHS nesting: site -> trust prefix

# ---- Trust-level attributes (constant within a trust) -----------------------
trusts <- sites %>%
  distinct(trust) %>%
  mutate(
    Trust_Name = sprintf("Synthetic NHS Trust %03d", row_number()),
    # ~12% of trusts get an excluded colour, to exercise the Stata colour drop.
    # Set prob of the last three to 0 if you want NO colour exclusions.
    Trust_Name_colour = sample(
      c("Blank", "Yellow", "Green", "Grey", "Light Red", "Pink Red", "Orange"),
      n(), replace = TRUE, prob = c(.45, .20, .15, .08, .05, .04, .03)),
    Staff_engagement = round(rnorm(n(), 6.90, 0.12), 4),   # ~1-10 scale
    Moral            = round(rnorm(n(), 5.95, 0.15), 4),
    mean             = round(pmin(pmax(rnorm(n(), 0.93, 0.025), 0.82), 0.99), 4) # bed occ
  )

# ---- Assemble: site-level fields + trust attributes -------------------------
provider <- sites %>%
  left_join(trusts, by = "trust") %>%
  mutate(
    Bowel_ca_surgery     = sample(c(1, 0, NA), n(), replace = TRUE, prob = c(.50, .40, .10)),
    Oesophageal_Surgery  = rbinom(n(), 1, 0.20),   # extra (OG-relevant); script can use this
    Comprehensive_centre = rbinom(n(), 1, 0.15),
    Teaching_hospitals   = rbinom(n(), 1, 0.25),
    Latest_Rating = sample(
      c("Outstanding", "Good", "Requires Improvement", "Inadequate", "Not rated"),
      n(), replace = TRUE, prob = c(.05, .45, .35, .08, .07)),
    # a few missing bed-occupancy values, as in the real file
    mean = if_else(runif(n()) < 0.05, NA_real_, mean)
  ) %>%
  select(Trust_Name, Trust_Name_colour, Hospital_site_code,
         Bowel_ca_surgery, Oesophageal_Surgery,
         Comprehensive_centre, Teaching_hospitals, Latest_Rating,
         Staff_engagement, Moral, mean)

# ---- Write Excel (default sheet name is "Sheet1", which the Stata script reads)
out_xlsx <- file.path(prov_dir, "NHSHospitals_services_SYNTH.xlsx")
write_xlsx(provider, out_xlsx)

cat("Wrote", nrow(provider), "site rows to:\n  ", out_xlsx, "\n")
cat("Distinct trusts:", n_distinct(provider$Trust_Name), "\n")
cat("Colour-excluded sites:",
    sum(provider$Trust_Name_colour %in% c("Light Red", "Pink Red", "Orange")), "\n")
print(count(provider, Latest_Rating))
