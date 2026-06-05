# =============================================================================
# Synthetic distance matrix for the provider win/loss (net patient flow) analysis
# -----------------------------------------------------------------------------
# Produces, from the synthetic cohort, the two inputs the Stata win/loss script
# needs, with the SAME column names and code system:
#
#   bowel_pairwise_distance_matrix_SYNTH.dta  -> lsoa11_code, sitecode, total_drive_time
#   NHSHospitals_valid_sites_SYNTH.dta        -> diag_hosp, valid
#
# Design (why this reproduces winners/losers without changing diag_hosp):
#   * Each hospital site gets a 2-D location and an "attractiveness" (pull).
#   * Each LSOA is placed at the location of the LOWEST-pull hospital among the
#     patients who live in it (its "local" district site).
#   * A patient diagnosed at their local site -> core. A co-resident patient
#     diagnosed at a higher-pull site -> arriver there, leaver at the local site.
#   => high-pull sites import, local sites export; total leavers == total arrivers.
#
# diag_hosp and lsoa11_code in the cohort are NOT modified. Tune REALISM KNOBS
# below to push the winner/loser split around.
# =============================================================================

library(tidyverse)
library(haven)

base_dir <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/"
set.seed(20260601)

# ---- REALISM KNOBS ----------------------------------------------------------
SPACE      <- 100      # size of the synthetic 2-D map (arbitrary units)
SCALE_MIN  <- 2.0      # minutes of drive time per map unit
HOME_SD    <- 1.0      # scatter of an LSOA around its local site (keep << site spacing)
NOISE_SD   <- 0.3      # random noise added to each drive time (minutes)
PULL_SPREAD<- 0.8     # >1 sharpens winner/loser contrast, <1 softens it

# ---- Load synthetic cohort (R side, before Stata export) --------------------
cohort <- readRDS(paste0(base_dir, "og_cohort_precwt_SYNTH.rds")) %>%
  select(pseudo_patientid, lsoa11_code, diag_hosp) %>%
  filter(!is.na(diag_hosp), !is.na(lsoa11_code))

# =============================================================================
# 1. Hospital geography: coordinates + attractiveness for each distinct site
# =============================================================================
sites <- tibble(sitecode = sort(unique(cohort$diag_hosp))) %>%
  mutate(
    sx   = runif(n(), 0, SPACE),
    sy   = runif(n(), 0, SPACE),
    # pull in (0,1): rank of a random draw gives an even spread; PULL_SPREAD
    # bends it so a few sites are strong attractors.
    pull = (rank(runif(n())) / (n() + 1))^(1 / PULL_SPREAD)
  )
n_site <- nrow(sites)
cat("Synthetic hospital sites:", n_site, "\n")

# =============================================================================
# 2. Place each LSOA at its lowest-pull ("local district") site + scatter
# =============================================================================
pat <- cohort %>%
  left_join(sites, by = c("diag_hosp" = "sitecode"))

lsoa_home <- pat %>%
  group_by(lsoa11_code) %>%
  slice_min(pull, n = 1, with_ties = FALSE) %>%   # local site = lowest pull here
  ungroup() %>%
  transmute(lsoa11_code,
            hx = sx + rnorm(n(), 0, HOME_SD),
            hy = sy + rnorm(n(), 0, HOME_SD))

cat("Distinct LSOAs:", nrow(lsoa_home),
    "| patients per LSOA:", round(nrow(pat) / nrow(lsoa_home), 2), "\n")

# =============================================================================
# 3. Build the LSOA x site drive-time matrix (full cross)
# =============================================================================
# NOTE: this is (n_LSOA x n_site) rows -- a few million. If memory is tight,
# replace the crossing() with a per-site loop binding rows.
dist_matrix <- tidyr::crossing(
  lsoa_home,
  sites %>% select(sitecode, sx, sy)
) %>%
  mutate(
    total_drive_time = round(
      sqrt((hx - sx)^2 + (hy - sy)^2) * SCALE_MIN + 1 + abs(rnorm(n(), 0, NOISE_SD)),
      2)
  ) %>%
  select(lsoa11_code, sitecode, total_drive_time)

cat("Distance matrix rows:", nrow(dist_matrix), "\n")

# =============================================================================
# 4. Valid-sites lookup (synthetic equivalent of the NHS provider Excel)
# =============================================================================
# All synthetic sites valid by default. To mimic the real exclusions, sample a
# fraction to mark invalid, e.g.:  filter(runif(n()) < 0.85)
valid_sites <- sites %>% transmute(diag_hosp = sitecode, valid = 1L) %>%
  filter(runif(n()) < 0.6)

# =============================================================================
# 5. Write the two Stata inputs
# =============================================================================
write_dta(dist_matrix, paste0(base_dir, "provider_level/OG_pairwise_distance_matrix_SYNTH.dta"))
write_dta(valid_sites, paste0(base_dir, "provider_level/NHSHospitals_valid_sites_SYNTH.dta"))
cat("Wrote distance matrix + valid-sites .dta files.\n")

# =============================================================================
# 6. R-side QC preview: reproduce the win/loss classification quickly
# =============================================================================
# nearest site per LSOA
nearest <- dist_matrix %>%
  group_by(lsoa11_code) %>%
  slice_min(total_drive_time, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(lsoa11_code, nearest_site = sitecode)

flows <- pat %>%
  left_join(nearest, by = "lsoa11_code") %>%
  mutate(core = diag_hosp == nearest_site)

leavers  <- flows %>% filter(!core) %>% count(site = nearest_site, name = "n_leavers")
arrivers <- flows %>% filter(!core) %>% count(site = diag_hosp,    name = "n_arrivers")

site_flow <- sites %>%
  transmute(site = sitecode, pull) %>%
  left_join(leavers,  by = "site") %>%
  left_join(arrivers, by = "site") %>%
  mutate(across(c(n_leavers, n_arrivers), ~replace_na(.x, 0L)),
         n_net_gain = n_arrivers - n_leavers,
         status = case_when(n_net_gain > 0 ~ "net importer",
                            n_net_gain < 0 ~ "net exporter",
                            TRUE           ~ "balanced"))

cat("\n--- QC preview (Poisson significance done in Stata) ---\n")
cat("Total leavers:", sum(site_flow$n_leavers),
    "| total arrivers:", sum(site_flow$n_arrivers),
    "(must be equal)\n")
cat("Patients involved in a flow:",
    round(100 * mean(!flows$core), 1), "% (real ~29%)\n")
print(count(site_flow, status))

# attractiveness should track net gain (high pull -> importer)
cat("\nCorrelation(pull, net_gain):",
    round(cor(site_flow$pull, site_flow$n_net_gain), 2), "\n")

