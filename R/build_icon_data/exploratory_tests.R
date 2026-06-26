# =============================================================================
# Exploratory waiting-times analysis
# -----------------------------------------------------------------------------
# Waiting times to decision-to-treat and from decision-to-treat to treatment,
# broken down by deprivation, year of diagnosis and treatment modality.
#
# Cohort: curative, non-emergency, stage 1-3 OG cancer diagnosed 2015-2022.
# Waiting times require dtt_valid == TRUE (a clean CWT decision-to-treat
# record anchored against a known curative treatment date).
#
# Waits reported:
#   wt_dx_to_dtt   : diagnosis to decision-to-treat
#   wt_dtt_to_tx   : decision-to-treat to first treatment (the 62-day clock)
#   wt_dx_to_tx    : diagnosis to first treatment (sum of the above)
#
# Breakdowns:
#   1. Deprivation (IMD quintile, reversed: 1 = most deprived)
#   2. Year of diagnosis (2015-2022)
#   3. Audit treatment modality (curative categories)
#
# Outputs: printed tables. Nothing saved unless write_outputs <- TRUE.
# =============================================================================

library(dplyr)
library(tidyr)

source("R/build_icon_data/01_define_parameters.R")
og <- readRDS(f_cohort_cwt)

# -----------------------------------------------------------------------------
# cohort restriction: curative, non-emergency, stage 1-3
# -----------------------------------------------------------------------------
curative <- og %>%
  filter(
    received_curative_tx_audit == TRUE,
    route_combined != "Emergency presentation",
    stage_clean %in% c("1", "2", "3")
  ) %>%
  mutate(
    imd  = factor(NHSE_reversed_imd_quintile_lsoas,
                  levels = c("1 - most deprived","2","3","4","5 - least deprived")),
    year = as.integer(ydiag),
    modality = tx_modality_audit,
    wt_dx_to_tx = wt_dx_to_dtt + wt_dtt_to_tx   # full pathway wait
  )

cat("==========  cohort  ==========\n")
cat("curative, non-emergency, stage 1-3 patients:", nrow(curative), "\n")
cat("  of which dtt_valid == TRUE (waits usable):",
    sum(curative$dtt_valid == TRUE, na.rm = TRUE), "\n\n")

# restrict to valid DTT records for the wait tables
wt <- curative %>% filter(dtt_valid == TRUE)

# helper: percentile summary, returns a named data-frame row
pct_row <- function(x) {
  x <- x[!is.na(x)]
  data.frame(
    n      = length(x),
    median = round(median(x), 0),
    p25    = round(quantile(x, 0.25), 0),
    p75    = round(quantile(x, 0.75), 0),
    p90    = round(quantile(x, 0.90), 0)
  )
}

# helper: apply pct_row to a grouped data-frame across the three waits
wait_summary <- function(df, grp_var) {
  df %>%
    group_by(across(all_of(grp_var))) %>%
    summarise(
      n           = n(),
      dx_dtt_med  = round(median(wt_dx_to_dtt,  na.rm = TRUE), 0),
      dx_dtt_p75  = round(quantile(wt_dx_to_dtt, 0.75, na.rm = TRUE), 0),
      dtt_tx_med  = round(median(wt_dtt_to_tx,  na.rm = TRUE), 0),
      dtt_tx_p75  = round(quantile(wt_dtt_to_tx, 0.75, na.rm = TRUE), 0),
      dx_tx_med   = round(median(wt_dx_to_tx,   na.rm = TRUE), 0),
      dx_tx_p75   = round(quantile(wt_dx_to_tx,  0.75, na.rm = TRUE), 0),
      .groups = "drop"
    )
}

cat("columns: n | dx->dtt median | dx->dtt p75 | dtt->tx median |",
    "dtt->tx p75 | dx->tx median | dx->tx p75  (all days)\n\n")

# =============================================================================
# 1. overall
# =============================================================================
cat("==========  overall  ==========\n")
ov <- wt %>%
  summarise(n = n(),
            dx_dtt_med = round(median(wt_dx_to_dtt, na.rm=TRUE),0),
            dx_dtt_p75 = round(quantile(wt_dx_to_dtt,.75,na.rm=TRUE),0),
            dtt_tx_med = round(median(wt_dtt_to_tx, na.rm=TRUE),0),
            dtt_tx_p75 = round(quantile(wt_dtt_to_tx,.75,na.rm=TRUE),0),
            dx_tx_med  = round(median(wt_dx_to_tx,  na.rm=TRUE),0),
            dx_tx_p75  = round(quantile(wt_dx_to_tx, .75,na.rm=TRUE),0))
print(as.data.frame(ov))

# =============================================================================
# 2. by deprivation (IMD quintile)
# =============================================================================
cat("\n==========  by deprivation (IMD quintile, 1 = most deprived)  ==========\n")
print(as.data.frame(wait_summary(wt, "imd")))

# gradient test: is there a monotone deprivation gradient?
imd_med <- wt %>%
  filter(!is.na(imd)) %>%
  group_by(imd) %>%
  summarise(dtt_tx_med = median(wt_dtt_to_tx, na.rm=TRUE), .groups="drop") %>%
  arrange(imd)
cat("\nDTT-to-tx medians across IMD quintiles (look for gradient):\n")
print(as.data.frame(imd_med))

# =============================================================================
# 3. by year of diagnosis
# =============================================================================
cat("\n==========  by year of diagnosis  ==========\n")
print(as.data.frame(wait_summary(wt, "year")))

# flag the covid dip years
cat("\n(note: 2020-2021 may reflect covid disruption to pathways)\n")

# =============================================================================
# 4. by audit treatment modality (curative categories only)
# =============================================================================
cat("\n==========  by treatment modality  ==========\n")
# keep only the curative modality categories for this table
curative_mods <- c("Surgery only","Surgery plus SACT/RT","Definitive chemoRT",
                   "EMR/ESD","Curative RT only")
wt_mod <- wt %>% filter(modality %in% curative_mods) %>%
  mutate(modality = factor(modality, levels = curative_mods))
print(as.data.frame(wait_summary(wt_mod, "modality")))

# =============================================================================
# 5. deprivation by modality - is the gradient different by treatment type?
# =============================================================================
cat("\n==========  DTT-to-tx median by modality and IMD  ==========\n")
wt_mod %>%
  filter(!is.na(imd)) %>%
  group_by(modality, imd) %>%
  summarise(n = n(), dtt_tx_med = round(median(wt_dtt_to_tx, na.rm=TRUE), 0),
            .groups = "drop") %>%
  pivot_wider(names_from = imd, values_from = c(n, dtt_tx_med),
              names_glue = "{imd}_{.value}") %>%
  as.data.frame() %>% print()

# =============================================================================
# 6. deprivation by year - does the gradient change over time?
# =============================================================================
cat("\n==========  DTT-to-tx median by year and IMD  ==========\n")
wt %>%
  filter(!is.na(imd)) %>%
  group_by(year, imd) %>%
  summarise(n = n(), dtt_tx_med = round(median(wt_dtt_to_tx, na.rm=TRUE), 0),
            .groups = "drop") %>%
  pivot_wider(names_from = imd, values_from = c(n, dtt_tx_med),
              names_glue = "{imd}_{.value}") %>%
  arrange(year) %>% as.data.frame() %>% print()

cat("\n==========  notes  ==========\n")
cat("Waits in days. dtt_valid == TRUE filter applied throughout.\n")
cat("Negative values arise where the CWT treatment date precedes the DTT date\n")
cat("(a known CWT data-quality artefact, flagged dtt_valid = FALSE; these are\n")
cat("excluded here).\n")
cat("IMD direction: 1 = most deprived, 5 = least deprived (reversed from NDRS default).\n")

# optional: save tables for further work
if (!exists("write_outputs")) write_outputs <- FALSE
if (isTRUE(write_outputs)) {
  out <- list(
    overall         = as.data.frame(ov),
    by_imd          = as.data.frame(wait_summary(wt, "imd")),
    by_year         = as.data.frame(wait_summary(wt, "year")),
    by_modality     = as.data.frame(wait_summary(wt_mod, "modality"))
  )
  saveRDS(out, file.path(dir_icon, "explore_waiting_times.rds"))
  cat("\ntables saved ->", file.path(dir_icon, "explore_waiting_times.rds"), "\n")
}
