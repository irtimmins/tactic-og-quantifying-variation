# =============================================================================
# Oesophageal (C15) waiting-times analysis - curative, non-emergency, stage 1-3
# -----------------------------------------------------------------------------
# Restricted to oesophageal cancer (cancer == "oesophagus", ICD-10 C15) per the
# TACTIC hospital-QI inclusion criterion. Curative, non-emergency, stage 1-3,
# diagnosed 2015-2022. Waiting times use dtt_valid == TRUE.
#
# Wait components covered (all in days), the endoscopy- and decision-to-treat
# pathway end to end:
#   endo_to_dtt : diagnostic endoscopy -> decision to treat   (audit clock start)
#   dx_to_dtt   : diagnosis            -> decision to treat
#   dtt_to_tx   : decision to treat    -> first treatment
#   endo_to_tx  : diagnostic endoscopy -> first treatment      (endo_to_dtt + dtt_to_tx)
#   dx_to_tx    : diagnosis            -> first treatment      (dx_to_dtt   + dtt_to_tx)
#
# All breakdowns are LONG: one row per group x wait component, with the grouping
# (IMD, year, modality) reading DOWN the page, never spread across columns.
#
# Outputs: printed tables. Nothing saved unless write_outputs <- TRUE.
# Needs: dplyr, tidyr.
# =============================================================================

library(dplyr)
library(tidyr)

source("R/build_icon_data/01_define_parameters.R")
og <- readRDS(f_cohort_cwt)

# -----------------------------------------------------------------------------
# cohort: oesophageal (C15), curative, non-emergency, stage 1-3
# -----------------------------------------------------------------------------
coh <- og %>%
  filter(received_curative_tx_audit == TRUE,
         route_combined != "Emergency presentation",
         cancer == "oesophagus",
         stage_clean %in% c("1","2","3")) %>%
  mutate(
    imd      = factor(NHSE_reversed_imd_quintile_lsoas,
                      levels = c("1 - most deprived","2","3","4","5 - least deprived")),
    year     = as.integer(ydiag),
    modality = factor(tx_modality_audit),
    wt_endo_to_tx = wt_endo_to_dtt + wt_dtt_to_tx,
    wt_dx_to_tx   = wt_dx_to_dtt   + wt_dtt_to_tx
  )

cat("==========  cohort  ==========\n")
cat("oesophageal (C15), curative, non-emergency, stage 1-3:", nrow(coh), "\n")
cat("  of which dtt_valid == TRUE (waits usable):",
    sum(coh$dtt_valid == TRUE, na.rm = TRUE), "\n")

wt <- coh %>% filter(dtt_valid == TRUE)

wait_cols <- c(endo_to_dtt = "wt_endo_to_dtt",
               dx_to_dtt   = "wt_dx_to_dtt",
               dtt_to_tx   = "wt_dtt_to_tx",
               endo_to_tx  = "wt_endo_to_tx",
               dx_to_tx    = "wt_dx_to_tx")
wait_levels <- names(wait_cols)

# long summary of every wait component within a grouping. grp = NULL -> overall.
# Output is blocked BY WAIT COMPONENT: all rows for one wait (across the grouping)
# sit together, so a gradient reads straight down the block.
wait_long <- function(df, grp = NULL) {
  d <- df %>% rename(all_of(wait_cols))
  d <- if (is.null(grp)) mutate(d, .grp = "all") else rename(d, .grp = all_of(grp))
  long <- d %>%
    select(.grp, all_of(wait_levels)) %>%
    pivot_longer(all_of(wait_levels), names_to = "wait", values_to = "days") %>%
    filter(!is.na(days)) %>%
    group_by(wait, .grp) %>%
    summarise(n = n(),
              median = round(median(days), 0),
              p25 = round(quantile(days, .25), 0),
              p75 = round(quantile(days, .75), 0),
              p90 = round(quantile(days, .90), 0),
              .groups = "drop") %>%
    mutate(wait = factor(wait, levels = wait_levels)) %>%
    arrange(wait, .grp) %>%
    select(wait, .grp, n, median, p25, p75, p90)
  names(long)[2] <- if (is.null(grp)) "cohort" else grp
  as.data.frame(long)
}

# 1. overall
cat("\n==========  1. overall, by wait component  ==========\n")
print(wait_long(wt))

# 2. by deprivation
cat("\n==========  2. by deprivation (IMD quintile, 1 = most deprived)  ==========\n")
print(wait_long(wt %>% filter(!is.na(imd)), "imd"))

# 3. by year
cat("\n==========  3. by year of diagnosis  ==========\n")
print(wait_long(wt, "year"))
cat("(2020-2021 may reflect covid disruption to pathways)\n")

# 4. by modality
cat("\n==========  4. by treatment modality  ==========\n")
curative_mods <- c("Surgery only","Surgery plus SACT/RT","Definitive chemoRT",
                   "EMR/ESD","Curative RT only")
wt_mod <- wt %>% filter(modality %in% curative_mods) %>%
  mutate(modality = factor(modality, levels = curative_mods))
print(wait_long(wt_mod, "modality"))

# 5. key-leg medians by deprivation (long)
cat("\n==========  5. key-leg medians by deprivation (long)  ==========\n")
wt %>%
  filter(!is.na(imd)) %>%
  select(imd, endo_to_dtt = wt_endo_to_dtt, dtt_to_tx = wt_dtt_to_tx) %>%
  pivot_longer(c(endo_to_dtt, dtt_to_tx), names_to = "leg", values_to = "days") %>%
  filter(!is.na(days)) %>%
  group_by(leg, imd) %>%
  summarise(n = n(), median = round(median(days), 0),
            p75 = round(quantile(days, .75), 0), .groups = "drop") %>%
  mutate(leg = factor(leg, levels = c("endo_to_dtt","dtt_to_tx"))) %>%
  arrange(leg, imd) %>% as.data.frame() %>% print()

# 6. decision-to-treatment wait by year and deprivation (long)
cat("\n==========  6. decision-to-treatment wait by year and deprivation (long)  ==========\n")
wt %>%
  filter(!is.na(imd)) %>%
  group_by(year, imd) %>%
  summarise(n = n(), dtt_to_tx_med = round(median(wt_dtt_to_tx, na.rm = TRUE), 0),
            .groups = "drop") %>%
  arrange(year, imd) %>% as.data.frame() %>% print()

cat("\n==========  notes  ==========\n")
cat("Oesophageal (C15) only. Waits in days; dtt_valid == TRUE throughout.\n")
cat("endo_to_dtt is the audit-aligned clock start (preferred over dx_to_dtt).\n")
cat("Negative waits (CWT treatment date before the decision-to-treat date) are a\n")
cat("known CWT data-quality artefact and are already excluded via dtt_valid.\n")
cat("IMD direction: 1 = most deprived, 5 = least deprived (reversed from NDRS).\n")

if (!exists("write_outputs")) write_outputs <- FALSE
if (isTRUE(write_outputs)) {
  out <- list(overall    = wait_long(wt),
              by_imd      = wait_long(wt %>% filter(!is.na(imd)), "imd"),
              by_year     = wait_long(wt, "year"),
              by_modality = wait_long(wt_mod, "modality"))
  saveRDS(out, file.path(dir_icon, "explore_waiting_times_c15.rds"))
  cat("\ntables saved ->", file.path(dir_icon, "explore_waiting_times_c15.rds"), "\n")
}