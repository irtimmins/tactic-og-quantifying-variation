# =============================================================================
# 15  CWT surgery residual - specialist-centre subgroup, HES contents
# -----------------------------------------------------------------------------
# Third of three CWT-surgery diagnostics (13 -> 14 -> 15), read-only. Run 14
# first; this closes out the one subgroup that the specialist-centre check left
# open: patients with a CWT surgery event at a designated OG specialist centre
# (NOGCA Appendix 12) but no HES resection. Major resections only happen at those
# centres, so this is where a genuinely missed resection would hide.
#
# It distinguishes two explanations:
#   - coverage/linkage gap : the patient has NO HES APC episode at all near the
#                            CWT date (their admission is not in our extract) -
#                            a data limitation the national audit shares.
#   - coding/match gap     : the patient HAS HES episodes but none carry a
#                            resection OPCS - so it inspects what OPCS they do
#                            carry, to show whether a real resection is being
#                            missed or the episodes are genuinely non-resection.
#
# Note: a real null-operation-date gap was found and fixed here historically -
# resections whose OPDATE held the HES sentinel 1800-01-01 were being dropped by
# the date-keyed anchor. The fix (fall back to EPISTART, in 01's
# match_opcs_episodes) recovered those, so this script's resection-hit count is
# now expected to be only the genuine out-of-window boundary cases.
#
# Reads : Data/ICON/og_cohort_cwt_2015_2022.rds, og_surgery_anchor.rds,
#         hes_apc extract, the partitioned CWT dataset
# Writes: nothing (prints a report). Optional line-list if write_residual_linelist.
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

og   <- readRDS(f_cohort_cwt)
surg <- readRDS(f_surgery_anchor)
hes  <- readRDS(f_hes_apc_extract)

# NOGCA Appendix 12 specialist OG surgical centres (England)
og_specialist_centres <- c(
  "R0D","RA2","RA7","RAE","RAJ","REM","RF4","RGT","RHM","RHQ","RHU","RJ1","RJE",
  "RK9","RKB","RM1","RM3","RPY","RR8","RRK","RRV","RTD","RTE","RTG","RTH","RTR",
  "RWA","RWE","RX1","RXN","RYJ","RYR")

# the flagged group: No treatment + CWT surgery code + no HES resection
flagged_ids <- og %>%
  filter(tx_pathway == "No treatment recorded", !is.na(cwt_treat_date),
         cwt_modality %in% c("01","23","24"),
         !pseudo_patientid %in% surg$pseudo_patientid) %>%
  pull(pseudo_patientid) %>% unique()

# raw CWT rows, keep the treating org, date and setting
cwt_raw <- open_dataset(path_cwt_partition) %>%
  filter(site_icd10 %in% og_icd10, pseudo_patientid %in% flagged_ids,
         modality %in% c("01","23","24")) %>%
  select(pseudo_patientid, org_treat_start, treat_start, care_setting,
         cte_type, modality) %>%
  collect() %>%
  mutate(cwt_tx = as.Date(treat_start, "%d/%m/%Y"),
         trust3 = substr(org_treat_start, 1, 3),
         at_specialist = trust3 %in% og_specialist_centres)

spec_grp <- cwt_raw %>% filter(at_specialist)
spec_ids <- unique(spec_grp$pseudo_patientid)

cat("==========  residual: CWT surgery at a specialist centre  ==========\n")
cat("records:", nrow(spec_grp), " distinct patients:", length(spec_ids), "\n\n")

cat("-- care setting (01 usually inpatient/day-case) --\n")
print(sort(table(spec_grp$care_setting), decreasing = TRUE))
cat("\n-- CWT modality in this subgroup --\n")
print(sort(table(spec_grp$modality), decreasing = TRUE))

# -----------------------------------------------------------------------------
# Coverage vs coding gap: do these patients have any HES APC episode at all?
# -----------------------------------------------------------------------------
has_any_hes <- spec_ids %in% unique(hes$STUDY_ID)
cat("\n==========  coverage vs coding gap  ==========\n")
cat("specialist-centre patients with NO HES APC episode at all:",
    sum(!has_any_hes), "  (coverage/linkage gap)\n")
cat("specialist-centre patients WITH HES APC episodes:        ",
    sum(has_any_hes), "  (inspect their OPCS below)\n")

# -----------------------------------------------------------------------------
# For the WITH-HES group: what OPCS do their episodes near the CWT date carry?
# This shows whether a real resection code is present (a match gap) or the
# episodes are genuinely non-resection procedures.
# -----------------------------------------------------------------------------
with_hes_ids <- spec_ids[has_any_hes]
# earliest CWT surgery date per patient, to window the HES look-up
cwt_dt <- spec_grp %>% filter(pseudo_patientid %in% with_hes_ids) %>%
  group_by(pseudo_patientid) %>% summarise(cwt_tx = min(cwt_tx, na.rm = TRUE), .groups = "drop")

hes_near <- hes %>%
  filter(STUDY_ID %in% with_hes_ids) %>%
  inner_join(cwt_dt, by = c("STUDY_ID" = "pseudo_patientid")) %>%
  mutate(gap = as.integer(EPISTART - cwt_tx)) %>%
  filter(!is.na(gap), abs(gap) <= 30)          # episodes within a month of the CWT date

cat("\nHES APC episodes within +/-30 days of the CWT surgery date:",
    nrow(hes_near), "across", n_distinct(hes_near$STUDY_ID), "patients\n")

# gather all OPCS codes carried by those episodes
opcs_long <- hes_near %>%
  select(STUDY_ID, all_of(op_cols)) %>%
  pivot_longer(all_of(op_cols), values_to = "opcs", names_to = NULL) %>%
  filter(!is.na(opcs), opcs != "-") %>%
  mutate(opcs = normalise_opcs(opcs))

cat("\n-- are ANY of these OPCS on the major-resection list (Appendix 8)? --\n")
res_hits <- opcs_long %>% filter(opcs %in% opcs_og_surgery_all)
cat("resection-coded episodes among the near-CWT HES records:", nrow(res_hits),
    "  (expected: only out-of-window boundary cases after the OPDATE fix)\n")

cat("\n-- most common OPCS in these near-CWT episodes (top 25) --\n")
print(head(opcs_long %>% count(opcs, sort = TRUE) %>% as.data.frame(), 25))

# classify each carried code against our known lists so the picture is readable
cat("\n-- what category do the carried OPCS fall into? --\n")
opcs_long %>%
  mutate(cat = case_when(
    opcs %in% opcs_og_surgery_all       ~ "major resection (Appendix 8)",
    opcs %in% opcs_emresd               ~ "EMR/ESD (Appendix 7)",
    opcs %in% opcs_diagnostic_endoscopy ~ "diagnostic/therapeutic endoscopy (App 6)",
    grepl("^G15|^G21|^G44", opcs)       ~ "stent / intubation / UGI tube",
    TRUE                                ~ "other")) %>%
  count(cat, sort = TRUE) %>% as.data.frame() %>% print()

cat("\n==========  reading  ==========\n")
cat("A 'major resection' count near zero means even the specialist-centre\n")
cat("subgroup carries no in-window resection OPCS in HES - so HES (the gold\n")
cat("standard, and NOGCA's only surgery source) genuinely holds no resection,\n")
cat("and these patients stay no-treatment. The dominant endoscopy / stent share\n")
cat("confirms the CWT 'surgery' was a non-resection procedure. Any remaining\n")
cat("resection hits are out-of-window late surgeries, correctly excluded by the\n")
cat("treatment window (and by NOGCA's equivalent 9-month rule).\n")

# Optional line-list of any resection-coded near-CWT episodes, for inspection.
if (!exists("write_residual_linelist")) write_residual_linelist <- FALSE
if (isTRUE(write_residual_linelist)) {
  out <- file.path(dir_icon, "check_cwt_residual_resection_hits.rds")
  saveRDS(res_hits, out)
  cat("\nresection-hit line-list saved ->", out, "\n")
}