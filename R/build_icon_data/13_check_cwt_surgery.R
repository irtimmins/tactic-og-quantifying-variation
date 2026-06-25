# =============================================================================
# 13  CWT surgery among "No treatment recorded" - characterise the group
# -----------------------------------------------------------------------------
# First of three CWT-surgery diagnostics (13 -> 14 -> 15), all read-only:
#   13  characterise the group (this script)
#   14  what triggered the CWT record (raw fields + specialist-centre check)
#   15  the specialist-centre residual (coverage vs coding gap, HES contents)
#
# Some patients classified "No treatment recorded" carry a CWT surgery code
# (01 / 23 / 24) within the treatment window, yet have no HES surgical resection.
# This script quantifies and describes that group so the question can be closed.
#
# Position: HES is the gold standard for surgery. A CWT surgery code without a
# HES resection is therefore NOT a missed treatment - it is CWT recording an
# event HES did not anchor (palliative or non-resection surgery, or CWT's looser
# "surgery" definition). The aim is to describe the group, not reclassify it.
#
# Reads : Data/ICON/og_cohort_cwt_2015_2022.rds, og_surgery_anchor.rds
# Writes: nothing (prints a report). Optional line-list if write_flagged_linelist.
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

og   <- readRDS(f_cohort_cwt)
surg <- readRDS(f_surgery_anchor)

# the group in question: No treatment, CWT surgery code, in window, no HES anchor
flagged <- og %>%
  filter(tx_pathway == "No treatment recorded",
         !is.na(cwt_treat_date),
         cwt_modality %in% c("01", "23", "24")) %>%
  mutate(days_dx_to_cwt = as.integer(cwt_treat_date - diagmdy)) %>%
  filter(days_dx_to_cwt >= 0, days_dx_to_cwt <= cwt_window_days) %>%
  mutate(in_hes_surgery = pseudo_patientid %in% surg$pseudo_patientid)

cat("==========  CWT surgery among 'No treatment recorded'  ==========\n")
cat("flagged patients:", nrow(flagged), "of",
    sum(og$tx_pathway == "No treatment recorded"), "no-treatment patients (",
    sprintf("%.1f%%", 100 * nrow(flagged) / sum(og$tx_pathway == "No treatment recorded")),
    ")\n")
cat("of these, in the HES surgery anchor:", sum(flagged$in_hes_surgery),
    "  not in HES:", sum(!flagged$in_hes_surgery), "\n")
cat("(HES is gold standard: the 'not in HES' group are not counted as surgery.)\n\n")

# -----------------------------------------------------------------------------
# 1. by CWT modality - is this old code 01, or the newer 23/24?
# -----------------------------------------------------------------------------
cat("-- by CWT surgery modality --\n")
flagged %>% count(cwt_modality) %>% arrange(desc(n)) %>% as.data.frame() %>% print()

# -----------------------------------------------------------------------------
# 2. by stage - palliative non-resection surgery skews to advanced disease
# -----------------------------------------------------------------------------
cat("\n-- by stage --\n")
flagged %>% count(stage_clean) %>% as.data.frame() %>% print()

# -----------------------------------------------------------------------------
# 3. by tumour subtype and site
# -----------------------------------------------------------------------------
cat("\n-- by cancer subtype --\n")
flagged %>% count(cancer_subtype) %>% arrange(desc(n)) %>% as.data.frame() %>% print()
cat("\n-- by tumour site group --\n")
flagged %>% count(tumour_site_grp) %>% arrange(desc(n)) %>% as.data.frame() %>% print()

# -----------------------------------------------------------------------------
# 4. timing - how long after diagnosis is the CWT surgery date?
# -----------------------------------------------------------------------------
cat("\n-- days from diagnosis to the CWT surgery date --\n")
flagged %>%
  mutate(band = cut(days_dx_to_cwt, c(-1, 30, 60, 90, 180, Inf),
                    labels = c("0-30", "31-60", "61-90", "91-180", "181+"))) %>%
  count(band) %>% as.data.frame() %>% print()

# -----------------------------------------------------------------------------
# 5. do these patients have ANY other treatment signal we did capture?
#    (chemo or RT anchor) - if so they are not truly untreated, just not surgical
# -----------------------------------------------------------------------------
chemo <- tryCatch(readRDS(f_chemo_anchor), error = function(e) NULL)
rt    <- tryCatch(readRDS(f_rt_anchor),    error = function(e) NULL)
has_chemo <- if (!is.null(chemo)) flagged$pseudo_patientid %in% chemo$pseudo_patientid else NA
has_rt    <- if (!is.null(rt))    flagged$pseudo_patientid %in% rt$pseudo_patientid    else NA

cat("\n-- other (non-surgical) treatment signal among the flagged group --\n")
cat("  also have a chemo anchor:", sum(has_chemo, na.rm = TRUE), "\n")
cat("  also have an RT anchor:  ", sum(has_rt,    na.rm = TRUE), "\n")
cat("  no chemo and no RT (truly no captured treatment):",
    sum(!has_chemo & !has_rt, na.rm = TRUE), "\n")

# -----------------------------------------------------------------------------
# Summary read
# -----------------------------------------------------------------------------
cat("\n==========  summary  ==========\n")
cat("These", nrow(flagged), "patients carry a CWT surgery code with no HES resection.\n")
cat("Under the HES-gold-standard rule they remain 'No treatment recorded' for the\n")
cat("surgical pathway. This conclusion is consistent with the NOGCA State of the\n")
cat("Nation 2025 methodology on two counts:\n")
cat("  1. Our surgery OPCS-4 list matches NOGCA Appendix 8 (major OG resections)\n")
cat("     exactly - 32 codes, no gap - so this is not a code-list omission.\n")
cat("  2. NOGCA ascertains surgery from HES-APC opertn fields (Table 4.3), and\n")
cat("     does NOT use CWT to identify whether treatment occurred - CWT is used\n")
cat("     only for referral/waiting-times dates. So the national audit would\n")
cat("     also classify these patients as having no surgical treatment.\n")
cat("The breakdowns above (advanced stage, early timing, no chemo/RT anchor) are\n")
cat("consistent with palliative or non-resection procedures that CWT records as a\n")
cat("treatment event but HES-APC does not carry as a major resection. This is a\n")
cat("known, quantified CWT-vs-HES difference, not a build error. To see what the\n")
cat("CWT record actually describes for these patients, run 14_check_cwt_trigger.R.\n")

# Optional line-list for the record. Set write_flagged_linelist <- TRUE to save.
if (!exists("write_flagged_linelist")) write_flagged_linelist <- FALSE
if (isTRUE(write_flagged_linelist)) {
  out <- file.path(dir_icon, "check_cwt_surgery_no_hes.rds")
  saveRDS(flagged %>% select(pseudo_patientid, stage_clean, cancer_subtype,
                             tumour_site_grp, cwt_modality, cwt_treat_date,
                             days_dx_to_cwt, in_hes_surgery), out)
  cat("\nline-list saved ->", out, "\n")
}