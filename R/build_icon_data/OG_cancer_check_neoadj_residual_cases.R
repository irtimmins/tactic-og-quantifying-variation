# =============================================================================
# OG cancer - neoadjuvant residual check (the ~360 surgery-anchored cases)
# -----------------------------------------------------------------------------
# The neoadjuvant check showed 95-98% of neoadjuvant patients anchor on the
# earlier chemo/RT event, exactly per guidance. The small residual (~4%) anchor
# on the later surgery instead. This asks why: did those patients have a chemo
# (or chemoRT/RT) CWT row that the merge passed over because it fell outside the
# DTT window or lost the earliest-consistent-row selection, or do they simply
# have no systemic CWT row, so surgery was the only thing to anchor on?
#
# This is a diminishing-returns check (~360 of ~9,000 neoadjuvant patients, and
# it does not change the treatment classification), run once to characterise the
# residual and close the investigation.
# =============================================================================

library(tidyverse)
library(arrow)

base_dir <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/"
cwt_path <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/CWT/11_CWT_data_partitioned"

tx_window_days <- 270L
dtt_min_offset <- -30L
og_icd10 <- c("C150","C151","C152","C153","C154","C155","C158","C159","C15",
              "C160","C161","C162","C163","C164","C165","C166","C168","C169","C16")

og <- readRDS(paste0(base_dir, "og_cohort_cwt_2015_2022.rds"))

neoadj_pathways <- c("Surgery + neoadjuvant chemo",
                     "Surgery + neoadjuvant chemoRT",
                     "Surgery + neoadjuvant RT")

# the residual: neoadjuvant pathway but CWT anchored on a surgery modality
resid <- og %>%
  filter(tx_pathway %in% neoadj_pathways,
         !is.na(cwt_treat_date),
         cwt_modality %in% c("01","23","24")) %>%
  select(pseudo_patientid, diagmdy, tx_pathway, cwt_modality,
         cwt_treat_date, sact_date, rt_date, surgery_date)

cat("Surgery-anchored neoadjuvant patients (the residual):", nrow(resid), "\n")
cat("by pathway:\n")
resid %>% count(tx_pathway) %>% print()

# pull ALL their CWT rows (not just the anchored one) to see what was available
resid_ids <- unique(resid$pseudo_patientid)
cwt_all <- open_dataset(cwt_path) %>%
  filter(site_icd10 %in% og_icd10) %>%
  collect() %>%
  mutate(pseudo_patientid = as.character(pseudo_patientid),
         cwt_treat = as.Date(treat_start,        format = "%d/%m/%Y"),
         cwt_dtt   = as.Date(treat_period_start, format = "%d/%m/%Y")) %>%
  filter(pseudo_patientid %in% resid_ids)

# classify each CWT row's modality into a broad group
grp_of <- function(m) case_when(
  m %in% c("01","23","24") ~ "surgery",
  m %in% c("02","14","15") ~ "chemo",
  m == "04"                ~ "chemoRT",
  m %in% c("05","06","13") ~ "radiotherapy",
  m %in% c("07","08","09") ~ "palliative",
  TRUE                     ~ "other")

cwt_all <- cwt_all %>%
  left_join(resid %>% select(pseudo_patientid, diagmdy), by = "pseudo_patientid") %>%
  mutate(grp        = grp_of(modality),
         dtt_offset = as.integer(cwt_dtt - diagmdy),
         in_window  = !is.na(dtt_offset) &
           dtt_offset >= dtt_min_offset & dtt_offset <= tx_window_days)

# per patient: did a systemic (chemo/chemoRT/RT) CWT row exist at all, and was
# it in the merge's DTT window?
per <- cwt_all %>%
  group_by(pseudo_patientid) %>%
  summarise(
    has_systemic_row        = any(grp %in% c("chemo","chemoRT","radiotherapy")),
    has_systemic_in_window  = any(grp %in% c("chemo","chemoRT","radiotherapy") & in_window),
    n_cwt_rows              = n(),
    .groups = "drop"
  )

cat("\nOf the residual, what systemic CWT rows did they have?\n")
per %>% summarise(
  n                       = n(),
  with_systemic_row       = sum(has_systemic_row),
  with_systemic_in_window = sum(has_systemic_in_window),
  no_systemic_row         = sum(!has_systemic_row)
) %>% print(width = Inf)

cat("\nAttribution of the residual:\n")
per %>%
  mutate(reason = case_when(
    has_systemic_in_window  ~ "1 systemic CWT row WAS in window (merge selection)",
    has_systemic_row        ~ "2 systemic CWT row existed but OUT of window",
    TRUE                    ~ "3 no systemic CWT row at all (surgery only in CWT)"
  )) %>%
  count(reason) %>% mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(reason) %>% print(width = Inf)

cat("\nReading:\n",
    " reason 1 -> a chemo/RT row was available and in-window but the merge took\n",
    "   the surgery row instead; a merge tie-break worth a look IF the count is\n",
    "   material (it governs only the CWT date, not the pathway/curative flag).\n",
    " reason 2 -> the systemic row sat outside -30..", tx_window_days, "d, so the merge\n",
    "   correctly could not use it; the surgery row was the only in-window option.\n",
    " reason 3 -> CWT genuinely holds only a surgery row for these patients; the\n",
    "   neoadjuvant chemo was not recorded in CWT. Nothing to fix.\n",
    "\nEither way this is ~", nrow(resid), " patients and does not change tx_pathway,\n",
    " first_tx_date, or the curative/any-treatment flags.\n")