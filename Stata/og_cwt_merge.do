*  OG cancer - CWT merge in Stata
*  ---------------------------------------------------------------------------
*  Reproduces the R minimal merge (og_cwt_merge) on the two synthetic .dta
*  files: the registry + treatment cohort (Table A) and the raw CWT records
*  (Table B). For each patient it picks the single CWT record that matches the
*  derived treatment pathway, takes the earliest valid decision-to-treat (DTT)
*  among the matching records, and derives the waiting-time intervals and the
*  audit categories.
*
*  The logic, in order:
*    1. read and parse the CWT records, assign each a broad modality group
*    2. apply the surgery-01 / 23-24 date-split rule
*    3. keep records within the DTT window of each patient's diagnosis
*    4. keep records whose modality group is consistent with the pathway
*       (only where at least one consistent record exists for that patient)
*    5. anchor on the pathway's primary modality, then the earliest DTT
*    6. attach the anchor to the cohort and derive intervals + audit fields
*
*  Plain Stata only - no Mata, no user-written packages. Set the path below.
*  ---------------------------------------------------------------------------

clear all
set more off
version 14

*  ---------------------------------------------------------------------------
*  Paths and merge constants (kept together so the rules are easy to audit)
*  ---------------------------------------------------------------------------
local base "D:/Projects/#2045_ICON_TACTIC/Project4_OG_variation_deviants/tactic-og-quantifying-variation/Data/synthetic"

local reg_dta "`base'/og_ncras_treatment_synthetic.dta"
local cwt_dta "`base'/og_cwt_records_synthetic.dta"
local out_dta "`base'/og_cohort_synthetic_stata.dta"

local tx_window_days = 270        // DTT must sit within this many days of diagnosis
local dtt_min_offset = -30        // earliest DTT relative to diagnosis (days)
local treat_tol_days = 14         // treatment may precede the DTT by up to this
local surg_switch    = td(01oct2020)   // CWT modality code 01 retired, 23/24 introduced on this date


*  ===========================================================================
*  Step 1-2.  CWT records: parse dates, assign modality group, date-split rule
*  ===========================================================================
use "`cwt_dta'", clear

*  the CWT date fields arrive as "dd/mm/yyyy" strings (as in the raw extract),
*  so parse them to Stata daily dates
gen dtt_date   = date(treat_period_start, "DMY")
gen treat_date = date(treat_start,        "DMY")
gen mdt_d      = date(mdt_date,           "DMY")
format dtt_date treat_date mdt_d %td

*  broad modality group, following the CWT data dictionary
gen str12 mod_group = ""
replace mod_group = "surgery"      if inlist(modality, "01", "23", "24")
replace mod_group = "chemo"        if inlist(modality, "02", "14", "15")
replace mod_group = "hormone"      if modality == "03"
replace mod_group = "chemort"      if modality == "04"
replace mod_group = "radiotherapy" if inlist(modality, "05", "06", "13")
replace mod_group = "palliative"   if inlist(modality, "07", "08", "09")
replace mod_group = "other"        if modality == "97"
replace mod_group = "declined"     if modality == "98"

*  surgery-01 date split: code 01 counts as surgery only before the switch
*  date; 23/24 count as surgery only on or after it. Records that fall on the
*  wrong side of the switch lose their surgery grouping.
replace mod_group = "" if modality == "01" & treat_date >= `surg_switch'
replace mod_group = "" if inlist(modality, "23", "24") & treat_date < `surg_switch'

*  drop records with no usable group, declined records, and records with no DTT
drop if mod_group == "" | mod_group == "declined" | missing(dtt_date)

keep pseudo_patientid dtt_date treat_date mdt_d modality mod_group
tempfile cwt_clean
save "`cwt_clean'"


*  ===========================================================================
*  Step 3.  Bring in the cohort's diagnosis date and pathway, keep in-window
*  ===========================================================================
*  We need, per patient: diagnosis date (for the window) and tx_pathway (for the
*  consistency check). Take them from Table A.
use pseudo_patientid diagmdy tx_pathway first_tx_date using "`reg_dta'", clear
tempfile cohort_keys
save "`cohort_keys'"

use "`cwt_clean'", clear
merge m:1 pseudo_patientid using "`cohort_keys'", keep(match) nogenerate

*  days from diagnosis to the decision-to-treat, then keep the in-window records
gen days_dx_to_dtt = dtt_date - diagmdy
keep if days_dx_to_dtt >= `dtt_min_offset' & days_dx_to_dtt <= `tx_window_days'


*  ===========================================================================
*  Step 4.  Flag records whose modality group is consistent with the pathway
*  ===========================================================================
*  Each pathway has a set of plausible clock-stop modality groups. group_ok = 1
*  marks a record whose group is one of them.
gen byte group_ok = 0

replace group_ok = 1 if tx_pathway == "EMR/ESD only"                  & inlist(mod_group, "surgery", "other")
replace group_ok = 1 if tx_pathway == "EMR/ESD then surgery"          & mod_group == "surgery"
replace group_ok = 1 if tx_pathway == "Surgery + neoadjuvant chemoRT" & inlist(mod_group, "surgery", "chemort", "chemo", "radiotherapy")
replace group_ok = 1 if tx_pathway == "Surgery + neoadjuvant chemo"   & inlist(mod_group, "surgery", "chemo")
replace group_ok = 1 if tx_pathway == "Surgery + neoadjuvant RT"      & inlist(mod_group, "surgery", "radiotherapy", "chemort")
replace group_ok = 1 if tx_pathway == "Surgery + adjuvant chemo"      & inlist(mod_group, "surgery", "chemo")
replace group_ok = 1 if tx_pathway == "Surgery only"                  & mod_group == "surgery"
replace group_ok = 1 if tx_pathway == "Surgery + other"               & inlist(mod_group, "surgery", "other")
replace group_ok = 1 if tx_pathway == "Definitive chemoRT"            & inlist(mod_group, "chemort", "chemo", "radiotherapy")
replace group_ok = 1 if tx_pathway == "Curative RT only"              & inlist(mod_group, "radiotherapy", "chemort")
replace group_ok = 1 if tx_pathway == "Palliative chemo + RT"         & inlist(mod_group, "chemo", "radiotherapy", "chemort", "palliative")
replace group_ok = 1 if tx_pathway == "SACT only"                     & inlist(mod_group, "chemo", "hormone", "palliative")
replace group_ok = 1 if tx_pathway == "Palliative RT only"            & inlist(mod_group, "radiotherapy", "palliative")
replace group_ok = 1 if tx_pathway == "No treatment recorded"         & inlist(mod_group, "palliative", "other")

*  is this record's group the pathway's single PRIMARY (defining) clock-stop?
*  The primary modality wins the tie-break below, so a neoadjuvant patient
*  anchors on their chemo/RT rather than an earlier surgical record.
gen byte is_primary = 0
replace is_primary = 1 if tx_pathway == "EMR/ESD only"                  & mod_group == "surgery"
replace is_primary = 1 if tx_pathway == "EMR/ESD then surgery"          & mod_group == "surgery"
replace is_primary = 1 if tx_pathway == "Surgery + neoadjuvant chemoRT" & mod_group == "chemort"
replace is_primary = 1 if tx_pathway == "Surgery + neoadjuvant chemo"   & mod_group == "chemo"
replace is_primary = 1 if tx_pathway == "Surgery + neoadjuvant RT"      & mod_group == "radiotherapy"
replace is_primary = 1 if tx_pathway == "Surgery + adjuvant chemo"      & mod_group == "surgery"
replace is_primary = 1 if tx_pathway == "Surgery only"                  & mod_group == "surgery"
replace is_primary = 1 if tx_pathway == "Surgery + other"               & mod_group == "surgery"
replace is_primary = 1 if tx_pathway == "Definitive chemoRT"            & mod_group == "chemort"
replace is_primary = 1 if tx_pathway == "Curative RT only"              & mod_group == "radiotherapy"
replace is_primary = 1 if tx_pathway == "Palliative chemo + RT"         & mod_group == "chemo"
replace is_primary = 1 if tx_pathway == "SACT only"                     & mod_group == "chemo"
replace is_primary = 1 if tx_pathway == "Palliative RT only"            & mod_group == "radiotherapy"
replace is_primary = 1 if tx_pathway == "No treatment recorded"         & mod_group == "palliative"

*  Where a patient has at least one pathway-consistent record, keep only those.
*  Where none are consistent, keep all in-window records (the patient still gets
*  an anchor, just not a pathway-matched one).
bysort pseudo_patientid (group_ok): gen byte any_match = group_ok[_N]
drop if any_match == 1 & group_ok == 0


*  ===========================================================================
*  Step 5.  Anchor: primary modality first, then earliest DTT
*  ===========================================================================
*  Sort so that, within each patient, primary records come first (is_primary
*  descending) and then the earliest DTT, then keep the first record per
*  patient. gsort handles the mixed sort directions and leaves the data sorted
*  by pseudo_patientid, so the by-group count below is valid.
gsort pseudo_patientid -is_primary dtt_date
by pseudo_patientid: gen byte pick = (_n == 1)
keep if pick == 1
drop pick

*  rename to the cohort's CWT field names and keep just the anchor fields
rename dtt_date   cwt_dtt_date
rename treat_date cwt_treat_date
rename mdt_d      cwt_mdt_date
rename modality   cwt_modality
keep pseudo_patientid cwt_dtt_date cwt_treat_date cwt_mdt_date cwt_modality
tempfile cwt_anchor
save "`cwt_anchor'"


*  ===========================================================================
*  Step 6.  Attach the anchor to the full cohort and derive the outputs
*  ===========================================================================
use "`reg_dta'", clear
merge 1:1 pseudo_patientid using "`cwt_anchor'", keep(master match) nogenerate

*  --- the six core waiting-time intervals --------------------------------
*  each of diagnosis / endoscopy to each of decision-to-treat / first treatment,
*  plus the DTT -> treatment link
gen wt_endo_to_dx  = diagmdy       - endoscopy_date
gen wt_dx_to_dtt   = cwt_dtt_date  - diagmdy
gen wt_endo_to_dtt = cwt_dtt_date  - endoscopy_date
gen wt_dx_to_tx    = first_tx_date - diagmdy
gen wt_endo_to_tx  = first_tx_date - endoscopy_date
gen wt_dtt_to_tx   = first_tx_date - cwt_dtt_date

*  --- per-modality component intervals (NA where that arm did not occur) --
gen wt_dx_to_surg   = surgery_date - diagmdy
gen wt_endo_to_surg = surgery_date - endoscopy_date
gen wt_dx_to_sact   = sact_date    - diagmdy
gen wt_endo_to_sact = sact_date    - endoscopy_date
gen wt_dx_to_rt     = rt_date      - diagmdy
gen wt_endo_to_rt   = rt_date      - endoscopy_date

*  --- treatment-sequencing gaps (neoadjuvant / adjuvant) -----------------
gen wt_sact_to_surg = surgery_date - sact_date
gen wt_surg_to_sact = sact_date    - surgery_date
gen wt_rt_to_surg   = surgery_date - rt_date

*  --- DTT validity -------------------------------------------------------
*  valid if the DTT is on or after diagnosis and the treatment does not precede
*  it by more than the tolerance. Stata treats a missing value as larger than
*  any number, so the tolerance test is guarded with explicit !missing() checks
*  to stop a missing interval passing by default. EMR/ESD pathways are set to
*  missing because a decision-to-treat is less meaningful there.
gen byte dtt_valid = !missing(cwt_dtt_date) & !missing(wt_dx_to_dtt) ///
    & !missing(wt_dtt_to_tx) & wt_dx_to_dtt >= 0 ///
    & wt_dtt_to_tx >= -`treat_tol_days'
replace dtt_valid = . if inlist(tx_pathway, "EMR/ESD only", "EMR/ESD then surgery")

*  --- audit modality category (NOGCA Tables 3 and 4 groupings) -----------
gen str30 tx_modality_audit = ""
replace tx_modality_audit = "Surgery only"                 if tx_pathway == "Surgery only"
replace tx_modality_audit = "Surgery plus SACT/RT"         if inlist(tx_pathway, ///
    "Surgery + neoadjuvant chemo", "Surgery + neoadjuvant chemoRT", ///
    "Surgery + neoadjuvant RT", "Surgery + adjuvant chemo", "Surgery + other")
replace tx_modality_audit = "EMR/ESD"                      if inlist(tx_pathway, ///
    "EMR/ESD only", "EMR/ESD then surgery")
replace tx_modality_audit = "Definitive chemoRT"           if tx_pathway == "Definitive chemoRT"
replace tx_modality_audit = "Curative RT only"             if tx_pathway == "Curative RT only"
replace tx_modality_audit = "Chemo/RT only (non-curative)" if inlist(tx_pathway, ///
    "Palliative chemo + RT", "SACT only", "Palliative RT only")
replace tx_modality_audit = "No treatment recorded"        if tx_pathway == "No treatment recorded"

*  --- audit intent -------------------------------------------------------
gen str12 tx_intent_audit = ""
replace tx_intent_audit = "Curative"     if inlist(tx_modality_audit, ///
    "Surgery only", "Surgery plus SACT/RT", "EMR/ESD", ///
    "Definitive chemoRT", "Curative RT only")
replace tx_intent_audit = "Non-curative" if tx_modality_audit == "Chemo/RT only (non-curative)"
replace tx_intent_audit = "No treatment" if tx_modality_audit == "No treatment recorded"

*  --- treatment received within the window -------------------------------
*  any treatment counts the CWT treatment date as well as first_tx_date, because
*  first_tx_date is the curative clock-stop (missing for palliative patients): a
*  palliative patient with an in-window CWT treatment still received treatment.
gen days_dx_to_cwt_treat = cwt_treat_date - diagmdy

gen byte received_any_tx = (tx_pathway != "No treatment recorded") & ( ///
       (!missing(first_tx_date) & wt_dx_to_tx <= `tx_window_days') | ///
       (!missing(cwt_treat_date) & days_dx_to_cwt_treat >= 0 ///
            & days_dx_to_cwt_treat <= `tx_window_days') )

gen byte received_curative_tx_audit = (tx_intent_audit == "Curative") ///
    & !missing(first_tx_date) & wt_dx_to_tx <= `tx_window_days'

drop days_dx_to_cwt_treat


*  ===========================================================================
*  Quick checks and save
*  ===========================================================================
display _n "Pathway mix:"
tabulate tx_pathway, missing

display _n "Audit Table 4 (stage 1-3): % curative and % any treatment"
preserve
    keep if inlist(stage_clean, "1", "2", "3")
    summarize received_curative_tx_audit received_any_tx
restore

display _n "CWT coverage (share with a DTT):"
count if !missing(cwt_dtt_date)
display "  " r(N) " of " _N

save "`out_dta'", replace
display _n "Saved `out_dta'"
