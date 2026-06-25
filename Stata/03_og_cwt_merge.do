* =============================================================================
* OG cancer - step 3: merge the CWT records onto the derived cohort
* -----------------------------------------------------------------------------
* Reads the derived cohort from step 1 and the raw CWT records, then:
*   1. assign each CWT row a broad modality group, applying the surgery
*      transition-window rule (01 / 23 / 24)
*   2. keep in-window records whose modality group is consistent with the
*      patient's pathway
*   3. anchor each patient on the pathway's primary modality, then the earliest
*      decision-to-treat (DTT) date
*   4. attach the DTT node and derive the waiting-time family and dtt_valid
*   5. derive the broader treatment categories if helpful
*
* Run 02_og_derive_pathway.do first.
*
* Input : Test_data/og_derived_synthetic_stata.dta
*         Test_data/og_cwt_records_synthetic_stata.dta
* Output: Test_data/og_cohort_synthetic_stata.dta
* =============================================================================

clear all
set more off

* -----------------------------------------------------------------------------
* Data directory - hard-coded. Edit this one line if the data moves.
* -----------------------------------------------------------------------------
local in_dir  "D:/Projects/#2045_ICON_TACTIC/Project4_OG_variation_deviants/tactic-og-quantifying-variation/Test_data"
local out_dir "`in_dir'"
local derived "`in_dir'/og_derived_synthetic_stata.dta"
local cwtfile "`in_dir'/og_cwt_records_synthetic_stata.dta"

* the derived cohort is produced by 02
capture confirm file "`derived'"
if _rc di as error "Note: run 02_og_derive_pathway.do first if this stops - derived cohort not found."

* interval / window constants
local tx_window_days  = 270
local dtt_min_offset  = -30
local treat_tol_days  = 14

* In CWT the 01 surgical modality code was retired in 2020.
* Hence for surgery there is a coding changeover: 
* within the transition window 01/23/24 all count as surgery; 
* before it only 01, after it only 23/24. 
local surg_start = td(01jan2020)
local surg_end   = td(30jun2021)

* =============================================================================
* 1.  CWT records: parse, assign modality group, apply the surgery rule
* =============================================================================
use "`cwtfile'", clear

* The CWT .dta stores its dates as dd/mm/YYYY strings; parse them into numeric
* Stata dates under the names the merge uses.
gen double cwt_dtt_date   = date(treat_period_start, "DMY")
gen double cwt_treat_date = date(treat_start,        "DMY")
gen double cwt_mdt_date   = date(mdt_date,           "DMY")
format cwt_dtt_date cwt_treat_date cwt_mdt_date %td

* broad modality group, following the CWT data dictionary
gen str12 mod_group = ""
replace mod_group = "surgery"      if inlist(modality, "01", "23", "24")
replace mod_group = "chemo"        if inlist(modality, "02", "14", "15")
replace mod_group = "hormone"      if modality == "03"
replace mod_group = "chemort"      if modality == "04"
replace mod_group = "radiotherapy" if inlist(modality, "05", "06", "13")
replace mod_group = "palliative"   if inlist(modality, "07", "08", "09")
replace mod_group = "other"        if modality == "97"
replace mod_group = "declined"     if modality == "98"

* transition-window rule: drop 01 after the window, and 23/24 before it, so each
* surgery code only counts in the period where it was the live coding
replace mod_group = "" if modality == "01" & cwt_treat_date > `surg_end'
replace mod_group = "" if inlist(modality, "23", "24") & cwt_treat_date < `surg_start'

* keep usable, dated, non-declined records
drop if mod_group == "" | mod_group == "declined" | missing(cwt_dtt_date)

tempfile cwt_grouped
save `cwt_grouped'

* =============================================================================
* 2.  Identify candidate dtt rows: in-window, and pathway-consistent
* =============================================================================
use "`derived'", clear
keep pseudo_patientid diagmdy tx_pathway first_tx_date
tempfile pw
save `pw'

use `cwt_grouped', clear
merge m:1 pseudo_patientid using `pw', keep(match) nogenerate

gen days_dx_to_dtt = cwt_dtt_date - diagmdy
keep if days_dx_to_dtt >= `dtt_min_offset' & days_dx_to_dtt <= `tx_window_days'

* is this record's modality group a plausible clock-stop for the pathway?
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

* the single defining ("primary") modality group for each pathway, used to break
* ties when a patient has several pathway-consistent records
gen str12 primary_group = ""
replace primary_group = "surgery"      if inlist(tx_pathway, "EMR/ESD only", "EMR/ESD then surgery", "Surgery + adjuvant chemo", "Surgery only", "Surgery + other")
replace primary_group = "chemort"      if inlist(tx_pathway, "Surgery + neoadjuvant chemoRT", "Definitive chemoRT")
replace primary_group = "chemo"        if inlist(tx_pathway, "Surgery + neoadjuvant chemo", "Palliative chemo + RT", "SACT only")
replace primary_group = "radiotherapy" if inlist(tx_pathway, "Surgery + neoadjuvant RT", "Curative RT only", "Palliative RT only")
replace primary_group = "palliative"   if tx_pathway == "No treatment recorded"
gen byte is_primary = mod_group == primary_group

* if a patient has any pathway-consistent record, keep only those; otherwise keep
* all their in-window records (so the earliest DTT can still anchor)
bysort pseudo_patientid: egen byte any_match = max(group_ok)
keep if (any_match == 1 & group_ok == 1) | (any_match == 0)

* =============================================================================
* 3.  Get down to one row per patient: 
* primary modality wins, then the earliest DTT if still multiple records.
* =============================================================================
gsort pseudo_patientid -is_primary cwt_dtt_date
by pseudo_patientid: keep if _n == 1

rename modality cwt_modality
keep pseudo_patientid cwt_dtt_date cwt_mdt_date cwt_treat_date cwt_modality
tempfile anchor
save `anchor'

* =============================================================================
* 4.  Now attach to the derived nrcas + treatment cohort;
* derive all the waiting time variables, run a few checks
* =============================================================================
use "`derived'", clear
merge 1:1 pseudo_patientid using `anchor', keep(master match) nogenerate

* the core intervals: diagnosis / endoscopy to each of DTT and first treatment,
* plus the DTT -> treatment link
gen wt_endo_to_dx  = diagmdy      - endoscopy_date
gen wt_dx_to_dtt   = cwt_dtt_date  - diagmdy
gen wt_endo_to_dtt = cwt_dtt_date  - endoscopy_date
gen wt_dx_to_tx    = first_tx_date - diagmdy
gen wt_endo_to_tx  = first_tx_date - endoscopy_date
gen wt_dtt_to_tx   = first_tx_date - cwt_dtt_date

* per-modality component intervals (which arm the clock stopped on)
gen wt_dx_to_surg   = surgery_date - diagmdy
gen wt_endo_to_surg = surgery_date - endoscopy_date
gen wt_dx_to_sact   = sact_date    - diagmdy
gen wt_endo_to_sact = sact_date    - endoscopy_date
gen wt_dx_to_rt     = rt_date      - diagmdy
gen wt_endo_to_rt   = rt_date      - endoscopy_date

* treatment-sequencing gaps (neoadjuvant / adjuvant)
gen wt_sact_to_surg = surgery_date - sact_date
gen wt_surg_to_sact = sact_date    - surgery_date
gen wt_rt_to_surg   = surgery_date - rt_date

* dtt_valid: DTT on/after diagnosis and treatment not before it beyond tolerance.
* The !missing(first_tx_date) guard matters in Stata: a missing wt_dtt_to_tx would
* otherwise satisfy ">= -tol" (Stata treats missing as +infinity), wrongly scoring
* the non-curative pathways valid. 
gen byte dtt_valid = !missing(cwt_dtt_date) & !missing(first_tx_date) & ///
                     wt_dx_to_dtt >= 0 & wt_dtt_to_tx >= -`treat_tol_days'
replace dtt_valid = . if missing(cwt_dtt_date)
* (double check how DTT works for minimally invasive) replace dtt_valid = . if inlist(tx_pathway, "EMR/ESD only", "EMR/ESD then surgery")

* =============================================================================
* 5.  May be helpful to present broader audit categories
* can decide later whether to have neoadjuvant and adjuvant SACT/RT seperately. 
* =============================================================================
gen str30 tx_modality_audit = ""
replace tx_modality_audit = "Surgery only"          if tx_pathway == "Surgery only"
replace tx_modality_audit = "Surgery plus SACT/RT"  if inlist(tx_pathway, "Surgery + neoadjuvant chemo", "Surgery + neoadjuvant chemoRT", "Surgery + neoadjuvant RT", "Surgery + adjuvant chemo", "Surgery + other")
replace tx_modality_audit = "EMR/ESD"               if inlist(tx_pathway, "EMR/ESD only", "EMR/ESD then surgery")
replace tx_modality_audit = "Definitive chemoRT"    if tx_pathway == "Definitive chemoRT"
replace tx_modality_audit = "Curative RT only"      if tx_pathway == "Curative RT only"
replace tx_modality_audit = "Chemo/RT only (non-curative)" if inlist(tx_pathway, "Palliative chemo + RT", "SACT only", "Palliative RT only")
replace tx_modality_audit = "No treatment recorded" if tx_pathway == "No treatment recorded"

gen str12 tx_intent_audit = ""
replace tx_intent_audit = "Curative" if inlist(tx_modality_audit, "Surgery only", "Surgery plus SACT/RT", "EMR/ESD", "Definitive chemoRT", "Curative RT only")
replace tx_intent_audit = "Non-curative" if tx_modality_audit == "Chemo/RT only (non-curative)"
replace tx_intent_audit = "No treatment" if tx_modality_audit == "No treatment recorded"

* received any treatment within nine months: counts first_tx_date, or the CWT
* treatment date (so palliative patients with an in-window CWT treatment count)
gen byte received_any_tx = 0
replace received_any_tx = 1 if tx_pathway != "No treatment recorded" & ///
    ( (!missing(first_tx_date) & wt_dx_to_tx <= `tx_window_days') | ///
      (!missing(cwt_treat_date) & (cwt_treat_date - diagmdy) >= 0 & ///
                                  (cwt_treat_date - diagmdy) <= `tx_window_days') )

gen byte received_curative_tx_audit = 0
replace received_curative_tx_audit = 1 if tx_intent_audit == "Curative" & ///
    !missing(first_tx_date) & wt_dx_to_tx <= `tx_window_days'

* -----------------------------------------------------------------------------
* Report and save
* -----------------------------------------------------------------------------
di as text _n "Merged cohort: " _N " patients"
quietly count if !missing(cwt_dtt_date)
di as text "CWT coverage (with a DTT): " %4.1f 100*r(N)/_N "%"

di as text _n "Audit Table 4 (stage 1-3): % curative and % any treatment"
preserve
    keep if inlist(stage_clean, "1", "2", "3")
    quietly summarize received_curative_tx_audit
    local pc_cur = 100*r(mean)
    quietly summarize received_any_tx
    local pc_any = 100*r(mean)
    di as text "  curative:  " %3.0f `pc_cur' "%"
    di as text "  any tx:    " %3.0f `pc_any' "%"
restore

save "`out_dir'/og_cohort_synthetic_stata.dta", replace
di as result _n "Saved og_cohort_synthetic_stata.dta (" _N " patients)."
