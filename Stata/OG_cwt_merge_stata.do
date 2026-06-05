*==============================================================================
* OG Cancer Waiting Times -- CWT merge (Stata translation)
*------------------------------------------------------------------------------
* Joins the Cancer Waiting Times (CWT) decision-to-treat (DTT) node onto the
* pre-CWT cohort, derives the DTT interval metrics, flags validity, and prints 
* the same validation tables.
*
* Works on either the SYNTHETIC or the REAL data -- just point the file names
* below at whichever pair you want to test.
*
*------------------------------------------------------------------------------
* Notes on variable types after that export:
*   - Cohort DATE columns (diagmdy, endoscopy_date, surgery_date, sact_date,
*     rt_date, first_tx_date, finmdy) arrive as Stata daily dates (%td).
*   - CWT DATE columns (treat_period_start, treat_start, crtp_date,
*     date_first_seen, mdt_date) arrive as STRINGS "DD/MM/YYYY" and are parsed
*     below with date(...,"DMY").
*   - tx_pathway, route_combined, stage_clean, IMD arrive as STRINGS.
*==============================================================================

clear all
set more off

*------------------------------------------------------------------------------
* 0. Settings -- swap these three filenames for the real data
*------------------------------------------------------------------------------

global base_dir "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/"

local cohort_file "${base_dir}og_cohort_precwt_SYNTH.dta" // Table A: one row per patient
local cwt_file    "${base_dir}cwt_records_SYNTH.dta" // Table B: many rows per patient
local out_file    "${base_dir}og_cohort_cwt_SYNTH.dta"  // merged output

local txwin = 270        // treatment window in days (9 months post-diagnosis)

tempfile cwt_filtered anchor    // temporary datasets used along the way


*==============================================================================
* 1. Read the CWT records, parse dates, keep OG sites, drop excluded modalities
*==============================================================================
use "`cwt_file'", clear

* --- Parse the string dates into Stata daily dates ---------------------------
* date("31/12/2020","DMY") returns the number of days since 01jan1960.
* The %td format just makes it print as a calendar date; the underlying value
* is an integer, so subtracting two dates gives a number of days.
gen cwt_dtt_date   = date(treat_period_start, "DMY")   // decision-to-treat (DTT)
gen cwt_treat_date = date(treat_start,        "DMY")   // first treatment
gen cwt_mdt_date   = date(mdt_date,           "DMY")   // MDT meeting
format cwt_dtt_date cwt_treat_date cwt_mdt_date %td

* --- Keep only oesophago-gastric sites (ICD-10 C15x or C16x) -----------------
keep if inlist(substr(site_icd10, 1, 3), "C15", "C16")

* --- Drop modalities that should never anchor a clock ------------------------
* 97/98/99 are non-treatment / unknown codes.
drop if inlist(modality, "97", "98", "99")
* My understanding is that endoscopic codes 23/24 were only introduced mid-2020; 
* before then they are unreliable, so drop any 23/24 record with a treatment 
* date before 01jun2020.
drop if inlist(modality, "23", "24") & cwt_treat_date < td(01jun2020)

* --- Diagnostic: how many CWT records per patient? ---------------------------
preserve
    bysort pseudo_patientid: gen _k = _N      // _N = number of rows in the group
    bysort pseudo_patientid: keep if _n == 1  // one row per patient
    di _n(1) "--- CWT records per patient ---"
    tab _k
restore

* --- Diagnostic: modality distribution ---------------------------------------
di _n(1) "--- Modality distribution ---"
tab modality

save "`cwt_filtered'"


*==============================================================================
* 2. Build the CWT anchor: earliest valid DTT per patient within the window
*==============================================================================
* Start from the cohort (one row per patient) so we only keep cohort patients,
* then attach their CWT records (a 1-to-many join: one cohort row can match
* several CWT rows).
use "`cohort_file'", clear
keep pseudo_patientid diagmdy
merge 1:m pseudo_patientid using "`cwt_filtered'"

* _merge==1 : cohort patient with NO matching CWT record  -> drop
* _merge==2 : CWT record for a non-cohort patient          -> drop
* _merge==3 : matched                                       -> keep
drop if _merge == 2
drop _merge

* --- Days from diagnosis to DTT ----------------------------------------------
gen days_dx_to_dtt = cwt_dtt_date - diagmdy

* --- Keep records inside the window ------------------------------------------
keep if !missing(days_dx_to_dtt) & days_dx_to_dtt >= -30 & days_dx_to_dtt <= `txwin'

* --- Keep the EARLIEST DTT per patient ---------------------------------------
bysort pseudo_patientid (cwt_dtt_date): keep if _n == 1

keep pseudo_patientid cwt_dtt_date cwt_treat_date cwt_mdt_date modality days_dx_to_dtt
save "`anchor'"

* --- Report anchor counts / completeness -------------------------------------
count
di "CWT anchor patients (real target for ~36,197): " r(N)
quietly count if !missing(cwt_dtt_date)
di "DTT completeness: " %4.1f 100*r(N)/_N "%"
quietly count if !missing(cwt_mdt_date)
di "MDT completeness: " %4.1f 100*r(N)/_N "%  (real ~40%)"
di _n(1) "days_dx_to_dtt (real: median 39, IQR 24-56):"
summarize days_dx_to_dtt, detail


*==============================================================================
* 3. Validation against the pipeline treatment dates
*==============================================================================
* Attach the cohort's own treatment dates to the anchor (one row per patient).
use "`cohort_file'", clear
keep pseudo_patientid diagmdy first_tx_date surgery_date sact_date rt_date tx_pathway
merge 1:1 pseudo_patientid using "`anchor'", keep(match) nogen

* Three comparison intervals:
gen dtt_to_cwt_treat = cwt_treat_date - cwt_dtt_date   // internal CWT: DTT -> treat
gen dtt_to_tx        = first_tx_date  - cwt_dtt_date   // DTT -> our curative tx
gen cwt_vs_first_tx  = cwt_treat_date - first_tx_date  // CWT treat vs our tx

* --- DTT to CWT treat date (internal consistency) ----------------------------
di _n(1) "--- DTT to CWT treat date (real: median 11, IQR 3-18, no negatives) ---"
summarize dtt_to_cwt_treat, detail
quietly count if dtt_to_cwt_treat < 0 & !missing(dtt_to_cwt_treat)
di "negatives: " r(N)

* --- DTT to first_tx_date -----------------------------------------------------
di _n(1) "--- DTT to first_tx_date (real: median 14, IQR 7-22, 5.3% negative) ---"
summarize dtt_to_tx, detail
quietly count if dtt_to_tx < 0 & !missing(dtt_to_tx)
quietly count if !missing(dtt_to_tx)
local denom = r(N)
quietly count if dtt_to_tx < 0 & !missing(dtt_to_tx)
di "negative dtt_to_tx: " %4.1f 100*r(N)/`denom' "%"

* --- CWT treat date vs first_tx_date -----------------------------------------
* Build 0/100 indicators so their MEAN is directly a percentage.
gen exact100 = 100*(cwt_vs_first_tx == 0)            if !missing(cwt_vs_first_tx)
gen w5_100   = 100*(abs(cwt_vs_first_tx) <= 5)       if !missing(cwt_vs_first_tx)
gen w14_100  = 100*(abs(cwt_vs_first_tx) <= 14)      if !missing(cwt_vs_first_tx)
di _n(1) "--- CWT treat vs first_tx (real: 71.1% exact, 85.6% within 14d) ---"
tabstat exact100 w5_100 w14_100, stat(mean n) columns(statistics)
di "median / p25 / p75 of cwt_vs_first_tx:"
summarize cwt_vs_first_tx, detail

* --- Negative dtt_to_tx by pathway -------------------------------------------
* mean of this 0/100 indicator within each pathway = % negative in that pathway.
gen negtx100 = 100*(dtt_to_tx < 0) if !missing(dtt_to_tx)
di _n(1) "--- Negative dtt_to_tx by pathway (real: EMR/ESD then surgery ~50%) ---"
tabstat negtx100, by(tx_pathway) stat(mean n)


*==============================================================================
* 4. Merge the DTT node onto the FULL cohort + derive intervals & validity
*==============================================================================
use "`cohort_file'", clear
merge 1:1 pseudo_patientid using "`anchor'", ///
    keepusing(cwt_dtt_date cwt_mdt_date cwt_treat_date) nogen

* --- DTT-based waiting-time intervals ----------------------------------------
gen wt_endo_to_dtt = cwt_dtt_date  - endoscopy_date   // staging interval
gen wt_dtt_to_tx   = first_tx_date - cwt_dtt_date     // scheduling interval
gen wt_dx_to_dtt   = cwt_dtt_date  - diagmdy

*------------------------------------------------------------------------------
* --- dtt_valid flag ----------------------
*------------------------------------------------------------------------------
*
* Target three-valued logic:
*   - no CWT DTT                                   -> 0  (FALSE)
*   - CWT DTT present but before diagnosis         -> 0  (FALSE)
*   - CWT DTT present, on/after dx, tx gap KNOWN   -> 1 if gap>=-14 else 0
*   - CWT DTT present, on/after dx, tx gap MISSING -> . (missing, = R's NA)
*   - EMR/ESD pathways (DTT not meaningful)        -> . (missing, = R's NA)

gen byte dtt_valid = .       // start everything as missing

* (a) No CWT DTT at all -> not valid (in R, FALSE & anything = FALSE)
replace dtt_valid = 0 if missing(cwt_dtt_date)

* (b) CWT DTT present but it falls BEFORE diagnosis -> not valid
replace dtt_valid = 0 if !missing(cwt_dtt_date) & wt_dx_to_dtt < 0

* (c) CWT DTT present, on/after diagnosis, AND we know the DTT->treatment gap:
*     valid if treatment is on/after DTT (allowing a 14-day tolerance).
*     The !missing(wt_dtt_to_tx) guard is what reproduces R's NA: if the gap is
*     unknown (no curative treatment) we leave dtt_valid missing.
replace dtt_valid = (wt_dtt_to_tx >= -14)               ///
    if !missing(cwt_dtt_date) & wt_dx_to_dtt >= 0 & !missing(wt_dtt_to_tx)

* (d) EMR/ESD pathways: DTT is not a meaningful clock start -> force missing.
replace dtt_valid = . if inlist(tx_pathway, "EMR/ESD only", "EMR/ESD then surgery")

label define validlbl 0 "not valid" 1 "valid"
label values dtt_valid validlbl

save "`out_file'", replace
di _n(1) "Saved `out_file'"


*==============================================================================
* 5. Post-merge summaries
*==============================================================================

* --- dtt_valid by pathway (row % = % valid within each pathway) --------------
* "if !missing(cwt_dtt_date)" matches the R filter; tab drops the missing
* dtt_valid rows (the EMR/ESD pathways) automatically, as R's !is.na() did.
di _n(1) "--- dtt_valid by pathway (real: curative RT ~99.6%, neoadj chemo ~95.7%) ---"
tab tx_pathway dtt_valid if !missing(cwt_dtt_date), row

* --- Staging (endo->DTT) and scheduling (DTT->tx) intervals ------------------
* Again, guard EVERY ">=0" with !missing() so missing intervals are excluded
* rather than treated as large positive numbers.
di _n(1) "--- Intervals where dtt_valid==1 (real: endo->DTT med 44, DTT->tx med 15) ---"
tabstat wt_endo_to_dtt wt_dtt_to_tx                                          ///
    if dtt_valid == 1                                                        ///
    & !missing(wt_endo_to_dtt) & wt_endo_to_dtt >= 0                         ///
    & !missing(wt_dtt_to_tx)   & wt_dtt_to_tx   >= 0,                        ///
    stat(n p25 p50 p75) columns(statistics)

* --- Intervals by pathway ----------------------------------------------------
di _n(1) "--- Intervals by pathway ---"
tabstat wt_endo_to_dtt wt_dtt_to_tx                                          ///
    if dtt_valid == 1                                                        ///
    & !missing(wt_endo_to_dtt) & wt_endo_to_dtt >= 0                         ///
    & !missing(wt_dtt_to_tx)   & wt_dtt_to_tx   >= 0,                        ///
    by(tx_pathway) stat(n p50) columns(statistics)

* --- Deprivation gradient (staging vs scheduling component) ------------------
di _n(1) "--- Deprivation gradient (real: endo->tx ~60 across IMD, very flat) ---"
tabstat wt_endo_to_tx wt_endo_to_dtt wt_dtt_to_tx                            ///
    if dtt_valid == 1                                                        ///
    & !missing(wt_endo_to_dtt) & wt_endo_to_dtt >= 0                         ///
    & !missing(wt_dtt_to_tx)   & wt_dtt_to_tx   >= 0                         ///
    & !missing(NHSE_reversed_imd_quintile_lsoas),                           ///
    by(NHSE_reversed_imd_quintile_lsoas) stat(n p50) columns(statistics)

* --- Missing-reason summary --------------------------------------------------
* The order matters: each line only fills
* rows still left blank ("" ), so earlier rules take precedence.
gen str40 missing_reason = ""
replace missing_reason = "no endoscopy, no treatment"  ///
    if missing(endoscopy_date) & missing(first_tx_date)
replace missing_reason = "no endoscopy, has treatment" ///
    if missing(endoscopy_date) & !missing(first_tx_date)
replace missing_reason = "has endoscopy, no treatment" ///
    if !missing(endoscopy_date) & missing(first_tx_date)
replace missing_reason = "negative wait"               ///
    if missing_reason == "" & !missing(wt_endo_to_tx) & wt_endo_to_tx < 0
replace missing_reason = "complete"                    ///
    if missing_reason == ""
di _n(1) "--- Missing reason summary ---"
* NOTE: 'has endoscopy, no treatment' is large by design -- first_tx_date is
* curative-only, so all palliative / SACT-only / no-treatment patients land here.
tab missing_reason
