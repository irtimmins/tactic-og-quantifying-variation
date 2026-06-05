********************************************************************************
* Win/loss analysis on OG cancer cohort
* Inputs:  OG_cohort.dta   (patients + demographics + diag_hosp)
*          pairwise_distance_matrix.dta (lsoa x site x drive time)
* Output:  net_gain_loss.dta  (site-level win/loss classification)
* Should produce identical classification to R script
********************************************************************************

********************************************************************************
* Load OG cancer patient cohort with lsoa11_code and diag_hosp.
********************************************************************************

clear
use "E:\Data_PHE\Extracts\#2045_ICON_TACTIC\Derived\og_cohort_precwt_SYNTH.dta"
* Keep only variables containing "hosp" or "diag"
* Can also keep waiting times variables e.g. "wt_" if needed
keep pseudo_patientid *hosp* *diag* wt_*  lsoa11_code

* Restrict to valid hospital codes.
merge m:1 diag_hosp using "E:\Data_PHE\Extracts\#2045_ICON_TACTIC\Derived\Provider_level\NHSHospitals_valid_sites_SYNTH.dta", keep(match) nogenerate

* Sense checks
distinct diag_hosp
count
di "N patients after provider exclusion: " r(N)
tab diag_hosp, missing

********************************************************************************
* Join patients to distance matrix
* Creates long format: one row per patient x site combination
********************************************************************************

* 1. Note which hospitals are actually in the cohort
preserve
    keep diag_hosp
    duplicates drop
    rename diag_hosp sitecode          // match the matrix's key name
    tempfile cohort_sites
    save `cohort_sites'
restore

* 2. Load the matrix, keep only those sites, then join to patients
preserve
    use "E:\Data_PHE\Extracts\#2045_ICON_TACTIC\Derived\provider_level\OG_pairwise_distance_matrix_SYNTH.dta", clear
    merge m:1 sitecode using `cohort_sites', keep(match) nogenerate   // inner join on sitecode
    tempfile dm
    save `dm'
restore

joinby lsoa11_code using `dm'
count
* Expect rows = (patients) x (distinct cohort sites)  

********************************************************************************
* Identify each patient's nearest hospital
* Tiebreaker 1: prefer diagnosis hospital
* Tiebreaker 2: prefer highest volume site (by unique patients diagnosed)
********************************************************************************

* Flag nearest hospital per patient (exact ties handled in steps below)
bysort pseudo_patientid: egen min_drive_time = min(total_drive_time)
gen byte nearest_hospital = (abs(total_drive_time - min_drive_time) < 0.00005)
label var nearest_hospital "Is this the patient's nearest hospital? (1=Yes, 0=No)"

bysort pseudo_patientid: egen n_nearest = total(nearest_hospital)
tab n_nearest

*-------------------------------------------------------------------------------
* Tiebreaker 1: where two hospitals are equidistant, prefer diag_hosp
*-------------------------------------------------------------------------------

bysort pseudo_patientid: replace nearest_hospital = 0 ///
    if n_nearest == 2 & sitecode != diag_hosp & nearest_hospital == 1

drop n_nearest
bysort pseudo_patientid: egen n_nearest = total(nearest_hospital)
tab n_nearest

*-------------------------------------------------------------------------------
* Tiebreaker 2: for remaining unresolved patients (diag_hosp not among tied
* sites), prefer the site with most unique patients diagnosed there
*-------------------------------------------------------------------------------

bysort pseudo_patientid: gen tag = (_n == 1)
bysort diag_hosp: egen site_volume = total(tag)
drop tag

* Pick row with shortest drive time, breaking ties by highest site volume
bysort pseudo_patientid (total_drive_time -site_volume): replace nearest_hospital = 1 ///
    if n_nearest == 0 & _n == 1

drop site_volume

*-------------------------------------------------------------------------------
* Final check: every patient must have exactly 1 nearest hospital
*-------------------------------------------------------------------------------

drop n_nearest
bysort pseudo_patientid: egen n_nearest = total(nearest_hospital)
tab n_nearest

count if n_nearest != 1
assert n_nearest == 1

********************************************************************************
* Verify diagnosing hospital is in the distance matrix
* Otherwise make further exclusions.
********************************************************************************

gen byte diaghosp_in_dm = (sitecode == diag_hosp)
bysort pseudo_patientid: egen pat_diaghosp_in_dm = max(diaghosp_in_dm)
label var pat_diaghosp_in_dm "Is patient's diagnosing hospital in distance matrix? (1=Yes, 0=No)"

* Exclude patients whose diagnosis hospital is not in the distance matrix
* Hopefully this isn't loads of patients
drop if pat_diaghosp_in_dm == 0

* Check
bysort pseudo_patientid: gen tag = (_n == 1)
count if tag == 1
drop tag


********************************************************************************
* Classify each patient-site row
********************************************************************************

gen byte diagnosed_here = (sitecode == diag_hosp)
label var diagnosed_here "Was patient diagnosed at this site? (1=Yes, 0=No)"

* Drop rows irrelevant to classification
* Keep only: nearest site row AND diagnosing site row per patient
drop if diagnosed_here == 0 & nearest_hospital == 0

* Sanity check: each patient should have 1 or 2 rows only
bysort pseudo_patientid: gen n_rows = _N
tab n_rows
* 1 = core patient (diagnosed at nearest site)
* 2 = leaver or arriver (bypassed or attracted)
assert n_rows <= 2
drop n_rows

gen byte leaver      = (diagnosed_here == 0 & nearest_hospital == 1)
gen byte arriver     = (diagnosed_here == 1 & nearest_hospital == 0)
gen byte core_patient = (diagnosed_here == 1 & nearest_hospital == 1)

label var leaver       "Did patient bypass this (nearest) site? (1=Yes, 0=No)"
label var arriver      "Did patient travel past a closer site to come here? (1=Yes, 0=No)"
label var core_patient "Was patient diagnosed at their nearest site? (1=Yes, 0=No)"

********************************************************************************
* Aggregate to site level
********************************************************************************

bysort sitecode: egen n_leavers  = total(leaver)
bysort sitecode: egen n_arrivers = total(arriver)
bysort sitecode: egen n_core     = total(core_patient)
gen n_net_gain = n_arrivers - n_leavers

label var n_leavers  "No. patients who bypassed this site for diagnosis elsewhere"
label var n_arrivers "No. patients attracted here despite a closer site existing"
label var n_core     "No. patients diagnosed at their nearest site"
label var n_net_gain "Net patient gain (n_arrivers - n_leavers)"

duplicates drop sitecode, force
keep sitecode n_leavers n_arrivers n_core n_net_gain
sort sitecode

count
* Should equal number of unique sites 

* Sanity check: total leavers must equal total arrivers
quietly summarize n_leavers
scalar total_leave = r(sum)
quietly summarize n_arrivers
scalar total_arrive = r(sum)
di "Total leavers: "  total_leave  " | Total arrivers: " total_arrive
assert total_leave == total_arrive

********************************************************************************
* Test for significant difference between arrivers and leavers
* Conditional Poisson test
* iri always tests first > second so swap arguments when leavers > arrivers
********************************************************************************

gen p_value = .
label var p_value "One-sided p-value: significant imbalance between arrivers and leavers?"

gen N = _N
forvalues i = 1/`=N' {
    if n_arrivers[`i'] >= n_leavers[`i'] {
        iri `=n_arrivers[`i']' `=n_leavers[`i']' 1 1
    }
    else {
        iri `=n_leavers[`i']' `=n_arrivers[`i']' 1 1
    }
    replace p_value = r(p) in `i'
}
drop N

********************************************************************************
* Classify each site's competitive status
********************************************************************************

gen byte competitor_status = .
replace  competitor_status = 1 if n_net_gain >  0 & p_value <= 0.05
replace  competitor_status = 2 if n_net_gain <  0 & p_value <= 0.05
replace  competitor_status = 3 if p_value >  0.05

label define competitor_status_lbl 1 "Winner" 2 "Loser" 3 "Insignificant diff."
label values competitor_status competitor_status_lbl
label var competitor_status "Site competition classification based on net patient gain"

* Final output 
tab competitor_status, missing

* Summary to compare directly with R output
tabstat n_leavers n_arrivers n_net_gain, ///
    stats(sum mean median min max) columns(statistics)

********************************************************************************
* Save the provider-level results
********************************************************************************

save "E:\Data_PHE\Extracts\#2045_ICON_TACTIC\Derived\provider_level\og_cohort_net_gain_loss_2015_2022_SYNTH.dta", replace

