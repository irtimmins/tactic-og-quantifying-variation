********************************************************************************
* Provider-level analysis: mean waiting times by site characteristics
* Replicates R script in Stata
********************************************************************************

********************************************************************************
* Load data
********************************************************************************

* Load cohort
clear
use "E:\Data_PHE\Extracts\#2045_ICON_TACTIC\Derived\og_cohort_cwt_SYNTH.dta"

* Keep only valid waiting time records
keep if !missing(wt_dx_to_tx) & !missing(wt_dx_to_dtt) & !missing(wt_dtt_to_tx)
keep if wt_dx_to_tx >= 0 & wt_dx_to_dtt >= 0 & wt_dtt_to_tx >= 0

********************************************************************************
* Load and prepare provider reference file
********************************************************************************

preserve

    import excel "E:\Data_PHE\Extracts\#2045_ICON_TACTIC\Derived\provider_level\NHSHospitals_services_SYNTH.xlsx", ///
        sheet("Sheet1") firstrow clear

    * Drop rows missing key identifiers
    drop if missing(Trust_Name)
    drop if missing(Hospital_site_code)

    * Exclude flagged trust colours
    drop if inlist(Trust_Name_colour, "Light Red", "Pink Red", "Orange")

    * Within each site code, keep single row or row with bowel surgery info
    destring Bowel_ca_surgery, replace force
    bysort Hospital_site_code: gen n_rows = _N
    drop if n_rows > 1 & missing(Bowel_ca_surgery)
    drop n_rows

    * Rename for merging
    rename Hospital_site_code diag_hosp

    * Binary factors - destring first as Excel imports as string
    destring Comprehensive_centre Teaching_hospitals, replace force

    gen comprehensive = (Comprehensive_centre == 1)
    label define comp 0 "Non-comprehensive" 1 "Comprehensive"
    label values comprehensive comp

    gen teaching = (Teaching_hospitals == 1)
    label define teach 0 "Non-teaching" 1 "Teaching"
    label values teaching teach

    * CQC rating - ordered numeric
    gen cqc_rating = .
    replace cqc_rating = 1 if Latest_Rating == "Inadequate"
    replace cqc_rating = 2 if Latest_Rating == "Requires Improvement"
    replace cqc_rating = 3 if Latest_Rating == "Good"
    replace cqc_rating = 4 if Latest_Rating == "Outstanding"
    label define cqc 1 "Inadequate" 2 "Requires Improvement" 3 "Good" 4 "Outstanding"
    label values cqc_rating cqc

    * Staff engagement and morale quintiles
    destring Staff_engagement Moral, replace force
    xtile staff_eng_cat = Staff_engagement, nquantiles(5)
    xtile moral_cat     = Moral,            nquantiles(5)
    label define quint 1 "Q1 (lowest)" 2 "Q2" 3 "Q3" 4 "Q4" 5 "Q5 (highest)"
    label values staff_eng_cat quint
    label values moral_cat     quint

    * Bed occupancy (mean = bed occupancy rate from Excel)
    destring mean, replace force
    gen bed_occ_cat = .
    replace bed_occ_cat = 0 if mean <  0.95 & !missing(mean)
    replace bed_occ_cat = 1 if mean >= 0.95 & !missing(mean)
    label define bedocc 0 "Normal (<95%)" 1 "High (>=95%)"
    label values bed_occ_cat bedocc

    * Keep only variables needed for merge
    keep diag_hosp comprehensive teaching cqc_rating staff_eng_cat moral_cat bed_occ_cat


    save "E:\Data_PHE\Extracts\#2045_ICON_TACTIC\Derived\provider_level\temp_OG_synth_provider_covariates.dta", replace

restore

********************************************************************************
* Load net gain / competitor status and prepare
********************************************************************************

preserve

    use "E:\Data_PHE\Extracts\#2045_ICON_TACTIC\Derived\provider_level\og_cohort_net_gain_loss_2015_2022_SYNTH.dta", clear

    * competitor_status is already labelled numeric in .dta
    * 1=Winner 2=Loser 3=Insignificant diff.
    rename sitecode diag_hosp
    keep diag_hosp competitor_status

    save "E:\Data_PHE\Extracts\#2045_ICON_TACTIC\Derived\provider_level\temp_og_net_gain_covariates_SYNTH.dta", replace

restore

********************************************************************************
* Merge provider covariates and competitor status onto cohort
********************************************************************************

merge m:1 diag_hosp using "E:\Data_PHE\Extracts\#2045_ICON_TACTIC\Derived\provider_level\temp_OG_synth_provider_covariates.dta", keep(master match) nogenerate

merge m:1 diag_hosp using "E:\Data_PHE\Extracts\#2045_ICON_TACTIC\Derived\provider_level\temp_og_net_gain_covariates_SYNTH.dta", keep(master match) nogenerate
* Check merge
count
tab competitor_status, missing

********************************************************************************
* Mean waiting times by covariate group
*  Produces: mean, SE, N patients and N hospitals per group for each waiting time outcome
********************************************************************************


local wt_vars   wt_dx_to_tx wt_dx_to_dtt wt_dtt_to_tx
local wt_labels "Dx to Treatment" "Dx to DTT" "DTT to Treatment"

local covariates comprehensive teaching cqc_rating staff_eng_cat moral_cat bed_occ_cat competitor_status

foreach wt of local wt_vars {
    foreach cov of local covariates {

        di "*** Outcome: `wt' | Covariate: `cov' ***"

        * Mean, SE and patient N per group
        tabstat `wt', by(`cov') stats(mean semean n) nototal

        * Number of distinct hospitals per group
        di "*** Hospital N: `wt' | Covariate: `cov' ***"
        preserve
            bysort `cov' diag_hosp: keep if _n == 1
            * Create a numeric dummy to count rows with tabstat
            gen byte hosp_count = 1
            tabstat hosp_count, by(`cov') stats(n) nototal
        restore

    }
}

