/*=============================================================================
  01 - Check input data formats (OG)
  -----------------------------------------------------------------------------
  Check that the two synthetic dataset inputs are in the shape the rest of the
  OG Stata scripts expect. Reports every problem it finds (it does not stop at
  the first), then prints a final summary. Run this before the merge/analysis.

  Checks:
    - both files exist
    - all variables the pipeline reads are present, with the right storage type
    - the patient id is a string in both files and unique in the cohort
    - the registry date fields are numeric Stata dates in a plausible range
    - the CWT date strings parse as dd/mm/YYYY
    - the CWT modality codes are strings holding the expected values
    - the CWT site codes look like oesophago-gastric ICD-10 (C15x / C16x)
    - the merge key overlaps between the two files

  Inputs (in syn):
    og_cohort_precwt_SYNTH.dta    Table A, one row per patient (pre-CWT cohort)
    cwt_records_SYNTH.dta         Table B, raw CWT records (long)
=============================================================================*/

clear all
set more off

* set syn if running this file on its own
if "$syn" == "" global syn "Data/synthetic"

global n_err  = 0
global n_warn = 0

* helper: confirm a variable exists and (optionally) has the expected type
capture program drop chkvar
program define chkvar
    args fname vname vtype
    capture confirm variable `vname'
    if _rc {
        display as error "  FAIL [`fname']: required variable not found - `vname'"
        global n_err = $n_err + 1
        exit
    }
    if "`vtype'" == "string" {
        capture confirm string variable `vname'
        if _rc {
            display as error "  FAIL [`fname']: `vname' should be string but is numeric"
            global n_err = $n_err + 1
        }
    }
    if "`vtype'" == "numeric" {
        capture confirm numeric variable `vname'
        if _rc {
            display as error "  FAIL [`fname']: `vname' should be numeric but is string"
            global n_err = $n_err + 1
        }
    }
    display "  ok [`fname']: `vname'"
end

*==============================================================================
* Table A: pre-CWT cohort (one row per patient)
*==============================================================================
display _n "{hline 70}"
display "Checking Table A: og_cohort_precwt_SYNTH.dta"
display "{hline 70}"

capture confirm file "$syn/og_cohort_precwt_SYNTH.dta"
if _rc {
    display as error "  FAIL: file not found in $syn"
    global n_err = $n_err + 1
}
else {
    use "$syn/og_cohort_precwt_SYNTH.dta", clear

    * merge key and the dates the merge / waiting-time logic reads
    chkvar "A" pseudo_patientid string
    chkvar "A" diagmdy          numeric
    chkvar "A" first_tx_date    numeric
    chkvar "A" endoscopy_date   numeric
    chkvar "A" surgery_date     numeric
    chkvar "A" sact_date        numeric
    chkvar "A" rt_date          numeric

    * cohort / tumour descriptors
    chkvar "A" pseudo_tumourid                  string
    chkvar "A" ydiag                            numeric
    chkvar "A" agediag                          numeric
    chkvar "A" sex                              numeric
    chkvar "A" tumour_site_grp                  string
    chkvar "A" cancer_subtype                   string
    chkvar "A" stage_clean                      string
    chkvar "A" route_combined                   string
    chkvar "A" emergency_admission              numeric
    chkvar "A" ps_num                           numeric
    chkvar "A" ethnicity_group_broad            string
    chkvar "A" NHSE_reversed_imd_quintile_lsoas string

    * organisation / trust fields used by the provider analysis
    chkvar "A" diag_hosp                        string
    chkvar "A" diag_trust                       string
    chkvar "A" first_trust                      string
    chkvar "A" tx_trust                         string
    chkvar "A" change_trust                     numeric
    chkvar "A" lsoa11_code                       string

    * pathway and the headline waiting times
    chkvar "A" tx_pathway                       string
    chkvar "A" wt_dx_to_tx                       numeric
    chkvar "A" wt_endo_to_tx                     numeric

    * survival fields kept as raw registry inputs
    chkvar "A" finmdy                            numeric
    chkvar "A" died                              numeric

    * one row per patient
    capture isid pseudo_patientid
    if _rc {
        display as error "  FAIL [A]: pseudo_patientid is not unique (expected one row per patient)"
        global n_err = $n_err + 1
    }

    * registry dates should be real Stata dates in a sensible range
    capture confirm numeric variable diagmdy
    if !_rc {
        quietly summarize diagmdy
        local dmin = r(min)
        local dmax = r(max)
        display "  diagmdy range: " %td `dmin' " to " %td `dmax'
        if `dmin' < td(01jan2010) | `dmax' > td(31dec2030) {
            display as error "  WARN [A]: diagmdy outside 2010-2030 - is it a true Stata date?"
            global n_warn = $n_warn + 1
        }
    }
    capture confirm numeric variable first_tx_date
    if !_rc {
        quietly summarize first_tx_date
        display "  first_tx_date range: " %td r(min) " to " %td r(max)
    }

    * change_trust should be a 0/1 flag (missing allowed where no curative tx)
    capture confirm numeric variable change_trust
    if !_rc {
        quietly tab change_trust, missing
        quietly count if !inlist(change_trust, 0, 1) & !missing(change_trust)
        if r(N) > 0 {
            display as error "  WARN [A]: change_trust has values other than 0/1/missing"
            global n_warn = $n_warn + 1
        }
    }
}

*==============================================================================
* Table B: raw CWT records
*==============================================================================
display _n "{hline 70}"
display "Checking Table B: cwt_records_SYNTH.dta"
display "{hline 70}"

capture confirm file "$syn/cwt_records_SYNTH.dta"
if _rc {
    display as error "  FAIL: file not found in $syn"
    global n_err = $n_err + 1
}
else {
    use "$syn/cwt_records_SYNTH.dta", clear

    chkvar "B" pseudo_patientid   string
    chkvar "B" site_icd10         string
    chkvar "B" modality           string
    chkvar "B" crtp_date          string
    chkvar "B" date_first_seen    string
    chkvar "B" mdt_date           string
    chkvar "B" treat_period_start string
    chkvar "B" treat_start        string

    * all CWT date strings must parse as dd/mm/YYYY
    foreach d in crtp_date date_first_seen mdt_date treat_period_start treat_start {
        capture confirm string variable `d'
        if !_rc {
            quietly count if !missing(`d')
            local n_nonmiss = r(N)
            quietly gen double _chk = date(`d', "DMY")
            quietly count if missing(_chk) & !missing(`d')
            local n_bad = r(N)
            drop _chk
            if `n_bad' > 0 {
                display as error "  FAIL [B]: `n_bad' of `n_nonmiss' `d' values do not parse as dd/mm/YYYY"
                global n_err = $n_err + 1
            }
            else {
                display "  `d': all non-missing values parse as dd/mm/YYYY"
            }
        }
    }

    * modality should be a string holding the OG treatment codes
    capture confirm string variable modality
    if !_rc {
        display "  modality values (counts):"
        tab modality

        * at least some of the main treatment codes should be present
        local found_keep = 0
        foreach m in 01 02 03 04 {
            quietly count if modality == "`m'"
            if r(N) > 0 local found_keep = 1
        }
        if !`found_keep' {
            display as error "  WARN [B]: none of the main modalities (01-04) found - check coding"
            global n_warn = $n_warn + 1
        }
    }

    * site codes should be oesophago-gastric (C15x / C16x)
    capture confirm string variable site_icd10
    if !_rc {
        quietly count if !missing(site_icd10)
        local n_site = r(N)
        quietly count if inlist(substr(site_icd10, 1, 3), "C15", "C16")
        local n_og = r(N)
        display "  site_icd10: `n_og' of `n_site' records are C15x / C16x"
        if `n_og' == 0 {
            display as error "  WARN [B]: no C15x / C16x site codes found - check site_icd10 coding"
            global n_warn = $n_warn + 1
        }
    }
}

*==============================================================================
* Cross-file: merge key compatibility and overlap
*==============================================================================
display _n "{hline 70}"
display "Checking the merge key across the two files"
display "{hline 70}"

capture confirm file "$syn/og_cohort_precwt_SYNTH.dta"
local haveA = (_rc == 0)
capture confirm file "$syn/cwt_records_SYNTH.dta"
local haveB = (_rc == 0)

if `haveA' & `haveB' {
    use pseudo_patientid using "$syn/og_cohort_precwt_SYNTH.dta", clear
    quietly count
    local nA = r(N)
    tempfile akeys
    quietly save `akeys'

    use pseudo_patientid using "$syn/cwt_records_SYNTH.dta", clear
    quietly bysort pseudo_patientid: keep if _n == 1
    quietly count
    local nB = r(N)

    merge 1:1 pseudo_patientid using `akeys'
    quietly count if _merge == 3
    local n_match = r(N)
    quietly count if _merge == 2
    local n_aonly = r(N)

    display "  Table A patients: `nA'"
    display "  Table B distinct patients: `nB'"
    display "  patients in both files: `n_match'"
    display "  Table A patients with no CWT record: `n_aonly'"

    if `n_match' == 0 {
        display as error "  FAIL: no patient ids overlap - the merge would return nothing"
        global n_err = $n_err + 1
    }
}
else {
    display as error "  skipped (one or both files missing)"
}

*==============================================================================
* Summary
*==============================================================================
display _n "{hline 70}"
if $n_err == 0 & $n_warn == 0 {
    display "Input check passed: no problems found. Safe to run the merge."
}
else {
    display "Input check finished with $n_err error(s) and $n_warn warning(s)."
    if $n_err > 0 display as error "Fix the errors above before running the merge."
}
display "{hline 70}"
