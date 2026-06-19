* =============================================================================
* OG cancer - step 2: derive the treatment pathway
* -----------------------------------------------------------------------------
* Reads the raw synthetic cohort (the per-patient treatment dates, curative
* descriptors, chemo provenance and provider codes) and derives, in order:
*   - the treatment-presence flags    (had_surgery, had_sact, ...)
*   - the sequencing flags            (sact_before_surgery, concurrent_chemo_rt)
*   - tx_pathway                       the treatment-pathway classification
*   - first_tx_date                    the clock-stop date for that pathway
*   - tx_trust                         the provider of the clock-stop treatment
*
* This mirrors og_derive_pathway() in the R pipeline so the same classification
* runs on either platform. The pathway is a function of the flags and dates
* alone; nothing is pre-supplied. Step 3 (03_og_cwt_merge.do) merges the CWT
* records onto the result.
*
* Input : Test_data/og_ncras_treatment_synthetic_stata.dta
* Output: Test_data/og_derived_synthetic_stata.dta
* =============================================================================

clear all
set more off

* -----------------------------------------------------------------------------
* Data directory - hard-coded. Edit this one line if the data moves.
* -----------------------------------------------------------------------------
local in_dir  "D:/Projects/#2045_ICON_TACTIC/Project4_OG_variation_deviants/tactic-og-quantifying-variation/Test_data"
local out_dir "`in_dir'"

use "`in_dir'/og_ncras_treatment_synthetic_stata.dta", clear

* -----------------------------------------------------------------------------
* Treatment-presence flags
* -----------------------------------------------------------------------------
gen byte had_emresd           = !missing(emresd_date)
gen byte had_surgery          = !missing(surgery_date)
gen byte had_curative_surgery = !missing(surgery_date) & curative_surgery == 1
gen byte had_sact             = !missing(sact_date)
gen byte had_rt               = !missing(rt_date)
gen byte had_curative_rt      = !missing(rt_date) & rt_curative == 1
gen byte had_palliative_rt    = !missing(rt_date) & rt_curative == 0

* chemo provenance is optional (the HES-supplement guard). If a column is absent,
* create a neutral placeholder so the chemoRT guard below simply falls back to
* "SACT chemo always counts" - matching the R function's defaulting.
capture confirm variable chemo_source
if _rc gen str4 chemo_source = ""
capture confirm variable hes_chemo_date
if _rc gen double hes_chemo_date = .

* chemo eligible to define a non-surgical definitive-chemoRT pathway: SACT chemo
* always counts; HES-only chemo only when within 28 days of the curative RT, so
* a separate HES chemo episode cannot manufacture chemoRT.
gen byte had_chemo_for_chemort = had_sact & ///
    ( chemo_source != "hes" | ///
      ( !missing(hes_chemo_date) & !missing(rt_date) & ///
        abs(hes_chemo_date - rt_date) <= 28 ) )

* -----------------------------------------------------------------------------
* Sequencing flags
* -----------------------------------------------------------------------------
gen byte sact_before_surgery = had_sact & had_surgery & sact_date < surgery_date
gen byte sact_after_surgery  = had_sact & had_surgery & sact_date > surgery_date
gen byte rt_before_surgery   = had_rt   & had_surgery & rt_date   < surgery_date
gen byte rt_after_surgery    = had_rt   & had_surgery & rt_date   > surgery_date
gen byte concurrent_chemo_rt = had_sact & had_curative_rt & ///
                               abs(sact_date - rt_date) <= 14

gen byte received_curative_tx = had_emresd | had_curative_surgery | had_curative_rt

* -----------------------------------------------------------------------------
* tx_pathway - first matching rule wins (mirrors the R case_when ladder)
* -----------------------------------------------------------------------------
gen str40 tx_pathway = ""
replace tx_pathway = "EMR/ESD only" ///
    if tx_pathway == "" & had_emresd & !had_surgery & !had_sact & !concurrent_chemo_rt
replace tx_pathway = "EMR/ESD then surgery" ///
    if tx_pathway == "" & had_emresd & had_surgery
replace tx_pathway = "Surgery + neoadjuvant chemoRT" ///
    if tx_pathway == "" & had_surgery & sact_before_surgery & rt_before_surgery
replace tx_pathway = "Surgery + neoadjuvant chemo" ///
    if tx_pathway == "" & had_surgery & sact_before_surgery & !rt_before_surgery
replace tx_pathway = "Surgery + neoadjuvant RT" ///
    if tx_pathway == "" & had_surgery & rt_before_surgery & !sact_before_surgery
replace tx_pathway = "Surgery + adjuvant chemo" ///
    if tx_pathway == "" & had_surgery & sact_after_surgery & !sact_before_surgery
replace tx_pathway = "Surgery only" ///
    if tx_pathway == "" & had_surgery & !had_sact & !concurrent_chemo_rt
replace tx_pathway = "Surgery + other" ///
    if tx_pathway == "" & had_surgery
replace tx_pathway = "Definitive chemoRT" ///
    if tx_pathway == "" & !had_surgery & had_curative_rt & had_chemo_for_chemort
replace tx_pathway = "Curative RT only" ///
    if tx_pathway == "" & !had_surgery & had_curative_rt & !had_chemo_for_chemort
replace tx_pathway = "Palliative chemo + RT" ///
    if tx_pathway == "" & !had_surgery & had_palliative_rt & had_sact
replace tx_pathway = "SACT only" ///
    if tx_pathway == "" & !had_surgery & had_sact & !had_curative_rt
replace tx_pathway = "Palliative RT only" ///
    if tx_pathway == "" & !had_surgery & had_palliative_rt & !had_sact
replace tx_pathway = "No treatment recorded" ///
    if tx_pathway == ""

* -----------------------------------------------------------------------------
* first_tx_date - the clock-stop date for the pathway
* neoadjuvant RT / chemoRT set the date; neoadjuvant chemo alone does not (the
* surgery does). min() over the relevant dates where two can apply.
* -----------------------------------------------------------------------------
gen first_tx_date = .
format first_tx_date %td

replace first_tx_date = emresd_date ///
    if inlist(tx_pathway, "EMR/ESD only", "EMR/ESD then surgery")
replace first_tx_date = min(sact_date, rt_date) ///
    if tx_pathway == "Surgery + neoadjuvant chemoRT"
replace first_tx_date = rt_date    if tx_pathway == "Surgery + neoadjuvant RT"
replace first_tx_date = sact_date  if tx_pathway == "Surgery + neoadjuvant chemo"
replace first_tx_date = surgery_date ///
    if inlist(tx_pathway, "Surgery + adjuvant chemo", "Surgery only", "Surgery + other")
replace first_tx_date = min(sact_date, rt_date) ///
    if tx_pathway == "Definitive chemoRT"
replace first_tx_date = rt_date    if tx_pathway == "Curative RT only"
* all other (non-curative) pathways leave first_tx_date missing

* -----------------------------------------------------------------------------
* tx_trust - provider of the clock-stop treatment (first 3 chars of the code)
* surgical / EMR pathways take the surgery provider; RT-anchored pathways take
* the RT provider. SACT's own provider is never the trust source - neoadjuvant
* chemo still takes surgery's, because the curative act is the surgery.
* -----------------------------------------------------------------------------
gen str3 tx_trust = ""
replace tx_trust = substr(surgery_provider, 1, 3) ///
    if inlist(tx_pathway, "EMR/ESD only", "EMR/ESD then surgery", ///
                          "Surgery + neoadjuvant chemo", "Surgery + adjuvant chemo", ///
                          "Surgery only", "Surgery + other")
replace tx_trust = substr(rt_provider, 1, 3) ///
    if inlist(tx_pathway, "Surgery + neoadjuvant chemoRT", "Surgery + neoadjuvant RT", ///
                          "Definitive chemoRT", "Curative RT only")
* non-curative pathways leave tx_trust empty

* -----------------------------------------------------------------------------
* Report the pathway mix and save
* -----------------------------------------------------------------------------
di as text _n "Derived pathway mix:"
tabulate tx_pathway, sort

save "`out_dir'/og_derived_synthetic_stata.dta", replace
di as result _n "Saved og_derived_synthetic_stata.dta (" _N " patients). Next: 03_og_cwt_merge.do"
