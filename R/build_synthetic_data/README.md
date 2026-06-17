# OG minimal synthetic pipeline

A compact synthetic-data pipeline for the oesophago-gastric (OG) cancer
waiting-times work. It builds a synthetic registry + treatment cohort and a
synthetic Cancer Waiting Times (CWT) records table, then merges them on exactly
the same logic the real pipeline uses, reproducing the NOGCA audit pathways and
categorisations. The same minimal merge also runs on the real ICON cohort once
it is condensed to the minimal column set, so the merge code is written once and
shared between synthetic and real data.

The design follows the lighter colon example rather than the full 128-column OG
cohort: a minimal Table A carries identity, a handful of patient/tumour
descriptors, the treatment anchor dates, the derived pathway, and survival -
enough to drive the merge and the audit tables, nothing more.

## Scripts and run order

The scripts are numbered in dependency order.

1. `01_og_minimal_merge.R` - the shared engine. Defines the merge constants,
   the modality/pathway lookups, `og_cwt_merge()` (the linkage + audit
   categorisation), the minimal Table A column contract (`og_minimal_cols`), and
   `condense_icon_to_minimal()` for reducing the real cohort. This script is
   sourced by 03 and 04; it is not run on its own.

2. `02_og_profile.R` - runs on the secure server against the real post-merge
   cohort (`og_cohort_cwt_2015_2022.rds`). Extracts disclosure-safe aggregate
   distributions - patient/tumour marginals, the pathway mix by stage x subtype,
   treatment-timing intervals, trust/hospital volume, and the CWT per-record
   structure - and writes `og_profile_for_synthetic.rds` plus `og_minimal_spec.rds`.
   The profile is what makes the synthetic data resemble the real cohort; without
   it the generator falls back to built-in defaults.

3. `03_og_generate.R` - runs anywhere, no real data needed. Sources 01, reads the
   profile/spec if present, and builds the synthetic Table A (with pathway
   assigned conditional on stage x subtype) and Table B (CWT records with
   modality consistent with each patient's pathway), then runs `og_cwt_merge()`.
   Saves the three outputs as `.rds` and `.dta`.

4. `04_og_validate.R` - sources 01, reads the synthetic outputs and the profile,
   and prints conformance, internal-consistency, pathway-mix and audit-target
   checks.

Typical flow: run 02 on the server to refresh the profile, take the profile out
through output checking, then run 03 and 04 anywhere. If you only need a
plausible synthetic cohort and do not care about matching the real marginals,
run 03 and 04 with no profile present (the defaults are clinically shaped).

## Running the minimal merge on the real data

`condense_icon_to_minimal()` in script 01 reduces the full ICON
`og_cohort_cwt_2015_2022.rds` to the minimal Table A, deriving `cci_group` from
`rcs_ch_score` and `tumour_site_grp` from `cancer_subtype` if they are absent.
You can then re-run `og_cwt_merge()` on the condensed real cohort and confirm it
reproduces the audit tables, which is the cross-check that the minimal merge is
faithful to the full pipeline:

```r
source("01_og_minimal_merge.R")
full <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/og_cohort_cwt_2015_2022.rds")
minimal <- condense_icon_to_minimal(full)
# minimal already carries the merge outputs from the full pipeline; to re-run the
# merge itself you also need the raw CWT records as Table B.
```

## Outputs (written to `Data/synthetic/`)

- `og_ncras_treatment_synthetic.rds/.dta` - Table A, the registry + treatment
  cohort (left side of the merge), realistic on its own.
- `og_cwt_records_synthetic.rds/.dta` - Table B, the raw CWT records (right
  side), one row per recorded treatment event, dates as dd/mm/yyyy.
- `og_cohort_synthetic.rds/.dta` - the merged analysis cohort, with the CWT
  decision-to-treat node attached and the audit categories derived.
- `og_profile_for_synthetic.rds`, `og_minimal_spec.rds` - from script 02.

## What is reproduced, and one known simplification

The pathway model reproduces the audit numbers: stage-1-3 curative ~50%,
any-treatment ~75%, with SCC highest on definitive chemoRT, gastric highest on
surgery-only, matching the real cohort. The merge reproduces the modality /
pathway matching, the neoadjuvant primary-modality tie-break, and the
`dtt_valid` logic.

One deliberate simplification: the generator places the chemo and RT arms of a
definitive-chemoRT course close together, so synthetic `dtt_valid` for that
pathway is ~100%, where the real cohort sits at ~84% because the two arms start
weeks apart. The synthetic is cleaner than reality here; this is flagged rather
than forced, since a minimal example does not need to reproduce every timing
quirk. If you need that fidelity, widen the sact-to-rt gap for the definitive
chemoRT branch in script 03.
