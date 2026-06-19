# OG cancer treatment pathways - Stata workflow


## The two steps

There are two scripts, run in order:

**`02_og_derive_pathway.do` - work out each patient's treatment pathway.**
It looks at the dates of the treatments a patient received (endoscopy, surgery,
chemotherapy, radiotherapy, and so on) and decides which *pathway* they followed
- for example "surgery only", "chemotherapy before surgery", or "no treatment
recorded". It also works out the single date that counts as the start of
definitive treatment (the "clock-stop"), and which hospital trust delivered it.

**`03_og_cwt_merge.do` - add the waiting-times record**
The NHS keeps a separate Cancer Waiting Times (CWT) record of treatment events.
This step matches the right CWT record to each patient, works out the waiting
times between the key milestones (diagnosis, decision to treat, first treatment),
and groups everyone into the standard national audit categories. The result is
the final patient table.

Run `02` first, then `03`. Each reads from and writes to the data folder named at
the very top of the file (one line, easy to change if the data moves).


## What goes in: the two raw datasets

### 1. The cohort - one row per patient

This is the registry record of who the patients are and the dates of the
treatments they received. The pathway logic is built entirely from these.

| Variable            | What it is                                            |
|---------------------|-------------------------------------------------------|
| `pseudo_patientid`  | Anonymous patient identifier (links the two datasets) |
| `diagmdy`           | Date of cancer diagnosis                              |
| `stage_clean`       | Cancer stage (1-4)                                    |
| `endoscopy_date`    | Date of the diagnostic endoscopy                     |
| `emresd_date`       | Date of endoscopic resection (EMR/ESD), if any       |
| `surgery_date`      | Date of major surgery, if any                        |
| `sact_date`         | Date chemotherapy started, if any                    |
| `rt_date`           | Date radiotherapy started, if any                    |
| `curative_surgery`  | Was the surgery done with curative intent? (yes/no)  |
| `rt_curative`       | Was the radiotherapy curative rather than palliative?|
| `surgery_provider`  | Code for the trust that did the surgery              |
| `rt_provider`       | Code for the trust that delivered the radiotherapy   |

A few extra fields (chemotherapy source, tumour subtype, age, sex, deprivation,
ethnicity, survival) can be included for later analysis but are **not required**
to work out the pathway.

### 2. The CWT records - one row per treatment event

This is the NHS Cancer Waiting Times feed. A patient can appear on several rows
(one per recorded event), and this step picks the one that matches their pathway.

| Variable             | What it is                                           |
|----------------------|------------------------------------------------------|
| `pseudo_patientid`   | Anonymous patient identifier (links to the cohort)   |
| `modality`           | Treatment type code (e.g. 01 = surgery, 02 = chemo)  |
| `site_icd10`         | Cancer site code (should be C15x / C16x for OG)      |
| `treat_period_start` | Decision-to-treat date                               |
| `treat_start`        | Date treatment actually started                      |
| `mdt_date`           | Date discussed at the multidisciplinary team meeting |

The date fields arrive as text (e.g. `04/03/2018`) and the script converts them
to proper dates before using them.


## What comes out: the final patient table

The final table (`og_cohort_synthetic_stata.dta`) has one row per patient with,
among others:

| Variable                     | What it tells you                                  |
|------------------------------|----------------------------------------------------|
| `tx_pathway`                 | The treatment pathway the patient followed         |
| `first_tx_date`              | The clock-stop date (start of definitive treatment)|
| `tx_trust`                   | The trust that delivered the clock-stop treatment  |
| `tx_modality_audit`          | The national audit treatment group                 |
| `tx_intent_audit`            | Curative / non-curative / no treatment             |
| `received_curative_tx_audit` | Did the patient get curative treatment? (yes/no)   |
| `received_any_tx`            | Did the patient get any treatment? (yes/no)        |
| `wt_*` (several)             | Waiting times between the key milestones, in days   |
| `dtt_valid`                  | Is the decision-to-treat date usable for a wait?    |


## The treatment pathways

Every patient is placed in exactly one of these:

- EMR/ESD only
- EMR/ESD then surgery
- Surgery only
- Surgery + neoadjuvant chemo (chemo before surgery)
- Surgery + neoadjuvant chemoRT (chemo and radiotherapy before surgery)
- Surgery + neoadjuvant RT (radiotherapy before surgery)
- Surgery + adjuvant chemo (chemo after surgery)
- Surgery + other
- Definitive chemoRT (curative chemo + radiotherapy, no surgery)
- Curative RT only
- Palliative chemo + RT
- SACT only (systemic anti-cancer therapy / chemo only)
- Palliative RT only
- No treatment recorded


## A note on the surgery codes

Around late 2020 the NHS changed the CWT code for surgery: the old code `01` was
replaced by new codes `23`/`24`. For most of 2020 both were in use. The script
handles this with a short *transition window* (Jan 2020 - Jun 2021) during which
all three codes count as surgery; before it only `01` counts, after it only
`23`/`24`. This stops genuine 2020 operations being missed because of the coding
change. The window dates are set near the top of `03` if they ever need adjusting.


## Running it

1. Put the two raw `.dta` files in the data folder (the path at the top of each
   script - currently `Test_data`).
2. Run `02_og_derive_pathway.do`. It prints the pathway mix and saves the derived
   cohort.
3. Run `03_og_cwt_merge.do`. It prints the CWT coverage and the headline audit
   figures (% curative, % any treatment) and saves the final cohort.

If a file is missing the script stops and tells you which one. The scripts do not
change the raw data - they only read it and write new output files.
