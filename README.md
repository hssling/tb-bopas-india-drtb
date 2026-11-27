# TB-BOPAS — India's drug-resistant tuberculosis surveillance, 2017–2024, and an evidence map of tuberculosis interventions

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24167197.svg)](https://doi.org/10.5281/zenodo.24167197)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

Analysis code, input data, outputs and verification scripts for the manuscript
*"Applying a Bayesian priority-action score to India's drug-resistant tuberculosis
surveillance, 2017–2024: no worsening axis, and the hazard of uncorrected reporting
discontinuities."*

## What the study found

Five national drug-resistant TB (DR-TB) axes were reconstructed from World Health
Organization Global TB Programme country data for India.

| Axis | Slope | Classification |
|---|---:|---|
| MDR/RR-TB treatment success (2014–2022 cohorts) | +4.23 pp/yr (+3.01 to +5.49) | improving |
| XDR-TB treatment success (2014–2022 cohorts) | +6.38 pp/yr (+5.47 to +7.90) | improving |
| RR/MDR-TB notifications per 100 000 (2017–2024) | +0.09 pp/yr (−0.11 to +0.30) | flat/uncertain |
| Fluoroquinolone-resistant RR-TB, pre-XDR (2019–2024) | −0.77 pp/yr (−5.82 to +4.50) | flat/uncertain |
| Rifampicin resistance, new cases — within era (2021–2024) | −0.36 pp/yr (−0.56 to −0.22) | improving |

**No axis is significantly worsening**, and the score returns no urgent priority
(11 non-zero cells of 60, maximum 0.090).

The methodological result is the transferable one. Rifampicin resistance among new
cases rises at +0.25 pp/yr when fitted naively across 2017–2024, but that is a
**+1.58-point level shift at 2021** — the WHO definitional revision plus a 3.7-fold
growth in the testing denominator — with declining slopes inside each era. Fitting
across that discontinuity is by itself enough to make the score nominate a confident
top-ranked national priority, 2.7× the corrected maximum, on an axis that is
improving.

## Reproduce

Python 3.11. The analysis has **no third-party numerical dependency** — the
incomplete beta function is evaluated by continued fraction with bisection for the
credible bound. Only the figures need `matplotlib`.

```bash
python code/build_tb_bopas_ijtb.py      # -> outputs/   (trajectory fits, scores, 6 sensitivity scenarios)
python code/build_figures_ijtb.py       # -> figures/   (needs matplotlib; writes 600 dpi PNG + TIFF)
python code/build_supplement_ijtb.py    # -> outputs/supplementary_tables.md
```

Verified from a clean checkout: `build_tb_bopas_ijtb.py` reproduces every number
reported in the manuscript, including all five slopes, the +1.58 level shift, the
0.090 maximum score and the count of 11 non-zero cells.

## Verify

```bash
python verification/test_beta_bound.py       # hand-rolled special function vs SciPy
python verification/recompute_sensitivity.py # sensitivity probes on the pipeline
```

`test_beta_bound.py` checks the continued-fraction incomplete beta against
`scipy.special.betainc` over 156 points, the credible bound against
`scipy.stats.beta.ppf` over 20 cases, boundary and monotonicity properties, and
**every one of the twelve favourability bounds published in the manuscript**. Worst
deviation ~1e-9. SciPy is the test oracle only; the analysis itself does not use it.

## Contents

| Path | What it is |
|---|---|
| `data/ntep_trajectories.csv` | The five surveillance axes. One row per datapoint, each carrying its numerator, denominator, WHO field names and source dataset, so every value is traceable to a published field rather than to a summary table. |
| `data/tb_interventions_unified.csv` | The 3,589-record open-access literature harvest (20 May 2026). Note the column named `abstract_present` holds abstract *text*, not a flag, and is empty for 1,586 records. |
| `code/build_tb_bopas_ijtb.py` | Trajectory fitting with bootstrap and *t*-based intervals, segmentation at documented reporting breaks, corpus coding, scoring, six sensitivity scenarios. |
| `code/harvest_tb_interventions.py` | The literature harvest. |
| `outputs/` | Every result the manuscript cites: slopes, node counts, the full 60-cell score grid, sensitivity scenarios, corpus diagnostics. |
| `figures/` | Figures 1–3. TIFF print masters are git-ignored and regenerate at 600 dpi. |
| `verification/` | Pre-submission checks, including the SciPy reference test. |
| `validation/` | A stratified 150-record worksheet and analysis script for validating the automated direction-of-effect classifier against independent human coding. **Not yet completed** — classifier accuracy is unmeasured, and the manuscript states this as a limitation. |

## Two things to read before trusting the literature stream

Both are stated in the manuscript and repeated here because they bound what the
evidence stream can support.

1. **Only 1,831 of the 3,589 harvested records (51.0%) match an intervention node**
   and contribute to the analysis; 1,758 match none. 44% of records carry no abstract
   text, so their coding rests on the title alone.
2. **The direction-of-effect classifier was corrected** relative to the parent
   framework: contrary cues now take precedence over favour cues (28.1% of records
   scored favourable carried both), and prevent-family nomenclature was removed from
   the evaluative dictionary (it was the sole trigger for 13.6%, and it names two of
   the twelve nodes). The uncorrected classifier is retained as a sensitivity scenario.

## Data sources

- WHO Global TB Programme country data for India, retrieved 8 September 2026 from
  <https://extranet.who.int/tme/> (notifications, outcomes, drug-resistance
  surveillance, burden estimates). Re-checked 12 November 2026: India's latest data
  year was unchanged at 2024.
- Central TB Division, Ministry of Health and Family Welfare, Government of India.
  *India TB Report 2024*. Consulted as the national primary source and cross-checked
  against the WHO fields; total MDR/RR-TB notified for 2023 agrees to within 1.0%
  (63,929 against 64,595).
- PubMed E-utilities and Europe PMC (`OPEN_ACCESS:Y`) for the intervention literature.

## Version 1.3.0: evidence map of the intervention literature

**What it is.** Record-level coding of the 3589 record–intervention pairs in the literature harvest
(12 interventions; 3572 records after merging 17 duplicate titles). Each record is coded,
for the intervention it was retrieved for, as an evaluation (with direction, design and setting), as relevant
without evaluating it (tier 2), or as off-topic, and carries the verbatim span of title or abstract that supports
the code. A 60-record random sample was coded independently by two human coders. Fixed, ordered rules assign each
intervention an evidence-readiness tier, with sensitivity analyses.

**Headline counts** (`outputs/full_bopas/results.json`): 746 evaluations, 1295 tier-2 and
1531 off-topic records; human coders agreed at kappa 0.81; final machine codes against
the coders' final codes kappa 0.64.

| Path | What it is |
|---|---|
| `data/full_bopas/pairs_no_abstract.jsonl` | One line per record–intervention pair: identifiers, intervention, year, title, `title_only`, `abstract_source`, and the SHA-256 and length of the abstract text that was coded. The abstract text itself is not redistributed. |
| `data/full_bopas/codes/*.codes.jsonl` | Every coding layer, unchanged. Applied in the order `batch_*` → `recode_*` → `zscope_*` → `ztier_*` → `zz*`; a later layer overrides an earlier one for the same `pair_id`. |
| `data/full_bopas/codebook.md` | The instructions given to the model, including the scope of each intervention and the tier-2 rule. |
| `data/full_bopas/spot_check_codes.csv` | The 60-record human check: item, pair_id, each coder's independent code, the agreed code, the final code after adjudication, the recorded reason, and the machine code at sampling. No names. |
| `code/full_bopas/paper/analysis.py` | Evidence counts, Wilson intervals, lower credible bounds, readiness tiers and sensitivity analyses, agreement statistics. Writes `outputs/full_bopas/results.json`, `readiness.csv`, `agreement.json`. |
| `code/full_bopas/paper/build_figures.py` | Figures from `results.json` into `outputs/full_bopas/figures/` (PNG). |
| `outputs/full_bopas/evidence_by_intervention.csv` | Per-intervention summary from the coding pipeline (`code/full_bopas/analyse_codes.py`), used as a cross-check. |

**Reproduce.**

```
python -I code/full_bopas/paper/analysis.py
python -I code/full_bopas/paper/build_figures.py
```

Two fields need the abstract text and are written as `null` without it: `quotes_verbatim` and `off_topic_terms`.
To compute them, rebuild the abstracts with `code/full_bopas/prepare_batches.py` (from
`data/tb_interventions_unified.csv`), `fetch_abstracts.py` and `merge_fetched.py` (NCBI E-utilities), then run
`analysis.py --pairs-with-abstracts <pairs_all.jsonl>`; every abstract is checked against the stored SHA-256 first.
Run from the repository with these files, `analysis.py` reproduced the committed `results.json` field for field,
except those two fields; with the abstracts supplied, it reproduced every field.

**AI coding disclosure.** The records were coded by a large language model, Claude Opus 5.5 (Anthropic), in
November 2026, in batches, following `codebook.md`. Part of the work for the bovine TB, milk and wildlife/occupational interventions
(840 final codes) was done by a Claude Sonnet model with the identical codebook, because the provider's automated
safety filter halted several Opus runs. Codebook wording changed between rounds and earlier codes were
re-examined after each change; every round is kept as a separate layer. Each code carries a verbatim supporting quote, matched to the coded
text by script (3572/3572 matched). Two human coders checked a random sample; see
`outputs/full_bopas/agreement.json`. The analysis code, figures and this documentation were written with Claude
under the authors' direction.

**Data provenance.** Pairs come from `data/tb_interventions_unified.csv` (searched 20 May 2026). Records that had no
abstract were given the PubMed abstract fetched by PMID in November 2026, attached only when its title agreed with
the record's title; suspected mismatched abstracts were removed. Quotes in the code files are short verbatim spans of
titles and abstracts. The repository holds code, data and checks only.

## Version 1.1.0: trend component and classifier validation

- `code/trend_component.py` recomputes the surveillance-trend term of the score (without the literature and feasibility terms) under five analytic choices and writes `outputs/trend_component.json`. Run: `python -I code/trend_component.py outputs/trend_component.json`.
- `validation/validation_codes_final.csv` holds the completed validation of the direction-of-effect classifier: 150 record-intervention pairs with the classifier code, two coders' independent codes and the consensus code (identifiers and codes only; titles and abstracts are in the corpus). `python -I validation/summarise_validation_final.py` writes `outputs/validation_final.json`.
- Result: inter-coder kappa 0.90; 95 of 150 pairs did not evaluate the assigned intervention; on the 55 genuine evaluations the classifier's direction matched the consensus in 29. The literature-evidence term is not supported by the automated classifier.

## Version 1.1.1: zoonotic tuberculosis reporting

- `data/who_zoonotic_tb_2026-11-26.csv`: WHO Global Tuberculosis Programme zoonotic TB dataset (https://extranet.who.int/tme/generateCSV.asp?ds=zoonotic, retrieved 26 November 2026; years 2018 and 2023).
- `python -I code/zoonotic_reporting.py` writes `outputs/zoonotic_reporting.json`: India reported no speciation results in 2018 and no zoonotic TB diagnoses in 2023; 34 of 204 responding countries reported zoonotic diagnoses in 2023.

## Version 1.2.3: animal-health reports

- `code/animal_sector_search.py` searches the Department of Animal Husbandry and Dairying annual reports 2023-24, 2024-25 and 2025-26 and the ICAR-NIVEDI annual reports 2024 and 2025 (URLs and SHA-256 in the output) for bovine tuberculosis; brucellosis and foot-and-mouth disease are counted as controls. Output: `outputs/animal_sector_search.json`. The reports are not redistributed.

## Version 1.2.2: per-intervention validation shares

- `outputs/revision_analyses.json` now also reports, for each intervention, the share of sampled pairs that are not evaluations under the authors' consensus (primary reference) and under the final codes after re-examination.

## Version 1.2.1: new versus previously treated patients

- `data/who_dr_surveillance_2026-11-26.csv`: WHO drug-resistance surveillance dataset (https://extranet.who.int/tme/generateCSV.asp?ds=dr_surveillance, retrieved 26 November 2026). `code/revision_analyses.py` now tests whether the 2021 rise in rifampicin resistance among new patients reflects a shift from previously treated patients: it does not; the rise appears in both groups and among all tested patients.

## Version 1.2.0: analyses added after peer review

- `python -I code/revision_analyses.py` (needs numpy and scipy) writes `outputs/revision_analyses.json`: segmented models with t intervals at each series' reporting break; rifampicin-testing coverage and the share of notified RR-TB classed as new; treatment-cohort coverage and outcome composition; validation metrics against the final codes, the pre-check consensus and each coder, with Wilson intervals; bedaquiline intervals; and the distribution of WHO zoonotic answers.
- `data/who_2026-08-25/`: the WHO country files used (notifications, outcomes, laboratory services, data dictionary; release of 25 August 2026).
- `validation/validation_codes_final.csv` now also holds the consensus before the AI-assisted check and the AI recoding; `validation/ai_blind_recode_prompt.md` gives the model, inputs and instructions.

## Citation

Cite the archive by DOI.

- **Version 1.3.0:** DOI assigned on release
- **Version 1.2.3:** 10.5281/zenodo.24566970
- **Version 1.3.0:** 10.5281/zenodo.24594999
- **Version 1.2.2:** 10.5281/zenodo.24566563
- **Version 1.2.1:** 10.5281/zenodo.24566120
- **Version 1.2.0:** 10.5281/zenodo.24566030
- **Version 1.1.1:** 10.5281/zenodo.24565313
- **Version 1.1.0:** 10.5281/zenodo.24564979
- **Version 1.0.1:** 10.5281/zenodo.24169718
- **Version 1.0.0:** 10.5281/zenodo.24167198
- **All versions (concept DOI, resolves to the latest):** 10.5281/zenodo.24167197

Machine-readable metadata is in `CITATION.cff` and `.zenodo.json`.

## Licence

CC BY 4.0. See `LICENSE`, including a note on the terms attaching to the underlying
WHO and bibliographic sources.
