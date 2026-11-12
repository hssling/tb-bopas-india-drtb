# TB-BOPAS — a Bayesian priority-action score applied to India's drug-resistant tuberculosis surveillance, 2017–2024

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

## Citation

See `CITATION.cff`. Once this repository has a Zenodo DOI, cite the DOI.

## Licence

CC BY 4.0. See `LICENSE`, including a note on the terms attaching to the underlying
WHO and bibliographic sources.
