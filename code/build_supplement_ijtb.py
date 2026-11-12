"""Generate the supplementary material from pipeline outputs.

Every table here is emitted from the deposited CSVs. Nothing is typed by hand:
in the predecessor submission the hand-written supplementary tables drifted
from the pipeline and carried stale values.
"""
from __future__ import annotations
import csv
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = ROOT / "outputs"
DEST = ROOT / "outputs" / "supplementary_tables.md"

NODES = {
    "1": ("Shorter all-oral regimen rollout (BPaL/BPaL-M)",
          "BPaL, BPaL-M and the all-oral combination regimens from the Nix-TB, ZeNix and TB-PRACTECAL programmes, including paediatric and pregnancy adaptations"),
    "2": ("Molecular diagnostic scale-up",
          "Xpert MTB/RIF, Xpert MTB/XDR, line probe assays, Truenat"),
    "3": ("Whole-genome sequencing for drug-susceptibility prediction",
          "Sequencing-based resistance prediction and mutation catalogues"),
    "4": ("Active case finding",
          "Door-to-door, contact investigation, high-risk-population screening"),
    "5": ("TB preventive therapy", "3HP, 1HP, 6H, levofloxacin regimens for DR-TB contacts"),
    "6": ("Airborne infection prevention and control",
          "Ventilation, ultraviolet germicidal irradiation, respirator programmes"),
    "7": ("Nutritional support", "Patient and household-contact nutritional support"),
    "8": ("Private-sector engagement",
          "Informant networks, mandatory notification, voucher schemes"),
    "9": ("Integrated TB-HIV-diabetes care", "Joint and comorbidity care models"),
    "10": ("Bovine tuberculosis control", "Test-and-cull, herd vaccination, dairy-herd control"),
    "11": ("Milk pasteurisation / dairy food chain",
           "Pasteurisation enforcement and dairy food-chain safety"),
    "12": ("Wildlife reservoir surveillance / occupational respiratory protection",
           "Wildlife reservoir surveillance and slaughterhouse occupational protection"),
}

TARGET_ORDER = ["RR_new_pct", "RRMDR_notif_100k", "preXDR_fqr_pct", "MDR_TSR_pct", "XDR_TSR_pct"]
TARGET_SHORT = {"RR_new_pct": "RR, new cases", "RRMDR_notif_100k": "RR/MDR notif.",
                "preXDR_fqr_pct": "pre-XDR", "MDR_TSR_pct": "MDR/RR success",
                "XDR_TSR_pct": "XDR success"}


def read(name):
    with (OUT / name).open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def main():
    traj = read("../source/ntep_trajectories.csv") if False else list(
        csv.DictReader((ROOT / "data" / "ntep_trajectories.csv").open(encoding="utf-8", newline="")))
    slopes = {r["series"]: r for r in read("ijtb_slopes.csv")}
    nodes = read("ijtb_node_counts.csv")
    grid = read("ijtb_bopas_long.csv")
    sens = read("ijtb_sensitivity.csv")
    diag = (OUT / "ijtb_diagnostics.txt").read_text(encoding="utf-8")

    L = []
    A = L.append
    A("# Supplementary Material")
    A("")
    A("**Manuscript:** Applying a Bayesian priority-action score to India's drug-resistant "
      "tuberculosis surveillance, 2017–2024: no worsening axis, and the hazard of uncorrected "
      "reporting discontinuities.")
    A("")
    A("*All tables in this document are generated directly from the deposited analysis outputs "
      "by `build_supplement_ijtb.py`. No value is transcribed by hand.*")
    A("")
    A("## Contents")
    A("")
    A("- Supplementary Table S1. Surveillance dataset with per-datapoint provenance.")
    A("- Supplementary Table S2. Trajectory fits, including segmented fits at the documented reporting break.")
    A("- Supplementary Table S3. The twelve pre-specified intervention nodes.")
    A("- Supplementary Table S4. Per-node corpus counts, favourability bounds, and coding diagnostics.")
    A("- Supplementary Table S5. Mechanistic-coupling matrix.")
    A("- Supplementary Table S6. Complete priority-score grid (all 60 cells).")
    A("- Supplementary Table S7. Sensitivity scenarios.")
    A("- Supplementary Material S8. Corpus coverage and classifier-bias diagnostics.")
    A("- Supplementary Table S9. Cross-check of the surveillance inputs against the national programme report.")
    A("")
    A("---")
    A("")

    # ---- S1 -------------------------------------------------------------
    A("## Supplementary Table S1. Surveillance dataset with per-datapoint provenance")
    A("")
    A("Every value used in the analysis, with the numerator, denominator and World Health "
      "Organization field names from which it was derived. Source datasets are the WHO Global TB "
      "Programme country data downloads (notifications, treatment outcomes, drug-resistance "
      "surveillance, burden estimates), retrieved 8 September 2026.")
    A("")
    A("| Axis | Year | Value | Numerator | Denominator | WHO field(s) | Dataset |")
    A("|---|:---:|---:|---:|---:|---|---|")
    for r in traj:
        A(f"| {r['label']} | {r['year']} | {float(r['value']):.3f} | {r['numerator'] or '—'} | "
          f"{r['denominator'] or '—'} | `{r['who_fields']}` | {r['who_dataset']} |")
    A("")
    notes = {}
    for r in traj:
        if r["note"]:
            notes.setdefault(r["note"], set()).add(r["label"])
    A("Notes:")
    for note, labels in notes.items():
        A(f"- *{'; '.join(sorted(labels))}:* {note}")
    A("")
    A("---")
    A("")

    # ---- S2 -------------------------------------------------------------
    A("## Supplementary Table S2. Trajectory fits, with segmented fits at the reporting break")
    A("")
    A("Slopes by ordinary least squares in percentage points per year. Bootstrap intervals are "
      "percentile intervals from 4 000 case resamples; OLS intervals are *t*-based on the "
      "regression standard error with *n*−2 degrees of freedom.")
    A("")
    A("| Axis | Years | *n* | Slope | Bootstrap 95% CI | OLS 95% CI | Classification | Basis of primary analysis |")
    A("|---|:---:|:---:|---:|---|---|---|---|")
    for key in TARGET_ORDER:
        r = slopes[key]
        A(f"| {r['label']} | {r['years']} | {r['n']} | {float(r['slope_full']):+.3f} | "
          f"{float(r['boot_lo']):+.3f} to {float(r['boot_hi']):+.3f} | "
          f"{float(r['ols_lo']):+.3f} to {float(r['ols_hi']):+.3f} | "
          f"{r['classification_full']} | {r['primary_basis']} |")
    A("")
    brk = slopes["RR_new_pct"]
    A(f"**Segmented fit, {brk['label']}** (reporting break at {brk['break_year']}):")
    A("")
    A("| Era | Slope | Bootstrap 95% CI | Classification |")
    A("|---|---:|---|---|")
    A(f"| Pre-break (2017–2020) | {float(brk['slope_pre']):+.3f} | "
      f"{float(brk['pre_boot_lo']):+.3f} to {float(brk['pre_boot_hi']):+.3f} | {brk['classification_pre']} |")
    A(f"| Post-break (2021–2024) | {float(brk['slope_post']):+.3f} | "
      f"{float(brk['post_boot_lo']):+.3f} to {float(brk['post_boot_hi']):+.3f} | {brk['classification_post']} |")
    A(f"| Naive full period (2017–2024) | {float(brk['slope_full']):+.3f} | "
      f"{float(brk['boot_lo']):+.3f} to {float(brk['boot_hi']):+.3f} | {brk['classification_full']} |")
    A("")
    A(f"Level shift between era means: **{float(brk['level_shift']):+.2f} percentage points**. "
      "The primary analysis uses the post-break slope.")
    A("")
    A("---")
    A("")

    # ---- S3 -------------------------------------------------------------
    A("## Supplementary Table S3. The twelve pre-specified intervention nodes")
    A("")
    A("| Node | Label | Scope |")
    A("|:---:|---|---|")
    for nid, (lbl, scope) in NODES.items():
        A(f"| {nid} | {lbl} | {scope} |")
    A("")
    A("---")
    A("")

    # ---- S4 -------------------------------------------------------------
    A("## Supplementary Table S4. Per-node corpus counts and favourability bounds")
    A("")
    A("*k* is the number of records coded as favouring the intervention and *n* the number coded "
      "directionally (favourable or contrary); *L* is the 5th percentile of the Beta(*k*+1, "
      "*n*−*k*+1) posterior. The corrected classifier gives contrary cues precedence over favour "
      "cues and excludes prevent-family nomenclature from the evaluative dictionary; the "
      "uncorrected classifier is the one used in the parent framework.")
    A("")
    A("| Node | Records matched | *k*/*n* corrected | *L* corrected | *k*/*n* uncorrected | *L* uncorrected |")
    A("|:---:|---:|:---:|---:|:---:|---:|")
    for r in nodes:
        A(f"| {r['node_id']} | {r['n_total']} | {r['k_corrected']}/{r['n_dir_corrected']} | "
          f"{r['L_corrected']} | {r['k_uncorrected']}/{r['n_dir_uncorrected']} | {r['L_uncorrected']} |")
    A("")
    A("Note that *n* (directionally coded) rather than the matched-record count is the denominator "
      "of *L*: records with no evaluative cue, or with an explicitly unclear cue, are excluded "
      "from the estimate but retained in the matched count.")
    A("")
    A("---")
    A("")

    # ---- S5 -------------------------------------------------------------
    A("## Supplementary Table S5. Mechanistic-coupling matrix")
    A("")
    A("Values are 0 (no plausible mechanism), 0.5 (partial or indirect) and 1.0 (direct). "
      "Re-specified for the revised target set on the parent framework's three-level scheme. "
      "Because only one axis retains a positive point slope, the ranked output is insensitive to "
      "this matrix in the present application (Supplementary Table S7, M-strict scenario).")
    A("")
    M = defaultdict(dict)
    for r in grid:
        M[r["intervention"]][r["series"]] = r["M"]
    A("| Intervention node | " + " | ".join(TARGET_SHORT[k] for k in TARGET_ORDER) + " |")
    A("|---|" + ":---:|" * len(TARGET_ORDER))
    for nid, (lbl, _) in NODES.items():
        match = next((k for k in M if k.startswith(lbl[:28])), None)
        if match is None:
            continue
        A(f"| {nid} — {lbl} | " + " | ".join(M[match][k] for k in TARGET_ORDER) + " |")
    A("")
    A("---")
    A("")

    # ---- S6 -------------------------------------------------------------
    A("## Supplementary Table S6. Complete priority-score grid")
    A("")
    A("All 60 candidate cells under the primary analysis, ordered by score. `floor_binding` "
      "records whether the interval half-width floor of 0.5 was active for that cell.")
    A("")
    A("| Score | Intervention node | Target axis | *L* | *T* | *M* | √-term *C* | Half-width used | Floor binding | Axis classification |")
    A("|---:|---|---|---:|---:|:---:|:---:|---:|:---:|---|")
    for r in grid:
        A(f"| {r['BOPAS']} | {r['intervention']} | {r['target']} | {r['L']} | {r['T']} | "
          f"{r['M']} | {r['C']} | {r['half_width_used']} | {r['floor_binding']} | {r['trajectory_class']} |")
    A("")
    A("---")
    A("")

    # ---- S7 -------------------------------------------------------------
    A("## Supplementary Table S7. Sensitivity scenarios")
    A("")
    A("| Scenario | Non-zero cells | Rank 1 | Rank 2 | Rank 3 |")
    A("|---|:---:|---|---|---|")
    for r in sens:
        A(f"| {r['scenario']} | {r['non_zero_cells']} | {r['rank1']} | {r['rank2']} | {r['rank3']} |")
    A("")
    A("---")
    A("")

    # ---- S8 -------------------------------------------------------------
    A("## Supplementary Material S8. Corpus coverage and classifier-bias diagnostics")
    A("")
    A("Emitted by the analysis script as `ijtb_diagnostics.txt`.")
    A("")
    A("```")
    A(diag.rstrip())
    A("```")
    A("")
    A("---")
    A("")

    # ---- S9 --------------------------------------------------------------
    # The only values in this document not emitted by the pipeline. Each is
    # transcribed from the cited national report with its source identified.
    A("## Supplementary Table S9. Cross-check against the national programme report")
    A("")
    A("The surveillance inputs are taken from the World Health Organization country data because "
      "those fields are machine-readable and carry explicit numerators and denominators. The "
      "India TB Report 2024 was consulted as the national primary source for the same period. "
      "The figures below are transcribed from that report with the section or table identified, "
      "and are the only values in this supplement not emitted by the analysis pipeline.")
    A("")
    A("**Agreement on total drug-resistant TB notified, 2023**")
    A("")
    A("| Source | MDR/RR-TB notified, 2023 |")
    A("|---|---:|")
    A("| India TB Report 2024, \"Performance in 2023\" | 63 929 |")
    A("| WHO country data (`conf_rr_nfqr` + `conf_rr_fqr`) | 64 595 |")
    A("| Difference | +666 (+1.0%) |")
    A("")
    A("**Disagreement on the pre-XDR proportion, 2023.** The proportion depends entirely on the "
      "denominator, and all three readings below are defensible for the same year:")
    A("")
    A("| Denominator | Numerator / denominator | Proportion | Source |")
    A("|---|---:|---:|---|")
    A("| All notified MDR/RR-TB | 11 749 / 63 929 | 18.4% | India TB Report 2024, \"Performance in 2023\" |")
    A("| Second-line culture DST conducted | 2 667 / 18 505 | 14.4% | India TB Report 2024, Table 4.4 |")
    A("| RR-TB with a fluoroquinolone result | 11 134 / 24 866 | 44.8% | WHO `rr_fqr` / `rr_dst_rlt_fq` — **used in this paper** |")
    A("")
    A("The spread is 14.4% to 44.8%, a 3.1-fold range. The axis is therefore labelled throughout "
      "as a proportion of RR-TB tested for fluoroquinolone susceptibility, and is never compared "
      "against a proportion computed on a different base.")
    A("")
    A("**Two XDR-TB counts for the same year.** The national report gives 114 XDR-TB cases in its "
      "programme narrative for 2023 and 354 in its second-line culture-DST table (Table 4.4), and "
      "describes a programme data-validation exercise initiated to investigate such "
      "discrepancies. This is one reason no XDR-TB prevalence axis is fitted here, alongside the "
      "January 2021 change in the XDR-TB case definition.")
    A("")

    DEST.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"wrote {DEST} ({len(L)} lines)")


if __name__ == "__main__":
    main()
