"""Analyses added after peer review: reporting breaks in every series, testing and cohort coverage,
outcome composition, urgency-term scenarios, validation sensitivity analyses and zoonotic answers.

Inputs (all in this repository):
  data/ntep_trajectories.csv                      indicator values with numerators and denominators
  data/who_2026-08-25/TB_notifications_2026-08-25.csv, TB_outcomes_2026-08-25.csv,
  TB_laboratory_diagnostic_services_2026-08-25.csv WHO Global TB Programme country data (release 25 Aug 2026)
  outputs/ijtb_slopes.csv, code/build_tb_bopas_ijtb.py  (bootstrap half-widths and the causal-link matrix)
  validation/validation_codes_final.csv           classifier, coder, pre-check consensus, AI recode, final codes
  data/who_zoonotic_tb_2026-11-26.csv             WHO zoonotic TB dataset
Output: outputs/revision_analyses.json.   Run: python -I code/revision_analyses.py
Requires numpy and scipy (the original pipeline is standard-library only).
"""
from __future__ import annotations

import ast
import csv
import json
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
WHO = ROOT / "data" / "who_2026-08-25"
OUT = ROOT / "outputs" / "revision_analyses.json"
FLOOR = 0.5


def series():
    s: dict[str, dict[int, dict]] = {}
    for r in csv.DictReader((ROOT / "data" / "ntep_trajectories.csv").open(encoding="utf-8")):
        s.setdefault(r["series"], {})[int(r["year"])] = {
            "value": float(r["value"]), "num": float(r["numerator"]), "den": float(r["denominator"])}
    return s


def india(fn: str) -> dict[int, dict]:
    return {int(r["year"]): r for r in csv.DictReader((WHO / fn).open(encoding="utf-8", errors="replace"))
            if r["iso3"] == "IND"}


def num(x):
    return float(x) if x not in (None, "", "NA") else None


def ols(x, y):
    """Slope with 95% t interval (df = n - 2)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    res = stats.linregress(x, y)
    df = len(x) - 2
    h = stats.t.ppf(0.975, df) * res.stderr if df > 0 else float("nan")
    return {"n": len(x), "slope": res.slope, "lo": res.slope - h, "hi": res.slope + h}


def segmented(years, values, brk):
    """y = b0 + b1*(t-brk) + b2*post + b3*(t-brk)*post; 95% t intervals, df = n - 4."""
    t = np.asarray(years, float) - brk
    post = (np.asarray(years) >= brk).astype(float)
    X = np.column_stack([np.ones_like(t), t, post, t * post])
    y = np.asarray(values, float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    df = len(y) - 4
    s2 = resid @ resid / df
    cov = s2 * np.linalg.inv(X.T @ X)
    q = stats.t.ppf(0.975, df)

    def ci(c):
        est = float(c @ beta)
        se = math.sqrt(float(c @ cov @ c))
        return {"est": est, "lo": est - q * se, "hi": est + q * se}
    return {"break": brk, "n": len(y), "df": df,
            "pre_slope": ci(np.array([0, 1, 0, 0])),
            "level_change_at_break": ci(np.array([0, 0, 1, 0])),
            "post_slope": ci(np.array([0, 1, 0, 1])),
            "slope_change": ci(np.array([0, 0, 0, 1]))}


def wilson(k, n):
    if n == 0:
        return None
    z = 1.959964
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return {"k": k, "n": n, "p": p, "lo": c - h, "hi": c + h}


CATS = ["favours intervention", "no effect or favours control", "unclear",
        "not an evaluation of this intervention"]
FAV, NE = CATS[0], CATS[3]


def kappa(x, y, cats=CATS):
    n = len(x)
    po = sum(a == b for a, b in zip(x, y)) / n
    cx, cy = Counter(x), Counter(y)
    pe = sum(cx[c] * cy[c] for c in cats) / (n * n)
    k = (po - pe) / (1 - pe) if pe < 1 else float("nan")
    se = math.sqrt(po * (1 - po) / (n * (1 - pe) ** 2)) if pe < 1 else float("nan")
    return {"agree": round(po * n), "n": n, "kappa": k, "lo": k - 1.96 * se, "hi": k + 1.96 * se}


def classifier_metrics(rows, ref):
    ev = [r for r in rows if r[ref] != NE]
    tp = sum(r["machine_code"] == FAV and r[ref] == FAV for r in ev)
    fp = sum(r["machine_code"] == FAV and r[ref] != FAV for r in ev)
    fn = sum(r["machine_code"] != FAV and r[ref] == FAV for r in ev)
    tn = len(ev) - tp - fp - fn
    bin_m = ["F" if r["machine_code"] == FAV else "N" for r in ev]
    bin_r = ["F" if r[ref] == FAV else "N" for r in ev]
    return {"reference": ref,
            "not_an_evaluation": wilson(len(rows) - len(ev), len(rows)),
            "assignment_precision": wilson(len(ev), len(rows)),
            "direction_accuracy": wilson(tp + tn, len(ev)),
            "ppv_favourable_grade": wilson(tp, tp + fp),
            "sensitivity_conditional": wilson(tp, tp + fn),
            "specificity_conditional": wilson(tn, tn + fp),
            "favourable_grades_genuine": wilson(tp, sum(r["machine_code"] == FAV for r in rows)),
            "kappa_2x2_on_evaluations": kappa(bin_m, bin_r, ["F", "N"]),
            "TP": tp, "FP": fp, "FN": fn, "TN": tn}


def main():
    S = series()
    notif = india("TB_notifications_2026-08-25.csv")
    outc = india("TB_outcomes_2026-08-25.csv")
    lab = india("TB_laboratory_diagnostic_services_2026-08-25.csv")
    res: dict = {}

    # 1. reporting breaks in every series (segmented OLS, t intervals)
    def yv(name, y0=None, y1=None):
        d = S[name]
        ys = [y for y in sorted(d) if (y0 is None or y >= y0) and (y1 is None or y <= y1)]
        return ys, [d[y]["value"] for y in ys]
    rr_y, rr_v = yv("RR_new_pct")
    res["rr_new"] = {
        "one_line": ols(rr_y, rr_v),
        "pre_2021": ols(*yv("RR_new_pct", None, 2020)),
        "post_2021": ols(*yv("RR_new_pct", 2021, None)),
        "segmented": {str(b): segmented(rr_y, rr_v, b) for b in (2019, 2020, 2021)},
        "absolute_resistant_new": {y: S["RR_new_pct"][y]["num"] for y in rr_y},
    }
    for name, brk in (("MDR_TSR_pct", 2020), ("XDR_TSR_pct", 2020), ("RRMDR_notif_100k", 2020)):
        ys, vs = yv(name)
        res[name] = {"one_line": ols(ys, vs), "pre": ols(*yv(name, None, brk - 1)),
                     "post": ols(*yv(name, brk, None)) if len([y for y in ys if y >= brk]) > 2 else None,
                     "segmented": segmented(ys, vs, brk)}
    res["preXDR_fqr_pct"] = {"one_line": ols(*yv("preXDR_fqr_pct"))}

    # 2. coverage
    cov = {}
    for y in range(2017, 2025):
        n_rr = S["RRMDR_notif_100k"][y]["num"]
        lab_new = num(notif[y].get("pulm_labconf_new"))
        cov[y] = {
            "rif_tested_new": S["RR_new_pct"][y]["den"],
            "pulm_labconf_new": lab_new,
            "rif_testing_coverage_new": S["RR_new_pct"][y]["den"] / lab_new if lab_new else None,
            "new_share_of_notified_rr": S["RR_new_pct"][y]["num"] / n_rr,
            "fq_tested_share_of_notified_rr": (S["preXDR_fqr_pct"][y]["den"] / n_rr
                                               if y in S["preXDR_fqr_pct"] else None),
            "bdq_tested_share_of_notified_rr": (S["BDQ_res_pct"][y]["den"] / n_rr
                                                if y in S["BDQ_res_pct"] else None),
            "labs_bdq_dst": num(lab[y].get("dst_bdq")) if y in lab else None,
        }
    for y in range(2017, 2023):
        o = outc[y]
        mdr, xdr = num(o["mdr_coh"]), num(o["xdr_coh"])
        n_rr = S["RRMDR_notif_100k"][y]["num"]
        cov[y].update({"mdr_cohort_share_of_notified": mdr / n_rr,
                       "mdr_plus_xdr_cohort_share_of_notified": (mdr + xdr) / n_rr})
    res["coverage"] = cov

    # 3. outcome composition by cohort year (percent of cohort)
    comp = {}
    for y in range(2014, 2023):
        o = outc[y]
        comp[y] = {}
        for g in ("mdr", "xdr"):
            c = num(o[f"{g}_coh"])
            comp[y][g] = {"cohort": c, **{k: 100 * num(o[f"{g}_{k}"]) / c
                                           for k in ("succ", "fail", "died", "lost", "neval")}}
    res["outcome_composition"] = comp

    # 4. urgency-term scenarios
    src = (ROOT / "code" / "build_tb_bopas_ijtb.py").read_text(encoding="utf-8")
    M = ast.literal_eval(re.search(r"^M_TB\s*=\s*(\{.*?^\})", src, re.S | re.M).group(1))
    slopes = list(csv.DictReader((ROOT / "outputs" / "ijtb_slopes.csv").open(encoding="utf-8")))
    res["slopes_csv_columns"] = list(slopes[0].keys())
    res["link_matrix_rows"] = len(M)

    # 5. validation sensitivity
    V = list(csv.DictReader((ROOT / "validation" / "validation_codes_final.csv").open(encoding="utf-8")))
    res["validation"] = {
        "coders_independent": kappa([r["coder_A_independent"] for r in V],
                                    [r["coder_B_independent"] for r in V]),
        "coder_A_vs_final": kappa([r["coder_A_independent"] for r in V], [r["consensus_code"] for r in V]),
        "coder_B_vs_final": kappa([r["coder_B_independent"] for r in V], [r["consensus_code"] for r in V]),
        "ai_vs_consensus_before_check": kappa([r["ai_blind_recode"] for r in V],
                                              [r["consensus_before_ai_check"] for r in V]),
        "ai_vs_final": kappa([r["ai_blind_recode"] for r in V], [r["consensus_code"] for r in V]),
        "pairs_changed_after_ai_check": sum(r["consensus_before_ai_check"] != r["consensus_code"] for r in V),
        "unanimous_pairs_changed": sum(r["coder_A_independent"] == r["coder_B_independent"]
                                       != r["consensus_code"] for r in V),
        "classifier_vs": {ref: classifier_metrics(V, ref) for ref in
                          ("consensus_code", "consensus_before_ai_check",
                           "coder_A_independent", "coder_B_independent")},
        "classifier_four_code_vs_final": kappa([r["machine_code"] for r in V], [r["consensus_code"] for r in V]),
        "classifier_possible_codes": sorted({r["machine_code"] for r in V}),
    }

    # 6. bedaquiline
    res["bedaquiline"] = {y: wilson(int(S["BDQ_res_pct"][y]["num"]), int(S["BDQ_res_pct"][y]["den"]))
                          for y in sorted(S["BDQ_res_pct"])}

    # 7. zoonotic answers
    Z = list(csv.DictReader((ROOT / "data" / "who_zoonotic_tb_2026-11-26.csv").open(encoding="utf-8")))
    z18 = [r for r in Z if r["year"] == "2018"]
    z23 = [r for r in Z if r["year"] == "2023"]
    res["zoonotic"] = {
        "2018_pulm_spec_available": dict(Counter(r["pulm_spec_available"] or "blank" for r in z18)),
        "2018_ep_spec_available": dict(Counter(r["ep_spec_available"] or "blank" for r in z18)),
        "2023_ident_zoonotic": dict(Counter(r["ident_zoonotic"] or "blank" for r in z23)),
    }

    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    print(json.dumps(res, indent=1, default=float)[:6000])


if __name__ == "__main__":
    main()
