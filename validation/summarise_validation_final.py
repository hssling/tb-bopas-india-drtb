"""Summarise the completed validation of the direction-of-effect classifier.

Reads validation/validation_codes_final.csv (150 record-intervention pairs: the classifier's code,
two coders' independent codes, and the final consensus code) and writes
outputs/validation_final.json: inter-coder agreement and Cohen's kappa (independent codes), the share
of pairs judged not to evaluate the intervention (consensus), and the classifier's agreement with the
consensus. Standard library only. Run: python -I validation/summarise_validation_final.py
"""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "validation" / "validation_codes_final.csv"
OUT = ROOT / "outputs" / "validation_final.json"
CATS = ["favours intervention", "no effect or favours control", "unclear",
        "not an evaluation of this intervention"]
FAV, NOT_EVAL = CATS[0], CATS[3]


def kappa(x: list[str], y: list[str]) -> tuple[float, float]:
    n = len(x)
    po = sum(a == b for a, b in zip(x, y)) / n
    cx, cy = Counter(x), Counter(y)
    pe = sum(cx[c] * cy[c] for c in CATS) / (n * n)
    return po, (po - pe) / (1 - pe)


def main() -> None:
    rows = list(csv.DictReader(SRC.open(encoding="utf-8", newline="")))
    a = [r["coder_A_independent"] for r in rows]
    b = [r["coder_B_independent"] for r in rows]
    cons = [r["consensus_code"] for r in rows]
    mach = [r["machine_code"] for r in rows]
    po_ab, k_ab = kappa(a, b)
    ev = [r for r in rows if r["consensus_code"] != NOT_EVAL]
    tp = sum(r["machine_code"] == FAV and r["consensus_code"] == FAV for r in ev)
    fp = sum(r["machine_code"] == FAV and r["consensus_code"] != FAV for r in ev)
    fn = sum(r["machine_code"] != FAV and r["consensus_code"] == FAV for r in ev)
    tn = len(ev) - tp - fp - fn
    by_node: dict[str, list[int]] = {}
    for r in rows:
        v = by_node.setdefault(r["intervention_node"], [0, 0])
        v[0] += r["consensus_code"] == NOT_EVAL
        v[1] += 1
    res = {
        "pairs": len(rows),
        "title_only": sum(r["title_only"] == "yes" for r in rows),
        "independent_agreement": f"{round(po_ab * len(rows))}/{len(rows)}",
        "independent_kappa": round(k_ab, 2),
        "consensus_counts": dict(Counter(cons)),
        "not_an_evaluation": f"{len(rows) - len(ev)}/{len(rows)}",
        "classifier_exact_agreement": f"{sum(m == c for m, c in zip(mach, cons))}/{len(rows)}",
        "classifier_kappa": round(kappa(mach, cons)[1], 2),
        "classifier_on_evaluations": {"n": len(ev), "correct": tp + tn, "TP": tp, "FP": fp,
                                      "FN": fn, "TN": tn, "PPV": round(tp / (tp + fp), 2),
                                      "sensitivity": round(tp / (tp + fn), 2),
                                      "specificity": round(tn / (tn + fp), 2)},
        "classifier_favours_calls_correct": f"{tp}/{mach.count(FAV)}",
        "not_an_evaluation_by_node": {k: f"{v[0]}/{v[1]}" for k, v in
                                      sorted(by_node.items(), key=lambda kv: int(kv[0]))},
    }
    OUT.write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
