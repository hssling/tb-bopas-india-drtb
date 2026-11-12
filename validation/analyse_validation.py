"""Analyse the filled validation worksheet.

Run this only after both coders have filled their own column independently and
recorded a consensus label for the rows where they differ.

Reports:
  1. inter-coder agreement (raw and Cohen's kappa, 4 categories)
  2. classifier agreement with the human consensus (raw, kappa, confusion)
  3. the binary question that drives the score: among records the humans judge
     to be directional evaluations, how often does the classifier call
     "favours intervention" correctly (sensitivity, specificity, PPV)
  4. what the favourability bound L would be per node on human labels, and how
     far the published bound sits from it

Reads validation_worksheet.xlsx if openpyxl can open it, else the CSV.
"""
from __future__ import annotations

import csv
import importlib.util
import math
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "code" / "build_tb_bopas_ijtb.py"
spec = importlib.util.spec_from_file_location("tb_ijtb", SRC)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

FAVOURS = "favours intervention"
AGAINST = "no effect or favours control"
UNCLEAR = "unclear"
NOT_EVAL = "not an evaluation of this intervention"
NO_CUE = "no evaluative cue"


def load_worksheet() -> list[dict]:
    xlsx = HERE / "validation_worksheet.xlsx"
    if xlsx.exists():
        try:
            from openpyxl import load_workbook
            wb = load_workbook(xlsx, data_only=True)
            sh = wb["worksheet"]
            header_row = None
            for i, row in enumerate(sh.iter_rows(values_only=True), start=1):
                if row and row[0] == "item":
                    header_row = i
                    header = [str(c) if c is not None else "" for c in row]
                    break
            if header_row:
                out = []
                for row in sh.iter_rows(min_row=header_row + 1, values_only=True):
                    if row[0] is None:
                        continue
                    out.append({h: ("" if v is None else str(v))
                                for h, v in zip(header, row)})
                return out
        except Exception as exc:            # fall back to the CSV
            print(f"(could not read the xlsx: {exc}; using the CSV)")
    with (HERE / "validation_worksheet.csv").open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def kappa(pairs: list[tuple[str, str]]) -> tuple[float, float]:
    """Cohen's kappa and raw agreement for a list of (rater1, rater2) labels."""
    n = len(pairs)
    if n == 0:
        return float("nan"), float("nan")
    po = sum(1 for a, b in pairs if a == b) / n
    c1, c2 = Counter(a for a, _ in pairs), Counter(b for _, b in pairs)
    pe = sum((c1[k] / n) * (c2[k] / n) for k in set(c1) | set(c2))
    k = (po - pe) / (1 - pe) if pe < 1 else float("nan")
    return k, po


def interpret(k: float) -> str:
    if k != k:
        return "undefined"
    for lo, label in ((0.81, "almost perfect"), (0.61, "substantial"),
                      (0.41, "moderate"), (0.21, "fair"), (0.0, "slight")):
        if k >= lo:
            return label
    return "poor (worse than chance)"


def main() -> None:
    rows = load_worksheet()
    key = {r["item"]: r for r in
           csv.DictReader((HERE / "validation_key.csv").open(encoding="utf-8", newline=""))}

    filled_a = [r for r in rows if r.get("coder_A_code", "").strip()]
    filled_b = [r for r in rows if r.get("coder_B_code", "").strip()]
    filled_c = [r for r in rows if r.get("consensus_code", "").strip()]
    print(f"worksheet rows: {len(rows)}")
    print(f"  coder A filled : {len(filled_a)}")
    print(f"  coder B filled : {len(filled_b)}")
    print(f"  consensus filled: {len(filled_c)}")
    if not filled_a or not filled_b:
        print("\nNothing to analyse yet. Both coders need to fill their own column in")
        print("validation_worksheet.xlsx (independently), then record a consensus label")
        print("for the rows where they differ, then re-run this script.")
        return

    # ---- 1. inter-coder agreement ---------------------------------------
    both = [(r["coder_A_code"].strip(), r["coder_B_code"].strip()) for r in rows
            if r.get("coder_A_code", "").strip() and r.get("coder_B_code", "").strip()]
    k, po = kappa(both)
    print(f"\n1. INTER-CODER AGREEMENT  (n = {len(both)})")
    print(f"   raw agreement  {po:.1%}")
    print(f"   Cohen's kappa  {k:.3f}  ({interpret(k)})")
    disagreements = [(a, b) for a, b in both if a != b]
    if disagreements:
        print(f"   {len(disagreements)} disagreement(s), most common patterns:")
        for (a, b), n in Counter(disagreements).most_common(5):
            print(f"      A={a[:34]:<36} B={b[:34]:<36} x{n}")

    # ---- 2. classifier vs human consensus -------------------------------
    def human(r):
        return (r.get("consensus_code", "").strip()
                or r.get("coder_A_code", "").strip())

    paired = [(key[r["item"]]["machine_code"], human(r)) for r in rows
              if human(r) and r["item"] in key]
    if not paired:
        print("\nNo rows carry a human label yet.")
        return
    km, pom = kappa([(mc, hc) for mc, hc in paired])
    print(f"\n2. CLASSIFIER vs HUMAN CONSENSUS  (n = {len(paired)})")
    print(f"   raw agreement  {pom:.1%}")
    print(f"   Cohen's kappa  {km:.3f}  ({interpret(km)})")
    conf = Counter(paired)
    labels = sorted({l for pair in paired for l in pair})
    w = max(len(l) for l in labels) + 2
    print(f"\n   confusion matrix (rows = classifier, columns = human)")
    print("   " + " " * w + "".join(f"{l[:14]:>16}" for l in labels))
    for mc in labels:
        line = f"   {mc[:w-2]:<{w}}"
        for hc in labels:
            line += f"{conf.get((mc, hc), 0):>16}"
        print(line)

    # ---- 3. the binary decision the score depends on --------------------
    directional = [(mc, hc) for mc, hc in paired if hc in (FAVOURS, AGAINST)]
    tp = sum(1 for mc, hc in directional if mc == FAVOURS and hc == FAVOURS)
    fp = sum(1 for mc, hc in directional if mc == FAVOURS and hc == AGAINST)
    fn = sum(1 for mc, hc in directional if mc != FAVOURS and hc == FAVOURS)
    tn = sum(1 for mc, hc in directional if mc != FAVOURS and hc == AGAINST)
    print(f"\n3. BINARY 'FAVOURS INTERVENTION' DECISION  (n = {len(directional)} "
          f"records the humans judged directional)")
    if tp + fn:
        print(f"   sensitivity {tp/(tp+fn):.1%}   ({tp}/{tp+fn} human-favourable records caught)")
    if tn + fp:
        print(f"   specificity {tn/(tn+fp):.1%}   ({tn}/{tn+fp} human-non-favourable records rejected)")
    if tp + fp:
        print(f"   PPV         {tp/(tp+fp):.1%}   ({tp}/{tp+fp} machine-favourable calls confirmed)")
    if fp > tp * 0.25 and tp:
        print("   NOTE: false positives remain material; the favourability bounds are")
        print("         still likely to be optimistic even after the correction.")

    # ---- 4. what L would be on human labels -----------------------------
    print("\n4. FAVOURABILITY BOUND ON HUMAN LABELS, PER NODE")
    per_node = defaultdict(lambda: {"k": 0, "n": 0})
    for r in rows:
        hc = human(r)
        if not hc or r["item"] not in key:
            continue
        nid = key[r["item"]]["intervention_node"]
        if hc == FAVOURS:
            per_node[nid]["k"] += 1
            per_node[nid]["n"] += 1
        elif hc == AGAINST:
            per_node[nid]["n"] += 1
    published = {r["node_id"]: float(r["L_corrected"]) for r in
                 csv.DictReader((HERE.parent / "outputs" / "ijtb_node_counts.csv")
                                .open(encoding="utf-8"))}
    print(f"   {'node':>5}{'k/n (human)':>14}{'L human':>10}{'L published':>13}"
          f"{'difference':>12}   label")
    diffs = []
    for nid in sorted(per_node, key=int):
        d = per_node[nid]
        if d["n"] == 0:
            continue
        lh = m.lower_beta_bound(d["k"], d["n"])
        lp = published.get(nid, float("nan"))
        diffs.append(lp - lh)
        print(f"   {nid:>5}{f'{d[chr(107)]}/{d[chr(110)]}':>14}{lh:>10.3f}{lp:>13.3f}"
              f"{lp-lh:>+12.3f}   {m.NODE_LABELS[nid][:34]}")
    if diffs:
        mean_d = sum(diffs) / len(diffs)
        print(f"\n   mean signed difference (published - human): {mean_d:+.3f}")
        print("   Positive means the published bounds are optimistic relative to human")
        print("   coding; the sampled n per node is small, so read this as a direction")
        print("   and rough magnitude, not a corrected estimate.")

    print("\nText for the manuscript can be drafted from the figures above; the")
    print("relevant limitation paragraph is Section 4.4.")


if __name__ == "__main__":
    main()
