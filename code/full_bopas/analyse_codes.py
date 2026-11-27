"""Merge and check the AI codes of the TB-BOPAS corpus, summarise the evidence per intervention and
recompute the score with the evidence component L taken from the AI codes.

Usage: python -I code/full_bopas/analyse_codes.py <exploration_dir>
Reads <dir>/pairs_all.jsonl and <dir>/codes/*.codes.jsonl; writes <dir>/results/.
"""
import csv
import json
import math
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import build_tb_bopas_ijtb as B  # noqa: E402  (score functions; main() is guarded)

CODES = ("favours intervention", "no effect or favours control", "unclear", "not an evaluation of this intervention")
FAV, NULL, UNC, NOT = CODES
TOP_DESIGNS = ("systematic review/meta-analysis", "randomised trial")
# strict: only designs with a comparator or a reference standard count towards direction
COMPARATIVE = TOP_DESIGNS + ("non-randomised comparative (cohort, quasi-experimental, before-after, interrupted time series)",
                             "diagnostic accuracy", "economic evaluation")


def norm(s: str) -> str:
    return " ".join(s.lower().split())


def wilson(k: int, n: int, z: float = 1.96):
    if n == 0:
        return (float("nan"),) * 3
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, c - h, c + h


def load(d: Path):
    pairs = {p["pair_id"]: p for p in map(json.loads, (d / "pairs_all.jsonl").read_text(encoding="utf-8").splitlines())}
    codes, problems = {}, Counter()
    # batch_* first; recode_* (full text of abstracts the first run cut at 3,500 characters) override them
    # zscope_* (first-round evaluations rechecked against the written scope rules) override both
    files = (sorted((d / "codes").glob("batch_*.codes.jsonl")) + sorted((d / "codes").glob("recode_*.codes.jsonl"))
             + sorted((d / "codes").glob("zscope_*.codes.jsonl")))
    for f in files:
        recode = f.name.startswith(("recode_", "zscope_"))
        for line in f.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            try:
                c = json.loads(line)
            except json.JSONDecodeError:
                problems["bad json line"] += 1
                continue
            pid = int(c.get("pair_id", -1))
            if pid not in pairs:
                problems["unknown pair_id"] += 1
                continue
            if c.get("code") not in CODES:
                problems["invalid code"] += 1
                continue
            if pid in codes:
                problems["recoded with full text" if recode else "duplicate pair_id"] += 1
            codes[pid] = c
    for pid, c in codes.items():
        p = pairs[pid]
        q = norm(c.get("quote", ""))
        c["quote_verbatim"] = bool(q) and q in norm(p["title"] + " " + p["abstract"])
    return pairs, codes, problems


def main() -> None:
    d = Path(sys.argv[1])
    out = d / "results"
    out.mkdir(exist_ok=True)
    pairs, codes, problems = load(d)
    coded_all = [(pairs[i], c) for i, c in sorted(codes.items())]
    # one study once per intervention: merge records whose normalised titles match within a node
    seen, coded, dups = set(), [], 0
    for p, c in coded_all:
        key = (p["node_id"], "".join(ch for ch in p["title"].lower() if ch.isalnum())[:90])
        if key in seen:
            dups += 1
            continue
        seen.add(key)
        coded.append((p, c))
    lines = [f"pairs {len(pairs)}; coded {len(codes)}; duplicate titles merged {dups}; analysed {len(coded)}; with abstract "
             f"{sum(not p['title_only'] for p, _ in coded)}; problems {dict(problems)}",
             f"quotes verbatim: {sum(c['quote_verbatim'] for _, c in coded)}/{len(coded)}"]

    by = defaultdict(list)
    for p, c in coded:
        by[p["node_id"]].append((p, c))
    node_counts, strict_counts, rows = {}, {}, []
    for nid, lbl in B.NODE_LABELS.items():
        items = by.get(nid, [])
        cc = Counter(c["code"] for _, c in items)
        n = len(items)
        ev = n - cc[NOT]
        k, nd = cc[FAV], cc[FAV] + cc[NULL]
        top = [c for _, c in items if c["code"] != NOT and c.get("design") in TOP_DESIGNS]
        ktop = sum(c["code"] == FAV for c in top)
        ndtop = sum(c["code"] in (FAV, NULL) for c in top)
        strict = [c for _, c in items if c["code"] != NOT and c.get("design") in COMPARATIVE]
        kst = sum(c["code"] == FAV for c in strict)
        ndst = sum(c["code"] in (FAV, NULL) for c in strict)
        strict_counts[nid] = {"k_favours": kst, "n_directional": ndst, "n_total": n}
        title_only_ev = sum(p["title_only"] and c["code"] != NOT for p, c in items)
        india = sum(c["code"] != NOT and c.get("setting") == "India" for _, c in items)
        share, lo, hi = wilson(ev, n)
        node_counts[nid] = {"k_favours": k, "n_directional": nd, "n_total": n}
        rows.append({"node_id": nid, "intervention": lbl, "pairs_coded": n, "evaluations": ev,
                     "evaluation_share": f"{share:.3f}", "share_lo": f"{lo:.3f}", "share_hi": f"{hi:.3f}",
                     "favours": k, "no_effect_or_control": cc[NULL], "unclear": cc[UNC],
                     "favourable_share": f"{(k / nd if nd else float('nan')):.3f}",
                     "L_ai": f"{B.lower_beta_bound(k, nd):.3f}",
                     "comparative_evaluations": len(strict), "L_ai_strict": f"{B.lower_beta_bound(kst, ndst):.3f}",
                     "SR_or_RCT_evaluations": len(top), "L_ai_SR_RCT": f"{B.lower_beta_bound(ktop, ndtop):.3f}",
                     "india_evaluations": india, "title_only_evaluations": title_only_ev})
    with (out / "evidence_by_intervention.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    by_series, _ = B.load_series()
    traj = B.fit_trajectories(by_series)
    grid, _ = B.compute_bopas(node_counts, traj)
    with (out / "bopas_ai_evidence.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(grid[0]))
        w.writeheader()
        w.writerows(grid)
    nz = [g for g in grid if g["BOPAS"] > 0]
    lines.append(f"score: {len(nz)} non-zero cells of {len(grid)}; maximum {max(g['BOPAS'] for g in grid):.3f}")
    lines += [f"  {g['intervention']} x {g['target']}: {g['BOPAS']:.3f} (L {g['L']}, T {g['T']})" for g in nz[:10]]

    gs, _ = B.compute_bopas(strict_counts, traj)
    with (out / "bopas_ai_evidence_strict.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(gs[0]))
        w.writeheader()
        w.writerows(gs)
    lines.append(f"strict score: {sum(g['BOPAS'] > 0 for g in gs)} non-zero; maximum {max(g['BOPAS'] for g in gs):.3f}")
    # abstracts that may belong to another record: few title words in the abstract
    def words(s):
        return {w for w in __import__("re").findall(r"[a-z]{5,}", s.lower())}
    mism = [p["pair_id"] for p, c in coded if not p["title_only"] and len(words(p["title"])) >= 4
            and len(words(p["title"]) & words(p["abstract"])) / len(words(p["title"])) < 0.3]
    lines.append(f"possible title-abstract mismatches: {len(mism)} {mism[:30]}")
    designs = Counter(c.get("design") for _, c in coded if c["code"] != NOT)
    lines.append("designs among evaluations: " + json.dumps(designs.most_common()))
    lines.append("codes overall: " + json.dumps(Counter(c["code"] for _, c in coded).most_common()))
    (out / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # random spot-check sample for the authors: stratified by code, fixed seed, AI code withheld
    rng = random.Random(20261127)
    quota = {FAV: 20, NULL: 10, UNC: 10, NOT: 20}
    sample = []
    for code, q in quota.items():
        pool = [(p, c) for p, c in coded if c["code"] == code and not p["title_only"]]
        sample += rng.sample(pool, min(q, len(pool)))
    rng.shuffle(sample)
    with (out / "spot_check_key.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["item", "pair_id", "ai_code"])
        for i, (p, c) in enumerate(sample, 1):
            w.writerow([i, p["pair_id"], c["code"]])
    import openpyxl
    from openpyxl.worksheet.datavalidation import DataValidation
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "spot_check"
    ws.append(["item", "pair_id", "intervention", "title", "abstract", "coder_A", "coder_B", "agreed_code"])
    for i, (p, c) in enumerate(sample, 1):
        ws.append([i, p["pair_id"], p["intervention"], p["title"], p["abstract"], "", "", ""])
        for cell in ws[ws.max_row][2:5]:          # harvested text is data, never a formula
            cell.data_type = "s"
    dv = DataValidation(type="list", formula1='"' + ",".join(CODES) + '"', allow_blank=True)
    ws.add_data_validation(dv)
    dv.add(f"F2:H{len(sample) + 1}")
    for col, wdt in zip("ABCDEFGH", (6, 8, 30, 50, 90, 28, 28, 28)):
        ws.column_dimensions[col].width = wdt
    wb.save(out / "spot_check_worksheet.xlsx")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
