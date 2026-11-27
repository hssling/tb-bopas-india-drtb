"""Agreement between the final AI codes and the authors' adjudicated spot-check codes (60 random pairs).

Reads <dir>/human_checks/ (filled worksheet, key, adjudication) and the final AI codes (same override
order as analyse_codes.py). Reports Cohen's kappa for the four codes and for evaluation-or-not, and lists
pairs where the authors coded an evaluation but the AI coded a tier-2 (relevant, not evaluating) record.
Usage: python -I code/full_bopas/spot_check_agreement.py <exploration_dir>
"""
import csv
import json
import sys
from collections import Counter
from pathlib import Path

import openpyxl

NOT = "not an evaluation of this intervention"


def kappa(a, b):
    n = len(a)
    po = sum(x == y for x, y in zip(a, b)) / n
    ca, cb = Counter(a), Counter(b)
    pe = sum(ca[k] * cb[k] for k in set(a) | set(b)) / (n * n)
    return (po - pe) / (1 - pe), po


def final_codes(d: Path) -> dict:
    codes = {}
    for pat in ("batch_*", "recode_*", "zscope_*", "ztier_*", "zz*"):
        for f in sorted((d / "codes").glob(pat + ".codes.jsonl")):
            for line in f.read_text(encoding="utf-8").splitlines():
                if line.strip():
                    c = json.loads(line)
                    codes[int(c["pair_id"])] = c
    return codes


def main() -> None:
    d = Path(sys.argv[1])
    h = d / "human_checks"
    ws = openpyxl.load_workbook(h / "spot_check_FILLED_2026-11-27.xlsx").active
    hdr = [c.value for c in ws[1]]
    rows = [dict(zip(hdr, [c.value for c in r])) for r in ws.iter_rows(min_row=2)]
    wa = openpyxl.load_workbook(h / "spot_check_adjudication_FILLED_2026-11-27.xlsx").active
    ha = [c.value for c in wa[1]]
    adj = {int(r[0].value): (r[ha.index("FINAL_code")].value, r[ha.index("why_they_differ")].value)
           for r in wa.iter_rows(min_row=2)}
    codes = final_codes(d)
    human = [(adj.get(int(r["item"]), (r["agreed_code"],))[0] or "").strip() for r in rows]
    ai = [codes[int(r["pair_id"])]["code"] for r in rows]
    t2 = [codes[int(r["pair_id"])].get("tier2") is True for r in rows]
    k4, p4 = kappa(ai, human)
    ke, pe = kappa([c != NOT for c in ai], [c != NOT for c in human])
    out = [f"60 spot-check pairs, final AI codes vs authors' adjudicated codes",
           f"  four codes: kappa {k4:.2f}, agreement {p4:.0%}",
           f"  evaluation or not: kappa {ke:.2f}, agreement {pe:.0%}"]
    rel = [r for r, a, hc, t in zip(rows, ai, human, t2) if a == NOT and t and hc != NOT]
    out.append(f"  authors coded an evaluation where the AI coded tier 2 (relevant, not evaluating): {len(rel)}")
    out += [f"    item {r['item']} pair {r['pair_id']}: {str(r['title'])[:110]}" for r in rel]
    keep = [i for i, r in enumerate(rows) if r not in rel]
    k2, p2 = kappa([ai[i] for i in keep], [human[i] for i in keep])
    out.append(f"  excluding those (n={len(keep)}): kappa {k2:.2f}, agreement {p2:.0%}")
    out.append("  confusion (AI -> authors): " + json.dumps(Counter(zip(ai, human)).most_common()))
    (d / "results" / "spot_check_agreement.txt").write_text("\n".join(out) + "\n", encoding="utf-8")
    print("\n".join(out))


if __name__ == "__main__":
    main()
