#!/usr/bin/env python3
"""Trend component of the priority score under alternative analytic choices.

The literature term L failed validation, so the paper reports only the parts
of the score that come from surveillance: the trend term
    T = max(signed slope / max(half-width, floor), 0)
(signed so that positive means worsening; half-width = half the bootstrap
95% range), and the number of intervention-indicator cells that can score
above zero, i.e. cells with T > 0 and a non-zero causal link M. This is the
rule in compute_bopas() of the pipeline, without L and C.

Reads ijtb_slopes.csv and the M_TB literal of the pipeline script (read-only).
Usage: python -I code/trend_component.py outputs/trend_component.json
"""
from __future__ import annotations

import ast
import csv
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent                                     # repository root
SLOPES = ROOT / "outputs" / "ijtb_slopes.csv"
CODE = ROOT / "code" / "build_tb_bopas_ijtb.py"
FLOOR = 0.5


def m_matrix():
    src = CODE.read_text(encoding="utf-8")
    return ast.literal_eval(re.search(r"^M_TB\s*=\s*(\{.*?^\})", src, re.S | re.M).group(1))


def t_value(slope, lo, hi, sign, floor):
    half = abs(hi - lo) / 2
    used = max(half, floor) if floor is not None else half
    return max(slope * sign / used, 0.0), half, used


def scenario(rows, M, *, naive=False, floor=FLOOR, strict=False):
    T = {}
    detail = {}
    for r in rows:
        s = r["series"]; sign = int(r["sign"])
        if naive:
            slope, lo, hi = float(r["slope_full"]), float(r["boot_lo"]), float(r["boot_hi"])
        else:
            slope, lo, hi = float(r["slope_primary"]), float(r["primary_lo"]), float(r["primary_hi"])
        t, half, used = t_value(slope, lo, hi, sign, floor)
        T[s] = t
        detail[s] = dict(slope=slope, half_width=round(half, 4), half_used=round(used, 4),
                         floor_binding=floor is not None and half < floor, T=round(t, 3))
    cells = 0
    for node, links in M.items():
        for s, m in links.items():
            if strict and m == 0.5:
                m = 0.0
            if T.get(s, 0) > 0 and m > 0:
                cells += 1
    top = max(T, key=T.get)
    return dict(T=detail, cells=cells, top_series=top, top_T=round(T[top], 3))


def main():
    rows = list(csv.DictReader(SLOPES.open(encoding="utf-8")))
    M = m_matrix()
    res = {
        "primary": scenario(rows, M),
        "no_floor": scenario(rows, M, floor=None),
        "one_line": scenario(rows, M, naive=True),
        "one_line_no_floor": scenario(rows, M, naive=True, floor=None),
        "partial_links_zero": scenario(rows, M, strict=True),
    }
    res["ratio_one_line_vs_primary"] = round(res["one_line"]["top_T"] / res["primary"]["top_T"], 2)
    res["ratio_no_floor"] = round(res["one_line_no_floor"]["top_T"] / res["no_floor"]["top_T"], 2)
    Path(sys.argv[1]).write_text(json.dumps(res, indent=1), encoding="utf-8")
    for k, v in res.items():
        if isinstance(v, dict):
            print(f"{k:<20} cells {v['cells']:>2}  top {v['top_series']:<18} T {v['top_T']}")
        else:
            print(k, v)


if __name__ == "__main__":
    main()
