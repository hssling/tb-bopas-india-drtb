"""India's answers in WHO's zoonotic tuberculosis dataset, with the global context.

Reads data/who_zoonotic_tb_2026-11-26.csv (WHO Global Tuberculosis Programme, dataset "zoonotic",
https://extranet.who.int/tme/generateCSV.asp?ds=zoonotic, retrieved 26 November 2026; collected
for 2018 and 2023) and writes outputs/zoonotic_reporting.json. Standard library only.
Run: python -I code/zoonotic_reporting.py
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "data" / "who_zoonotic_tb_2026-11-26.csv"
OUT = ROOT / "outputs" / "zoonotic_reporting.json"
FIELDS = ["pulm_spec_available", "newrel_pulm_spec", "newrel_pulm_mbovis", "ep_spec_available",
          "newrel_ep_spec", "newrel_ep_mbovis", "ident_zoonotic", "zoonotic", "disagg_zoonotic",
          "m_bovis", "m_orygis", "m_caprae", "m_pinnipedii"]


def has(v: str | None) -> bool:
    return v not in (None, "", "NA")


def main() -> None:
    rows = list(csv.DictReader(SRC.open(encoding="utf-8", newline="")))
    india = {r["year"]: {f: r[f] for f in FIELDS if has(r[f])} for r in rows if r["iso3"] == "IND"}
    by_year = {}
    for y in sorted({r["year"] for r in rows}):
        yr = [r for r in rows if r["year"] == y]
        by_year[y] = {
            "countries_in_file": len(yr),
            "answered_ident_zoonotic": sum(has(r["ident_zoonotic"]) for r in yr),
            "ident_zoonotic_yes": sum(r["ident_zoonotic"] == "1" for r in yr),
            "reported_zoonotic_count_above_zero": sum(has(r["zoonotic"]) and r["zoonotic"] != "0" for r in yr),
            "pulm_speciation_available": sum(r["pulm_spec_available"] == "1" for r in yr),
        }
    res = {
        "source": SRC.name,
        "sha256": hashlib.sha256(SRC.read_bytes()).hexdigest(),
        "years_collected": sorted({r["year"] for r in rows}),
        "india": india,
        "all_countries_by_year": by_year,
        "country_years_answering_ident_zoonotic": sum(has(r["ident_zoonotic"]) for r in rows),
        "country_years_ident_zoonotic_yes": sum(r["ident_zoonotic"] == "1" for r in rows),
        "countries_reporting_zoonotic_count_above_zero": len({r["iso3"] for r in rows
                                                              if has(r["zoonotic"]) and r["zoonotic"] != "0"}),
    }
    OUT.write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
