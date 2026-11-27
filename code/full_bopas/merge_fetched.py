"""Attach fetched PubMed abstracts to corpus pairs and batch the pairs that need (re)coding.

A fetched abstract is used only when its PMID matches and its title agrees with the corpus title (token
overlap >= 0.6), so an abstract can never attach to another record. Pairs that gain an abstract (title-only
records and suspected mismatches) go to batches F01..; pairs still without an abstract go to T_rest
(coded from the title). Outputs are named recode_* when coded so they override earlier codes.
Usage: python -I code/full_bopas/merge_fetched.py <exploration_dir> <mismatch_pair_ids_comma_separated>
"""
import json
import re
import sys
from pathlib import Path


def toks(s: str) -> set[str]:
    return set(re.findall(r"[a-z0-9]{3,}", s.lower()))


def main() -> None:
    d = Path(sys.argv[1])
    mism = {int(x) for x in sys.argv[2].split(",") if x}
    pairs = [json.loads(l) for l in (d / "pairs_all.jsonl").read_text(encoding="utf-8").splitlines()]
    fetched = {}
    for g in map(json.loads, (d / "fetched" / "pubmed_abstracts.jsonl").read_text(encoding="utf-8").splitlines()):
        fetched[g["pmid"]] = g
    gained, still, rejected = [], [], 0
    for p in pairs:
        if not (p["title_only"] or p["pair_id"] in mism):
            continue
        g = fetched.get(p["pmid"].strip())
        ok = False
        if g and g["abstract"]:
            a, b = toks(p["title"]), toks(g["title"])
            ok = bool(a) and len(a & b) / len(a) >= 0.6
            rejected += not ok
        if ok:
            p["abstract"], p["title_only"], p["abstract_source"] = g["abstract"][:12000], False, "pubmed_efetch_xml"
            gained.append(p)
        elif p["title_only"]:
            still.append(p)
        else:                                  # suspected mismatch without a verified replacement: title only
            p["abstract"], p["title_only"], p["abstract_source"] = "", True, "dropped_mismatch"
            still.append(p)
    (d / "pairs_all.jsonl").write_text("".join(json.dumps(p, ensure_ascii=False) + "\n" for p in pairs), encoding="utf-8")
    for j in range(0, len(gained), 150):
        (d / "batches" / f"batch_F{j // 150 + 1:02d}.jsonl").write_text(
            "".join(json.dumps(p, ensure_ascii=False) + "\n" for p in gained[j:j + 150]), encoding="utf-8")
    (d / "batches" / "batch_T_rest.jsonl").write_text(
        "".join(json.dumps(p, ensure_ascii=False) + "\n" for p in still), encoding="utf-8")
    print(f"gained an abstract {len(gained)} (of which suspected mismatches {sum(p['pair_id'] in mism for p in gained)}); "
          f"title only {len(still)}; title disagreement rejected {rejected}; F batches {-(-len(gained) // 150)}")


if __name__ == "__main__":
    main()
