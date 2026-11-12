"""Build a 150-record adjudication worksheet for validating the automated
direction-of-effect classifier.

Reviewer 1, Major 3: "Machine coding of direction of effect is defensible;
machine coding that is not validated against a human-coded sample is not."

Design decisions that matter for the validity of the exercise:

  * The machine's label is NOT shown in the worksheet. If it were, both
    coders would anchor on it and the agreement statistic would be
    meaningless. It is written to a separate key file that the analysis
    script reads after the humans have finished.
  * Sampling is stratified by intervention node, proportional to the number
    of directionally coded records in that node, with a floor of 4 per node
    so that every node contributes something. Node 9 has only 3 directional
    records, so it contributes all of them.
  * Records are presented in a fixed shuffled order (seeded), so the two
    coders see the same sequence but node membership is not clustered.
  * Both coders fill their own column independently. A third column records
    the consensus reached on disagreements.

Outputs, in this directory:
  validation_worksheet.xlsx   what the two authors fill in
  validation_worksheet.csv    same content, for anyone who prefers plain text
  validation_key.csv          machine labels, withheld from the coders
"""
from __future__ import annotations

import csv
import importlib.util
import random
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "code" / "build_tb_bopas_ijtb.py"
ROOT = HERE.parent
CORPUS = ROOT / "data" / "tb_interventions_unified.csv"

TARGET_N = 150
FLOOR_PER_NODE = 4
SEED = 20260908

spec = importlib.util.spec_from_file_location("tb_ijtb", SRC)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

CODES = ["favours intervention", "no effect or favours control", "unclear",
         "not an evaluation of this intervention"]


def machine_label(text: str) -> str:
    """Reproduce the corrected classifier's label for one record."""
    if m.RE_UNCLEAR.search(text):
        return "unclear"
    if m.RE_CONTRARY.search(text):
        return "no effect or favours control"
    if m.RE_FAVOURS.search(text):
        return "favours intervention"
    return "no evaluative cue"


def main() -> None:
    rows = list(csv.DictReader(CORPUS.open(encoding="utf-8", newline="")))

    # every (record, node) pair the classifier codes directionally
    by_node: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        text = m._record_text(r)
        if m.RE_UNCLEAR.search(text):
            continue
        if not (m.RE_CONTRARY.search(text) or m.RE_FAVOURS.search(text)):
            continue
        for nid, rx in m.NODE_REGEX.items():
            if rx.search(text):
                by_node[nid].append(r)

    total = sum(len(v) for v in by_node.values())
    print(f"directionally coded (record, node) pairs available: {total}")

    # proportional allocation with a floor, then trim/top-up to TARGET_N
    alloc = {nid: max(FLOOR_PER_NODE, round(TARGET_N * len(v) / total))
             for nid, v in by_node.items()}
    for nid in alloc:                       # cannot sample more than exist
        alloc[nid] = min(alloc[nid], len(by_node[nid]))
    while sum(alloc.values()) > TARGET_N:   # trim the largest strata first
        nid = max(alloc, key=lambda k: (alloc[k], len(by_node[k])))
        if alloc[nid] <= FLOOR_PER_NODE:
            break
        alloc[nid] -= 1
    while sum(alloc.values()) < TARGET_N:
        nid = max(alloc, key=lambda k: len(by_node[k]) - alloc[k])
        if alloc[nid] >= len(by_node[nid]):
            break
        alloc[nid] += 1

    rng = random.Random(SEED)
    sample = []
    for nid, n_take in sorted(alloc.items(), key=lambda kv: int(kv[0])):
        for r in rng.sample(by_node[nid], n_take):
            sample.append((nid, r))
        print(f"  node {nid:>2}  {m.NODE_LABELS[nid][:46]:<48} "
              f"available {len(by_node[nid]):>4}  sampled {n_take}")
    rng.shuffle(sample)
    print(f"total sampled: {len(sample)}")

    # ---- worksheet ------------------------------------------------------
    header = ["item", "record_id", "intervention_node", "title", "abstract",
              "coder_A_code", "coder_B_code", "consensus_code", "notes"]
    ws_rows, key_rows = [], []
    for i, (nid, r) in enumerate(sample, start=1):
        rid = (r.get("doi") or "").strip() or ("PMID:" + (r.get("pmid") or "").strip())
        abstract = (r.get("abstract_present") or "").strip()
        ws_rows.append({
            "item": i,
            "record_id": rid or f"row{i}",
            "intervention_node": f"{nid} - {m.NODE_LABELS[nid]}",
            "title": (r.get("title") or "").strip(),
            "abstract": abstract if abstract else
                        "[no abstract in the deposited corpus - code from the title alone]",
            "coder_A_code": "", "coder_B_code": "", "consensus_code": "", "notes": "",
        })
        key_rows.append({"item": i, "record_id": rid or f"row{i}",
                         "intervention_node": nid,
                         "machine_code": machine_label(m._record_text(r)),
                         "title_only": "yes" if not abstract else "no"})

    with (HERE / "validation_worksheet.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=header)
        w.writeheader()
        w.writerows(ws_rows)
    with (HERE / "validation_key.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(key_rows[0].keys()))
        w.writeheader()
        w.writerows(key_rows)

    # ---- xlsx, with dropdowns so the codes stay clean -------------------
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter
    from openpyxl.worksheet.datavalidation import DataValidation

    wb = Workbook()
    sh = wb.active
    sh.title = "worksheet"

    instructions = [
        "Direction-of-effect validation worksheet",
        "",
        "One row per (record, intervention node) pair. For each row, read the title and abstract "
        "and decide what the study reports ABOUT THAT INTERVENTION, using the dropdown in your "
        "own column:",
        "",
        "    favours intervention                    - the study reports a benefit for this intervention",
        "    no effect or favours control            - no benefit, a null result, or harm",
        "    unclear                                 - mixed, heterogeneous, or genuinely undecidable",
        "    not an evaluation of this intervention  - the record does not evaluate this intervention at all",
        "",
        "Coder A and Coder B fill their own columns INDEPENDENTLY, without discussion and without "
        "looking at each other's column. Only afterwards, work through the rows where you differ "
        "and record the agreed label in consensus_code.",
        "",
        "The classifier's own label is deliberately not shown here - it is held in "
        "validation_key.csv so that it cannot anchor your judgment. Once both coder columns and "
        "the consensus column are filled, run:",
        "",
        "    python analyse_validation.py",
        "",
        "which reports inter-coder agreement, classifier accuracy against the human consensus, "
        "and what the favourability bounds would be if the human labels were used.",
    ]
    for line in instructions:
        sh.append([line])
    sh.append([])
    start = sh.max_row + 1

    sh.append(header)
    hdr_font = Font(bold=True, color="FFFFFF")
    hdr_fill = PatternFill("solid", fgColor="1F6F8B")
    for c in range(1, len(header) + 1):
        cell = sh.cell(row=start, column=c)
        cell.font = hdr_font
        cell.fill = hdr_fill
        cell.alignment = Alignment(vertical="center")
    for r in ws_rows:
        sh.append([r[h] for h in header])

    widths = {"item": 6, "record_id": 26, "intervention_node": 34, "title": 60,
              "abstract": 90, "coder_A_code": 30, "coder_B_code": 30,
              "consensus_code": 30, "notes": 30}
    for i, h in enumerate(header, start=1):
        sh.column_dimensions[get_column_letter(i)].width = widths[h]
    for row in sh.iter_rows(min_row=start + 1, max_row=sh.max_row):
        for cell in row:
            cell.alignment = Alignment(wrap_text=True, vertical="top")

    dv = DataValidation(type="list", formula1='"' + ",".join(CODES) + '"',
                        allow_blank=True)
    dv.error = "Choose one of the four listed codes."
    dv.prompt = "Pick the code that matches what the study reports about this intervention."
    sh.add_data_validation(dv)
    first, last = start + 1, start + len(ws_rows)
    for col in ("F", "G", "H"):
        dv.add(f"{col}{first}:{col}{last}")
    sh.freeze_panes = f"A{start + 1}"

    key = wb.create_sheet("codes")
    key.append(["code", "meaning"])
    meanings = ["the study reports a benefit for this intervention",
                "no benefit, a null result, or harm",
                "mixed, heterogeneous, or genuinely undecidable",
                "the record does not evaluate this intervention at all"]
    for c, mean in zip(CODES, meanings):
        key.append([c, mean])
    key.column_dimensions["A"].width = 40
    key.column_dimensions["B"].width = 60

    wb.save(HERE / "validation_worksheet.xlsx")
    print(f"\nwrote validation_worksheet.xlsx ({len(ws_rows)} rows), "
          f"validation_worksheet.csv and validation_key.csv")
    n_title_only = sum(1 for k in key_rows if k["title_only"] == "yes")
    print(f"records in the sample with no abstract (title-only coding): {n_title_only}")


if __name__ == "__main__":
    main()
