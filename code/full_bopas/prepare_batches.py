"""Prepare (record, intervention) pairs of the TB-BOPAS corpus for one round of AI coding.

Reads data/tb_interventions_unified.csv and outputs/ijtb_node_counts.csv (intervention labels), strips
PubMed citation headers and author-affiliation blocks, and writes JSONL batches to the exploration
folder given as the first argument. Records without an abstract are coded from the title.
Usage: python -I code/full_bopas/prepare_batches.py <out_dir> [batch_size]
"""
import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CAP = 12000  # characters; 3,500 (first run) cut the results off 165 long abstracts
HEAD = r"(BACKGROUND|INTRODUCTION|OBJECTIVES?|AIMS?|PURPOSE|METHODS?|SETTING|RATIONALE|CONTEXT|IMPORTANCE|SUMMARY|DESIGN)\s*:"


def clean(title: str, text: str) -> str:
    t = text.strip()
    if not t:
        return ""
    i = t.find(title[:60]) if title else -1
    if i >= 0:                                   # drop the citation header before the title
        t = t[i + len(title):]
    m = re.search(r"Author information:", t)
    if m:                                        # drop affiliations up to the first section heading
        h = re.search(HEAD, t[m.end():])
        t = t[:m.start()] + (t[m.end() + h.start():] if h else t[m.end():])
    t = re.split(r"\s(?:©|Copyright ©|DOI: 10\.|PMID: \d)", t)[0]
    return re.sub(r"\s+", " ", t).strip()[:CAP]


def main() -> None:
    out = Path(sys.argv[1])
    size = int(sys.argv[2]) if len(sys.argv) > 2 else 150
    labels = {r["node_id"]: r["label"] for r in csv.DictReader(open(ROOT / "outputs/ijtb_node_counts.csv", encoding="utf-8"))}
    rows = list(csv.DictReader(open(ROOT / "data/tb_interventions_unified.csv", encoding="utf-8")))
    pairs = []
    for k, r in enumerate(rows, 1):
        text = clean(r["title"], r["abstract_present"])
        pairs.append({"pair_id": k, "source": r["source"], "node_id": r["node_id"], "intervention": labels[r["node_id"]],
                      "pmid": r["pmid"], "doi": r["doi"], "year": r["year"], "title": r["title"].strip(),
                      "abstract": text, "title_only": not bool(text)})
    if len(sys.argv) > 3 and sys.argv[3] == "--recode-capped":
        # rewrite pairs_all.jsonl with full text and batch the pairs the 3,500 cap had cut; batches untouched
        (out / "pairs_all.jsonl").write_text("".join(json.dumps(p, ensure_ascii=False) + "\n" for p in pairs), encoding="utf-8")
        cut = [p for p in pairs if len(p["abstract"]) > 3500]
        for j in range(0, len(cut), 85):
            (out / "batches" / f"batch_R{j // 85 + 1:02d}.jsonl").write_text(
                "".join(json.dumps(p, ensure_ascii=False) + "\n" for p in cut[j:j + 85]), encoding="utf-8")
        print(f"{len(cut)} pairs longer than 3,500 characters re-batched")
        return
    with_ab = [p for p in pairs if not p["title_only"]]
    title_only = [p for p in pairs if p["title_only"]]
    (out / "pairs_all.jsonl").write_text("".join(json.dumps(p, ensure_ascii=False) + "\n" for p in pairs), encoding="utf-8")
    n = 0
    for group, tag in ((with_ab, "A"), (title_only, "T")):
        bs = size if tag == "A" else size * 4      # title-only items are short
        for j in range(0, len(group), bs):
            n += 1
            f = out / "batches" / f"batch_{tag}{j // bs + 1:02d}.jsonl"
            f.write_text("".join(json.dumps(p, ensure_ascii=False) + "\n" for p in group[j:j + bs]), encoding="utf-8")
    print(f"{len(pairs)} pairs: {len(with_ab)} with abstract, {len(title_only)} title only; {n} batch files")


if __name__ == "__main__":
    main()
