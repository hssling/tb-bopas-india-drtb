"""Term search of India's animal-health annual reports for bovine tuberculosis.

The reports are public documents and are not redistributed here. Download them from the URLs in
SOURCES (retrieved 26 November 2026), place them in a folder, and run:
    python -I code/animal_sector_search.py FOLDER
Requires PyMuPDF (fitz). Writes outputs/animal_sector_search.json with page counts, SHA-256 digests,
term counts and the context of every tuberculosis mention. Brucellosis and foot-and-mouth disease,
both under national control programmes, are counted as controls showing that the text was searchable.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import fitz

ROOT = Path(__file__).resolve().parent.parent
SOURCES = {
    "DAHD_2023-24.pdf": "https://www.dahd.gov.in/static/uploads/2025/08/62243e99836ac752f30c37a792bc997f.pdf",
    "DAHD_2024-25.pdf": "https://www.dahd.gov.in/static/uploads/2025/08/e62dde0085200270e6dd8bc455528ffe.pdf",
    "DAHD_2025-26.pdf": "https://www.dahd.gov.in/static/uploads/2026/04/29232c53442fadce97015a62fd715e74.pdf",
    "NIVEDI_2024.pdf": "https://nivedi.res.in/pdf/Reports/NIVEDI_Annual_Report_2024.pdf",
    "NIVEDI_2025.pdf": "https://nivedi.res.in/pdf/Reports/Annual%20Report%202025.pdf",
}
TERMS = {
    "tuberculosis": r"(?<!para)tubercul(?:osis|ar)",
    "bovine_tb": r"bovine\s+t(?:uberculosis|b)\b|\bbTB\b",
    "m_bovis": r"m(?:ycobacterium|\.)\s*bovis",
    "tuberculin": r"tuberculin",
    "brucellosis_control": r"brucell",
    "fmd_control": r"foot[\s-]*and[\s-]*mouth|\bFMD\b",
}


def main(folder: Path) -> None:
    res = {}
    for name, url in SOURCES.items():
        p = folder / name
        doc = fitz.open(p)
        pages = [pg.get_text() for pg in doc]
        text = "\n".join(pages)
        res[name] = {
            "url": url,
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
            "pages": len(pages),
            "pages_with_text": sum(len(x.strip()) > 200 for x in pages),
            "counts": {k: len(re.findall(v, text, re.I)) for k, v in TERMS.items()},
            "tuberculosis_contexts": [re.sub(r"\s+", " ", text[max(0, m.start() - 160):m.end() + 160])
                                      for m in re.finditer(TERMS["tuberculosis"], text, re.I)],
        }
    out = ROOT / "outputs" / "animal_sector_search.json"
    out.write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")
    for k, v in res.items():
        print(k, v["pages"], v["counts"])


if __name__ == "__main__":
    main(Path(sys.argv[1]))
