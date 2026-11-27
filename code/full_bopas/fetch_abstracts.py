"""Fetch PubMed abstracts for corpus pairs that have no abstract, or whose stored abstract may belong to
another record (few title words in the abstract), so they can be coded from full text.

Uses NCBI E-utilities efetch (retmode=xml) and keys every abstract by its own PMID and title, which avoids
the order-based mapping that produced mismatches in the original text harvest. Generic tool name; no
e-mail or personal data is sent. Writes <dir>/fetched/pubmed_abstracts.jsonl.
Usage: python -I code/full_bopas/fetch_abstracts.py <exploration_dir> [extra_pair_ids_comma_separated]
"""
import json
import sys
import time
import urllib.parse
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UA = "tb-bopas-research/1.0 (python urllib)"


def text_of(el) -> str:
    return " ".join("".join(el.itertext()).split()) if el is not None else ""


def fetch(pmids: list[str], tries: int = 4) -> list[dict]:
    for k in range(tries):                   # NCBI returns transient 502/429; back off and retry
        try:
            return _fetch(pmids)
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
            if k == tries - 1:
                raise
            time.sleep(3 * (k + 1))


def _fetch(pmids: list[str]) -> list[dict]:
    data = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(pmids), "retmode": "xml",
                                   "tool": "tb_bopas_research"}).encode()
    req = urllib.request.Request(EFETCH, data=data, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=120) as r:
        root = ET.fromstring(r.read())
    out = []
    for art in root.iter("PubmedArticle"):
        pmid = text_of(art.find(".//MedlineCitation/PMID"))
        title = text_of(art.find(".//Article/ArticleTitle"))
        parts = []
        for ab in art.findall(".//Article/Abstract/AbstractText"):
            label = ab.get("Label")
            parts.append((label + ": " if label else "") + text_of(ab))
        out.append({"pmid": pmid, "title": title, "abstract": " ".join(parts)})
    return out


def main() -> None:
    d = Path(sys.argv[1])
    extra = {int(x) for x in sys.argv[2].split(",")} if len(sys.argv) > 2 and sys.argv[2] else set()
    pairs = [json.loads(l) for l in (d / "pairs_all.jsonl").read_text(encoding="utf-8").splitlines()]
    want = sorted({p["pmid"].strip() for p in pairs if p["pmid"].strip() and (p["title_only"] or p["pair_id"] in extra)})
    (d / "fetched").mkdir(exist_ok=True)
    got = []
    for i in range(0, len(want), 100):
        got += fetch(want[i:i + 100])
        time.sleep(0.5)                      # stay under NCBI's 3 requests per second
    path = d / "fetched" / "pubmed_abstracts.jsonl"
    path.write_text("".join(json.dumps(g, ensure_ascii=False) + "\n" for g in got), encoding="utf-8")
    print(f"requested {len(want)} PMIDs; received {len(got)} records; with abstract "
          f"{sum(bool(g['abstract']) for g in got)}; {path.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
