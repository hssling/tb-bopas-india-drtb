"""
harvest_tb_interventions.py
===========================

Intervention-focused open-access harvest for the TB-BOPAS manuscript.

Why this exists
---------------
The existing mycobacterial One Health corpus at D:/One Health TB/...
was built with *organism-focused* queries (M. bovis, M. orygis, M. avium,
M. abscessus, MTBC-animal). That corpus underweights *intervention-focused*
literature (BPaL trials, Xpert deployment studies, airborne IPC trials,
TPT trials), which the TB-BOPAS score needs to populate its per-node
favours-intervention counts.

This script runs 12 TB-intervention-specific Boolean queries against
PubMed E-utilities and Europe PMC (OPEN_ACCESS:Y), deduplicates the
unioned set by DOI then by normalised title, and writes a single CSV
ready to feed into build_tb_bopas.py via the --corpus override.

Outputs
-------
data_raw/tb_interventions_pubmed.csv     per-node PubMed hits (raw)
data_raw/tb_interventions_europepmc.csv  per-node Europe PMC hits (raw)
data_processed/tb_interventions_unified.csv   deduplicated unified corpus
logs/tb_interventions_harvest.log         per-node hit counts

Usage
-----
    python code/harvest_tb_interventions.py

Then re-run build_tb_bopas.py with:
    python code/build_tb_bopas.py --corpus data_processed/tb_interventions_unified.csv

Polite to the APIs
------------------
PubMed E-utilities is throttled to 3 requests/second per NCBI guidance;
Europe PMC has no published rate limit but we sleep 0.4 s between
requests. Each node fetches up to retmax=200 records, so the maximum
total volume is 12 nodes x 200 = 2 400 records per source.
"""
from __future__ import annotations
import csv, json, time, re, hashlib
from pathlib import Path
from urllib.parse import quote
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

ROOT = Path(__file__).resolve().parent.parent
RAW  = ROOT / "data_raw";         RAW.mkdir(exist_ok=True)
PROC = ROOT / "data_processed";   PROC.mkdir(exist_ok=True)
LOGS = ROOT / "logs";              LOGS.mkdir(exist_ok=True)

UA = "TB-BOPAS-harvest/1.0 (research; contact: hssling@yahoo.com)"

# Date window: BPaL approval (2019) to current; reduces irrelevant historical hits.
DATE_FROM = "2017"
DATE_TO   = "2026"

# 12 TB-intervention Boolean queries, paired with the BOPAS node id.
# Each query combines: pathogen scope (tuberculosis) AND intervention vocabulary
# AND comparative-design language. Queries are tuned per source syntax.
TB_QUERIES_PUBMED: list[tuple[str, str, str]] = [
    ("1",  "shorter_regimen",
     '(tuberculosis[tiab] AND ('
     'BPaL[tiab] OR "BPaL-M"[tiab] OR bedaquiline[tiab] OR pretomanid[tiab] OR '
     '"all-oral regimen"[tiab] OR "shorter regimen"[tiab] OR "Nix-TB"[tiab] OR '
     'ZeNix[tiab] OR PRACTECAL[tiab])) AND ('
     'randomized[tiab] OR trial[tiab] OR "controlled study"[tiab] OR cohort[tiab] OR '
     '"interrupted time series"[tiab])'),
    ("2",  "molecular_dx",
     '("Xpert MTB/RIF"[tiab] OR "Xpert MTB/XDR"[tiab] OR GeneXpert[tiab] OR '
     '"line probe assay"[tiab] OR Truenat[tiab]) AND tuberculosis[tiab] AND ('
     'sensitivity[tiab] OR specificity[tiab] OR rollout[tiab] OR implementation[tiab] OR '
     'scale-up[tiab] OR yield[tiab])'),
    ("3",  "wgs_dst",
     '("whole genome sequencing"[tiab] OR WGS[tiab] OR "next-generation sequencing"[tiab]) AND '
     'tuberculosis[tiab] AND ("drug susceptibility"[tiab] OR resistance prediction[tiab] OR '
     '"catalogue of mutations"[tiab])'),
    ("4",  "active_case_finding",
     '("active case finding"[tiab] OR "systematic screening"[tiab] OR '
     '"contact investigation"[tiab] OR "house-to-house"[tiab]) AND tuberculosis[tiab] AND ('
     'cluster[tiab] OR controlled[tiab] OR "interrupted time series"[tiab] OR "before-after"[tiab])'),
    ("5",  "tpt",
     '("TB preventive therapy"[tiab] OR "TB preventive treatment"[tiab] OR TPT[tiab] OR '
     '3HP[tiab] OR 1HP[tiab] OR "isoniazid preventive"[tiab] OR "levofloxacin prophylaxis"[tiab]) AND ('
     'tuberculosis[tiab] AND (randomized[tiab] OR trial[tiab] OR "controlled study"[tiab]))'),
    ("6",  "airborne_ipc",
     '("airborne infection control"[tiab] OR "airborne infection prevention"[tiab] OR '
     '"respiratory protection"[tiab] OR N95[tiab] OR "ultraviolet germicidal"[tiab] OR UVGI[tiab] OR '
     'ventilation[tiab]) AND tuberculosis[tiab]'),
    ("7",  "nutrition_support",
     '("nutritional support"[tiab] OR "food support"[tiab] OR undernutrition[tiab] OR '
     '"calorie supplementation"[tiab]) AND tuberculosis[tiab] AND ('
     'randomized[tiab] OR trial[tiab] OR cohort[tiab])'),
    ("8",  "private_sector",
     '("private sector engagement"[tiab] OR "JEET"[tiab] OR "joint effort"[tiab] OR '
     '"notification scheme"[tiab] OR "informant network"[tiab] OR "voucher scheme"[tiab]) AND '
     'tuberculosis[tiab]'),
    ("9",  "integrated_tbhivdm",
     '(("TB-HIV"[tiab] OR "TB and HIV"[tiab] OR "TB-diabetes"[tiab] OR "TB and diabetes"[tiab]) AND '
     '("integrated care"[tiab] OR collaboration[tiab] OR joint[tiab] OR coordinated[tiab]))'),
    ("10", "bovine_tb_control",
     '("bovine tuberculosis"[tiab] OR "Mycobacterium bovis"[tiab]) AND ('
     '"test-and-cull"[tiab] OR "test and slaughter"[tiab] OR vaccination[tiab] OR '
     '"BCG cattle"[tiab] OR "herd vaccination"[tiab] OR control[tiab]) AND ('
     'cattle[tiab] OR dairy[tiab] OR herd[tiab])'),
    ("11", "milk_pasteurisation",
     '("milk pasteurization"[tiab] OR "milk pasteurisation"[tiab] OR "raw milk"[tiab] OR '
     '"unpasteurised milk"[tiab] OR "unpasteurized milk"[tiab]) AND ('
     '"Mycobacterium bovis"[tiab] OR "bovine tuberculosis"[tiab] OR zoonotic[tiab])'),
    ("12", "wildlife_occupational",
     '("wildlife reservoir"[tiab] OR badger[tiab] OR "white-tailed deer"[tiab] OR '
     '"abattoir worker"[tiab] OR slaughterhouse[tiab] OR "occupational tuberculosis"[tiab]) AND ('
     '"Mycobacterium bovis"[tiab] OR "Mycobacterium orygis"[tiab] OR zoonotic[tiab])'),
]

# Europe PMC equivalent (different syntax: simpler boolean, supports OPEN_ACCESS:Y)
def _epmc_query_for(node_id: str, slug: str, pubmed_q: str) -> str:
    # strip pubmed [tiab] tags, keep the boolean expressions
    q = re.sub(r"\[tiab\]", "", pubmed_q)
    return f"({q}) AND OPEN_ACCESS:Y AND FIRST_PDATE:[{DATE_FROM}-01-01 TO {DATE_TO}-12-31]"


def http_json(url: str, retries: int = 3, sleep: float = 0.4):
    last = None
    for k in range(retries):
        try:
            req = Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
            with urlopen(req, timeout=45) as r:
                return json.loads(r.read().decode("utf-8", errors="replace"))
        except (HTTPError, URLError, TimeoutError) as e:
            last = e
            time.sleep(sleep * (2 ** k))
    raise last


# ------------------------------------------------------------------ PubMed
def fetch_pubmed_node(node_id: str, slug: str, query: str, retmax: int = 200) -> list[dict]:
    esearch = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
               f"?db=pubmed&term={quote(query)}&retmax={retmax}&retmode=json")
    js = http_json(esearch); pmids = js["esearchresult"]["idlist"]
    count = int(js["esearchresult"]["count"])
    print(f"[pubmed][{slug}] hits={count}; retrieving {len(pmids)}")
    rows: list[dict] = []
    for i in range(0, len(pmids), 100):
        chunk = pmids[i:i + 100]
        esum = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                f"?db=pubmed&id={','.join(chunk)}&retmode=json")
        sj = http_json(esum); time.sleep(0.34)
        efetch = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                  f"?db=pubmed&id={','.join(chunk)}&rettype=abstract&retmode=text")
        try:
            req = Request(efetch, headers={"User-Agent": UA})
            with urlopen(req, timeout=60) as r:
                blob = r.read().decode("utf-8", errors="replace")
        except Exception:
            blob = ""
        ab_map: dict[str, str] = {}
        for piece in re.split(r"\n\n(?=\d+\. )", blob):
            m = re.search(r"PMID:\s*(\d+)", piece)
            if m: ab_map[m.group(1)] = re.sub(r"\s+", " ", piece).strip()
        for pmid in chunk:
            d = sj.get("result", {}).get(pmid, {})
            if not d: continue
            doi = ""
            for art_id in d.get("articleids", []):
                if art_id.get("idtype") == "doi":
                    doi = art_id.get("value", ""); break
            rows.append({
                "source": "pubmed", "node_id": node_id, "node_slug": slug,
                "pmid": pmid, "doi": doi.lower(),
                "title": d.get("title", "").strip(),
                "year": (d.get("pubdate", "") or "")[:4],
                "journal": d.get("fulljournalname", "") or d.get("source", ""),
                "abstract": ab_map.get(pmid, ""),
            })
        time.sleep(0.34)
    return rows


# --------------------------------------------------------------- Europe PMC
def fetch_europepmc_node(node_id: str, slug: str, query: str, page_size: int = 100,
                          max_pages: int = 2) -> list[dict]:
    rows: list[dict] = []
    cursor = "*"
    page = 0
    while page < max_pages:
        url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search"
               f"?query={quote(query)}&format=json&pageSize={page_size}&cursorMark={quote(cursor)}")
        js = http_json(url); time.sleep(0.4)
        hits = js.get("resultList", {}).get("result", [])
        if not hits: break
        for h in hits:
            rows.append({
                "source": "europepmc", "node_id": node_id, "node_slug": slug,
                "pmid": h.get("pmid", "") or "", "doi": (h.get("doi", "") or "").lower(),
                "title": (h.get("title", "") or "").strip(),
                "year": (h.get("firstPublicationDate", "") or "")[:4],
                "journal": h.get("journalTitle", "") or "",
                "abstract": (h.get("abstractText", "") or "").replace("\n", " ").strip(),
            })
        cursor = js.get("nextCursorMark", cursor)
        page += 1
        if len(hits) < page_size: break
    print(f"[europepmc][{slug}] retrieved {len(rows)} OA records")
    return rows


# --------------------------------------------------------------- deduplication
def _norm_title(s: str) -> str:
    s = re.sub(r"[^a-z0-9 ]", " ", s.lower())
    return re.sub(r"\s+", " ", s).strip()


def deduplicate(rows: list[dict]) -> list[dict]:
    seen_doi: set[str] = set(); seen_title: set[str] = set(); out: list[dict] = []
    for r in rows:
        doi = (r.get("doi") or "").strip().lower()
        title = _norm_title(r.get("title", ""))
        if doi and doi in seen_doi: continue
        # key on (normalised title + year) to absorb DOI-less duplicates
        ty = f"{title}|{r.get('year', '')}"
        if ty in seen_title and not doi: continue
        if doi: seen_doi.add(doi)
        if ty:  seen_title.add(ty)
        out.append(r)
    return out


# ------------------------------------------------------------------- main
def main() -> None:
    log_path = LOGS / "tb_interventions_harvest.log"
    log = log_path.open("w", encoding="utf-8")
    log.write(f"TB intervention-focused OA harvest, window {DATE_FROM}-{DATE_TO}\n\n")

    all_pubmed: list[dict] = []; all_epmc: list[dict] = []
    for node_id, slug, q_pubmed in TB_QUERIES_PUBMED:
        try:
            pm_rows = fetch_pubmed_node(node_id, slug, q_pubmed)
        except Exception as e:
            print(f"[pubmed][{slug}] ERROR: {e}"); pm_rows = []
        all_pubmed.extend(pm_rows)
        log.write(f"[pubmed][{slug}] {len(pm_rows)} records retrieved\n")

        q_epmc = _epmc_query_for(node_id, slug, q_pubmed)
        try:
            ep_rows = fetch_europepmc_node(node_id, slug, q_epmc)
        except Exception as e:
            print(f"[europepmc][{slug}] ERROR: {e}"); ep_rows = []
        all_epmc.extend(ep_rows)
        log.write(f"[europepmc][{slug}] {len(ep_rows)} records retrieved\n")

    # write per-source raw CSVs
    def _dump(path: Path, rows: list[dict]) -> None:
        if not rows: return
        cols = ["source", "node_id", "node_slug", "pmid", "doi", "title", "year", "journal", "abstract"]
        with path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)

    _dump(RAW / "tb_interventions_pubmed.csv", all_pubmed)
    _dump(RAW / "tb_interventions_europepmc.csv", all_epmc)

    # union + dedup
    unified = deduplicate(all_pubmed + all_epmc)
    log.write(f"\nRAW union: {len(all_pubmed) + len(all_epmc)}\n")
    log.write(f"AFTER DEDUP: {len(unified)}\n\n")

    # node-level counts in the deduplicated set
    from collections import Counter
    node_counter = Counter(r["node_slug"] for r in unified)
    log.write("PER-NODE COUNTS (deduplicated, may appear under multiple nodes if a record matched several queries on the source side):\n")
    for slug, n in sorted(node_counter.items(), key=lambda kv: -kv[1]):
        log.write(f"  {slug:<24} {n}\n")

    # write the unified CSV with extra columns expected by build_tb_bopas.py
    # (it expects 'title' and 'abstract_present' — we map abstract -> abstract_present)
    cols = ["source", "node_id", "node_slug", "pmid", "doi", "title", "year", "journal",
            "abstract_present"]
    out_csv = PROC / "tb_interventions_unified.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(cols)
        for r in unified:
            w.writerow([r.get("source", ""), r.get("node_id", ""), r.get("node_slug", ""),
                        r.get("pmid", ""), r.get("doi", ""), r.get("title", ""),
                        r.get("year", ""), r.get("journal", ""), r.get("abstract", "")])

    log.close()
    print(f"\nOK")
    print(f"  raw pubmed:    {len(all_pubmed)} -> {RAW/'tb_interventions_pubmed.csv'}")
    print(f"  raw epmc:      {len(all_epmc)} -> {RAW/'tb_interventions_europepmc.csv'}")
    print(f"  deduplicated:  {len(unified)} -> {out_csv}")
    print(f"  log:           {log_path}")


if __name__ == "__main__":
    main()
