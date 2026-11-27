#!/usr/bin/env python3
"""Evidence-map analysis of the TB intervention corpus (record-level AI codes with human validation).

Reads the repository copies of the coded corpus:
  data/full_bopas/pairs_no_abstract.jsonl   one line per record-intervention pair (abstract text removed;
                                            its SHA-256 and length kept)
  data/full_bopas/codes/*.codes.jsonl       code layers, applied in the order batch_* -> recode_* -> zscope_*
                                            -> ztier_* -> zz* (a later layer overrides an earlier one)
  data/full_bopas/spot_check_codes.csv      60-record human spot check (independent codes, agreed, adjudicated)
  validation/validation_codes_final.csv     earlier 150-pair validation (reported for reference only)
and writes outputs/full_bopas/results.json, readiness.csv and agreement.json.

Two fields need the abstract text: quotes_verbatim (every supporting quote found in title + abstract) and
off_topic_terms (term counts in off-topic records). Without abstracts they are written as null and a message
is printed. To fill them, rebuild the abstracts (code/full_bopas/prepare_batches.py on data/
tb_interventions_unified.csv, then fetch_abstracts.py and merge_fetched.py, which call NCBI E-utilities) and pass
the resulting pairs_all.jsonl with --pairs-with-abstracts PATH; each abstract is checked against the stored
SHA-256 before use.

Tier rules, first match wins: (a) no favourable trial and at least one null trial -> evaluate before (further)
scale-up; (b) <= 5 evaluations -> research gap; ready to scale with monitoring
(>= 10 SR/MA or RCT and comparative-design lower bound >= 0.5); implement and evaluate in India
(>= 10 SR/MA, RCT or non-randomised comparative evaluations and bound >= 0.5); otherwise evaluate before
(further) scale-up. Sensitivity analyses: all-design bound, trial-only bound for 'ready', count thresholds
5 and 20, and the tiers without the null-trial rule.

Usage: python -I code/full_bopas/paper/analysis.py [--pairs-with-abstracts PATH]
"""
from __future__ import annotations

import csv
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

from scipy.stats import beta as _beta

import argparse
import hashlib

REPO = Path(__file__).resolve().parents[3]
DATA = REPO / "data" / "full_bopas"
BUILD = REPO / "outputs" / "full_bopas"
BUILD.mkdir(parents=True, exist_ok=True)
ABSTRACTS: dict | None = None      # pair_id -> abstract text, only when --pairs-with-abstracts is given

FAV, NULL, UNC, NOT = ("favours intervention", "no effect or favours control", "unclear",
                       "not an evaluation of this intervention")
SR, RCT = "systematic review/meta-analysis", "randomised trial"
NRS = "non-randomised comparative (cohort, quasi-experimental, before-after, interrupted time series)"
ACC, ECON = "diagnostic accuracy", "economic evaluation"
PROG, MOD, OTHER = "implementation/programme report (no comparator)", "modelling", "other"
STRONG = (SR, RCT)
OUTCOME_COMP = (SR, RCT, NRS)                    # comparative designs measuring effects on outcomes
COMPARATIVE = (SR, RCT, NRS, ACC, ECON)          # pipeline's strict set: comparator or reference standard
LABELS = {
    "1": "Shorter all-oral regimens for RR-TB",
    "2": "Rapid molecular diagnostics (Xpert, Truenat)",
    "3": "Sequencing for drug-susceptibility prediction",
    "4": "Active case finding",
    "5": "TB preventive treatment",
    "6": "Airborne infection prevention and control",
    "7": "Nutritional support",
    "8": "Private-sector engagement",
    "9": "Integrated TB-HIV-diabetes care",
    "10": "Bovine TB control",
    "11": "Milk pasteurisation and dairy-chain controls",
    "12": "Wildlife surveillance and occupational protection",
}
ONE_HEALTH = {"10", "11", "12"}

# ---- readiness rules (see the module docstring) ----
GAP_MAX_EVAL = 5        # research gap: at most 5 evaluations
STRONG_MIN = 10         # ready: at least 10 SR/MA or RCT evaluations ...
LOWER_MIN = 0.5         # ... and lower credible bound of the favourable share >= 0.5
OUTCOME_MIN = 10        # implement and evaluate: at least 10 outcome-comparative evaluations (+ bound >= 0.5)
INDIA_MIN = 10          # India evidence flagged scarce below 10 India-setting evaluations
TIERS = ("Ready to scale with monitoring", "Implement and evaluate in India",
         "Evaluate before (further) scale-up", "Research gap")


def lower_bound(k: int, n: int) -> float:
    """5th percentile of Beta(k+1, n-k+1) (uniform prior); 0 when n = 0."""
    return 0.0 if n == 0 else float(_beta.ppf(0.05, k + 1, n - k + 1))


def wilson(k: int, n: int, z: float = 1.959964):
    if n == 0:
        return (float("nan"),) * 3
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, c - h, c + h


def kappa(a, b):
    n = len(a)
    po = sum(x == y for x, y in zip(a, b)) / n
    ca, cb = Counter(a), Counter(b)
    pe = sum(ca[k] * cb[k] for k in set(a) | set(b)) / (n * n)
    return (po - pe) / (1 - pe), po


def load():
    pairs = {p["pair_id"]: p for p in map(json.loads, (DATA / "pairs_no_abstract.jsonl").read_text(encoding="utf-8").splitlines())}
    for pid, p in pairs.items():
        p["abstract"] = ABSTRACTS.get(pid, "") if ABSTRACTS is not None else None
    codes, src, first = {}, {}, {}
    files_by_round = Counter()
    for pat in ("batch_*", "recode_*", "zscope_*", "ztier_*", "zz*"):
        for f in sorted((DATA / "codes").glob(pat + ".codes.jsonl")):
            files_by_round[pat] += 1
            for line in f.read_text(encoding="utf-8").splitlines():
                if line.strip():
                    c = json.loads(line)
                    pid = int(c["pair_id"])
                    codes[pid] = c
                    src[pid] = f.name
                    if pat in ("batch_*", "recode_*"):
                        first[pid] = c["code"]
    return pairs, codes, src, first, files_by_round


def norm(s: str) -> str:
    return " ".join(s.lower().split())


def tier_for(r: dict, lb_key: str = "L_strict", strong_min=STRONG_MIN, outcome_min=OUTCOME_MIN,
             ready_key: str | None = None, null_rule: bool = True) -> str:
    """Readiness rules, applied in order (null-trial rule first). ready_key: bound used for the 'ready' test only
    (sensitivity: the bound from SR/MA and RCT evaluations alone)."""
    if null_rule and r["rct_favours"] == 0 and r["rct_null"] >= 1:     # rule 1: every directional trial null
        return TIERS[2]
    if r["evaluations"] <= GAP_MAX_EVAL:
        return TIERS[3]
    if r["strong"] >= strong_min and r[ready_key or lb_key] >= LOWER_MIN:
        return TIERS[0]
    if r["outcome_comparative"] >= outcome_min and r[lb_key] >= LOWER_MIN:
        return TIERS[1]
    return TIERS[2]


def main() -> None:
    pairs, codes, src, first, files_by_round = load()
    res: dict = {"files_by_round": dict(files_by_round)}
    res["pairs"] = len(pairs)
    res["coded"] = len(codes)
    coded_all = [(pairs[i], c, src[i], first.get(i)) for i, c in sorted(codes.items())]
    seen, coded, dups = set(), [], 0
    for p, c, s, f in coded_all:
        key = (p["node_id"], "".join(ch for ch in p["title"].lower() if ch.isalnum())[:90])
        if key in seen:
            dups += 1
            continue
        seen.add(key)
        coded.append((p, c, s, f))
    N = len(coded)
    res.update(merged=dups, analysed=N)
    res["sources_pairs"] = dict(Counter(p["source"] for p in pairs.values()))
    res["abstract_source_pairs"] = dict(Counter(p.get("abstract_source", "original harvest") for p in pairs.values()))
    res["title_only_pairs"] = sum(p["title_only"] for p in pairs.values())
    res["title_only_analysed"] = sum(p["title_only"] for p, *_ in coded)
    if ABSTRACTS is None:
        res["quotes_verbatim"] = None
        print("quotes_verbatim and off_topic_terms need abstract text: written as null "
              "(pass --pairs-with-abstracts PATH to compute them)")
    else:
        res["quotes_verbatim"] = sum(bool(norm(c.get("quote", ""))) and norm(c.get("quote", "")) in norm(p["title"] + " " + p["abstract"])
                                     for p, c, *_ in coded)
    # per-node retrieval by source
    res["records_per_node"] = {nid: dict(Counter(p["source"] for p in pairs.values() if p["node_id"] == nid)) for nid in LABELS}

    # years
    yrs = [int(p["year"]) for p, *_ in coded if str(p["year"]).isdigit()]
    yc = Counter(yrs)
    res["years"] = {"n_with_year": len(yrs), "min": min(yrs), "max": max(yrs), "median": statistics.median(yrs),
                    "share_2025_2026": round((yc[2025] + yc[2026]) / len(yrs), 3), "n_2025_2026": yc[2025] + yc[2026],
                    "share_2023_2026": round(sum(yc[y] for y in (2023, 2024, 2025, 2026)) / len(yrs), 3),
                    "before_2017": sum(v for y, v in yc.items() if y < 2017),
                    "by_year": {str(y): yc[y] for y in sorted(yc)}}
    by_src_year = {}
    for s_ in ("pubmed", "europepmc"):
        ys = [int(p["year"]) for p, *_ in coded if p["source"] == s_ and str(p["year"]).isdigit()]
        by_src_year[s_] = {"n": len(ys), "median": statistics.median(ys),
                           "share_2025_2026": round(sum(y >= 2025 for y in ys) / len(ys), 3)}
    res["years_by_source"] = by_src_year
    ev_years = [int(p["year"]) for p, c, *_ in coded if c["code"] != NOT and str(p["year"]).isdigit()]
    res["years"]["evaluations_share_2025_2026"] = round(sum(y >= 2025 for y in ev_years) / len(ev_years), 3)

    # model provenance of final codes
    res["final_codes_from_sonnet"] = sum("sonnet" in s for _, _, s, _ in coded)
    res["final_codes_from_sonnet_by_node"] = dict(Counter(p["node_id"] for p, c, s, _ in coded if "sonnet" in s))
    res["final_code_round"] = dict(Counter(s.split("_")[0] for _, _, s, _ in coded))
    changed = [(f, c["code"]) for p, c, s, f in coded if f is not None and f != c["code"]]
    res["first_round_vs_final_changed"] = len(changed)
    res["first_round_vs_final_n"] = sum(f is not None for *_, f in coded)
    res["first_round_to_final_transitions"] = Counter(changed).most_common(8)
    first_eval = sum(f is not None and f != NOT for *_, f in coded)
    res["first_round_evaluations"] = first_eval

    # overall
    cc = Counter(c["code"] for _, c, *_ in coded)
    ev = N - cc[NOT]
    tier2 = sum(c["code"] == NOT and c.get("tier2") is True for _, c, *_ in coded)
    off = cc[NOT] - tier2
    res["overall"] = {"evaluations": ev, "eval_share_ci": [round(x, 3) for x in wilson(ev, N)],
                      "favours": cc[FAV], "no_effect": cc[NULL], "unclear": cc[UNC],
                      "tier2": tier2, "tier2_share": round(tier2 / N, 3), "off_topic": off, "off_share": round(off / N, 3),
                      "fav_share_directional": [round(x, 3) for x in wilson(cc[FAV], cc[FAV] + cc[NULL])],
                      "eval_plus_tier2_share": round((ev + tier2) / N, 3),
                      "tier2_missing_flag": sum(c["code"] == NOT and c.get("tier2") not in (True, False) for _, c, *_ in coded)}
    # direction by design (all interventions)
    dd = {}
    for d in (SR, RCT, NRS, ACC, ECON, MOD, PROG, OTHER):
        items = [c for _, c, *_ in coded if c["code"] != NOT and c.get("design") == d]
        k = sum(c["code"] == FAV for c in items)
        n0 = sum(c["code"] == NULL for c in items)
        u = sum(c["code"] == UNC for c in items)
        p, lo, hi = wilson(k, k + n0)
        dd[d] = {"evaluations": len(items), "favours": k, "no_effect": n0, "unclear": u, "directional": k + n0,
                 "fav_share": round(p, 3), "lo": round(lo, 3), "hi": round(hi, 3)}
    res["direction_by_design"] = dd
    res["designs_other_than_listed"] = dict(Counter(c.get("design") for _, c, *_ in coded if c["code"] != NOT
                                                    and c.get("design") not in dd))
    strong_all = [c for _, c, *_ in coded if c["code"] != NOT and c.get("design") in STRONG]
    k_s = sum(c["code"] == FAV for c in strong_all); n_s = sum(c["code"] == NULL for c in strong_all)
    nonc = [c for _, c, *_ in coded if c["code"] != NOT and c.get("design") in (PROG,)]
    res["fav_share_strong"] = [round(x, 3) for x in wilson(k_s, k_s + n_s)] + [k_s, n_s]
    k_p = sum(c["code"] == FAV for c in nonc); n_p = sum(c["code"] == NULL for c in nonc)
    res["fav_share_programme"] = [round(x, 3) for x in wilson(k_p, k_p + n_p)] + [k_p, n_p]

    # settings
    res["setting_evaluations"] = dict(Counter(c.get("setting") for _, c, *_ in coded if c["code"] != NOT))
    india_ev = sum(c["code"] != NOT and c.get("setting") == "India" for _, c, *_ in coded)
    res["india_evaluations"] = india_ev
    res["india_share_of_evaluations"] = round(india_ev / ev, 3)
    india_strong = [c for _, c, *_ in coded if c["code"] != NOT and c.get("setting") == "India" and c.get("design") in STRONG]
    res["india_strong"] = len(india_strong)
    res["india_strong_by_node"] = dict(Counter(p["node_id"] for p, c, *_ in coded if c["code"] != NOT
                                               and c.get("setting") == "India" and c.get("design") in STRONG))
    res["india_designs"] = dict(Counter(c.get("design") for _, c, *_ in coded if c["code"] != NOT and c.get("setting") == "India"))

    # per intervention
    by = defaultdict(list)
    for p, c, *_ in coded:
        by[p["node_id"]].append((p, c))
    rows = []
    for nid, lbl in LABELS.items():
        items = by[nid]
        n = len(items)
        cn = Counter(c["code"] for _, c in items)
        e = n - cn[NOT]
        evs = [c for _, c in items if c["code"] != NOT]
        dz = Counter(c.get("design") for c in evs)
        st = [c for c in evs if c.get("design") in COMPARATIVE]
        sg = [c for c in evs if c.get("design") in STRONG]
        k_st = sum(c["code"] == FAV for c in st); nd_st = sum(c["code"] in (FAV, NULL) for c in st)
        k_sg = sum(c["code"] == FAV for c in sg); nd_sg = sum(c["code"] in (FAV, NULL) for c in sg)
        share, lo, hi = wilson(e, n)
        t2 = sum(c["code"] == NOT and c.get("tier2") is True for _, c in items)
        r = {"node_id": nid, "intervention": lbl, "group": "One Health" if nid in ONE_HEALTH else "Human health",
             "records": n, "evaluations": e, "eval_share": round(share, 3), "eval_lo": round(lo, 3), "eval_hi": round(hi, 3),
             "tier2": t2, "off_topic": n - e - t2,
             "favours": cn[FAV], "no_effect": cn[NULL], "unclear": cn[UNC],
             "fav_share": round(cn[FAV] / (cn[FAV] + cn[NULL]), 3) if cn[FAV] + cn[NULL] else None,
             "L_all": round(lower_bound(cn[FAV], cn[FAV] + cn[NULL]), 3),
             "comparative": len(st), "L_strict": round(lower_bound(k_st, nd_st), 3),
             "rct_favours": sum(c["code"] == FAV for c in evs if c.get("design") == RCT),
             "rct_null": sum(c["code"] == NULL for c in evs if c.get("design") == RCT),
             "strong": len(sg), "SR": dz[SR], "RCT": dz[RCT], "NRS": dz[NRS], "accuracy": dz[ACC],
             "economic": dz[ECON], "modelling": dz[MOD], "programme": dz[PROG], "other_design": dz[OTHER],
             "outcome_comparative": dz[SR] + dz[RCT] + dz[NRS],
             "L_strong": round(lower_bound(k_sg, nd_sg), 3), "strong_directional": nd_sg,
             "k_strict": k_st, "n_strict_directional": nd_st, "k_strong": k_sg,
             "settings": "{} / {} / {} / {} / {}".format(*[sum(c.get("setting") == x for c in evs) for x in (
                 "India", "other low/middle-income", "high-income", "multi-country", "not stated")]),
             "fav_share_strict": round(k_st / nd_st, 3) if nd_st else None,
             "india": sum(c.get("setting") == "India" for c in evs),
             "title_only_evals": sum(p["title_only"] and c["code"] != NOT for p, c in items)}
        r["india_flag"] = "scarce" if r["india"] < INDIA_MIN else "available"
        rows.append(r)
    # readiness and sensitivity
    for r in rows:
        r["tier"] = tier_for(r)
        r["tier_all_designs"] = tier_for(r, "L_all")
        r["tier_strong_bound"] = tier_for(r, ready_key="L_strong")
        r["tier_thresholds_5"] = tier_for(r, strong_min=5, outcome_min=5)
        r["tier_thresholds_20"] = tier_for(r, strong_min=20, outcome_min=20)
        r["tier_without_null_trial_rule"] = tier_for(r, null_rule=False)
    res["interventions"] = rows
    res["tier_changes"] = {s: [(r["intervention"], r["tier"], r[s]) for r in rows if r[s] != r["tier"]]
                           for s in ("tier_all_designs", "tier_strong_bound", "tier_thresholds_5", "tier_thresholds_20")}
    res["tier_counts"] = dict(Counter(r["tier"] for r in rows))
    res["null_trial_rule_moved"] = [(r["intervention"], r["tier_without_null_trial_rule"], r["tier"], f"{r['rct_null']} null / {r['rct_favours']} favourable trials")
                                    for r in rows if r["tier_without_null_trial_rule"] != r["tier"]]
    res["null_trial_rule_applies"] = [r["intervention"] for r in rows if r["rct_favours"] == 0 and r["rct_null"] >= 1]
    sens = ("tier_all_designs", "tier_strong_bound", "tier_thresholds_5", "tier_thresholds_20")
    res["tier_changed_any"] = sum(any(r[k] != r["tier"] for k in sens) for r in rows)
    res["strong_top3"] = sum(r["strong"] for r in rows if r["node_id"] in ("1", "4", "5"))
    res["strong_total"] = sum(r["strong"] for r in rows)
    res["india_top3"] = sum(r["india"] for r in rows if r["node_id"] in ("2", "7", "8"))
    res["derived_shares"] = {
        "strong_top3_of_strong": round(res["strong_top3"] / res["strong_total"], 3),
        "sonnet_of_analysed": round(res["final_codes_from_sonnet"] / N, 3),
        "changed_of_first_round": round(res["first_round_vs_final_changed"] / res["first_round_vs_final_n"], 3),
        "rct_share_of_evaluations": round(dd[RCT]["evaluations"] / ev, 3),
        "rct_nulls_of_all_nulls": f"{dd[RCT]['no_effect']}/{cc[NULL]}",
        "pasteurisation_eval_share": [r["eval_share"] for r in rows if r["node_id"] == "11"][0],
        "accuracy_xpert": [f"{r['accuracy']}/{r['evaluations']}" for r in rows if r["node_id"] == "2"][0],
        "accuracy_seq": [f"{r['accuracy']}/{r['evaluations']}" for r in rows if r["node_id"] == "3"][0],
        "bovine_high_income": int([r["settings"] for r in rows if r["node_id"] == "10"][0].split(" / ")[2]),
        "bovine_trials_directional_fav": [f"{r['k_strong']}/{r['strong_directional']}" for r in rows if r["node_id"] == "10"][0],
        "india_one_health": sum(r["india"] for r in rows if r["group"] == "One Health"),
        "india_strong_one_health": sum(v for k, v in res["india_strong_by_node"].items() if k in ONE_HEALTH),
    }

    # group totals
    grp = {}
    for g in ("Human health", "One Health"):
        rr = [r for r in rows if r["group"] == g]
        n = sum(r["records"] for r in rr); e = sum(r["evaluations"] for r in rr)
        grp[g] = {"records": n, "evaluations": e, "share": [round(x, 3) for x in wilson(e, n)],
                  "tier2": sum(r["tier2"] for r in rr), "off_topic": sum(r["off_topic"] for r in rr),
                  "india": sum(r["india"] for r in rr), "strong": sum(r["strong"] for r in rr)}
    oh_excl_bovine = [r for r in rows if r["node_id"] in ("11", "12")]
    n = sum(r["records"] for r in oh_excl_bovine); e = sum(r["evaluations"] for r in oh_excl_bovine)
    grp["Pasteurisation + wildlife/occupational"] = {"records": n, "evaluations": e, "share": [round(x, 3) for x in wilson(e, n)],
                                                     "tier2": sum(r["tier2"] for r in oh_excl_bovine)}
    res["groups"] = grp

    # what fills the off-topic share (term counts in off-topic records)
    def off_terms(nid, terms):
        off = [(p, c) for p, c in by[nid] if c["code"] == NOT and c.get("tier2") is not True]
        return {"off_topic": len(off), **{t: sum(t in (p["title"] + " " + p["abstract"]).lower() for p, _ in off) for t in terms}}
    res["off_topic_terms"] = None if ABSTRACTS is None else {
        "6": off_terms("6", ["mechanical ventilation", "covid"]),
        "11": off_terms("11", ["brucell", "raw milk"]),
        "12": off_terms("12", ["brucell", "covid"])}
    # cross-check with pipeline results
    pipe = {r["node_id"]: r for r in csv.DictReader(open(BUILD / "evidence_by_intervention.csv", encoding="utf-8"))}
    mism = []
    for r in rows:
        p = pipe[r["node_id"]]
        for a, b in (("records", "pairs_coded"), ("evaluations", "evaluations"), ("favours", "favours"),
                     ("no_effect", "no_effect_or_control"), ("unclear", "unclear"), ("comparative", "comparative_evaluations"),
                     ("strong", "SR_or_RCT_evaluations"), ("india", "india_evaluations"), ("tier2", "tier2_relevant"),
                     ("off_topic", "off_topic")):
            if int(p[b]) != r[a]:
                mism.append((r["node_id"], a, r[a], p[b]))
        for a, b in (("L_all", "L_ai"), ("L_strict", "L_ai_strict"), ("L_strong", "L_ai_SR_RCT")):
            if abs(float(p[b]) - r[a]) > 0.0015:
                mism.append((r["node_id"], a, r[a], p[b]))
    res["pipeline_crosscheck_mismatches"] = mism

    # ---- spot check (repository CSV; no names) ----
    srows = list(csv.DictReader(open(DATA / "spot_check_codes.csv", encoding="utf-8")))
    for r in srows:
        r["intervention"] = pairs[int(r["pair_id"])]["intervention"]
    adj = {int(r["item"]): (r["final_code"], r["reason"]) for r in srows if r["reason"]}
    key = {int(r["item"]): r["ai_code"] for r in srows}
    A = [str(r["coder_A"]).strip() for r in srows]
    Bc = [str(r["coder_B"]).strip() for r in srows]
    final_h = [r["final_code"].strip() for r in srows]
    ai = [codes[int(r["pair_id"])]["code"] for r in srows]
    t2 = [codes[int(r["pair_id"])].get("tier2") is True for r in srows]
    kab, pab = kappa(A, Bc)
    k4, p4 = kappa(ai, final_h)
    kb, pb = kappa([c != NOT for c in ai], [c != NOT for c in final_h])
    rel = [i for i, (a, hc, t) in enumerate(zip(ai, final_h, t2)) if a == NOT and t and hc != NOT]
    keep = [i for i in range(len(srows)) if i not in rel]
    kx, px = kappa([ai[i] for i in keep], [final_h[i] for i in keep])
    ai_at_sampling = [key[int(r["item"])] for r in srows]
    kfirst, pfirst = kappa(ai_at_sampling, final_h)
    agreed = [str(r["agreed_code"]).strip() for r in srows]
    kag, pag = kappa(ai, agreed)
    kag0, _ = kappa(ai_at_sampling, agreed)
    kagb, _ = kappa([c != NOT for c in ai], [c != NOT for c in agreed])
    res["spot_check"] = {
        "ai_final_vs_agreed_before_adjudication_kappa4": round(kag, 2),
        "ai_final_vs_agreed_before_adjudication_kappa_eval": round(kagb, 2),
        "ai_at_sampling_vs_agreed_kappa4": round(kag0, 2),
        "adjudication_changed_authors_code": sum(a_ != f_ for a_, f_ in zip(agreed, final_h)),
        "n": len(srows), "A_vs_B_kappa": round(kab, 2), "A_vs_B_agree": f"{sum(x == y for x, y in zip(A, Bc))}/{len(A)}",
        "adjudicated_items": len(adj), "why_counts": dict(Counter(v[1] for v in adj.values())),
        "ai_final_vs_human_kappa4": round(k4, 2), "agree4": f"{sum(x == y for x, y in zip(ai, final_h))}/{len(ai)}",
        "ai_final_vs_human_kappa_eval": round(kb, 2),
        "agree_eval": f"{sum((x != NOT) == (y != NOT) for x, y in zip(ai, final_h))}/{len(ai)}",
        "ai_tier2_human_eval": len(rel), "kappa4_excluding": round(kx, 2), "n_excluding": len(keep),
        "agree4_excluding": f"{sum(ai[i] == final_h[i] for i in keep)}/{len(keep)}",
        "ai_at_sampling_vs_human_kappa4": round(kfirst, 2),
        "sample_ai_at_sampling": dict(Counter(ai_at_sampling)), "sample_ai_final": dict(Counter(ai)),
        "human_final": dict(Counter(final_h)),
        "ai_not_eval_sampled": sum(a == NOT for a in ai),
        "ai_not_eval_tier2_sampled": sum(a == NOT and t for a, t in zip(ai, t2)),
        "ai_not_eval_offtopic_sampled": sum(a == NOT and not t for a, t in zip(ai, t2)),
        "offtopic_human_eval": sum(a == NOT and not t and hc != NOT for a, t, hc in zip(ai, t2, final_h)),
        "ai_eval_human_not": sum(a != NOT and hc == NOT for a, hc in zip(ai, final_h)),
        "confusion": [[list(k), v] for k, v in Counter(zip(ai, final_h)).most_common()],
        "direction_agree_both_eval": f"{sum(a == hc for a, hc in zip(ai, final_h) if a != NOT and hc != NOT)}/"
                                     f"{sum(a != NOT and hc != NOT for a, hc in zip(ai, final_h))}",
        "sample_by_node": dict(Counter(str(r["intervention"]) for r in srows)),
        "sample_title_only": sum(pairs[int(r["pair_id"])]["title_only"] for r in srows),
    }
    # ---- earlier 150-pair validation (for reference) ----
    v = list(csv.DictReader(open(REPO / "validation" / "validation_codes_final.csv", encoding="utf-8")))
    res["pilot150"] = {
        "n": len(v),
        "human_human_kappa": round(kappa([r["coder_A_independent"] for r in v], [r["coder_B_independent"] for r in v])[0], 2),
        "ai_vs_consensus_before": round(kappa([r["ai_blind_recode"] for r in v], [r["consensus_before_ai_check"] for r in v])[0], 2),
        "ai_vs_final_consensus": round(kappa([r["ai_blind_recode"] for r in v], [r["consensus_code"] for r in v])[0], 2),
    }
    (BUILD / "results.json").write_text(json.dumps(res, indent=1, ensure_ascii=False, default=str), encoding="utf-8")
    (BUILD / "agreement.json").write_text(json.dumps({"spot_check_60": res["spot_check"], "validation_150": res["pilot150"]},
                                                     indent=1, ensure_ascii=False), encoding="utf-8")
    with (BUILD / "readiness.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(json.dumps({k: res[k] for k in ("pairs", "merged", "analysed", "overall", "groups", "tier_counts", "tier_changes",
                                          "pipeline_crosscheck_mismatches", "spot_check", "pilot150", "years",
                                          "years_by_source")}, indent=1, default=str))
    for r in rows:
        print(f"{r['node_id']:>2} {r['intervention'][:42]:42} n{r['records']:4} E{r['evaluations']:4} S{r['strong']:3} "
              f"O{r['outcome_comparative']:3} C{r['comparative']:3} Ls{r['L_strict']:.3f} La{r['L_all']:.3f} "
              f"Lsg{r['L_strong']:.3f}({r['strong_directional']}) I{r['india']:3} -> {r['tier']} | {r['tier_strong_bound']}")


def _load_abstracts(path: Path) -> dict:
    stored = {p["pair_id"]: p for p in map(json.loads, (DATA / "pairs_no_abstract.jsonl").read_text(encoding="utf-8").splitlines())}
    out, bad = {}, 0
    for p in map(json.loads, path.read_text(encoding="utf-8").splitlines()):
        ab = p.get("abstract", "") or ""
        want = stored[p["pair_id"]]["abstract_sha256"]
        got = hashlib.sha256(ab.encode("utf-8")).hexdigest() if ab else None
        if got != want:
            bad += 1
            continue
        out[p["pair_id"]] = ab
    print(f"abstracts: {len(out)} pairs match the stored SHA-256; {bad} do not (treated as missing)")
    if bad:
        raise SystemExit("abstracts differ from the coded text; quotes_verbatim would not be comparable")
    return out


if __name__ == "__main__":
    sys.dont_write_bytecode = True
    ap = argparse.ArgumentParser(description="Evidence-map analysis")
    ap.add_argument("--pairs-with-abstracts", type=Path, default=None)
    args = ap.parse_args()
    if args.pairs_with_abstracts:
        ABSTRACTS = _load_abstracts(args.pairs_with_abstracts)
    main()
