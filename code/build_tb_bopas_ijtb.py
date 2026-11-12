"""
build_tb_bopas_ijtb.py
======================

TB-BOPAS, rebuilt for the Indian Journal of Tuberculosis submission.

What changed relative to `code/build_tb_bopas.py` (and why):

  1. SURVEILLANCE DATA. The previous script carried five hard-coded
     `STAND_IN_SERIES` literals and never read the CSV it documented as their
     source (that CSV did not exist). This script reads
     `ntep_trajectories.csv`, in which every datapoint carries its numerator,
     denominator, WHO field names and source dataset. Values are the WHO
     Global TB Programme country data for India (extranet.who.int/tme,
     downloaded 8 September 2026).

  2. DIRECTION-OF-EFFECT CODING. The previous classifier tested favour cues
     before contrary cues (so "reduced notifications but not significantly"
     scored as favourable) and treated `prevent|prevention` as an evaluative
     favour cue, although "prevention" is the *name* of two intervention
     nodes. Both defects are corrected here; the uncorrected coding is
     retained as `code_corpus_uncorrected` for the sensitivity analysis.

  3. UNCERTAINTY NORMALISATION. The half-width floor of 0.5 is retained but
     is now explicit in the formula, reported, and accompanied by an
     unfloored sensitivity analysis. In the previous version the floor bound
     for all five series, so the normalisation never operated.

  4. REPORTING DISCONTINUITIES. Series with a documented reporting break are
     fitted both naively over the full period and segmented at the break.
     For the primary analysis the within-era slope is used.

Outputs (all under this directory's ../outputs/):
  ijtb_slopes.csv          per-series slopes, CIs, segmented fits, classification
  ijtb_node_counts.csv     per-node k/n and L, corrected and uncorrected
  ijtb_bopas_long.csv      full BOPAS grid
  ijtb_bopas_top.csv       ranked non-zero cells
  ijtb_sensitivity.csv     every sensitivity scenario reported in the paper
  ijtb_diagnostics.txt     corpus coverage and classifier-bias diagnostics
"""
from __future__ import annotations
import csv, math, random, re
from collections import Counter, defaultdict
from pathlib import Path

random.seed(20260908)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent                                     # repository root
OUT = ROOT / "outputs"; OUT.mkdir(exist_ok=True)
TRAJ_CSV = ROOT / "data" / "ntep_trajectories.csv"
CORPUS = ROOT / "data" / "tb_interventions_unified.csv"

# --------------------------------------------------------------------------
# Trajectory targets
# --------------------------------------------------------------------------
# sign = +1 when rising values mean a worsening burden, -1 when rising values
# mean an improving burden (treatment-success series).
# break_year = first year of the post-discontinuity era, where a reporting
# break is documented in the source data; None if no break.
TARGETS = [
    ("RR_new_pct",       "Rifampicin resistance among new pulmonary TB cases (%)", +1, 2021),
    ("RRMDR_notif_100k", "RR/MDR-TB notifications per 100 000 population",         +1, None),
    ("preXDR_fqr_pct",   "Fluoroquinolone-resistant RR-TB (pre-XDR) fraction (%)", +1, None),
    ("MDR_TSR_pct",      "MDR/RR-TB treatment success rate (%)",                   -1, None),
    ("XDR_TSR_pct",      "XDR-TB treatment success rate (%)",                      -1, None),
]
DESCRIPTIVE_ONLY = {"BDQ_res_pct"}      # too few reporting years to fit

HALF_WIDTH_FLOOR = 0.5   # see Methods: caps the influence of implausibly
                         # tight intervals on 6-9 annual aggregates


# --------------------------------------------------------------------------
# Stage 1: intervention-node coding of the open-access corpus
# --------------------------------------------------------------------------
NODE_LABELS = {
    "1":  "Shorter all-oral regimen rollout (BPaL/BPaL-M)",
    "2":  "Molecular diagnostic scale-up (Xpert MTB/RIF, Xpert MTB/XDR)",
    "3":  "Whole-genome sequencing for drug-susceptibility prediction",
    "4":  "Active case finding",
    "5":  "TB preventive therapy",
    "6":  "Airborne infection prevention and control",
    "7":  "Nutritional support",
    "8":  "Private-sector engagement",
    "9":  "Integrated TB-HIV-diabetes care",
    "10": "Bovine tuberculosis control",
    "11": "Milk pasteurisation / dairy food chain",
    "12": "Wildlife reservoir surveillance / occupational respiratory protection",
}

NODE_REGEX = {
    "1":  re.compile(r"\b(bpal|bedaquiline|pretomanid|nix[- ]?tb|zenix|practecal|all[- ]?oral.*regimen|shorter.*regimen|6[- ]?month.*regimen)\b", re.I),
    "2":  re.compile(r"\b(xpert|genexpert|cepheid|mtb/rif|mtb/xdr|line[- ]?probe|truenat|molecular diagnostic)\b", re.I),
    "3":  re.compile(r"\b(whole[- ]genome sequenc|wgs|genom.*surveillance|drug[- ]susceptibility prediction|catalogue of mutations)\b", re.I),
    "4":  re.compile(r"\b(active case[- ]finding|door[- ]to[- ]door|systematic screening|contact (investigation|tracing))\b", re.I),
    "5":  re.compile(r"\b(preventive therapy|tb preventive|tpt|3hp|1hp|6h|isoniazid prophylaxis|levofloxacin.*contact)\b", re.I),
    "6":  re.compile(r"\b(airborne (infection )?(prevention|control)|ipc|n95|respirator program|ventilation|ultraviolet germicidal|uvgi)\b", re.I),
    "7":  re.compile(r"\b(nutritional support|food support|nutrition.*tuberculosis|undernutrition.*tb)\b", re.I),
    "8":  re.compile(r"\b(private[- ]sector engagement|joint effort.*tb|jeet|informant network|notification scheme|voucher)\b", re.I),
    "9":  re.compile(r"\b(integrated tb[- ]?hiv|tb[- ]?diabetes|comorbidity care|joint tb hiv)\b", re.I),
    "10": re.compile(r"\b(bovine tuberculosis|m\.? bovis|test[- ]and[- ]cull|cattle.*tuberculosis|herd.*vaccination)\b", re.I),
    "11": re.compile(r"\b(milk pasteurisation|dairy.*food chain|raw milk|unpasteurised milk|pasteuri[sz]ation enforcement)\b", re.I),
    "12": re.compile(r"\b(wildlife reservoir|slaughterhouse.*protection|occupational.*tuberculosis|abattoir.*respiratory)\b", re.I),
}

# Evaluative cues. Note what is NOT here: prevent / prevention / prevented.
# Those words name interventions ("TB preventive therapy", "infection
# prevention and control") rather than reporting a favourable result.
RE_FAVOURS = re.compile(r"\b(reduce[ds]?|decreased?|lower(?:ed)?|improve[ds]?|effective|"
                        r"successful|favou?rable|protective)\b", re.I)
RE_FAVOURS_UNCORRECTED = re.compile(r"\b(reduce[ds]?|decreased?|lower(?:ed)?|improve[ds]?|effective|"
                                    r"successful|favou?rable|protective|prevent(?:ed|ion)?)\b", re.I)
RE_CONTRARY = re.compile(r"\b(increase[ds]?|rose|rising|no effect|did not differ|not significant|"
                         r"inconclusive|ineffective|fail(?:ed)?)\b", re.I)
RE_UNCLEAR = re.compile(r"\b(mixed (results|evidence)|heterogen|unclear|further (research|study) needed)\b", re.I)


def _read_corpus() -> list[dict[str, str]]:
    with CORPUS.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def _record_text(row: dict[str, str]) -> str:
    # NB: the column named 'abstract_present' holds abstract text, not a flag.
    return " ".join([row.get("title", ""), row.get("abstract_present", "")])


def code_corpus(rows, *, corrected: bool = True) -> dict[str, dict[str, int]]:
    """Per-node {k_favours, n_directional, n_total}.

    corrected=True  : contrary cues take precedence; nomenclature words excluded.
    corrected=False : reproduces the original (favour-first, nomenclature-inclusive)
                      classifier for the sensitivity analysis.
    """
    favour_re = RE_FAVOURS if corrected else RE_FAVOURS_UNCORRECTED
    counts = {nid: Counter() for nid in NODE_LABELS}
    for row in rows:
        text = _record_text(row)
        for nid, rx in NODE_REGEX.items():
            if not rx.search(text):
                continue
            counts[nid]["n_total"] += 1
            if RE_UNCLEAR.search(text):
                continue
            if corrected:
                if RE_CONTRARY.search(text):            # contrary wins
                    counts[nid]["n_directional"] += 1
                elif favour_re.search(text):
                    counts[nid]["k_favours"] += 1
                    counts[nid]["n_directional"] += 1
            else:
                if favour_re.search(text):              # original precedence
                    counts[nid]["k_favours"] += 1
                    counts[nid]["n_directional"] += 1
                elif RE_CONTRARY.search(text):
                    counts[nid]["n_directional"] += 1
    return {nid: dict(c) for nid, c in counts.items()}


def corpus_diagnostics(rows) -> str:
    lines, matched, multi, titleonly = [], 0, 0, 0
    per_node_titleonly = Counter()
    for row in rows:
        text = _record_text(row)
        hits = [n for n, rx in NODE_REGEX.items() if rx.search(text)]
        if hits:
            matched += 1
            if not row.get("abstract_present", "").strip():
                titleonly += 1
                for n in hits:
                    per_node_titleonly[n] += 1
        if len(hits) > 1:
            multi += 1
    total = len(rows)
    lines.append(f"harvested records in corpus            : {total}")
    lines.append(f"records matching >=1 intervention node : {matched} ({matched/total:.1%})")
    lines.append(f"records matching no node (unused)      : {total-matched} ({(total-matched)/total:.1%})")
    lines.append(f"records matching >1 node               : {multi}")
    no_abs = sum(1 for r in rows if not r.get("abstract_present", "").strip())
    lines.append(f"records with no abstract text          : {no_abs} ({no_abs/total:.1%})")
    lines.append(f"coded records relying on title only    : {titleonly}")
    lines.append("")
    lines.append("title-only coded records per node:")
    for nid, lbl in NODE_LABELS.items():
        lines.append(f"  node {nid:>2}  {lbl[:56]:<58} {per_node_titleonly[nid]}")

    # classifier bias diagnostics, measured on the uncorrected classifier
    prevent = re.compile(r"\b(prevent(?:ed|ion)?)\b", re.I)
    kf = both = only_prevent = 0
    triggers = Counter()
    for row in rows:
        text = _record_text(row)
        if not any(rx.search(text) for rx in NODE_REGEX.values()):
            continue
        if RE_UNCLEAR.search(text):
            continue
        hits = RE_FAVOURS_UNCORRECTED.findall(text)
        if not hits:
            continue
        kf += 1
        words = {w.lower() for w in hits}
        for w in words:
            triggers[w] += 1
        if all(prevent.fullmatch(w) for w in words):
            only_prevent += 1
        if RE_CONTRARY.search(text):
            both += 1
    lines.append("")
    lines.append("bias in the original (uncorrected) direction-of-effect classifier:")
    lines.append(f"  records scored favours-intervention                    : {kf}")
    lines.append(f"  ... containing contrary cues too (resolved as favours) : {both} ({both/kf:.1%})")
    lines.append(f"  ... whose only favour cue is a prevent-family word     : {only_prevent} ({only_prevent/kf:.1%})")
    lines.append(f"  favour-cue frequency: {triggers.most_common()}")
    return "\n".join(lines)


# --------------------------------------------------------------------------
# Stage 2: trajectory fitting
# --------------------------------------------------------------------------
def _ols(xs, ys):
    n = len(xs); xb = sum(xs)/n; yb = sum(ys)/n
    sxx = sum((x-xb)**2 for x in xs)
    return sum((x-xb)*(y-yb) for x, y in zip(xs, ys))/sxx if sxx else 0.0


def _bootstrap_ci(xs, ys, B=4000):
    pts = list(zip(xs, ys)); slopes = []
    for _ in range(B):
        smp = [pts[random.randrange(len(pts))] for _ in pts]
        sx, sy = zip(*smp); slopes.append(_ols(list(sx), list(sy)))
    slopes.sort()
    return slopes[int(0.025*B)], slopes[int(0.975*B)-1]


# two-sided t critical values, df = n-2
_T95 = {2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365, 8: 2.306, 9: 2.262}


def _ols_ci(xs, ys):
    """Conventional OLS slope CI. Reported alongside the bootstrap because a
    case bootstrap on 6-9 annual aggregates understates uncertainty."""
    n = len(xs)
    if n <= 2:
        return None, None
    b = _ols(xs, ys); xb = sum(xs)/n; yb = sum(ys)/n; a = yb - b*xb
    resid = [y-(a+b*x) for x, y in zip(xs, ys)]
    s2 = sum(r*r for r in resid)/(n-2)
    sxx = sum((x-xb)**2 for x in xs)
    se = math.sqrt(s2/sxx) if sxx else 0.0
    t = _T95.get(n-2, 1.96)
    return b - t*se, b + t*se


def _classify(lo, hi, sign):
    if lo is None:
        return "not fitted"
    if lo > 0:
        return "worsening" if sign == +1 else "improving"
    if hi < 0:
        return "improving" if sign == +1 else "worsening"
    return "flat/uncertain"


def load_series():
    by_series = defaultdict(list)
    labels = {}
    with TRAJ_CSV.open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh):
            by_series[r["series"]].append((int(r["year"]), float(r["value"])))
            labels[r["series"]] = r["label"]
    for k in by_series:
        by_series[k].sort()
    return by_series, labels


def fit_trajectories(by_series):
    out = []
    for key, label, sign, brk in TARGETS:
        pts = by_series[key]
        xs = [float(y) for y, _ in pts]; ys = [v for _, v in pts]
        slope = _ols(xs, ys)
        blo, bhi = _bootstrap_ci(xs, ys)
        olo, ohi = _ols_ci(xs, ys)
        row = dict(series=key, label=label, sign=sign, n=len(pts),
                   years=f"{int(xs[0])}-{int(xs[-1])}",
                   slope_full=round(slope, 4),
                   boot_lo=round(blo, 4), boot_hi=round(bhi, 4),
                   ols_lo=(None if olo is None else round(olo, 4)),
                   ols_hi=(None if ohi is None else round(ohi, 4)),
                   classification_full=_classify(blo, bhi, sign),
                   break_year=brk or "")
        if brk:
            pre = [(y, v) for y, v in pts if y < brk]
            post = [(y, v) for y, v in pts if y >= brk]
            for tag, seg in (("pre", pre), ("post", post)):
                sx = [float(y) for y, _ in seg]; sy = [v for _, v in seg]
                s = _ols(sx, sy)
                lo, hi = _bootstrap_ci(sx, sy)
                row[f"slope_{tag}"] = round(s, 4)
                row[f"{tag}_boot_lo"] = round(lo, 4)
                row[f"{tag}_boot_hi"] = round(hi, 4)
                row[f"classification_{tag}"] = _classify(lo, hi, sign)
            row["level_shift"] = round(sum(v for _, v in post)/len(post)
                                       - sum(v for _, v in pre)/len(pre), 4)
            # primary analysis uses the within-era (post-break) slope
            row["slope_primary"] = row["slope_post"]
            row["primary_lo"] = row["post_boot_lo"]
            row["primary_hi"] = row["post_boot_hi"]
            row["classification_primary"] = row["classification_post"]
            row["primary_basis"] = f"within-era slope, {brk}-{int(xs[-1])}"
        else:
            row["slope_primary"] = row["slope_full"]
            row["primary_lo"] = row["boot_lo"]
            row["primary_hi"] = row["boot_hi"]
            row["classification_primary"] = row["classification_full"]
            row["primary_basis"] = f"full period, {row['years']}"
        out.append(row)
    return out


# --------------------------------------------------------------------------
# Stage 3: BOPAS
# --------------------------------------------------------------------------
# Mechanistic coupling M[node][series], re-specified for the revised target
# set on the same three-level scheme as the parent framework
# (0 = no plausible mechanism, 0.5 = partial/indirect, 1.0 = direct).
M_TB = {
    "1":  {"RR_new_pct": 0.0, "RRMDR_notif_100k": 0.5, "preXDR_fqr_pct": 1.0, "MDR_TSR_pct": 1.0, "XDR_TSR_pct": 1.0},
    "2":  {"RR_new_pct": 1.0, "RRMDR_notif_100k": 1.0, "preXDR_fqr_pct": 1.0, "MDR_TSR_pct": 1.0, "XDR_TSR_pct": 1.0},
    "3":  {"RR_new_pct": 0.5, "RRMDR_notif_100k": 0.5, "preXDR_fqr_pct": 1.0, "MDR_TSR_pct": 0.5, "XDR_TSR_pct": 1.0},
    "4":  {"RR_new_pct": 1.0, "RRMDR_notif_100k": 1.0, "preXDR_fqr_pct": 0.5, "MDR_TSR_pct": 0.5, "XDR_TSR_pct": 0.5},
    "5":  {"RR_new_pct": 1.0, "RRMDR_notif_100k": 1.0, "preXDR_fqr_pct": 0.5, "MDR_TSR_pct": 0.0, "XDR_TSR_pct": 0.0},
    "6":  {"RR_new_pct": 1.0, "RRMDR_notif_100k": 1.0, "preXDR_fqr_pct": 1.0, "MDR_TSR_pct": 0.5, "XDR_TSR_pct": 0.5},
    "7":  {"RR_new_pct": 0.0, "RRMDR_notif_100k": 0.5, "preXDR_fqr_pct": 0.0, "MDR_TSR_pct": 1.0, "XDR_TSR_pct": 1.0},
    "8":  {"RR_new_pct": 0.5, "RRMDR_notif_100k": 1.0, "preXDR_fqr_pct": 0.5, "MDR_TSR_pct": 1.0, "XDR_TSR_pct": 0.5},
    "9":  {"RR_new_pct": 0.0, "RRMDR_notif_100k": 0.5, "preXDR_fqr_pct": 0.5, "MDR_TSR_pct": 1.0, "XDR_TSR_pct": 1.0},
    "10": {"RR_new_pct": 0.5, "RRMDR_notif_100k": 0.5, "preXDR_fqr_pct": 0.0, "MDR_TSR_pct": 0.0, "XDR_TSR_pct": 0.0},
    "11": {"RR_new_pct": 0.5, "RRMDR_notif_100k": 0.5, "preXDR_fqr_pct": 0.0, "MDR_TSR_pct": 0.0, "XDR_TSR_pct": 0.0},
    "12": {"RR_new_pct": 0.5, "RRMDR_notif_100k": 0.0, "preXDR_fqr_pct": 0.0, "MDR_TSR_pct": 0.0, "XDR_TSR_pct": 0.0},
}

# Open-access full-text coverage per node, carried unchanged from the parent
# framework. Reported as a limitation: these are the parent study's proxy
# values, not re-derived here, which is why a C-removed sensitivity analysis
# is reported.
C_TB = {"1": 0.71, "2": 0.65, "3": 0.78, "4": 0.62, "5": 0.66, "6": 0.69,
        "7": 0.58, "8": 0.55, "9": 0.61, "10": 0.66, "11": 0.59, "12": 0.52}


def _betacf(x, a, b, max_iter=200, eps=3e-7):
    qab, qap, qam = a+b, a+1, a-1
    c = 1.0; d = 1.0 - qab*x/qap
    if abs(d) < 1e-30: d = 1e-30
    d = 1.0/d; h = d
    for mm in range(1, max_iter+1):
        m2 = 2*mm
        aa = mm*(b-mm)*x/((qam+m2)*(a+m2))
        d = 1.0 + aa*d
        if abs(d) < 1e-30: d = 1e-30
        c = 1.0 + aa/c
        if abs(c) < 1e-30: c = 1e-30
        d = 1.0/d; h *= d*c
        aa = -(a+mm)*(qab+mm)*x/((a+m2)*(qap+m2))
        d = 1.0 + aa*d
        if abs(d) < 1e-30: d = 1e-30
        c = 1.0 + aa/c
        if abs(c) < 1e-30: c = 1e-30
        d = 1.0/d; delta = d*c; h *= delta
        if abs(delta-1.0) < eps: break
    return h


def _ibeta(x, a, b):
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    lbt = (math.lgamma(a+b) - math.lgamma(a) - math.lgamma(b)
           + a*math.log(x) + b*math.log(1-x))
    bt = math.exp(lbt)
    if x < (a+1)/(a+b+2):
        return bt*_betacf(x, a, b)/a
    return 1.0 - bt*_betacf(1-x, b, a)/b


def lower_beta_bound(k, n, alpha=0.10):
    """5th percentile of Beta(k+1, n-k+1): the lower bound of the central
    90% posterior credible interval for P(favours intervention)."""
    if n == 0:
        return 0.0
    a, b = k+1, n-k+1
    target = alpha/2.0
    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = 0.5*(lo+hi)
        if _ibeta(mid, a, b) < target:
            lo = mid
        else:
            hi = mid
    return 0.5*(lo+hi)


def compute_bopas(node_counts, traj_rows, *, half_floor=HALF_WIDTH_FLOOR,
                  use_C=True, slope_field="slope_primary",
                  lo_field="primary_lo", hi_field="primary_hi"):
    L = {nid: lower_beta_bound(c.get("k_favours", 0), c.get("n_directional", 0))
         for nid, c in node_counts.items()}
    rows = []
    for nid in NODE_LABELS:
        for t in traj_rows:
            slope_nat = t[slope_field] * t["sign"]      # +ve = worsening
            half = abs(t[hi_field] - t[lo_field]) / 2.0
            half_used = max(half, half_floor) if half_floor is not None else half
            T = max(slope_nat/half_used, 0.0) if half_used > 0 else 0.0
            M = M_TB[nid][t["series"]]
            C = C_TB[nid] if use_C else 1.0
            score = L[nid] * T * M * math.sqrt(C)
            rows.append(dict(intervention=NODE_LABELS[nid], target=t["label"],
                             series=t["series"], L=round(L[nid], 3), T=round(T, 3),
                             M=f"{M:.2f}", C=f"{C:.2f}", BOPAS=round(score, 3),
                             slope_used=f"{slope_nat:+.3f}",
                             half_width_used=round(half_used, 3),
                             floor_binding=("yes" if half_floor is not None and half < half_floor else "no"),
                             trajectory_class=t["classification_primary"]))
    rows.sort(key=lambda r: -r["BOPAS"])
    return rows, L


def main() -> None:
    rows = _read_corpus()
    by_series, _ = load_series()

    print("Stage 1: coding the open-access corpus (corrected classifier) ...")
    corrected = code_corpus(rows, corrected=True)
    uncorrected = code_corpus(rows, corrected=False)
    (OUT / "ijtb_diagnostics.txt").write_text(corpus_diagnostics(rows), encoding="utf-8")

    print("Stage 2: fitting trajectories from ntep_trajectories.csv ...")
    traj = fit_trajectories(by_series)
    keys = ["series", "label", "sign", "n", "years", "slope_full", "boot_lo", "boot_hi",
            "ols_lo", "ols_hi", "classification_full", "break_year", "slope_pre",
            "pre_boot_lo", "pre_boot_hi", "classification_pre", "slope_post",
            "post_boot_lo", "post_boot_hi", "classification_post", "level_shift",
            "slope_primary", "primary_lo", "primary_hi", "classification_primary",
            "primary_basis"]
    with (OUT / "ijtb_slopes.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in traj:
            w.writerow(r)

    with (OUT / "ijtb_node_counts.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["node_id", "label", "n_total", "k_corrected", "n_dir_corrected",
                    "L_corrected", "k_uncorrected", "n_dir_uncorrected", "L_uncorrected"])
        for nid, lbl in NODE_LABELS.items():
            c, u = corrected[nid], uncorrected[nid]
            w.writerow([nid, lbl, c.get("n_total", 0),
                        c.get("k_favours", 0), c.get("n_directional", 0),
                        f"{lower_beta_bound(c.get('k_favours',0), c.get('n_directional',0)):.3f}",
                        u.get("k_favours", 0), u.get("n_directional", 0),
                        f"{lower_beta_bound(u.get('k_favours',0), u.get('n_directional',0)):.3f}"])

    print("Stage 3: computing BOPAS ...")
    grid, _ = compute_bopas(corrected, traj)
    with (OUT / "ijtb_bopas_long.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(grid[0].keys())); w.writeheader(); w.writerows(grid)
    nonzero = [r for r in grid if r["BOPAS"] > 0]
    with (OUT / "ijtb_bopas_top.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(grid[0].keys())); w.writeheader(); w.writerows(nonzero)

    print("Stage 4: sensitivity analyses ...")
    scenarios = []

    def record(name, g, note=""):
        top = [r for r in g if r["BOPAS"] > 0][:5]
        scenarios.append(dict(scenario=name, non_zero_cells=sum(1 for r in g if r["BOPAS"] > 0),
                              rank1=(f"{top[0]['intervention']} x {top[0]['target']} ({top[0]['BOPAS']})" if top else "none"),
                              rank2=(f"{top[1]['intervention']} x {top[1]['target']} ({top[1]['BOPAS']})" if len(top) > 1 else "none"),
                              rank3=(f"{top[2]['intervention']} x {top[2]['target']} ({top[2]['BOPAS']})" if len(top) > 2 else "none"),
                              note=note))

    record("Primary (corrected coding, floor 0.5, within-era slopes)", grid)
    record("Uncorrected direction-of-effect classifier",
           compute_bopas(uncorrected, traj)[0],
           "original favour-first, nomenclature-inclusive classifier")
    record("Half-width floor removed", compute_bopas(corrected, traj, half_floor=None)[0],
           "normalisation exactly as described in the parent framework")
    record("Feasibility factor C removed", compute_bopas(corrected, traj, use_C=False)[0],
           "sqrt(C) set to 1 for all nodes")
    record("Naive full-period slopes (no discontinuity correction)",
           compute_bopas(corrected, traj, slope_field="slope_full",
                         lo_field="boot_lo", hi_field="boot_hi")[0],
           "fits RR_new_pct across the 2021 reporting break")
    strict = {n: {k: (0.0 if v == 0.5 else v) for k, v in d.items()} for n, d in M_TB.items()}
    saved = dict(M_TB)
    try:
        M_TB.update(strict)
        record("M-strict (partial couplings set to 0)", compute_bopas(corrected, traj)[0])
    finally:
        M_TB.update(saved)
    with (OUT / "ijtb_sensitivity.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(scenarios[0].keys())); w.writeheader(); w.writerows(scenarios)

    print("\n--- trajectory classifications (primary) ---")
    for t in traj:
        print(f"  {t['label'][:56]:<58} {t['slope_primary']:+.3f}/yr "
              f"[{t['primary_lo']:+.3f},{t['primary_hi']:+.3f}]  {t['classification_primary']}"
              f"   ({t['primary_basis']})")
    print(f"\n--- BOPAS: {len(nonzero)} non-zero cells of {len(grid)} ---")
    for r in nonzero[:10]:
        print(f"  {r['BOPAS']:.3f}  {r['intervention'][:44]:<46} x {r['target'][:46]}")
    if not nonzero:
        print("  (none: no trajectory has a positive worsening slope)")
    print(f"\noutputs written to {OUT}")


if __name__ == "__main__":
    main()
