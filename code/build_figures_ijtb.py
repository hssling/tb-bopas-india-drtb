"""Figures for the IJTB submission. Emits 400-dpi PNG + TIFF."""
from __future__ import annotations
import csv, math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = ROOT / "figures"; OUT.mkdir(exist_ok=True)
OUTPUTS = ROOT / "outputs"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.linewidth": 0.7, "axes.edgecolor": "#333333",
    "axes.labelsize": 8.5, "axes.titlesize": 9,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5, "legend.frameon": False,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.08,
})
INK, ACCENT, WARN, MUTED = "#1c1c1c", "#1f6f8b", "#b3452f", "#9a9a9a"


def save(fig, stem):
    for ext in ("png", "tiff"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=400)
    plt.close(fig)
    print("  wrote", stem)


def load_series():
    d = defaultdict(list); labels = {}
    for r in csv.DictReader((ROOT / "data" / "ntep_trajectories.csv").open(encoding="utf-8")):
        d[r["series"]].append((int(r["year"]), float(r["value"])))
        labels[r["series"]] = r["label"]
    for k in d:
        d[k].sort()
    return d, labels


def load_slopes():
    return {r["series"]: r for r in csv.DictReader((OUTPUTS / "ijtb_slopes.csv").open(encoding="utf-8"))}


def ols(xs, ys):
    n = len(xs); xb = sum(xs)/n; yb = sum(ys)/n
    sxx = sum((x-xb)**2 for x in xs)
    b = sum((x-xb)*(y-yb) for x, y in zip(xs, ys))/sxx if sxx else 0.0
    return b, yb - b*xb


# --------------------------------------------------------------- Figure 1
def figure1(series, labels, slopes):
    order = ["RR_new_pct", "RRMDR_notif_100k", "preXDR_fqr_pct", "MDR_TSR_pct", "XDR_TSR_pct"]
    short = {"RR_new_pct": "Rifampicin resistance,\nnew pulmonary cases (%)",
             "RRMDR_notif_100k": "RR/MDR-TB notifications\nper 100 000",
             "preXDR_fqr_pct": "Fluoroquinolone-resistant\nRR-TB (pre-XDR, %)",
             "MDR_TSR_pct": "MDR/RR-TB treatment\nsuccess (%)",
             "XDR_TSR_pct": "XDR-TB treatment\nsuccess (%)"}
    fig, axes = plt.subplots(1, 5, figsize=(11.6, 2.6))
    for ax, key in zip(axes, order):
        pts = series[key]; s = slopes[key]
        xs = [y for y, _ in pts]; ys = [v for _, v in pts]
        cls = s["classification_primary"]
        col = WARN if cls == "worsening" else (ACCENT if cls == "improving" else MUTED)
        ax.plot(xs, ys, "o-", color=col, lw=1.3, ms=3.6, zorder=3)
        if key == "RR_new_pct":
            brk = int(s["break_year"])
            ax.axvline(brk - 0.5, color=INK, lw=0.8, ls=(0, (3, 2)), zorder=1)
            for seg in ([p for p in pts if p[0] < brk], [p for p in pts if p[0] >= brk]):
                sx = [y for y, _ in seg]; sy = [v for _, v in seg]
                b, a = ols([float(x) for x in sx], sy)
                ax.plot(sx, [a + b*x for x in sx], color=INK, lw=1.0, zorder=2)
            ax.annotate("reporting\nbreak", xy=(brk - 0.5, max(ys)), xytext=(-2, -2),
                        textcoords="offset points", ha="right", va="top",
                        fontsize=6.8, color=INK)
        else:
            b, a = ols([float(x) for x in xs], ys)
            ax.plot(xs, [a + b*x for x in xs], color=INK, lw=1.0, zorder=2)
        ax.set_title(short[key], fontsize=8, pad=6)
        ax.xaxis.set_major_locator(MaxNLocator(4, integer=True))
        ax.tick_params(length=2.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        pos = (0.98, 0.03, "right") if key == "RR_new_pct" else (0.02, 0.03, "left")
        ax.text(pos[0], pos[1], f"{float(s['slope_primary']):+.2f} pp/yr\n{cls}",
                transform=ax.transAxes, fontsize=6.8, color=col, va="bottom",
                ha=pos[2], bbox=dict(facecolor="white", edgecolor="none",
                                     alpha=0.85, pad=1.5))
    axes[0].set_ylabel("value")
    fig.suptitle("India, WHO Global TB Programme country data: five drug-resistant tuberculosis "
                 "surveillance axes with fitted trends", fontsize=9, y=1.06)
    save(fig, "fig_ijtb_1_trajectories")


# --------------------------------------------------------------- Figure 2
def figure2(series, slopes):
    s = slopes["RR_new_pct"]; pts = series["RR_new_pct"]; brk = int(s["break_year"])
    xs = [y for y, _ in pts]; ys = [v for _, v in pts]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(8.6, 3.1),
                                  gridspec_kw={"width_ratios": [1.45, 1]})
    ax.plot(xs, ys, "o", color=INK, ms=4.5, zorder=4)
    bf, af = ols([float(x) for x in xs], ys)
    ax.plot(xs, [af + bf*x for x in xs], color=WARN, lw=1.6, zorder=3,
            label=f"naive full-period fit: {bf:+.2f} pp/yr (classified worsening)")
    for tag, seg, off, ha in (("2017-2020", [p for p in pts if p[0] < brk], (4, 8), "left"),
                              ("2021-2024", [p for p in pts if p[0] >= brk], (-4, 10), "right")):
        sx = [y for y, _ in seg]; sy = [v for _, v in seg]
        b, a = ols([float(x) for x in sx], sy)
        ax.plot(sx, [a + b*x for x in sx], color=ACCENT, lw=1.8, zorder=3)
        ax.annotate(f"{tag}: {b:+.2f} pp/yr", xy=(sx[-1], a + b*sx[-1]),
                    xytext=off, textcoords="offset points", fontsize=7,
                    color=ACCENT, ha=ha)
    ax.axvline(brk - 0.5, color=MUTED, lw=0.9, ls=(0, (3, 2)))
    ax.annotate(f"level shift {float(s['level_shift']):+.2f} pp\n(WHO 2021 definition revision;\nuniversal DST scale-up)",
                xy=(brk - 0.5, 3.4), xytext=(6, 0), textcoords="offset points",
                fontsize=7, color=MUTED, va="center")
    ax.set_xlabel("year"); ax.set_ylabel("rifampicin resistance among new\npulmonary TB cases (%)")
    ax.set_title("A. One reporting discontinuity, two opposite conclusions", fontsize=8.5, loc="left")
    ax.legend(loc="lower right", fontsize=7)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    rows = {r["scenario"]: r for r in csv.DictReader((OUTPUTS / "ijtb_sensitivity.csv").open(encoding="utf-8"))}
    naive = next(r for k, r in rows.items() if k.startswith("Naive"))
    prim = next(r for k, r in rows.items() if k.startswith("Primary"))
    vals = [float(prim["rank1"].rsplit("(", 1)[1].rstrip(")")),
            float(naive["rank1"].rsplit("(", 1)[1].rstrip(")"))]
    bars = ax2.bar(["discontinuity\ncorrected", "uncorrected"], vals,
                   color=[ACCENT, WARN], width=0.55)
    for b, v in zip(bars, vals):
        ax2.text(b.get_x() + b.get_width()/2, v, f"{v:.3f}", ha="center", va="bottom", fontsize=8)
    ax2.set_ylabel("highest priority score returned")
    ax2.set_title("B. Effect on the top-ranked priority", fontsize=8.5, loc="left")
    ax2.set_ylim(0, max(vals)*1.28)
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)
    ax2.text(0.5, 0.60, "spurious priority on an axis that is\nimproving within the current era",
             transform=ax2.transAxes, ha="center", va="center",
             fontsize=7, color=WARN)
    save(fig, "fig_ijtb_2_discontinuity")


# --------------------------------------------------------------- Figure 3
def figure3():
    nodes = list(csv.DictReader((OUTPUTS / "ijtb_node_counts.csv").open(encoding="utf-8")))
    nodes.sort(key=lambda r: -float(r["L_uncorrected"]))
    labels = [r["label"].replace(" (Xpert MTB/RIF, Xpert MTB/XDR)", "")
                        .replace(" / occupational respiratory protection", "")
                        .replace(" for drug-susceptibility prediction", " for DST") for r in nodes]
    unc = [float(r["L_uncorrected"]) for r in nodes]
    cor = [float(r["L_corrected"]) for r in nodes]
    y = range(len(nodes))
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    for i, (u, c) in enumerate(zip(unc, cor)):
        ax.plot([c, u], [i, i], color="#cccccc", lw=1.4, zorder=1, solid_capstyle="round")
    ax.scatter(unc, list(y), s=34, color=WARN, zorder=3, label="original classifier")
    ax.scatter(cor, list(y), s=34, color=ACCENT, zorder=3, label="corrected classifier")
    ax.set_yticks(list(y)); ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("lower bound of the 90% credible interval for P(favours intervention)")
    ax.set_xlim(0, 1)
    ax.set_title("Effect of correcting the direction-of-effect classifier on every\n"
                 "intervention node's favourability bound", fontsize=9, loc="left", pad=8)
    ax.legend(loc="lower right")
    ax.grid(axis="x", color="#eeeeee", lw=0.7)
    ax.set_axisbelow(True)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="y", length=0)
    save(fig, "fig_ijtb_3_classifier")


if __name__ == "__main__":
    series, labels = load_series()
    slopes = load_slopes()
    print("building figures ...")
    figure1(series, labels, slopes)
    figure2(series, slopes)
    figure3()
    print("done ->", OUT)
