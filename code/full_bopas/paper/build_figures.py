#!/usr/bin/env python3
"""Figures for the TB intervention evidence map (matplotlib). Reads outputs/full_bopas/results.json only.

Reads outputs/full_bopas/results.json; writes outputs/full_bopas/figures/*.png (300 dpi).
Okabe-Ito colour-blind-safe palette; drawn at final print width (174 mm); no text below 7 pt.
A layout check fails the build if any text leaves the figure or overlaps another text.
Usage: python -I code/full_bopas/paper/build_figures.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.text import Text

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "outputs" / "full_bopas"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
R = json.loads((OUT / "results.json").read_text(encoding="utf-8"))
ROWS = R["interventions"]
DPI = 600
W = 174 / 25.4
# Okabe-Ito
BLUE, ORANGE, SKY, GREEN, YELLOW, VERM, PURPLE, GREY = ("#0072B2", "#E69F00", "#56B4E9", "#009E73",
                                                         "#F0E442", "#D55E00", "#CC79A7", "#BBBBBB")
INK = "#1c1c1c"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.linewidth": 0.6,
                     "axes.titlesize": 8, "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.fontsize": 7})

ORDER = ["5", "2", "1", "3", "4", "9", "8", "7", "6", "10", "12", "11"]   # human health, then One Health


def ordered():
    by = {r["node_id"]: r for r in ROWS}
    return [by[i] for i in ORDER]


def label(r):
    return r["intervention"] + ("*" if r["group"] == "One Health" else "")


def check_layout(fig, name):
    fig.canvas.draw()
    rnd = fig.canvas.get_renderer()
    fb = fig.bbox
    texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    bad = []
    ext = []
    for t in texts:
        e = t.get_window_extent(rnd)
        if e.x0 < fb.x0 - 1 or e.x1 > fb.x1 + 1 or e.y0 < fb.y0 - 1 or e.y1 > fb.y1 + 1:
            bad.append(f"leaves figure: {t.get_text()!r}")
        if t.get_fontsize() < 7:
            bad.append(f"below 7 pt: {t.get_text()!r}")
        ext.append((t, e))
    for i in range(len(ext)):
        for j in range(i + 1, len(ext)):
            a, b = ext[i][1], ext[j][1]
            if a.overlaps(b) and a.width > 0 and b.width > 0:
                ov = (min(a.x1, b.x1) - max(a.x0, b.x0)) * (min(a.y1, b.y1) - max(a.y0, b.y0))
                if ov > 0.15 * min(a.width * a.height, b.width * b.height):
                    bad.append(f"overlap: {ext[i][0].get_text()!r} / {ext[j][0].get_text()!r}")
    if bad:
        raise SystemExit(f"{name} layout check failed:\n  " + "\n  ".join(bad))
    print(f"  {name}: {len(texts)} texts inside the figure, none overlapping, none below 7 pt")


def save(fig, stem):
    check_layout(fig, stem)
    fig.savefig(FIG / f"{stem}.png", dpi=300, facecolor="white")
    print(f"  {stem}.png written")
    plt.close(fig)


def figure1():
    rows = ordered()[::-1]
    fig, ax = plt.subplots(figsize=(W, 4.2))
    fig.subplots_adjust(left=0.37, right=0.90, top=0.90, bottom=0.12)
    y = range(len(rows))
    segs = [("evaluations", "Evaluates the intervention (tier 1)", BLUE),
            ("tier2", "Relevant, does not evaluate it (tier 2)", SKY),
            ("off_topic", "Off-topic", GREY)]
    for key, lab, col in segs:
        left = [sum(r[k] for k, _, _ in segs[:[s[0] for s in segs].index(key)]) / r["records"] for r in rows]
        w = [r[key] / r["records"] for r in rows]
        ax.barh(list(y), w, left=left, color=col, edgecolor="white", linewidth=0.4, height=0.72, label=lab)
        for yi, (l, wi, r) in enumerate(zip(left, w, rows)):
            if wi >= 0.06:
                ax.text(l + wi / 2, yi, f"{r[key]}", ha="center", va="center", fontsize=7,
                        color="white" if col in (BLUE,) else INK)
    for yi, r in enumerate(rows):
        ax.text(1.01, yi, f"n = {r['records']}", va="center", ha="left", fontsize=7)
    ax.set_yticks(list(y))
    ax.set_yticklabels([label(r) for r in rows])
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax.set_xlabel("Share of records retrieved for the intervention")
    ax.axhline(2.5, color=INK, lw=0.6, ls=(0, (3, 2)))
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, loc="upper left", bbox_to_anchor=(0.01, 0.995), ncol=3, frameon=False, handlelength=1.2,
               columnspacing=1.0)
    fig.text(0.01, 0.02, "* One Health interventions (below the dashed line)", fontsize=7)
    save(fig, "Figure_1")


def figure2():
    rows = ordered()[::-1]
    fig = plt.figure(figsize=(W, 4.6))
    ax = fig.add_axes([0.37, 0.14, 0.33, 0.72])
    bx = fig.add_axes([0.80, 0.14, 0.18, 0.72])
    groups = [("SR/MA or RCT", lambda r: r["SR"] + r["RCT"], BLUE),
              ("Non-randomised comparative", lambda r: r["NRS"], GREEN),
              ("Diagnostic accuracy", lambda r: r["accuracy"], ORANGE),
              ("Programme report", lambda r: r["programme"], VERM),
              ("Economic, modelling, other", lambda r: r["economic"] + r["modelling"] + r["other_design"], GREY)]
    left = [0] * len(rows)
    for lab, f, col in groups:
        w = [f(r) for r in rows]
        ax.barh(range(len(rows)), w, left=left, color=col, edgecolor="white", linewidth=0.4, height=0.72, label=lab)
        left = [a + b for a, b in zip(left, w)]
    for yi, r in enumerate(rows):
        ax.text(r["evaluations"] + 1.5, yi, str(r["evaluations"]), va="center", fontsize=7)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([label(r) for r in rows])
    ax.set_xlim(0, 135)
    ax.set_xlabel("Evaluations (number)")
    ax.axhline(2.5, color=INK, lw=0.6, ls=(0, (3, 2)))
    ax.set_title("A  Design of the evaluations", loc="left", x=-1.05, fontweight="bold")
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, loc="upper left", bbox_to_anchor=(0.01, 0.995), ncol=3, frameon=False, handlelength=1.2,
               columnspacing=0.9)
    dd = R["direction_by_design"]
    order = [("randomised trial", "RCT"), ("systematic review/meta-analysis", "SR/MA"),
             ("non-randomised comparative (cohort, quasi-experimental, before-after, interrupted time series)", "Non-rand."),
             ("diagnostic accuracy", "Accuracy"), ("economic evaluation", "Economic"),
             ("implementation/programme report (no comparator)", "Programme"), ("modelling", "Modelling")]
    for yi, (k, lab) in enumerate(order[::-1]):
        d = dd[k]
        n = d["favours"] + d["no_effect"]
        bx.errorbar(d["fav_share"], yi, xerr=[[d["fav_share"] - d["lo"]], [d["hi"] - d["fav_share"]]],
                    fmt="o", color=BLUE if lab in ("RCT", "SR/MA") else INK, ms=3.5, elinewidth=0.9, capsize=1.8)
        bx.text(0.02, yi + 0.33, f"{d['favours']}/{n}", fontsize=7, va="center", ha="left", color=INK)
    bx.set_yticks(range(len(order)))
    bx.set_yticklabels([lab for _, lab in order[::-1]])
    bx.set_xlim(0, 1.04)
    bx.set_ylim(-0.6, len(order) - 0.2)
    bx.set_xticks([0, 0.5, 1])
    bx.set_xticklabels(["0%", "50%", "100%"])
    bx.set_xlabel("Favourable share (95% CI)")
    bx.set_title("B  Direction by design", loc="left", x=-0.62, fontweight="bold")
    save(fig, "Figure_2")


def figure3():
    rows = ordered()[::-1]
    fig = plt.figure(figsize=(W, 4.3))
    ax = fig.add_axes([0.37, 0.13, 0.30, 0.74])
    bx = fig.add_axes([0.73, 0.13, 0.25, 0.74])
    for yi, r in enumerate(rows):
        col = ORANGE if r["group"] == "One Health" else BLUE
        ax.errorbar(r["eval_share"], yi, xerr=[[r["eval_share"] - r["eval_lo"]], [r["eval_hi"] - r["eval_share"]]],
                    fmt="o", color=col, ms=3.5, elinewidth=0.9, capsize=1.8)
    g = R["groups"]
    for name, col in (("Human health", BLUE), ("One Health", ORANGE)):
        ax.axvline(g[name]["share"][0], color=col, lw=0.7, ls=(0, (2, 2)))
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([label(r) for r in rows])
    ax.set_xlim(0, 0.5)
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.set_xticklabels(["0%", "10%", "20%", "30%", "40%", "50%"])
    ax.set_xlabel("Records that evaluate the intervention (95% CI)")
    ax.axhline(2.5, color=INK, lw=0.6, ls=(0, (3, 2)))
    ax.set_title("A  Evaluation share", loc="left", x=-1.2, fontweight="bold")
    ax.text(0.5, 1.035, f"dotted lines: human health {g['Human health']['share'][0]:.0%}; One Health "
            f"{g['One Health']['share'][0]:.0%}", transform=ax.transAxes, fontsize=7, ha="center")
    ind = [r["india"] for r in rows]
    oth = [r["evaluations"] - r["india"] for r in rows]
    bx.barh(range(len(rows)), ind, color=VERM, height=0.72, label="India setting")
    bx.barh(range(len(rows)), oth, left=ind, color=GREY, height=0.72, label="Other or not stated")
    for yi, r in enumerate(rows):
        bx.text(r["evaluations"] + 2, yi, f"{r['india']}/{r['evaluations']}", va="center", fontsize=7)
    bx.set_yticks(range(len(rows)))
    bx.set_yticklabels([])
    bx.set_xlim(0, 150)
    bx.set_xlabel("Evaluations (number)")
    bx.axhline(2.5, color=INK, lw=0.6, ls=(0, (3, 2)))
    bx.set_title("B  India-setting evaluations", loc="left", fontweight="bold")
    for axx in (ax, bx):
        axx.set_ylim(-0.6, len(rows) - 0.4)
    bx.legend(loc="lower right", frameon=False, handlelength=1.2)
    fig.text(0.01, 0.015, "* One Health interventions (below the dashed lines)", fontsize=7)
    save(fig, "Figure_3")


def graphical_abstract():
    o = R["overall"]
    dd = R["direction_by_design"]
    strong = sum(r["strong"] for r in ROWS)
    fig = plt.figure(figsize=(W, 4.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 74)
    ax.axis("off")
    ax.text(2, 70, "Tuberculosis interventions: what the retrieved literature actually evaluates",
            fontsize=9.5, fontweight="bold", va="center")
    steps = [(f"{R['analysed']:,} records retrieved".replace(",", " "), 58, GREY),
             (f"{o['evaluations']} evaluate the intervention ({o['eval_share_ci'][0]:.0%})", 49, BLUE),
             (f"{strong} systematic reviews or trials", 40, GREEN),
             (f"{R['india_strong']} of these set in India", 31, VERM)]
    for i, (txt, yy, col) in enumerate(steps):
        w = 40 - i * 6
        ax.add_patch(plt.Rectangle((2 + (40 - w) / 2, yy - 3.6), w, 7.2, color=col, alpha=0.9 if i else 0.6))
        ax.text(22, yy, txt, ha="center", va="center", fontsize=7.5, color="white" if i else INK)
    ax.text(22, 20, f"Favourable share: trials {dd['randomised trial']['fav_share']:.0%},",
            ha="center", va="center", fontsize=7.5)
    ax.text(22, 15.5, f"programme reports {dd['implementation/programme report (no comparator)']['fav_share']:.0%}",
            ha="center", va="center", fontsize=7.5)
    ax.text(22, 9, "Direction counts need design restriction", ha="center", va="center", fontsize=7.5, style="italic")
    tiers = [("Ready to scale with monitoring", BLUE), ("Implement and evaluate in India", GREEN),
             ("Evaluate before (further) scale-up", ORANGE), ("Research gap", VERM)]
    ax.text(47, 63, "Evidence-readiness tier (rule-based)", fontsize=8, fontweight="bold", va="center")
    yy = 57
    for t, col in tiers:
        names = [r["intervention"] for r in ROWS if r["tier"] == t]
        lines, cur = [], ""
        for nm in names:                       # greedy wrap at about 62 characters per line
            nxt = (cur + "; " + nm) if cur else nm
            if len(nxt) > 62 and cur:
                lines.append(cur + ";")
                cur = nm
            else:
                cur = nxt
        lines.append(cur)
        ax.add_patch(plt.Rectangle((47, yy - 1.8 - 3.2 * len(lines)), 1.6, 3.6 + 3.2 * len(lines), color=col))
        ax.text(50, yy, t, fontsize=7.5, fontweight="bold", va="center")
        for li, ln in enumerate(lines):
            ax.text(50, yy - 3.2 - 3.2 * li, ln, fontsize=7, va="center")
        yy -= 5.8 + 3.2 * len(lines)
    save(fig, "Graphical_Abstract")


if __name__ == "__main__":
    figure1()
    figure2()
    figure3()
    graphical_abstract()
