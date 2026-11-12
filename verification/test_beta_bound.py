"""Unit test for the hand-rolled incomplete beta function and credible bound.

Reviewer 1, Minor 6: the analysis script implements the regularised incomplete
beta function by continued fraction and inverts it by bisection, so that the
pipeline has no third-party numerical dependency. A hand-rolled special
function must be tested against a reference implementation before anyone
trusts a published number that depends on it.

This test compares against scipy. scipy is NOT a dependency of the analysis
itself -- it is used here, and only here, as the reference oracle.

Run:  python test_beta_bound.py
"""
import importlib.util
import sys
from pathlib import Path

from scipy.special import betainc          # regularised incomplete beta
from scipy.stats import beta as beta_dist  # for the quantile reference

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "code" / "build_tb_bopas_ijtb.py"
spec = importlib.util.spec_from_file_location("tb_ijtb", SRC)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

TOL_IBETA = 1e-6      # continued fraction is documented to ~3e-7
TOL_BOUND = 1e-5      # bisection runs 80 iterations, so this is generous
failures = []


def check(label, got, want, tol):
    ok = abs(got - want) <= tol
    if not ok:
        failures.append(f"{label}: got {got!r}, want {want!r}, diff {abs(got-want):.3e}")
    return ok


# ---------------------------------------------------------------------------
# 1. regularised incomplete beta against scipy.special.betainc
# ---------------------------------------------------------------------------
print("1. regularised incomplete beta I_x(a,b) vs scipy.special.betainc")
grid_x = [1e-6, 0.001, 0.01, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 0.99, 0.999]
grid_ab = [(1, 1), (0.5, 0.5), (2, 5), (5, 2), (10, 10), (1, 100), (100, 1),
           (68, 22), (138, 52), (102, 24), (50, 63), (96, 50), (15, 21)]
n_checked = worst = 0
worst_at = None
for a, b in grid_ab:
    for x in grid_x:
        got = m._ibeta(x, a, b)
        want = float(betainc(a, b, x))
        n_checked += 1
        d = abs(got - want)
        if d > worst:
            worst, worst_at = d, (x, a, b)
        check(f"I_{x}({a},{b})", got, want, TOL_IBETA)
print(f"   {n_checked} points checked over {len(grid_ab)} (a,b) pairs")
print(f"   worst absolute deviation {worst:.3e} at x={worst_at[0]}, a={worst_at[1]}, b={worst_at[2]}")

# boundary behaviour
print("\n2. boundary conditions")
for a, b in [(1, 1), (3, 7)]:
    check(f"I_0({a},{b}) == 0", m._ibeta(0.0, a, b), 0.0, 0.0)
    check(f"I_1({a},{b}) == 1", m._ibeta(1.0, a, b), 1.0, 0.0)
    check(f"I_-0.5({a},{b}) == 0", m._ibeta(-0.5, a, b), 0.0, 0.0)
    check(f"I_1.5({a},{b}) == 1", m._ibeta(1.5, a, b), 1.0, 0.0)
print("   0, 1 and out-of-range arguments clamp correctly")

# ---------------------------------------------------------------------------
# 3. the credible bound itself against the exact Beta quantile
# ---------------------------------------------------------------------------
# lower_beta_bound(k, n, alpha) returns the alpha/2 quantile of
# Beta(k+1, n-k+1) -- the lower limit of the central 100(1-alpha)% interval.
print("\n3. lower_beta_bound(k, n) vs scipy.stats.beta.ppf(alpha/2, k+1, n-k+1)")
cases = [(0, 0), (0, 1), (1, 1), (0, 10), (5, 10), (10, 10), (1, 3),
         (67, 89), (137, 207), (129, 188), (87, 172), (75, 128), (49, 111),
         (14, 34), (1, 3), (95, 144), (45, 82), (24, 41), (43, 76), (999, 1000)]
print(f"   {'k':>5}{'n':>6}{'ours':>12}{'scipy':>12}{'diff':>11}")
for k, n in cases:
    got = m.lower_beta_bound(k, n)
    if n == 0:
        want = 0.0                      # documented special case
    else:
        want = float(beta_dist.ppf(0.05, k + 1, n - k + 1))
    check(f"bound({k},{n})", got, want, TOL_BOUND)
    print(f"   {k:>5}{n:>6}{got:>12.6f}{want:>12.6f}{abs(got-want):>11.2e}")

# ---------------------------------------------------------------------------
# 4. the twelve values the paper actually publishes
# ---------------------------------------------------------------------------
print("\n4. every node bound reported in the manuscript, against scipy")
import csv
OUT = HERE.parent / "outputs" / "ijtb_node_counts.csv"
rows = list(csv.DictReader(OUT.open(encoding="utf-8")))
for r in rows:
    k, n = int(r["k_corrected"]), int(r["n_dir_corrected"])
    published = float(r["L_corrected"])
    want = float(beta_dist.ppf(0.05, k + 1, n - k + 1))
    ok = check(f"node {r['node_id']} L", published, want, 5e-4)  # published to 3 dp
    print(f"   node {r['node_id']:>2}  k/n={k}/{n:<4} published {published:.3f}  "
          f"scipy {want:.6f}  {'OK' if ok else 'MISMATCH'}")

# ---------------------------------------------------------------------------
# 5. monotonicity and range sanity
# ---------------------------------------------------------------------------
print("\n5. structural properties")
prev = -1.0
mono = True
for k in range(0, 51):
    v = m.lower_beta_bound(k, 50)
    if v < prev:
        mono = False
    prev = v
    if not (0.0 <= v <= 1.0):
        failures.append(f"bound({k},50) = {v} outside [0,1]")
print(f"   bound is non-decreasing in k at fixed n : {mono}")
if not mono:
    failures.append("bound not monotone in k")
print("   all bounds within [0, 1]                : "
      f"{not any('outside' in f for f in failures)}")
# tighter data should raise the lower bound
a1, a2 = m.lower_beta_bound(9, 10), m.lower_beta_bound(90, 100)
print(f"   bound(9/10)={a1:.4f} < bound(90/100)={a2:.4f}    : {a1 < a2}")
if not a1 < a2:
    failures.append("bound did not tighten with more data at the same proportion")

# ---------------------------------------------------------------------------
print("\n" + "=" * 62)
if failures:
    print(f"FAILED: {len(failures)} check(s)")
    for f in failures:
        print("  -", f)
    sys.exit(1)
print("ALL CHECKS PASSED")
print("The pure-stdlib incomplete beta and credible bound agree with scipy to")
print(f"within {TOL_IBETA:.0e} (function) and {TOL_BOUND:.0e} (quantile), and every node")
print("bound published in the manuscript is reproduced by the reference implementation.")
