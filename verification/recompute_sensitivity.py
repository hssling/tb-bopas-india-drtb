"""Read-only sensitivity probes on the published TB-BOPAS pipeline."""
import importlib.util, csv, math, re
from pathlib import Path
spec = importlib.util.spec_from_file_location('tb', 'code/build_tb_bopas_ijtb.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

CORPUS = Path('data/tb_interventions_unified.csv')
rows = list(csv.DictReader(open(CORPUS, encoding='utf-8', newline='')))
traj = m.fit_trajectories()


def bopas_grid(node_counts, half_floor):
    L = {n: m._lower_beta_bound(c.get('k_favours', 0), c.get('n_directional', 0))
         for n, c in node_counts.items()}
    out = []
    for nid in m.NODE_LABELS:
        for t in traj:
            slope = float(t['slope']) * int(t['sign'])
            half = abs(float(t['ci_hi']) - float(t['ci_lo'])) / 2.0
            if half_floor is not None:
                half = max(half, half_floor)
            T = max(slope / half, 0.0)
            v = L[nid] * T * m.M_TB[nid][t['key']] * math.sqrt(m.C_TB[nid])
            out.append((v, m.NODE_LABELS[nid], t['label'], L[nid], T))
    out.sort(reverse=True)
    return out, L


# corrected direction-of-effect coder:
#  1. contrary / no-effect cues take precedence over favour cues
#  2. intervention-nomenclature words ("prevent", "prevention") no longer count
#     as evaluative favour cues
FAV = re.compile(r"\b(reduce[ds]?|decreased?|lower(?:ed)?|improve[ds]?|effective|"
                 r"successful|favou?rable|protective)\b", re.I)


def corrected_counts():
    counts = {n: {'k_favours': 0, 'n_directional': 0, 'n_total': 0} for n in m.NODE_LABELS}
    for r in rows:
        text = ' '.join([r.get('title', ''), r.get('abstract_present', '')])
        for nid, rx in m.NODE_REGEX.items():
            if not rx.search(text):
                continue
            counts[nid]['n_total'] += 1
            if m.RE_UNCLEAR.search(text):
                continue
            if m.RE_CONTRARY.search(text):        # precedence: contrary wins
                counts[nid]['n_directional'] += 1
            elif FAV.search(text):
                counts[nid]['k_favours'] += 1
                counts[nid]['n_directional'] += 1
    return counts


published = m.code_corpus(CORPUS)
corrected = corrected_counts()

print('=== A. is the CI half-width floor of 0.5 active? ===')
for t in traj:
    half = abs(float(t['ci_hi']) - float(t['ci_lo'])) / 2
    flag = 'FLOOR ACTIVE' if half < 0.5 else 'floor inactive'
    print(f"  {t['label']:38} half-width={half:.4f} -> used {max(half,0.5):.3f}  {flag}")

base, Lp = bopas_grid(published, 0.5)
nofloor, _ = bopas_grid(published, None)
print('\n=== B. top 8: as published (floor 0.5) vs floor removed ===')
for i in range(8):
    a = f"{base[i][0]:.3f} {base[i][1][:30]} x {base[i][2][:20]}"
    b = f"{nofloor[i][0]:.3f} {nofloor[i][1][:30]} x {nofloor[i][2][:20]}"
    print(f"{i+1:>3}  {a:<58}|  {b}")

corr, Lc = bopas_grid(corrected, 0.5)
print('\n=== C. corrected direction-of-effect coding: k/n and L per node ===')
print(f"{'node':<50}{'k/n pub':>14}{'L pub':>8}{'k/n corr':>14}{'L corr':>8}")
for nid, lbl in m.NODE_LABELS.items():
    p, c = published[nid], corrected[nid]
    kp = f"{p.get('k_favours',0)}/{p.get('n_directional',0)}"
    kc = f"{c['k_favours']}/{c['n_directional']}"
    print(f"{lbl[:48]:<50}{kp:>14}{Lp[nid]:>8.3f}{kc:>14}{Lc[nid]:>8.3f}")

print('\n=== D. top 8 under corrected coding (floor kept) ===')
for i in range(8):
    print(f"  {i+1}. {corr[i][0]:.3f}  {corr[i][1][:40]} x {corr[i][2][:26]}")

print('\n=== E. non-zero cell counts ===')
for name, g in (('as published', base), ('corrected coding', corr), ('floor removed', nofloor)):
    print(f"  {name:<18} {sum(1 for x in g if x[0] > 0)} of {len(g)}")
