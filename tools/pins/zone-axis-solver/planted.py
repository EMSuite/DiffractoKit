"""ZoneAxisSolver planted set (DiffractoKit S5): spot patterns at known orientations, planted with py4DSTEM.

  ../References/python/py4dstem/bin/python -I -B tools/pins/zone-axis-solver/planted.py

Writes Tests/DiffractoGeometryTests/Fixtures/zone-axis-planted.json. No Python function solves one pattern alone (py4DSTEM's
match_single_pattern is a template match on a 2 degree grid), so the pin is the planted truth:
  * 8 phases: Al, Au, theta' (Al2Cu), beta'' (Mg5Si6), Mg (cell and sites read from tools/pins/phase-vectors' fixture), Si,
    alpha-Fe and Mg2Si (written below);
  * per phase the first 14 zone axes [u v w] with |u|, |v|, |w| <= 3, one per class of the structure's proper rotations
    and the sign (spglib), by u^2 + v^2 + w^2;
  * per zone axis 2 orientations, py4DSTEM's `parse_orientation(zone_axis_lattice=uvw) @ m3z(gamma)` at a seeded random
    gamma, and per orientation 5 conditions;
  * the pattern: `generate_diffraction_pattern(orientation_matrix=M, sigma_excitation_error=0.02,
    tol_excitation_error_mult=3, tol_intensity=1e-5, k_max=1.5)` at 200 kV, 000 dropped, spots below 1e-3 of the strongest
    dropped, the strongest 40 kept;
  * conditions (1 px = 0.00125 1/Angstrom, a 2048 px detector across +-1.28 1/Angstrom): clean (Gaussian noise 0.5 px),
    scale (positions times 1.02, noise 0.5 px), missing (noise 1 px, 30 % of the spots dropped, 3 false spots at
    0.25 to 1.2 1/Angstrom), click3 and click6 (the 3 or 6 shortest zero-layer spots, h u + k v + l w = 0, of at least 10 %
    of the strongest, one per direction, noise 0.5 px).
"""
import importlib.metadata as md
import itertools
import json
import math
from pathlib import Path

import numpy as np
import spglib
from py4DSTEM.process.diffraction import Crystal

KIT = Path(__file__).resolve().parents[3]
SRC = KIT / 'Tests/DiffractoOrientationTests/Fixtures/py4dstem-phase-vector-pins.json'
OUT = KIT / 'Tests/DiffractoGeometryTests/Fixtures/zone-axis-planted.json'
SEED = 20261123
K_MAX, SIGMA, MULT, TOL_I, TOL_F, VOLTS = 1.5, 0.02, 3.0, 1e-5, 1e-4, 200e3
DQ, ZONES, REPEATS, STRONGEST = 0.00125, 14, 2, 40
assert md.version('py4DSTEM') == '0.14.19'

FCC = [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)]
EXTRA = [
    {'name': 'Si', 'cell': dict(a=5.4309, b=5.4309, c=5.4309, alpha=90, beta=90, gamma=90),
     'sites': [{'z': 14, 'frac': list(p), 'occupancy': 1.0} for p in FCC + [tuple(x + .25 for x in f) for f in FCC]]},
    {'name': 'Fe', 'cell': dict(a=2.8665, b=2.8665, c=2.8665, alpha=90, beta=90, gamma=90),
     'sites': [{'z': 26, 'frac': [0, 0, 0], 'occupancy': 1.0}, {'z': 26, 'frac': [.5, .5, .5], 'occupancy': 1.0}]},
    {'name': 'Mg2Si', 'cell': dict(a=6.351, b=6.351, c=6.351, alpha=90, beta=90, gamma=90),
     'sites': [{'z': 14, 'frac': list(p), 'occupancy': 1.0} for p in FCC]
            + [{'z': 12, 'frac': [(x + d) % 1 for x in f], 'occupancy': 1.0} for f in FCC for d in (.25, .75)]},
]


def m3z(g):
    c, s = np.cos(g), np.sin(g)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])


def crystal(ph):
    c = ph['cell']
    pos = np.array([s['frac'] for s in ph['sites']], dtype=np.float64)
    num = np.array([s['z'] for s in ph['sites']], dtype=np.intp)
    occ = np.array([s['occupancy'] for s in ph['sites']], dtype=np.float64)
    cr = Crystal(pos, num, np.array([c['a'], c['b'], c['c'], c['alpha'], c['beta'], c['gamma']]), occupancy=occ)
    cr.calculate_structure_factors(k_max=K_MAX, tol_structure_factor=TOL_F)
    cr.setup_diffraction(VOLTS)
    return cr, pos, num


def zone_axes(cr, pos, num):
    lattice = cr.lat_real
    sym = spglib.get_symmetry((lattice, pos, num), symprec=1e-3)
    rots = [r for r in sym['rotations'] if round(np.linalg.det(r)) == 1]
    seen, out = set(), []
    for t in sorted(itertools.product(range(-3, 4), repeat=3), key=lambda t: (sum(x * x for x in t), t)):
        if t == (0, 0, 0) or math.gcd(math.gcd(abs(t[0]), abs(t[1])), abs(t[2])) != 1:
            continue
        key = max(tuple(int(x) for x in s * (r @ np.array(t))) for r in rots for s in (1, -1))
        if key in seen:
            continue
        seen.add(key)
        out.append(t)
    return out[:ZONES]


def main():
    rng = np.random.default_rng(SEED)
    src = json.load(open(SRC))
    phases = [p for p in src['phases'] if p['name'] in ('Al', 'Au', 'theta_prime_Al2Cu', 'beta_double_prime_Mg5Si6', 'Mg')]
    phases += EXTRA
    out = {'generator': 'tools/pins/zone-axis-solver/planted.py', 'seed': SEED,
           'versions': {'py4DSTEM': md.version('py4DSTEM'), 'numpy': np.__version__, 'spglib': spglib.__version__},
           'commit': 'f050d207851f8accd3aa550365baaf1f38230428', 'volts': VOLTS, 'pixel': DQ,
           'matrix_layout': 'row-major 9 numbers; columns are lab x, y, z in the crystal Cartesian frame',
           'spots_layout': 'qx0, qy0, qx1, qy1, ... (1/Angstrom)', 'phases': []}
    for ph in phases:
        cr, pos, num = crystal(ph)
        rec = {'name': ph['name'], 'cell': ph['cell'], 'sites': ph['sites'], 'cases': []}
        for uvw in zone_axes(cr, pos, num):
            for _ in range(REPEATS):
                M = cr.parse_orientation(zone_axis_lattice=np.array(uvw, float)) @ m3z(rng.uniform(0, 2 * math.pi))
                d = cr.generate_diffraction_pattern(orientation_matrix=M, sigma_excitation_error=SIGMA,
                                                    tol_excitation_error_mult=MULT, tol_intensity=TOL_I, k_max=K_MAX).data
                hkl = np.stack([d['h'], d['k'], d['l']], 1).astype(int)
                xy = np.stack([d['qx'], d['qy']], 1).astype(float)
                inten = np.asarray(d['intensity'], float)
                keep = np.any(hkl != 0, axis=1)
                hkl, xy, inten = hkl[keep], xy[keep], inten[keep]
                keep = inten >= 1e-3 * inten.max()
                hkl, xy, inten = hkl[keep], xy[keep], inten[keep]
                o = np.argsort(-inten, kind='stable')[:STRONGEST]
                hkl, xy, inten = hkl[o], xy[o], inten[o]
                for cond in ('clean', 'scale', 'missing', 'click3', 'click6'):
                    q, noise = xy.copy(), 0.5
                    if cond == 'scale':
                        q = q * 1.02
                    if cond == 'missing':
                        noise = 1.0
                        q = q[rng.random(len(q)) >= 0.3]
                    if cond.startswith('click'):
                        want = int(cond[5:])
                        zero = (hkl @ np.array(uvw) == 0) & (inten >= 0.1 * inten.max())
                        cand = q[zero]
                        r = np.linalg.norm(cand, axis=1)
                        picked = []
                        for i in np.argsort(r, kind='stable'):
                            same_row = any(abs(cand[i, 0] * cand[k, 1] - cand[i, 1] * cand[k, 0]) < 0.02 * r[i] * r[k]
                                           and r[i] > 1.5 * r[k] for k in picked)
                            if not same_row:
                                picked.append(i)
                            if len(picked) >= want:
                                break
                        q = cand[picked]
                    q = q + rng.normal(scale=noise * DQ, size=q.shape)
                    if cond == 'missing':
                        r = rng.uniform(0.25, 1.2, 3)
                        a = rng.uniform(0, 2 * math.pi, 3)
                        q = np.vstack([q, np.stack([r * np.cos(a), r * np.sin(a)], 1)])
                    rec['cases'].append({'zone': list(uvw), 'condition': cond, 'matrix': M.reshape(-1).tolist(),
                                         'scale': 1.02 if cond == 'scale' else 1.0,
                                         'spots': [round(float(v), 7) for v in q.reshape(-1)]})
        out['phases'].append(rec)
        print(ph['name'], len(rec['cases']), flush=True)
    OUT.write_text(json.dumps(out, allow_nan=False, separators=(',', ':')) + '\n')
    print(OUT, OUT.stat().st_size)


if __name__ == '__main__':
    main()
