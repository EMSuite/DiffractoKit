"""Double diffraction between a precipitate and its matrix (DiffractoKit S8, #163): the recorded counts.

  ../References/python/py4dstem/bin/python -I -B tools/pins/double-diffraction/pins.py

No Python package computes double diffraction between two crystals, so this records the closed rule of
`SpotPattern.doubleDiffraction(with:routeWeightCut:)` on py4DSTEM f050d207 (0.14.19) spot patterns, the pinned reference of
`SpotPattern` (P42): `generate_diffraction_pattern(orientation_matrix=M, sigma_excitation_error=0.02,
tol_excitation_error_mult=3, tol_intensity=1e-5, k_max=1.5)` at 200 kV.
  * Matrix: Al (Fm-3m, a 4.0495 A, `CrystalPresets.aluminum`) down [001], [110] and [112]: `parse_orientation(zone_axis_lattice)`,
    no in-plane turn (`CrystalOrientation(zoneAxis:inPlaneAngle: 0)`).
  * Precipitates: beta'' (`CrystalPresets.betaDoublePrime`), planted theta' and T1 (the cells and sites of
    PrecipitateFixtures), each in the first variant of its relationship ([100] || [320], [001] || [1 -3 0]; (001) || (001),
    [100] || [100]; (0001) || (111), [2 1 0] || [1 -1 0]): R by the triad construction of `OrientationRelationship.rotation`,
    the precipitate's orientation R^T M.
  * Rule: every pair (matrix spot a, precipitate spot b) with (I_a / max I_matrix)(I_b / max I_precipitate) >= 0.1, at
    a + b in the pattern plane, unless within 1e-6 A^-1 of the direct beam or of a spot of either pattern; sums within
    1e-6 of each other are one position.
Written per case: the position count, every route (matrix hkl, precipitate hkl, weight) and the number of pairs, on either
side of the cut, whose weight lies within 1e-3 relative of it (py4DSTEM's complex64 structure factors put intensities 5e-5 relative from the kit's, P42).
"""
import importlib.metadata as md
import json
from pathlib import Path

import numpy as np
from py4DSTEM.process.diffraction import Crystal

KIT = Path(__file__).resolve().parents[3]
OUT = KIT / 'Tests/DiffractoGeometryTests/Fixtures/double-diffraction-pins.json'
VOLTS, K_MAX, CUT, EPS = 200e3, 1.5, 0.1, 1e-6
assert md.version('py4DSTEM') == '0.14.19'


def lattice(a, b, c, al=90.0, be=90.0, ga=90.0):
    al, be, ga = np.radians([al, be, ga])
    cx = c * np.cos(be); cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    real = np.array([[a, 0, 0], [b * np.cos(ga), b * np.sin(ga), 0], [cx, cy, np.sqrt(c * c - cx * cx - cy * cy)]])
    return real, np.linalg.inv(real).T


def four_i(z, x, zc):
    return [(z, (x, 0, zc)), (z, (-x, 0, -zc)), (z, (x + .5, .5, zc)), (z, (-x + .5, .5, -zc))]


BETA = [(12, (0, 0, 0)), (12, (.5, .5, 0))] + four_i(12, .3459, .089) + four_i(12, .430, .652) + four_i(14, .0565, .649) \
    + four_i(14, .1885, .224) + four_i(14, .2171, .617)
PHASES = {
    'Al': ((4.0495, 4.0495, 4.0495, 90, 90, 90), [(13, p) for p in [(0, 0, 0), (.5, .5, 0), (.5, 0, .5), (0, .5, .5)]]),
    "beta''": ((15.16, 4.05, 6.74, 90, 105.3, 90), BETA),
    "theta'": ((4.04, 4.04, 5.80, 90, 90, 90), [(29, (0, 0, 0)), (29, (.5, .5, .5)), (13, (0, .5, .25)), (13, (.5, 0, .25)),
                                               (13, (0, .5, .75)), (13, (.5, 0, .75))]),
    'T1': ((4.965, 4.965, 9.345, 90, 90, 120), [(29, (0, 0, 0)), (13, (1 / 3, 2 / 3, .5)), (13, (2 / 3, 1 / 3, .5)), (3, (0, 0, .5))]),
}
RELATIONS = {
    "beta''": ((('direction', [1, 0, 0]), ('direction', [3, 2, 0])), (('direction', [0, 0, 1]), ('direction', [1, -3, 0]))),
    "theta'": ((('plane', [0, 0, 1]), ('plane', [0, 0, 1])), (('direction', [1, 0, 0]), ('direction', [1, 0, 0]))),
    'T1': ((('plane', [0, 0, 1]), ('plane', [1, 1, 1])), (('direction', [2, 1, 0]), ('direction', [1, -1, 0]))),
}
ZONES = [[0, 0, 1], [1, 1, 0], [1, 1, 2]]


def crystal(name):
    cell, sites = PHASES[name]
    cr = Crystal(np.array([s[1] for s in sites], float), np.array([s[0] for s in sites], np.intp), np.array(cell, float))
    cr.calculate_structure_factors(k_max=K_MAX, tol_structure_factor=1e-4)
    cr.setup_diffraction(VOLTS)
    return cr


def vec(rows, kind, idx):
    real, rec = rows
    return np.asarray(idx, float) @ (real if kind == 'direction' else rec)


def triad(a, b):
    e1 = a / np.linalg.norm(a); r = b - (b @ e1) * e1; e2 = r / np.linalg.norm(r)
    return np.column_stack([e1, e2, np.cross(e1, e2)])


def spots(cr, M):
    d = cr.generate_diffraction_pattern(orientation_matrix=M, sigma_excitation_error=0.02, tol_excitation_error_mult=3,
                                        tol_intensity=1e-5, k_max=K_MAX).data
    return [((int(h), int(k), int(l)), np.array([x, y]), float(i)) for h, k, l, x, y, i in
            zip(d['h'], d['k'], d['l'], d['qx'], d['qy'], d['intensity']) if (h, k, l) != (0, 0, 0)]


def double(A, B):
    ma = max(s[2] for s in A); mb = max(s[2] for s in B)
    kin = [s[1] for s in A + B]
    pos, routes, near = [], [], 0
    for ha, pa, ia in A:
        for hb, pb, ib in B:
            w = (ia / ma) * (ib / mb)
            near += abs(w - CUT) <= 1e-3 * CUT
            if w < CUT:
                continue
            p = pa + pb
            if np.linalg.norm(p) <= EPS or any(np.linalg.norm(p - q) <= EPS for q in kin):
                continue
            routes.append([*ha, *hb, w])
            if not any(np.linalg.norm(p - q) <= EPS for q in pos):
                pos.append(p)
    return len(pos), routes, near


def main():
    al = crystal('Al')
    al_rows = lattice(*PHASES['Al'][0][:3])
    out = {'generator': 'tools/pins/double-diffraction/pins.py', 'py4DSTEM': md.version('py4DSTEM'),
           'commit': 'f050d207851f8accd3aa550365baaf1f38230428', 'volts': VOLTS, 'kMax': K_MAX, 'cut': CUT, 'cases': []}
    for name in ("beta''", "theta'", 'T1'):
        cr = crystal(name)
        c = PHASES[name][0]
        rows = lattice(c[0], c[1], c[2], c[3], c[4], c[5])
        (f, s) = RELATIONS[name]
        R = triad(vec(al_rows, *f[1]), vec(al_rows, *s[1])) @ triad(vec(rows, *f[0]), vec(rows, *s[0])).T
        for z in ZONES:
            M = al.parse_orientation(zone_axis_lattice=np.array(z))
            A, B = spots(al, M), spots(cr, R.T @ M)
            n, routes, near = double(A, B)
            out['cases'].append({'precipitate': name, 'zone': z, 'matrixSpots': len(A), 'precipitateSpots': len(B),
                                 'positions': n, 'routes': routes, 'nearCut': near, 'matrix': M.reshape(-1).tolist()})
            print(name, z, len(A), len(B), n, len(routes), near)
    OUT.write_text(json.dumps(out, allow_nan=False, separators=(',', ':')) + '\n')
    print(OUT, OUT.stat().st_size)


if __name__ == '__main__':
    main()
