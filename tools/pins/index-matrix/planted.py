"""The Index reference matrix: `PatternOrientation.index` on py4DSTEM-planted spot patterns of eight crystals of five
systems, with py4DSTEM's own answer beside each.

  ../References/python/py4dstem/bin/python -I -B tools/pins/index-matrix/planted.py

Writes Tests/DiffractoOrientationTests/Fixtures/index-matrix-planted.json. It widens tools/pins/pattern-orientation
(Al and Si) to crystals the suite's users index: fcc Au, rock-salt MgO, bcc Fe, perovskite SrTiO3, hcp Mg and Zr,
rutile TiO2 and 2H-WS2. The pattern rule and py4DSTEM's call are that script's, unchanged:
  * per phase the first 8 zone-axis classes with |u|, |v|, |w| <= 3 (spglib proper rotations and the sign), each at a
    tilt of 0, 1 and 2 degrees, 2 orientations each (`parse_orientation(zone_axis_lattice=uvw) @ m3z(gamma)`, tilted about
    a seeded in-plane axis), and 10 seeded random orientations (scipy `Rotation.random`);
  * the pattern: `generate_diffraction_pattern(orientation_matrix=M, sigma_excitation_error=0.02)` at 200 kV, k_max 1.5,
    000 dropped, spots below 1e-3 of the strongest dropped, the strongest 40 kept; positions with Gaussian noise of 0.5 px
    (1 px = 0.00125 1/Angstrom), intensities times exp(N(0, 0.2));
  * py4DSTEM's answer: `orientation_plan(accel_voltage=200e3)` at its defaults (2 degree steps) and
    `match_single_pattern(num_matches_return=5)`.
Cells follow py4DSTEM's convention (a along x, b in the xy plane), which is the kit's `UnitCell`.
"""
import importlib.metadata as md
import itertools
import json
import math
import time
import warnings
from pathlib import Path

import numpy as np
import spglib
from scipy.spatial.transform import Rotation as SR

warnings.filterwarnings('ignore')
import py4DSTEM  # noqa: E402
from py4DSTEM.process.diffraction import Crystal  # noqa: E402

KIT = Path(__file__).resolve().parents[3]
OUT = KIT / 'Tests/DiffractoOrientationTests/Fixtures/index-matrix-planted.json'
SEED = 20261127
K_MAX, SIGMA, VOLTS, DQ, STRONGEST = 1.5, 0.02, 200e3, 0.00125, 40
ZONES, REPEATS, TILTS, RANDOM = 8, 2, (0.0, 1.0, 2.0), 10
MATCHES = 5
assert md.version('py4DSTEM') == '0.14.19'

FCC = [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)]
HCP = [(1 / 3, 2 / 3, .25), (2 / 3, 1 / 3, .75)]
U = 0.3053      # rutile O
ZS = 0.621      # 2H-WS2 S (4f)
PHASES = [
    {'name': 'Au', 'cell': [4.0782, 4.0782, 4.0782, 90, 90, 90], 'sites': [(79, f) for f in FCC]},
    {'name': 'MgO', 'cell': [4.2117, 4.2117, 4.2117, 90, 90, 90],
     'sites': [(12, f) for f in FCC] + [(8, ((f[0] + .5) % 1, f[1], f[2])) for f in FCC]},
    {'name': 'Fe', 'cell': [2.8665, 2.8665, 2.8665, 90, 90, 90], 'sites': [(26, (0, 0, 0)), (26, (.5, .5, .5))]},
    {'name': 'SrTiO3', 'cell': [3.905, 3.905, 3.905, 90, 90, 90],
     'sites': [(38, (0, 0, 0)), (22, (.5, .5, .5)), (8, (.5, .5, 0)), (8, (.5, 0, .5)), (8, (0, .5, .5))]},
    {'name': 'Mg', 'cell': [3.2094, 3.2094, 5.2105, 90, 90, 120], 'sites': [(12, f) for f in HCP]},
    {'name': 'Zr', 'cell': [3.232, 3.232, 5.147, 90, 90, 120], 'sites': [(40, f) for f in HCP]},
    {'name': 'TiO2', 'cell': [4.594, 4.594, 2.959, 90, 90, 90],
     'sites': [(22, (0, 0, 0)), (22, (.5, .5, .5)), (8, (U, U, 0)), (8, (1 - U, 1 - U, 0)),
               (8, (.5 + U, .5 - U, .5)), (8, (.5 - U, .5 + U, .5))]},
    {'name': 'WS2', 'cell': [3.153, 3.153, 12.323, 90, 90, 120],
     'sites': [(74, HCP[0]), (74, HCP[1]), (16, (1 / 3, 2 / 3, ZS)), (16, (2 / 3, 1 / 3, ZS - .5)),
               (16, (2 / 3, 1 / 3, 1 - ZS)), (16, (1 / 3, 2 / 3, 1.5 - ZS))]},
]


def m3z(g):
    c, s = np.cos(g), np.sin(g)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])


def tilted(M, ax, ay):
    rx = SR.from_rotvec([ax, 0, 0]).as_matrix()
    ry = SR.from_rotvec([0, ay, 0]).as_matrix()
    return M @ (ry @ rx).T


def pattern(cr, M, rng):
    d = cr.generate_diffraction_pattern(orientation_matrix=M, sigma_excitation_error=SIGMA).data
    r = np.hypot(d['qx'], d['qy'])
    keep = r > 1e-6
    if keep.sum() < 2:
        return None
    keep &= d['intensity'] >= 1e-3 * d['intensity'][keep].max()
    data = d[keep]
    data = data[np.argsort(-data['intensity'], kind='stable')[:STRONGEST]].copy()
    if len(data) < 4:
        return None
    data['qx'] += rng.normal(scale=0.5 * DQ, size=len(data))
    data['qy'] += rng.normal(scale=0.5 * DQ, size=len(data))
    data['intensity'] *= np.exp(rng.normal(scale=0.2, size=len(data)))
    return data


def case(cr, M, data, **extra):
    ori = cr.match_single_pattern(py4DSTEM.PointList(data=data), num_matches_return=MATCHES, verbose=False)
    spots = np.column_stack([data['qx'], data['qy'], data['intensity']]).reshape(-1)
    return dict(extra, matrix=np.asarray(M).reshape(-1).tolist(), spots=[float('%.8g' % v) for v in spots],
                py4dstem=[np.asarray(m).reshape(-1).tolist() for m in np.asarray(ori.matrix)[:MATCHES]],
                py4dstem_corr=[float(c) for c in np.asarray(ori.corr)[:MATCHES]])


def zone_classes(cr, pos, num):
    rots = [r for r in spglib.get_symmetry((cr.lat_real, pos, num), symprec=1e-3)['rotations'] if round(np.linalg.det(r)) == 1]
    seen, zones = set(), []
    for t in sorted(itertools.product(range(-3, 4), repeat=3), key=lambda t: (sum(x * x for x in t), t)):
        if t == (0, 0, 0) or math.gcd(math.gcd(abs(t[0]), abs(t[1])), abs(t[2])) != 1:
            continue
        key = max(tuple(int(x) for x in s * (r @ np.array(t))) for r in rots for s in (1, -1))
        if key not in seen:
            seen.add(key)
            zones.append(t)
    return zones, len(rots)


def main():
    out = {'generator': 'tools/pins/index-matrix/planted.py', 'seed': SEED,
           'versions': {'py4DSTEM': md.version('py4DSTEM'), 'numpy': np.__version__, 'spglib': spglib.__version__},
           'commit': 'f050d207851f8accd3aa550365baaf1f38230428', 'volts': VOLTS, 'pixel': DQ, 'kmax': K_MAX,
           'matrix_layout': 'row-major 9 numbers; columns are lab x, y, z in the crystal Cartesian frame (py4DSTEM)',
           'spots_layout': 'qx, qy, intensity per spot, flat (1/Angstrom, Angstrom^-4)', 'phases': []}
    rng = np.random.default_rng(SEED)
    for ph in PHASES:
        t0 = time.time()
        pos = np.array([s[1] for s in ph['sites']], float) % 1
        num = np.array([s[0] for s in ph['sites']], np.intp)
        cr = Crystal(pos, num, ph['cell'])
        cr.calculate_structure_factors(K_MAX)
        cr.orientation_plan(accel_voltage=VOLTS, progress_bar=False)
        zones, nrot = zone_classes(cr, pos, num)
        rec = {'name': ph['name'], 'cell': ph['cell'], 'z': num.tolist(), 'frac': pos.tolist(),
               'proper_rotations': nrot, 'cases': []}
        for uvw in zones[:ZONES]:
            for tilt in TILTS:
                for _ in range(REPEATS):
                    while True:
                        M0 = cr.parse_orientation(zone_axis_lattice=np.array(uvw, float)) @ m3z(rng.uniform(0, 2 * math.pi))
                        phi = rng.uniform(0, 2 * math.pi)
                        t = math.radians(tilt)
                        M = tilted(M0, t * math.cos(phi), t * math.sin(phi))
                        data = pattern(cr, M, rng)
                        if data is not None:
                            rec['cases'].append(case(cr, M, data, set='zones', zone=list(uvw), tilt=tilt))
                            break
        n = 0
        while n < RANDOM:
            M = SR.random(random_state=int(rng.integers(1 << 31))).as_matrix().T
            data = pattern(cr, M, rng)
            if data is None:
                continue
            rec['cases'].append(case(cr, M, data, set='random', zone=None, tilt=None))
            n += 1
        out['phases'].append(rec)
        print(ph['name'], len(rec['cases']), 'cases', f'{time.time() - t0:.0f} s', 'zones', zones[:ZONES], flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, allow_nan=False, separators=(',', ':')) + '\n')
    print(OUT, OUT.stat().st_size)


if __name__ == '__main__':
    main()
