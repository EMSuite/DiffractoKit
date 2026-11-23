"""CrystalOrientation pins (DiffractoKit S1): the orientation matrix of a zone axis and where g appears, from py4DSTEM.

  ../References/python/py4dstem/bin/python -I -B tools/pins/crystal-orientation/pins-py4dstem.py

Writes Tests/DiffractoCoreTests/Fixtures/crystal-orientation-pins.json. The oracle is py4DSTEM f050d207 (0.14.19), only called:
  * the 7 phases and 37 zone axes of tools/pins/phase-vectors (read from that fixture: cell, sites, uvw);
  * per axis `Crystal.parse_orientation(zone_axis_lattice=uvw)` (columns x', y', z' in the crystal's Cartesian frame);
  * per axis, at in-plane angles 0, 0.7 and 2.5 rad, the matrix `parse_orientation(uvw) @ m3z(gamma)` (crystal_ACOM.py:
    m3z = [[c, s, 0], [-s, c, 0], [0, 0, 1]]) and the peaks `generate_diffraction_pattern(orientation_matrix=...)` at 200 kV
    (h, k, l, qx, qy): where g appears for that matrix;
  * 20 seeded random rotations (SciPy) of the first and last phase with their peaks the same way.
py4DSTEM states no tolerance for these functions (it ships no tests): numpy.testing.assert_allclose's float64 defaults apply.
"""
import importlib.metadata as md
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.spatial.transform import Rotation as SR
from py4DSTEM.process.diffraction import Crystal

KIT = Path(__file__).resolve().parents[3]
SRC = KIT / 'Tests/DiffractoOrientationTests/Fixtures/py4dstem-phase-vector-pins.json'
OUT = KIT / 'Tests/DiffractoCoreTests/Fixtures/crystal-orientation-pins.json'
SEED = 20261123
GAMMAS = [0.0, 0.7, 2.5]
K_MAX = 1.2
assert md.version('py4DSTEM') == '0.14.19'


def m3z(g):
    c, s = np.cos(g), np.sin(g)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])


def crystal(ph):
    c = ph['cell']
    pos = np.array([s['frac'] for s in ph['sites']], dtype=np.float64)
    num = np.array([s['z'] for s in ph['sites']], dtype=np.intp)
    occ = np.array([s['occupancy'] for s in ph['sites']], dtype=np.float64)
    cr = Crystal(pos, num, np.array([c['a'], c['b'], c['c'], c['alpha'], c['beta'], c['gamma']]), occupancy=occ)
    cr.calculate_structure_factors(k_max=K_MAX, tol_structure_factor=1e-4)
    cr.setup_diffraction(200e3)
    return cr


def peaks(cr, M):
    p = cr.generate_diffraction_pattern(orientation_matrix=M, k_max=K_MAX, tol_intensity=1e-6)
    d = p.data
    keep = (d['h'] != 0) | (d['k'] != 0) | (d['l'] != 0)
    return [[int(h), int(k), int(l), float(x), float(y)]
            for h, k, l, x, y in zip(d['h'][keep], d['k'][keep], d['l'][keep], d['qx'][keep], d['qy'][keep])]


def main():
    src = json.load(open(SRC))
    rng = np.random.default_rng(SEED)
    out = {'generator': 'tools/pins/crystal-orientation/pins-py4dstem.py', 'seed': SEED,
           'versions': {'py4DSTEM': md.version('py4DSTEM'), 'numpy': np.__version__, 'scipy': scipy.__version__},
           'commit': 'f050d207851f8accd3aa550365baaf1f38230428', 'k_max': K_MAX, 'voltage_V': 200e3,
           'm3z': '[[c, s, 0], [-s, c, 0], [0, 0, 1]] on the right of parse_orientation (crystal_ACOM.py)',
           'matrix_layout': 'row-major 9 numbers of the 3x3 matrix whose columns are lab x, y, z in the crystal Cartesian frame',
           'phases': [], 'random': []}
    for ph in src['phases']:
        cr = crystal(ph)
        rec = {'name': ph['name'], 'cell': ph['cell'], 'lat_real': cr.lat_real.tolist(), 'axes': []}
        for z in ph['zone_axes']:
            uvw = [int(x) for x in z['uvw']]
            O = cr.parse_orientation(zone_axis_lattice=np.array(uvw))
            rows = []
            for g in GAMMAS:
                M = O @ m3z(g)
                rows.append({'gamma': g, 'matrix': M.reshape(-1).tolist(), 'peaks': peaks(cr, M)})
            rec['axes'].append({'uvw': uvw, 'matrix': O.reshape(-1).tolist(), 'in_plane': rows})
        out['phases'].append(rec)
    for ph in (src['phases'][0], src['phases'][-1]):
        cr = crystal(ph)
        for M in SR.random(10, random_state=rng.integers(2**31)).as_matrix():
            out['random'].append({'phase': ph['name'], 'matrix': M.reshape(-1).tolist(), 'peaks': peaks(cr, M)})
    OUT.write_text(json.dumps(out, allow_nan=False, separators=(',', ':')) + '\n')
    print(OUT, OUT.stat().st_size, sum(len(p['axes']) for p in out['phases']))


if __name__ == '__main__':
    main()
