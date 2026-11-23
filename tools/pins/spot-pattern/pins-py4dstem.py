"""SpotPattern pins (DiffractoKit S2): the kinematic spot pattern of a structure at an orientation, from py4DSTEM.

  ../References/python/py4dstem/bin/python -I -B tools/pins/spot-pattern/pins-py4dstem.py

Writes Tests/DiffractoGeometryTests/Fixtures/spot-pattern-pins.json. The oracle is py4DSTEM f050d207 (0.14.19), only called:
  * the 7 phases of tools/pins/phase-vectors (cell, sites; read from that fixture), `Crystal.calculate_structure_factors(
    k_max=1.5, tol_structure_factor=1e-4)`;
  * at 80, 200 and 300 kV (`setup_diffraction`), per phase: its first three zone axes (`parse_orientation`) turned by
    0.4 rad (`@ m3z`), the same three tilted by 1 degree about lab x and -2 degrees about lab y (M (R_y R_x)^T, the kit's
    `tilted(by:)`), and three seeded random rotations;
  * per orientation `generate_diffraction_pattern(orientation_matrix=M, sigma_excitation_error=0.02,
    tol_excitation_error_mult=3, tol_intensity=1e-5, k_max=1.5)` (h, k, l, qx, qy, intensity; 000 dropped) and, for those
    reflections, `excitation_errors(M.T @ g)`.
py4DSTEM keeps its structure factors in complex64 and ships no tests: positions at assert_allclose's float64 default
(rtol 1e-7), intensities at 5e-5 relative (the float32 rounding of F, plan row S2).
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
OUT = KIT / 'Tests/DiffractoGeometryTests/Fixtures/spot-pattern-pins.json'
SEED = 20261124
K_MAX = 1.5
SIGMA, MULT, TOL_I, TOL_F = 0.02, 3.0, 1e-5, 1e-4
VOLTS = [80e3, 200e3, 300e3]
GAMMA = 0.4
TILT = (np.deg2rad(1.0), np.deg2rad(-2.0))
assert md.version('py4DSTEM') == '0.14.19'


def m3z(g):
    c, s = np.cos(g), np.sin(g)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])


def tilted(M, ax, ay):
    rx = SR.from_rotvec([ax, 0, 0]).as_matrix()
    ry = SR.from_rotvec([0, ay, 0]).as_matrix()
    return M @ (ry @ rx).T


def crystal(ph):
    c = ph['cell']
    pos = np.array([s['frac'] for s in ph['sites']], dtype=np.float64)
    num = np.array([s['z'] for s in ph['sites']], dtype=np.intp)
    occ = np.array([s['occupancy'] for s in ph['sites']], dtype=np.float64)
    cr = Crystal(pos, num, np.array([c['a'], c['b'], c['c'], c['alpha'], c['beta'], c['gamma']]), occupancy=occ)
    cr.calculate_structure_factors(k_max=K_MAX, tol_structure_factor=TOL_F)
    return cr


def pattern(cr, M):
    p = cr.generate_diffraction_pattern(orientation_matrix=M, sigma_excitation_error=SIGMA, tol_excitation_error_mult=MULT,
                                        tol_intensity=TOL_I, k_max=K_MAX)
    d = p.data
    sg_all = cr.excitation_errors(M.T @ cr.g_vec_all)
    index = {tuple(int(x) for x in h): i for i, h in enumerate(cr.hkl.T)}
    rows = []
    for h, k, l, x, y, i in zip(d['h'], d['k'], d['l'], d['qx'], d['qy'], d['intensity']):
        if h == 0 and k == 0 and l == 0:
            continue
        j = index[(int(h), int(k), int(l))]
        rows.append([int(h), int(k), int(l), float(x), float(y), float(i), float(sg_all[j])])
    return rows


def main():
    src = json.load(open(SRC))
    rng = np.random.default_rng(SEED)
    out = {'generator': 'tools/pins/spot-pattern/pins-py4dstem.py', 'seed': SEED,
           'versions': {'py4DSTEM': md.version('py4DSTEM'), 'numpy': np.__version__, 'scipy': scipy.__version__},
           'commit': 'f050d207851f8accd3aa550365baaf1f38230428',
           'settings': {'k_max': K_MAX, 'sigma_excitation_error': SIGMA, 'tol_excitation_error_mult': MULT,
                        'tol_intensity': TOL_I, 'tol_structure_factor': TOL_F},
           'matrix_layout': 'row-major 9 numbers; columns are lab x, y, z in the crystal Cartesian frame',
           'row_layout': 'h, k, l, qx, qy, intensity, excitation error', 'phases': []}
    for ph in src['phases']:
        cr = crystal(ph)
        mats = []
        for z in ph['zone_axes'][:3]:
            M = cr.parse_orientation(zone_axis_lattice=np.array(z['uvw'])) @ m3z(GAMMA)
            mats.append(('zone %s' % z['uvw'], M))
            mats.append(('zone %s tilted' % z['uvw'], tilted(M, *TILT)))
        for M in SR.random(3, random_state=rng.integers(2**31)).as_matrix():
            mats.append(('random', M))
        rec = {'name': ph['name'], 'cell': ph['cell'], 'sites': ph['sites'], 'patterns': []}
        for volts in VOLTS:
            cr.setup_diffraction(volts)
            for label, M in mats:
                rec['patterns'].append({'label': label, 'volts': volts, 'wavelength': cr.wavelength,
                                        'matrix': M.reshape(-1).tolist(), 'spots': pattern(cr, M)})
        out['phases'].append(rec)
    OUT.write_text(json.dumps(out, allow_nan=False, separators=(',', ':')) + '\n')
    n = sum(len(p['spots']) for ph in out['phases'] for p in ph['patterns'])
    print(OUT, OUT.stat().st_size, n)


if __name__ == '__main__':
    main()
