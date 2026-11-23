"""Tilt refinement planted set (DiffractoKit S5b): spot patterns with intensities at known tilted orientations, from py4DSTEM.

  ../References/python/py4dstem/bin/python -I -B tools/pins/tilt-refinement/planted.py

Writes Tests/DiffractoGeometryTests/Fixtures/tilt-planted.json. No Python function refines the tilt of one pattern on its spot
intensities, so the pin is the planted truth:
  * Al and Si (written below), the first 13 zone-axis classes to index 3 (spglib proper rotations and the sign);
  * per zone axis and tilt 2 orientations: `parse_orientation(zone_axis_lattice=uvw) @ m3z(gamma)` at a seeded gamma, then
    tilted by t = 0, 0.5, 1, 2 or 3 degrees about a seeded in-plane axis (M (R_y R_x)^T with (x, y) = t (cos phi, sin phi),
    the kit's `tilted(by:)`);
  * the pattern: `generate_diffraction_pattern(orientation_matrix=M, sigma_excitation_error=0.02,
    tol_excitation_error_mult=3, tol_intensity=1e-5, k_max=1.5)` at 200 kV, 000 dropped, spots below 1e-3 of the strongest
    dropped, the strongest 40 kept; positions with Gaussian noise of 0.5 px (1 px = 0.00125 1/Angstrom), intensities times
    exp(N(0, 0.2)).
"""
import importlib.metadata as md
import itertools
import json
import math
from pathlib import Path

import numpy as np
import spglib
from scipy.spatial.transform import Rotation as SR
from py4DSTEM.process.diffraction import Crystal

KIT = Path(__file__).resolve().parents[3]
OUT = KIT / 'Tests/DiffractoGeometryTests/Fixtures/tilt-planted.json'
SEED = 20261126
K_MAX, SIGMA, MULT, TOL_I, TOL_F, VOLTS = 1.5, 0.02, 3.0, 1e-5, 1e-4, 200e3
DQ, ZONES, REPEATS, STRONGEST, TILTS = 0.00125, 13, 2, 40, (0.0, 0.5, 1.0, 2.0, 3.0)
assert md.version('py4DSTEM') == '0.14.19'

FCC = [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)]
PHASES = [
    {'name': 'Al', 'cell': dict(a=4.0495, b=4.0495, c=4.0495, alpha=90, beta=90, gamma=90),
     'sites': [{'z': 13, 'frac': list(p), 'occupancy': 1.0} for p in FCC]},
    {'name': 'Si', 'cell': dict(a=5.4309, b=5.4309, c=5.4309, alpha=90, beta=90, gamma=90),
     'sites': [{'z': 14, 'frac': list(p), 'occupancy': 1.0} for p in FCC + [tuple(x + .25 for x in f) for f in FCC]]},
]


def m3z(g):
    c, s = np.cos(g), np.sin(g)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])


def tilted(M, ax, ay):
    rx = SR.from_rotvec([ax, 0, 0]).as_matrix()
    ry = SR.from_rotvec([0, ay, 0]).as_matrix()
    return M @ (ry @ rx).T


def main():
    rng = np.random.default_rng(SEED)
    out = {'generator': 'tools/pins/tilt-refinement/planted.py', 'seed': SEED,
           'versions': {'py4DSTEM': md.version('py4DSTEM'), 'numpy': np.__version__, 'spglib': spglib.__version__},
           'commit': 'f050d207851f8accd3aa550365baaf1f38230428', 'volts': VOLTS, 'pixel': DQ,
           'matrix_layout': 'row-major 9 numbers; columns are lab x, y, z in the crystal Cartesian frame',
           'spots_layout': 'qx, qy, intensity per spot, flat (1/Angstrom, Angstrom^-4)', 'phases': []}
    for ph in PHASES:
        c = ph['cell']
        pos = np.array([s['frac'] for s in ph['sites']], float)
        num = np.array([s['z'] for s in ph['sites']], np.intp)
        cr = Crystal(pos, num, np.array([c['a'], c['b'], c['c'], c['alpha'], c['beta'], c['gamma']]))
        cr.calculate_structure_factors(k_max=K_MAX, tol_structure_factor=TOL_F)
        cr.setup_diffraction(VOLTS)
        rots = [r for r in spglib.get_symmetry((cr.lat_real, pos, num), symprec=1e-3)['rotations'] if round(np.linalg.det(r)) == 1]
        seen, zones = set(), []
        for t in sorted(itertools.product(range(-3, 4), repeat=3), key=lambda t: (sum(x * x for x in t), t)):
            if t == (0, 0, 0) or math.gcd(math.gcd(abs(t[0]), abs(t[1])), abs(t[2])) != 1:
                continue
            key = max(tuple(int(x) for x in s * (r @ np.array(t))) for r in rots for s in (1, -1))
            if key not in seen:
                seen.add(key)
                zones.append(t)
        rec = {'name': ph['name'], 'cell': c, 'sites': ph['sites'], 'cases': []}
        for uvw in zones[:ZONES]:
            for tilt in TILTS:
                for _ in range(REPEATS):
                    M0 = cr.parse_orientation(zone_axis_lattice=np.array(uvw, float)) @ m3z(rng.uniform(0, 2 * math.pi))
                    phi = rng.uniform(0, 2 * math.pi)
                    t = math.radians(tilt)
                    M = tilted(M0, t * math.cos(phi), t * math.sin(phi))
                    d = cr.generate_diffraction_pattern(orientation_matrix=M, sigma_excitation_error=SIGMA,
                                                        tol_excitation_error_mult=MULT, tol_intensity=TOL_I, k_max=K_MAX).data
                    hkl = np.stack([d['h'], d['k'], d['l']], 1).astype(int)
                    xy = np.stack([d['qx'], d['qy']], 1).astype(float)
                    inten = np.asarray(d['intensity'], float)
                    keep = np.any(hkl != 0, axis=1)
                    xy, inten = xy[keep], inten[keep]
                    keep = inten >= 1e-3 * inten.max()
                    xy, inten = xy[keep], inten[keep]
                    o = np.argsort(-inten, kind='stable')[:STRONGEST]
                    xy, inten = xy[o], inten[o]
                    q = xy + rng.normal(scale=0.5 * DQ, size=xy.shape)
                    a = inten * np.exp(rng.normal(scale=0.2, size=len(inten)))
                    spots = np.column_stack([q, a]).reshape(-1)
                    rec['cases'].append({'zone': list(uvw), 'tilt': tilt, 'matrix': M.reshape(-1).tolist(),
                                         'spots': [float('%.8g' % v) for v in spots]})
        out['phases'].append(rec)
        print(ph['name'], len(rec['cases']), flush=True)
    OUT.write_text(json.dumps(out, allow_nan=False, separators=(',', ':')) + '\n')
    print(OUT, OUT.stat().st_size)


if __name__ == '__main__':
    main()
