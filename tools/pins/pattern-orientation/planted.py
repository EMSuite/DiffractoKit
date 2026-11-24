"""Full orientation of one spot pattern (DiffractoKit S6): py4DSTEM-planted patterns and py4DSTEM's own answer.

  ../References/python/py4dstem/bin/python -I -B tools/pins/pattern-orientation/planted.py

Writes Tests/DiffractoOrientationTests/Fixtures/pattern-orientation-planted.json. Three sets, each for Al and Si (written below):
  * zones: the first 13 zone-axis classes to index 3 (spglib proper rotations and the sign), each at a tilt of 0, 1, 2 and
    3 degrees, 2 orientations each: `parse_orientation(zone_axis_lattice=uvw) @ m3z(gamma)` at a seeded gamma, tilted about a
    seeded in-plane axis (M (R_y R_x)^T with (x, y) = t (cos phi, sin phi), the kit's `tilting(_:by:)`);
  * random: the spike's 40 random orientations per crystal (part6c_random.py: default_rng(9), scipy `Rotation.random`), the
    three its combined route got wrong marked `spike_wrong`;
  * edge: DiffractoKit #65, Al [001] and [011] at tilts 0, 0.5, 1 and 2 degrees, 3 orientations each, the spots cut at
    |q| <= 0.75 1/Angstrom (the AlMgSi demo's detector edge, 63.5 px x 0.012) while the template library reaches 1.5;
  * the pattern: `generate_diffraction_pattern(orientation_matrix=M, sigma_excitation_error=0.02)` at 200 kV, k_max 1.5,
    000 dropped, spots below 1e-3 of the strongest dropped, the strongest 40 kept; positions with Gaussian noise of 0.5 px
    (1 px = 0.00125 1/Angstrom), intensities times exp(N(0, 0.2));
  * py4DSTEM's answer: `orientation_plan(accel_voltage=200e3)` at its defaults (2 degree steps, corr_kernel_size 0.08,
    power_intensity 0.25) and `match_single_pattern(num_matches_return=5)`, its matrices and correlations.
"""
import importlib.metadata as md
import itertools
import json
import math
import warnings
from pathlib import Path

import numpy as np
import spglib
from scipy.spatial.transform import Rotation as SR

warnings.filterwarnings('ignore')
import py4DSTEM  # noqa: E402
from py4DSTEM.process.diffraction import Crystal  # noqa: E402

KIT = Path(__file__).resolve().parents[3]
OUT = KIT / 'Tests/DiffractoOrientationTests/Fixtures/pattern-orientation-planted.json'
SEED = 20261127
K_MAX, SIGMA, VOLTS, DQ, STRONGEST = 1.5, 0.02, 200e3, 0.00125, 40
ZONES, REPEATS, TILTS = 13, 2, (0.0, 1.0, 2.0, 3.0)
EDGE, EDGE_TILTS, EDGE_REPEATS = 0.75, (0.0, 0.5, 1.0, 2.0), 3
MATCHES = 5
# the random orientations the spike's combined route got wrong (part6f.json: Al 16 and 20, Si 2; two leave a measured spot
# unexplained, one is a tie of a sparse pattern)
SPIKE_WRONG = {'Al': {16, 20}, 'Si': {2}}
assert md.version('py4DSTEM') == '0.14.19'

FCC = [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)]
PHASES = [
    {'name': 'Al', 'a': 4.0495, 'z': 13, 'frac': FCC},
    {'name': 'Si', 'a': 5.4309, 'z': 14, 'frac': FCC + [tuple(x + .25 for x in f) for f in FCC]},
]


def m3z(g):
    c, s = np.cos(g), np.sin(g)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])


def tilted(M, ax, ay):
    rx = SR.from_rotvec([ax, 0, 0]).as_matrix()
    ry = SR.from_rotvec([0, ay, 0]).as_matrix()
    return M @ (ry @ rx).T


def pattern(cr, M, rng, cut=None):
    """The spike's pattern rule (part6a/part6c): 000 and spots below 1e-3 of the strongest dropped, the strongest 40."""
    d = cr.generate_diffraction_pattern(orientation_matrix=M, sigma_excitation_error=SIGMA).data
    r = np.hypot(d['qx'], d['qy'])
    keep = r > 1e-6
    if cut is not None:
        keep &= r <= cut
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


def main():
    out = {'generator': 'tools/pins/pattern-orientation/planted.py', 'seed': SEED,
           'versions': {'py4DSTEM': md.version('py4DSTEM'), 'numpy': np.__version__, 'spglib': spglib.__version__},
           'commit': 'f050d207851f8accd3aa550365baaf1f38230428', 'volts': VOLTS, 'pixel': DQ, 'kmax': K_MAX,
           'matrix_layout': 'row-major 9 numbers; columns are lab x, y, z in the crystal Cartesian frame (py4DSTEM)',
           'spots_layout': 'qx, qy, intensity per spot, flat (1/Angstrom, Angstrom^-4)', 'phases': []}
    rng = np.random.default_rng(SEED)
    for ph in PHASES:
        pos = np.array(ph['frac'], float)
        num = np.full(len(pos), ph['z'], np.intp)
        cr = Crystal(pos, num, ph['a'])
        cr.calculate_structure_factors(K_MAX)
        cr.orientation_plan(accel_voltage=VOLTS, progress_bar=False)
        rots = [r for r in spglib.get_symmetry((cr.lat_real, pos, num), symprec=1e-3)['rotations'] if round(np.linalg.det(r)) == 1]
        seen, zones = set(), []
        for t in sorted(itertools.product(range(-3, 4), repeat=3), key=lambda t: (sum(x * x for x in t), t)):
            if t == (0, 0, 0) or math.gcd(math.gcd(abs(t[0]), abs(t[1])), abs(t[2])) != 1:
                continue
            key = max(tuple(int(x) for x in s * (r @ np.array(t))) for r in rots for s in (1, -1))
            if key not in seen:
                seen.add(key)
                zones.append(t)
        rec = {'name': ph['name'], 'a': ph['a'], 'z': ph['z'], 'frac': [list(f) for f in ph['frac']], 'cases': []}

        def planted(uvw, tilt, cut, set_name):
            while True:
                M0 = cr.parse_orientation(zone_axis_lattice=np.array(uvw, float)) @ m3z(rng.uniform(0, 2 * math.pi))
                phi = rng.uniform(0, 2 * math.pi)
                t = math.radians(tilt)
                M = tilted(M0, t * math.cos(phi), t * math.sin(phi))
                data = pattern(cr, M, rng, cut)
                if data is not None:
                    rec['cases'].append(case(cr, M, data, set=set_name, zone=list(uvw), tilt=tilt))
                    return

        for uvw in zones[:ZONES]:
            for tilt in TILTS:
                for _ in range(REPEATS):
                    planted(uvw, tilt, None, 'zones')
        if ph['name'] == 'Al':
            for uvw in ((0, 0, 1), (0, 1, 1)):
                for tilt in EDGE_TILTS:
                    for _ in range(EDGE_REPEATS):
                        planted(uvw, tilt, EDGE, 'edge')
        out['phases'].append(rec)

    # the spike's random set, in its own order and with its own generator (both crystals, as part6c_random.py)
    spike = np.random.default_rng(9)
    for rec in out['phases']:
        pos = np.array(rec['frac'], float)
        cr = Crystal(pos, np.full(len(pos), rec['z'], np.intp), rec['a'])
        cr.calculate_structure_factors(K_MAX)
        cr.orientation_plan(accel_voltage=VOLTS, progress_bar=False)
        n = 0
        while n < 40:
            Mrows = SR.random(random_state=spike.integers(1 << 31)).as_matrix()
            bp = cr.generate_diffraction_pattern(orientation_matrix=Mrows.T, sigma_excitation_error=SIGMA)
            r = np.hypot(bp.data['qx'], bp.data['qy'])
            if (r > 1e-6).sum() < 4:
                continue
            keep = (r > 1e-6) & (bp.data['intensity'] >= 1e-3 * bp.data['intensity'][r > 1e-6].max())
            data = bp.data[keep]
            data = data[np.argsort(-data['intensity'])[:STRONGEST]].copy()
            if len(data) < 4:
                continue
            data['qx'] += spike.normal(scale=0.5 * DQ, size=len(data))
            data['qy'] += spike.normal(scale=0.5 * DQ, size=len(data))
            data['intensity'] *= np.exp(spike.normal(scale=0.2, size=len(data)))
            rec['cases'].append(case(cr, Mrows.T, data, set='random', zone=None, tilt=None,
                                     spike_wrong=n in SPIKE_WRONG[rec['name']]))
            n += 1
        print(rec['name'], {s: sum(c['set'] == s for c in rec['cases']) for s in ('zones', 'edge', 'random')}, flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, allow_nan=False, separators=(',', ':')) + '\n')
    print(OUT, OUT.stat().st_size)


if __name__ == '__main__':
    main()
