"""PatternSpots pins (DiffractoKit S4): py4DSTEM's get_maxima_2D on planted spot patterns.

  ../References/python/py4dstem/bin/python -I -B tools/pins/pattern-spots/pins.py

Writes Tests/DiffractoGeometryTests/Fixtures/pattern-spots-pins.json. The oracle is py4DSTEM f050d207 (0.14.19),
`py4DSTEM.preprocess.utils.get_maxima_2D(ar, subpixel="poly", sigma, minAbsoluteIntensity, minRelativeIntensity,
relativeToPeak, minSpacing, edgeBoundary, maxNumPeaks)`, only called. py4DSTEM ships no test for it, so numpy.testing's
assert_allclose float64 default (rtol 1e-7) is the tolerance of the positions and intensities; the maxima and their order
are exact.

The planted images (256 x 256, counts as uint16): kinematic spot patterns from py4DSTEM's own generate_diffraction_pattern
(Al [001], Al [011], Al [112] off by 1 degree, Si [011], Si [001] turned 17 degrees, Mg [2 -1 -1 0]) at 200 kV, drawn as
Gaussian spots (sigma 1.6 px) on a falling background with a direct beam, then Poisson noise (seeded). The truth list is
every drawn spot with its exact centre and its peak height over the background.

Settings sets: "py4dstem_defaults" (get_maxima_2D's own: sigma 0, max 1 peak), "kit_defaults" (PatternSpots.Settings()),
"strict" and "unsmoothed" (sigma 0: integer counts, so ties and NumPy's order are exercised).
"""
import base64
import importlib.metadata as md
import json
from pathlib import Path

import numpy as np
import scipy
from py4DSTEM.preprocess.utils import get_maxima_2D
from py4DSTEM.process.diffraction import Crystal

KIT = Path(__file__).resolve().parents[3]
OUT = KIT / 'Tests/DiffractoGeometryTests/Fixtures/pattern-spots-pins.json'
SEED = 20261125
N = 256
SPOT_SIGMA = 1.6
assert md.version('py4DSTEM') == '0.14.19'

SETTINGS = {
    'py4dstem_defaults': dict(sigma=0, minAbsoluteIntensity=0, minRelativeIntensity=0, relativeToPeak=0, minSpacing=0,
                              edgeBoundary=1, maxNumPeaks=1),
    'kit_defaults': dict(sigma=2, minAbsoluteIntensity=0, minRelativeIntensity=0.005, relativeToPeak=0, minSpacing=10,
                         edgeBoundary=20, maxNumPeaks=70),
    'strict': dict(sigma=1, minAbsoluteIntensity=40, minRelativeIntensity=0.02, relativeToPeak=1, minSpacing=6,
                   edgeBoundary=5, maxNumPeaks=200),
    'unsmoothed': dict(sigma=0, minAbsoluteIntensity=0, minRelativeIntensity=0.003, relativeToPeak=0, minSpacing=4,
                       edgeBoundary=3, maxNumPeaks=300),
}

PHASES = {
    'Al': (np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]]), [13] * 4, [4.0495, 4.0495, 4.0495, 90, 90, 90]),
    'Si': (np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5], [.25, .25, .25], [.75, .75, .25], [.75, .25, .75],
                     [.25, .75, .75]]), [14] * 8, [5.431, 5.431, 5.431, 90, 90, 90]),
    'Mg': (np.array([[1 / 3, 2 / 3, .25], [2 / 3, 1 / 3, .75]]), [12] * 2, [3.209, 3.209, 5.211, 90, 90, 120]),
}

CASES = [('Al', [0, 0, 1], 0.0, 0.0, 58.0), ('Al', [0, 1, 1], 0.4, 0.0, 52.0), ('Al', [1, 1, 2], 1.1, 1.0, 64.0),
         ('Si', [0, 1, 1], -0.3, 0.0, 70.0), ('Si', [0, 0, 1], 17 * np.pi / 180, 0.0, 66.0),
         ('Mg', [2, -1, -1, 0], 0.9, 0.0, 74.0)]


def m3z(g):
    c, s = np.cos(g), np.sin(g)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])


def tilt_x(deg):
    t = np.deg2rad(deg)
    return np.array([[1, 0, 0], [0, np.cos(t), -np.sin(t)], [0, np.sin(t), np.cos(t)]])


def plant(rng, phase, uvw, gamma, tilt_deg, px_per_inv_a):
    pos, num, cell = PHASES[phase]
    cr = Crystal(pos, np.array(num), np.array(cell, dtype=float))
    cr.calculate_structure_factors(k_max=2.2, tol_structure_factor=1e-4)
    cr.setup_diffraction(200e3)
    z = np.array(uvw)
    if z.shape[0] == 4:
        z = cr.hexagonal_to_lattice(z)
    M = cr.parse_orientation(zone_axis_lattice=z) @ m3z(gamma) @ tilt_x(tilt_deg)
    p = cr.generate_diffraction_pattern(orientation_matrix=M, sigma_excitation_error=0.02, k_max=2.2, tol_intensity=1e-4)
    d = p.data
    centre = np.array([N / 2 + rng.uniform(-3, 3), N / 2 + rng.uniform(-3, 3)])
    rows = centre[0] + d['qx'] * px_per_inv_a
    cols = centre[1] + d['qy'] * px_per_inv_a
    inten = d['intensity'] / d['intensity'][(d['h'] != 0) | (d['k'] != 0) | (d['l'] != 0)].max()
    r, c = np.mgrid[0:N, 0:N].astype(float)
    rad = np.hypot(r - centre[0], c - centre[1])
    background = 30 + 400 * np.exp(-rad / 25.0)
    image = background.copy()
    truth = []
    for rr, cc, ii, h, k, l in zip(rows, cols, inten, d['h'], d['k'], d['l']):
        beam = h == 0 and k == 0 and l == 0
        height = 6000.0 if beam else 1800.0 * ii
        if height < 5 or not (2 <= rr < N - 2 and 2 <= cc < N - 2):
            continue
        image += height * np.exp(-((r - rr) ** 2 + (c - cc) ** 2) / (2 * SPOT_SIGMA ** 2))
        truth.append([float(rr), float(cc), float(height), int(h), int(k), int(l)])
    counts = np.minimum(rng.poisson(image), 65535).astype('<u2')
    return counts, centre, truth


def main():
    rng = np.random.default_rng(SEED)
    out = {'generator': 'tools/pins/pattern-spots/pins.py', 'seed': SEED,
           'versions': {'py4DSTEM': md.version('py4DSTEM'), 'numpy': np.__version__, 'scipy': scipy.__version__},
           'commit': 'f050d207851f8accd3aa550365baaf1f38230428', 'size': N, 'spot_sigma': SPOT_SIGMA,
           'call': 'get_maxima_2D(ar.astype(float64), subpixel="poly", **settings); x = row (axis 0), y = column',
           'settings': SETTINGS, 'images': []}
    for phase, uvw, gamma, tilt, scale in CASES:
        counts, centre, truth = plant(rng, phase, uvw, gamma, tilt, scale)
        ar = counts.astype(np.float64)
        maxima = {}
        for name, s in SETTINGS.items():
            m = get_maxima_2D(ar, subpixel='poly', **s)
            maxima[name] = [[float(a), float(b), float(c)] for a, b, c in zip(m['x'], m['y'], m['intensity'])]
        out['images'].append({'label': f'{phase} {uvw} gamma {gamma:.3f} tilt {tilt}', 'centre': centre.tolist(),
                              'counts_u16le_b64': base64.b64encode(counts.tobytes()).decode(), 'truth': truth,
                              'maxima': maxima})
    OUT.write_text(json.dumps(out, allow_nan=False, separators=(',', ':')) + '\n')
    print(OUT, OUT.stat().st_size, [len(i['truth']) for i in out['images']],
          [{k: len(v) for k, v in i['maxima'].items()} for i in out['images']])


if __name__ == '__main__':
    main()
