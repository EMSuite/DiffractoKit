"""Index on recorded data: py4DSTEM's spots and its own orientation match on patterns of a real 4D-STEM scan, so
`PatternOrientation.index` is held on measured spots (detector noise, a twin, an elliptical projector distortion), not
only on planted ones (tools/pins/index-matrix).

  ../References/python/py4dstem/bin/python -I -B tools/pins/index-real/nanowire.py

The data: pyxem's demo dataset, Zenodo 14113591 (Carter Francis, 2024, CC-BY-4.0), `twinned_nanowire.hdf5`: a
precession (2 deg) scan of a twinned GaAs nanowire seen down <110>, 100 x 30 positions, 144 x 144 uint8 at 300 kV.
Downloaded once to $INDEX_REAL_DATA (default ../References/CuboKit/zenodo, CuboKit's tools/bragg-real cache) and
checked against Zenodo's MD5. Only spot lists derived from it are committed.

The recipe, step by step as its users run it:
  1. the disks: CuboKit's tools/bragg-real recipe (py4DSTEM tutorial orient_strain_01_WS2: probe radius of the mean
     pattern, synthetic probe width 0.7, sigmoid kernel radii 1.2 r and 4 r, `find_Bragg_disks` with the tutorial's
     keywords, minAbsoluteIntensity 8 and minPeakSpacing 5), on the whole scan; `measure_origin`, `fit_origin`;
  2. the distortion: pyxem-demos notebook 01 corrects this dataset with
     `apply_affine_transformation([[0.99, 0, 0], [0, 0.69, 0], [0, 0, 1]])` (skimage x = column, y = row, about the
     pattern centre); on the spots that is the row offset from the origin times 0.69 and the column offset times 0.99;
  3. the calibration: pyxem-demos' recip_cal (d111 of GaAs, a = 5.6535 Angstrom, at 11.4 px) as the start, then
     py4DSTEM's `Crystal.calibrate_pixel_size(bragg_k_power=2)` against GaAs;
  4. py4DSTEM's answer: `orientation_plan(accel_voltage=300e3)` at its defaults (2 degree steps) and
     `match_single_pattern(num_matches_return=5)` on every sixth scan row in three columns of the wire, patterns with
     at least 6 spots besides the direct beam.

Writes Tests/DiffractoOrientationTests/Fixtures/index-real-nanowire.json: per case the scan position, the spots (qx, qy
in 1/Angstrom about the origin, intensity; the direct beam dropped) and py4DSTEM's five matrices and correlations.
"""
import hashlib
import importlib.metadata as md
import json
import os
import platform
import urllib.request
import warnings
from pathlib import Path

import h5py
import numpy as np

warnings.filterwarnings('ignore')
os.environ.setdefault('TQDM_DISABLE', '1')
import py4DSTEM  # noqa: E402
from emdfile import PointList  # noqa: E402
from py4DSTEM.process.calibration import get_probe_size  # noqa: E402
from py4DSTEM.process.diffraction import Crystal  # noqa: E402

assert md.version('py4DSTEM') == '0.14.19'

KIT = Path(__file__).resolve().parents[3]
OUT = KIT / 'Tests/DiffractoOrientationTests/Fixtures/index-real-nanowire.json'
DATA = Path(os.environ.get('INDEX_REAL_DATA', KIT.parent / 'References' / 'CuboKit' / 'zenodo'))
SOURCE = {'record': 14113591, 'file': 'twinned_nanowire.hdf5', 'md5': '2765fa855db8bc011252dbb425facc81'}
DISKS = dict(corrPower=1.0, sigma_dp=0.0, sigma_cc=0.0, edgeBoundary=4, minRelativeIntensity=0.0, relativeToPeak=0,
             subpixel='poly', upsample_factor=8, maxNumPeaks=1000, minAbsoluteIntensity=8.0, minPeakSpacing=5.0)
AFFINE_ROW, AFFINE_COLUMN = 0.69, 0.99
A_GAAS = 5.6535
RECIP_CAL = np.sqrt(3 / A_GAAS ** 2) / 11.4
VOLTS, MATCHES, MIN_SPOTS = 300e3, 5, 6
POSITIONS = [(i, j) for i in range(2, 100, 6) for j in (8, 15, 22)]
FCC = [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)]
SITES = [(31, f) for f in FCC] + [(33, tuple(x + .25 for x in f)) for f in FCC]


def fetch():
    path = DATA / SOURCE['file']
    if not path.exists():
        DATA.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(f"https://zenodo.org/records/{SOURCE['record']}/files/{SOURCE['file']}?download=1", path)
    digest = hashlib.md5(path.read_bytes()).hexdigest()
    assert digest == SOURCE['md5'], f'{path}: md5 {digest}, Zenodo says {SOURCE["md5"]}'
    return path


def main():
    with h5py.File(fetch(), 'r') as h:
        data = np.asarray(h['Experiments/__unnamed__/data'])
    cube = py4DSTEM.DataCube(data=data)
    dp_mean = data.astype(np.float64).mean(axis=(0, 1))
    r = float(get_probe_size(dp_mean)[0])
    q = cube.Qshape
    probe = py4DSTEM.Probe.generate_synthetic_probe(radius=r, width=0.7, Qshape=q)
    kernel = probe.get_kernel(mode='sigmoid', origin=(q[0] / 2, q[0] / 2), radii=(r * 1.2, r * 4))
    bp = cube.find_Bragg_disks(template=kernel, **DISKS)
    bp.measure_origin()
    bp.fit_origin()
    ox, oy = bp.calibration.get_origin()
    for i in range(bp.shape[0]):
        for j in range(bp.shape[1]):
            p = bp.raw[i, j].data
            p['qx'] = ox[i, j] + (p['qx'] - ox[i, j]) * AFFINE_ROW
            p['qy'] = oy[i, j] + (p['qy'] - oy[i, j]) * AFFINE_COLUMN
    bp.calibration.set_Q_pixel_size(RECIP_CAL)
    bp.calibration.set_Q_pixel_units('A^-1')
    cr = Crystal(positions=[s[1] for s in SITES], numbers=[s[0] for s in SITES], cell=A_GAAS)
    cr.calculate_structure_factors(1.5)
    bp = cr.calibrate_pixel_size(bragg_peaks=bp, bragg_k_power=2.0, plot_result=False)
    pixel = float(bp.calibration.get_Q_pixel_size())
    cr.orientation_plan(accel_voltage=VOLTS, progress_bar=False)
    cases = []
    for i, j in POSITIONS:
        d = bp.cal[i, j].data
        keep = np.hypot(d['qx'], d['qy']) > 1e-3
        if keep.sum() < MIN_SPOTS:
            continue
        pl = PointList(d[keep].copy())
        o = cr.match_single_pattern(pl, num_matches_return=MATCHES, verbose=False)
        spots = np.column_stack([d['qx'][keep], d['qy'][keep], d['intensity'][keep]]).reshape(-1)
        cases.append(dict(position=[i, j], spots=[float('%.10g' % v) for v in spots],
                          py4dstem=[np.asarray(m).reshape(-1).tolist() for m in np.asarray(o.matrix)[:MATCHES]],
                          py4dstem_corr=[float(c) for c in np.asarray(o.corr)[:MATCHES]]))
        print(i, j, int(keep.sum()), np.round(o.corr, 2), flush=True)
    out = {'generator': 'tools/pins/index-real/nanowire.py',
           'source': dict(SOURCE, doi='10.5281/zenodo.14113591', licence='CC-BY-4.0'),
           'versions': {'py4DSTEM': md.version('py4DSTEM'), 'numpy': np.__version__, 'h5py': h5py.__version__,
                        'python': platform.python_version()},
           'volts': VOLTS, 'cell': [A_GAAS] * 3 + [90, 90, 90], 'z': [s[0] for s in SITES],
           'frac': [list(s[1]) for s in SITES], 'probe_radius_px': r, 'affine': [AFFINE_ROW, AFFINE_COLUMN],
           'pixel_start': RECIP_CAL, 'pixel': pixel,
           'matrix_layout': 'row-major 9 numbers; columns are lab x, y, z in the crystal Cartesian frame (py4DSTEM)',
           'spots_layout': 'qx, qy, intensity per spot, flat (1/Angstrom about the fitted origin)', 'cases': cases}
    OUT.write_text(json.dumps(out, allow_nan=False, separators=(',', ':')) + '\n')
    print(OUT, len(cases), 'cases, pixel', pixel)


if __name__ == '__main__':
    main()
