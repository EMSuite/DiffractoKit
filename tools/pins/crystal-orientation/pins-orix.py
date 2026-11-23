"""CrystalOrientation pins (DiffractoKit S1): the conversion to and from orix quaternions.

  ../References/python/orix/bin/python -I -B tools/pins/crystal-orientation/pins-orix.py

Writes Tests/DiffractoCoreTests/Fixtures/crystal-orientation-orix-pins.json. The oracle is orix 0.15.0:
  * `Orientation.from_matrix(M).data` for 40 seeded random rotations and six edges (identity, half turns about x, y, z and
    [111], a turn of 1e-9 rad), M the kit's matrix (columns lab x, y, z in crystal coordinates; orix's `to_matrix` is the
    same matrix, see Misorientation.swift);
  * `Orientation(q).to_matrix()` for 20 seeded random unit quaternions (a, b, c, d), some with a < 0.
Tolerance: np.allclose defaults (rtol 1e-5, atol 1e-8), the tolerance of orix's own quaternion tests.
"""
import importlib.metadata as md
import json
from pathlib import Path

import numpy as np
import scipy
from orix.quaternion import Orientation
from scipy.spatial.transform import Rotation as SR

KIT = Path(__file__).resolve().parents[3]
OUT = KIT / 'Tests/DiffractoCoreTests/Fixtures/crystal-orientation-orix-pins.json'
SEED = 20261124
assert md.version('orix') == '0.15.0'


def main():
    rng = np.random.default_rng(SEED)
    edges = [np.eye(3), np.diag([1.0, -1, -1]), np.diag([-1.0, 1, -1]), np.diag([-1.0, -1, 1]),
             SR.from_rotvec(np.pi * np.ones(3) / np.sqrt(3)).as_matrix(), SR.from_rotvec([0, 0, 1e-9]).as_matrix()]
    mats = edges + list(SR.random(40, random_state=rng.integers(2**31)).as_matrix())
    from_matrix = [{'matrix': M.reshape(-1).tolist(), 'q': Orientation.from_matrix(M).data.reshape(-1).tolist()} for M in mats]
    qs = rng.normal(size=(20, 4))
    qs /= np.linalg.norm(qs, axis=1, keepdims=True)
    to_matrix = [{'q': q.tolist(), 'matrix': Orientation(q).to_matrix().reshape(-1).tolist()} for q in qs]
    out = {'generator': 'tools/pins/crystal-orientation/pins-orix.py', 'seed': SEED,
           'versions': {'orix': md.version('orix'), 'numpy': np.__version__, 'scipy': scipy.__version__},
           'matrix_layout': 'row-major 9 numbers', 'quaternion_layout': '(a, b, c, d), a the real part',
           'from_matrix': from_matrix, 'to_matrix': to_matrix}
    OUT.write_text(json.dumps(out, allow_nan=False, indent=0) + '\n')
    print(OUT, OUT.stat().st_size)


if __name__ == '__main__':
    main()
