"""orix StereographicProjection pins (DiffractoKit S12), the reference for `StereographicProjection.project` and `direction`.

  ../References/python/orix/bin/python -I -B tools/pins/orix/stereographic-pins.py

orix 0.15.0, the default projection (pole = -1: the projection point is the south pole [00-1], the upper hemisphere is kept).
Forward: 1000 seeded directions of random length over the whole sphere, plus edge cases (the poles, the equator, z just
below and above orix's region cut of -1e-9); for each, whether orix keeps it (`v <= projection.region`, a z of the vector AS
GIVEN above -1e-9) and the (X, Y) of `vector2xy`. Inverse: 1000 seeded (X, Y) in [-2, 2]^2 and `InverseStereographicProjection
.xy2vector`. Writes Tests/DiffractoGeometryTests/Fixtures/orix-stereographic-pins.json. Deterministic.
"""
import importlib.metadata as md
import json
from pathlib import Path

import numpy as np
from orix.projections import InverseStereographicProjection, StereographicProjection
from orix.vector import Vector3d

KIT = Path(__file__).resolve().parents[3]
OUT = KIT / 'Tests/DiffractoGeometryTests/Fixtures/orix-stereographic-pins.json'
SEED = 20261123
assert md.version('orix') == '0.15.0' and md.version('numpy') == '2.5.4'


def main():
    rng = np.random.default_rng(SEED)
    v = rng.normal(size=(1000, 3))
    v /= np.linalg.norm(v, axis=1)[:, None]
    v *= rng.uniform(0.2, 5.0, size=(1000, 1))
    edges = np.array([[0, 0, 1], [0, 0, -1], [1, 0, 0], [0, 1, 0], [0.6, -0.8, 0], [1, 0, -5e-10], [1, 0, -2e-9],
                      [0.3, 0.4, 1e-12], [0, 0, 3.0], [2.0, -1.0, 0.5]], dtype=float)
    v = np.vstack([edges, v])
    proj = StereographicProjection()
    forward = []
    for row in v:
        vec = Vector3d(row)
        kept = bool(np.asarray(vec <= proj.region).ravel()[0])
        x, y = proj.vector2xy(vec)
        assert (len(np.ravel(x)) == 1) == kept
        forward.append({'v': [float(c) for c in row], 'kept': kept,
                        'xy': [float(np.ravel(x)[0]), float(np.ravel(y)[0])] if kept else None})
    xy = rng.uniform(-2, 2, size=(1000, 2))
    xy = np.vstack([[[0.0, 0.0], [1.0, 0.0], [0.0, -1.0]], xy])
    back = InverseStereographicProjection().xy2vector(xy[:, 0], xy[:, 1]).data
    inverse = [{'xy': [float(a), float(b)], 'v': [float(c) for c in r]} for (a, b), r in zip(xy, back)]
    kept = sum(f['kept'] for f in forward)
    out = {'versions': {'orix': md.version('orix'), 'numpy': np.__version__}, 'seed': SEED,
           'convention': 'pole -1: kept when z (as given) > -1e-9; X = x/(z+1), Y = y/(z+1) of the unit vector; '
                         'inverse (2X, 2Y, 1 - X^2 - Y^2) / (1 + X^2 + Y^2)',
           'forward': forward, 'inverse': inverse}
    OUT.write_text(json.dumps(out, separators=(',', ':')) + '\n')
    print(f'wrote {OUT}: {len(forward)} forward ({kept} kept), {len(inverse)} inverse')


if __name__ == '__main__':
    main()
