"""SpotPattern cross-check (DiffractoKit S2): diffsims 0.7.0 `SimulationGenerator.calculate_diffraction2d`.

  ../References/python/diffsims/bin/python -I -B tools/pins/spot-pattern/pins-diffsims.py

Writes Tests/DiffractoGeometryTests/Fixtures/spot-pattern-diffsims-pins.json. The second check of plan row S2 (py4DSTEM is the
reference). Per phase of tools/pins/phase-vectors (sites already expanded, given to diffsims as P1) at 200 kV, eight seeded
random orix rotations, `calculate_diffraction2d(reciprocal_radius=1.5, with_direct_beam=False,
max_excitation_error=0.03)` (lobato scattering, lorentzian shape factor: diffsims' defaults):
  * per spot h, k, l, diffsims' x, y, z (Å⁻¹, lab), its intensity and |F|² = intensity / its lorentzian of its own s;
  * the orientation in the kit's convention: R (det +1) with xyz = R g, g the kit's Cartesian g (py4DSTEM's lat_inv), fitted
    by Kabsch on all spots (rms recorded); the kit's matrix is R^T (both put the beam along -z).
diffsims states np.allclose (rtol 1e-5, atol 1e-8) in its tests; the plan row sets position 1e-7, |F|^2 1e-6 relative and
s 1e-5 Å⁻¹. Its two documented deviations are left out of the comparison: (1) diffsims keeps reflections the kit's |F| cut
(1e-4) drops, (2) its excitation error is the exact distance to the sphere along z, s = 1/λ − sqrt(1/λ² − x² − y²) − z,
and its envelope a lorentzian, so intensities are compared before the envelope.
"""
import importlib.metadata as md
import json
import math
from pathlib import Path

import numpy as np
from diffpy.structure import Atom, Lattice, Structure
from diffsims.generators.simulation_generator import SimulationGenerator
from orix.crystal_map import Phase
from orix.quaternion import Rotation

KIT = Path(__file__).resolve().parents[3]
SRC = KIT / 'Tests/DiffractoOrientationTests/Fixtures/py4dstem-phase-vector-pins.json'
OUT = KIT / 'Tests/DiffractoGeometryTests/Fixtures/spot-pattern-diffsims-pins.json'
SEED = 20261125
RADIUS, SMAX, KV = 1.5, 0.03, 200.0
SYMBOL = {1: 'H', 6: 'C', 12: 'Mg', 13: 'Al', 14: 'Si', 26: 'Fe', 29: 'Cu', 79: 'Au', 3: 'Li', 8: 'O', 22: 'Ti', 28: 'Ni'}
assert md.version('diffsims') == '0.7.0'


def lat_inv(cell):
    """py4DSTEM's lat_real (a along x, b in the xy plane) and its inverse transpose (rows are b1, b2, b3)."""
    a, b, c = cell['a'], cell['b'], cell['c']
    al, be, ga = (math.radians(cell[k]) for k in ('alpha', 'beta', 'gamma'))
    cx = c * math.cos(be)
    cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    real = np.array([[a, 0, 0], [b * math.cos(ga), b * math.sin(ga), 0], [cx, cy, math.sqrt(c * c - cx * cx - cy * cy)]])
    return np.linalg.inv(real).T


def kabsch(P, Q):
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, float(np.sqrt(np.mean(np.sum((Q - P @ R.T) ** 2, axis=1))))


def main():
    src = json.load(open(SRC))
    rng = np.random.default_rng(SEED)
    gen = SimulationGenerator(accelerating_voltage=KV)
    lam = gen.wavelength
    out = {'generator': 'tools/pins/spot-pattern/pins-diffsims.py', 'seed': SEED,
           'versions': {p: md.version(p) for p in ('diffsims', 'orix', 'diffpy.structure', 'numpy')},
           'reciprocal_radius': RADIUS, 'max_excitation_error': SMAX, 'volts': KV * 1e3, 'wavelength': lam,
           'row_layout': 'h, k, l, x, y, z, intensity, |F|^2 (diffsims)', 'phases': []}
    for ph in src['phases']:
        c = ph['cell']
        latt = Lattice(c['a'], c['b'], c['c'], c['alpha'], c['beta'], c['gamma'])
        atoms = [Atom(SYMBOL[s['z']], xyz=list(s['frac']), lattice=latt, occupancy=s['occupancy']) for s in ph['sites']]
        phase = Phase(name=ph['name'], space_group=1, structure=Structure(atoms=atoms, lattice=latt))
        B = lat_inv(c)
        rec = {'name': ph['name'], 'patterns': []}
        for _ in range(8):
            q = rng.normal(size=4)
            rot = Rotation(q / np.linalg.norm(q))
            sim = gen.calculate_diffraction2d(phase, rotation=rot, reciprocal_radius=RADIUS, with_direct_beam=False,
                                              max_excitation_error=SMAX)
            v = sim.coordinates
            hkl = np.rint(v.hkl).astype(int)
            xyz = np.asarray(v.data, float)
            inten = np.asarray(v.intensity, float)
            if len(hkl) < 3:
                continue
            R, rms = kabsch(hkl @ B, xyz)
            r = 1 / lam
            s = (r - np.sqrt(r * r - xyz[:, 0] ** 2 - xyz[:, 1] ** 2)) - xyz[:, 2]
            sigma = math.pi / SMAX
            shape = sigma / (math.pi * (sigma ** 2 * s ** 2 + 1)) * SMAX
            rows = [[*map(int, h), *map(float, p), float(i), float(i / f)] for h, p, i, f in zip(hkl, xyz, inten, shape)]
            rec['patterns'].append({'matrix': R.T.reshape(-1).tolist(), 'kabsch_rms': rms, 'spots': rows})
        out['phases'].append(rec)
    OUT.write_text(json.dumps(out, allow_nan=False, separators=(',', ':')) + '\n')
    print(OUT, OUT.stat().st_size, sum(len(p['spots']) for ph in out['phases'] for p in ph['patterns']))


if __name__ == '__main__':
    main()
