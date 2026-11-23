"""Planted truth for the double-tilt holder (DiffractoKit S12): `DoubleTiltHolder.tilts`, `orientation(_:at:movedTo:)` and
`calibrateAxis`.

  ../References/python/pymatgen/bin/python -I -B tools/pins/tilt-planner/planted.py

No Python package computes holder tilts (ALPHABETA, tompas and pycotem expose no callable function at a pinned version), so
the pins are planted from the model in NumPy, written from the definitions below and not from any Swift:

  holder: a sample vector s goes to the pattern frame as  p = Rz(phi) Rx(alpha) Ry(beta) s,
          phi the angle of the alpha axis in the pattern (from lab x toward lab y), beta about the axis alpha carries;
  pattern orientation (the kit's CrystalOrientation, py4DSTEM's matrix): M^T = Rz(phi) Rx(alpha) Ry(beta) S, S fixed;
  the tilts that bring a crystal direction t (sample frame at zero tilt: S t) onto the beam (+z):
          alpha = asin(t_y), beta = atan2(-t_x, t_z).

Cases (seeded):
  tilts:        500 planted (S, phi, current tilts, target direction): the closed-form tilts and the moved orientation.
  calibration:  300 trials in each band of tilt change, 10 to 20 and 20 to 30 degrees: two orientations of one crystal at
                two random tilts within +-30 degrees, each with an independent 0.2 degree error (a rotation about a random
                axis), the planted phi; half of the trials store the second orientation turned by a half turn about the beam
                (the ambiguity of one spot pattern), flagged.
The NumPy calibration below (phi by a 0.5 degree scan and a golden refinement of the misorientation between the measured
rotation and Rz(phi) Q Rz(phi)^T, both half turns of the second pattern tried) is run on the same trials as a self-check of
the bounds the Swift test holds: the 95th percentile of the axis error at most 1.5 degrees from 10 degrees of tilt change,
and the half turn named right in every trial.
Writes Tests/DiffractoGeometryTests/Fixtures/tilt-planner-planted.json.
"""
import json
import math
from pathlib import Path

import numpy as np

KIT = Path(__file__).resolve().parents[3]
OUT = KIT / 'Tests/DiffractoGeometryTests/Fixtures/tilt-planner-planted.json'
SEED = 20261124
NOISE_DEG = 0.2


def rx(a):
    c, s = math.cos(a), math.sin(a)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def ry(a):
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def rz(a):
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def axis_angle(axis, angle):
    k = np.asarray(axis, float) / np.linalg.norm(axis)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + math.sin(angle) * K + (1 - math.cos(angle)) * K @ K


def random_rotation(rng):
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])


def angle_of(R):
    return math.acos(max(-1.0, min(1.0, (np.trace(R) - 1) / 2)))


def holder(phi, a, b):
    return rz(phi) @ rx(a) @ ry(b)


def calibrate(Ma, ta, Mb, tb):
    """phi (rad), misfit (rad), half turn (bool), the other half turn's misfit"""
    Q = rx(tb[0]) @ ry(tb[1]) @ ry(ta[1]).T @ rx(ta[0]).T
    best = []
    for k in (0, 1):
        R = rz(k * math.pi).T @ Mb.T @ Ma
        f = lambda p: angle_of(R @ (rz(p) @ Q @ rz(p).T).T)
        grid = [math.radians(0.5 * i) for i in range(720)]
        p0 = min(grid, key=f)
        lo, hi = p0 - math.radians(0.5), p0 + math.radians(0.5)
        g = (math.sqrt(5) - 1) / 2
        for _ in range(60):
            c, d = hi - g * (hi - lo), lo + g * (hi - lo)
            if f(c) < f(d):
                hi = d
            else:
                lo = c
        p = (lo + hi) / 2
        best.append((f(p), p % (2 * math.pi), k == 1))
    best.sort()
    return best[0][1], best[0][0], best[0][2], best[1][0]


def M_of(phi, t, S):
    return (holder(phi, *t) @ S).T


def main():
    rng = np.random.default_rng(SEED)
    lim = math.radians(30)
    tilts = []
    worst = 0.0
    for _ in range(500):
        S = random_rotation(rng)
        phi = rng.uniform(0, 2 * math.pi)
        cur = rng.uniform(-lim, lim, size=2)
        t = rng.normal(size=3)
        t /= np.linalg.norm(t)
        M0 = M_of(phi, cur, S)
        ts = S @ t
        alpha, beta = math.asin(ts[1]), math.atan2(-ts[0], ts[2])
        M1 = M_of(phi, (alpha, beta), S)
        on_beam = M1.T @ t
        worst = max(worst, math.degrees(math.atan2(math.hypot(on_beam[0], on_beam[1]), on_beam[2])))
        tilts.append({'M0': M0.tolist(), 'phi_deg': math.degrees(phi), 'current_deg': np.degrees(cur).tolist(),
                      'direction': t.tolist(), 'target_deg': [math.degrees(alpha), math.degrees(beta)], 'M1': M1.tolist()})
    assert worst < 1e-9, worst
    print(f'tilts: 500 planted, target on the beam within {worst:.1e} degrees')

    calibration = []
    for band in ((10.0, 20.0), (20.0, 30.0)):
        errs, wrong = [], 0
        n = 0
        while n < 300:
            S = random_rotation(rng)
            phi = rng.uniform(0, 2 * math.pi)
            ta, tb = rng.uniform(-lim, lim, size=2), rng.uniform(-lim, lim, size=2)
            change = math.degrees(angle_of(rx(tb[0]) @ ry(tb[1]) @ ry(ta[1]).T @ rx(ta[0]).T))
            if not band[0] <= change < band[1]:
                continue
            noise = lambda: axis_angle(rng.normal(size=3), math.radians(NOISE_DEG))
            Ma = (noise() @ holder(phi, *ta) @ S).T
            Mb = (noise() @ holder(phi, *tb) @ S).T
            flipped = bool(n % 2)
            if flipped:
                Mb = Mb @ rz(math.pi)
            p, misfit, half, other = calibrate(Ma, ta, Mb, tb)
            err = abs((math.degrees(p - phi) + 180) % 360 - 180)
            errs.append(err)
            wrong += half != flipped
            calibration.append({'band_deg': list(band), 'tilt_change_deg': change, 'phi_deg': math.degrees(phi),
                                'tilts_a_deg': np.degrees(ta).tolist(), 'tilts_b_deg': np.degrees(tb).tolist(),
                                'Ma': Ma.tolist(), 'Mb': Mb.tolist(), 'b_half_turned': flipped,
                                'numpy_phi_deg': math.degrees(p), 'numpy_misfit_deg': math.degrees(misfit)})
            n += 1
        p95 = float(np.percentile(errs, 95))
        print(f'calibration {band}: median {np.median(errs):.2f}, p95 {p95:.2f} degrees; half turn wrong in {wrong}')
        assert p95 <= 1.5 and wrong == 0
    out = {'seed': SEED, 'noise_deg': NOISE_DEG, 'model': __doc__.split('Cases')[0].strip(),
           'matrices': 'M as rows: the kit CrystalOrientation.matrix (columns = lab x, y, z in the crystal frame)',
           'tilts': tilts, 'calibration': calibration}
    OUT.write_text(json.dumps(out, separators=(',', ':')) + '\n')
    print(f'wrote {OUT}')


if __name__ == '__main__':
    main()
