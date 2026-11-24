"""The simulated two-beam CBED discs for TwoBeamThickness (DiffractoGeometry): Si 220 at 200 kV, static lattice, potential
slices from abTEM 1.0.10 (Lobato, infinite projection, four slices of 1.358 Å per cubic cell), float64.

Run (an abTEM 1.0.10 venv: pip install abtem==1.0.10 numpy==2.5.4 scipy==1.18.1; about 3 min per tilt):
    <venv>/bin/python -I -B tools/pins/two-beam-thickness/disc.py [out.json]
Writes Tests/DiffractoGeometryTests/Fixtures/two-beam-thickness-disc.json. Not in tools/pins/regenerate.sh: the reference
job has no abTEM venv, and the fixture is a simulation the test measures, not a function value it matches bit for bit.

The disc along g is one plane wave per point: each sample is a beam tilted by dth from the exact Bragg position of 220
(transverse wave vector -g/2 + dth/lambda along g), plus a fixed tilt tau g perpendicular to g in the (001) plane, which
takes the beam off the [001] zone axis so that the other beams stay weak. Every wave is propagated with the multislice
propagator exp(-i pi lambda dz |k + kt|^2) through the same transmission functions, and the intensity of beam +g is
recorded every 4 cells (21.72 Å). This is the cloud experiment of 23 Nov 2026 (issue #170,
"Thickness from a two-beam CBED disc"), the same code with its outputs reduced to the profiles the tests read:

  tau 8.3  (108 mrad off [001], a clean two-beam orientation): 60, 100, 150, 200, 250 and 300 nm, each within 3 % expected
  tau 3.2  (42 mrad, a strong line through the disc): 250 nm, where two thicknesses fit equally

The deviation parameter of sample i is s = g * dth (Å^-1), exact for the profile across the disc to first order in dth.
"""
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
import abtem  # noqa: E402
import dask  # noqa: E402
from abtem.atoms import orthogonalize_cell, wrap_with_tolerance  # noqa: E402
from abtem.core.energy import energy2sigma, energy2wavelength  # noqa: E402
from ase.build import bulk  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "Tests/DiffractoGeometryTests/Fixtures/two-beam-thickness-disc.json"
dask.config.set(scheduler="synchronous")
abtem.config.set({"local_diagnostics.progress_bar": False, "precision": "float64"})
E0 = 200e3
CASES = [(8.3, "clean", [600, 1000, 1500, 2000, 2500, 3000]), (3.2, "line", [2500])]


def cell():
    a = bulk("Si", cubic=True)
    a.rotate(45, "z", rotate_cell=True)
    u = wrap_with_tolerance(orthogonalize_cell(a), tol=1e-6)
    u.positions[np.abs(u.positions) < 1e-9] = 0
    box = np.diag(u.cell)
    u.positions[:, 2] += box[2] / 8
    pot = abtem.Potential(u, sampling=0.05, slice_thickness=box[2] / 4, projection="infinite")
    v = np.asarray(pot.build().compute().array)
    n = v.shape[1] // 2
    tol = 1e-9 * v.max()
    assert np.allclose(v[:, :n, :n], v[:, n:, n:], atol=tol) and np.allclose(v[:, :n, :n], v[:, n:, :n], atol=tol)
    return v[:, :n, :n], box[0] / 2, box[2] / 4, n  # the 220 period: half the rotated cell


def run(tau, v, length, dz, n, ncell=560, step=4):
    lam = energy2wavelength(E0)
    sig = energy2sigma(E0)
    g = 2 / length
    theta_b = lam * g / 2
    trans = np.exp(1j * sig * v)
    k = np.fft.fftfreq(n, length / n)
    kx, ky = np.meshgrid(k, k, indexing="ij")
    mask = (kx ** 2 + ky ** 2) < (2 / 3 * k.max()) ** 2
    dth = np.linspace(-6.5e-3, 6.5e-3, 261)
    ktx = (-theta_b + dth) / lam
    kty = tau * g
    prop = np.exp(-1j * np.pi * lam * dz * ((kx[None] + ktx[:, None, None]) ** 2 + (ky[None] + kty) ** 2)) * mask[None]
    psi = np.ones((len(dth), n, n), complex)
    out = np.zeros((len(dth), ncell // step))
    for c in range(ncell):
        for s in range(4):
            psi = np.fft.ifft2(np.fft.fft2(psi * trans[s]) * prop)
        if (c + 1) % step == 0:
            out[:, (c + 1) // step - 1] = np.abs(np.fft.fft2(psi)[:, 2, 0] / n ** 2) ** 2
    thickness = (np.arange(ncell // step) + 1) * step * 4 * dz
    return dth, thickness, out, g, lam


def main():
    v, length, dz, n = cell()
    cases = []
    for tau, name, wanted in CASES:
        t0 = time.perf_counter()
        dth, thickness, ig, g, lam = run(tau, v, length, dz, n)
        print(f"tau {tau}: {time.perf_counter() - t0:.0f} s", flush=True)
        for w in wanted:
            j = int(np.argmin(np.abs(thickness - w)))
            cases.append({
                "name": f"{name}-{round(thickness[j] / 10)}nm",
                "orientation": name,
                "tilt_off_001_mrad": tau * g * lam * 1e3,
                "thickness": float(thickness[j]),
                "deviation": [float(x) for x in g * dth],
                "intensity": [float(x) for x in ig[:, j]],
            })
    doc = {
        "source": "tools/pins/two-beam-thickness/disc.py: abTEM 1.0.10 potential (Lobato, infinite projection), "
                  "plane-wave multislice per point of the 220 disc, float64",
        "reflection": [2, 2, 0], "phase": "Si (Fd-3m, a = 5.43 Å, ASE bulk)", "volts": E0,
        "wavelength": float(lam), "g": float(g), "slice_thickness": float(dz),
        "pendelloesung_period_at_bragg": 956.0,
        "cases": cases,
    }
    OUT.write_text(json.dumps(doc, indent=1) + "\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
