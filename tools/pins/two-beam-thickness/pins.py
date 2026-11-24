"""Closed-form pins for TwoBeamThickness (DiffractoGeometry): the two-beam intensity of the diffracted beam, the positions
of its minima and the fringe line through them (Kelly, Jostsons, Blake and Napier, phys. stat. sol. (a) 31, 771 (1975);
Allen, Philos. Mag. A 43, 325 (1981)), in NumPy, and the two-beam extinction distance of Si 220 at 200 kV from the Lobato
table (py4DSTEM's, the kit's ElectronScatteringFactors) with the kit's constants.

Run (the pymatgen venv, see tools/pins/requirements/pymatgen.txt; only NumPy is used):
    ../References/python/pymatgen/bin/python -I -B tools/pins/two-beam-thickness/pins.py [out.json]
Writes Tests/DiffractoGeometryTests/Fixtures/two-beam-thickness-pins.json (deterministic).

    I(s) = sin^2(pi t s_eff) / (xi s_eff)^2,   s_eff^2 = s^2 + 1/xi^2
    s_n  = sqrt(n^2/t^2 - 1/xi^2),  n > t/xi
    (s_n/n)^2 = 1/t^2 - (1/xi^2)(1/n^2)       numpy.polyfit(1/n^2, (s_n/n)^2, 1)

No Python package does this measurement, so the closed form is the reference. Tolerances (the issue): the intensity and the
fitted thickness and extinction distance on noiseless profiles rtol 1e-6, the fringe line through exact minima rtol 1e-7
(assert_allclose's float64 default).
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "Tests/DiffractoGeometryTests/Fixtures/two-beam-thickness-pins.json"

# (thickness, extinction distance) in Å: 40 to 400 nm, 50 to 150 nm
CASES = [(400.0, 500.0), (400.0, 1500.0), (1000.0, 956.0), (1500.0, 1200.0), (2500.0, 700.0), (4000.0, 500.0),
         (4000.0, 1500.0)]
# Si: the Lobato row of Z = 14 (the kit's table, py4DSTEM's scattering_factors.txt), diamond cubic, a = 5.431 Å (ASE)
SI_LOBATO = [2.87189142611612350e+00, 5.08487103642989610e+00, -2.06173501195173530e+00, 4.29178185305126190e-01,
             2.17114024204478720e+00, 3.66485434192162230e-01, -6.63073633058801900e-02, 1.19710611296903400e-01,
             3.01070709670513740e-03, 1.43994536128397540e-02]
HC_EV_A = 6.626_070_15e-34 * 299_792_458.0 / 1.602_176_634e-19 * 1e10
MC2_EV = 9.109_383_7015e-31 * 299_792_458.0 ** 2 / 1.602_176_634e-19


def intensity(s, t, xi):
    se = np.sqrt(s * s + 1 / (xi * xi))
    return np.sin(np.pi * t * se) ** 2 / (xi * se) ** 2


def minimum(n, t, xi):
    return math.sqrt(n * n / (t * t) - 1 / (xi * xi))


def si220(volts):
    a = 5.431
    g2 = 8 / a ** 2
    f = sum(SI_LOBATO[2 * i] * (2 + SI_LOBATO[2 * i + 1] * g2) / (1 + SI_LOBATO[2 * i + 1] * g2) ** 2 for i in range(5))
    big_f = 8 * f / a ** 3  # |F_220| of diamond cubic: 8 f, over V
    lam = HC_EV_A / math.sqrt(volts * (2 * MC2_EV + volts))
    gamma = 1 + volts / MC2_EV
    sin_b = lam * math.sqrt(g2) / 2
    return math.pi * math.sqrt(1 - sin_b ** 2) / (lam * gamma * big_f)


def main():
    cases = []
    for t, xi in CASES:
        n0 = math.floor(t / xi) + 1
        orders = list(range(n0, n0 + 5))
        mins = [minimum(n, t, xi) for n in orders]
        # a profile reaching half a fringe past the fourth minimum, 241 samples, off-centre by 1e-5 Å^-1
        smax = 0.5 * (mins[3] + mins[4])
        s = np.linspace(-smax, smax, 241)
        slope, icpt = np.polyfit(1 / np.asarray(orders, float) ** 2, (np.asarray(mins) / np.asarray(orders, float)) ** 2, 1)
        cases.append({
            "thickness": t, "extinction_distance": xi,
            "deviation": [float(x) for x in s[::20]],
            "intensity": [float(x) for x in intensity(s[::20], t, xi)],
            "orders": orders, "minima": mins,
            "profile_half_width": smax, "profile_samples": 241, "profile_centre": 1e-5,
            "fringe_line_thickness": 1 / math.sqrt(icpt), "fringe_line_extinction_distance": 1 / math.sqrt(-slope),
        })
    doc = {
        "source": "tools/pins/two-beam-thickness/pins.py: the closed form in NumPy " + np.__version__,
        "cases": cases,
        "si220_200kV_extinction_distance": si220(200e3),
        "si220_pendelloesung_period_abtem": 956.0,
    }
    OUT.write_text(json.dumps(doc, indent=1) + "\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
