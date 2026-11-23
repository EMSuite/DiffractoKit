"""Camera-length pins (DiffractoKit S2): R = L tan 2θ, sin θ = λ/2d, in 60-digit arithmetic (mpmath).

  ../References/python/py4dstem/bin/python -I -B tools/pins/spot-pattern/radius-pins.py

Writes Tests/DiffractoGeometryTests/Fixtures/spot-radius-pins.json: 540 cases, 3 voltages (80, 200, 300 kV) x 5 camera
lengths (100, 250, 500, 1000, 2000 mm) x 36 spacings (0.4 to 12 Å, log spaced). λ = hc / sqrt(eV (2 m0c² + eV)) with
the CODATA 2018 constants of `ElectronWavelength` (h, c, e exact; m0 = 9.1093837015e-31 kg). No Python package states a
tolerance for the form: assert_allclose's float64 default (rtol 1e-7) applies.
"""
import importlib.metadata as md
import json
from pathlib import Path

import mpmath as mp

mp.mp.dps = 60
KIT = Path(__file__).resolve().parents[3]
OUT = KIT / 'Tests/DiffractoGeometryTests/Fixtures/spot-radius-pins.json'
h, c, e, m0 = mp.mpf('6.62607015e-34'), mp.mpf(299792458), mp.mpf('1.602176634e-19'), mp.mpf('9.1093837015e-31')


def wavelength(volts):
    v = mp.mpf(volts)
    hc = h * c / e * mp.mpf(10) ** 10
    mc2 = m0 * c * c / e
    return hc / mp.sqrt(v * (2 * mc2 + v))


def main():
    cases = []
    ds = [mp.mpf('0.4') * (mp.mpf(30) ** (mp.mpf(i) / 35)) for i in range(36)]
    for volts in (80_000, 200_000, 300_000):
        lam = wavelength(volts)
        for L in (100, 250, 500, 1000, 2000):
            for d in ds:
                d = mp.mpf(float(d))       # the double the test passes
                theta = mp.asin(lam / (2 * d))
                cases.append([volts, L, float(d), float(lam), float(mp.mpf(L) * mp.tan(2 * theta))])
    OUT.write_text(json.dumps({'generator': 'tools/pins/spot-pattern/radius-pins.py', 'versions': {'mpmath': md.version('mpmath')},
                               'digits': 60, 'layout': 'volts, camera length (mm), d (Å), λ (Å), R (mm)', 'cases': cases},
                              separators=(',', ':')) + '\n')
    print(OUT, len(cases))


if __name__ == '__main__':
    main()
