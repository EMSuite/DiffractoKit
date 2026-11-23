"""DiffractoKit: the pin generator for powder X-ray patterns (PowderPattern.xray) and pymatgen's X-ray scattering table.

Run, from any directory, with the pinned oracle:

    ../References/python/pymatgen/bin/python -I -B tools/pins/powder-xrd/pins.py

It writes Tests/DiffractoScatteringTests/Fixtures/pymatgen-powder-xrd-pins.json. Deterministic: no clock, no random numbers;
two runs give the same bytes.

The oracle is pymatgen's own `XRDCalculator.get_pattern` (pymatgen 2026.9.24, the code in pymatgen-core 2026.10.16), called
on a `Structure` built from the cell (`Lattice.from_parameters(..., vesta=True)`, the kit's convention) and the sites of
each structure exactly as the committed structure-factor pins hold them (Tests/DiffractoScatteringTests/Fixtures/
numpy-structure-factor-pins.json: the presets' exact sites, the CIFs' sites as written). Every structure of that file but
the synthetic `Weak-cut` (a unit form factor construct of the |F| cut) is taken: the 19 committed structures. Two planted
structures the 2.4 design names are added: Fe (bcc, the preset's a = 2.8665 A) and Mg2Si (antifluorite, Fm-3m, a = 6.351 A).

What is pinned, per case:
  - the peaks of `get_pattern(scaled=True)`: 2 theta, the scaled intensity (strongest 100), d, and pymatgen's hkl labels
    with their multiplicities (`get_unique_families`: permutations of |h|, |k|, |l|; four indices h k i l on a lattice
    pymatgen calls hexagonal);
  - the unscaled intensities of `get_pattern(scaled=False)` (2 sum |F|^2 LP over the merged reflections);
  - the case's settings: wavelength, 2 theta range, Debye-Waller factors.
Cases: every structure with Cu K-alpha (1.54184 A) at 0 to 90 degrees (pymatgen's defaults); Si with Mo K-alpha; Al and
calcite at 20 to 60 degrees (the lower bound filters); Mg at 0 to 180 degrees (the whole limiting sphere); Au and WS2 with
Debye-Waller factors.

The table: pymatgen's `atomic_scattering_params.json` (99 elements, H to Cf, four (a, b) pairs each), written out with its
sha256 so the kit's copy is held to it bit for bit.

Exit status: 0 all self-checks hold; 1 a self-check failed (assert, loud).
"""
import hashlib
import json
import math
import os
import sys
import warnings
from importlib.metadata import version as dist_version

import numpy as np

import pymatgen.analysis.diffraction.xrd as xrd_module
from pymatgen.analysis.diffraction.xrd import ATOMIC_SCATTERING_PARAMS, WAVELENGTHS, XRDCalculator
from pymatgen.core import Element, Lattice, Structure

HERE = os.path.dirname(os.path.abspath(__file__))
KIT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
FIXTURES = os.path.join(KIT, "Tests", "DiffractoScatteringTests", "Fixtures")
SOURCE = os.path.join(FIXTURES, "numpy-structure-factor-pins.json")
OUT = os.path.join(FIXTURES, "pymatgen-powder-xrd-pins.json")

LEFT_OUT = {"Weak-cut"}
PLANTED = [
    ("Fe", (2.8665, 2.8665, 2.8665, 90.0, 90.0, 90.0), [(26, (0.0, 0.0, 0.0), 1.0), (26, (0.5, 0.5, 0.5), 1.0)]),
    ("Mg2Si", (6.351, 6.351, 6.351, 90.0, 90.0, 90.0),
     [(14, p, 1.0) for p in ((0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))]
     + [(12, (x, y, z), 1.0) for x in (0.25, 0.75) for y in (0.25, 0.75) for z in (0.25, 0.75)]),
]
EXTRA_CASES = [  # (structure, radiation key or None, wavelength, two-theta range, Debye-Waller factors by symbol)
    ("Si", "MoKa", None, (0.0, 90.0), {}),
    ("Al", "CuKa", None, (20.0, 60.0), {}),
    ("calcite_R-3c_hex", "CuKa", None, (20.0, 60.0), {}),
    ("Mg", "CuKa", None, (0.0, 180.0), {}),
    ("Au", "CuKa", None, (0.0, 90.0), {"Au": 0.62}),
    ("WS2", "CuKa", None, (0.0, 90.0), {"W": 0.25, "S": 0.5}),
]

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append(name)
    assert ok, f"SELF-CHECK FAILED: {name} {detail}"


def structure_of(cell, sites):
    lattice = Lattice.from_parameters(*cell, vesta=True)
    species = [{Element.from_Z(z).symbol: occ} for (z, _, occ) in sites]
    coords = [list(p) for (_, p, _) in sites]
    return Structure(lattice, species, coords)


def families_json(fams):
    return [{"indices": [int(v) for v in e["hkl"]], "multiplicity": int(e["multiplicity"])} for e in fams]


def case_block(name, cell, sites, radiation, wavelength, tt_range, dw):
    st = structure_of(cell, sites)
    lam = WAVELENGTHS[radiation] if radiation else wavelength
    calc = XRDCalculator(wavelength=radiation if radiation else wavelength, debye_waller_factors=dw or None)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        scaled = calc.get_pattern(st, scaled=True, two_theta_range=tt_range)
        raw = calc.get_pattern(st, scaled=False, two_theta_range=tt_range)
    check(f"{name}: scaled and unscaled patterns have the same peaks", list(scaled.x) == list(raw.x))
    check(f"{name}: the strongest scaled peak is 100", abs(max(scaled.y) - 100) <= 1e-12)
    peaks = []
    for x, y, yr, hk, d in zip(scaled.x, scaled.y, raw.y, scaled.hkls, scaled.d_hkls):
        peaks.append({"two_theta": float(x), "intensity": float(y), "absolute_intensity": float(yr), "d": float(d),
                      "families": families_json(hk)})
    tt = [p["two_theta"] for p in peaks]
    check(f"{name}: peaks in ascending 2 theta", tt == sorted(tt))
    print(f"  {name}: {radiation or wavelength} {tt_range} dw={dw}: {len(peaks)} peaks"
          + (f", strongest at {peaks[int(np.argmax([p['intensity'] for p in peaks]))]['two_theta']:.4f}" if peaks else ""))
    return {"structure": name, "radiation": radiation, "wavelength": float(lam),
            "two_theta_range": [float(tt_range[0]), float(tt_range[1])],
            "debye_waller_factors": {str(Element(s).Z): float(b) for s, b in sorted(dw.items())},
            "is_hexagonal_lattice": bool(st.lattice.is_hexagonal()),
            "peak_count": len(peaks), "peaks": peaks}


def main():
    src = json.load(open(SOURCE, encoding="utf-8"))
    structures = []
    for s in src["structures"]:
        if s["name"] in LEFT_OUT:
            continue
        c = s["cell"]
        cell = (c["a"], c["b"], c["c"], c["alpha"], c["beta"], c["gamma"])
        sites = [(x["z"], tuple(x["frac"]), x["occupancy"]) for x in s["sites"]]
        structures.append((s["name"], cell, sites, "structure-factor pins (" + s["build"]["how"] + ")"))
    check("19 committed structures", len(structures) == 19, str(len(structures)))
    for name, cell, sites in PLANTED:
        structures.append((name, cell, sites, "planted"))
    by_name = {n: (c, s) for n, c, s, _ in structures}

    print("--- pymatgen XRDCalculator.get_pattern ---")
    cases = [case_block(n, c, s, "CuKa", None, (0.0, 90.0), {}) for n, c, s, _ in structures]
    for (n, rad, lam, rng, dw) in EXTRA_CASES:
        c, s = by_name[n]
        cases.append(case_block(n, c, s, rad, lam, rng, dw))

    al = next(x for x in cases if x["structure"] == "Al" and x["two_theta_range"] == [0.0, 90.0] and not x["debye_waller_factors"])
    check("Al with Cu K-alpha: 5 peaks, 111 at 38.5 degrees and strongest",
          al["peak_count"] == 5 and abs(al["peaks"][0]["two_theta"] - 38.51) < 0.01 and al["peaks"][0]["intensity"] == 100,
          str([(p["two_theta"], p["intensity"]) for p in al["peaks"]]))

    table_path = os.path.join(os.path.dirname(xrd_module.__file__), "atomic_scattering_params.json")
    table_bytes = open(table_path, "rb").read()
    table = [{"symbol": sym, "z": Element(sym).Z, "coefficients": [[float(a), float(b)] for a, b in ATOMIC_SCATTERING_PARAMS[sym]]}
             for sym in ATOMIC_SCATTERING_PARAMS]
    table.sort(key=lambda e: e["z"])
    zs = [e["z"] for e in table]
    check("the table: 99 symbols, Z 1 to 98 without a gap; D is a second Z = 1 row equal to H",
          len(table) == 99 and sorted(set(zs)) == list(range(1, 99))
          and [e["coefficients"] for e in table if e["z"] == 1][0] == [e["coefficients"] for e in table if e["z"] == 1][1],
          str(len(table)))
    check("four (a, b) pairs per element", all(len(e["coefficients"]) == 4 for e in table))

    out = {
        "versions": {"pymatgen": dist_version("pymatgen"), "pymatgen-core": dist_version("pymatgen-core"), "numpy": np.__version__},
        "python": sys.version.split()[0],
        "definitions": {
            "f_x": "f(s) = Z - 41.78214 s^2 sum_i a_i exp(-b_i s^2), s = sin(theta)/lambda = |g|/2, |g| = 1/d (no 2 pi)",
            "lp": "(1 + cos^2 2theta) / (sin^2 theta cos theta)",
            "intensity": "2 |F|^2 LP over the Friedel half-space (F = sum_e f_e exp(-B_e s^2) sum_j occ_j exp(2 pi i h.r_j)), "
                         "merged where round(2theta / 1e-5) is equal, peaks at or below 1e-3 percent of the strongest dropped, "
                         "scaled so the strongest is 100",
            "families": "pymatgen get_unique_families: hkl and -hkl of every merged reflection, grouped by sorted |indices|, "
                        "each named by its largest member; four indices (h, k, -h-k, l) when Lattice.is_hexagonal()",
        },
        "table_file": "pymatgen/analysis/diffraction/atomic_scattering_params.json",
        "table_sha256_of_file": hashlib.sha256(table_bytes).hexdigest(),
        "wavelengths": {k: float(v) for k, v in WAVELENGTHS.items()},
        "structures": [{"name": n, "source": src_, "cell": list(c), "sites": [{"z": z, "frac": list(p), "occupancy": o} for z, p, o in s]}
                       for n, c, s, src_ in structures],
        "table": table,
        "cases": cases,
    }
    dump = lambda o: json.dumps(o, separators=(",", ":"), allow_nan=False)
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("{")
        keys = list(out)
        for i, k in enumerate(keys):
            if k in ("structures", "table", "cases"):
                f.write(f'"{k}":[\n' + ",\n".join(dump(x) for x in out[k]) + "]")
            else:
                f.write(f'"{k}":' + dump(out[k]))
            f.write(",\n" if i < len(keys) - 1 else "\n")
        f.write("}\n")
    json.load(open(OUT, encoding="utf-8"))
    print(f"\nwrote {OUT} ({os.path.getsize(OUT) / 1e6:.2f} MB, {len(cases)} cases); {len(CHECKS)} self-checks passed")


if __name__ == "__main__":
    main()
