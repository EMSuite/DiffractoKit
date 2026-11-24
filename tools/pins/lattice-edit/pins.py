"""DiffractoKit: the pin generator for lattice edits within the crystal system (DiffractoKit #149, plan row K2).

Run, from any directory, with the pinned oracle:

    ../References/python/pymatgen/bin/python -I -B tools/pins/lattice-edit/pins.py

It writes Tests/DiffractoCoreTests/Fixtures/lattice-edit-pins.json. No Python function edits a cell within its crystal
system, so the edited cells are STATED here, case by case, from the rule (the free parameters of each system, the tied
ones following, the fixed ones kept): cubic a; tetragonal and hexagonal a, c; trigonal in hexagonal axes a, c, in
rhombohedral axes a, alpha; orthorhombic a, b, c; monoclinic (unique axis b) a, b, c, beta; triclinic all six; a
primitive cell of a centred cubic lattice is scaled. pymatgen 2026.9.24 (pymatgen-core 2026.10.16) is the oracle for
the d-spacings of each edited cell: `Lattice.from_parameters(..., vesta=True).d_hkl` for every hkl in -3...3 but 000.
The Swift test edits the kit's structure, checks the free parameters and the edited cell against this file, and its d
against pymatgen's at rtol 1e-7 (assert_allclose's float64 default; pymatgen ships no test for d_hkl).

File:line pins into the pymatgen-core clone (References/pymatgen-core-ae9a834/src/pymatgen/):
  core/lattice.py       Lattice.from_parameters (vesta=True branch) 389-445, d_hkl 258
"""
import datetime
import itertools
import json
import math
import os
from importlib.metadata import version as dist_version

import numpy as np
from pymatgen.core import Lattice

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "..", "Tests", "DiffractoCoreTests", "Fixtures", "lattice-edit-pins.json")

# name, where the Swift test takes the structure from, the system, the input cell, the free parameters, the edit, the
# edited cell (the rule applied by hand).
A_PRIM = 4.0495 / math.sqrt(2)
CASES = [
    ("silicon", "preset:silicon", "cubic", (5.4309, 5.4309, 5.4309, 90, 90, 90), ["a"],
     {"a": 5.5}, (5.5, 5.5, 5.5, 90, 90, 90)),
    ("magnesium", "preset:magnesium", "hexagonal", (3.2094, 3.2094, 5.2108, 90, 90, 120), ["a", "c"],
     {"a": 3.25, "c": 5.25}, (3.25, 3.25, 5.25, 90, 90, 120)),
    ("alpha quartz", "cif:alpha_quartz_P3221.cif", "trigonal", (4.9134, 4.9134, 5.4052, 90, 90, 120), ["a", "c"],
     {"a": 4.95, "c": 5.45}, (4.95, 4.95, 5.45, 90, 90, 120)),
    ("rutile", "code:rutile", "tetragonal", (4.594, 4.594, 2.959, 90, 90, 90), ["a", "c"],
     {"a": 4.6, "c": 3.0}, (4.6, 4.6, 3.0, 90, 90, 90)),
    ("orthorhombic", "code:orthorhombic", "orthorhombic", (3.0, 4.0, 5.0, 90, 90, 90), ["a", "b", "c"],
     {"a": 3.1, "b": 4.2, "c": 4.9}, (3.1, 4.2, 4.9, 90, 90, 90)),
    ("beta''", "preset:betaDoublePrime", "monoclinic", (15.16, 4.05, 6.74, 90, 105.3, 90), ["a", "b", "c", "beta"],
     {"a": 15.2, "b": 4.1, "c": 6.8, "beta": 106.0}, (15.2, 4.1, 6.8, 90, 106.0, 90)),
    ("triclinic", "code:triclinic", "triclinic", (5.0, 6.0, 7.0, 80, 85, 95),
     ["a", "b", "c", "alpha", "beta", "gamma"],
     {"a": 5.1, "b": 6.1, "c": 7.1, "alpha": 81.0, "beta": 86.0, "gamma": 96.0}, (5.1, 6.1, 7.1, 81, 86, 96)),
    ("rhombohedral axes", "code:rhombohedral", "trigonal", (6.375, 6.375, 6.375, 46.08, 46.08, 46.08), ["a", "alpha"],
     {"a": 6.4, "alpha": 47.0}, (6.4, 6.4, 6.4, 47, 47, 47)),
    ("aluminium, primitive cell", "code:aluminiumPrimitive", "cubic", (A_PRIM, A_PRIM, A_PRIM, 60, 60, 60), ["a"],
     {"a": 2.9}, (2.9, 2.9, 2.9, 60, 60, 60)),
]

HKL = [hkl for hkl in itertools.product(range(-3, 4), repeat=3) if hkl != (0, 0, 0)]


def main():
    cases = []
    for name, source, system, cell, free, edits, edited in CASES:
        lattice = Lattice.from_parameters(*edited, vesta=True)
        d = [float(lattice.d_hkl(hkl)) for hkl in HKL]
        assert all(np.isfinite(d)) and min(d) > 0
        # The edit keeps what the system fixes: the unedited input parameters of a free group stay, the fixed ones stay.
        cases.append({"name": name, "source": source, "system": system, "inputCell": list(map(float, cell)),
                      "free": free, "edits": {k: float(v) for k, v in edits.items()},
                      "editedCell": list(map(float, edited)), "d": d})
    out = {
        "generator": "tools/pins/lattice-edit/pins.py",
        "date": datetime.date.today().isoformat(),
        "versions": {p: dist_version(p) for p in ("pymatgen", "pymatgen-core", "numpy")},
        "tolerance": {"d_rtol": 1e-7},
        "hkl": [list(h) for h in HKL],
        "cases": cases,
    }
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
        f.write("\n")
    print(f"wrote {os.path.normpath(OUT)}: {len(cases)} cases x {len(HKL)} hkl")


if __name__ == "__main__":
    main()
