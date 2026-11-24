"""Write Sources/DiffractoScattering/XRayScatteringTable.swift from pymatgen's atomic_scattering_params.json.

    ../References/python/pymatgen/bin/python -I -B tools/pins/powder-xrd/make-table.py

The kit's copy is carried by this script, never by hand: every number is written as Python's shortest round-trip repr of
the double the JSON decodes to, which Swift parses to the same double. Rows are by Z (1 to 98); D, the table's second
Z = 1 entry, equals H and is not repeated. PowderXRDPinsTests holds the copy to the pinned table bit for bit.
"""
import hashlib
import json
import os

import pymatgen.analysis.diffraction.xrd as xrd_module
from pymatgen.core import Element

HERE = os.path.dirname(os.path.abspath(__file__))
KIT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.join(KIT, "Sources", "DiffractoScattering", "XRayScatteringTable.swift")

path = os.path.join(os.path.dirname(xrd_module.__file__), "atomic_scattering_params.json")
raw = open(path, "rb").read()
table = json.loads(raw)
rows = {}
for sym, pairs in table.items():
    z = Element(sym).Z
    flat = [float(v) for pair in pairs for v in pair]
    assert len(flat) == 8, sym
    if z in rows:
        assert rows[z][1] == flat, f"{sym} differs from {rows[z][0]}"
        continue
    rows[z] = (sym, flat)
assert sorted(rows) == list(range(1, 99))

lines = [
    "//",
    "//  XRayScatteringTable.swift",
    "//  Role: pymatgen's X-ray atomic scattering parameters (pymatgen/analysis/diffraction/atomic_scattering_params.json,",
    f"//        sha256 {hashlib.sha256(raw).hexdigest()[:16]}…), four (a, b) pairs per element, Z 1 to 98.",
    "//        WRITTEN BY tools/pins/powder-xrd/make-table.py; do not edit by hand.",
    "//",
    "",
    "extension XRayScatteringFactors {",
    "",
    "    /// Row Z − 1 = [a0, b0, a1, b1, a2, b2, a3, b3] (s in Å⁻¹, s = sin θ / λ).",
    "    public static let parameters: [[Double]] = [",
]
for z in range(1, 99):
    sym, flat = rows[z]
    lines.append("        [" + ", ".join(repr(v) for v in flat) + f"],  // {z:3d} {sym}")
lines += ["    ]", "}", ""]
open(OUT, "w", encoding="utf-8").write("\n".join(lines))
print(f"wrote {OUT}")
