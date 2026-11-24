"""Pins for the CIF that names its space group by Hermann-Mauguin symbol or IT number and lists no operations (K6).

Run: ../References/python/pymatgen/bin/python -I -B tools/pins/cif-symmetry/space-group-pins.py
Versions: spglib 2.8.0, pymatgen-core 2026.10.16 (tools/pins/requirements/pymatgen.txt).

Writes Tests/DiffractoCoreTests/Fixtures/space-group-symbol-pins.json:
  symbols: every symbol of Sources/DiffractoCore/SpaceGroupTable.swift with the operations pymatgen gives it (SpaceGroup,
           or the first symm_ops.json entry, as CifParser.get_symops), or for a symbol pymatgen does not read, spglib's
           operations of the setting the table names; numbers: SpaceGroup.from_int_number(n). Each set is the first 16 hex
           digits of the SHA-256 of its sorted 12-character codes (hall-symbol-pins.json's encoding), concatenated.
  cifs:    planted CIFs with no operation loop: pymatgen CifParser.parse_structures(primitive=False) cell, sites and
           warnings; where pymatgen does not read the symbol (it then builds P1), the atom count spglib's operations of
           the table's setting give, and the kit's documented answer.
"""
import hashlib, json, pathlib, re, warnings
from importlib.metadata import version
import numpy as np
import spglib
import pymatgen.symmetry as pymatgen_symmetry
from monty.serialization import loadfn
from pymatgen.core.operations import SymmOp
from pymatgen.io.cif import CifParser, str2float
from pymatgen.symmetry.groups import SYMM_DATA, SpaceGroup

assert version("pymatgen-core") == "2026.10.16"
assert spglib.__version__ == "2.8.0"
root = pathlib.Path(__file__).resolve().parents[3]
out = root / "Tests/DiffractoCoreTests/Fixtures/space-group-symbol-pins.json"


def encode(r, t):
    rot = "".join("-0+"[int(round(v)) + 1] for v in np.asarray(r).ravel())
    scaled = np.asarray(t, float) % 1.0 * 12
    assert np.allclose(scaled, np.rint(scaled), atol=1e-6), t
    return rot + "".join("0123456789ab"[int(v) % 12] for v in np.rint(scaled))


def digest(pairs):
    codes = sorted(encode(r, t) for r, t in pairs)
    assert len(set(codes)) == len(codes)
    return hashlib.sha256("".join(codes).encode()).hexdigest()[:16]


def of_ops(ops):
    return digest((o.rotation_matrix, o.translation_vector) for o in ops)


def of_hall(hall):
    d = spglib.get_symmetry_from_database(hall)
    return digest(zip(d["rotations"], d["translations"]))


def sanitized(s):
    return re.sub(r"[\s_]", "", s)


# pymatgen's names, enumerated here from pymatgen itself (not from the kit's table): the encoding first, then the first
# symm_ops.json entry of each other name. A name whose operations are no spglib setting is listed as left out.
symbols = {}
left_out = []
encoding = {sanitized(k): k for k in SYMM_DATA["space_group_encoding"]}
spglib_sets = {digest(zip(*(lambda d: (d["rotations"], d["translations"]))(spglib.get_symmetry_from_database(h))))
               for h in range(1, 531)}
for name, key in encoding.items():
    symbols[name] = of_ops(SpaceGroup(key).symmetry_ops)
for entry in loadfn(pathlib.Path(pymatgen_symmetry.__file__).parent / "symm_ops.json"):
    name = re.sub(r"\s+", "", entry["hermann_mauguin"])
    if name in symbols or name in left_out:
        continue
    d = of_ops([SymmOp.from_xyz_str(s) for s in entry["symops"]])
    if d in spglib_sets:
        symbols[name] = d
    else:
        left_out.append(name)
# spglib's names pymatgen does not read: the setting the kit's table names (read from the generated source).
table = (root / "Sources/DiffractoCore/SpaceGroupTable.swift").read_text()
section = table.split("static let spglibSymbolSettings")[1].split("]\n")[0]
for name, hall in re.findall(r'"((?:[^"\\]|\\.)*)": (\d+)', section):
    name = name.replace('\\"', '"')
    assert name not in symbols, name
    symbols[name] = of_hall(int(hall))
numbers = [of_ops(SpaceGroup.from_int_number(n).symmetry_ops) for n in range(1, 231)]


def cif(tags, cell, sites):
    lines = ["data_t"] + [f"{t} {v}" for t, v in tags]
    for k, v in zip(["length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"], cell):
        lines.append(f"_cell_{k} {v}")
    lines += ["loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y",
              "_atom_site_fract_z"]
    lines += [f"{el}{i + 1} {el} {f[0]} {f[1]} {f[2]}" for i, (el, f) in enumerate(sites)]
    return "\n".join(lines) + "\n"


def spglib_atoms(hall, sites):
    d = spglib.get_symmetry_from_database(hall)
    pts = []
    for el, f in sites:
        for r, t in zip(d["rotations"], d["translations"]):
            p = (r @ np.array(f, float) + t) % 1.0
            if not any(e == el and np.allclose(((q - p + 0.5) % 1.0) - 0.5, 0, atol=1e-4) for e, q in pts):
                pts.append((el, p))
    return len(pts)


HM = "_symmetry_space_group_name_H-M"
NUM = "_symmetry_Int_Tables_number"
specs = [
    # (name, tags, cell, sites, kit: "pymatgen" | ("spglib", hall) | "refused")
    ("rock salt by symbol", [(HM, "'F m -3 m'")], (5.6402,) * 3 + (90, 90, 90),
     [("Na", (0, 0, 0)), ("Cl", (0.5, 0.5, 0.5))], "pymatgen"),
    ("silicon by number", [(NUM, "227")], (5.4309,) * 3 + (90, 90, 90), [("Si", (0, 0, 0))], "pymatgen"),
    ("magnesium by symbol and number", [(HM, "'P 63/m m c'"), (NUM, "194")], (3.2094, 3.2094, 5.2108, 90, 90, 120),
     [("Mg", (0.3333333, 0.6666667, 0.25))], "pymatgen"),
    ("rutile by symbol with underscore", [(HM, "'P 4_2/m n m'")], (4.594, 4.594, 2.959, 90, 90, 90),
     [("Ti", (0, 0, 0)), ("O", (0.3053, 0.3053, 0))], "pymatgen"),
    ("baddeleyite by full symbol", [(HM, "'P 1 21/c 1'")], (5.1505, 5.2116, 5.3173, 90, 99.23, 90),
     [("Zr", (0.2754, 0.0395, 0.2083)), ("O", (0.07, 0.3317, 0.3447)), ("O", (0.4496, 0.7569, 0.4792))], "pymatgen"),
    ("P 32 2 1 general positions by number", [(NUM, "154")], (4.9134, 4.9134, 5.4052, 90, 90, 120),
     [("Si", (0.4697, 0, 0)), ("O", (0.4135, 0.2669, 0.1191))], "pymatgen"),
    ("corundum by symbol, hexagonal axes", [(HM, "'R -3 c'")], (4.7589, 4.7589, 12.991, 90, 90, 120),
     [("Al", (0, 0, 0.35216)), ("O", (0.30624, 0, 0.25))], "pymatgen"),
    ("baddeleyite by short symbol alone", [(HM, "'P 21/c'")], (5.1505, 5.2116, 5.3173, 90, 99.23, 90),
     [("Zr", (0.2754, 0.0395, 0.2083)), ("O", (0.07, 0.3317, 0.3447)), ("O", (0.4496, 0.7569, 0.4792))],
     ("spglib", 81)),
    ("silicon by symbol, origin choice 2", [(HM, "'F d -3 m :2'")], (5.4309,) * 3 + (90, 90, 90),
     [("Si", (0.125, 0.125, 0.125))], ("spglib", 526)),
    ("C 1, a setting the kit does not hold", [(HM, "'C 1'")], (5, 6, 7, 90, 100, 90), [("Si", (0.1, 0.2, 0.3))],
     "refused"),
    ("A -1 beside its number, a setting the kit does not hold", [(HM, "'A -1'"), (NUM, "2")], (5, 6, 7, 80, 100, 95),
     [("Si", (0.1, 0.2, 0.3))], "refused"),
]
cifs = []
for name, tags, cell, sites, kit in specs:
    text = cif(tags, cell, sites)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parser = CifParser.from_str(text)
        s = parser.parse_structures(primitive=False, symmetrized=False)[0]
    # The cell as pymatgen reads it from the text (str2float of the six values) before Lattice.from_parameters:
    # lattice.gamma comes back from the matrix 1 ulp off 120.
    block = parser._cif.data["t"]
    written = [str2float(block[f"_cell_{k}"]) for k in
               ["length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"]]
    assert np.allclose(written, [s.lattice.a, s.lattice.b, s.lattice.c, s.lattice.alpha, s.lattice.beta, s.lattice.gamma],
                       rtol=1e-12, atol=0)
    entry = {"name": name, "text": text, "pymatgenWarnings": list(parser.warnings),
             "pymatgen": {"cell": written,
                          "sites": [{"z": site.specie.Z, "fractional": (site.frac_coords % 1.0).tolist()} for site in s]}}
    if kit == "pymatgen":
        entry["kit"] = "pymatgen"
    elif kit == "refused":
        entry["kit"] = "refused"
        entry["refusedAs"] = dict(tags)[HM].strip("'")
    else:
        entry["kit"] = "spglib"
        entry["hall"] = spglib.get_spacegroup_type(kit[1]).hall_symbol
        entry["spglibAtoms"] = spglib_atoms(kit[1], sites)
    cifs.append(entry)

json.dump({"spglib": spglib.__version__, "pymatgen-core": version("pymatgen-core"), "symbols": symbols, "leftOut": left_out,
           "numbers": numbers, "cifs": cifs}, open(out, "w"), indent=0, sort_keys=True)
out.write_text(out.read_text() + "\n")
print(out.name, len(symbols), "symbols", len(numbers), "numbers")
for c in cifs:
    print(c["name"], len(c["pymatgen"]["sites"]), c["kit"], c.get("spglibAtoms"), c["pymatgenWarnings"])
