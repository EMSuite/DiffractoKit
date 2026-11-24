"""Writes Sources/DiffractoCore/SpaceGroupTable.swift: the settings pymatgen's CifParser expands a CIF with when the file
names its space group by Hermann-Mauguin symbol or IT number and lists no operations, as Hall numbers of spglib.

Run: ../References/python/pymatgen/bin/python -I -B tools/pins/cif-symmetry/space-group-table.py
Versions: spglib 2.8.0, pymatgen-core 2026.10.16 (tools/pins/requirements/pymatgen.txt).

pymatgen-core 2026.10.16 `CifParser.get_symops`, with no operation loop, removes white space and underscores from the
symbol, looks it up among `SYMM_DATA["space_group_encoding"]` (`SpaceGroup(symbol)`), else takes the first entry of
`symmetry/symm_ops.json` whose `hermann_mauguin` (white space removed) is equal; with no symbol it reads the IT number
(`SpaceGroup.from_int_number`). Each of those operation sets is matched here to the one spglib Hall setting with exactly
the same operations (rotations, and translations mod 1 in twelfths); the kit reads that Hall symbol (`HallSymbol`, exact
against all 530). A symbol whose first symm_ops.json entry is no spglib setting (a change of basis such as "C 1",
"A -1", "B 1 21 1") is left out and listed: the kit refuses such a file. No number here is typed.
"""
import json, pathlib, re, warnings
from importlib.metadata import version
import numpy as np
import spglib
import pymatgen.symmetry as pymatgen_symmetry
from monty.serialization import loadfn
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.groups import SYMM_DATA, SpaceGroup

assert version("pymatgen-core") == "2026.10.16"
assert spglib.__version__ == "2.8.0"
root = pathlib.Path(__file__).resolve().parents[3]


def code(rotation, translation):
    scaled = np.asarray(translation, float) % 1.0 * 12
    assert np.allclose(scaled, np.rint(scaled), atol=1e-6), translation
    return tuple(int(v) for v in np.rint(np.asarray(rotation)).ravel()) + tuple(int(v) % 12 for v in np.rint(scaled))


def key(pairs):
    codes = [code(r, t) for r, t in pairs]
    assert len(set(codes)) == len(codes)
    return frozenset(codes)


by_operations = {}
for hall in range(1, 531):
    d = spglib.get_symmetry_from_database(hall)
    by_operations.setdefault(key(zip(d["rotations"], d["translations"])), hall)  # the smallest of equal sets


def hall_of(ops):
    return by_operations.get(key((o.rotation_matrix, o.translation_vector) for o in ops))


def sanitized(symbol):
    return re.sub(r"[\s_]", "", symbol)


symbols = {}
for name in SYMM_DATA["space_group_encoding"]:
    hall = hall_of(SpaceGroup(name).symmetry_ops)
    assert hall is not None, name
    symbols[sanitized(name)] = hall
left_out = []
for entry in loadfn(pathlib.Path(pymatgen_symmetry.__file__).parent / "symm_ops.json"):
    name = re.sub(r"\s+", "", entry["hermann_mauguin"])
    if name in symbols or name in left_out:
        continue  # pymatgen's encoding wins, and the first symm_ops.json entry of a name wins
    hall = hall_of([SymmOp.from_xyz_str(s) for s in entry["symops"]])
    if hall is None:
        left_out.append(name)
    else:
        symbols[name] = hall
numbers = [hall_of(SpaceGroup.from_int_number(n).symmetry_ops) for n in range(1, 231)]
assert None not in numbers
assert all(spglib.get_spacegroup_type(h).number == n for n, h in zip(range(1, 231), numbers))
hall_symbols = [spglib.get_spacegroup_type(h).hall_symbol for h in range(1, 531)]

# pymatgen reads no short symbol ("P 21/c", "Fd-3m" is its own name but "P 42/m n m" is not): with no IT number it
# expands such a file as P1. The kit reads spglib's short and full international symbols after pymatgen's names and the
# number: the setting pymatgen gives the symbol's IT number when it is one of the symbol's settings, else the first.
candidates = {}
for hall in range(1, 531):
    t = spglib.get_spacegroup_type(hall)
    names = {t.international_short, t.international_full, t.international}
    if t.choice:  # "Fd-3m:2", "R-3m:R", "P121/c1:b1": the setting spelled out
        names |= {f"{n}:{t.choice}" for n in list(names)}
    for name in names:
        name = sanitized(name)
        if name not in symbols:
            candidates.setdefault(name, []).append(hall)
spglib_symbols = {}
for name, halls in candidates.items():
    number = spglib.get_spacegroup_type(halls[0]).number
    assert all(spglib.get_spacegroup_type(h).number == number for h in halls), name
    standard = numbers[number - 1]
    spglib_symbols[name] = standard if standard in halls else min(halls)


def swift_string(s):
    assert "\\" not in s
    return '"' + s.replace('"', '\\"') + '"'


def rows(items, per_line):
    return "\n".join("        " + " ".join(items[i:i + per_line]) for i in range(0, len(items), per_line))


source = f'''//
//  SpaceGroupTable.swift
//  Role: the space-group settings pymatgen-core {version("pymatgen-core")} `CifParser` expands a CIF with when the file names its
//        group by Hermann-Mauguin symbol or IT number and lists no operations, as spglib {spglib.__version__} Hall numbers.
//        Generated by tools/pins/cif-symmetry/space-group-table.py; do not edit. {len(symbols)} pymatgen symbols and
//        {len(spglib_symbols)} more of spglib's; the {len(left_out)} symbols whose pymatgen setting is no spglib setting are left out: {", ".join(left_out)}.
//

enum SpaceGroupTable {{
    /// spglib's Hall symbol of each Hall number: index 0 is Hall number 1, index 529 is 530.
    static let hallSymbols: [String] = [
{rows([swift_string(s) + "," for s in hall_symbols], 8)}
    ]

    /// The Hall number of the setting `SpaceGroup.from_int_number(n)` gives: index 0 is IT number 1.
    static let numberSettings: [Int] = [
{rows([f"{h}," for h in numbers], 20)}
    ]

    /// The Hall number of each Hermann-Mauguin symbol pymatgen reads, white space and underscores removed.
    static let symbolSettings: [String: Int] = [
{rows([f"{swift_string(s)}: {h}," for s, h in sorted(symbols.items())], 6)}
    ]

    /// The symbols pymatgen reads in a change of basis that is no spglib setting: a file naming one is refused.
    static let leftOutSymbols: Set<String> = [{", ".join(swift_string(s) for s in left_out)}]

    /// spglib's short and full international symbols pymatgen does not read (white space and underscores removed): the
    /// Hall number of `numberSettings` when it is one of the symbol's settings, else the symbol's first setting.
    static let spglibSymbolSettings: [String: Int] = [
{rows([f"{swift_string(s)}: {h}," for s, h in sorted(spglib_symbols.items())], 6)}
    ]
}}
'''
(root / "Sources/DiffractoCore/SpaceGroupTable.swift").write_text(source)
print(len(symbols), "symbols;", "left out:", left_out, ";", len(spglib_symbols), "spglib symbols")
