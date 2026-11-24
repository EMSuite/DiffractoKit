"""Pins for naming the space-group type from operations (K7).

Run: ../References/python/pymatgen/bin/python -I -B tools/pins/cif-symmetry/space-group-type-pins.py
Versions: spglib 2.8.0 (tools/pins/requirements/pymatgen.txt).

Writes Tests/DiffractoCoreTests/Fixtures/space-group-type-pins.json: for each of the 530 Hall settings spglib lists, a
seeded random origin shift p, and what spglib.get_spacegroup_type_from_symmetry names the setting's operations shifted
by it (t' = t + (I - R) p mod 1; lattice the identity, symprec 1e-5): IT number, Hall number, short international
symbol. The kit rebuilds the same operations from the Hall symbol (exact, P37) and the shift.
"""
import json, pathlib
import numpy as np
import spglib

assert spglib.__version__ == "2.8.0"
out = pathlib.Path(__file__).resolve().parents[3] / "Tests/DiffractoCoreTests/Fixtures/space-group-type-pins.json"
rng = np.random.default_rng(20261123)
settings = []
for hall in range(1, 531):
    d = spglib.get_symmetry_from_database(hall)
    shift = rng.random(3)
    translations = np.array([(t + (np.eye(3) - r) @ shift) % 1.0 for r, t in zip(d["rotations"], d["translations"])])
    named = spglib.get_spacegroup_type_from_symmetry(d["rotations"], translations, lattice=np.eye(3), symprec=1e-5)
    settings.append({"hall": hall, "hallSymbol": spglib.get_spacegroup_type(hall).hall_symbol, "shift": shift.tolist(),
                     "number": named.number, "hallNumber": named.hall_number,
                     "international": named.international_short})
json.dump({"spglib": spglib.__version__, "seed": 20261123, "settings": settings}, open(out, "w"), indent=0)
out.write_text(out.read_text() + "\n")
print(out.name, len(settings), len({s["number"] for s in settings}), "types")
