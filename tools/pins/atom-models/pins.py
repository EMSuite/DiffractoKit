"""DiffractoKit: the pin generator for the atom-model readers, extended XYZ and POSCAR (DiffractoKit #169, plan row MD1).

Run, from any directory, with the pinned oracle:

    ../References/python/ase/bin/python -I -B tools/pins/atom-models/pins.py

It writes the corpus Tests/DiffractoCoreTests/Fixtures/atom-models/*.{xyz,vasp} and the pins file
Tests/DiffractoCoreTests/Fixtures/atom-models-pins.json: for every file what `ase.io.read` of ASE 3.29.0 returns (atomic
numbers, Cartesian positions, the cell's three rows, pbc, and the number of frames for XYZ), or that the kit refuses it
and by which name.

The corpus:
  ase-*      ASE 3.29.0's own test inputs, verbatim (ase/test/fio/test_extxyz.py, ase/test/fio/vasp/test_poscar.py,
             ase/test/fio/vasp/test_vasp_poscar.py, ase/test/testdata/vasp/POSCAR_example_1, doc/ase/io/172_ext.xyz).
  export-*   written here by ASE 3.29.0's own writers (`ase.io.write`, formats extxyz, xyz and vasp) from `ase.build`.
  layout-*   written here by hand in the layout VESTA 3 and atomsk write (VESTA: File > Export Data > VASP; atomsk: -o
             POSCAR and .exyz). Neither program runs where this script runs; their real exports are DiffractoKit #212.

The numbers: positions and cells the file states in Cartesian Angstrom are pinned exact (the JSON float is the double ASE
holds); converted ones (Direct coordinates, a scale factor) at rtol 1e-7 (assert_allclose's float64 default; ASE's own
reader tests use np.testing.assert_allclose with its defaults).

File:line pins into ase 3.29.0 (site-packages/ase/):
  io/vasp.py      parse_poscar_scaling_factor 35, read_vasp 150, read_vasp_configuration 167
  io/extxyz.py    key_val_str_to_dict 67, parse_properties 338, _read_xyz_frame 390, read_xyz 569
  cell.py         Cell.cartesian_positions 276
"""
import datetime
import io
import json
import os
import textwrap
from importlib.metadata import version as dist_version

import numpy as np
import ase
import ase.io
from ase import Atoms
from ase.build import bulk, molecule
from ase.constraints import FixAtoms

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURES = os.path.join(HERE, "..", "..", "..", "Tests", "DiffractoCoreTests", "Fixtures")
CORPUS = os.path.join(FIXTURES, "atom-models")
OUT = os.path.join(FIXTURES, "atom-models-pins.json")

# --- ASE's own inputs, verbatim -------------------------------------------------------------------------------------

ASE_FILES = {}

# ase/test/testdata/vasp/POSCAR_example_1: a VASP 4 file (no element line); ASE takes the species from an OUTCAR next
# to it. The kit reads no other file, so it refuses this by name.
ASE_FILES["ase-POSCAR_example_1.vasp"] = """\
from OUTCAR_example_1, but removing species
 1.0000000000000000
    17.9343499999999985    0.0000000000000000    0.0000000000000000
     0.0000000000000000   17.9343499999999985    0.0000000000000000
     0.0000000000000000    0.0000000000000000   17.9343499999999985
  18
Cartesian
  0.0000000000000000  0.0000000000000000  0.0000000000000000
 16.8215199999999996 15.7196300000000004 17.6689499999999988
  0.3892300000000000  0.0195000000000000 15.5347100000000005
 16.7979799999999990 15.7855600000000003 15.1218400000000006
 15.8288600000000006  0.1312800000000000 16.6206700000000005
  0.9646700000000000 15.6925000000000008 16.9140900000000016
  2.2704800000000001  0.0633500000000000 16.5629199999999983
 16.3836000000000013  1.9928800000000000  0.0423800000000000
  0.8269100000000000  2.0872000000000002 16.8188999999999993
 16.5645300000000013  1.7942199999999999 15.3298100000000002
 15.6910299999999996 17.8037099999999988  1.3030100000000000
  2.2413099999999999 17.9060200000000016  1.2370300000000001
  2.9509699999999999  2.2455900000000000 17.9213600000000000
 17.4524700000000017  2.1773400000000001  1.0751900000000001
 15.9658200000000008  2.8128700000000002 16.0000900000000002
  1.6050700000000000 15.5837799999999991  1.5149800000000000
 17.0127200000000016 16.5097800000000014  2.6060500000000002
  1.4612400000000001  1.6224600000000000  0.7436700000000000
"""

# ase/test/fio/vasp/test_poscar.py BUF_WITH_HASHES: POTCAR labels with hashes on the element line.
ASE_FILES["ase-poscar-hashes.vasp"] = """\
Na1
   1.0000000000000000
     5.9524044248000987   -0.0007582713987269   -0.0012396575089365
    -0.0014963523464428    6.3481246432503076   -0.0209332994306307
    -0.0014359618019133   -3.1516690207714437    5.5105006246735595
  Na_pv/6a2f546d
               1
Direct
  0.9711350580843501  0.0342911094149940  0.0033526422550807
"""

# ase/test/fio/vasp/test_vasp_poscar.py BUF_NEGATIVE and BUF_MULTIPLE.
ASE_FILES["ase-poscar-negative-scale.vasp"] = """\
 H
-120.00000000000000
     1.0000000000000000    0.0000000000000000    0.0000000000000000
     2.0000000000000000    3.0000000000000000    0.0000000000000000
    -3.0000000000000000    1.0000000000000000    5.0000000000000000
 H
   1
Cartesian
  0.1000000000000000  0.2000000000000000  0.3000000000000000
"""
ASE_FILES["ase-poscar-three-scales.vasp"] = """\
 H
 3.0000000000000000 2.0000000000000000 1.0000000000000000
     1.0000000000000000    0.0000000000000000    0.0000000000000000
     2.0000000000000000    3.0000000000000000    0.0000000000000000
    -3.0000000000000000    1.0000000000000000    5.0000000000000000
 H
   1
Cartesian
  0.1000000000000000  0.2000000000000000  0.3000000000000000
"""

# doc/ase/io/172_ext.xyz (extra integer and string columns, two trailing blank lines).
ASE_FILES["ase-172_ext.xyz"] = """\
7
Lattice="6.0 0.0 0.0 0.0 6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3:molid:I:1:type:S:1 pbc="T T T"
O       1.613900000000000     -0.762100000000000     -0.000000000000000 1 O
C      -0.327900000000000      0.522700000000000      0.000000000000000 1 CT
C       0.392400000000000     -0.722900000000000     -0.000000000000000 1 C
H      -0.960000000000000      0.580900000000000      0.887500000000000 1 HC
H      -0.960000000000000      0.580900000000000     -0.887500000000000 1 HC
H       0.346400000000000      1.382000000000000      0.000000000000000 1 HC
H      -0.104900000000000     -1.581400000000000     -0.000000000000000 1 H1


"""

# ase/test/fio/test_extxyz.py: test_vec_cell (VEC lines, whitespace), the upper-case and trailing-blank-line file,
# test_read_slash, test_blank_comment, test_constraints_int and three of test_pbc_property's cases.
ASE_FILES["ase-vec-cell.xyz"] = """1
    Coordinates
    C         -7.28250        4.71303       -3.82016
      VEC1 1.0 0.1 1.1
    1

    C         -7.28250        4.71303       -3.82016
    VEC1 1.0 0.1 1.1
    """
ASE_FILES["ase-uppercase-blank-line.xyz"] = """4
    Coordinates
    MG        -4.25650        3.79180       -2.54123
    C         -1.15405        2.86652       -1.26699
    C         -5.53758        3.70936        0.63504
    C         -7.28250        4.71303       -3.82016

    """
ASE_FILES["ase-slash.xyz"] = """4
    key1=a key2=a/b key3=a@b key4="a@b"
    Mg        -4.25650        3.79180       -2.54123
    C         -1.15405        2.86652       -1.26699
    C         -5.53758        3.70936        0.63504
    C         -7.28250        4.71303       -3.82016
    """
ASE_FILES["ase-blank-comment.xyz"] = """4

    Mg        -4.25650        3.79180       -2.54123
    C         -1.15405        2.86652       -1.26699
    C         -5.53758        3.70936        0.63504
    C         -7.28250        4.71303       -3.82016
    """
ASE_FILES["ase-move-mask.xyz"] = """3
Properties=species:S:1:pos:R:3:move_mask:I:1 pbc="F F F"
O        0.00000000       0.00000000       0.11926200  1
H        0.00000000       0.76323900      -0.47704700  0
H        0.00000000      -0.76323900      -0.47704700  0"""
for i, pbc in enumerate(["False false T", "True", "T T F"]):
    ASE_FILES[f"ase-pbc-{i + 1}.xyz"] = f"""2
Lattice="3.608 0.0 0.0 -1.804 3.125 0.0 0.0 0.0 21.3114930844" pbc="{pbc}"
As           1.8043384632       1.0417352974      11.3518747709
As          -0.0000000002       2.0834705948       9.9596183135"""

# --- exports by ASE's own writers -----------------------------------------------------------------------------------


def export(atoms_or_images, fmt, **kw):
    buf = io.StringIO()
    ase.io.write(buf, atoms_or_images, format=fmt, **kw)
    return buf.getvalue()


def exports():
    files = {}
    si = bulk("Si", "diamond", a=5.4309, cubic=True).repeat((2, 2, 2))
    files["export-si-cubic-222.xyz"] = export(si, "extxyz")
    mg = bulk("Mg", "hcp", a=3.2094, c=5.2108).repeat((3, 3, 2))
    files["export-mg-hcp-332.xyz"] = export(mg, "extxyz")
    water = molecule("H2O")
    files["export-h2o-no-cell.xyz"] = export(water, "extxyz")
    files["export-si-plain.xyz"] = export(bulk("Si", "diamond", a=5.4309, cubic=True), "xyz")
    frames = []
    for k in range(3):
        cu = bulk("Cu", "fcc", a=3.615, cubic=True).repeat((k + 1, 1, 1))
        cu.rattle(stdev=0.01, seed=k)
        frames.append(cu)
    files["export-cu-three-frames.xyz"] = export(frames, "extxyz")
    sto = Atoms("SrTiO3", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)],
                cell=[3.905, 3.905, 3.905], pbc=True).repeat((2, 2, 2))
    files["export-srtio3-direct.vasp"] = export(sto, "vasp", direct=True)
    files["export-srtio3-cartesian.vasp"] = export(sto, "vasp", direct=False, sort=True)
    al = bulk("Al", "fcc", a=4.0495, cubic=True).repeat((3, 3, 3))
    del al[13]
    al.set_constraint(FixAtoms(indices=range(0, 20)))
    files["export-al-vacancy-selective.vasp"] = export(al, "vasp", direct=True)
    tri = Atoms("GaAs", scaled_positions=[(0, 0, 0), (0.25, 0.25, 0.25)],
                cell=[[0, 2.8265, 2.8265], [2.8265, 0, 2.8265], [2.8265, 2.8265, 0]], pbc=True)
    files["export-gaas-primitive.vasp"] = export(tri, "vasp", direct=True)
    return files


# --- the layouts of VESTA and atomsk, by hand ------------------------------------------------------------------------

LAYOUT_FILES = {
    # VESTA 3, File > Export Data > VASP (POSCAR), "Fractional coordinates".
    "layout-vesta-si.vasp": textwrap.dedent("""\
        New structure
           1.0
                5.4309000969         0.0000000000         0.0000000000
                0.0000000000         5.4309000969         0.0000000000
                0.0000000000         0.0000000000         5.4309000969
           Si
            8
        Direct
             0.000000000         0.000000000         0.000000000
             0.000000000         0.500000000         0.500000000
             0.500000000         0.000000000         0.500000000
             0.500000000         0.500000000         0.000000000
             0.750000000         0.250000000         0.750000000
             0.250000000         0.250000000         0.250000000
             0.250000000         0.750000000         0.750000000
             0.750000000         0.750000000         0.250000000
        """),
    # atomsk -o POSCAR: a comment line, the scale 1.0, Cartesian.
    "layout-atomsk-al.vasp": textwrap.dedent("""\
        # Al fcc
             1.00000000
             4.04950000     0.00000000     0.00000000
             0.00000000     4.04950000     0.00000000
             0.00000000     0.00000000     4.04950000
          Al
             4
        Cartesian
             0.00000000     0.00000000     0.00000000
             2.02475000     2.02475000     0.00000000
             0.00000000     2.02475000     2.02475000
             2.02475000     0.00000000     2.02475000
        """),
    # atomsk -o exyz: extended XYZ with Lattice, Properties and pbc.
    "layout-atomsk-al.xyz": textwrap.dedent("""\
        4
        Lattice="4.04950000 0.00000000 0.00000000 0.00000000 4.04950000 0.00000000 0.00000000 0.00000000 4.04950000" Properties=species:S:1:pos:R:3 pbc="T T T"
        Al     0.00000000     0.00000000     0.00000000
        Al     2.02475000     2.02475000     0.00000000
        Al     0.00000000     2.02475000     2.02475000
        Al     2.02475000     0.00000000     2.02475000
        """),
}

# What the kit refuses although ASE reads it, by the refusal's name (AtomModel.ReadError.name in Swift).
KIT_REFUSES = {
    "ase-POSCAR_example_1.vasp": "noElementLine",
}

# Whether a number of the file is stated in Cartesian Angstrom (exact) or converted (rtol 1e-7).
CONVERTED_POSITIONS = {"ase-poscar-hashes.vasp", "ase-poscar-negative-scale.vasp", "ase-poscar-three-scales.vasp",
                       "export-srtio3-direct.vasp", "export-al-vacancy-selective.vasp", "export-gaas-primitive.vasp",
                       "layout-vesta-si.vasp"}
CONVERTED_CELLS = {"ase-poscar-negative-scale.vasp", "ase-poscar-three-scales.vasp"}


def frame_count(text):
    frames = list(ase.io.iread(io.StringIO(text), index=":", format="extxyz"))
    return len(frames)


def main():
    os.makedirs(CORPUS, exist_ok=True)
    files = dict(ASE_FILES)
    files.update(exports())
    files.update(LAYOUT_FILES)
    cases = []
    for name in sorted(files):
        text = files[name]
        with open(os.path.join(CORPUS, name), "w", newline="\n") as f:
            f.write(text)
        fmt = "extxyz" if name.endswith(".xyz") else "vasp"
        case = {"file": name, "format": "xyz" if fmt == "extxyz" else "poscar"}
        if name in KIT_REFUSES:
            case["refusal"] = KIT_REFUSES[name]
            cases.append(case)
            continue
        atoms = ase.io.read(io.StringIO(text), format=fmt)
        case["numbers"] = [int(z) for z in atoms.numbers]
        case["positions"] = [[float(v) for v in p] for p in atoms.positions]
        case["cell"] = [[float(v) for v in row] for row in np.array(atoms.cell)]
        case["pbc"] = [bool(b) for b in atoms.pbc]
        case["frames"] = frame_count(text) if fmt == "extxyz" else 1
        case["positionsExact"] = name not in CONVERTED_POSITIONS
        case["cellExact"] = name not in CONVERTED_CELLS
        cases.append(case)
    out = {
        "generator": "tools/pins/atom-models/pins.py",
        "issue": "EMSuite/DiffractoKit#169",
        "ase": dist_version("ase"),
        "numpy": dist_version("numpy"),
        "rtol": 1e-7,
        "date": datetime.date.today().isoformat(),
        "cases": cases,
    }
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
        f.write("\n")
    print(f"wrote {len(cases)} cases to {os.path.relpath(OUT)}")


if __name__ == "__main__":
    main()
