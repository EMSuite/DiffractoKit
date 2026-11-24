"""orix pins for the variant orientations of an orientation relationship (DiffractoKit S7, #162).

  ../References/python/orix/bin/python -I -B tools/pins/orix/variant-pins.py

Three planted relations of a precipitate in aluminium, each stated as two parallel pairs the way a paper writes it:
  * beta'' (Mg5Si6, C2/m, a 15.16, b 4.05, c 6.74 A, beta 105.3 deg; Andersen et al., Acta Mater. 46 (1998) 3283):
    [100] || [3 2 0]Al first, [001] || [1 -3 0]Al second. The second pair is 0.045 deg from parallel (105.3 against
    arccos(-3/sqrt(130)) = 105.2551 deg), which the kit returns as the misfit.
  * theta' (Al2Cu, tetragonal, a 4.04, c 5.80 A; Silcock, Heal and Hardy, J. Inst. Met. 82 (1953) 239): (001) || (001)Al,
    [100] || [100]Al.
  * T1 (Al2CuLi, hexagonal, a 4.965, c 9.345 A; Huang and Ardell, Mater. Sci. Technol. 3 (1987) 176): (0001) || (111)Al,
    [2 1 0] (= <10-10> in four indices) || [1 -1 0]Al.
The rotation R (precipitate Cartesian to Al Cartesian, the kit's UnitCell frames: a along x, b in the xy plane) is the
triad construction of `OrientationRelationship.rotation`, rewritten here in numpy. orix 0.15.0 `Misorientation(R,
symmetry=(Gp, Gm)).equivalent()` gives the set Gm R Gp; the classes under Gp on the right are the variants. Written per
case: R, the distinct equivalents, the class count, the misfit of the second pair. Tolerance atol 1e-4
(orix test_misorientation.py); the class count and the misfit are exact up to round-off.
"""
import importlib.metadata as md
import importlib.util
import json
from pathlib import Path
import numpy as np
from orix.quaternion import Misorientation

KIT = Path(__file__).resolve().parents[3]
FIX = KIT / 'Tests/DiffractoOrientationTests/Fixtures'
assert md.version('orix') == '0.15.0'
spec = importlib.util.spec_from_file_location('mis', Path(__file__).with_name('misorientation-pins.py'))
mis = importlib.util.module_from_spec(spec); spec.loader.exec_module(mis)


def lattice(a, b, c, al=90.0, be=90.0, ga=90.0):
    """Rows a, b, c (real) and a*, b*, c* (reciprocal, no 2 pi): the kit's UnitCell frame."""
    al, be, ga = np.radians([al, be, ga])
    cx = c * np.cos(be); cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    real = np.array([[a, 0, 0], [b * np.cos(ga), b * np.sin(ga), 0], [cx, cy, np.sqrt(c * c - cx * cx - cy * cy)]])
    return real, np.linalg.inv(real).T


def vec(rows, kind, idx):
    real, rec = rows
    return np.asarray(idx, float) @ (real if kind == 'direction' else rec)


def triad(a, b):
    e1 = a / np.linalg.norm(a); r = b - (b @ e1) * e1; e2 = r / np.linalg.norm(r)
    return np.column_stack([e1, e2, np.cross(e1, e2)])


def angle(a, b):
    return float(np.arctan2(np.linalg.norm(np.cross(a, b)), a @ b))


AL = lattice(4.05, 4.05, 4.05)
CASES = [
    ('beta-double-prime', '2/m', lattice(15.16, 4.05, 6.74, be=105.3),
     (('direction', [1, 0, 0]), ('direction', [3, 2, 0])), (('direction', [0, 0, 1]), ('direction', [1, -3, 0])), 12),
    ('theta-prime', '4/mmm', lattice(4.04, 4.04, 5.80),
     (('plane', [0, 0, 1]), ('plane', [0, 0, 1])), (('direction', [1, 0, 0]), ('direction', [1, 0, 0])), 3),
    ('t1', '6/mmm', lattice(4.965, 4.965, 9.345, ga=120.0),
     (('plane', [0, 0, 1]), ('plane', [1, 1, 1])), (('direction', [2, 1, 0]), ('direction', [1, -1, 0])), 4),
]


def main():
    out = {'orix': md.version('orix'), 'cases': []}
    sm, mops = mis.symmetry('m-3m')
    for name, psym, prows, first, second, expected in CASES:
        p1, m1 = vec(prows, *first[0]), vec(AL, *first[1])
        p2, m2 = vec(prows, *second[0]), vec(AL, *second[1])
        R = triad(m1, m2) @ triad(p1, p2).T
        sp, pops = mis.symmetry(psym)
        m = Misorientation(mis.rot(R).data); m.symmetry = (sp.proper_subgroup, sm.proper_subgroup)
        eq = mis.mats(m.equivalent())
        uniq = []
        for e in eq:
            if not any(np.allclose(e, u, atol=1e-9) for u in uniq):
                uniq.append(e)
        # orix's set is Gm R Gp: assert it from the kit's operators, then count the classes under Gp.
        brute = [g @ R @ s for g in mops for s in pops]
        assert all(any(np.abs(e - b).max() < 1e-9 for b in brute) for e in uniq) and len(uniq) == len({tuple(np.round(b, 9).ravel()) for b in brute})
        classes = []
        for e in uniq:
            if not any(any(np.abs(c @ s - e).max() < 1e-9 for s in pops) for c in classes):
                classes.append(e)
        assert len(classes) == expected, (name, len(classes))
        out['cases'].append({'name': name, 'precipitateSymmetry': psym, 'R': R.tolist(), 'equivalents': [u.tolist() for u in uniq],
                             'variantCount': len(classes), 'secondPairMisfit': angle(R @ p2, m2)})
    (FIX / 'orix-variant-pins.json').write_text(json.dumps(out, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n')
    print('ok', [(c['name'], len(c['equivalents']), c['variantCount'], np.degrees(c['secondPairMisfit'])) for c in out['cases']])


if __name__ == '__main__':
    main()
