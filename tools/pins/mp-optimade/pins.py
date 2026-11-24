"""pymatgen's reading of the recorded keyless Materials Project entries (DiffractoKit #150, plan row K3a).

Run with the pymatgen environment: python -I -B tools/pins/mp-optimade/pins.py (tools/pins/regenerate.sh does).
Input: the entries tools/materials-project/record.py recorded byte for byte from optimade.materialsproject.org into
Tests/DiffractoCoreTests/Fixtures/mp-optimade/ (no network here). For each entry:

  served       pymatgen `OptimadeStructureAdapter.get_structure`: the served cell's parameters and every site's species
               and fractional coordinates (the Cartesian positions in the served lattice)
  conventional `SpacegroupAnalyzer(symprec=1e-3, angle_tolerance=-1).get_conventional_standard_structure()`: the
               Bravais lattice, the Laue class, the space group, the cell parameters, the volume and the site count

symprec 1e-3 Å is a thousand times the server's rounding (1e-6 Å on every coordinate) and far below any real
displacement. The recorder kept only entries whose space group is the same at 1e-4, 1e-3 and 1e-2 Å, and this script
asserts that pymatgen agrees, so the pins do not hang on the tolerance. angle_tolerance -1 leaves the metric to symprec
(spglib's own default, and the recorder's): pymatgen's default of 5 degrees reads mp-623511 (Sn, an I cell with edges
3.45 and 3.55 Å) as tetragonal, which the kit's metric tolerance of 1e-5 does not. No expected value is calculated by the
Swift port.
"""
from importlib.metadata import version
from pathlib import Path
import json
import platform

import numpy as np
import spglib
from pymatgen.io.optimade import OptimadeStructureAdapter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

assert version('pymatgen') == '2026.9.24'
assert version('pymatgen-core') == '2026.10.16'
assert np.__version__ == '2.5.4'
assert spglib.__version__ == '2.8.0'

ROOT = Path(__file__).resolve().parents[3] / 'Tests/DiffractoCoreTests/Fixtures/mp-optimade'
SYMPREC = 1e-3
ANGLE_TOLERANCE = -1
LAUE = {'1': '-1', '-1': '-1', '2': '2/m', 'm': '2/m', '2/m': '2/m', '222': 'mmm', 'mm2': 'mmm', 'mmm': 'mmm',
        '4': '4/m', '-4': '4/m', '4/m': '4/m', '422': '4/mmm', '4mm': '4/mmm', '-42m': '4/mmm', '4/mmm': '4/mmm',
        '3': '-3', '-3': '-3', '32': '-3m', '3m': '-3m', '-3m': '-3m', '6': '6/m', '-6': '6/m', '6/m': '6/m',
        '622': '6/mmm', '6mm': '6/mmm', '-6m2': '6/mmm', '6/mmm': '6/mmm', '23': 'm-3', 'm-3': 'm-3',
        '432': 'm-3m', '-43m': 'm-3m', 'm-3m': 'm-3m'}
FAMILY = {'triclinic': 'a', 'monoclinic': 'm', 'orthorhombic': 'o', 'tetragonal': 't', 'trigonal': 'h',
          'hexagonal': 'h', 'cubic': 'c'}


def parameters(lattice):
    return [lattice.a, lattice.b, lattice.c, lattice.alpha, lattice.beta, lattice.gamma]


def bravais(analyzer):
    centring = analyzer.get_space_group_symbol()[0]
    centring = 'C' if centring in 'ABC' else centring
    family = FAMILY[analyzer.get_crystal_system()]
    if family == 'h' and centring == 'R':
        return 'hR'
    return family + centring


cases = []
for path in sorted(ROOT.glob('entry-*.json')):
    resource = json.loads(path.read_bytes())['data']
    served = OptimadeStructureAdapter.get_structure(resource)
    analyzer = SpacegroupAnalyzer(served, symprec=SYMPREC, angle_tolerance=ANGLE_TOLERANCE)
    conventional = analyzer.get_conventional_standard_structure()
    # The recorder kept only entries whose space group does not depend on the tolerance; hold pymatgen to that too.
    for symprec in (1e-4, 1e-2):
        other = SpacegroupAnalyzer(served, symprec=symprec, angle_tolerance=ANGLE_TOLERANCE).get_space_group_number()
        assert other == analyzer.get_space_group_number(), f'{path.name}: space group {other} at symprec {symprec}'
    cases.append({
        'file': path.name,
        'id': resource['id'],
        'served': {'parameters': parameters(served.lattice),
                   'sites': [{'element': site.specie.symbol, 'fractional': site.frac_coords.tolist()} for site in served]},
        'bravais': bravais(analyzer),
        'laue': LAUE[analyzer.get_point_group_symbol()],
        'spaceGroupNumber': analyzer.get_space_group_number(),
        'spaceGroupSymbol': analyzer.get_space_group_symbol(),
        'conventional': {'parameters': parameters(conventional.lattice), 'volume': conventional.volume,
                         'siteCount': len(conventional)},
    })

out = {'generator': 'tools/pins/mp-optimade/pins.py', 'python': platform.python_version(), 'pymatgen': version('pymatgen'),
       'pymatgenCore': version('pymatgen-core'), 'spglib': spglib.__version__, 'numpy': np.__version__,
       'symprec': SYMPREC, 'angleTolerance': ANGLE_TOLERANCE, 'cases': cases}
(ROOT / 'pymatgen-pins.json').write_text(json.dumps(out, sort_keys=True, indent=1) + '\n')
print(f'{len(cases)} entries pinned')
for case in cases:
    print(case['id'], case['bravais'], case['spaceGroupSymbol'], case['laue'], case['conventional']['siteCount'])
