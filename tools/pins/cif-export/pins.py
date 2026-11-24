"""Pin C1: pymatgen-core 2026.10.16 `CifParser` reads each CIF DiffractoKit's `CIFExport` wrote (the corpus of
Tests/DiffractoKitTestSupport/CIFExportCorpus.swift, committed byte for byte under Fixtures/cif-export/written/): the cell
and every expanded site as pymatgen builds them. Run with References/python/pymatgen/bin/python -I -B from the repo root.
The kit's texts come from `CIFExportPymatgenTests` (it prints a text it does not find committed); nothing here is
computed by the Swift port. `parse_structures(primitive=False)`: the cell as written, no reduction.
"""
from pathlib import Path
from importlib.metadata import version
import hashlib, json, platform, warnings
import numpy as np
from pymatgen.io.cif import CifParser
import spglib
assert version('pymatgen-core') == '2026.10.16'
assert spglib.__version__ == '2.8.0'
root = Path(__file__).resolve().parents[3] / 'Tests/DiffractoCoreTests/Fixtures/cif-export'

structures = []
for path in sorted((root / 'written').glob('*.cif')):
    text = path.read_text(encoding='utf-8')
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        s = CifParser.from_str(text).parse_structures(primitive=False)[0]
    sites = []
    for site in s:
        (element, occupancy), = site.species.items()
        sites.append(dict(z=element.Z, frac=[float(v) for v in site.frac_coords], occupancy=float(occupancy)))
    a, b, c = s.lattice.abc
    alpha, beta, gamma = s.lattice.angles
    structures.append(dict(name=path.stem, cifSha256=hashlib.sha256(text.encode('utf-8')).hexdigest(),
                           cell=dict(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma), sites=sites,
                           warnings=sorted({str(w.message) for w in caught})))
out = dict(python=platform.python_version(), pymatgenCore=version('pymatgen-core'), spglib=spglib.__version__,
           numpy=np.__version__, structures=structures)
(root / 'pymatgen-read-pins.json').write_text(json.dumps(out, sort_keys=True, indent=2) + '\n')
print(f'{len(structures)} CIFs read by pymatgen')
for st in structures:
    print(st['name'], len(st['sites']), st['warnings'])
