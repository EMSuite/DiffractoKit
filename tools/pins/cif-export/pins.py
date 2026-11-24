"""Pin C1: pymatgen-core 2026.10.16 `CifParser` reads each CIF DiffractoKit's `CIFExport` wrote (the corpus of
Tests/DiffractoKitTestSupport/CIFExportCorpus.swift, committed byte for byte under Fixtures/cif-export/written/): the cell
and every expanded site as pymatgen builds them. Run with References/python/pymatgen/bin/python -I -B from the repo root.
The kit's texts come from `CIFExportPymatgenTests` (it prints a text it does not find committed); nothing here is
computed by the Swift port. `parse_structures(primitive=False)`: the cell as written, no reduction. Each text is read twice:
with pymatgen's defaults (`frac`; `frac_tolerance` 1e-4 rounds a coordinate near 1/3, 2/3, 1/6, ... to the ideal value) and
with `frac_tolerance=0` (`fracAsWritten`, the same expansion without that rounding). The kit is compared with the second at
1e-7 and with the first within the 1e-4 of pymatgen's own rounding.
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
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        written = CifParser.from_str(text, frac_tolerance=0).parse_structures(primitive=False)[0]
    assert len(written) == len(s), path.name
    sites = []
    for site, unrounded in zip(s, written):
        (element, occupancy), = site.species.items()
        assert unrounded.species == site.species, path.name
        frac = [float(v) for v in site.frac_coords]
        as_written = [float(v) for v in unrounded.frac_coords]
        assert max(abs(a - b) for a, b in zip(frac, as_written)) <= 1e-4, path.name
        sites.append(dict(z=element.Z, frac=frac, fracAsWritten=as_written, occupancy=float(occupancy)))
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
