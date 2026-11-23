"""pymatgen's VESTA bonds (CutOffDictNN.from_preset("vesta_2019").get_nn_info) of six structures.
Run with References/python/pymatgen/bin/python -I -B (tools/pins/requirements/pymatgen.txt).
Every site is written as pymatgen holds it (full cell, coordinates as given, not wrapped), so the kit builds the same
cell from the same numbers; each bond is site, neighbour, image (neighbour's coordinate minus the stored one,
pymatgen's jimage) and length. No expected value is computed by the Swift port.
"""
from pathlib import Path
from importlib.metadata import version
import json, platform
import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.analysis.local_env import CutOffDictNN
assert version('pymatgen') == '2026.9.24'
assert version('pymatgen-core') == '2026.10.16'
assert np.__version__ == '2.5.4'

def cubic(a):
    return Lattice.from_parameters(a, a, a, 90, 90, 90, vesta=True)

structures = {
    # Silicon: the committed silicon_Fd-3m.cif cell, origin choice 1.
    'silicon': Structure.from_spacegroup('Fd-3m', cubic(5.4309), ['Si'], [[0, 0, 0]]),
    # The same eight atoms moved by -1/8 along every axis: coordinates below 0, so the images are tested as stored.
    'silicon shifted': None,
    'rocksalt NaCl': Structure.from_spacegroup('Fm-3m', cubic(5.6402), ['Na', 'Cl'], [[0, 0, 0], [0.5, 0.5, 0.5]]),
    'aluminium': Structure.from_spacegroup('Fm-3m', cubic(4.0495), ['Al'], [[0, 0, 0]]),
    'strontium titanate': Structure.from_spacegroup('Pm-3m', cubic(3.905), ['Sr', 'Ti', 'O'],
                                                    [[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0]]),
    # Baddeleyite ZrO2, P2_1/c, close to the room-temperature refinement (Zr 7-coordinate).
    'baddeleyite ZrO2': Structure.from_spacegroup(
        'P2_1/c', Lattice.from_parameters(5.1505, 5.2116, 5.3173, 90, 99.23, 90, vesta=True), ['Zr', 'O', 'O'],
        [[0.2754, 0.0395, 0.2083], [0.0700, 0.3317, 0.3447], [0.4496, 0.7569, 0.4792]]),
}
si = structures['silicon']
structures['silicon shifted'] = Structure(si.lattice, si.species, si.frac_coords - 0.125)

nn = CutOffDictNN.from_preset('vesta_2019')
cases = []
for name, s in structures.items():
    bonds = []
    for i in range(len(s)):
        for info in nn.get_nn_info(s, i):
            j = int(info['site_index'])
            image = [int(v) for v in info['image']]
            # The image is the neighbour's coordinate minus the stored one, exactly.
            assert np.allclose(info['site'].frac_coords, s[j].frac_coords + image, atol=1e-12, rtol=0), (name, i, j)
            bonds.append([i, j, *image, float(info['weight'])])
    bonds.sort(key=lambda b: tuple(b[:5]))
    l = s.lattice
    cases.append(dict(name=name, a=l.a, b=l.b, c=l.c, alpha=l.alpha, beta=l.beta, gamma=l.gamma,
                      sites=[dict(z=site.specie.Z, fractional=site.frac_coords.tolist()) for site in s],
                      bonds=bonds))
out = dict(python=platform.python_version(), pymatgen=version('pymatgen'), pymatgenCore=version('pymatgen-core'),
           numpy=np.__version__, preset='vesta_2019', cases=cases)
path = Path(__file__).resolve().parents[3] / 'Tests/DiffractoCoreTests/Fixtures/vesta-bond-pins.json'
path.write_text(json.dumps(out, indent=1) + '\n')
for c in cases:
    per = {}
    for b in c['bonds']:
        per.setdefault(b[0], 0); per[b[0]] += 1
    zs = [s['z'] for s in c['sites']]
    print(c['name'], len(c['sites']), 'atoms', len(c['bonds']), 'bonds; per site', sorted(set((zs[k], v) for k, v in per.items())))
