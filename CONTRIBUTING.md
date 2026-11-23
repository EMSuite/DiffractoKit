# Contributing to DiffractoKit

Thank you for your interest. Bug reports, questions and pull requests are welcome.

## Build and test

You need a Swift 6.2 toolchain on macOS 15 or later.

```sh
swift test -c release
```

The package builds optimised code in release mode; a Debug build of the reflection loops and the ring fit is slow. The tests
read committed fixtures and need no Python.

## How numbers are checked

Every number a scientist reads is covered by a test against a reference: pymatgen and spglib for cells, symmetry and
d-spacings, py4DSTEM for scattering factors and diffraction geometry, SciPy and NumPy for numerical results, or a closed form
or planted truth where no package applies. The reference outputs are committed as JSON fixtures under `Tests/*/Fixtures/`,
with the generator scripts under `tools/pins/`. `docs/parity.md` lists what is held against which reference.

A test's tolerance is the one the Python package's own tests use for that function; a result outside it is a bug. A change
that moves a fixture says why in the pull request and in the CHANGELOG; regenerate a fixture with its script and the pinned
versions of the reference.

## Regenerating fixtures

Every fixture generator lives under `tools/pins/<subject>/`, and each script's header names the command that runs it. They
run in three Python 3.12 virtual environments outside the repository, beside the checkout under `../References/python/`,
one per oracle, each built from its requirements file:

| Venv | Requirements | Runs |
|---|---|---|
| `../References/python/pymatgen` | `tools/pins/requirements/pymatgen.txt` (pymatgen-core 2026.10.16, spglib 2.8.0, NumPy 2.5.4, SciPy 1.18.1) | every generator except the two below |
| `../References/python/orix` | `tools/pins/requirements/orix.txt` (orix 0.15.0) | `tools/pins/orix` |
| `../References/python/py4dstem` | `tools/pins/requirements/py4dstem.txt` (py4DSTEM f050d207) | `tools/pins/acom-matcher`, and the py4DSTEM halves the pymatgen scripts of `phase-vectors` and `orientation-symmetry` start themselves |

```sh
for oracle in pymatgen orix py4dstem; do
  python3.12 -m venv ../References/python/$oracle
  ../References/python/$oracle/bin/pip install -r tools/pins/requirements/$oracle.txt
done
tools/pins/structure-factors/fetch-table.sh                  # py4DSTEM's scattering table, sha256-locked
python3 tools/crystal-structures/make_beta_double_prime.py   # the beta'' CIF the cell pins read
```

The last two write into the checkout's own gitignored `References/`.

`tools/pins/regenerate.sh --setup` does all of this (and the two pinned source clones some generators check against);
`tools/pins/regenerate.sh` then runs every generator in order and `tools/pins/compare.py` judges what moved: equal, float
noise, or DIFFERS. The Linux reference job (`.github/workflows/references.yml`, weekly and on every PR that touches the
generators or the fixtures) runs exactly these on GitHub's Ubuntu runner and fails on any DIFFERS. The committed fixtures
were written on a Mac, and Linux reproduces them to float noise (`|a - b| <= 1e-12 + 1e-9 max(|a|, |b|)`); what does not
reproduce yet is listed, with its reason, in `tools/pins/known-differences.txt`. `orientation-symmetry/pins.py` and
`acom-matcher/pins.py` read the 4DSTEM Studio source and stay Mac-only. Put the committed files back after a check with
`git checkout -- Tests Sources`.

Every fixture that records a pymatgen version records pymatgen-core 2026.10.16 (the `pymatgen` metapackage is 2026.9.24).
A regenerated fixture that moves is its own commit with the reason and a CHANGELOG line.

## Pull requests

- Keep a change focused, with a test that fails before the change and passes after it.
- Tests use Swift Testing (`import Testing`, `@Test`, `#expect`).
- Public API only: new types and functions are `public`, and the separate package in `tools/consumer-check` must still build.
- Add a line to `CHANGELOG.md` for anything a user can see.
- By contributing you agree that your contribution is licensed under GPL-3.0-only, like the rest of the library.

## Reporting issues

Open an issue with the DiffractoKit version, your platform and Swift version, a small CIF or code sample that reproduces the
problem, what you expected and what you got. For security problems see `SECURITY.md`.
