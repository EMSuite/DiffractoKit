"""DiffractoKit: record the keyless Materials Project responses the tests decode (DiffractoKit #150, plan row K3a).

Run by .github/workflows/record-materials-project.yml on GitHub's Linux runner, which reaches the service (the cloud
sessions and the Mac's sandbox do not). It needs Python 3.12, numpy and spglib 2.8.0 (tools/pins/requirements/pymatgen.txt
pins both); it sends no key and nothing about anyone.

It writes the response bodies byte for byte as the server sent them, and a manifest with each URL, the date and the
SHA-256 of the bytes, to Tests/DiffractoCoreTests/Fixtures/mp-optimade/:

  info.json                      GET /v1/info
  search-<name>.json             GET /v1/structures?filter=...&response_fields=chemical_formula_reduced,nsites,_mp_stability
  entry-<id>.json                GET /v1/structures/<id>, one entry for each Bravais lattice it could fill
  manifest.json

The entries are chosen, not typed: for each candidate formula the search is recorded, its rows are sorted by the lowest
energy above the hull over the functionals listed, the five lowest are fetched, and spglib names each one's space group
at symprec 1e-4, 1e-3 and 1e-2 Å (a hundred to ten thousand times the server's rounding of 1e-6 Å). An entry is kept only
when all three agree: a calculated cell that is nearly of a higher symmetry (a monoclinic angle 0.03 degrees off 90, a
bcc cell strained by 0.2 %) reads differently at each tolerance, and no reference can pin it. The first stable entry of
each of the 14 Bravais lattices is kept; lattices no candidate fills are listed in the manifest as missing.
"""
import datetime
import hashlib
import json
import os
import sys
import time
import urllib.parse
import urllib.request

import numpy as np
import spglib

BASE = "https://optimade.materialsproject.org/v1"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "..", "Tests", "DiffractoCoreTests", "Fixtures", "mp-optimade"))
FIELDS = "chemical_formula_reduced,nsites,_mp_stability"
# OPTIMADE reduced formulas list the elements alphabetically.
CANDIDATES = ["Al", "Fe", "Mg", "Si", "ClCs", "Po", "In", "Sn", "O2Ti", "Al2O3", "Bi", "Ga", "U", "S", "I", "CuO",
              "O2Zr", "CFe3", "Al2Cu", "AlB2", "B", "Se", "Te", "P", "Mn", "Li2MnO3", "Ca", "Na", "Hg", "Pa", "CaCO3",
              "K2O4S", "Al2O5Si", "MgO3Si", "Cu", "Zn", "Sb", "As", "Ge"]
# The version of what is chosen: 1 took the first entry of each lattice, 2 only entries whose space group does not
# depend on the tolerance. Raise it when the selection changes.
RECORDING = 2
# The tolerances an entry's space group must not depend on (Å).
SYMPRECS = [1e-4, 1e-3, 1e-2]
LATTICES = ["aP", "mP", "mC", "oP", "oC", "oI", "oF", "tP", "tI", "hR", "hP", "cP", "cI", "cF"]


def get(url):
    request = urllib.request.Request(url, headers={"User-Agent": "DiffractoKit-fixture-recorder"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return response.read()
        except Exception as error:  # noqa: BLE001 - a flaky network gets three retries
            if attempt == 3:
                raise
            print(f"retry {url}: {error}", file=sys.stderr)
            time.sleep(5 * (attempt + 1))


def search_url(formula):
    query = urllib.parse.urlencode({"filter": f'chemical_formula_reduced="{formula}"', "response_fields": FIELDS,
                                    "page_limit": 500})
    return f"{BASE}/structures?{query}"


def bravais(attributes):
    lattice = np.array(attributes["lattice_vectors"], dtype=float)
    cartesian = np.array(attributes["cartesian_site_positions"], dtype=float)
    fractional = cartesian @ np.linalg.inv(lattice)
    names = attributes["species_at_sites"]
    kinds = {name: i + 1 for i, name in enumerate(sorted(set(names)))}
    cell = (lattice, fractional, [kinds[n] for n in names])
    datasets = [spglib.get_symmetry_dataset(cell, symprec=symprec) for symprec in SYMPRECS]
    if any(dataset is None for dataset in datasets) or len({dataset.number for dataset in datasets}) != 1:
        return None
    symbol = datasets[0].international
    number = datasets[0].number
    centring = symbol[0]
    if number <= 2:
        family = "a"
    elif number <= 15:
        family = "m"
    elif number <= 74:
        family = "o"
    elif number <= 142:
        family = "t"
    elif number <= 194:
        family = "h"
    else:
        family = "c"
    if centring in "ABC":
        centring = "C"
    if family == "h" and centring == "R":
        return "hR"
    return family + centring


def main():
    # The workflow runs on every push to a pull request that holds this file. A recording made with this selection
    # stands (the service's answers differ from day to day only in their metadata), unless --force asks for a new one;
    # a change to what is chosen raises RECORDING, which records again on the next push.
    manifest_path = os.path.join(OUT, "manifest.json")
    if "--force" not in sys.argv and os.path.exists(manifest_path):
        with open(manifest_path) as f:
            if json.load(f).get("recording") == RECORDING:
                print(f"the committed recording is recording {RECORDING}; nothing to do (--force records again)")
                return
    os.makedirs(OUT, exist_ok=True)
    # A new recording replaces the old one whole (the pins file beside it is tools/pins/mp-optimade/pins.py's).
    for name in os.listdir(OUT):
        if name.endswith(".json") and name != "pymatgen-pins.json":
            os.remove(os.path.join(OUT, name))
    manifest = {"recorder": "tools/materials-project/record.py", "date": datetime.date.today().isoformat(),
                "recording": RECORDING, "spglib": spglib.__version__, "symprecs": SYMPRECS, "files": [], "lattices": {}}

    def keep(name, url, body, note):
        with open(os.path.join(OUT, name), "wb") as f:
            f.write(body)
        manifest["files"].append({"file": name, "url": url, "sha256": hashlib.sha256(body).hexdigest(),
                                  "bytes": len(body), "note": note})

    url = f"{BASE}/info"
    keep("info.json", url, get(url), "the service's info")
    url = search_url("Xe3Kr5")
    keep("search-none.json", url, get(url), "a formula with no entry: the service answers, with no rows")

    for formula in CANDIDATES:
        if len(manifest["lattices"]) == len(LATTICES):
            break
        url = search_url(formula)
        body = get(url)
        rows = json.loads(body)["data"]

        def lowest(row):
            stability = (row.get("attributes") or {}).get("_mp_stability") or {}
            values = [v.get("energy_above_hull") for v in stability.values() if isinstance(v, dict)]
            values = [v for v in values if isinstance(v, (int, float))]
            return min(values) if values else float("inf")

        kept_search = False
        for row in sorted(rows, key=lambda r: (lowest(r), r["id"]))[:5]:
            entry_url = f"{BASE}/structures/{urllib.parse.quote(row['id'])}"
            entry = get(entry_url)
            lattice = bravais(json.loads(entry)["data"]["attributes"])
            time.sleep(0.5)
            if lattice is None or lattice in manifest["lattices"]:
                continue
            if not kept_search:
                keep(f"search-{formula}.json", url, body, f"the search for {formula}")
                kept_search = True
            keep(f"entry-{row['id']}.json", entry_url, entry, f"{formula}, {lattice} (spglib, the same at symprec 1e-4, 1e-3 and 1e-2)")
            manifest["lattices"][lattice] = row["id"]
        time.sleep(0.5)

    manifest["missing"] = [name for name in LATTICES if name not in manifest["lattices"]]
    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1, sort_keys=True)
        f.write("\n")
    print(json.dumps({"lattices": manifest["lattices"], "missing": manifest["missing"]}, indent=1))


if __name__ == "__main__":
    main()
