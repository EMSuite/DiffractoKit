"""DiffractoKit: record what a keyless Materials Project search can say about each row's symmetry (DiffractoKit #263).

Run by .github/workflows/record-materials-project.yml on GitHub's Linux runner, which reaches the service (the cloud
sessions and the Mac's sandbox do not). It needs numpy and spglib 2.8.0, sends no key and nothing about anyone.

It writes the response bodies byte for byte, and a manifest with each URL, the date and the SHA-256 of the bytes, to
Tests/DiffractoCoreTests/Fixtures/mp-optimade-search/:

  info-structures.json      GET /v1/info/structures: every property the service lists for a structure
  search-GaAs-all.json      the GaAs search asking for every listed property (which ones name a symmetry, filled?)
  search-GaAs-sites.json    the GaAs search asking for the row fields and the structure fields (lattice, sites, species)
  manifest.json             also spglib's space group of each GaAs row at symprec 1e-4, 1e-3 and 1e-2 Å

The second search is the one the kit sends when no field names the symmetry; spglib's numbers are the reference its
own detection on each row is checked against.
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
OUT = os.path.normpath(os.path.join(HERE, "..", "..", "Tests", "DiffractoCoreTests", "Fixtures", "mp-optimade-search"))
RECORDING = 1
ROW_FIELDS = "chemical_formula_reduced,nsites,_mp_stability"
SITE_FIELDS = ROW_FIELDS + ",lattice_vectors,cartesian_site_positions,species_at_sites,species"
SYMPRECS = [1e-4, 1e-3, 1e-2]


def get(url):
    request = urllib.request.Request(url, headers={"User-Agent": "DiffractoKit-fixture-recorder"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read()
        except Exception as error:  # noqa: BLE001 - a flaky network gets three retries
            if attempt == 3:
                raise
            print(f"retry {url}: {error}", file=sys.stderr)
            time.sleep(5 * (attempt + 1))


def search_url(formula, fields):
    query = urllib.parse.urlencode({"filter": f'chemical_formula_reduced="{formula}"', "response_fields": fields,
                                    "page_limit": 500})
    return f"{BASE}/structures?{query}"


def space_groups(attributes):
    lattice = np.array(attributes["lattice_vectors"], dtype=float)
    fractional = np.array(attributes["cartesian_site_positions"], dtype=float) @ np.linalg.inv(lattice)
    names = attributes["species_at_sites"]
    kinds = {name: i + 1 for i, name in enumerate(sorted(set(names)))}
    out = []
    for symprec in SYMPRECS:
        dataset = spglib.get_symmetry_dataset((lattice, fractional, [kinds[n] for n in names]), symprec=symprec)
        out.append(None if dataset is None else {"number": dataset.number, "symbol": dataset.international})
    return out


def main():
    manifest_path = os.path.join(OUT, "manifest.json")
    if "--force" not in sys.argv and os.path.exists(manifest_path):
        with open(manifest_path) as f:
            if json.load(f).get("recording") == RECORDING:
                print(f"the committed recording is recording {RECORDING}; nothing to do (--force records again)")
                return
    os.makedirs(OUT, exist_ok=True)
    manifest = {"recorder": "tools/materials-project/record-search-symmetry.py",
                "date": datetime.date.today().isoformat(), "recording": RECORDING, "spglib": spglib.__version__,
                "symprecs": SYMPRECS, "files": [], "spglib_space_groups": {}}

    def keep(name, url, body, note):
        with open(os.path.join(OUT, name), "wb") as f:
            f.write(body)
        manifest["files"].append({"file": name, "url": url, "sha256": hashlib.sha256(body).hexdigest(),
                                  "bytes": len(body), "note": note})

    url = f"{BASE}/info/structures"
    body = get(url)
    keep("info-structures.json", url, body, "every property the service lists for a structure")
    properties = sorted(json.loads(body)["data"]["properties"])

    url = search_url("AsGa", ",".join(properties))
    keep("search-GaAs-all.json", url, get(url), "the GaAs search asking for every listed property")

    url = search_url("AsGa", SITE_FIELDS)
    body = get(url)
    keep("search-GaAs-sites.json", url, body, "the GaAs search with the row fields and the structure fields")
    for row in json.loads(body)["data"]:
        manifest["spglib_space_groups"][row["id"]] = space_groups(row["attributes"])

    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=1, sort_keys=True)
        f.write("\n")
    print(json.dumps(manifest["spglib_space_groups"], indent=1))


if __name__ == "__main__":
    main()
