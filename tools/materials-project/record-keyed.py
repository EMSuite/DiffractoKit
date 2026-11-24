"""DiffractoKit: record the keyed Materials Project answers the tests decode (DiffractoKit #151, plan row K3b).

Run once on the owner's Mac, from the DiffractoKit checkout, in Terminal (not in a sandboxed session, which cannot reach
the service):

    python3 -I tools/materials-project/record-keyed.py

Standard library only. The key is the owner's own and never leaves the Mac's keychain except in the X-API-KEY header
of these requests: it is read with `security find-generic-password` (service org.emsuite.materials-project-recorder,
account api-key). The first run, with no such item, asks `security add-generic-password` to prompt for it; the key is
typed into that prompt, never onto a command line, and this script never prints or writes it. Every body is checked
for the key before it is written, and the run stops if one holds it.

It writes the answers byte for byte, and a manifest with each URL (which never holds the key), the date, the SHA-256,
and the rows each search must decode to, to Tests/DiffractoCoreTests/Fixtures/mp-keyed/:

  heartbeat.json                 GET /heartbeat (sent without the key): the database version
  search-<name>.json             GET /materials/summary/?<query>&deprecated=false&batch_id_not_eq=...&_fields=...
  summary-<id>.json              GET /materials/summary/?material_ids=<id>&_fields=..., for each id of the keyless
                                 recording (mp-optimade/manifest.json), so that both routes can be compared
  manifest.json

The requests are spelled as `MaterialsProjectImport.searchRequest` and `request(materialID:apiKey:)` build them.
"""
import datetime
import hashlib
import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

BASE = "https://api.materialsproject.org"
HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURES = os.path.normpath(os.path.join(HERE, "..", "..", "Tests", "DiffractoCoreTests", "Fixtures"))
OUT = os.path.join(FIXTURES, "mp-keyed")
SERVICE = "org.emsuite.materials-project-recorder"
ACCOUNT = "api-key"
SEARCH_FIELDS = "material_id,formula_pretty,nsites,energy_above_hull,is_stable,theoretical,deprecated,symmetry"
SUMMARY_FIELDS = "material_id,formula_pretty,structure,symmetry,deprecated,last_updated"
EXCLUDED_BATCH = "gnome_r2scan_statics"
RECORDING = 1
# name -> the query part of MaterialsProjectImport.searchRequest, in its order.
SEARCHES = {
    "chemsys-Mg-Si": [("chemsys", "Mg-Si")],
    "elements-Al": [("elements", "Al")],
    "formula-Al2Cu": [("formula", "Al2Cu")],
    "formula-Si": [("formula", "Si")],
    "none": [("formula", "Kr5Xe3")],
}


def read_key():
    def find():
        result = subprocess.run(["security", "find-generic-password", "-s", SERVICE, "-a", ACCOUNT, "-w"],
                                capture_output=True, text=True)
        return result.stdout.strip() if result.returncode == 0 else ""

    key = find()
    if not key:
        print(f"No key in the keychain yet. Paste your Materials Project API key at the prompt "
              f"(it is stored as {SERVICE} / {ACCOUNT} and not shown).")
        # -w as the last option makes `security` prompt for the value, so it never sits on a command line.
        subprocess.run(["security", "add-generic-password", "-s", SERVICE, "-a", ACCOUNT, "-U", "-w"], check=True)
        key = find()
    if not key:
        sys.exit("The key could not be read from the keychain.")
    return key


def url_of(path, items):
    query = urllib.parse.urlencode(items, safe=",")
    return f"{BASE}/{path}" + (f"?{query}" if query else "")


def get(url, key):
    headers = {"User-Agent": "DiffractoKit-fixture-recorder", "Accept": "application/json"}
    if key:
        headers["X-API-KEY"] = key
    request = urllib.request.Request(url, headers=headers)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return response.read()
        except urllib.error.HTTPError as error:
            if key and error.code in (401, 403):
                sys.exit(f"The Materials Project did not accept the key (HTTP {error.code}). Delete the keychain item "
                         f"(security delete-generic-password -s {SERVICE} -a {ACCOUNT}) and run again.")
            if attempt == 3:
                raise
        except Exception:  # noqa: BLE001 - a flaky network gets three retries
            if attempt == 3:
                raise
        print(f"retry {url}", file=sys.stderr)
        time.sleep(5 * (attempt + 1))


def expected_rows(body):
    """The rows MaterialsProjectImport.decodeSearch gives: the most stable first, no energy last, then the id."""
    answer = json.loads(body)
    rows = []
    for item in answer.get("data") or []:
        symmetry = item.get("symmetry")
        rows.append({
            "materialID": item["material_id"], "formula": item.get("formula_pretty"), "siteCount": item.get("nsites"),
            "energyAboveHull": item.get("energy_above_hull"), "isStable": item.get("is_stable"),
            "theoretical": item.get("theoretical"), "deprecated": item.get("deprecated"),
            "symmetry": None if symmetry is None else {
                "crystalSystem": symmetry.get("crystal_system"), "symbol": symmetry.get("symbol"),
                "number": symmetry.get("number"), "pointGroup": symmetry.get("point_group"),
                "hall": symmetry.get("hall")},
        })
    rows.sort(key=lambda r: (float("inf") if r["energyAboveHull"] is None else r["energyAboveHull"], r["materialID"]))
    total = (answer.get("meta") or {}).get("total_doc")
    return rows, (total if total is not None else len(rows)) > len(rows)


def main():
    key = read_key()
    os.makedirs(OUT, exist_ok=True)
    for name in os.listdir(OUT):
        if name.endswith(".json"):
            os.remove(os.path.join(OUT, name))
    manifest = {"recorder": "tools/materials-project/record-keyed.py", "date": datetime.date.today().isoformat(),
                "recording": RECORDING, "files": [], "searches": [], "summaries": []}

    def keep(name, url, body, note):
        if key.encode() in body or key in url:
            sys.exit(f"{name} holds the key; nothing more is written. Delete {OUT} before committing.")
        with open(os.path.join(OUT, name), "wb") as f:
            f.write(body)
        manifest["files"].append({"file": name, "url": url, "sha256": hashlib.sha256(body).hexdigest(),
                                  "bytes": len(body), "note": note})

    url = url_of("heartbeat", [])
    body = get(url, None)
    keep("heartbeat.json", url, body, "the service's state and database version")
    manifest["databaseVersion"] = json.loads(body)["db_version"]

    for name, query in SEARCHES.items():
        url = url_of("materials/summary/", query + [
            ("deprecated", "false"), ("batch_id_not_eq", EXCLUDED_BATCH), ("_fields", SEARCH_FIELDS),
            ("_sort_fields", "energy_above_hull"), ("_limit", "500")])
        body = get(url, key)
        keep(f"search-{name}.json", url, body, f"the search {name}")
        rows, more = expected_rows(body)
        manifest["searches"].append({"file": f"search-{name}.json", "more": more, "rows": rows})
        time.sleep(0.5)

    with open(os.path.join(FIXTURES, "mp-optimade", "manifest.json")) as f:
        ids = sorted(set(json.load(f)["lattices"].values()))
    for material_id in ids:
        url = url_of("materials/summary/", [("material_ids", material_id), ("_fields", SUMMARY_FIELDS)])
        keep(f"summary-{material_id}.json", url, get(url, key), f"the structure of {material_id}")
        manifest["summaries"].append(material_id)
        time.sleep(0.5)

    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1, sort_keys=True)
        f.write("\n")
    print(f"Recorded {len(manifest['files'])} answers (database {manifest['databaseVersion']}) in {OUT}.")


if __name__ == "__main__":
    main()
