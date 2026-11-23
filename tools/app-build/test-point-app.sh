#!/bin/bash
# Tests of point-app.sh on planted manifests (run by .github/workflows/app-build.yml on every PR; needs no Swift).
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
fails=0
check() { if eval "$2"; then echo "ok    $1"; else echo "FAIL  $1"; fails=1; fi; }
setup() {
  T=$(mktemp -d); mkdir -p "$T/app" "$T/candidate"; touch "$T/candidate/Package.swift"
  printf '%s\n' 'let package = Package(' '    dependencies: [' \
    '        .package(url: "https://github.com/EMSuite/ImagoKit", from: "0.11.0"),' "$1" '    ],' ')' > "$T/app/Package.swift"
}

setup '        .package(url: "https://github.com/mac4DSTEM/DiffractoKit", from: "0.11.0"),'
"$here/point-app.sh" "$T/app" "$T/candidate" "DiffractoKit" > /dev/null
check "old owner: path dependency under the URL's identity" 'grep -q "\.package(path: \"\.\./DiffractoKit\")," "$T/app/Package.swift"'
check "candidate moved beside the app" '[ -f "$T/DiffractoKit/Package.swift" ] && [ ! -e "$T/candidate" ]'
check "other dependencies untouched" 'grep -q "github.com/EMSuite/ImagoKit\", from: \"0.11.0\"" "$T/app/Package.swift"'
rm -rf "$T"

setup '        .package(url: "https://github.com/EMSuite/DiffractoKit", from: "0.12.0"),'
"$here/point-app.sh" "$T/app" "$T/candidate" "DiffractoKit" > /dev/null
check "new owner" 'grep -q "\.package(path: \"\.\./DiffractoKit\")," "$T/app/Package.swift" && [ -d "$T/DiffractoKit" ]'
rm -rf "$T"

setup '        .package(url: "https://github.com/EMSuite/SpectroKit", from: "0.16.0"),'
check "refuses a manifest without the package" '! "$here/point-app.sh" "$T/app" "$T/candidate" "DiffractoKit" > /dev/null'
rm -rf "$T"
exit $fails
