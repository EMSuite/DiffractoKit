//
//  StructureSnapshotTests.swift
//  Role: the Codable structure snapshot and the content fingerprint (DiffractoKit #148, plan row K1). Every committed
//        structure survives encode and decode bit for bit and rebuilds to an equal structure with the same fingerprint;
//        the fingerprint is equal for equal content and changes with any cell or site value (a position by a step of its
//        1e-6 grid; one ulp of a position, as a hexagonal CIF read back gives it, is the same content, #220); the schema
//        version is written and a newer one is refused.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import DiffractoCore

struct StructureSnapshotTests {

    typealias Snapshot = CrystalStructure.Snapshot

    /// Every committed CIF that reads, by file name, sorted.
    static func committedStructures() throws -> [(name: String, structure: CrystalStructure)] {
        let root = FixtureDirs.tests
        guard let walker = FileManager.default.enumerator(at: root, includingPropertiesForKeys: nil) else { return [] }
        var files: [URL] = []
        for case let url as URL in walker where url.pathExtension.lowercased() == "cif" { files.append(url) }
        var loaded: [(name: String, structure: CrystalStructure)] = []
        var refused: [String] = []
        for url in files.sorted(by: { $0.path < $1.path }) {
            let bytes = try Data(contentsOf: url)
            if let crystal = try? CrystalLoading.load(bytes: bytes, fileName: url.lastPathComponent) {
                loaded.append((url.lastPathComponent, crystal.structure))
            } else {
                refused.append(url.lastPathComponent)
            }
        }
        // The edge cases refused by design (three_structures.cif names no block of its three); any other CIF that stops
        // reading fails here.
        #expect(refused == ["macromolecular.cif", "magnetic.cif", "nan_coordinate.cif", "three_structures.cif"], "\(refused)")
        return loaded
    }

    static let presets: [(name: String, structure: CrystalStructure)] = [
        ("Al", CrystalPresets.aluminum), ("Au", CrystalPresets.gold), ("Ni", CrystalPresets.nickel),
        ("Cu", CrystalPresets.copper), ("Fe", CrystalPresets.iron), ("Si", CrystalPresets.silicon),
        ("Mg", CrystalPresets.magnesium), ("beta''", CrystalPresets.betaDoublePrime),
        ("WS2", CrystalPresets.tungstenDisulfide),
    ]

    static func encoder() -> JSONEncoder {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys]
        return encoder
    }

    /// The bits of a value with a negative zero read as zero (the snapshot stores zero for it).
    static func bits(_ value: Double) -> UInt64 { (value == 0 ? 0 : value).bitPattern }

    static func expectSameBits(_ lhs: CrystalStructure, _ rhs: CrystalStructure, _ name: String) {
        let l = lhs.cell, r = rhs.cell
        for (x, y) in zip([l.a, l.b, l.c, l.alphaDeg, l.betaDeg, l.gammaDeg], [r.a, r.b, r.c, r.alphaDeg, r.betaDeg, r.gammaDeg]) {
            #expect(bits(x) == bits(y), "\(name): cell \(x) vs \(y)")
        }
        #expect(lhs.sites.count == rhs.sites.count, "\(name): site count")
        for (s, t) in zip(lhs.sites, rhs.sites) {
            #expect(s.z == t.z && s.label == t.label, "\(name): site")
            for k in 0..<3 { #expect(bits(s.fractional[k]) == bits(t.fractional[k]), "\(name): coordinate") }
            #expect(bits(s.occupancy) == bits(t.occupancy), "\(name): occupancy")
            #expect(s.isotropicDisplacement.map(bits) == t.isotropicDisplacement.map(bits), "\(name): displacement")
        }
    }

    /// The round trip: snapshot, JSON, snapshot; the same bytes again; the rebuilt structure equal, bit for bit, with
    /// the same fingerprint.
    static func expectRoundTrip(_ structure: CrystalStructure, _ name: String) throws {
        let snapshot = Snapshot(structure)
        let data = try encoder().encode(snapshot)
        let decoded = try JSONDecoder().decode(Snapshot.self, from: data)
        #expect(decoded == snapshot, "\(name): decoded snapshot")
        #expect(try encoder().encode(decoded) == data, "\(name): encoded again")
        let rebuilt = decoded.structure()
        #expect(rebuilt == structure, "\(name): rebuilt structure")
        expectSameBits(rebuilt, structure, name)
        #expect(rebuilt.fingerprint == structure.fingerprint, "\(name): fingerprint of the rebuilt structure")
        #expect(decoded.fingerprint == structure.fingerprint, "\(name): stored fingerprint")
    }

    @Test func everyCommittedStructureSurvivesEncodeAndDecodeBitForBit() throws {
        let committed = try Self.committedStructures()
        // 64 committed CIFs; the four edge cases refused by design (mmCIF, magnetic, a NaN coordinate, three structures
        // with no block named) do not read. 28 of them are the kit's own CIFExport output (cif-export/written/).
        #expect(committed.count == 60, "\(committed.count) committed structures read")
        for (name, structure) in committed + Self.presets {
            try Self.expectRoundTrip(structure, name)
        }
    }

    @Test func theConventionalSettingAndGivenOperationsRoundTripToo() throws {
        for (name, structure) in try Self.committedStructures() + Self.presets {
            if let conventional = structure.conventionalSetting?.structure {
                #expect(conventional.setting == .conventional)
                try Self.expectRoundTrip(conventional, "\(name) conventional")
            }
        }
        let given = CrystalStructure(cell: UnitCell.cubic(a: 4.0782),
                                     sites: [AtomSite(z: 79, fractional: [0, 0, 0], label: "Au1")],
                                     operations: [.identity], spaceGroupNumber: 1, name: "gold, P1 as given")
        try Self.expectRoundTrip(given, "given")
        #expect(given.detectionSettings == nil)
    }

    /// The detection's settings travel with the snapshot: a CIF read at a Cartesian tolerance and a looser metric
    /// rebuilds the same structure, not the default read.
    @Test func theDetectionSettingsTravelWithTheSnapshot() throws {
        var policy = CIFImport.Policy.default
        policy.positionTolerance = .cartesian(1e-3)
        policy.metricLengthTolerance = 1e-5
        policy.symmetry.gridSnapFactor = 2
        let doc = try CIFImport.read(try Fixtures.text("alpha_quartz_P3221.cif"), policy: policy)
        #expect(doc.structure.detectionSettings?.tolerance == .cartesian(1e-3))
        #expect(doc.structure.detectionSettings?.metricLengthTolerance == 1e-5)
        #expect(doc.structure.detectionSettings?.policy.gridSnapFactor == 2)
        #expect((doc.structure.detectionSettings?.coarsestHalfStep ?? 0) > 0)
        try Self.expectRoundTrip(doc.structure, "quartz, own policy")
        let unconventional = CrystalStructure.detectingSymmetry(
            cell: CrystalPresets.magnesium.cell, sites: CrystalPresets.magnesium.sites, conventional: false)
        try Self.expectRoundTrip(unconventional, "Mg, file cell only")
    }

    @Test func theSchemaVersionIsWrittenAndANewerOneIsRefused() throws {
        let data = try Self.encoder().encode(Snapshot(CrystalPresets.gold))
        let object = try #require(try JSONSerialization.jsonObject(with: data) as? [String: Any])
        #expect(object["schemaVersion"] as? Int == Snapshot.currentSchemaVersion)
        #expect(Snapshot.currentSchemaVersion == 1)
        let symmetry = try #require(object["symmetry"] as? [String: Any])
        #expect(symmetry["source"] as? String == "detected")

        func decode(_ edit: (inout [String: Any]) -> Void) -> Error? {
            var copy = object
            edit(&copy)
            do {
                let edited = try JSONSerialization.data(withJSONObject: copy)
                _ = try JSONDecoder().decode(Snapshot.self, from: edited)
                return nil
            } catch {
                return error
            }
        }
        #expect(decode { _ in } == nil)
        let newer = decode { $0["schemaVersion"] = 2 }
        #expect(newer is DecodingError)
        #expect("\(String(describing: newer))".contains("schema 2"))
        #expect(decode { $0["schemaVersion"] = 0 } is DecodingError)
        #expect(decode { $0["symmetry"] = ["source": "guessed"] } is DecodingError)
        #expect(decode {
            var settings = symmetry["settings"] as? [String: Any] ?? [:]
            settings["toleranceKind"] = "spherical"
            var edited = symmetry
            edited["settings"] = settings
            $0["symmetry"] = edited
        } is DecodingError)
        #expect(decode {
            $0["symmetry"] = ["source": "given", "operations": [["rotation": [1, 0, 0, 0, 1, 0, 0, 0], "translation": [0, 0, 0]] as [String: Any]]] as [String: Any]
        } is DecodingError)
    }

    @Test func theFingerprintIsEqualForEqualContent() throws {
        let gold = CrystalPresets.gold
        #expect(gold.fingerprint == CrystalPresets.fcc(a: 4.0782, z: 79).fingerprint)
        #expect(gold.fingerprint.count == 64)
        #expect(gold.fingerprint.allSatisfy { $0.isHexDigit })

        // The order of the sites and of the operations, the name, the labels and the stated number are not content.
        let sites = gold.sites
        let relabelled = sites.enumerated().map { i, s in
            AtomSite(z: s.z, fractional: s.fractional, occupancy: s.occupancy, label: "Au\(i)")
        }
        let a = CrystalStructure(cell: gold.cell, sites: sites, operations: gold.operations)
        let b = CrystalStructure(cell: gold.cell, sites: relabelled.reversed(), operations: gold.operations.reversed(),
                                 spaceGroupNumber: 225, name: "gold")
        #expect(a.fingerprint == b.fingerprint)
        #expect(a.fingerprint == gold.fingerprint)

        // Two reads that find the same operations: the settings they were given are not content.
        let cartesian = CrystalStructure.detectingSymmetry(cell: gold.cell, sites: sites, tolerance: .cartesian(1e-3))
        #expect(cartesian.operations.count == gold.operations.count)
        #expect(cartesian.fingerprint == gold.fingerprint)

        // A negative zero counts as zero.
        let negativeZero = CrystalStructure(cell: gold.cell, sites: [AtomSite(z: 79, fractional: [-0.0, 0, 0])],
                                            operations: [.identity])
        let positiveZero = CrystalStructure(cell: gold.cell, sites: [AtomSite(z: 79, fractional: [0, 0, 0])],
                                            operations: [.identity])
        #expect(negativeZero.fingerprint == positiveZero.fingerprint)

        // A position one ulp off, as a CIF's operations put the image of a third (1 − 1/3 against 2/3), and a position
        // across the cell edge (1 − 1e-12 against 0), are the same content (#220).
        func one(_ x: Double) -> CrystalStructure {
            CrystalStructure(cell: gold.cell, sites: [AtomSite(z: 79, fractional: [x, 0.25, 0])], operations: [.identity])
        }
        #expect((1 - 1.0 / 3) != 2.0 / 3)
        #expect(one(1 - 1.0 / 3).fingerprint == one(2.0 / 3).fingerprint)
        #expect(one(0.6666667).fingerprint == one(0.6666667.nextUp).fingerprint)
        #expect(one(1 - 1e-12).fingerprint == one(0).fingerprint)
        let stored = try JSONDecoder().decode(Snapshot.self, from: Self.encoder().encode(Snapshot(negativeZero)))
        #expect(stored.sites[0].fractional.x.sign == .plus)
    }

    @Test func theFingerprintChangesWithAnyCellOrSiteValue() {
        let cell = UnitCell(a: 5.1, b: 6.2, c: 7.3, alphaDeg: 88, betaDeg: 97, gammaDeg: 101)
        let site = AtomSite(z: 14, fractional: [0.1, 0.2, 0.3], occupancy: 0.9, isotropicDisplacement: 0.01)
        let other = AtomSite(z: 8, fractional: [0.6, 0.7, 0.8])
        func structure(_ cell: UnitCell, _ first: AtomSite, operations: [SymmetryOperation] = [.identity]) -> CrystalStructure {
            CrystalStructure(cell: cell, sites: [first, other], operations: operations)
        }
        let base = structure(cell, site)
        var variants: [(String, CrystalStructure)] = []
        let p = [cell.a, cell.b, cell.c, cell.alphaDeg, cell.betaDeg, cell.gammaDeg]
        for k in 0..<6 {
            var q = p
            q[k] = q[k].nextUp
            variants.append(("cell \(k)", structure(UnitCell(a: q[0], b: q[1], c: q[2], alphaDeg: q[3], betaDeg: q[4],
                                                              gammaDeg: q[5]), site)))
        }
        for k in 0..<3 {
            var moved = site
            moved.fractional[k] += 1e-6
            variants.append(("coordinate \(k)", structure(cell, moved)))
        }
        var changed = site
        changed.occupancy = site.occupancy.nextUp
        variants.append(("occupancy", structure(cell, changed)))
        changed = site
        changed.z = 15
        variants.append(("element", structure(cell, changed)))
        changed = site
        changed.isotropicDisplacement = 0.01.nextUp
        variants.append(("displacement", structure(cell, changed)))
        changed = site
        changed.isotropicDisplacement = nil
        variants.append(("no displacement", structure(cell, changed)))
        variants.append(("one site fewer", CrystalStructure(cell: cell, sites: [site], operations: [.identity])))
        variants.append(("an inversion", structure(cell, site, operations: [.identity, SymmetryOperation(rotation: .inversion)])))

        var seen: [String: String] = [base.fingerprint: "base"]
        for (name, variant) in variants {
            #expect(seen[variant.fingerprint] == nil, "\(name) has the fingerprint of \(seen[variant.fingerprint] ?? "")")
            seen[variant.fingerprint] = name
        }
        #expect(seen.count == variants.count + 1)
    }

    /// Settings the detection reads as 0 are stored as 0: a NaN metric tolerance or policy value gives a structure equal
    /// to the one read at 0, equal to itself, and written.
    @Test func settingsAreStoredAsTheDetectionReadsThem() throws {
        let mg = CrystalPresets.magnesium
        var nanPolicy = SymmetryDetection.Policy.default
        nanPolicy.gridSnapFactor = .nan
        var zeroPolicy = SymmetryDetection.Policy.default
        zeroPolicy.gridSnapFactor = 0
        let nan = CrystalStructure.detectingSymmetry(cell: mg.cell, sites: mg.sites, metricLengthTolerance: .nan,
                                                     metricAngleToleranceDeg: -1, policy: nanPolicy)
        let zero = CrystalStructure.detectingSymmetry(cell: mg.cell, sites: mg.sites, metricLengthTolerance: 0,
                                                      metricAngleToleranceDeg: 0, policy: zeroPolicy)
        #expect(nan == zero)
        #expect(nan == nan)
        try Self.expectRoundTrip(nan, "Mg, NaN settings")
        let refused = CrystalStructure.detectingSymmetry(cell: mg.cell, sites: mg.sites, tolerance: .fractional(.nan))
        #expect(refused.symmetryRefusal != nil)
        #expect(refused == refused)
    }

    /// A snapshot of a structure with a non-finite value is not written: JSON has no NaN.
    @Test func aNonFiniteValueIsNotWritten() {
        let structure = CrystalStructure(cell: UnitCell.cubic(a: 4), sites: [AtomSite(z: 13, fractional: [0, 0, 0], occupancy: .nan)],
                                         operations: [.identity])
        #expect(!structure.fingerprint.isEmpty)
        #expect(throws: EncodingError.self) { try Self.encoder().encode(Snapshot(structure)) }
    }
}
