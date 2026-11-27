//
//  MaterialsProjectSearchSymmetryTests.swift
//  Role: each row of a keyless Materials Project search names its space group (DiffractoKit #263, the kit half of
//        Electron-Microscopy-Studio#1115), so zincblende and wurtzite GaAs can be told apart before one is chosen.
//        Recorded byte for byte in `Fixtures/mp-optimade-search/` by `tools/materials-project/record-search-symmetry.py`
//        on GitHub's Linux runner, with spglib 2.8.0's space group of every row at symprec 1e-4, 1e-3 and 1e-2 Å.
//  Done when: the service states no symmetry (every `space_group_*` field of every row is null, so asking for them
//  does not help); the search the kit sends asks for the structure; every row read from it names spglib's group (number
//  and short symbol) or, where the kit's conventional cell is chosen by the metric and reads no group
//  (`CrystalStructure.spaceGroupType`, docs/parity.md), at least its crystal system; GaAs mp-2534 is F-43m (216) and
//  mp-8883 is P6_3mc (186).
//

import CryptoKit
import DiffractoKitTestSupport
import Foundation
import Testing
import DiffractoCore

struct MaterialsProjectSearchSymmetryTests {

    static let directory = FixtureDirs.fixtures(of: "DiffractoCoreTests").appendingPathComponent("mp-optimade-search")

    static func data(_ name: String) throws -> Data { try Data(contentsOf: directory.appendingPathComponent(name)) }

    struct Manifest: Decodable {
        struct File: Decodable {
            let file: String
            let url: String
            let sha256: String
            let bytes: Int
        }
        struct Group: Decodable {
            let number: Int
            let symbol: String
        }
        let files: [File]
        let spglib_space_groups: [String: [Group?]]
    }

    static let manifest: Manifest = try! JSONDecoder().decode(Manifest.self, from: data("manifest.json"))

    @Test func theRecordingsAreTheBytesTheManifestNames() throws {
        #expect(Set(Self.manifest.files.map(\.file)) == ["info-structures.json", "search-GaAs-all.json", "search-GaAs-sites.json"])
        for file in Self.manifest.files {
            let body = try Self.data(file.file)
            #expect(body.count == file.bytes, "\(file.file) size")
            let digest = SHA256.hash(data: body).map { String(format: "%02x", $0) }.joined()
            #expect(digest == file.sha256, "\(file.file) SHA-256")
            #expect(file.url.hasPrefix("https://optimade.materialsproject.org/v1/"), "\(file.file) URL")
        }
        // The sites search is the one the kit sends.
        let sites = try #require(Self.manifest.files.first { $0.file == "search-GaAs-sites.json" })
        let request = try MaterialsProjectOptimade.searchRequest(.formula("GaAs"))
        #expect(sites.url == request.url?.absoluteString)
    }

    /// The diagnosis: asked for every property it lists, the service fills none that names a symmetry.
    @Test func theServiceStatesNoSymmetry() throws {
        struct Envelope: Decodable {
            struct Item: Decodable { let attributes: [String: JSONValue] }
            let data: [Item]
        }
        let all = try JSONDecoder().decode(Envelope.self, from: Self.data("search-GaAs-all.json"))
        #expect(all.data.count == Self.manifest.spglib_space_groups.count)
        for item in all.data {
            for field in ["space_group_it_number", "space_group_symbol_hermann_mauguin",
                          "space_group_symbol_hermann_mauguin_extended", "space_group_symbol_hall",
                          "space_group_symmetry_operations_xyz"] {
                #expect(item.attributes[field] == .null, "\(field)")
            }
        }
    }

    @Test func everyRowNamesSpglibsGroupOrItsCrystalSystem() throws {
        guard case .rows(let rows, let more) = try MaterialsProjectOptimade.decodeSearch(Self.data("search-GaAs-sites.json")) else {
            Issue.record("expected rows")
            return
        }
        #expect(!more)
        #expect(Set(rows.map(\.materialID)) == Set(Self.manifest.spglib_space_groups.keys))
        var named = 0
        for row in rows {
            let spglib = try #require(Self.manifest.spglib_space_groups[row.materialID]?.first ?? nil, "\(row.materialID)")
            let symmetry = try #require(row.symmetry, "\(row.materialID) has no symmetry")
            if let number = symmetry.number {
                named += 1
                #expect(number == spglib.number, "\(row.materialID): \(number), spglib \(spglib.number)")
                #expect(symmetry.symbol == spglib.symbol, "\(row.materialID): \(symmetry.symbol ?? "-"), spglib \(spglib.symbol)")
            } else {
                #expect(symmetry.crystalSystem != nil, "\(row.materialID): neither a group nor a crystal system")
            }
        }
        print("GaAs search: \(named) of \(rows.count) rows name spglib's group, the rest their crystal system")

        let zincblende = try #require(rows.first { $0.materialID == "mp-2534" }?.symmetry)
        #expect(zincblende.number == 216)
        #expect(zincblende.symbol == "F-43m")
        #expect(zincblende.crystalSystem == "Cubic")
        let wurtzite = try #require(rows.first { $0.materialID == "mp-8883" }?.symmetry)
        #expect(wurtzite.number == 186)
        #expect(wurtzite.symbol == "P6_3mc")
        #expect(wurtzite.crystalSystem == "Hexagonal")
    }

    /// An answer without the structure fields (a recording made with the three small fields) has rows without symmetry.
    @Test func aRowWithoutItsStructureHasNoSymmetry() throws {
        let body = #"{"data":[{"id":"mp-7","type":"structures","attributes":{"nsites":1,"chemical_formula_reduced":"Al","_mp_stability":{"r2scan":{"energy_above_hull":0}}}}],"meta":{"data_available":1,"more_data_available":false}}"#
        guard case .rows(let rows, _) = try MaterialsProjectOptimade.decodeSearch(Data(body.utf8)) else {
            Issue.record("expected rows")
            return
        }
        #expect(rows.map(\.symmetry) == [nil])
    }
}

/// Enough of JSON to tell null from a value.
private enum JSONValue: Decodable, Equatable {
    case null
    case value

    init(from decoder: any Decoder) throws {
        self = try decoder.singleValueContainer().decodeNil() ? .null : .value
    }
}
