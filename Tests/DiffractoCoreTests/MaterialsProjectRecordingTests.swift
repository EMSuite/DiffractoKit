//
//  MaterialsProjectRecordingTests.swift
//  Role: the keyless Materials Project route on the service's own answers (DiffractoKit #150, plan row K3a). The bodies
//        in `Fixtures/mp-optimade/` were recorded byte for byte by `tools/materials-project/record.py` on GitHub's Linux
//        runner (its manifest holds each URL, the date and the SHA-256); `pymatgen-pins.json` is pymatgen's reading of
//        the same entries (`tools/pins/mp-optimade/pins.py`, P48). No test touches the network.
//

import CryptoKit
import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore

struct MaterialsProjectRecordingTests {

    typealias Route = MaterialsProjectOptimade

    static let directory = FixtureDirs.fixtures(of: "DiffractoCoreTests").appendingPathComponent("mp-optimade")

    static func data(_ name: String) throws -> Data {
        try Data(contentsOf: directory.appendingPathComponent(name))
    }

    struct Manifest: Decodable {
        struct File: Decodable {
            let file: String
            let url: String
            let sha256: String
            let bytes: Int
        }
        let files: [File]
        let lattices: [String: String]
        let missing: [String]
    }

    struct Pins: Decodable {
        struct Served: Decodable {
            struct Site: Decodable {
                let element: String
                let fractional: [Double]
            }
            let parameters: [Double]
            let sites: [Site]
        }
        struct Conventional: Decodable {
            let parameters: [Double]
            let volume: Double
            let siteCount: Int
        }
        struct Case: Decodable, CustomTestStringConvertible, Sendable {
            let file: String
            let id: String
            let served: Served
            let bravais: String
            let laue: String
            let spaceGroupNumber: Int
            let conventional: Conventional
            var testDescription: String { "\(id) \(bravais)" }
        }
        let cases: [Case]
    }

    static let manifest: Manifest = try! JSONDecoder().decode(Manifest.self, from: data("manifest.json"))
    static let pins: Pins = try! JSONDecoder().decode(Pins.self, from: data("pymatgen-pins.json"))

    /// The recordings are the bytes the service sent: every file the manifest names has its size and SHA-256, and the
    /// pins cover every recorded entry.
    @Test func theRecordingsAreTheBytesTheManifestNames() throws {
        #expect(Self.manifest.files.count >= 3)
        for file in Self.manifest.files {
            let body = try Self.data(file.file)
            #expect(body.count == file.bytes, "\(file.file) size")
            let digest = SHA256.hash(data: body).map { String(format: "%02x", $0) }.joined()
            #expect(digest == file.sha256, "\(file.file) SHA-256")
            #expect(file.url.hasPrefix("https://optimade.materialsproject.org/v1/"), "\(file.file) URL")
        }
        let entries = Set(Self.manifest.files.map(\.file).filter { $0.hasPrefix("entry-") })
        #expect(Set(Self.pins.cases.map(\.file)) == entries)
        #expect(Set(Self.manifest.lattices.values.map { "entry-\($0).json" }) == entries)
        // At least ten of the 14 Bravais lattices are filled; the rest are listed as missing.
        #expect(Self.manifest.lattices.count >= 10)
        #expect(Self.manifest.lattices.count + Self.manifest.missing.count == BravaisLattice.allCases.count)
    }

    /// Every recorded search decodes, most stable first; the formula with no entry is "no results", not an empty service.
    @Test func theRecordedSearchesDecode() throws {
        let searches = Self.manifest.files.map(\.file).filter { $0.hasPrefix("search-") }
        #expect(searches.contains("search-none.json"))
        for name in searches {
            let answer = try Route.decodeSearch(Self.data(name))
            if name == "search-none.json" {
                #expect(answer == .noResults, "\(name)")
                continue
            }
            guard case .rows(let rows, _) = answer else {
                Issue.record("\(name): expected rows, got \(answer)")
                continue
            }
            #expect(!rows.isEmpty, "\(name)")
            let stability = rows.map { $0.lowestEnergyAboveHull ?? .infinity }
            #expect(zip(stability, stability.dropFirst()).allSatisfy { $0 <= $1 }, "\(name) order")
            let formula = String(name.dropFirst("search-".count).dropLast(".json".count))
            #expect(rows.allSatisfy { $0.formula == formula }, "\(name) formula")
        }
    }

    /// The Bravais lattice the structure holds: the family of its crystal system and the centring its detected pure
    /// translations make (1 P; 2 I for the body centre, else C, an A or B face named C as pymatgen's pins do; 3 R; 4 F).
    /// Read from the operations because the route's cell metric alone names a conventional centred cell P, and an
    /// entry served in hexagonal axes (mp-1094122) keeps its R centring only as operations.
    static func bravais(of structure: CrystalStructure) -> String? {
        guard let system = structure.crystalSystem else { return nil }
        let family: String
        switch system {
        case .triclinic: family = "a"
        case .monoclinic: family = "m"
        case .orthorhombic: family = "o"
        case .tetragonal: family = "t"
        case .trigonal, .hexagonal: family = "h"
        case .cubic: family = "c"
        }
        let translations = structure.operations.filter { $0.rotation.isIdentity }.map(\.translationMod1)
        let centring: String
        switch translations.count {
        case 1: centring = "P"
        case 2:
            let bodyCentre = translations.contains { t in
                [t.x, t.y, t.z].allSatisfy { abs($0 - 0.5) < 1e-3 }
            }
            centring = bodyCentre ? "I" : "C"
        case 3: centring = "R"
        case 4: centring = "F"
        default: return nil
        }
        return family + centring
    }

    /// Each recorded entry against pymatgen (P48). The served sites are linear algebra on the same served numbers
    /// (1e-9). The conventional cell is detected by the kit and by spglib through pymatgen: the Bravais lattice, the
    /// Laue class and the site count exact; the volume within 1e-5 relative (`metricLengthTolerance`, ten times the
    /// server's 1e-6 Å rounding); the cell parameters at the same tolerance where the setting is unique (cubic,
    /// tetragonal, hexagonal and rhombohedral in hexagonal axes), the sorted edges for orthorhombic (the two pick the
    /// axis order by different rules), and the volume alone for monoclinic and triclinic (the cell choice differs).
    @Test(arguments: pins.cases) func eachEntryReadsAsPymatgenReadsIt(_ pin: Pins.Case) throws {
        let entry = try Route.decodeEntry(Self.data(pin.file))
        #expect(entry.materialID == pin.id)
        #expect(entry.speciesAtSites.count == pin.served.sites.count)
        let rows = entry.latticeVectors
        let toFractional = simd_double3x3(columns: (rows[0], rows[1], rows[2])).inverse
        for (index, (position, site)) in zip(entry.cartesianPositions, pin.served.sites).enumerated() {
            let fractional = toFractional * position
            for k in 0..<3 {
                #expect(abs(fractional[k] - site.fractional[k]) <= 1e-9, "\(pin.id) site \(index) coordinate \(k)")
            }
        }

        let imported = try Route.structure(from: entry)
        let structure = imported.structure
        #expect(structure.name == pin.id)
        #expect(Self.bravais(of: structure) == pin.bravais, "\(pin.id) Bravais lattice")
        #expect(structure.laueClass?.rawValue == pin.laue, "\(pin.id) Laue class")
        #expect(structure.sites.count == pin.conventional.siteCount, "\(pin.id) site count")
        let tolerance = Route.metricLengthTolerance
        let volume = structure.cell.volume
        #expect(abs(volume - pin.conventional.volume) <= tolerance * pin.conventional.volume, "\(pin.id) volume \(volume)")

        let cell = structure.cell
        let got = [cell.a, cell.b, cell.c, cell.alphaDeg, cell.betaDeg, cell.gammaDeg]
        let want = pin.conventional.parameters
        switch pin.bravais.first {
        case "c", "t", "h":
            for k in 0..<6 {
                #expect(abs(got[k] - want[k]) <= tolerance * max(1, want[k]), "\(pin.id) parameter \(k): \(got[k]) vs \(want[k])")
            }
        case "o":
            for (g, w) in zip(got[0..<3].sorted(), want[0..<3].sorted()) {
                #expect(abs(g - w) <= tolerance * w, "\(pin.id) edge \(g) vs \(w)")
            }
            #expect(got[3...].allSatisfy { abs($0 - 90) <= 1e-9 }, "\(pin.id) angles")
        default:
            break
        }
    }
}
