//
//  SpaceGroupSymbolTests.swift
//  A CIF that names its space group by Hermann–Mauguin symbol or IT number and lists no operations (K6): expanded in
//  the setting pymatgen-core 2026.10.16 `CifParser` takes, from the generated `SpaceGroupTable` (P46,
//  tools/pins/cif-symmetry/space-group-pins.py). Operation sets exact; cell exact; sites 1e-7 against pymatgen.
//

import CryptoKit
import DiffractoKitTestSupport
import Foundation
import Testing
@testable import DiffractoCore

struct SpaceGroupSymbolTests {
    struct Pins: Decodable {
        struct Site: Decodable { let z: Int; let fractional: [Double] }
        struct Reading: Decodable { let cell: [Double]; let sites: [Site] }
        struct PlantedCIF: Decodable {
            let name: String
            let text: String
            let kit: String
            let pymatgen: Reading
            let pymatgenWarnings: [String]
            let hall: String?
            let spglibAtoms: Int?
            let refusedAs: String?
        }
        let spglib: String
        let symbols: [String: String]
        let leftOut: [String]
        let numbers: [String]
        let cifs: [PlantedCIF]
    }

    static func pins() throws -> Pins {
        try JSONDecoder().decode(Pins.self, from: Data(contentsOf: Fixtures.directory.appendingPathComponent("space-group-symbol-pins.json")))
    }

    /// The pins' digest: SHA-256 of the sorted 12-character codes of `HallSymbolTests.encode`, first 16 hex digits.
    static func digest(_ ops: [SymmetryOperation]) -> String {
        let codes = ops.map(HallSymbolTests.encode).sorted().joined()
        return SHA256.hash(data: Data(codes.utf8)).map { String(format: "%02x", $0) }.joined().prefix(16).description
    }

    /// Every symbol and number of the table reads to exactly the operations pymatgen gives it (spglib's for a symbol
    /// pymatgen does not read). RED UNDER: no table.
    @Test func everySymbolAndNumberReadsToPymatgensOperations() throws {
        let pins = try Self.pins()
        #expect(pins.spglib == "2.8.0")
        #expect(pins.symbols.count == SpaceGroupTable.symbolSettings.count + SpaceGroupTable.spglibSymbolSettings.count)
        #expect(pins.numbers.count == 230)
        #expect(Set(pins.leftOut) == SpaceGroupTable.leftOutSymbols)
        var failures: [String] = []
        for (symbol, want) in pins.symbols.sorted(by: { $0.key < $1.key }) {
            guard let setting = SpaceGroupTable.setting(symbols: [symbol], number: nil) else {
                failures.append("\(symbol): no setting"); continue
            }
            let ops = try HallSymbol.operations(setting.hallSymbol)
            if Self.digest(ops) != want { failures.append("\(symbol) -> \(setting.hallSymbol)") }
        }
        for (index, want) in pins.numbers.enumerated() {
            let setting = try #require(SpaceGroupTable.setting(symbols: [], number: index + 1))
            #expect(setting.symbol == nil)
            if Self.digest(try HallSymbol.operations(setting.hallSymbol)) != want {
                failures.append("No. \(index + 1) -> \(setting.hallSymbol)")
            }
        }
        #expect(failures.isEmpty, "\(failures.count) differ: \(failures.prefix(20))")
    }

    @Test func theLookupIsPymatgensOrder() throws {
        // White space and underscores do not count; case does.
        #expect(SpaceGroupTable.setting(symbols: ["P 6_3/m m c"], number: nil)?.hallSymbol == "-P 6c 2c")
        #expect(SpaceGroupTable.setting(symbols: ["p 63/m m c"], number: nil) == nil)
        // A symbol pymatgen reads wins over the number beside it, as in CifParser.get_symops.
        #expect(SpaceGroupTable.setting(symbols: ["F m -3 m"], number: 194)?.hallSymbol == "-F 4 2 3")
        // A short symbol pymatgen does not read comes after the number: P 21/c beside No. 14 is the number's setting.
        let shortWithNumber = try #require(SpaceGroupTable.setting(symbols: ["P 21/c"], number: 14))
        #expect(shortWithNumber.symbol == nil && shortWithNumber.hallSymbol == "-P 2ybc")
        #expect(SpaceGroupTable.setting(symbols: ["P 21/c"], number: nil)?.hallSymbol == "-P 2ybc")
        #expect(SpaceGroupTable.setting(symbols: ["R -3 m :R"], number: nil)?.hallSymbol == "-P 3* 2")
        #expect(SpaceGroupTable.setting(symbols: ["C 1"], number: nil) == nil)
        // A left-out symbol is never replaced by the number beside it (pymatgen reads the symbol first).
        #expect(SpaceGroupTable.setting(symbols: ["A -1"], number: 2) == nil)
        // A symbol neither package reads falls to the number.
        #expect(SpaceGroupTable.setting(symbols: ["F banana"], number: 225)?.hallSymbol == "-F 4 2 3")
        #expect(SpaceGroupTable.setting(symbols: [], number: 0) == nil)
        #expect(SpaceGroupTable.setting(symbols: [], number: 231) == nil)
    }

    /// Planted CIFs with no loop: the cell exact and the sites at 1e-7 as pymatgen reads them; where pymatgen does not
    /// read the symbol (it builds P1 with a warning), the atom count of spglib's setting; "C 1" refused.
    /// RED UNDER: no table (every file refused as missing operations).
    @Test func aCIFNamedBySymbolOrNumberAloneIsExpanded() throws {
        let pins = try Self.pins()
        #expect(pins.cifs.count == 11)
        for planted in pins.cifs {
            switch planted.kit {
            case "pymatgen":
                let doc = try CIFImport.read(planted.text)
                let cell = doc.structure.cell
                #expect([cell.a, cell.b, cell.c, cell.alphaDeg, cell.betaDeg, cell.gammaDeg] == planted.pymatgen.cell,
                        "\(planted.name)")
                #expect(Self.sameSites(doc.structure.sites, planted.pymatgen.sites, tolerance: 1e-7), "\(planted.name)")
                #expect(doc.hadSymmetryOperations, "\(planted.name)")
                #expect(doc.warnings.contains { warning in
                    switch warning {
                    case .operationsFromSpaceGroupSymbol, .operationsFromSpaceGroupNumber: return true
                    default: return false
                    }
                }, "\(planted.name)")
            case "spglib":
                let doc = try CIFImport.read(planted.text)
                #expect(doc.structure.sites.count == planted.spglibAtoms, "\(planted.name)")
                #expect(planted.pymatgen.sites.count < (planted.spglibAtoms ?? 0), "\(planted.name): pymatgen builds P1")
                #expect(planted.pymatgenWarnings.contains { $0.contains("Defaulting to P1") }, "\(planted.name)")
                let symbol = try #require(doc.declaredSpaceGroupSymbol)
                #expect(doc.warnings.contains(.operationsFromSpaceGroupSymbol(symbol: symbol, hallSymbol: planted.hall ?? "")),
                        "\(planted.name)")
            default:
                // pymatgen reads these symbols (before any number) in a change of basis the kit does not hold.
                let symbol = try #require(planted.refusedAs)
                #expect(throws: CIFImportError.missingSymmetryOperations(spaceGroup: symbol), "\(planted.name)") {
                    try CIFImport.read(planted.text)
                }
            }
        }
    }

    @Test func theWarningsNameTheSettingUsed() throws {
        let silicon = try #require(try Self.pins().cifs.first { $0.name == "silicon by number" })
        let doc = try CIFImport.read(silicon.text)
        #expect(doc.warnings == [.operationsFromSpaceGroupNumber(number: 227, hallSymbol: "F 4d 2 3 -1d")])
        #expect(doc.warnings[0].description.contains("standard setting F 4d 2 3 -1d"))
        let rockSalt = try #require(try Self.pins().cifs.first { $0.name == "rock salt by symbol" })
        #expect(try CIFImport.read(rockSalt.text).warnings == [
            .operationsFromSpaceGroupSymbol(symbol: "F m -3 m", hallSymbol: "-F 4 2 3"),
        ])
        // A loop, even the identity alone, is what the file lists: the symbol is not used (pymatgen uses the loop too).
        let identityLoop = rockSalt.text + "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n"
        #expect(throws: CIFImportError.missingSymmetryOperations(spaceGroup: "F m -3 m")) { try CIFImport.read(identityLoop) }
    }

    /// The same atoms as a multiset: equal atomic numbers, positions equal modulo 1 within `tolerance`.
    static func sameSites(_ got: [AtomSite], _ want: [Pins.Site], tolerance: Double) -> Bool {
        guard got.count == want.count else { return false }
        var unused = Array(want.indices)
        for site in got {
            guard let match = unused.firstIndex(where: { index in
                let w = want[index]
                guard w.z == site.z else { return false }
                return (0..<3).allSatisfy { axis in
                    let d = site.fractional[axis] - w.fractional[axis]
                    return abs(d - d.rounded()) <= tolerance
                }
            }) else { return false }
            unused.remove(at: match)
        }
        return true
    }
}
