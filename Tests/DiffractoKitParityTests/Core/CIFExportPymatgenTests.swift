//
//  CIFExportPymatgenTests.swift
//  Pin C1: the CIF `CIFExport` writes for each structure of the corpus (`CIFExportCorpus`: P2's 19, six with an edited
//  lattice, three Materials Project imports) is the text committed under Fixtures/cif-export/written/, byte for byte, and
//  pymatgen-core 2026.10.16 `CifParser` reads that text (tools/pins/cif-export/pins.py -> pymatgen-read-pins.json) to the
//  structure's cell, relative 1e-7, and its expanded sites as a multiset: element and occupancy equal, coordinates 1e-7
//  modulo 1 against pymatgen's reading with `frac_tolerance=0`, and within 1e-4 (that tolerance's default) of its default
//  reading, which rounds 0.333333 to 1/3 (suite design "Files, figures and exports", pin C1). A text not committed, or committed differently, is printed
//  between `CIFExport C1 begin <name>` and `CIFExport C1 end` lines, so a recording can be taken from the CI log.
//  Red before: `CIFExport` did not exist.
//

import CryptoKit
import DiffractoKitTestSupport
import Foundation
import Testing
import DiffractoCore

struct CIFExportPymatgenTests {
    struct Pins: Decodable {
        struct Site: Decodable { let z: Int; let frac: [Double]; let fracAsWritten: [Double]; let occupancy: Double }
        struct Cell: Decodable { let a, b, c, alpha, beta, gamma: Double }
        struct Structure: Decodable { let name: String; let cifSha256: String; let cell: Cell; let sites: [Site] }
        let pymatgenCore: String
        let structures: [Structure]
    }

    static let directory = Fixtures.directory.appendingPathComponent("cif-export")

    @Test func theKitWritesTheCommittedTexts() throws {
        var differing: [String] = []
        for (name, structure) in try CIFExportCorpus.structures() {
            let text = try CIFExport.text(of: structure).text
            let committed = try? String(contentsOf: Self.directory.appendingPathComponent("written/\(name).cif"), encoding: .utf8)
            if committed != text {
                differing.append(name)
                print("CIFExport C1 begin \(name)\n\(text)CIFExport C1 end")
            }
        }
        #expect(differing.isEmpty, "texts not as committed: \(differing)")
    }

    @Test func pymatgenReadsEachTextToTheStructure() throws {
        let data = try Data(contentsOf: Self.directory.appendingPathComponent("pymatgen-read-pins.json"))
        let pins = try JSONDecoder().decode(Pins.self, from: data)
        #expect(pins.pymatgenCore == "2026.10.16")
        let corpus = try CIFExportCorpus.structures()
        #expect(pins.structures.count == corpus.count)
        for (name, structure) in corpus {
            guard let pin = pins.structures.first(where: { $0.name == name }) else { Issue.record("\(name) not pinned"); continue }
            let text = try String(contentsOf: Self.directory.appendingPathComponent("written/\(name).cif"), encoding: .utf8)
            let sha = SHA256.hash(data: Data(text.utf8)).map { String(format: "%02x", $0) }.joined()
            #expect(sha == pin.cifSha256, "\(name): the pins were not made from the committed text")
            let cell = structure.cell
            for (mine, theirs, what) in [(cell.a, pin.cell.a, "a"), (cell.b, pin.cell.b, "b"), (cell.c, pin.cell.c, "c"),
                                         (cell.alphaDeg, pin.cell.alpha, "alpha"), (cell.betaDeg, pin.cell.beta, "beta"),
                                         (cell.gammaDeg, pin.cell.gamma, "gamma")] {
                #expect(abs(mine - theirs) <= 1e-7 * abs(theirs), "\(name) \(what): \(mine) against pymatgen's \(theirs)")
            }
            #expect(structure.sites.count == pin.sites.count, "\(name) site count")
            var unmatched = Array(structure.sites.indices)
            for site in pin.sites {
                let found = unmatched.firstIndex { i in
                    let mine = structure.sites[i]
                    return mine.z == site.z && mine.occupancy == site.occupancy
                        && (0..<3).allSatisfy { fractionalAxisDelta(mine.fractional[$0], site.fracAsWritten[$0]) <= 1e-7 }
                        && (0..<3).allSatisfy { fractionalAxisDelta(mine.fractional[$0], site.frac[$0]) <= 1e-4 }
                }
                if let found { unmatched.remove(at: found) } else { Issue.record("\(name): pymatgen's Z\(site.z) \(site.frac) is not a site") }
            }
            #expect(unmatched.isEmpty, "\(name): sites pymatgen did not read \(unmatched)")
        }
    }
}
