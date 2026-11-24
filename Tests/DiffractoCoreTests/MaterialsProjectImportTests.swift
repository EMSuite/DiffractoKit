//
//  MaterialsProjectImportTests.swift
//  The Materials Project JSON models, the document-to-structure route and the expectation verdicts. Numerics against
//  pymatgen SpacegroupAnalyzer (fixtures/mp-import/pymatgen-pins.json, tools/pins/mp-import/pins.py; pymatgen 2026.9.24,
//  spglib 2.8.0): the conventional cell's volume and site count (volume 1e-7 relative (the reference's tolerance; was 1e-9): the MP inputs carry their cell to 3-4
//  digits and pymatgen's standard cell is rebuilt from lengths and angles), the Laue class and the metric family exact.
//  Verdicts, refusals and the formula parser are exact logic. Nothing touches the network.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore

struct MaterialsProjectImportTests {
    private struct Pins: Decodable {
        struct Lat: Decodable { let lattice: [[Double]]; let sites: [Site] }
        struct Site: Decodable { let abc: [Double]; let species: [String: Double] }
        struct Case: Decodable { let input: String; let spaceGroupNumber: Int; let spaceGroupSymbol: String; let conventional: Lat; let primitive: Lat }
        let cases: [Case]
    }
    private static let dir = URL(fileURLWithPath: #filePath).deletingLastPathComponent().appendingPathComponent("Fixtures/mp-import")
    private static let pins = try! JSONDecoder().decode(Pins.self, from: Data(contentsOf: dir.appendingPathComponent("pymatgen-pins.json")))

    private func det(_ m: [[Double]]) -> Double {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    @Test func servedInputsGiveTheConventionalStructurePymatgenFinds() throws {
        for c in Self.pins.cases {
            let data = try Data(contentsOf: Self.dir.appendingPathComponent(c.input))
            let doc = try MaterialsProjectImport.firstDocument(in: MaterialsProjectImport.decode(data))
            let imported = try MaterialsProjectImport.structure(from: doc)
            #expect(imported.structure.sites.count == c.conventional.sites.count, "\(c.input)")
            let volume = abs(det(c.conventional.lattice))
            #expect(abs(imported.structure.cell.volume - volume) <= 1e-7 * volume, "\(c.input) volume")
            let expectedFamily: UnitCell.MetricFamily = c.spaceGroupNumber >= 195 ? .cubic : .hexagonal
            #expect(imported.structure.cell.metricFamily(lengthTolerance: 1e-4, angleToleranceDeg: 1e-2) == expectedFamily, "\(c.input)")
            #expect(imported.structure.laueClass == (c.spaceGroupNumber >= 195 ? .mThreeBarM : .sixOverMmm), "\(c.input)")
            #expect(imported.structure.name == doc.materialID)
        }
    }

    @Test func aPrimitiveCentredCellIsConvertedToPymatgensConventionalCell() throws {
        for c in Self.pins.cases {
            // A document built from pymatgen's primitive standard cell, as MP serves centred lattices.
            let p = c.primitive
            let sites = p.sites.map { site in
                MaterialsProjectSite(species: site.species.map { MaterialsProjectSpecies(element: $0.key, occu: $0.value) }, abc: site.abc)
            }
            let l = p.lattice
            func len(_ r: [Double]) -> Double { (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).squareRoot() }
            func ang(_ u: [Double], _ v: [Double]) -> Double {
                acos((u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (len(u) * len(v))) * 180 / .pi
            }
            let lattice = MaterialsProjectLattice(a: len(l[0]), b: len(l[1]), c: len(l[2]), alpha: ang(l[1], l[2]),
                                                  beta: ang(l[0], l[2]), gamma: ang(l[0], l[1]), matrix: l)
            let doc = MaterialsProjectDocument(
                materialID: "mp-1", structure: MaterialsProjectStructure(lattice: lattice, sites: sites),
                symmetry: MaterialsProjectSymmetry(crystalSystem: c.spaceGroupNumber >= 195 ? "cubic" : "hexagonal",
                                                   symbol: c.spaceGroupSymbol, number: c.spaceGroupNumber))
            let imported = try MaterialsProjectImport.structure(from: doc)
            #expect(imported.structure.sites.count == c.conventional.sites.count, "\(c.input)")
            let volume = abs(det(c.conventional.lattice))
            #expect(abs(imported.structure.cell.volume - volume) <= 1e-7 * volume, "\(c.input) volume")
            #expect(imported.structure.laueClass == (c.spaceGroupNumber >= 195 ? .mThreeBarM : .sixOverMmm), "\(c.input)")
            // Fm-3m and Fd-3m start with F: converted; P6_3/mmc starts with P: as served.
            if c.spaceGroupSymbol.hasPrefix("F") { #expect(imported.cell == .conventional(centring: "F"), "\(c.input)") }
            else { #expect(imported.cell == .asServed, "\(c.input)") }
        }
    }

    @Test func envelopeDecodingAndRefusalsAreNamed() throws {
        // `errors` absent on a success; unknown keys ignored.
        let ok = try MaterialsProjectImport.decode(Data(#"{"data":[{"material_id":"mp-1","nsites":2,"formula_pretty":"Al"}],"meta":{"total_doc":1}}"#.utf8))
        #expect(ok.errors.isEmpty && ok.data.count == 1)
        // errors are checked before data.
        let failed = MaterialsProjectSummaryResponse(data: [], errors: [MaterialsProjectError(code: 401, message: "bad key")])
        #expect(throws: MaterialsProjectImport.Failure.apiError(code: 401, message: "bad key")) {
            _ = try MaterialsProjectImport.firstDocument(in: failed)
        }
        #expect(throws: MaterialsProjectImport.Failure.emptyResult("the response's data array was empty.")) {
            _ = try MaterialsProjectImport.firstDocument(in: MaterialsProjectSummaryResponse(data: [], errors: []))
        }
        for bad in ["", "mp-", "mp-12a", "MP-12", "12", "mp- 12"] {
            #expect(throws: MaterialsProjectImport.Failure.invalidMaterialID(bad)) {
                _ = try MaterialsProjectImport.request(materialID: bad, apiKey: "k")
            }
        }
        let request = try MaterialsProjectImport.request(materialID: "mp-149", apiKey: "secret")
        #expect(request.value(forHTTPHeaderField: "X-API-KEY") == "secret" && request.httpMethod == "GET")
        #expect(request.url?.absoluteString.contains("material_ids=mp-149") == true)
        // Structural absence.
        #expect(throws: MaterialsProjectImport.Failure.malformed("mp-9 has no structure block.")) {
            _ = try MaterialsProjectImport.structure(from: MaterialsProjectDocument(materialID: "mp-9"))
        }
        let lattice = MaterialsProjectLattice(a: 3, b: 3, c: 3, alpha: 90, beta: 90, gamma: 90)
        #expect(throws: MaterialsProjectImport.Failure.malformed("mp-9's structure lists no atomic sites.")) {
            _ = try MaterialsProjectImport.structure(from: MaterialsProjectDocument(materialID: "mp-9", structure: MaterialsProjectStructure(lattice: lattice, sites: [])))
        }
        let unknown = MaterialsProjectSite(species: [MaterialsProjectSpecies(element: "Xx", occu: 1)], abc: [0, 0, 0])
        #expect(throws: MaterialsProjectImport.Failure.malformed("Site 0 of mp-9 has an unrecognized element symbol \"Xx\".")) {
            _ = try MaterialsProjectImport.structure(from: MaterialsProjectDocument(materialID: "mp-9", structure: MaterialsProjectStructure(lattice: lattice, sites: [unknown])))
        }
        // A centred primitive cell whose stated system does not fit is read from its atoms, as spglib reads it (#218):
        // the fcc primitive cell stated hexagonal is cubic F, a = 2.9 √2.
        let primitive = MaterialsProjectLattice(a: 2.9, b: 2.9, c: 2.9, alpha: 60, beta: 60, gamma: 60)
        let site = MaterialsProjectSite(species: [MaterialsProjectSpecies(element: "Al", occu: 1)], abc: [0, 0, 0])
        let doc = MaterialsProjectDocument(materialID: "mp-9", structure: MaterialsProjectStructure(lattice: primitive, sites: [site]),
                                           symmetry: MaterialsProjectSymmetry(crystalSystem: "hexagonal", symbol: "Fm-3m", number: 225))
        let read = try MaterialsProjectImport.structure(from: doc)
        #expect(read.structure.crystalSystem == .cubic)
        #expect(read.structure.laueClass == .mThreeBarM)
        #expect(read.structure.sites.count == 4)
        #expect(abs(read.structure.cell.a - 2.9 * 2.0.squareRoot()) <= 1e-6 * 4.1)
        #expect(read.cell == .conventional(centring: "F"))
        #expect(read.symprec == 5e-4)  // the tight reading: hexagonal is reached at no step
    }

    @Test func expectationVerdicts() {
        let doc = MaterialsProjectDocument(materialID: "mp-1185307", formulaPretty: "Al2Cu",
                                           symmetry: MaterialsProjectSymmetry(symbol: "I4/mcm", number: 140))
        func verdict(_ f: String?, _ n: Int?) -> PhaseExpectationVerdict {
            MaterialsProjectImport.check(doc, against: PhaseExpectation(formula: f, spaceGroupNumber: n))
        }
        #expect(verdict(nil, nil) == .unchecked)
        #expect(verdict("Al2Cu", nil) == .matches)
        #expect(verdict("CuAl2", 140) == .matches)
        #expect(verdict("Al2Cu", 62) == .mismatch(reason: "mp-1185307 is Al2Cu in I4/mcm (140); expected space group 62."))
        #expect(verdict("AlCu", nil) == .mismatch(reason: "mp-1185307 is Al2Cu; expected AlCu."))
        #expect(verdict("al2cu", nil) == .mismatch(reason: "Could not parse expected formula \"al2cu\"."))
        #expect(verdict(nil, 140) == .matches)
        let bare = MaterialsProjectDocument(materialID: "mp-2")
        #expect(MaterialsProjectImport.check(bare, against: PhaseExpectation(formula: "Al")) == .mismatch(reason: "mp-2 has no usable formula to compare against \"Al\"."))
        #expect(MaterialsProjectImport.check(bare, against: PhaseExpectation(spaceGroupNumber: 225)) ==
                .mismatch(reason: "mp-2 is mp-2 in an unrecorded space group; expected space group 225."))
    }

    @Test func formulaParser() {
        #expect(MaterialsProjectImport.parseFormula("LiAl2Cu") == ["Li": 1, "Al": 2, "Cu": 1])
        #expect(MaterialsProjectImport.parseFormula("Mg5Si6") == ["Mg": 5, "Si": 6])
        #expect(MaterialsProjectImport.parseFormula("AlAl2") == ["Al": 3])
        for bad in ["", "al", "Al0", "Al-2", "(AlCu)2", "Al 2", "2Al", "É"] { #expect(MaterialsProjectImport.parseFormula(bad) == nil, "\(bad)") }
    }

    /// mp-1079297 (Si, C2/m) as the main API served it on 2026-11-24 (#151, the id put back in its legacy form):
    /// pymatgen's primitive cell, which the fixed C-centring matrix does not take to a conventional monoclinic cell. It
    /// is read as the keyless route reads it, to the conventional cell pymatgen gives for the same atoms (P48).
    static let monoclinicSilicon = Data(#"""
{"data":[{"formula_pretty":"Si","symmetry":{"crystal_system":"Monoclinic","symbol":"C2/m","hall":"-C 2y","number":12,"point_group":"2/m","symprec":0.1,"angle_tolerance":5.0,"version":"2.7.0"},"material_id":"mp-1079297","deprecated":false,"last_updated":"2021-07-15T02:29:39.366000","structure":{"@module":"pymatgen.core.structure","@class":"Structure","charge":0.0,"lattice":{"matrix":[[-0.0,3.84593083,0.0],[-0.37362086,0.0,6.32197728],[6.90216228,-1.92296642,-0.40892854]],"pbc":[true,true,true],"a":3.84593083,"b":6.333007916927259,"c":7.176689107292369,"alpha":96.52403067130003,"beta":105.54208897884945,"gamma":90.0,"volume":167.23079840226276},"properties":{},"sites":[{"species":[{"element":"Si","occu":1.0}],"abc":[0.05552936,0.87396351,0.11105773],"properties":{},"label":"Si","xyz":[0.4400074766936057,1.7921227421920184e-06,5.479762778384439]},{"species":[{"element":"Si","occu":1.0}],"abc":[0.94447064,0.12603649,0.88894227],"properties":{},"label":"Si","xyz":[6.088533943306394,1.9229626178772574,0.4332859616155613]},{"species":[{"element":"Si","occu":1.0}],"abc":[0.55899985,0.65339544,0.11799771],"properties":{},"label":"Si","xyz":[0.5703171768755004,1.9229691231134771,4.08249849526196]},{"species":[{"element":"Si","occu":1.0}],"abc":[0.44100015,0.34660456,0.88200229],"properties":{},"label":"Si","xyz":[5.9582242431245,-4.713113477446404e-06,1.83055024473804]},{"species":[{"element":"Si","occu":1.0}],"abc":[0.21274317,0.05730781,0.42548434],"properties":{},"label":"Si","xyz":[2.9153505690217787,3.4183190682939696e-06,0.1883059828374932]},{"species":[{"element":"Si","occu":1.0}],"abc":[0.78725683,0.94269219,0.57451566],"properties":{},"label":"Si","xyz":[3.613190850978221,1.9229609916809316,5.724742757162507]},{"species":[{"element":"Si","occu":1.0}],"abc":[0.72710984,0.58557765,0.45421767],"properties":{},"label":"Si","xyz":[2.9163000435937088,1.9229688236717255,3.5162660303404896]},{"species":[{"element":"Si","occu":1.0}],"abc":[0.27289016,0.41442235,0.54578233],"properties":{},"label":"Si","xyz":[3.612241376406291,-4.413671725876489e-06,2.39678270965951]}]}}]}
"""#.utf8)

    @Test func aPrimitiveCellTheFixedMatrixCannotConvertIsReadFromTheLattice() throws {
        let doc = try MaterialsProjectImport.firstDocument(in: MaterialsProjectImport.decode(Self.monoclinicSilicon))
        let imported = try MaterialsProjectImport.structure(from: doc)
        #expect(imported.cell == .conventional(centring: "C"))
        let cell = imported.structure.cell
        // pymatgen 2026.9.24 SpacegroupAnalyzer(symprec=1e-3).get_conventional_standard_structure of the same entry.
        let expected = [13.828530896705164, 3.84593083, 6.333007916927258, 90.0, 96.77278162787648, 90.0]
        let got = [cell.a, cell.b, cell.c, cell.alphaDeg, cell.betaDeg, cell.gammaDeg]
        for (g, e) in zip(got, expected) { #expect(abs(g - e) <= 1e-6 * e, "\(got)") }
        #expect(imported.structure.sites.count == 16)
        #expect(imported.structure.laueClass?.rawValue == "2/m")
    }
}
