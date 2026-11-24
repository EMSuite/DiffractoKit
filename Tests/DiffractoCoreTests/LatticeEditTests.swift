//
//  LatticeEditTests.swift
//  Role: lattice edits within the crystal system (DiffractoKit #149, plan row K2). The free parameters of each system,
//        the tied ones following, the fixed ones and the fractional sites kept, the d-spacings of the edited cells
//        against pymatgen 2026.9.24 (`tools/pins/lattice-edit/pins.py` -> `lattice-edit-pins.json`, rtol 1e-7), and
//        the refusals: a parameter that is not free, a value that is not a length or an angle, no cell, and an edit
//        that changes the Laue class.
//

import DiffractoKitTestSupport
import Foundation
import Testing
@testable import DiffractoCore

struct LatticeEditTests {

    struct Pins: Decodable {
        struct Case: Decodable {
            let name: String
            let source: String
            let system: String
            let inputCell: [Double]
            let free: [String]
            let edits: [String: Double]
            let editedCell: [Double]
            let d: [Double]
        }
        let hkl: [[Int]]
        let cases: [Case]
    }

    static func pins() throws -> Pins {
        let url = Fixtures.directory.appendingPathComponent("lattice-edit-pins.json")
        return try JSONDecoder().decode(Pins.self, from: Data(contentsOf: url))
    }

    static var rutile: CrystalStructure {
        let u = 0.3048
        return CrystalPresets.structure(a: 4.594, b: 4.594, c: 2.959, sites: [
            AtomSite(z: 22, fractional: [0, 0, 0]), AtomSite(z: 22, fractional: [0.5, 0.5, 0.5]),
            AtomSite(z: 8, fractional: [u, u, 0]), AtomSite(z: 8, fractional: [1 - u, 1 - u, 0]),
            AtomSite(z: 8, fractional: [0.5 + u, 0.5 - u, 0.5]), AtomSite(z: 8, fractional: [0.5 - u, 0.5 + u, 0.5]),
        ])
    }

    static func structure(_ source: String) throws -> CrystalStructure {
        switch source {
        case "preset:silicon": return CrystalPresets.silicon
        case "preset:magnesium": return CrystalPresets.magnesium
        case "preset:betaDoublePrime": return CrystalPresets.betaDoublePrime
        case "code:rutile": return rutile
        case "code:orthorhombic":
            return CrystalPresets.structure(a: 3, b: 4, c: 5, sites: [
                AtomSite(z: 13, fractional: [0, 0, 0]), AtomSite(z: 8, fractional: [0.5, 0.5, 0.5]),
            ])
        case "code:triclinic":
            return CrystalPresets.structure(a: 5, b: 6, c: 7, alphaDeg: 80, betaDeg: 85, gammaDeg: 95, sites: [
                AtomSite(z: 14, fractional: [0.1, 0.2, 0.3]), AtomSite(z: 8, fractional: [0.4, 0.1, 0.7]),
            ])
        case "code:rhombohedral":
            return CrystalPresets.structure(a: 6.375, b: 6.375, c: 6.375, alphaDeg: 46.08, betaDeg: 46.08, gammaDeg: 46.08,
                                            sites: [AtomSite(z: 20, fractional: [0, 0, 0]),
                                                    AtomSite(z: 6, fractional: [0.5, 0.5, 0.5])])
        case "code:aluminiumPrimitive":
            let a = 4.0495 / 2.0.squareRoot()
            return CrystalPresets.structure(a: a, b: a, c: a, alphaDeg: 60, betaDeg: 60, gammaDeg: 60,
                                            sites: [AtomSite(z: 13, fractional: [0, 0, 0])])
        default:
            let file = try #require(source.hasPrefix("cif:") ? String(source.dropFirst(4)) : nil, "unknown source \(source)")
            return try CIFImport.read(try Fixtures.text(file)).structure
        }
    }

    static func parameters(_ cell: UnitCell) -> [Double] {
        [cell.a, cell.b, cell.c, cell.alphaDeg, cell.betaDeg, cell.gammaDeg]
    }

    @Test func everySystemEditsItsFreeParametersAndMatchesPymatgensDSpacings() throws {
        let pins = try Self.pins()
        #expect(pins.cases.count == 9)
        #expect(pins.hkl.count == 342)
        var compared = 0
        for c in pins.cases {
            let structure = try Self.structure(c.source)
            #expect(structure.crystalSystem?.rawValue == c.system, "\(c.name): \(String(describing: structure.crystalSystem))")
            for (x, y) in zip(Self.parameters(structure.cell), c.inputCell) {
                #expect(abs(x - y) <= 1e-12 * max(1, abs(y)), "\(c.name): input cell")
            }
            #expect(structure.freeLatticeParameters.map(\.rawValue) == c.free, "\(c.name): free \(structure.freeLatticeParameters)")

            var edits: [LatticeParameter: Double] = [:]
            for (key, value) in c.edits { edits[try #require(LatticeParameter(rawValue: key))] = value }
            let edited = try structure.editingLattice(edits)
            #expect(Self.parameters(edited.cell) == c.editedCell, "\(c.name): edited cell \(Self.parameters(edited.cell))")
            #expect(edited.laueClass == structure.laueClass, "\(c.name): Laue class")
            #expect(edited.sites == structure.sites, "\(c.name): sites kept")
            #expect(edited.name == structure.name)
            for (hkl, d) in zip(pins.hkl, c.d) {
                let value = edited.cell.d(MillerIndex(hkl[0], hkl[1], hkl[2]))
                #expect(abs(value - d) <= 1e-7 * d, "\(c.name) \(hkl): \(value) vs pymatgen \(d)")
                compared += 1
            }
        }
        #expect(compared == 9 * 342)
    }

    /// Editing one parameter of a group moves the group; what is not edited keeps its bits.
    @Test func tiedParametersFollowAndTheRestKeepTheirBits() throws {
        let mg = CrystalPresets.magnesium
        let wider = try mg.editingLattice([.a: 3.3])
        #expect(Self.parameters(wider.cell) == [3.3, 3.3, mg.cell.c, 90, 90, 120])
        #expect(wider.cell.c.bitPattern == mg.cell.c.bitPattern)

        let rhombohedral = try Self.structure("code:rhombohedral")
        let sharper = try rhombohedral.editingLattice([.alpha: 45])
        #expect(Self.parameters(sharper.cell) == [6.375, 6.375, 6.375, 45, 45, 45])
        #expect(sharper.laueClass == .threeBarM)

        // An empty edit is the same structure.
        #expect(try mg.editingLattice([:]) == mg)
    }

    @Test func aParameterThatIsNotFreeIsRefusedWithTheFreeOnes() throws {
        let si = CrystalPresets.silicon
        for p in [LatticeParameter.b, .c, .alpha, .beta, .gamma] {
            #expect(throws: LatticeEditError.notFree(p, free: [.a], system: .cubic)) { try si.editingLattice([p: 5]) }
        }
        let error = LatticeEditError.notFree(.c, free: [.a], system: .cubic)
        #expect(error.description.contains("cubic"))
        #expect(error.description.contains("can be set are a."))
        #expect(throws: LatticeEditError.notFree(.gamma, free: [.a, .c], system: .hexagonal)) {
            try CrystalPresets.magnesium.editingLattice([.gamma: 110])
        }
        #expect(throws: LatticeEditError.notFree(.alpha, free: [.a, .b, .c, .beta], system: .monoclinic)) {
            try CrystalPresets.betaDoublePrime.editingLattice([.alpha: 95])
        }
    }

    @Test func aValueThatIsNotALengthOrAnAngleIsRefused() throws {
        let triclinic = try Self.structure("code:triclinic")
        for bad in [0, -1, Double.nan, .infinity] {
            #expect(throws: LatticeEditError.self) { try triclinic.editingLattice([.a: bad]) }
        }
        for bad in [0, 180, 200, -10, Double.nan] {
            #expect(throws: LatticeEditError.self) { try triclinic.editingLattice([.beta: bad]) }
        }
        #expect(throws: LatticeEditError.invalidValue(.a, -1)) { try triclinic.editingLattice([.a: -1]) }
        // Three angles of 120° make a flat cell.
        let flat = #expect(throws: LatticeEditError.self) {
            try triclinic.editingLattice([.alpha: 120, .beta: 120, .gamma: 120])
        }
        if case .invalidCell = flat {} else { Issue.record("expected invalidCell, got \(String(describing: flat))") }
    }

    /// Setting c equal to a on a tetragonal cell whose atoms are cubic gives a cubic crystal: refused, with both
    /// classes named. Rutile's atoms stay tetragonal on the same metric: allowed.
    @Test func anEditThatChangesTheLaueClassIsRefused() throws {
        let tetragonal = CrystalPresets.structure(a: 3, b: 3, c: 4, sites: [AtomSite(z: 26, fractional: [0, 0, 0])])
        #expect(tetragonal.laueClass == .fourOverMmm)
        #expect(throws: LatticeEditError.symmetryChanged(from: .fourOverMmm, to: .mThreeBarM)) {
            try tetragonal.editingLattice([.c: 3])
        }
        let message = LatticeEditError.symmetryChanged(from: .fourOverMmm, to: .mThreeBarM).description
        #expect(message.contains("m-3m") && message.contains("4/mmm"))

        let rutileOnACube = try Self.rutile.editingLattice([.c: 4.594])
        #expect(rutileOnACube.laueClass == .fourOverMmm)
    }

    /// Given operations are kept, and refused when the new metric does not keep them.
    @Test func givenOperationsAreKeptOrRefused() throws {
        let gold = CrystalPresets.gold
        let given = CrystalStructure(cell: gold.cell, sites: gold.sites, operations: gold.operations, name: "gold")
        #expect(given.freeLatticeParameters == [.a])
        let edited = try given.editingLattice([.a: 4.1])
        #expect(edited.operations == gold.operations)
        #expect(edited.symmetrySource == .operationsGiven)
        #expect(throws: LatticeEditError.operationsDoNotFit) {
            try given.replacingCell(UnitCell(a: 4, b: 4, c: 5, alphaDeg: 90, betaDeg: 90, gammaDeg: 90))
        }
    }

    /// The ties read from the rotations, on their own.
    @Test func tiesFollowTheRotations() {
        let free = { (s: CrystalStructure) in s.freeLatticeParameters }
        #expect(free(CrystalPresets.gold) == [.a])
        #expect(free(CrystalPresets.tungstenDisulfide) == [.a, .c])
        let cell = UnitCell(a: 5, b: 6, c: 7, alphaDeg: 90, betaDeg: 90, gammaDeg: 100)
        #expect(LatticeTies(rotations: [.identity], cell: cell).isRepresentative(.gamma))
        let mirrorC = IntegerMatrix([1, 0, 0], [0, 1, 0], [0, 0, -1])
        let ties = LatticeTies(rotations: [.identity, mirrorC], cell: cell)
        #expect(ties.isRepresentative(.gamma))
        #expect(!ties.isRepresentative(.alpha) && !ties.isRepresentative(.beta))
    }

    /// A primitive cell of a face-centred orthorhombic lattice: its rotations carry axes onto sums of axes, so a length
    /// moved alone breaks mmm. Only one scale of the three lengths is offered, and it keeps mmm.
    @Test func aPrimitiveCellOfACentredLatticeOffersOnlyWhatKeepsItsSymmetry() throws {
        let (a, b, c) = (3.0, 4.0, 5.0)
        let v1 = SIMD3<Double>(0, b, c) / 2, v2 = SIMD3<Double>(a, 0, c) / 2, v3 = SIMD3<Double>(a, b, 0) / 2
        func angle(_ x: SIMD3<Double>, _ y: SIMD3<Double>) -> Double {
            acos((x * y).sum() / ((x * x).sum() * (y * y).sum()).squareRoot()) * 180 / .pi
        }
        func length(_ x: SIMD3<Double>) -> Double { (x * x).sum().squareRoot() }
        let primitive = CrystalPresets.structure(
            a: length(v1), b: length(v2), c: length(v3),
            alphaDeg: angle(v2, v3), betaDeg: angle(v1, v3), gammaDeg: angle(v1, v2),
            sites: [AtomSite(z: 13, fractional: [0, 0, 0])])
        #expect(primitive.laueClass == .mmm)
        #expect(primitive.freeLatticeParameters == [.a])
        let scaled = try primitive.editingLattice([.a: length(v1) * 1.02])
        #expect(scaled.laueClass == .mmm)
        #expect(abs(scaled.cell.b / primitive.cell.b - 1.02) < 1e-12)
        #expect(scaled.cell.alphaDeg == primitive.cell.alphaDeg)
    }

    /// A hexagonal crystal in an orthohexagonal C cell. The conventional cell is read from the metric, which is oP for this
    /// cell (the C centring is in the atoms only), so the file and the conventional setting both read mmm and an edit of a
    /// alone is an edit within mmm: allowed, b and c kept (the known limit in docs/parity.md).
    @Test func anOrthohexagonalCellIsEditedWithinMmm() throws {
        let a = 3.2094
        let orthohexagonal = CrystalPresets.structure(a: a, b: a * 3.0.squareRoot(), c: 5.2108, sites: [
            AtomSite(z: 12, fractional: [0, 0, 0]), AtomSite(z: 12, fractional: [0.5, 0.5, 0]),
            AtomSite(z: 12, fractional: [0.5, 1.0 / 6.0, 0.5]), AtomSite(z: 12, fractional: [0, 2.0 / 3.0, 0.5]),
        ])
        #expect(orthohexagonal.laueClass == .mmm)
        #expect(orthohexagonal.conventionalSetting?.laueClass == .mmm)
        #expect(orthohexagonal.freeLatticeParameters == [.a, .b, .c])
        let edited = try orthohexagonal.editingLattice([.a: 3.3])
        #expect(edited.laueClass == .mmm)
        #expect(edited.cell.a == 3.3)
        #expect(edited.cell.b == orthohexagonal.cell.b && edited.cell.c == orthohexagonal.cell.c)
    }

    /// A structure whose symmetry was not read has no crystal system to edit within.
    @Test func aStructureWithoutASymmetryReadingIsNotEdited() {
        let mg = CrystalPresets.magnesium
        let refused = CrystalStructure.detectingSymmetry(cell: mg.cell, sites: mg.sites, tolerance: .fractional(.nan))
        #expect(refused.symmetryRefusal != nil)
        #expect(refused.freeLatticeParameters.isEmpty)
        #expect(throws: LatticeEditError.symmetryNotRead) { try refused.editingLattice([.a: 3.3]) }
    }
}
