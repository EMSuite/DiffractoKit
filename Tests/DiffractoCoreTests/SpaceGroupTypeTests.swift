//
//  SpaceGroupTypeTests.swift
//  The space-group type from operations (K7) against spglib 2.8.0 `get_spacegroup_type_from_symmetry` on all 530
//  settings with seeded random origin shifts (P47, tools/pins/cif-symmetry/space-group-type-pins.py), and against
//  pymatgen's SpacegroupAnalyzer number on the P2 structures. Exact: number, Hall number, symbol.
//

import DiffractoKitTestSupport
import Foundation
import simd
import Testing
import DiffractoCore

struct SpaceGroupTypeTests {
    struct Pins: Decodable {
        struct Setting: Decodable {
            let hall: Int
            let hallSymbol: String
            let shift: [Double]
            let number: Int
            let hallNumber: Int
            let international: String
        }
        let spglib: String
        let settings: [Setting]
    }

    static func pins() throws -> Pins {
        try JSONDecoder().decode(Pins.self, from: Data(contentsOf: Fixtures.directory.appendingPathComponent("space-group-type-pins.json")))
    }

    /// The operations of a Hall symbol with the origin moved by `shift`: t' = t + (I − R) p, mod 1.
    static func shifted(_ symbol: String, by shift: SIMD3<Double>) throws -> [SymmetryOperation] {
        try HallSymbol.operations(symbol).map { op in
            let r = op.rotation.rows.map { SIMD3<Double>(Double($0.x), Double($0.y), Double($0.z)) }
            let moved = SIMD3(shift.x - simd_dot(r[0], shift), shift.y - simd_dot(r[1], shift), shift.z - simd_dot(r[2], shift))
            let t = op.translation + moved
            return SymmetryOperation(rotation: op.rotation, translation: t - t.rounded(.down))
        }
    }

    /// RED UNDER: no naming (every setting nil).
    @Test func everySettingWithAnyOriginIsNamedAsSpglibNamesIt() throws {
        let pins = try Self.pins()
        #expect(pins.spglib == "2.8.0")
        #expect(pins.settings.count == 530)
        var failures: [String] = []
        for setting in pins.settings {
            let ops = try Self.shifted(setting.hallSymbol, by: SIMD3(setting.shift[0], setting.shift[1], setting.shift[2]))
            guard let type = SpaceGroupType.of(ops) else { failures.append("\(setting.hall): nil"); continue }
            if type.number != setting.number || type.hallNumber != setting.hallNumber
                || type.internationalSymbol != setting.international {
                failures.append("\(setting.hall): \(type.number) \(type.hallNumber) \(type.internationalSymbol)")
            }
            // The setting found is one with exactly these operations up to the origin (three pairs of No. 68 settings
            // list the same operations; the first is found).
            let same = try HallSymbol.operations(type.settingHallSymbol)
            if Set(same.map(\.rotation)) != Set(ops.map(\.rotation)) { failures.append("\(setting.hall): setting \(type.settingHallNumber)") }
        }
        #expect(failures.isEmpty, "\(failures.count) settings differ: \(failures.prefix(20))")
    }

    @Test func notAGroupInASettingIsNamedNothing() throws {
        let rockSalt = try HallSymbol.operations("-F 4 2 3")
        // A missing operation, a stray translation, a basis no setting uses, nothing at all.
        #expect(SpaceGroupType.of(rockSalt.filter { $0.rotation != .inversion }) == nil)
        #expect(SpaceGroupType.of(rockSalt.map { $0.rotation == .inversion ? SymmetryOperation(rotation: .inversion, translation: SIMD3(0.1, 0, 0)) : $0 }) == nil)
        // One centring copy moved off the lattice: the set is not closed.
        var broken = rockSalt
        let index = try #require(broken.lastIndex { $0.rotation == .inversion })
        broken[index] = SymmetryOperation(rotation: .inversion, translation: broken[index].translation + SIMD3(0.25, 0, 0))
        #expect(SpaceGroupType.of(broken) == nil)
        let skew = IntegerMatrix([0, 1, 0], [1, 1, 0], [0, 0, 1])
        #expect(SpaceGroupType.of([.identity, SymmetryOperation(rotation: skew)]) == nil)
        #expect(SpaceGroupType.of([]) == nil)
        #expect(SpaceGroupType.of([SymmetryOperation(rotation: .identity, translation: SIMD3(.nan, 0, 0))]) == nil)
        // A translation off by more than the tolerance is refused; within it, named.
        let p21c = try Self.shifted("-P 2ybc", by: SIMD3(0.1, 0.2, 0.3))
        let nudged = p21c.map { SymmetryOperation(rotation: $0.rotation, translation: $0.translation + SIMD3(0, 0.0004, 0)) }
        #expect(SpaceGroupType.of(nudged)?.number == 14)
        let far = p21c.map { $0.rotation == .identity ? $0 : SymmetryOperation(rotation: $0.rotation, translation: $0.translation + SIMD3(0, 0.05, 0)) }
        #expect(SpaceGroupType.of(far) == nil)
        // A centring copy left out: not a whole group.
        let rockSaltMinusOne = Array(rockSalt.filter { $0.rotation == .inversion }.dropFirst())
            + rockSalt.filter { $0.rotation != .inversion }
        #expect(SpaceGroupType.of(rockSaltMinusOne) == nil)
    }

    /// The structures of P2 (the committed CIFs and beta''): the type read from the atoms is pymatgen's
    /// SpacegroupAnalyzer number (symprec 1e-4). RED UNDER: no naming.
    @Test func theCIFStructuresAreNamedAsPymatgenNamesThem() throws {
        struct Cells: Decodable {
            struct Structure: Decodable {
                let name: String
                let cif_text: String?
                let cif_file: String?
                let space_group_number: Int
            }
            let structures: [Structure]
        }
        let pins = try JSONDecoder().decode(Cells.self, from: Data(contentsOf: Fixtures.directory.appendingPathComponent("pymatgen-cells-pins.json")))
        // The known limit: a cube holding an mm2 structure (two-fold along [110]) is in no setting's basis.
        let unnamed: Set<String> = ["tetragonal_a_equals_c"]
        var read = 0
        for entry in pins.structures {
            let text: String
            if let inline = entry.cif_text {
                text = inline
            } else if let file = entry.cif_file {
                text = try String(contentsOf: Fixtures.directory.appendingPathComponent(file), encoding: .utf8)
            } else {
                continue
            }
            let structure = try CIFImport.read(text).structure
            if unnamed.contains(entry.name) {
                #expect(structure.spaceGroupType == nil, "\(entry.name)")
            } else {
                #expect(structure.spaceGroupType?.number == entry.space_group_number, "\(entry.name)")
            }
            read += 1
        }
        #expect(read == 18)
        // The symbol is spglib's short one.
        let silicon = try CIFImport.read(try String(contentsOf: Fixtures.directory.appendingPathComponent("silicon_Fd-3m.cif"), encoding: .utf8))
        #expect(silicon.structure.spaceGroupType?.internationalSymbol == "Fd-3m")
    }

    /// A relaxed cell can read a lower group than its source states: both stay readable, neither replaces the other.
    /// The primitive cell of fcc strained 2 % along c, as the Materials Project serves a primitive cell, stated 225.
    /// The parameters are computed, not typed: rounded to seven digits they are no longer tetragonal at the kit's
    /// metric tolerance (1e-9) and read C2/m.
    @Test func aStatedNumberAndTheDetectedTypeAreBothKept() {
        let a = 4.0, c = 4.08
        let rows = [SIMD3(0, a / 2, c / 2), SIMD3(a / 2, 0, c / 2), SIMD3(a / 2, a / 2, 0)]
        func angle(_ i: Int, _ j: Int) -> Double {
            acos(simd_dot(rows[i], rows[j]) / (simd_length(rows[i]) * simd_length(rows[j]))) * 180 / .pi
        }
        let cell = UnitCell(a: simd_length(rows[0]), b: simd_length(rows[1]), c: simd_length(rows[2]),
                            alphaDeg: angle(1, 2), betaDeg: angle(0, 2), gammaDeg: angle(0, 1))
        let structure = CrystalStructure.detectingSymmetry(cell: cell, sites: [AtomSite(z: 13, fractional: .zero)],
                                                           spaceGroupNumber: 225, name: "relaxed")
        #expect(structure.spaceGroupNumber == 225)
        #expect(structure.spaceGroupType?.number == 139)
        #expect(structure.spaceGroupType?.internationalSymbol == "I4/mmm")
    }
}
