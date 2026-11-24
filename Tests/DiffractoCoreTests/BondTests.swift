//
//  BondTests.swift
//  VESTA's bonds against pymatgen-core 2026.10.16 CutOffDictNN.from_preset("vesta_2019") (P39,
//  tools/pins/vesta-bonds/pins.py): site indices and images exact, lengths at rtol 1e-7 (assert_allclose's
//  float64 default; pymatgen ships no tolerance for get_nn_info).
//

import DiffractoKitTestSupport
import Foundation
import simd
import Testing
import DiffractoCore

struct BondTests {
    struct Pins: Decodable {
        struct Site: Decodable { let z: Int; let fractional: [Double] }
        struct Case: Decodable {
            let name: String
            let a: Double, b: Double, c: Double, alpha: Double, beta: Double, gamma: Double
            let sites: [Site]
            let bonds: [[Double]]
        }
        let pymatgenCore: String
        let cases: [Case]
    }

    static func pins() throws -> Pins {
        try JSONDecoder().decode(Pins.self, from: Data(contentsOf: Fixtures.directory.appendingPathComponent("vesta-bond-pins.json")))
    }

    static func structure(_ pin: Pins.Case) -> CrystalStructure {
        let cell = UnitCell(a: pin.a, b: pin.b, c: pin.c, alphaDeg: pin.alpha, betaDeg: pin.beta, gammaDeg: pin.gamma)
        let sites = pin.sites.map { AtomSite(z: $0.z, fractional: SIMD3($0.fractional[0], $0.fractional[1], $0.fractional[2])) }
        return CrystalStructure(cell: cell, sites: sites, operations: [.identity], name: pin.name)
    }

    @Test func theBondsArePymatgensOnSixStructures() throws {
        let pins = try Self.pins()
        #expect(pins.pymatgenCore == "2026.10.16")
        #expect(pins.cases.map(\.name) == ["silicon", "silicon shifted", "rocksalt NaCl", "aluminium", "strontium titanate",
                                          "baddeleyite ZrO2"])
        for pin in pins.cases {
            let bonds = Self.structure(pin).bonds
            let expected = pin.bonds.map {
                Bond(site: Int($0[0]), neighbour: Int($0[1]), image: SIMD3(Int($0[2]), Int($0[3]), Int($0[4])), length: $0[5])
            }
            #expect(bonds.count == expected.count, "\(pin.name)")
            for (bond, want) in zip(bonds, expected) {
                #expect(bond.site == want.site && bond.neighbour == want.neighbour && bond.image == want.image,
                        "\(pin.name): \(bond) against \(want)")
                #expect(abs(bond.length - want.length) <= 1e-7 * abs(want.length), "\(pin.name): \(bond.length) against \(want.length)")
            }
        }
    }

    @Test func theCoordinationsAreTheIssuesCounts() throws {
        let pins = try Self.pins()
        func coordination(_ name: String) throws -> [Int: Set<Int>] {
            let pin = try #require(pins.cases.first { $0.name == name })
            let structure = Self.structure(pin)
            var count = [Int](repeating: 0, count: structure.sites.count)
            let bonds = structure.bonds
            for bond in bonds { count[bond.site] += 1 }
            var byElement: [Int: Set<Int>] = [:]
            for (index, site) in structure.sites.enumerated() { byElement[site.z, default: []].insert(count[index]) }
            return byElement
        }
        #expect(try coordination("silicon") == [14: [4]])
        #expect(try coordination("silicon shifted") == [14: [4]])
        #expect(try coordination("rocksalt NaCl") == [11: [6], 17: [6]])
        #expect(try coordination("aluminium") == [13: [0]])
        #expect(try coordination("strontium titanate") == [38: [12], 22: [6], 8: [6]])
        #expect(try coordination("baddeleyite ZrO2") == [40: [7], 8: [3, 4]])
    }

    @Test func everyBondIsListedFromBothAtoms() throws {
        for pin in try Self.pins().cases {
            let bonds = Self.structure(pin).bonds
            let set = Set(bonds.map { Key($0) })
            #expect(bonds.allSatisfy { set.contains(Key($0.reversed)) }, "\(pin.name)")
        }
    }

    private struct Key: Hashable {
        let site: Int, neighbour: Int, image: SIMD3<Int>
        init(_ bond: Bond) { site = bond.site; neighbour = bond.neighbour; image = bond.image }
    }

    @Test func theTableIsPymatgensFileWithItsLastEntryWinning() {
        // H-O and D-O are listed twice in the file (2.1, then 1.2); pymatgen keeps the last.
        #expect(BondCutoffs.vesta(1, 8) == 1.2)
        #expect(BondCutoffs.vesta(8, 1) == 1.2)
        #expect(BondCutoffs.vesta(14, 14) == 2.6)
        #expect(BondCutoffs.vesta(11, 17) == BondCutoffs.vesta(17, 11))
        // Metals: no pair, so no bond.
        #expect(BondCutoffs.vesta(13, 13) == nil)
        #expect(BondCutoffs.vesta(26, 26) == nil)
        #expect(BondCutoffs.vesta(0, 8) == nil)
        #expect(BondCutoffs.vesta(-1, 200) == nil)
    }

    @Test func aPartlyOccupiedSiteBondsByItsElement() {
        // VESTA draws it; pymatgen keys by species string ("Na:0.5") and gives none (docs/parity.md).
        let cell = UnitCell(a: 5.6402, b: 5.6402, c: 5.6402, alphaDeg: 90, betaDeg: 90, gammaDeg: 90)
        let structure = CrystalStructure(cell: cell, sites: [
            AtomSite(z: 11, fractional: SIMD3(0, 0, 0), occupancy: 0.5),
            AtomSite(z: 17, fractional: SIMD3(0.5, 0, 0)),
        ], operations: [.identity])
        let fromSodium = structure.bonds.filter { $0.site == 0 }
        #expect(fromSodium.count == 2)
        #expect(fromSodium.allSatisfy { abs($0.length - 2.8201) < 1e-12 })
    }

    @Test func twoElementsAtOnePositionAreAZeroLengthBondButOneElementIsOneAtom() {
        let cell = UnitCell(a: 10, b: 10, c: 10, alphaDeg: 90, betaDeg: 90, gammaDeg: 90)
        let structure = CrystalStructure(cell: cell, sites: [
            AtomSite(z: 14, fractional: SIMD3(0.1, 0.1, 0.1)),
            AtomSite(z: 8, fractional: SIMD3(0.1, 0.1, 0.1)),
        ], operations: [.identity])
        #expect(structure.bonds == [
            Bond(site: 0, neighbour: 1, image: .zero, length: 0),
            Bond(site: 1, neighbour: 0, image: .zero, length: 0),
        ])
        // Two silicon sites at one position: one atom to pymatgen (same species within 1e-8 Å), so no bond.
        let doubled = CrystalStructure(cell: cell, sites: [
            AtomSite(z: 14, fractional: SIMD3(0.1, 0.1, 0.1)),
            AtomSite(z: 14, fractional: SIMD3(0.1, 0.1, 0.1 + 1e-10)),
        ], operations: [.identity])
        #expect(doubled.bonds.isEmpty)
    }

    @Test func aCustomCutoffAndUnusableNumbers() {
        let cell = UnitCell(a: 4.0495, b: 4.0495, c: 4.0495, alphaDeg: 90, betaDeg: 90, gammaDeg: 90)
        let sites = [SIMD3<Double>(0, 0, 0), SIMD3(0, 0.5, 0.5), SIMD3(0.5, 0, 0.5), SIMD3(0.5, 0.5, 0)]
            .map { AtomSite(z: 13, fractional: $0) }
        let aluminium = CrystalStructure(cell: cell, sites: sites, operations: [.identity])
        #expect(aluminium.bonds.isEmpty)
        // A cut-off just above the nearest-neighbour distance a/sqrt(2) gives fcc's twelve.
        let twelve = aluminium.bonds(cutoff: { _, _ in 2.9 })
        #expect(twelve.count == 48)
        #expect(twelve.allSatisfy { abs($0.length - 4.0495 / 2.0.squareRoot()) < 1e-12 })
        #expect(aluminium.bonds(cutoff: { _, _ in .nan }).isEmpty)
        #expect(aluminium.bonds(cutoff: { _, _ in -1 }).isEmpty)

        let far = CrystalStructure(cell: cell, sites: [
            AtomSite(z: 14, fractional: SIMD3(0, 0, 0)), AtomSite(z: 8, fractional: SIMD3(1e12, 0, 0)),
            AtomSite(z: 8, fractional: SIMD3(.nan, 0, 0)),
        ], operations: [.identity])
        #expect(far.bonds.isEmpty)
    }
}
