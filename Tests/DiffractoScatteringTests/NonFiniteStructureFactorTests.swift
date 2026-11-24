//
//  NonFiniteStructureFactorTests.swift
//  Role: `StructureFactors.compute(_:in:)` refuses (nil) a structure whose volume or an occupancy or coordinate is not a
//        finite number (DiffractoKit #120, first block; plan row K9). On 0.18.1 it returned F = NaN or ±inf; only the
//        kinematic path validated the structure.
//

import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoScattering

struct NonFiniteStructureFactorTests {

    static func gold(occupancy: Double = 1, x: Double = 0, cell: UnitCell = .cubic(a: 4.0782)) -> CrystalStructure {
        let sites = [SIMD3<Double>(x, 0, 0), SIMD3(0.5, 0.5, 0), SIMD3(0.5, 0, 0.5), SIMD3(0, 0.5, 0.5)].enumerated().map {
            AtomSite(z: 79, fractional: $0.element, occupancy: $0.offset == 0 ? occupancy : 1)
        }
        return CrystalStructure(cell: cell, sites: sites, operations: [.identity])
    }

    @Test func aFiniteStructureHasItsF() throws {
        let f = try #require(StructureFactors.compute(MillerIndex(1, 1, 1), in: Self.gold()))
        #expect(f.real.isFinite && f.imaginary.isFinite && f.intensity > 0)
    }

    @Test(arguments: [Double.nan, .infinity, -.infinity])
    func aNonFiniteOccupancyOrCoordinateIsRefused(_ value: Double) {
        #expect(StructureFactors.compute(MillerIndex(1, 1, 1), in: Self.gold(occupancy: value)) == nil)
        #expect(StructureFactors.compute(MillerIndex(1, 1, 1), in: Self.gold(x: value)) == nil)
    }

    @Test func aVolumeThatIsNotAFinitePositiveNumberIsRefused() {
        let nanCell = UnitCell(a: .nan, b: 4, c: 4, alphaDeg: 90, betaDeg: 90, gammaDeg: 90)
        #expect(StructureFactors.compute(MillerIndex(1, 1, 1), in: Self.gold(cell: nanCell)) == nil)
        let flat = UnitCell(a: 0, b: 4, c: 4, alphaDeg: 90, betaDeg: 90, gammaDeg: 90)
        #expect(StructureFactors.compute(MillerIndex(1, 1, 1), in: Self.gold(cell: flat)) == nil)
        let huge = UnitCell.cubic(a: 1e200)
        #expect(StructureFactors.compute(MillerIndex(1, 1, 1), in: Self.gold(cell: huge)) == nil)   // V = inf
    }
}
