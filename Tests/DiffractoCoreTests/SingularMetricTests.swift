//
//  SingularMetricTests.swift
//  Role: silent wrong values refused (DiffractoKit #120, first block; plan row K9). On 0.18.1 a singular metric used the
//        identity as its inverse, so d, g and |g| came back as plausible numbers (d(100) = 1 Å for a cell with a = 0), and
//        the absolute 1e-15 bound on det G read a valid 1e-5 Å cubic cell as singular in the same way. A NaN metric passed
//        `preservesMetric`. Each case here failed on 0.18.1.
//

import Foundation
import Testing
import simd
import DiffractoCore

struct SingularMetricTests {

    @Test func aDegenerateCellHasNoReciprocalLattice() {
        // a = 0 (det G = 0) and three 120 degree angles (coplanar axes): both `degenerateVolume`, both NaN, never a number.
        for cell in [UnitCell(a: 0, b: 4, c: 4, alphaDeg: 90, betaDeg: 90, gammaDeg: 90),
                     UnitCell.rhombohedral(a: 5, alphaDeg: 120)] {
            #expect(cell.issues.contains(.degenerateVolume))
            #expect(cell.d(MillerIndex(1, 0, 0)).isNaN, "\(cell)")
            #expect(cell.gLength(MillerIndex(0, 1, 0)).isNaN, "\(cell)")
            #expect(cell.reciprocalMetricTensor.allSatisfy { $0.x.isNaN && $0.y.isNaN && $0.z.isNaN })
            #expect(cell.interplanarAngleDeg(MillerIndex(1, 0, 0), MillerIndex(0, 1, 0)).isNaN, "was 180")
        }
    }

    /// The threshold is relative to the cell: a 1e-5 Å cube is a valid (if absurd) cell and has its exact d.
    @Test(arguments: [1e-5, 1e-3, 4.0782, 1e5])
    func aSmallOrLargeCubeHasItsD(_ a: Double) {
        let cell = UnitCell.cubic(a: a)
        #expect(cell.isValid)
        #expect(abs(cell.d(MillerIndex(1, 0, 0)) - a) <= 1e-12 * a)
        #expect(abs(cell.d(MillerIndex(1, 1, 1)) - a / 3.0.squareRoot()) <= 1e-12 * a)
    }

    /// A smoke check that a valid monoclinic cell is inverted as before (the pymatgen parity pins hold the numbers).
    @Test func aValidCellIsStillInverted() {
        let cell = UnitCell.monoclinic(a: 14.83, b: 4.05, c: 6.74, betaDeg: 105.3)
        let g = cell.metricTensor
        let det = g[0].x * (g[1].y * g[2].z - g[1].z * g[2].y) - g[0].y * (g[1].x * g[2].z - g[1].z * g[2].x)
            + g[0].z * (g[1].x * g[2].y - g[1].y * g[2].x)
        let expected = (g[1].y * g[2].z - g[1].z * g[2].y) / det
        #expect(abs(cell.reciprocalMetricTensor[0].x - expected) <= 1e-14 * abs(expected))
        #expect(cell.d(MillerIndex(2, 0, 0)).isFinite)
    }

    /// A valid but extreme cell whose det G leaves Double's range has no reflection list.
    @Test func aCellBeyondDoublesRangeHasNoReflections() {
        let tiny = UnitCell.cubic(a: 1e-60)
        #expect(tiny.isValid)
        #expect(tiny.d(MillerIndex(1, 0, 0)).isNaN)
        #expect(throws: ReflectionListError.cellTooLargeOrInvalid(estimatedReflections: .infinity)) {
            try ReflectionList(cell: tiny, operations: [.identity], laueGroup: PointSymmetry(closureOf: []), laueClass: nil, kMax: 1)
        }
    }

    @Test func aNaNMetricPreservesNothing() {
        let nan = UnitCell(a: .nan, b: 4, c: 4, alphaDeg: 90, betaDeg: 90, gammaDeg: 90)
        // The identity is no exception: a NaN metric is no metric.
        #expect(!SymmetryDetection.preservesMetric(.identity, of: nan))
        #expect(!SymmetryDetection.preservesMetric(IntegerMatrix([0, -1, 0], [1, 0, 0], [0, 0, 1]), of: nan))
        #expect(SymmetryDetection.preservesMetric(.identity, of: .cubic(a: 4)))
    }
}
