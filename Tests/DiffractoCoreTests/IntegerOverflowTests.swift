//
//  IntegerOverflowTests.swift
//  Role: the integer sites of DiffractoKit #117 that trapped (a crashed test process on 0.18.1): the closure of a
//        generator that is not unimodular (its powers leave `Int` in about 63 steps), a Miller index or matrix product past
//        `Int`, and the orthorhombic C branch of the conventional cell, which force-unwrapped a face class a loose
//        tolerance need not give. Products now wrap (exact whenever the result fits, as for every symmetry rotation) and
//        the C branch falls back to the reduced cell instead of trapping.
//

import Foundation
import Testing
import simd
import DiffractoCore

struct IntegerOverflowTests {

    @Test func theClosureOfANonUnimodularGeneratorStops() {
        let doubling = IntegerMatrix([2, 0, 0], [0, 1, 0], [0, 0, 1])
        let group = PointSymmetry(closureOf: [doubling])
        #expect(group.order <= 200)
        #expect(doubling.order == nil)
    }

    @Test func productsPastIntWrapInsteadOfTrapping() {
        let big = IntegerMatrix([Int.max, 0, 0], [0, 1, 0], [0, 0, 1])
        let square = big.multiplied(by: big)
        #expect(square[0, 0] == Int.max &* Int.max)
        #expect(MillerIndex(Int.max, 1, 0).applying(big).h == Int.max &* Int.max)
        #expect(IntegerMatrix([Int.min, 0, 0], [0, 1, 0], [0, 0, 1]).negated[0, 0] == Int.min)
        #expect(IntegerMatrix([Int.max, 0, 0], [0, Int.max, 0], [0, 0, 1]).trace == Int.max &+ Int.max &+ 1)
        // Exact where it fits: the symmetry rotations are unchanged.
        let fourfold = IntegerMatrix([0, -1, 0], [1, 0, 0], [0, 0, 1])
        #expect(fourfold.order == 4)
        #expect(MillerIndex(1, 2, 3).applying(fourfold) == MillerIndex(2, -1, 3))
    }

    /// Near-orthorhombic C cells (a primitive setting, a close to b) under loose metric tolerances: the route answers a cell
    /// for every one, never traps, and two runs agree (the face is chosen from a sorted list, not a set's order).
    @Test func looseTolerancesOnCentredCellsNeverTrap() {
        for t in [1e-9, 1e-3, 0.02, 0.05, 0.1, 0.2, 0.3] {
            for (a, b, c) in [(3.0, 3.01, 5.0), (3.0, 4.0, 5.0), (4.0, 4.0, 4.02), (2.9, 5.1, 5.0), (5.0, 5.05, 5.1)] {
                // the primitive cell of the C lattice with conventional a, b, c: a' = (a - b)/2, b' = (a + b)/2, c' = c
                let ap = 0.5 * (a * a + b * b).squareRoot()
                let gamma = 2 * atan(b / a) * 180 / .pi
                let cell = UnitCell(a: ap, b: ap, c: c, alphaDeg: 90, betaDeg: 90, gammaDeg: gamma)
                let first = cell.conventionalCell(tolerance: t)
                let second = cell.conventionalCell(tolerance: t)
                #expect(first.bravais == second.bravais, "t \(t), \(a) \(b) \(c)")
                #expect(first.transformation == second.transformation, "t \(t), \(a) \(b) \(c)")
            }
        }
    }
}
