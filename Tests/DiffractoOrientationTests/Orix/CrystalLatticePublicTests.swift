//
//  CrystalLatticePublicTests.swift
//  Role: `CrystalLattice.standardBase` and `CrystalLattice.inverse` are reachable without `@testable` (DiffractoKit #267):
//        the app's Kikuchi structure factors call them. Red before: both were internal, so this file did not compile.
//

import DiffractoOrientation
import Foundation
import Testing

struct CrystalLatticePublicTests {

    private static func angle(_ u: SIMD3<Double>, _ v: SIMD3<Double>) -> Double {
        let d = (u * v).sum() / ((u * u).sum().squareRoot() * (v * v).sum().squareRoot())
        return acos(min(1, max(-1, d))) * 180 / .pi
    }

    /// The rows are the cell's edges: their lengths and the angles between them are the cell's.
    @Test(arguments: [[4.0, 4, 4, 90, 90, 90], [3.2094, 3.2094, 5.2105, 90, 90, 120], [5.1, 6.2, 7.3, 81, 97, 104]])
    func theStandardBaseRowsAreTheCellsEdges(_ p: [Double]) {
        let base = CrystalLattice.standardBase(a: p[0], b: p[1], c: p[2], alpha: p[3], beta: p[4], gamma: p[5])
        #expect(base.count == 3)
        for (row, length) in zip(base, p[0...2]) {
            #expect(abs((row * row).sum().squareRoot() - length) <= 1e-12 * length)
        }
        #expect(abs(Self.angle(base[1], base[2]) - p[3]) <= 1e-9)
        #expect(abs(Self.angle(base[0], base[2]) - p[4]) <= 1e-9)
        #expect(abs(Self.angle(base[0], base[1]) - p[5]) <= 1e-9)
    }

    @Test func theInverseUndoesTheBase() {
        let base = CrystalLattice.standardBase(a: 5.1, b: 6.2, c: 7.3, alpha: 81, beta: 97, gamma: 104)
        let inverse = CrystalLattice.inverse(base)
        // Rows times rows: (M M⁻¹)ᵢⱼ = Σₖ Mᵢₖ (M⁻¹)ₖⱼ.
        for i in 0..<3 {
            for j in 0..<3 {
                let entry = (0..<3).reduce(0.0) { $0 + base[i][$1] * inverse[$1][j] }
                #expect(abs(entry - (i == j ? 1 : 0)) <= 1e-12, "(\(i), \(j))")
            }
        }
    }
}
