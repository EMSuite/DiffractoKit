//
//  ZoneAxisSolverTests.swift
//  Role: `ZoneAxisSolver` (S5) on cases with a known answer: aluminium patterns planted with the kit's own `SpotPattern`
//        (solved to the planted orientation, the scale error found, the half turn flagged where the crystal lacks it), the
//        published austenite example (Weirich 2024), one row of spots and too few spots refused, the settings refused.
//        The planted set of eight phases is in DiffractoKitParityTests/Geometry/ZoneAxisSolverPlantedTests.swift.
//

import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoGeometry
import DiffractoScattering

@Suite struct ZoneAxisSolverTests {

    private static let fcc = [SIMD3(0.0, 0, 0), SIMD3(0, 0.5, 0.5), SIMD3(0.5, 0, 0.5), SIMD3(0.5, 0.5, 0)]
    private static let aluminium = CrystalStructure.detectingSymmetry(
        cell: .cubic(a: 4.0495), sites: fcc.map { AtomSite(z: 13, fractional: $0) }, conventional: false)
    private static let degree = Double.pi / 180

    /// The strongest `count` spots of the kit's pattern at `orientation`, 200 kV, |q| ≤ 1.5 Å⁻¹.
    private static func planted(_ structure: CrystalStructure, _ orientation: CrystalOrientation, count: Int = 40)
        throws -> [SIMD2<Double>] {
        let p = try SpotPattern.simulate(structure, orientation: orientation, volts: 200_000, kMax: 1.5)
        return p.spots.sorted { $0.intensity > $1.intensity }.prefix(count).map(\.position)
    }

    @Test func aluminiumDown001IsSolvedToThePlantedOrientation() throws {
        #expect(Self.aluminium.laueClass == .mThreeBarM)
        let truth = try CrystalOrientation(zoneAxis: SIMD3(0, 0, 1), inPlaneAngle: 0.7, cell: Self.aluminium.cell)
        let spots = try Self.planted(Self.aluminium, truth)
        let solver = try ZoneAxisSolver(structure: Self.aluminium, volts: 200_000, kMax: 1.6)
        let solution = solver.solve(spots)
        let best = try #require(solution.solved, "\(solution.verdict)")
        #expect(solver.angle(between: best.orientation, and: truth) <= 0.1 * Self.degree)
        #expect(abs(best.scale - 1) <= 5e-4)
        #expect(best.zoneAxisAngle <= 1e-9)
        #expect(Set([abs(best.zoneAxis.u), abs(best.zoneAxis.v), abs(best.zoneAxis.w)]) == Set([0, 1]))
        #expect(best.matched == spots.count && best.agreement > 0.95)
        // The half turn about [001] is a symmetry of the cube: nothing to flag.
        #expect(!best.halfTurnAmbiguous)
        // Every spot gets an hkl of the zero layer of the solved zone axis, at its own |g|.
        for (spot, hkl) in zip(spots, best.assignments) {
            let h = try #require(hkl)
            #expect(abs(Self.aluminium.cell.gLength(h) * best.scale - simd_length(spot)) <= 1e-9)
        }
    }

    @Test func aScaleErrorIsFound() throws {
        let truth = try CrystalOrientation(zoneAxis: SIMD3(0, 1, 1), inPlaneAngle: 2.1, cell: Self.aluminium.cell)
        let spots = try Self.planted(Self.aluminium, truth).map { $0 * 1.02 }
        let solver = try ZoneAxisSolver(structure: Self.aluminium, volts: 200_000, kMax: 1.6)
        let best = try #require(solver.solve(spots).solved)
        #expect(abs(best.scale - 1.02) <= 5e-4)
        #expect(solver.angle(between: best.orientation, and: truth) <= 0.1 * Self.degree)
    }

    @Test func theHalfTurnIsFlaggedWhenTheCrystalLacksIt() throws {
        // Down [1 1 2] of aluminium the zero layer is centrosymmetric and no higher layer is reached, but no two-fold axis of
        // m-3m lies along [1 1 2]: the half turn predicts the same spots and is a different orientation.
        let truth = try CrystalOrientation(zoneAxis: SIMD3(1, 1, 2), inPlaneAngle: 0.3, cell: Self.aluminium.cell)
        let spots = try Self.planted(Self.aluminium, truth)
        let solver = try ZoneAxisSolver(structure: Self.aluminium, volts: 200_000, kMax: 1.6)
        let best = try #require(solver.solve(spots).solved)
        #expect(best.halfTurnAmbiguous)
        let plain = solver.angle(between: best.orientation, and: truth, includingHalfTurn: false)
        let turned = solver.angle(between: best.orientation.rotatedInPlane(by: .pi), and: truth, includingHalfTurn: false)
        // One of the two is the planted orientation; the other is a different orientation (the half turn composed with the
        // cube's rotations is 60° away).
        #expect(min(plain, turned) <= 0.1 * Self.degree)
        #expect(max(plain, turned) >= 1 * Self.degree)
    }

    /// Weirich, J. Appl. Cryst. (2024), doi 10.1107/S1600576724010215, section 5.1: austenitic steel, an F cell
    /// with a = 3.704 Å, 200 kV; three measured spots d = 2.149, 1.304 and 1.116 Å at 0°, 89.03° and 59.25°; published
    /// answer 111, 220, 311 and zone axis [112]. The paper's own tolerance is 3 % on the lengths.
    @Test func thePublishedAusteniteExampleIsSolved() throws {
        let austenite = CrystalStructure.detectingSymmetry(
            cell: .cubic(a: 3.704), sites: Self.fcc.map { AtomSite(z: 26, fractional: $0) }, conventional: false)
        let spots = zip([2.149, 1.304, 1.116], [0.0, 89.03, 59.25]).map { d, a in
            SIMD2(cos(a * Self.degree), sin(a * Self.degree)) / d
        }
        let solver = try ZoneAxisSolver(structure: austenite, volts: 200_000, kMax: 1.6,
                                        settings: .init(tolerance: 0.03))
        let best = try #require(solver.solve(spots, clicked: true).solved)
        let zone = [abs(best.zoneAxis.u), abs(best.zoneAxis.v), abs(best.zoneAxis.w)].sorted()
        #expect(zone == [1, 1, 2])
        let families = try best.assignments.map { h -> [Int] in
            let m = try #require(h)
            return [abs(m.h), abs(m.k), abs(m.l)].sorted()
        }
        #expect(families == [[1, 1, 1], [0, 2, 2], [1, 1, 3]])
        #expect(best.halfTurnAmbiguous)
        #expect(abs(best.scale - 1) <= 0.03)
    }

    @Test func oneRowAndTooFewSpotsAreRefused() throws {
        let solver = try ZoneAxisSolver(structure: Self.aluminium, volts: 200_000, kMax: 1.6)
        let g = SIMD2(2.0 / 4.0495, 0).rotated(0.4)
        let row = (1...4).flatMap { [g * Double($0), -g * Double($0)] }
        #expect(solver.solve(row).verdict == .refused(.oneRow))
        #expect(solver.solve(row, clicked: true).verdict == .refused(.oneRow))
        #expect(solver.solve([g]).verdict == .refused(.tooFewSpots(1)))
        #expect(solver.solve([g, .zero, SIMD2(.nan, 1)]).verdict == .refused(.tooFewSpots(1)))
        // Two directions whose lengths no pair of aluminium's reflections has.
        #expect(solver.solve([SIMD2(0.11, 0), SIMD2(0, 0.13)]).verdict == .refused(.noCandidate))
    }

    @Test func badSettingsAreRefused() throws {
        for bad in [ZoneAxisSolver.Settings(tolerance: 0), .init(floor: .nan), .init(scaleWindow: -0.1),
                    .init(seeds: 1), .init(maxIndex: 0), .init(minimumAgreement: 1.5),
                    .init(pattern: .init(sigmaExcitationError: 0))] {
            #expect(throws: ZoneAxisSolverError.self) {
                try ZoneAxisSolver(structure: Self.aluminium, volts: 200_000, settings: bad)
            }
        }
        #expect(throws: ZoneAxisSolverError.pattern(.invalidVoltage(-1))) {
            try ZoneAxisSolver(structure: Self.aluminium, volts: -1)
        }
    }
}

private extension SIMD2 where Scalar == Double {
    func rotated(_ a: Double) -> SIMD2<Double> { SIMD2(cos(a) * x - sin(a) * y, sin(a) * x + cos(a) * y) }
}
