//
//  TiltRefinementTests.swift
//  Role: the tilt refinement on spot intensities (S5b) on patterns planted with the kit's own `SpotPattern`: a tilted aluminium
//        pattern is refined from the solver's zone-axis start to the planted orientation (or a look-alike flagged as an
//        alternative), the agreement rises, the closed form of the tilt, and the refusals. The planted set of Al and Si at
//        tilts of 0 to 3° is in DiffractoKitParityTests/Geometry/TiltRefinementPlantedTests.swift.
//

import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoGeometry
import DiffractoScattering

@Suite struct TiltRefinementTests {

    private static let fcc = [SIMD3(0.0, 0, 0), SIMD3(0, 0.5, 0.5), SIMD3(0.5, 0, 0.5), SIMD3(0.5, 0.5, 0)]
    private static let aluminium = CrystalStructure.detectingSymmetry(
        cell: .cubic(a: 4.0495), sites: fcc.map { AtomSite(z: 13, fractional: $0) }, conventional: false)
    private static let degree = Double.pi / 180

    @Test(arguments: [(SIMD3(0, 1, 1), 1.5, 0.4), (SIMD3(1, 1, 2), 2.0, 2.2), (SIMD3(0, 1, 3), 1.0, 4.0)])
    func aTiltedPatternIsRefinedToThePlantedOrientation(_ zone: SIMD3<Int>, _ tiltDegrees: Double, _ phi: Double) throws {
        let onZone = try CrystalOrientation(zoneAxis: zone, inPlaneAngle: 0.9, cell: Self.aluminium.cell)
        let t = tiltDegrees * Self.degree
        let truth = onZone.tilted(by: SIMD2(t * cos(phi), t * sin(phi)))
        let pattern = try SpotPattern.simulate(Self.aluminium, orientation: truth, volts: 200_000, kMax: 1.5)
        let strongest = pattern.spots.sorted { $0.intensity > $1.intensity }.prefix(40)
        let spots = strongest.map(\.position), intensities = strongest.map(\.intensity)
        let solver = try ZoneAxisSolver(structure: Self.aluminium, volts: 200_000, kMax: 1.6)
        let solution = solver.solve(spots)
        let result = try #require(solver.refineTilts(of: solution, spots: spots, intensities: intensities))
        let found = ([result.best] + result.alternatives).map {
            solver.angle(between: $0.orientation, and: truth, includingHalfTurn: false)
        }
        #expect(found.min()! <= 0.3 * Self.degree, "\(found.map { $0 / Self.degree })")
        #expect(result.best.agreement > result.best.startAgreement)
        #expect(result.best.agreement > 0.98 && result.best.agreement <= 1 + 1e-12)
        // The agreement of the planted orientation itself is as good as the refined one, within the grid's last step.
        let planted = solver.intensityAgreement(of: truth, scale: 1, spots: spots, intensities: intensities)
        #expect(planted >= result.best.agreement - 1e-3)
    }

    @Test func refusals() throws {
        let solver = try ZoneAxisSolver(structure: Self.aluminium, volts: 200_000, kMax: 1.6)
        let o = try CrystalOrientation(zoneAxis: SIMD3(0, 0, 1), cell: Self.aluminium.cell)
        let spots = [SIMD2(0.49, 0.0), SIMD2(0.0, 0.49)]
        #expect(solver.refineTilt(o, scale: 1, spots: spots, intensities: [1]) == nil)
        #expect(solver.refineTilt(o, scale: 0, spots: spots, intensities: [1, 1]) == nil)
        #expect(solver.refineTilt(o, scale: 1, spots: spots, intensities: [1, 1], settings: .init(step: 0)) == nil)
        #expect(solver.refineTilt(o, scale: 1, spots: [], intensities: []) == nil)
        #expect(solver.intensityAgreement(of: o, scale: 1, spots: spots, intensities: [1]) == -1)
        // Nothing to refine: a refused solution without candidates.
        #expect(solver.refineTilts(of: solver.solve([spots[0]]), spots: [spots[0]], intensities: [1]) == nil)
    }
}
