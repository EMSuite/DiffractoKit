//
//  PatternOrientationTests.swift
//  Role: `PatternOrientation.index` (S6) on aluminium patterns planted with the kit's own `SpotPattern` (positions and
//        intensities): a tilted zone axis comes back to its full orientation, the [001] and [011] zone axes of DiffractoKit #65
//        come back exact when the spots stop at the detector edge (0.75 Å⁻¹) while the templates reach 1.5 Å⁻¹, and the
//        refusals. The py4DSTEM-planted sets are in DiffractoKitParityTests/Orientation/PatternOrientationPlantedTests.swift.
//

import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoGeometry
import DiffractoOrientation

@Suite struct PatternOrientationTests {

    private static let aluminium: CrystalStructure = {
        let fcc = [SIMD3(0.0, 0, 0), SIMD3(0, 0.5, 0.5), SIMD3(0.5, 0, 0.5), SIMD3(0.5, 0.5, 0)]
        return CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0495), sites: fcc.map { AtomSite(z: 13, fractional: $0) },
                                                  conventional: false)
    }()
    private static let degree = Double.pi / 180
    private static let indexer: PatternOrientation = {
        do { return try PatternOrientation(structure: aluminium, volts: 200_000) } catch { fatalError("\(error)") }
    }()

    /// The strongest 40 spots of the kit's pattern at `orientation` (200 kV, |q| ≤ 1.5 Å⁻¹, or ≤ `cut`), with intensities.
    private static func planted(_ orientation: CrystalOrientation, cut: Double = .infinity)
        throws -> (spots: [SIMD2<Double>], intensities: [Double]) {
        let p = try SpotPattern.simulate(aluminium, orientation: orientation, volts: 200_000, kMax: 1.5)
        let kept = p.spots.filter { simd_length($0.position) <= cut }.sorted { $0.intensity > $1.intensity }.prefix(40)
        return (kept.map(\.position), kept.map(\.intensity))
    }

    /// The smallest angle from the truth over the best and its alternatives (radians; the half turn about the beam counts).
    private static func error(_ r: PatternOrientationResult, _ truth: CrystalOrientation) -> (best: Double, any: Double) {
        let solver = indexer.solver
        let best = solver.angle(between: r.best.orientation, and: truth, includingHalfTurn: false)
        let any = ([r.best] + r.alternatives).map { solver.angle(between: $0.orientation, and: truth, includingHalfTurn: false) }
            .min() ?? .infinity
        return (best, any)
    }

    @Test func aTiltedZoneAxisComesBackToItsFullOrientation() throws {
        let truth = try CrystalOrientation(zoneAxis: SIMD3(0, 1, 1), inPlaneAngle: 0.4, cell: Self.aluminium.cell)
            .tilted(by: SIMD2(1.2, -0.7) * Self.degree)
        let (spots, intensities) = try Self.planted(truth)
        let r = try #require(Self.indexer.index(spots: spots, intensities: intensities))
        let e = Self.error(r, truth)
        #expect(e.any <= 0.3 * Self.degree, "best \(e.best / Self.degree)°, any \(e.any / Self.degree)°")
        #expect(r.best.agreement > 0.98)
        #expect(r.best.positions.matched == spots.count)
        #expect(abs(r.best.scale - 1) <= 5e-4)
    }

    /// DiffractoKit #65: the template library reaches 1.5 Å⁻¹, the pattern stops at the detector edge (0.75 Å⁻¹), and the
    /// template step alone answers a beam direction a few degrees off [001] and [011]. The full route comes back to the zone axis.
    @Test(arguments: [(SIMD3(0, 0, 1), 0.3), (SIMD3(0, 1, 1), 1.1)])
    func zoneAxesCutAtTheDetectorEdgeComeBackExact(_ zone: SIMD3<Int>, _ inPlane: Double) throws {
        let truth = try CrystalOrientation(zoneAxis: zone, inPlaneAngle: inPlane, cell: Self.aluminium.cell)
        let (spots, intensities) = try Self.planted(truth, cut: 0.75)
        let r = try #require(Self.indexer.index(spots: spots, intensities: intensities))
        let e = Self.error(r, truth)
        #expect(e.any <= 0.5 * Self.degree, "\(zone): best \(e.best / Self.degree)°, any \(e.any / Self.degree)°")
        #expect(r.best.positions.matched == spots.count)
    }

    @Test func refusals() throws {
        let truth = try CrystalOrientation(zoneAxis: SIMD3(0, 0, 1), cell: Self.aluminium.cell)
        let (spots, intensities) = try Self.planted(truth)
        #expect(Self.indexer.index(spots: spots, intensities: Array(intensities.dropLast())) == nil)
        #expect(Self.indexer.index(spots: [spots[0], SIMD2(.nan, 0)], intensities: [1, 1]) == nil)
        for bad in [PatternOrientation.Settings(templateMatches: 0, solverCandidates: 0),
                    .init(templateMatches: -1), .init(tilt: .init(step: 0)), .init(lookAlikeMargin: .nan)] {
            #expect(throws: PatternOrientationError.self) {
                try PatternOrientation(structure: Self.aluminium, volts: 200_000, settings: bad)
            }
        }
        #expect(throws: PatternOrientationError.self) { try PatternOrientation(structure: Self.aluminium, volts: -1) }
    }
}
