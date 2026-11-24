//
//  TiltRefinementPlantedTests.swift
//  Role: the tilt refinement (S5b) on the planted set, read from Tests/DiffractoGeometryTests/Fixtures/tilt-planted.json
//        (tools/pins/tilt-refinement/planted.py: py4DSTEM f050d207's `generate_diffraction_pattern` at known tilted
//        orientations, positions and intensities with noise). No Python function refines the tilt of one pattern alone, so
//        the pin is the planted truth (docs/parity.md P44): Al and Si, 13 zone-axis classes to index 3, tilts 0, 0.5, 1, 2
//        and 3°, two orientations each, 260 patterns.
//  Done when (plan row S5b): from the solver's start the refined orientation (or a look-alike returned as an alternative) is
//  within 1° of the planted one, the median under 0.2°; a start that ends more than 1° off agrees less than the right one,
//  or agrees as well and is returned as an alternative. The worst error, 0.6°, is the regression bound (NumPy prototype of
//  the same arithmetic: medians 0.09 to 0.14° per tilt, worst 0.51°).
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoGeometry

private struct TiltCell: Decodable { let a, b, c, alpha, beta, gamma: Double }
private struct TiltSite: Decodable { let frac: [Double]; let occupancy: Double; let z: Int }
private struct TiltCase: Decodable {
    let zone: [Int]
    let tilt: Double
    let matrix: [Double]
    /// qx, qy, intensity per spot
    let spots: [Double]
}
private struct TiltPhase: Decodable { let name: String; let cell: TiltCell; let sites: [TiltSite]; let cases: [TiltCase] }
private struct TiltFile: Decodable { let versions: [String: String]; let volts: Double; let phases: [TiltPhase] }

private enum TiltPlanted {
    static let file: TiltFile = {
        let url = FixtureDirs.fixtures(of: "DiffractoGeometryTests").appendingPathComponent("tilt-planted.json")
        do { return try JSONDecoder().decode(TiltFile.self, from: Data(contentsOf: url)) } catch { fatalError("\(error)") }
    }()
    static var names: [String] { file.phases.map(\.name) }
}

@Suite("Tilt refinement on the planted set (py4DSTEM f050d207 patterns at known tilts)")
struct TiltRefinementPlantedTests {

    @Test(arguments: TiltPlanted.names)
    func refinedWithinOneDegreeAndRightStartsAgreeBest(_ name: String) throws {
        let phase = try #require(TiltPlanted.file.phases.first { $0.name == name })
        #expect(TiltPlanted.file.versions["py4DSTEM"] == "0.14.19")
        let c = phase.cell
        let structure = CrystalStructure.detectingSymmetry(
            cell: UnitCell(a: c.a, b: c.b, c: c.c, alphaDeg: c.alpha, betaDeg: c.beta, gammaDeg: c.gamma),
            sites: phase.sites.map { AtomSite(z: $0.z, fractional: SIMD3($0.frac[0], $0.frac[1], $0.frac[2]),
                                              occupancy: $0.occupancy) },
            conventional: false)
        let solver = try ZoneAxisSolver(structure: structure, volts: TiltPlanted.file.volts, kMax: 1.6)
        let degree = Double.pi / 180
        var errors: [Double: [Double]] = [:]
        for p in phase.cases {
            let m = p.matrix
            let truth = try CrystalOrientation(matrix: simd_double3x3(rows: [SIMD3(m[0], m[1], m[2]), SIMD3(m[3], m[4], m[5]),
                                                                              SIMD3(m[6], m[7], m[8])]))
            let rows = stride(from: 0, to: p.spots.count, by: 3)
            let spots = rows.map { SIMD2(p.spots[$0], p.spots[$0 + 1]) }, intensities = rows.map { p.spots[$0 + 2] }
            let label = "\(name) \(p.zone) tilt \(p.tilt)"
            guard let result = solver.refineTilts(of: solver.solve(spots), spots: spots, intensities: intensities) else {
                Issue.record("\(label): nothing refined"); continue
            }
            func off(_ r: RefinedOrientation) -> Double { solver.angle(between: r.orientation, and: truth, includingHalfTurn: false) }
            let reported = [result.best] + result.alternatives
            let error = reported.map(off).min()!
            #expect(error <= 1 * degree, "\(label): \(error / degree)°")
            errors[p.tilt, default: []].append(error)
            // Every start that ends more than 1° off agrees less than the best right one, or is reported (a look-alike).
            let right = result.refined.filter { off($0) <= 1 * degree }.map(\.agreement).max() ?? -2
            for r in result.refined where off(r) > 1 * degree && r.agreement >= right {
                let flagged = reported.contains { solver.angle(between: $0.orientation, and: r.orientation, includingHalfTurn: false) <= 1 * degree }
                #expect(flagged, "\(label): a wrong start agrees \(r.agreement) against \(right) and is not reported")
            }
        }
        for (tilt, e) in errors.sorted(by: { $0.key < $1.key }) {
            let sorted = e.sorted()
            let median = sorted[sorted.count / 2]
            print("TiltRefinement planted \(name) tilt \(tilt)°: median \(median / degree)°, worst \(sorted.last! / degree)°")
            #expect(median <= 0.2 * degree, "\(name) tilt \(tilt): median \(median / degree)°")
            #expect(sorted.last! <= 0.6 * degree, "\(name) tilt \(tilt): worst \(sorted.last! / degree)°")
        }
    }
}
