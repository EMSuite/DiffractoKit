//
//  ZoneAxisSolverPlantedTests.swift
//  Role: `ZoneAxisSolver` (S5) on the planted set, read from Tests/DiffractoGeometryTests/Fixtures/zone-axis-planted.json
//        (tools/pins/zone-axis-solver/planted.py: py4DSTEM f050d207's `generate_diffraction_pattern` at known orientations).
//        No Python function solves one pattern alone, so the pin is the planted truth (docs/parity.md P43): 8 phases, the
//        first 14 zone-axis classes to index 3, two in-plane turns each, five conditions (clean, scale +2 %, missing and
//        false spots, 3 and 6 clicked zero-layer spots), 1070 patterns.
//  Done when (plan row S5): every pattern is right (within 0.1° on a full pattern, 0.3° from clicked spots; the scale within
//  0.05 % on a full pattern; with 1 px noise, missing and false spots 0.2° and 0.15 %), or a tie that holds the planted
//  orientation, or refused; never a wrong answer. The counts of
//  right, tie and refused per phase are the kit's regression bound (measured on the prototype of the same arithmetic).
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoGeometry
import DiffractoScattering

private struct PlantedCell: Decodable { let a, b, c, alpha, beta, gamma: Double }
private struct PlantedSite: Decodable { let frac: [Double]; let occupancy: Double; let z: Int }
private struct PlantedCase: Decodable {
    let zone: [Int]
    let condition: String
    let matrix: [Double]
    let scale: Double
    /// qx0, qy0, qx1, qy1, ... (Å⁻¹)
    let spots: [Double]
}
private struct PlantedPhase: Decodable { let name: String; let cell: PlantedCell; let sites: [PlantedSite]; let cases: [PlantedCase] }
private struct PlantedFile: Decodable { let versions: [String: String]; let volts: Double; let phases: [PlantedPhase] }

private enum Planted {
    static let file: PlantedFile = {
        let url = FixtureDirs.fixtures(of: "DiffractoGeometryTests").appendingPathComponent("zone-axis-planted.json")
        do { return try JSONDecoder().decode(PlantedFile.self, from: Data(contentsOf: url)) } catch { fatalError("\(error)") }
    }()
    static var names: [String] { file.phases.map(\.name) }

    /// The least right patterns and the most refusals per phase over the 5 conditions (full, clicked): the regression bound.
    /// Measured on the NumPy prototype of the same arithmetic (right, tie, refused): Al 97/7/26, Au 100/12/18,
    /// theta' 58/79/3, beta'' 102/32/6, Mg 112/20/8, Si 99/11/20, Fe 89/10/31, Mg2Si 116/4/10; the bound leaves 3 either way
    /// for rounding at the margins.
    static let bounds: [String: (right: Int, refused: Int)] = [
        "Al": (94, 29), "Au": (97, 21), "theta_prime_Al2Cu": (55, 6), "beta_double_prime_Mg5Si6": (99, 9),
        "Mg": (109, 11), "Si": (96, 23), "Fe": (86, 34), "Mg2Si": (113, 13),
    ]
}

private enum Outcome: String { case right, tie, refused, wrong }

@Suite("ZoneAxisSolver on the planted set (py4DSTEM f050d207 patterns at known orientations)")
struct ZoneAxisSolverPlantedTests {

    @Test(arguments: Planted.names)
    func everyPatternIsRightATieOrRefused(_ name: String) throws {
        let phase = try #require(Planted.file.phases.first { $0.name == name })
        #expect(Planted.file.versions["py4DSTEM"] == "0.14.19")
        let c = phase.cell
        let structure = CrystalStructure.detectingSymmetry(
            cell: UnitCell(a: c.a, b: c.b, c: c.c, alphaDeg: c.alpha, betaDeg: c.beta, gammaDeg: c.gamma),
            sites: phase.sites.map { AtomSite(z: $0.z, fractional: SIMD3($0.frac[0], $0.frac[1], $0.frac[2]),
                                              occupancy: $0.occupancy) },
            conventional: false)
        let solver = try ZoneAxisSolver(structure: structure, volts: Planted.file.volts, kMax: 1.6)
        let degree = Double.pi / 180
        var counts: [Outcome: Int] = [:]
        var worst = (angle: 0.0, scale: 0.0)
        for p in phase.cases {
            let m = p.matrix
            let truth = try CrystalOrientation(matrix: simd_double3x3(rows: [SIMD3(m[0], m[1], m[2]), SIMD3(m[3], m[4], m[5]),
                                                                              SIMD3(m[6], m[7], m[8])]))
            let spots = stride(from: 0, to: p.spots.count, by: 2).map { SIMD2(p.spots[$0], p.spots[$0 + 1]) }
            let clicked = p.condition.hasPrefix("click")
            let solution = solver.solve(spots, clicked: clicked)
            // 0.1° on a full pattern at 0.5 px noise, 0.2° with 1 px, 30 % missing and 3 false spots, 0.3° from clicked spots;
            // the scale within 0.05 % at 0.5 px, 0.15 % at 1 px (clicked spots: not bounded).
            let missing = p.condition == "missing"
            let limit = (clicked ? 0.3 : missing ? 0.2 : 0.1) * degree
            let scaleLimit = missing ? 1.5e-3 : 5e-4
            func angle(_ x: ZoneAxisCandidate) -> Double {
                solver.angle(between: x.orientation, and: truth, includingHalfTurn: x.halfTurnAmbiguous)
            }
            let outcome: Outcome
            switch solution.verdict {
            case .refused: outcome = .refused
            case .solved:
                let best = solution.candidates[0]
                let a = angle(best)
                let scaleOff = abs(best.scale / p.scale - 1)
                if a <= limit, clicked || scaleOff <= scaleLimit {
                    outcome = .right
                    worst.angle = max(worst.angle, a)
                    if !clicked { worst.scale = max(worst.scale, scaleOff) }
                } else {
                    outcome = .wrong
                }
            case .tie(let n):
                outcome = solution.candidates.prefix(n).contains { angle($0) <= limit } ? .tie : .wrong
            }
            counts[outcome, default: 0] += 1
            #expect(outcome != .wrong, "\(name) \(p.zone) \(p.condition): \(solution.verdict) \(solution.candidates.prefix(3).map { ($0.zoneAxis, $0.agreement, angle($0) / degree) })")
        }
        print("ZoneAxisSolver planted \(name): right \(counts[.right, default: 0]) tie \(counts[.tie, default: 0]) refused \(counts[.refused, default: 0]) wrong \(counts[.wrong, default: 0]); worst angle \(worst.angle / degree)°, worst scale \(worst.scale)")
        if let bound = Planted.bounds[name] {
            #expect(counts[.right, default: 0] >= bound.right)
            #expect(counts[.refused, default: 0] <= bound.refused)
        }
    }
}
