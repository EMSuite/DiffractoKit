//
//  PatternOrientationPlantedTests.swift
//  Role: `PatternOrientation.index` (S6) on py4DSTEM-planted patterns with py4DSTEM's own answer, read from
//        Tests/DiffractoOrientationTests/Fixtures/pattern-orientation-planted.json (tools/pins/pattern-orientation/planted.py:
//        py4DSTEM f050d207 `generate_diffraction_pattern` at known orientations, `orientation_plan` at its defaults and
//        `match_single_pattern(num_matches_return=5)`). Al and Si: 13 zone-axis classes to index 3 at 0 to 3° tilt (104 each),
//        the spike's 40 random orientations each, and DiffractoKit #65's Al [001] and [011] cut at the detector edge (24).
//  Done when (plan row S6, docs/parity.md P45): wherever py4DSTEM's answer is within its 2° step of the truth, the kit's is
//  too (the best, or a flagged alternative); no answer more than 2° off goes unflagged; the spike's three wrong random
//  orientations come back right, tied or flagged. The counts and medians per set are the kit's regression bound (measured on
//  the NumPy prototype of the same arithmetic with py4DSTEM's matches as the template step).
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoGeometry
import DiffractoOrientation

private struct PlantedCase: Decodable {
    let set: String
    let zone: [Int]?
    let tilt: Double?
    let matrix: [Double]
    /// qx, qy, intensity per spot
    let spots: [Double]
    let py4dstem: [[Double]]
    let spikeWrong: Bool?

    enum CodingKeys: String, CodingKey { case set, zone, tilt, matrix, spots, py4dstem, spikeWrong = "spike_wrong" }
}
private struct PlantedPhase: Decodable { let name: String; let a: Double; let z: Int; let frac: [[Double]]; let cases: [PlantedCase] }
private struct PlantedFile: Decodable { let versions: [String: String]; let volts: Double; let phases: [PlantedPhase] }

private enum Planted {
    static let file: PlantedFile = {
        let url = FixtureDirs.fixtures(of: "DiffractoOrientationTests").appendingPathComponent("pattern-orientation-planted.json")
        do { return try JSONDecoder().decode(PlantedFile.self, from: Data(contentsOf: url)) } catch { fatalError("\(error)") }
    }()

    /// (phase, set): the least right answers (best within 2°) and the largest median error of the right ones (degrees).
    /// Prototype (right, median): Al zones 86, 0.14; Al edge 23, 0.78; Al random 39, 0.10; Si zones 97, 0.13; Si random 38,
    /// 0.11. The counts leave 3 (2 on the edge set) for the kit's own template step, the medians half again.
    static let bounds: [String: (right: Int, median: Double)] = [
        "Al zones": (83, 0.2), "Al edge": (21, 1.2), "Al random": (36, 0.2), "Si zones": (94, 0.2), "Si random": (35, 0.2),
    ]
    static let sets = ["Al zones", "Al edge", "Al random", "Si zones", "Si random"]

    static func matrix(_ m: [Double]) -> simd_double3x3 {
        simd_double3x3(rows: [SIMD3(m[0], m[1], m[2]), SIMD3(m[3], m[4], m[5]), SIMD3(m[6], m[7], m[8])])
    }
}

@Suite("PatternOrientation on py4DSTEM-planted patterns (py4DSTEM f050d207, its match_single_pattern alongside)")
struct PatternOrientationPlantedTests {

    @Test(arguments: Planted.sets)
    func rightWherePy4DSTEMIsRightAndNeverWrongUnflagged(_ key: String) throws {
        let parts = key.split(separator: " ").map(String.init)
        let phase = try #require(Planted.file.phases.first { $0.name == parts[0] })
        #expect(Planted.file.versions["py4DSTEM"] == "0.14.19")
        let structure = CrystalStructure.detectingSymmetry(
            cell: .cubic(a: phase.a), sites: phase.frac.map { AtomSite(z: phase.z, fractional: SIMD3($0[0], $0[1], $0[2])) },
            conventional: false)
        let indexer = try PatternOrientation(structure: structure, volts: Planted.file.volts)
        let solver = indexer.solver
        let degree = Double.pi / 180, right = 2 * degree
        var errors: [Double] = []
        var counts = (right: 0, tie: 0, flagged: 0, py4DSTEMRight: 0)
        for p in phase.cases where p.set == parts[1] {
            let truth = try CrystalOrientation(matrix: Planted.matrix(p.matrix))
            let xyI = stride(from: 0, to: p.spots.count, by: 3)
            let spots = xyI.map { SIMD2(p.spots[$0], p.spots[$0 + 1]) }, intensities = xyI.map { p.spots[$0 + 2] }
            let label = "\(key) \(p.zone.map { "\($0)" } ?? "random") tilt \(p.tilt ?? .nan)"
            guard let r = indexer.index(spots: spots, intensities: intensities) else {
                Issue.record("\(label): no answer")
                continue
            }
            func off(_ o: CrystalOrientation) -> Double { solver.angle(between: o, and: truth, includingHalfTurn: false) }
            let best = off(r.best.orientation)
            let anyRight = ([r.best] + r.alternatives).contains { off($0.orientation) <= right }
            let flagged = r.isFlagged()
            let py = try off(CrystalOrientation(matrix: Planted.matrix(p.py4dstem[0])))
            if best <= right { counts.right += 1; errors.append(best) } else if anyRight && !r.alternatives.isEmpty { counts.tie += 1 }
            if flagged { counts.flagged += 1 }
            if py <= right { counts.py4DSTEMRight += 1 }
            let message = "\(label): best \(best / degree)°, py4DSTEM \(py / degree)°, \(r.alternatives.count) alternatives, spread \(r.best.tiltSpread / degree)°"
            // within py4DSTEM's 2° step wherever py4DSTEM is right: the best, or a flagged alternative
            if py <= right { #expect(best <= right || (anyRight && flagged), "\(message)") }
            // never wrong without a flag
            #expect(best <= right || flagged, "\(message)")
            if p.spikeWrong == true { #expect(best <= right || anyRight || flagged, "spike's wrong case \(message)") }
        }
        errors.sort()
        let median = errors.isEmpty ? .nan : errors[errors.count / 2] / degree
        print("PatternOrientation planted \(key): right \(counts.right) tie \(counts.tie) flagged \(counts.flagged), py4DSTEM right \(counts.py4DSTEMRight); median \(median)°, worst right \((errors.last ?? .nan) / degree)°")
        let bound = try #require(Planted.bounds[key])
        #expect(counts.right >= bound.right)
        #expect(median <= bound.median)
    }
}
