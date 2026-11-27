//
//  IndexMatrixTests.swift
//  Role: the Index reference matrix. `PatternOrientation.index` on py4DSTEM-planted spot patterns of eight crystals of
//        five systems (fcc Au, rock-salt MgO, bcc Fe, perovskite SrTiO3, hcp Mg and Zr, rutile TiO2, 2H-WS2), read from
//        Tests/DiffractoOrientationTests/Fixtures/index-matrix-planted.json (tools/pins/index-matrix/planted.py: py4DSTEM
//        f050d207 `generate_diffraction_pattern` at known orientations, `orientation_plan` at its defaults and
//        `match_single_pattern(num_matches_return=5)`). Per crystal 8 zone-axis classes at 0, 1 and 2 degrees of tilt (48
//        cases) and 10 random orientations.
//  Done when (the rule of PatternOrientationPlantedTests, docs/parity.md P45): wherever py4DSTEM's answer is within its 2°
//  step of the truth, the kit's is too (the best, or a flagged alternative); no answer more than 2° off goes unflagged.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoGeometry
import DiffractoOrientation

private struct MatrixCase: Decodable {
    let set: String
    let zone: [Int]?
    let tilt: Double?
    let matrix: [Double]
    /// qx, qy, intensity per spot
    let spots: [Double]
    let py4dstem: [[Double]]
}
private struct MatrixPhase: Decodable { let name: String; let cell: [Double]; let z: [Int]; let frac: [[Double]]; let cases: [MatrixCase] }
private struct MatrixFile: Decodable { let versions: [String: String]; let volts: Double; let phases: [MatrixPhase] }

private enum IndexMatrix {
    static let file: MatrixFile = {
        let url = FixtureDirs.fixtures(of: "DiffractoOrientationTests").appendingPathComponent("index-matrix-planted.json")
        do { return try JSONDecoder().decode(MatrixFile.self, from: Data(contentsOf: url)) } catch { fatalError("\(error)") }
    }()
    static let phases = ["Au", "MgO", "Fe", "SrTiO3", "Mg", "Zr", "TiO2", "WS2"]

    static func matrix(_ m: [Double]) -> simd_double3x3 {
        simd_double3x3(rows: [SIMD3(m[0], m[1], m[2]), SIMD3(m[3], m[4], m[5]), SIMD3(m[6], m[7], m[8])])
    }
}

@Suite("Index matrix: PatternOrientation on eight crystals planted by py4DSTEM (f050d207), its match_single_pattern alongside")
struct IndexMatrixTests {

    @Test func fixtureCoversTheMatrix() {
        #expect(IndexMatrix.file.versions["py4DSTEM"] == "0.14.19")
        #expect(IndexMatrix.file.phases.map(\.name) == IndexMatrix.phases)
        for p in IndexMatrix.file.phases {
            #expect(p.cases.filter { $0.set == "zones" }.count == 48, "\(p.name) zones")
            #expect(p.cases.filter { $0.set == "random" }.count == 10, "\(p.name) random")
        }
    }

    @Test(arguments: IndexMatrix.phases)
    func rightWherePy4DSTEMIsRightAndNeverWrongUnflagged(_ name: String) throws {
        let phase = try #require(IndexMatrix.file.phases.first { $0.name == name })
        let c = phase.cell
        let cell = UnitCell(a: c[0], b: c[1], c: c[2], alphaDeg: c[3], betaDeg: c[4], gammaDeg: c[5])
        let sites = zip(phase.z, phase.frac).map { AtomSite(z: $0.0, fractional: SIMD3($0.1[0], $0.1[1], $0.1[2])) }
        let structure = CrystalStructure.detectingSymmetry(cell: cell, sites: sites, conventional: false)
        let indexer = try PatternOrientation(structure: structure, volts: IndexMatrix.file.volts)
        let solver = indexer.solver
        let degree = Double.pi / 180, right = 2 * degree
        var errors: [Double] = []
        var counts = (right: 0, flagged: 0, py4DSTEMRight: 0, cases: 0)
        for p in phase.cases {
            counts.cases += 1
            let truth = try CrystalOrientation(matrix: IndexMatrix.matrix(p.matrix))
            let xyI = stride(from: 0, to: p.spots.count, by: 3)
            let spots = xyI.map { SIMD2(p.spots[$0], p.spots[$0 + 1]) }, intensities = xyI.map { p.spots[$0 + 2] }
            let label = "\(name) \(p.zone.map { "\($0)" } ?? "random") tilt \(p.tilt ?? .nan)"
            guard let r = indexer.index(spots: spots, intensities: intensities) else {
                Issue.record("\(label): no answer")
                continue
            }
            func off(_ o: CrystalOrientation) -> Double { solver.angle(between: o, and: truth, includingHalfTurn: false) }
            let best = off(r.best.orientation)
            let anyRight = ([r.best] + r.alternatives).contains { off($0.orientation) <= right }
            let flagged = r.isFlagged()
            let py = try off(CrystalOrientation(matrix: IndexMatrix.matrix(p.py4dstem[0])))
            if best <= right { counts.right += 1; errors.append(best) }
            if flagged { counts.flagged += 1 }
            if py <= right { counts.py4DSTEMRight += 1 }
            let message = "\(label): best \(best / degree)°, py4DSTEM \(py / degree)°, \(r.alternatives.count) alternatives, spread \(r.best.tiltSpread / degree)°"
            if py <= right { #expect(best <= right || (anyRight && flagged), "\(message)") }
            #expect(best <= right || flagged, "\(message)")
        }
        errors.sort()
        let median = errors.isEmpty ? .nan : errors[errors.count / 2] / degree
        print("Index matrix \(name): \(counts.cases) cases, right \(counts.right) flagged \(counts.flagged), py4DSTEM right \(counts.py4DSTEMRight); median \(median)°, worst right \((errors.last ?? .nan) / degree)°")
        // the kit is right at least as often as py4DSTEM's 2° template match
        #expect(counts.right >= counts.py4DSTEMRight - counts.flagged)
    }
}
