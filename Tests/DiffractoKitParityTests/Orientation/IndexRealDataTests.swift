//
//  IndexRealDataTests.swift
//  Role: Index on recorded data. `PatternOrientation.index` on 50 patterns of a real 4D-STEM scan (pyxem's demo dataset,
//        Zenodo 14113591, a precession scan of a twinned GaAs nanowire seen down <110>, 300 kV; CC-BY-4.0), read from
//        Tests/DiffractoOrientationTests/Fixtures/index-real-nanowire.json (tools/pins/index-real/nanowire.py: py4DSTEM
//        0.14.19's disks, pyxem-demos' affine correction of the projector distortion, py4DSTEM's pixel calibration against
//        GaAs, then its `orientation_plan` at the defaults and `match_single_pattern(num_matches_return=5)`).
//  Done when (docs/parity.md P60): no truth is known for a recorded pattern, so each pattern is judged against py4DSTEM's
//  answer and against the spots. The kit's best lies within 3° of py4DSTEM's best match (its 2° template step and half a
//  cell either side) or of that match's Σ3 twin (a half turn about a <111>); or, where it lies further, it explains the
//  measured spots at least as well as py4DSTEM's match does (the intensity agreement C of each, py4DSTEM's with its in-plane
//  turn and scale fitted the same way); or the kit flags the pattern. A pattern of fewer than 8 spots may be refused.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoGeometry
import DiffractoOrientation

private struct RealCase: Decodable {
    let position: [Int]
    /// qx, qy, intensity per spot (1/Å about the fitted origin)
    let spots: [Double]
    let py4dstem: [[Double]]
    let py4dstem_corr: [Double]
}
private struct RealFile: Decodable {
    let versions: [String: String]
    let volts: Double
    let cell: [Double]
    let z: [Int]
    let frac: [[Double]]
    let cases: [RealCase]
}

private enum IndexReal {
    static let file: RealFile = {
        let url = FixtureDirs.fixtures(of: "DiffractoOrientationTests").appendingPathComponent("index-real-nanowire.json")
        do { return try JSONDecoder().decode(RealFile.self, from: Data(contentsOf: url)) } catch { fatalError("\(error)") }
    }()

    static func matrix(_ m: [Double]) -> simd_double3x3 {
        simd_double3x3(rows: [SIMD3(m[0], m[1], m[2]), SIMD3(m[3], m[4], m[5]), SIMD3(m[6], m[7], m[8])])
    }

    /// The identity and the half turns about the four <111> (crystal frame): a Σ3 twin of a sphalerite crystal.
    static let twins: [simd_double3x3] = [matrix_identity_double3x3] + [
        SIMD3<Double>(1, 1, 1), SIMD3(1, 1, -1), SIMD3(1, -1, 1), SIMD3(-1, 1, 1),
    ].map { v in
        let n = simd_normalize(v)
        return 2 * simd_double3x3(columns: (n * n.x, n * n.y, n * n.z)) - matrix_identity_double3x3
    }
}

@Suite("Index on recorded data: PatternOrientation on a twinned GaAs nanowire, py4DSTEM's match_single_pattern alongside")
struct IndexRealDataTests {

    @Test func fixtureCoversTheScan() {
        #expect(IndexReal.file.versions["py4DSTEM"] == "0.14.19")
        #expect(IndexReal.file.cases.count == 50)
        #expect(IndexReal.file.cases.allSatisfy { $0.spots.count >= 18 && $0.py4dstem.count == 5 })
    }

    @Test func agreesWithPy4DSTEMOrExplainsTheSpotsBetterOrFlags() throws {
        let f = IndexReal.file, c = f.cell
        let cell = UnitCell(a: c[0], b: c[1], c: c[2], alphaDeg: c[3], betaDeg: c[4], gammaDeg: c[5])
        let sites = zip(f.z, f.frac).map { AtomSite(z: $0.0, fractional: SIMD3($0.1[0], $0.1[1], $0.1[2])) }
        let structure = CrystalStructure.detectingSymmetry(cell: cell, sites: sites, conventional: false)
        let indexer = try PatternOrientation(structure: structure, volts: f.volts)
        let solver = indexer.solver, power = indexer.settings.tilt.intensityPower
        let degree = Double.pi / 180, near = 3 * degree
        var counts = (agree: 0, twin: 0, better: 0, flagged: 0, refused: 0, cases: 0)
        var gaps: [Double] = []
        for p in f.cases {
            counts.cases += 1
            let xyI = stride(from: 0, to: p.spots.count, by: 3)
            let spots = xyI.map { SIMD2(p.spots[$0], p.spots[$0 + 1]) }, intensities = xyI.map { p.spots[$0 + 2] }
            let label = "scan \(p.position) (\(spots.count) spots)"
            guard let r = indexer.index(spots: spots, intensities: intensities) else {
                counts.refused += 1
                #expect(spots.count < 8, "\(label): no answer")
                continue
            }
            let py = try CrystalOrientation(matrix: IndexReal.matrix(p.py4dstem[0]))
            let off = try IndexReal.twins.map { t in
                try solver.angle(between: r.best.orientation, and: CrystalOrientation(matrix: t * py.matrix))
            }
            let gap = off.min()!
            let flagged = r.isFlagged()
            if flagged { counts.flagged += 1 }
            if off[0] <= near { counts.agree += 1; gaps.append(off[0]); continue }
            if gap <= near { counts.twin += 1; continue }
            // further than py4DSTEM's step: the spots decide, py4DSTEM's match fitted (turn, scale) as the kit fits its own
            let fitted = solver.fit(py, spots: spots)
            let pyAgreement = fitted.map {
                solver.intensityAgreement(of: $0.orientation, scale: $0.scale, spots: spots, intensities: intensities,
                                          intensityPower: power)
            } ?? -1
            let better = r.best.agreement >= pyAgreement - 1e-9
            if better { counts.better += 1 }
            #expect(better || flagged,
                    "\(label): \(gap / degree)° from py4DSTEM's best or its twin, agreement \(r.best.agreement) against py4DSTEM's \(pyAgreement), matched \(r.best.positions.matched) against \(fitted?.matched ?? 0), \(r.alternatives.count) alternatives, spread \(r.best.tiltSpread / degree)°")
        }
        gaps.sort()
        let median = gaps.isEmpty ? .nan : gaps[gaps.count / 2] / degree
        print("Index on the nanowire: \(counts.cases) patterns, with py4DSTEM \(counts.agree), its twin \(counts.twin), further and explaining the spots better \(counts.better), flagged \(counts.flagged), refused \(counts.refused); median gap \(median)°")
        #expect(counts.refused <= 1)
    }
}
