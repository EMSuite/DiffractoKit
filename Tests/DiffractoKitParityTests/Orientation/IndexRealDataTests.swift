//
//  IndexRealDataTests.swift
//  Role: Index on recorded data. `PatternOrientation.index` on 50 patterns of a real 4D-STEM scan (pyxem's demo dataset,
//        Zenodo 14113591, a precession scan of a twinned GaAs nanowire seen down <110>, 300 kV; CC-BY-4.0), read from
//        Tests/DiffractoOrientationTests/Fixtures/index-real-nanowire.json (tools/pins/index-real/nanowire.py: py4DSTEM
//        0.14.19's disks, pyxem-demos' affine correction of the projector distortion, py4DSTEM's pixel calibration against
//        GaAs, then its `orientation_plan` at the defaults and `match_single_pattern(num_matches_return=5)`).
//  Done when (docs/parity.md P60): on every pattern the kit's best lies within 3° of py4DSTEM's best match (its 2° template
//  step and the half of a grid cell either side) or of that match's twin (the Σ3 twin of the wire: a half turn about a
//  <111>), or the kit flags the pattern. No truth is known for a recorded pattern; py4DSTEM's answer is the reference.
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

    @Test func agreesWithPy4DSTEMOrItsTwinOrFlags() throws {
        let f = IndexReal.file, c = f.cell
        let cell = UnitCell(a: c[0], b: c[1], c: c[2], alphaDeg: c[3], betaDeg: c[4], gammaDeg: c[5])
        let sites = zip(f.z, f.frac).map { AtomSite(z: $0.0, fractional: SIMD3($0.1[0], $0.1[1], $0.1[2])) }
        let structure = CrystalStructure.detectingSymmetry(cell: cell, sites: sites, conventional: false)
        let indexer = try PatternOrientation(structure: structure, volts: f.volts)
        let solver = indexer.solver
        let degree = Double.pi / 180, near = 3 * degree
        var counts = (agree: 0, twin: 0, flagged: 0, cases: 0)
        var gaps: [Double] = []
        for p in f.cases {
            counts.cases += 1
            let xyI = stride(from: 0, to: p.spots.count, by: 3)
            let spots = xyI.map { SIMD2(p.spots[$0], p.spots[$0 + 1]) }, intensities = xyI.map { p.spots[$0 + 2] }
            let label = "scan \(p.position)"
            guard let r = indexer.index(spots: spots, intensities: intensities) else {
                Issue.record("\(label): no answer")
                continue
            }
            let py = IndexReal.matrix(p.py4dstem[0])
            let off = try IndexReal.twins.map { t in
                try solver.angle(between: r.best.orientation, and: CrystalOrientation(matrix: t * py))
            }
            let gap = off.min()!
            let flagged = r.isFlagged()
            if off[0] <= near { counts.agree += 1; gaps.append(off[0]) } else if gap <= near { counts.twin += 1 }
            if flagged { counts.flagged += 1 }
            #expect(gap <= near || flagged,
                    "\(label): \(off[0] / degree)° from py4DSTEM's best, \(gap / degree)° from it or its twin, \(r.alternatives.count) alternatives, spread \(r.best.tiltSpread / degree)°")
        }
        gaps.sort()
        let median = gaps.isEmpty ? .nan : gaps[gaps.count / 2] / degree
        print("Index on the nanowire: \(counts.cases) patterns, with py4DSTEM \(counts.agree), its twin \(counts.twin), flagged \(counts.flagged); median gap \(median)°")
    }
}
