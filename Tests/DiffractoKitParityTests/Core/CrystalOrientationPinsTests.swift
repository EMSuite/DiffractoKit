//
//  CrystalOrientationPinsTests.swift
//  Role: `CrystalOrientation` (S1) against its references, read from Tests/DiffractoCoreTests/Fixtures:
//    crystal-orientation-pins.json (tools/pins/crystal-orientation/pins-py4dstem.py, py4DSTEM f050d207):
//      the zone-axis matrix of `parse_orientation` for 7 phases and 37 zone axes, the matrix at three in-plane angles
//      (`@ m3z`), and where every reflection of `generate_diffraction_pattern` appears for those matrices and 20 random ones;
//    crystal-orientation-orix-pins.json (tools/pins/crystal-orientation/pins-orix.py, orix 0.15.0):
//      `Orientation.from_matrix(M).data` and `Orientation(q).to_matrix()`.
//  py4DSTEM ships no tests: numpy.testing.assert_allclose's float64 default (rtol 1e-7) applies, with an absolute floor of
//  1e-12 for entries that are exact zeros up to rounding. orix: np.allclose defaults (rtol 1e-5, atol 1e-8), its own tests'.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore

private struct PinCell: Decodable { let a, b, c, alpha, beta, gamma: Double }
private struct PinInPlane: Decodable { let gamma: Double; let matrix: [Double]; let peaks: [[Double]] }
private struct PinAxis: Decodable { let uvw: [Int]; let matrix: [Double]; let inPlane: [PinInPlane] }
private struct PinPhase: Decodable { let name: String; let cell: PinCell; let latReal: [[Double]]; let axes: [PinAxis] }
private struct PinRandom: Decodable { let phase: String; let matrix: [Double]; let peaks: [[Double]] }
private struct Py4DSTEMPins: Decodable {
    let versions: [String: String]
    let commit: String
    let phases: [PinPhase]
    let random: [PinRandom]
}
private struct OrixRow: Decodable { let matrix: [Double]; let q: [Double] }
private struct OrixPins: Decodable { let versions: [String: String]; let fromMatrix: [OrixRow]; let toMatrix: [OrixRow] }

private enum Pins {
    static func load<T: Decodable>(_ name: String) -> T {
        let url = FixtureDirs.fixtures(of: "DiffractoCoreTests").appendingPathComponent(name)
        let decoder = JSONDecoder()
        decoder.keyDecodingStrategy = .convertFromSnakeCase
        do { return try decoder.decode(T.self, from: Data(contentsOf: url)) } catch { fatalError("\(name): \(error)") }
    }
    static let py4dstem: Py4DSTEMPins = load("crystal-orientation-pins.json")
    static let orix: OrixPins = load("crystal-orientation-orix-pins.json")
}

/// Row-major 9 numbers to a matrix.
private func matrix(_ r: [Double]) -> simd_double3x3 {
    simd_double3x3(rows: [SIMD3(r[0], r[1], r[2]), SIMD3(r[3], r[4], r[5]), SIMD3(r[6], r[7], r[8])])
}

private func cell(_ c: PinCell) -> UnitCell {
    UnitCell(a: c.a, b: c.b, c: c.c, alphaDeg: c.alpha, betaDeg: c.beta, gammaDeg: c.gamma)
}

/// The largest |got − want| − rtol |want| over the nine entries (≤ floor passes).
private func excess(_ got: simd_double3x3, _ want: simd_double3x3, rtol: Double) -> Double {
    var worst = -Double.infinity
    for c in 0..<3 { for r in 0..<3 { worst = max(worst, abs(got[c][r] - want[c][r]) - rtol * abs(want[c][r])) } }
    return worst
}

@Suite("CrystalOrientation against py4DSTEM f050d207 and orix 0.15.0")
struct CrystalOrientationPinsTests {

    @Test func pinsCarryTheirVersionsAndSizes() {
        #expect(Pins.py4dstem.versions["py4DSTEM"] == "0.14.19")
        #expect(Pins.py4dstem.commit.hasPrefix("f050d207"))
        #expect(Pins.py4dstem.phases.count == 7)
        #expect(Pins.py4dstem.phases.reduce(0) { $0 + $1.axes.count } == 37)
        #expect(Pins.py4dstem.random.count == 20)
        #expect(Pins.orix.versions["orix"] == "0.15.0")
        #expect(Pins.orix.fromMatrix.count == 46 && Pins.orix.toMatrix.count == 20)
    }

    /// The cell's real lattice is py4DSTEM's `lat_real` (the frame every number below is written in).
    @Test func theCellFrameIsPy4DSTEMs() {
        for p in Pins.py4dstem.phases {
            let c = cell(p.cell)
            for i in 0..<3 { for j in 0..<3 {
                #expect(abs(c.realLattice[i][j] - p.latReal[i][j]) <= 1e-12 * c.a + 1e-7 * abs(p.latReal[i][j]), "\(p.name) [\(i)][\(j)]")
            } }
        }
    }

    /// 7 phases × 37 zone axes: the orientation of [u v w] at in-plane angle 0 is `parse_orientation(zone_axis_lattice=uvw)`,
    /// and at 0.7 and 2.5 rad it is that matrix `@ m3z(gamma)`; the third column is the zone axis.
    @Test func zoneAxisMatricesEqualPy4DSTEMs() throws {
        var worst = -Double.infinity
        for p in Pins.py4dstem.phases {
            let c = cell(p.cell)
            for a in p.axes {
                let uvw = SIMD3(a.uvw[0], a.uvw[1], a.uvw[2])
                let o = try CrystalOrientation(zoneAxis: uvw, cell: c)
                let e = excess(o.matrix, matrix(a.matrix), rtol: 1e-7)
                worst = max(worst, e)
                #expect(e <= 1e-12, "\(p.name) \(a.uvw): \(e)")
                let direction = Double(uvw.x) * c.realLattice[0] + Double(uvw.y) * c.realLattice[1] + Double(uvw.z) * c.realLattice[2]
                #expect(simd_length(o.beamDirection - simd_normalize(direction)) <= 1e-12, "\(p.name) \(a.uvw): third column")
                for row in a.inPlane {
                    let turned = try CrystalOrientation(zoneAxis: uvw, inPlaneAngle: row.gamma, cell: c)
                    let et = excess(turned.matrix, matrix(row.matrix), rtol: 1e-7)
                    worst = max(worst, et)
                    #expect(et <= 1e-12, "\(p.name) \(a.uvw) at \(row.gamma): \(et)")
                }
            }
        }
        print("S1WORST zone-axis matrices \(worst)")
    }

    /// Where g appears: (Mᵀg).x, (Mᵀg).y equals py4DSTEM's (qx, qy) for every reflection of `generate_diffraction_pattern` at
    /// the 111 in-plane matrices and the 20 random rotations.
    @Test func reflectionsAppearWherePy4DSTEMPutsThem() throws {
        var worst = 0.0, count = 0
        func check(_ o: CrystalOrientation, _ c: UnitCell, _ peaks: [[Double]], _ label: String) {
            for pk in peaks {
                let g = c.gVector(MillerIndex(Int(pk[0]), Int(pk[1]), Int(pk[2])))
                let p = o.labCoordinates(of: g)
                let d = max(abs(p.x - pk[3]) - 1e-7 * abs(pk[3]), abs(p.y - pk[4]) - 1e-7 * abs(pk[4]))
                worst = max(worst, abs(p.x - pk[3]), abs(p.y - pk[4]))
                count += 1
                #expect(d <= 1e-12, "\(label) \(pk.prefix(3)): (\(p.x), \(p.y)) against (\(pk[3]), \(pk[4]))")
            }
        }
        var cells: [String: UnitCell] = [:]
        for p in Pins.py4dstem.phases {
            let c = cell(p.cell)
            cells[p.name] = c
            for a in p.axes {
                for row in a.inPlane {
                    let o = try CrystalOrientation(zoneAxis: SIMD3(a.uvw[0], a.uvw[1], a.uvw[2]), inPlaneAngle: row.gamma, cell: c)
                    check(o, c, row.peaks, "\(p.name) \(a.uvw) \(row.gamma)")
                }
            }
        }
        for r in Pins.py4dstem.random {
            let c = try #require(cells[r.phase])
            check(try CrystalOrientation(matrix: matrix(r.matrix)), c, r.peaks, "random \(r.phase)")
        }
        #expect(count > 1000)
        print("S1WORST pattern positions \(worst) over \(count) reflections")
    }

    /// orix `Orientation.from_matrix(M).data`: the same quaternion up to its sign, with a ≥ 0.
    @Test func quaternionsEqualOrixFromMatrix() throws {
        for row in Pins.orix.fromMatrix {
            let q = try CrystalOrientation(matrix: matrix(row.matrix)).orixQuaternion
            let want = SIMD4(row.q[0], row.q[1], row.q[2], row.q[3])
            let plus = simd_reduce_max(abs(q - want) - 1e-5 * abs(want)), minus = simd_reduce_max(abs(q + want) - 1e-5 * abs(want))
            #expect(min(plus, minus) <= 1e-8, "\(row.q) against \(q)")
            #expect(q.x >= 0)
            if want.x > 1e-6 { #expect(plus <= 1e-8, "the sign: \(row.q) against \(q)") }
        }
    }

    /// orix `Orientation(q).to_matrix()`: the same matrix, for unit quaternions of either sign of a.
    @Test func matricesEqualOrixToMatrix() throws {
        for row in Pins.orix.toMatrix {
            let o = try CrystalOrientation(orixQuaternion: SIMD4(row.q[0], row.q[1], row.q[2], row.q[3]))
            #expect(excess(o.matrix, matrix(row.matrix), rtol: 1e-5) <= 1e-8, "\(row.q)")
        }
    }
}
