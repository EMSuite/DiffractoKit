//
//  TiltPlannerPinsTests.swift
//  Role: the tilt planner's kit pieces (S12) against their references:
//    orix-stereographic-pins.json  tools/pins/orix/stereographic-pins.py: orix 0.15.0 StereographicProjection() and its inverse
//                                  on 1010 directions and 1003 points (pin P46), np.allclose defaults (rtol 1e-5, atol 1e-8),
//                                  which hemisphere exact;
//    tilt-planner-planted.json     tools/pins/tilt-planner/planted.py: planted holder tilts, moves and two-pattern
//                                  calibrations in NumPy (pin P47). Tilts and moves at 1e-9 (a closed form); the calibrated
//                                  axis within 1.5° at the 95th percentile for tilt changes from 10° with 0.2° error per
//                                  pattern, the half turn of the second pattern named right in every trial.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoGeometry

private struct StereoForward: Decodable { let v: [Double]; let kept: Bool; let xy: [Double]? }
private struct StereoInverse: Decodable { let xy: [Double]; let v: [Double] }
private struct StereoPins: Decodable { let forward: [StereoForward]; let inverse: [StereoInverse] }

private struct TiltCase: Decodable {
    let M0: [[Double]]
    let phiDeg: Double
    let currentDeg: [Double]
    let direction: [Double]
    let targetDeg: [Double]
    let M1: [[Double]]
}
private struct CalibrationCase: Decodable {
    let bandDeg: [Double]
    let tiltChangeDeg: Double
    let phiDeg: Double
    let tiltsADeg: [Double]
    let tiltsBDeg: [Double]
    let Ma: [[Double]]
    let Mb: [[Double]]
    let bHalfTurned: Bool
}
private struct PlantedPins: Decodable { let tilts: [TiltCase]; let calibration: [CalibrationCase] }

private enum TiltPins {
    static let dir = FixtureDirs.fixtures(of: "DiffractoGeometryTests")
    static let stereo: Result<StereoPins, any Error> = Result {
        try JSONDecoder().decode(StereoPins.self, from: Data(contentsOf: dir.appendingPathComponent("orix-stereographic-pins.json")))
    }
    static let planted: Result<PlantedPins, any Error> = Result {
        let d = JSONDecoder()
        d.keyDecodingStrategy = .convertFromSnakeCase
        return try d.decode(PlantedPins.self, from: Data(contentsOf: dir.appendingPathComponent("tilt-planner-planted.json")))
    }
    static func matrix(_ rows: [[Double]]) -> simd_double3x3 {
        simd_double3x3(rows: rows.map { SIMD3($0[0], $0[1], $0[2]) })
    }
}

private func allclose(_ a: Double, _ b: Double) -> Bool { abs(a - b) <= 1e-8 + 1e-5 * abs(b) }

private func maxDifference(_ a: simd_double3x3, _ b: simd_double3x3) -> Double {
    let e = a - b
    return [e.columns.0, e.columns.1, e.columns.2].map { simd_reduce_max(simd_abs($0)) }.max() ?? .infinity
}

struct TiltPlannerPinsTests {

    /// orix's forward projection: the same vectors kept (z as given above −1e-9) and (X, Y) within np.allclose.
    /// RED UNDER: (x, y)/(1 − z) (the other pole); the unit vector not taken; the hemisphere cut on the unit z or at 0.
    @Test func forwardMatchesOrix() throws {
        let pins = try TiltPins.stereo.get()
        #expect(pins.forward.count == 1010)
        var bad: [String] = []
        for f in pins.forward {
            let p = StereographicProjection.project(SIMD3(f.v[0], f.v[1], f.v[2]))
            switch (p, f.xy) {
            case (nil, nil): break
            case let (got?, want?):
                if !allclose(got.x, want[0]) || !allclose(got.y, want[1]) { bad.append("\(f.v): \(got) vs \(want)") }
            default: bad.append("\(f.v): kept \(p != nil), orix \(f.kept)")
            }
        }
        #expect(bad.isEmpty, "\(bad.count) differ: \(bad.prefix(3))")
    }

    /// orix's inverse projection within np.allclose, inside and outside the unit circle.
    @Test func inverseMatchesOrix() throws {
        let pins = try TiltPins.stereo.get()
        var bad = 0
        for i in pins.inverse {
            let v = StereographicProjection.direction(SIMD2(i.xy[0], i.xy[1]))
            if !(allclose(v.x, i.v[0]) && allclose(v.y, i.v[1]) && allclose(v.z, i.v[2])) { bad += 1 }
        }
        #expect(pins.inverse.count == 1003 && bad == 0)
    }

    /// The planted tilts and moves: the tilts that bring each direction onto the beam equal the closed form (1e-9°), and the
    /// orientation after the move equals the planted one (1e-9).
    /// RED UNDER: Ry Rx for Rx Ry (β about a fixed axis); φ ignored or of the other sign; β = atan2(t_x, t_z).
    @Test func plantedTiltsAndMoves() throws {
        let pins = try TiltPins.planted.get()
        #expect(pins.tilts.count == 500)
        var worstTilt = 0.0, worstMove = 0.0
        for c in pins.tilts {
            let holder = DoubleTiltHolder(axisAngleDeg: c.phiDeg)
            let m0 = try CrystalOrientation(matrix: TiltPins.matrix(c.M0))
            let current = DoubleTiltHolder.Tilts(alpha: c.currentDeg[0], beta: c.currentDeg[1])
            let target = try #require(holder.tilts(toBring: SIMD3(c.direction[0], c.direction[1], c.direction[2]),
                                                   orientation: m0, at: current))
            worstTilt = max(worstTilt, abs(target.tilts.alpha - c.targetDeg[0]), abs(target.tilts.beta - c.targetDeg[1]))
            let moved = try #require(holder.orientation(m0, at: current, movedTo: target.tilts))
            worstMove = max(worstMove, maxDifference(moved.matrix, TiltPins.matrix(c.M1)))
        }
        #expect(worstTilt <= 1e-9, "tilts differ by \(worstTilt)°")
        #expect(worstMove <= 1e-9, "moved orientation differs by \(worstMove)")
    }

    /// The planted calibrations: per band of tilt change, the axis within 1.5° at the 95th percentile, the half turn of the
    /// second pattern named right in every trial, and no trial refused.
    /// RED UNDER: the half turn not tried; Q composed in the wrong order; the refinement skipped (the scan's 0.5° steps alone
    /// still pass p95, so the median is held too: below 0.6°).
    @Test func plantedCalibrations() throws {
        let pins = try TiltPins.planted.get()
        #expect(pins.calibration.count == 600)
        for band in [10.0, 20.0] {
            let cases = pins.calibration.filter { $0.bandDeg[0] == band }
            var errors: [Double] = []
            var wrongHalf = 0, refused = 0
            for c in cases {
                let a = try CrystalOrientation(matrix: TiltPins.matrix(c.Ma))
                let b = try CrystalOrientation(matrix: TiltPins.matrix(c.Mb))
                do {
                    let cal = try DoubleTiltHolder.calibrateAxis(
                        a, at: .init(alpha: c.tiltsADeg[0], beta: c.tiltsADeg[1]),
                        b, at: .init(alpha: c.tiltsBDeg[0], beta: c.tiltsBDeg[1]))
                    var e = (cal.axisAngleDeg - c.phiDeg).truncatingRemainder(dividingBy: 360)
                    if e > 180 { e -= 360 }
                    if e < -180 { e += 360 }
                    errors.append(abs(e))
                    if cal.secondNeedsHalfTurn != c.bHalfTurned { wrongHalf += 1 }
                    #expect(abs(cal.tiltChangeDeg - c.tiltChangeDeg) <= 1e-6)
                } catch {
                    refused += 1
                }
            }
            errors.sort()
            let p95 = errors.isEmpty ? .infinity : errors[min(errors.count - 1, Int(0.95 * Double(errors.count)))]
            let median = errors.isEmpty ? .infinity : errors[errors.count / 2]
            #expect(cases.count == 300 && refused == 0 && wrongHalf == 0, "band \(band): refused \(refused), wrong half turn \(wrongHalf)")
            #expect(p95 <= 1.5 && median <= 0.6, "band \(band): median \(median)°, p95 \(p95)°")
        }
    }
}
