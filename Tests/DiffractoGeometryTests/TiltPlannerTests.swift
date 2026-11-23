//
//  TiltPlannerTests.swift
//  The stereographic projection, the zone-axis map and the double-tilt holder (S12) on closed forms, round trips and
//  refusals. orix's numbers and the planted calibrations are pins P46 and P47 (DiffractoKitParityTests/Geometry/
//  TiltPlannerPinsTests.swift).
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoGeometry

/// SplitMix64, so the planted round trips are the same on every run.
private struct SeededGenerator: RandomNumberGenerator {
    var state: UInt64
    mutating func next() -> UInt64 {
        state &+= 0x9E3779B97F4A7C15
        var z = state
        z = (z ^ (z >> 30)) &* 0xBF58476D1CE4E5B9
        z = (z ^ (z >> 27)) &* 0x94D049BB133111EB
        return z ^ (z >> 31)
    }
}

/// A uniform random rotation (Shoemake's quaternion).
private func randomRotation(_ rng: inout SeededGenerator) -> simd_double3x3 {
    let u1 = Double.random(in: 0..<1, using: &rng), u2 = Double.random(in: 0..<1, using: &rng)
    let u3 = Double.random(in: 0..<1, using: &rng)
    let a = (1 - u1).squareRoot(), b = u1.squareRoot()
    let q = simd_quatd(ix: a * sin(2 * .pi * u2), iy: a * cos(2 * .pi * u2), iz: b * sin(2 * .pi * u3), r: b * cos(2 * .pi * u3))
    return simd_double3x3(q)
}

private func angleFromBeamDeg(_ v: SIMD3<Double>) -> Double {
    atan2((v.x * v.x + v.y * v.y).squareRoot(), v.z) * 180 / .pi
}

struct TiltPlannerTests {

    /// Closed forms: the beam projects to the centre, the equator to the unit circle, 90° to radius 1 and 2 atan(r) back;
    /// the lower hemisphere is refused; direction(project(v)) = v for 1000 upper-hemisphere directions.
    @Test func projectionClosedForms() throws {
        #expect(StereographicProjection.project(SIMD3(0, 0, 1)) == SIMD2(0, 0))
        #expect(StereographicProjection.project(SIMD3(0, 0, -1)) == nil)
        #expect(StereographicProjection.project(.zero) == nil)
        let equator = try #require(StereographicProjection.project(SIMD3(0.6, 0.8, 0)))
        #expect(abs(length(equator) - 1) <= 1e-15)
        let tilted = try #require(StereographicProjection.project(SIMD3(sin(0.5), 0, cos(0.5))))
        #expect(abs(tilted.x - tan(0.25)) <= 1e-14)
        #expect(abs(StereographicProjection.angleFromCentreDeg(tilted) - 0.5 * 180 / .pi) <= 1e-12)
        var rng = SeededGenerator(state: 7)
        var worst = 0.0
        for _ in 0..<1000 {
            var v = SIMD3(Double.random(in: -1...1, using: &rng), Double.random(in: -1...1, using: &rng),
                          Double.random(in: 0.01...1, using: &rng))
            v = normalize(v)
            let back = StereographicProjection.direction(try #require(StereographicProjection.project(v)))
            worst = max(worst, length(back - v))
        }
        #expect(worst <= 1e-14)
    }

    /// The trace of a plane through the beam is a diameter: the plane x = 0 projects onto X = 0, |Y| ≤ 1; a band edge at
    /// 10° from it is the circle that crosses the X axis at tan(5°).
    @Test func circles() throws {
        let runs = StereographicProjection.circle(normal: SIMD3(1, 0, 0), samples: 720)
        let points = runs.flatMap { $0 }
        #expect(!points.isEmpty && points.allSatisfy { abs($0.x) <= 1e-12 && abs($0.y) <= 1 + 1e-12 })
        let edge = StereographicProjection.circle(normal: SIMD3(1, 0, 0), offsetDeg: 10, samples: 7200).flatMap { $0 }
        let nearAxis = try #require(edge.min { abs($0.y) < abs($1.y) })
        #expect(abs(nearAxis.x - tan(5 * .pi / 180)) <= 1e-3)
    }

    /// Al looking down [001]: [001] is the first pole at 0°, then the four ⟨103⟩ at 18.43° and the ⟨102⟩ at 26.57°; every
    /// pole lies within the asked angle and is projected where its direction says.
    @Test func aluminiumPoles() throws {
        let al = CrystalStructure(cell: .cubic(a: 4.0495), sites: [AtomSite(z: 13, fractional: .zero)], operations: [])
        let o = try CrystalOrientation(zoneAxis: SIMD3(0, 0, 1), cell: al.cell)
        let poles = ZoneAxisMap.poles(structure: al, orientation: o, maxIndex: 3, withinDeg: 30)
        #expect(poles.first?.uvw == SIMD3(0, 0, 1) && (poles.first?.angleFromBeamDeg ?? 1) <= 1e-12)
        let near = poles.dropFirst().prefix(4)
        #expect(near.allSatisfy { abs($0.angleFromBeamDeg - atan(1.0 / 3) * 180 / .pi) <= 1e-9 })
        #expect(near.allSatisfy { Set([abs($0.uvw.x), abs($0.uvw.y), $0.uvw.z]) == Set([0, 1, 3]) })
        #expect(poles.contains { $0.uvw == SIMD3(1, 0, 2) && abs($0.angleFromBeamDeg - atan(0.5) * 180 / .pi) <= 1e-9 })
        #expect(poles.allSatisfy { $0.angleFromBeamDeg <= 30 && abs(StereographicProjection.angleFromCentreDeg($0.position) - $0.angleFromBeamDeg) <= 1e-9 })
    }

    /// The traces of fcc Al's three lowest families ({111}, {200}, {220}): one per ±pair (4 + 3 + 6), each normal a unit
    /// vector, the Bragg angle asin(λ/2d).
    @Test func aluminiumTraces() throws {
        let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        let al = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0495), sites: fcc.map { AtomSite(z: 13, fractional: $0) })
        let o = try CrystalOrientation(zoneAxis: SIMD3(0, 0, 1), cell: al.cell)
        let lambda = ElectronWavelength.angstroms(volts: 200_000)
        let traces = ZoneAxisMap.planeTraces(structure: al, orientation: o, families: 3, wavelengthAngstroms: lambda)
        #expect(traces.count == 13)
        #expect(traces.allSatisfy { abs(length($0.normal) - 1) <= 1e-12 })
        let first = try #require(traces.first)
        #expect(abs(first.d - 4.0495 / 3.0.squareRoot()) <= 1e-9)
        #expect(abs((first.braggAngleDeg ?? 0) - asin(lambda / (2 * first.d)) * 180 / .pi) <= 1e-12)
    }

    /// 20 000 planted cases: random crystal, axis angle, current tilts within ±30° and target direction; after moving to the
    /// computed tilts the direction is on the beam within 1e-5°, and moving back restores the orientation.
    /// RED UNDER: α and β swapped; φ ignored; the move composed on the wrong side.
    @Test func roundTrips() throws {
        var rng = SeededGenerator(state: 20_261_124)
        var worstBeam = 0.0, worstBack = 0.0
        for _ in 0..<20_000 {
            let holder = DoubleTiltHolder(axisAngleDeg: Double.random(in: 0..<360, using: &rng))
            let o = try CrystalOrientation(matrix: randomRotation(&rng))
            let current = DoubleTiltHolder.Tilts(alpha: Double.random(in: -30...30, using: &rng),
                                                 beta: Double.random(in: -30...30, using: &rng))
            let t = normalize(SIMD3(Double.random(in: -1...1, using: &rng), Double.random(in: -1...1, using: &rng),
                                    Double.random(in: -1...1, using: &rng)))
            let target = try #require(holder.tilts(toBring: t, orientation: o, at: current))
            let moved = try #require(holder.orientation(o, at: current, movedTo: target.tilts))
            worstBeam = max(worstBeam, angleFromBeamDeg(moved.labCoordinates(of: t)))
            let back = try #require(holder.orientation(moved, at: target.tilts, movedTo: current))
            let e = back.matrix - o.matrix
            worstBack = max(worstBack, [e.columns.0, e.columns.1, e.columns.2].map { simd_reduce_max(simd_abs($0)) }.max() ?? 1)
        }
        #expect(worstBeam <= 1e-5, "worst \(worstBeam)°")
        #expect(worstBack <= 1e-12)
    }

    /// Reach: a direction within the limits is reachable, one needing α beyond them is not; the angle from the beam is
    /// reported whatever the holder. Without a calibrated axis there are no tilts.
    @Test func reachAndUncalibrated() throws {
        let o = try CrystalOrientation(matrix: matrix_identity_double3x3)
        let holder = DoubleTiltHolder(axisAngleDeg: 0)
        let near = try #require(holder.tilts(toBring: SIMD3(0, sin(0.3), cos(0.3)), orientation: o, at: .zero))
        #expect(near.reachable && abs(near.tilts.alpha - 0.3 * 180 / .pi) <= 1e-12 && abs(near.tilts.beta) <= 1e-12)
        #expect(abs(near.angleFromBeamDeg - 0.3 * 180 / .pi) <= 1e-12)
        let far = try #require(holder.tilts(toBring: SIMD3(0, sin(0.8), cos(0.8)), orientation: o, at: .zero))
        #expect(!far.reachable)
        #expect(DoubleTiltHolder().tilts(toBring: SIMD3(0, 0, 1), orientation: o, at: .zero) == nil)
        #expect(DoubleTiltHolder().orientation(o, at: .zero, movedTo: .zero) == nil)
    }

    /// Calibration refusals: a tilt change under 2°, a pair that fits no axis (the second pattern 40° off the stage), a
    /// tilt that is not finite. A clean planted pair gives the axis to 1e-6° and no half turn.
    @Test func calibration() throws {
        var rng = SeededGenerator(state: 99)
        let s = randomRotation(&rng)
        let phi = 37.0
        let holder = DoubleTiltHolder(axisAngleDeg: phi)
        let a = try CrystalOrientation(matrix: s)
        let ta = DoubleTiltHolder.Tilts(alpha: -5, beta: 3), tb = DoubleTiltHolder.Tilts(alpha: 12, beta: -8)
        let b = try #require(holder.orientation(a, at: ta, movedTo: tb))
        let cal = try DoubleTiltHolder.calibrateAxis(a, at: ta, b, at: tb)
        #expect(abs(cal.axisAngleDeg - phi) <= 1e-6 && !cal.secondNeedsHalfTurn && cal.misfitDeg <= 1e-6)
        let flipped = b.rotatedInPlane(by: .pi)
        #expect(try DoubleTiltHolder.calibrateAxis(a, at: ta, flipped, at: tb).secondNeedsHalfTurn)
        #expect(throws: DoubleTiltHolder.CalibrationError.self) {
            try DoubleTiltHolder.calibrateAxis(a, at: ta, b, at: DoubleTiltHolder.Tilts(alpha: -4, beta: 3))
        }
        let other = b.tilted(by: SIMD2(0.7, 0))      // 40° more than the stage says
        #expect(throws: DoubleTiltHolder.CalibrationError.self) { try DoubleTiltHolder.calibrateAxis(a, at: ta, other, at: tb) }
        #expect(throws: DoubleTiltHolder.CalibrationError.notFinite) {
            try DoubleTiltHolder.calibrateAxis(a, at: ta, b, at: DoubleTiltHolder.Tilts(alpha: .nan, beta: 0))
        }
    }

    /// The proper rotations of a cubic crystal are its 24 rotations, orthogonal in the Cartesian frame.
    @Test func cubicRotations() {
        let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        let al = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0495), sites: fcc.map { AtomSite(z: 13, fractional: $0) })
        let ops = DoubleTiltHolder.properRotations(of: al)
        #expect(ops.count == 24)
        #expect(ops.allSatisfy { abs($0.determinant - 1) <= 1e-12 })
    }
}
