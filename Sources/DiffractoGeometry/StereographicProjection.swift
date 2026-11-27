//
//  StereographicProjection.swift
//  Role: the stereographic projection of the upper hemisphere (z toward the electron source) onto the plane, orix 0.15.0's
//        `StereographicProjection()` (pole −1, the projection point at [0 0 −1]) and its inverse (pin P50):
//          kept when z (of the vector as given) > −1e-9; (X, Y) = (x, y) / (z + 1) of the unit vector;
//          inverse (x, y, z) = (2X, 2Y, 1 − X² − Y²) / (1 + X² + Y²).
//        In the pattern frame of a `CrystalOrientation` the beam (+z) is the centre and the primitive circle (radius 1) is
//        90° from it. Also the arcs of great and small circles the zone-axis map draws.
//

import Foundation
import simd

/// The stereographic projection of the upper hemisphere (orix's default projection): the static members. A value is orix's
/// `StereographicProjection(pole)` with either pole (the members after `init(pole:)`, from Electron Microscopy Studio's
/// orix port, DiffractoKit#179).
public struct StereographicProjection: Sendable, Equatable {

    /// orix's region cut: a vector is in the upper hemisphere when its z, as given, is above this.
    public static let hemisphereTolerance = -1e-9

    /// (X, Y) of the direction `v` (any length but zero), nil when it lies in the lower hemisphere (z ≤ −1e-9) or is not a
    /// finite non-zero vector.
    public static func project(_ v: SIMD3<Double>) -> SIMD2<Double>? {
        let n = length(v)
        guard n.isFinite, n > 0, v.z > hemisphereTolerance else { return nil }
        let u = v / n
        let denom = u.z + 1
        guard denom != 0 else { return SIMD2(0, 0) }
        return SIMD2(u.x / denom, u.y / denom)
    }

    /// The unit vector at (X, Y): inside the unit circle the upper hemisphere, outside it the lower.
    public static func direction(_ p: SIMD2<Double>) -> SIMD3<Double> {
        let r2 = p.x * p.x + p.y * p.y
        let denom = 1 + r2
        return SIMD3(2 * p.x / denom, 2 * p.y / denom, (1 - r2) / denom)
    }

    /// The angle (degrees) between the beam (+z) and the direction projected at (X, Y): 2 atan(r).
    public static func angleFromCentreDeg(_ p: SIMD2<Double>) -> Double {
        2 * atan((p.x * p.x + p.y * p.y).squareRoot()) * 180 / .pi
    }

    /// The points of the circle at `offsetDeg` from the great circle whose pole is `normal` (0: the great circle, the trace of
    /// the plane with that normal; ±θ: the two edges of a band of half-width θ), in the upper hemisphere, as runs of
    /// consecutive projected points (a run breaks where the circle leaves the hemisphere). `samples` points around the whole
    /// circle.
    public static func circle(normal: SIMD3<Double>, offsetDeg: Double = 0, samples: Int = 360) -> [[SIMD2<Double>]] {
        let n0 = length(normal)
        guard n0.isFinite, n0 > 0, samples >= 3, offsetDeg.isFinite, abs(offsetDeg) < 90 else { return [] }
        let n = normal / n0
        let helper: SIMD3<Double> = abs(n.z) < 0.9 ? SIMD3(0, 0, 1) : SIMD3(1, 0, 0)
        let e1 = normalize(cross(n, helper))
        let e2 = cross(n, e1)
        let off = offsetDeg * .pi / 180
        var runs: [[SIMD2<Double>]] = []
        var current: [SIMD2<Double>] = []
        for i in 0...samples {
            let t = 2 * Double.pi * Double(i) / Double(samples)
            let v = sin(off) * n + cos(off) * (cos(t) * e1 + sin(t) * e2)
            // drawn only where z ≥ 0 (stricter than project's −1e-9 cut), so a run ends on the primitive circle
            if let p = project(v), v.z >= 0 {
                current.append(p)
            } else if !current.isEmpty {
                runs.append(current)
                current = []
            }
        }
        if !current.isEmpty {
            // the run that crosses the start joins the first one
            if let first = runs.first, !runs.isEmpty, let head = first.first, let tail = current.last,
               length(head - tail) < 1e-12 {
                runs[0] = current + first.dropFirst()
            } else {
                runs.append(current)
            }
        }
        return runs
    }

    // MARK: - orix's `StereographicProjection(pole)` and `InverseStereographicProjection`

    /// −1 projects the upper hemisphere from the south pole, +1 the lower one from the north pole.
    public var pole: Int

    public init(pole: Int = -1) {
        precondition(pole == -1 || pole == 1, "the stereographic pole is -1 (upper hemisphere) or 1 (lower hemisphere)")
        self.pole = pole
    }

    /// The hemisphere this projection draws, orix's `SphericalRegion([0, 0, -pole])`: a vector is in it when its dot product
    /// with the region's normal is above −1e-9 (`v <= region`).
    public func isInHemisphere(_ v: SIMD3<Double>) -> Bool {
        Self.dot(SIMD3(0, 0, Double(-pole)), v) > -1e-9
    }

    /// `vector2xy`: the vectors of this hemisphere (the others are dropped, as orix drops them), each made unit, to
    /// (x, y) = −pole · (vx, vy) / (vz − pole); a zero denominator gives (0, 0).
    public func vectorToXY(_ vectors: [SIMD3<Double>]) -> [SIMD2<Double>] {
        vectors.filter(isInHemisphere).map { StereographicProjection.project($0, pole: pole) }
    }

    /// `_vector2xy` of one vector, whatever its hemisphere.
    public static func project(_ v: SIMD3<Double>, pole: Int) -> SIMD2<Double> {
        let u = unit(v)
        let p = Double(pole)
        let denom = u.z - p
        guard denom != 0 else { return SIMD2(0, 0) }
        return SIMD2(-p * u.x / denom, -p * u.y / denom)
    }

    /// `vector2xy_split`: the upper hemisphere's vectors projected with pole −1 and the lower one's with pole 1 (a vector on
    /// the equator is in both).
    public static func vectorToXYSplit(_ vectors: [SIMD3<Double>]) -> (upper: [SIMD2<Double>], lower: [SIMD2<Double>]) {
        (StereographicProjection(pole: -1).vectorToXY(vectors), StereographicProjection(pole: 1).vectorToXY(vectors))
    }

    /// `InverseStereographicProjection.xy2vector`: the unit vector of (x, y) in this projection's hemisphere.
    public func xyToVector(_ xy: SIMD2<Double>) -> SIMD3<Double> {
        let x = xy.x, y = xy.y
        let denom = 1 + x * x + y * y
        return SIMD3(2 * x / denom, 2 * y / denom, -Double(pole) * (1 - x * x - y * y) / denom)
    }

    /// orix's dot product, written out (the port's `Vector3.dot`).
    private static func dot(_ u: SIMD3<Double>, _ v: SIMD3<Double>) -> Double { u.x * v.x + u.y * v.y + u.z * v.z }

    /// `Vector3d.unit` (a zero vector stays zero), the port's `Vector3.unit`.
    private static func unit(_ v: SIMD3<Double>) -> SIMD3<Double> {
        let n = (v.x * v.x + v.y * v.y + v.z * v.z).squareRoot()
        return n == 0 ? SIMD3(0, 0, 0) : v / n
    }
}
