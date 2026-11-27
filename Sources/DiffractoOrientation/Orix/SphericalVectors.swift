//
//  SphericalVectors.swift
//  orix 0.15.0's `orix.vector.Vector3d` spherical coordinates and `orix.projections` (the stereographic projection and its
//  inverse, upper and lower hemisphere) on `SIMD3<Double>`, the vector type of the orientation core. GPL-3.0 (orix is
//  GPL-3.0-or-later), see NOTICE. Parity: `tools/parity/sem` bundle `m3/ipf-color-keys`, coordinates within 1e-7.
//

import Foundation

extension Vector3 {
    /// `Vector3d.azimuth`: atan2(y, x) in [0, 2π), each of x and y set to 0 first where it is within numpy's `isclose`
    /// of 0 (1e-8), as orix does against rounding.
    public static func azimuth(_ v: SIMD3<Double>) -> Double {
        let x = abs(v.x) <= 1e-8 ? 0 : v.x
        let y = abs(v.y) <= 1e-8 ? 0 : v.y
        let a = atan2(y, x)
        return a < 0 ? a + 2 * Double.pi : a
    }

    /// `Vector3d.polar`: arccos(z / |v|) in [0, π].
    public static func polar(_ v: SIMD3<Double>) -> Double {
        acos(v.z / RotationConversions.norm(v))
    }

    /// `Vector3d.get_circle`: `steps` vectors on the circle at `openingAngle` about `v` (a great circle by default), from
    /// rotating `v` by `openingAngle` about its perpendicular and then by numpy's `linspace(0, 2π, steps)` about itself.
    public static func circle(about v: SIMD3<Double>, openingAngle: Double = Double.pi / 2, steps: Int = 100) -> [SIMD3<Double>] {
        let r1 = Rotation.fromAxisAngle(perpendicular(v), openingAngle)
        return linspace(0, 2 * Double.pi, steps).map { w in (Rotation.fromAxisAngle(v, w) * r1) * v }
    }

    /// numpy's `linspace(start, stop, num)` (endpoint included): start + i · step, the last value exactly `stop`.
    public static func linspace(_ start: Double, _ stop: Double, _ num: Int) -> [Double] {
        guard num > 1 else { return num == 1 ? [start] : [] }
        let step = (stop - start) / Double(num - 1)
        var out = (0..<num).map { start + Double($0) * step }
        out[num - 1] = stop
        return out
    }

    /// numpy's `mod` on floats: the remainder with the sign of the divisor (a zero remainder is +0, NaN stays NaN).
    static func pythonMod(_ a: Double, _ b: Double) -> Double {
        var m = fmod(a, b)
        if m != 0 {
            if (b < 0) != (m < 0) { m += b }
        } else if !m.isNaN {
            m = 0
        }
        return m
    }
}

// orix's `StereographicProjection` and its inverse are DiffractoGeometry's `StereographicProjection` (one type).
