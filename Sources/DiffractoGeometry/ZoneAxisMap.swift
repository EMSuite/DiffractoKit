//
//  ZoneAxisMap.swift
//  Role: the zone-axis map of an indexed crystal (2.4 tilt planner, S12): its low-index zone axes and the traces of its
//        lowest plane families around the beam, in the pattern frame of a `CrystalOrientation` and projected with
//        `StereographicProjection` (beam at the centre). Angles between zone axes are exact whatever the holder; the
//        tilts to reach them need a calibrated `DoubleTiltHolder`.
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering

/// The zone axes and plane traces of one orientation, for the stereographic map.
public enum ZoneAxisMap {

    /// A zone axis [u v w] near the beam.
    public struct Pole: Sendable, Equatable {
        /// [u v w], no common factor.
        public let uvw: SIMD3<Int>
        /// Its unit direction in the pattern frame (lab x, y; z toward the source).
        public let direction: SIMD3<Double>
        /// Its stereographic position.
        public let position: SIMD2<Double>
        /// The angle between it and the beam (degrees).
        public let angleFromBeamDeg: Double
    }

    /// The trace of a lattice plane (hkl): the great circle of the directions in the plane.
    public struct Trace: Sendable, Equatable {
        /// One member of each ±(hkl) pair of the family.
        public let hkl: MillerIndex
        /// The family's name (the kit's orbit name, e.g. "{1 1 1}").
        public let family: String
        /// The plane normal (g) as a unit vector in the pattern frame.
        public let normal: SIMD3<Double>
        /// d (Å).
        public let d: Double
        /// The Bragg angle θ at the given wavelength (degrees): a Kikuchi band's edges lie at ±θ from the trace.
        public let braggAngleDeg: Double?
        /// The trace as runs of stereographic points (``StereographicProjection/circle(normal:offsetDeg:samples:)``).
        public let points: [[SIMD2<Double>]]
    }

    /// Every zone axis [u v w] with |u|, |v|, |w| ≤ `maxIndex` and no common factor whose direction u a₁ + v a₂ + w a₃ lies
    /// within `withinDeg` of the beam, ordered by the angle (ties by the smaller largest index, then the larger (u, v, w)).
    public static func poles(structure: CrystalStructure, orientation: CrystalOrientation, maxIndex: Int = 3,
                             withinDeg: Double = 60) -> [Pole] {
        let lattice = structure.cell.realLattice
        guard maxIndex >= 1, maxIndex <= 64, lattice.count == 3, structure.cell.isValid, withinDeg.isFinite else { return [] }
        var out: [(pole: Pole, key: Int)] = []
        for u in -maxIndex...maxIndex { for v in -maxIndex...maxIndex { for w in -maxIndex...maxIndex {
            guard gcd(gcd(abs(u), abs(v)), abs(w)) == 1 else { continue }
            let d = Double(u) * lattice[0] + Double(v) * lattice[1] + Double(w) * lattice[2]
            let lab = normalize(orientation.labCoordinates(of: d))
            let angle = atan2((lab.x * lab.x + lab.y * lab.y).squareRoot(), lab.z) * 180 / .pi
            guard angle <= withinDeg, let p = StereographicProjection.project(lab) else { continue }
            out.append((Pole(uvw: SIMD3(u, v, w), direction: lab, position: p, angleFromBeamDeg: angle),
                        max(abs(u), abs(v), abs(w))))
        } } }
        out.sort { x, y in
            if abs(x.pole.angleFromBeamDeg - y.pole.angleFromBeamDeg) > 1e-9 {
                return x.pole.angleFromBeamDeg < y.pole.angleFromBeamDeg
            }
            if x.key != y.key { return x.key < y.key }
            return (x.pole.uvw.x, x.pole.uvw.y, x.pole.uvw.z) > (y.pole.uvw.x, y.pole.uvw.y, y.pole.uvw.z)
        }
        return out.map { $0.pole }
    }

    /// The traces of the `families` plane families of largest d among the reflections the structure scatters into (the
    /// kinematic list after the |F| cut, electron form factor, out to `kMax`), one trace per ±(hkl) pair. With a wavelength
    /// each carries its Bragg angle. Empty when the list is refused.
    public static func planeTraces(structure: CrystalStructure, orientation: CrystalOrientation, families: Int = 3,
                                   kMax: Double = 1.0, wavelengthAngstroms: Double? = nil, samples: Int = 360) -> [Trace] {
        guard families >= 1, let list = try? structure.kinematicReflections(kMax: kMax, retainingAll: false) else { return [] }
        var out: [Trace] = []
        for kf in list.families.prefix(families) {
            var seen = Set<MillerIndex>()
            for hkl in kf.family.members {
                let negative = MillerIndex(-hkl.h, -hkl.k, -hkl.l)
                if seen.contains(negative) { continue }
                seen.insert(hkl)
                let g = structure.cell.gVector(hkl)
                let normal = normalize(orientation.labCoordinates(of: g))
                var bragg: Double?
                if let lambda = wavelengthAngstroms, lambda.isFinite, lambda > 0, lambda * kf.family.g / 2 <= 1 {
                    bragg = asin(lambda * kf.family.g / 2) * 180 / .pi
                }
                out.append(Trace(hkl: hkl, family: kf.family.name, normal: normal, d: kf.family.d, braggAngleDeg: bragg,
                                 points: StereographicProjection.circle(normal: normal, samples: samples)))
            }
        }
        return out
    }

    private static func gcd(_ a: Int, _ b: Int) -> Int { b == 0 ? a : gcd(b, a % b) }
}
