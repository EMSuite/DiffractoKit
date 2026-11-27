//
//  SlipTraces.swift
//  Moved from Electron Microscopy Studio's Sources/TEMSEM/Deformation (at 7b09bf5) by DiffractoKit#180; the bundles its moved
//  tests read are in Tests/DiffractoOrientationTests/Fixtures/orix-sem, the others stay with the app (docs/parity.md).
//  DefDAP 1.3.1's per-grain slip analysis (defdap.ebsd.Grain.calc_average_ori, calc_average_schmid_factors, calc_slip_traces):
//  the grain's average orientation, the Schmid factor of every slip system for a loading direction, and the angle and
//  inclination of every slip plane's trace on the map plane. Apache-2.0 (DefDAP, Copyright 2025 Mechanics of Microstructures
//  Group at The University of Manchester), see NOTICE. Method: Schmid & Boas (1935) for the factor |cos φ cos λ|; DefDAP
//  (Atkinson, Orozco-Caballero, da Fonseca et al.) for the trace convention.
//
//  Conventions as DefDAP's: the orientation is passive (sample to crystal), the load is transformed into the crystal as given
//  (DefDAP does not normalise it, so the factor scales with its squared length: pass a unit vector), a trace angle is measured
//  from the map's +y axis, counter-clockwise, in [0, π], and an inclination is the angle between the slip plane and the map
//  plane in [0, π/2]. A slip plane parallel to the map plane has no trace: its angle is NaN, as in DefDAP. Parity: bundles
//  `m5/slip-ti` (hcp) and `m5/slip-sdss` (fcc, bcc), Schmid factors and angles within 1e-6 (docs/parity.md).
//

import Foundation

/// The slip analysis of one orientation against a phase's slip systems.
public nonisolated enum SlipTraceAnalysis {
    /// One slip plane's trace on the map plane.
    public struct Trace: Sendable, Equatable {
        /// Radians from the map's +y axis, counter-clockwise, in [0, π]; NaN for a plane parallel to the map.
        public let angle: Double
        /// Radians between the slip plane and the map plane, in [0, π/2].
        public let inclination: Double
    }

    /// `Grain.calc_average_schmid_factors`: |(l·n)(l·d)| per system, grouped as the systems are, with l the load in crystal
    /// coordinates.
    static func schmidFactors(orientation q: SIMD4<Double>, load: SIMD3<Double>, systems: SlipSystemSet) -> [[Double]] {
        let l = DefDAPQuat.transform(q, load)
        return systems.groups.map { group in
            group.map { abs(Vector3.dot(l, $0.plane) * Vector3.dot(l, $0.direction)) }
        }
    }

    /// The largest Schmid factor over every system (the value DefDAP's Schmid-factor map shows when no plane is picked).
    static func maximumSchmidFactor(orientation q: SIMD4<Double>, load: SIMD3<Double>, systems: SlipSystemSet) -> Double {
        schmidFactors(orientation: q, load: load, systems: systems).joined().max() ?? .nan
    }

    /// `Grain.calc_slip_traces`: per slip-plane group (the first system's plane), the trace of the plane on the map plane.
    static func traces(orientation q: SIMD4<Double>, systems: SlipSystemSet) -> [Trace] {
        let screen = DefDAPQuat.transform(q, SIMD3(0, 0, 1))
        let back = DefDAPQuat.conjugate(q)
        return systems.groups.map { group in
            let n = group[0].plane
            let intersectionCrystal = Vector3.cross(screen, n)
            var inclination = acos(Vector3.dot(screen, n))
            if inclination > Double.pi / 2 { inclination = Double.pi - inclination }
            var intersection = DefDAPQuat.transform(back, intersectionCrystal)
            intersection = intersection / Vector3.dot(intersection, intersection).squareRoot()
            if intersection.x > 0 { intersection = -intersection }
            return Trace(angle: acos(intersection.y), inclination: inclination)
        }
    }

    /// `schmidFactors` of an orientation (lab to crystal, orix's and DefDAP's convention for the same Bunge angles; a grain's
    /// `GrainRow.orientation` is `Rotation(unit:)` of its quaternion).
    public static func schmidFactors(orientation: Rotation, load: SIMD3<Double>, systems: SlipSystemSet) -> [[Double]] {
        schmidFactors(orientation: orientation.q, load: load, systems: systems)
    }

    /// `maximumSchmidFactor` of an orientation.
    public static func maximumSchmidFactor(orientation: Rotation, load: SIMD3<Double>, systems: SlipSystemSet) -> Double {
        maximumSchmidFactor(orientation: orientation.q, load: load, systems: systems)
    }

    /// `traces` of an orientation.
    public static func traces(orientation: Rotation, systems: SlipSystemSet) -> [Trace] {
        traces(orientation: orientation.q, systems: systems)
    }

    /// `Grain.calc_average_ori` for a grain whose points (indices into `quaternions`) are listed in DefDAP's order.
    static func averageOrientation(of points: [Int], in quaternions: [SIMD4<Double>],
                                           structure: DefDAPCrystalStructure) -> SIMD4<Double> {
        DefDAPQuat.averageOrientation(points.map { quaternions[$0] }, structure: structure)
    }
}
