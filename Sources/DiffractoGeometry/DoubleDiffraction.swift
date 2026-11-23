//
//  DoubleDiffraction.swift
//  Role: the spots double diffraction adds to kinematic patterns (DiffractoKit S3). A beam diffracted by g₁ is diffracted
//        again by g₂ and leaves along g₁ + g₂. Where g₁ + g₂ is a reflection the structure forbids (|F| at or below the cut,
//        e.g. 002 in silicon down [110]) or a point on no lattice at all (g₁ of a matrix plus g₂ of a precipitate), a spot
//        appears that the kinematic pattern does not have. Each added spot names its routes (g₁ of which crystal, g₂ of which).
//
//  Rule (Williams & Carter, Transmission Electron Microscopy, 2nd ed. (2009), section 17.5): both g₁ and g₂ are spots of
//  their kinematic patterns (`SpotPatternSimulator.pattern`, the same settings), and the sum is excited by the same window,
//  |s(g₁ + g₂)| ≤ σ · multiple, with s from the lab vector of the sum (py4DSTEM's form, `SimulatedSpot.excitationError`).
//  Within one crystal the sum is added only when its hkl is not a kinematic reflection (|F| at or below the cut) and lies
//  within k_max; across crystals, only when it falls on none of the patterns' own spots and within the larger k_max of the two. No intensity is computed (kinematic
//  theory gives none): each route carries I(g₁)·I(g₂) to rank the routes.
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering

/// One crystal of a double-diffraction calculation: its reflections and its orientation.
public struct DiffractingCrystal: Sendable {
    public let simulator: SpotPatternSimulator
    public let orientation: CrystalOrientation

    public init(simulator: SpotPatternSimulator, orientation: CrystalOrientation) {
        self.simulator = simulator
        self.orientation = orientation
    }
}

/// g₁ of one crystal followed by g₂ of one crystal.
public struct DoubleDiffractionRoute: Sendable, Equatable {
    public struct Step: Sendable, Equatable, Hashable {
        /// The index of the crystal in the call's list.
        public let crystal: Int
        public let hkl: MillerIndex
    }
    public let first: Step
    public let second: Step
    /// I(g₁) · I(g₂), Å⁻⁸: ranks the routes, not an intensity.
    public let weight: Double
}

/// A spot that double diffraction adds.
public struct DoubleDiffractionSpot: Sendable, Equatable {
    /// (qx, qy) of g₁ + g₂, Å⁻¹, as `SimulatedSpot.position`.
    public let position: SIMD2<Double>
    /// The lab z of g₁ + g₂, Å⁻¹.
    public let qz: Double
    /// s of g₁ + g₂ (py4DSTEM's form), Å⁻¹.
    public let excitationError: Double
    /// The crystal and the hkl, when every route stays within one crystal; nil for a spot between crystals.
    public let crystal: Int?
    public let hkl: MillerIndex?
    /// Every route that lands here, the largest weight first.
    public let routes: [DoubleDiffractionRoute]
}

/// The kinematic patterns and the spots double diffraction adds to them.
public struct DoubleDiffractionPatterns: Sendable, Equatable {
    /// One kinematic pattern per crystal, in the call's order.
    public let patterns: [SpotPattern]
    /// The added spots, the strongest route first.
    public let added: [DoubleDiffractionSpot]
}

public enum DoubleDiffraction {

    /// Two spots closer than this (Å⁻¹, in all three lab coordinates) are one spot.
    public static let coincidence = 1e-6

    /// The added spots of one crystal at one orientation.
    public static func spots(simulator: SpotPatternSimulator, orientation: CrystalOrientation, volts: Double,
                             settings: SpotPatternSimulator.Settings = .init()) throws(SpotPatternError) -> DoubleDiffractionPatterns {
        try spots(crystals: [DiffractingCrystal(simulator: simulator, orientation: orientation)], volts: volts, settings: settings)
    }

    /// The added spots of several crystals in one beam (a precipitate in its matrix): routes within each crystal and between
    /// any two. Throws for an empty list and what `SpotPatternSimulator.pattern` refuses.
    public static func spots(crystals: [DiffractingCrystal], volts: Double,
                             settings: SpotPatternSimulator.Settings = .init()) throws(SpotPatternError) -> DoubleDiffractionPatterns {
        guard !crystals.isEmpty else { throw .invalidSettings("double diffraction needs at least one crystal") }
        var patterns: [SpotPattern] = []
        for c in crystals {
            patterns.append(try c.simulator.pattern(orientation: c.orientation, volts: volts, settings: settings))
        }
        let lambda = patterns[0].wavelength
        let sgMax = settings.sigmaExcitationError * settings.excitationMultiple
        // Each spot's lab vector and intensity.
        struct Lab { let crystal: Int; let hkl: MillerIndex; let p: SIMD3<Double>; let intensity: Double }
        var labs: [[Lab]] = []
        for (i, c) in crystals.enumerated() {
            let cell = c.simulator.structure.cell
            labs.append(patterns[i].spots.map {
                Lab(crystal: i, hkl: $0.hkl, p: c.orientation.labCoordinates(of: cell.gVector($0.hkl)), intensity: $0.intensity)
            })
        }
        let allowed = crystals.map { Set($0.simulator.reflections.map(\.structureFactor.hkl)) }
        let kinematic = labs.flatMap { $0 }.map(\.p)

        func excitation(_ p: SIMD3<Double>) -> Double {
            (2 * p.z - lambda * simd_length_squared(p)) / (2 - 2 * lambda * p.z)
        }
        struct Candidate { var p: SIMD3<Double>; var s: Double; var crystal: Int?; var hkl: MillerIndex?; var routes: [DoubleDiffractionRoute] }
        var candidates: [Candidate] = []
        func offer(_ a: Lab, _ b: Lab) {
            let route = DoubleDiffractionRoute(first: .init(crystal: a.crystal, hkl: a.hkl), second: .init(crystal: b.crystal, hkl: b.hkl),
                                               weight: a.intensity * b.intensity)
            var p: SIMD3<Double>, hkl: MillerIndex?, crystal: Int?
            if a.crystal == b.crystal {
                let sum = MillerIndex(a.hkl.h + b.hkl.h, a.hkl.k + b.hkl.k, a.hkl.l + b.hkl.l)
                guard !sum.isZero, !allowed[a.crystal].contains(sum) else { return }
                let c = crystals[a.crystal]
                let g = c.simulator.structure.cell.gVector(sum)
                guard simd_length(g) <= c.simulator.kMax else { return }
                p = c.orientation.labCoordinates(of: g)
                hkl = sum
                crystal = a.crystal
            } else {
                p = a.p + b.p
                let reach = max(crystals[a.crystal].simulator.kMax, crystals[b.crystal].simulator.kMax)
                guard simd_length(p) > coincidence, simd_length(p) <= reach else { return }
                guard !kinematic.contains(where: { simd_distance($0, p) <= coincidence }) else { return }
            }
            let s = excitation(p)
            guard s.isFinite, abs(s) <= sgMax else { return }
            if let k = candidates.firstIndex(where: { simd_distance($0.p, p) <= coincidence }) {
                candidates[k].routes.append(route)
                if candidates[k].crystal != crystal || candidates[k].hkl != hkl { candidates[k].crystal = nil; candidates[k].hkl = nil }
            } else {
                candidates.append(Candidate(p: p, s: s, crystal: crystal, hkl: hkl, routes: [route]))
            }
        }
        for i in labs.indices {
            for a in labs[i].indices { for b in a..<labs[i].count { offer(labs[i][a], labs[i][b]) } }
            for j in labs.indices where j > i {
                for a in labs[i] { for b in labs[j] { offer(a, b) } }
            }
        }
        let added = candidates.map { c in
            DoubleDiffractionSpot(position: SIMD2(c.p.x, c.p.y), qz: c.p.z, excitationError: c.s, crystal: c.crystal, hkl: c.hkl,
                                  routes: c.routes.sorted { $0.weight > $1.weight })
        }.sorted { ($0.routes.first?.weight ?? 0) > ($1.routes.first?.weight ?? 0) }
        return DoubleDiffractionPatterns(patterns: patterns, added: added)
    }
}
