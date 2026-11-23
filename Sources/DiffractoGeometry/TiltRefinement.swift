//
//  TiltRefinement.swift
//  Role: the tilt of an orientation away from its zone axis, refined on the measured spot intensities (DiffractoKit S5b).
//        `ZoneAxisSolver` names the zone from the spot positions; which reflections are excited (the Laue circle) says how
//        far and which way the crystal is tilted off it. The tilt is the one whose kinematic intensities agree best with the
//        measured ones.
//
//  Agreement (a normalised correlation, as the template-matching score of ASTAR, Rauch & Dupuy, Arch. Metall. Mater. 50
//  (2005) 87, and py4DSTEM's orientation matching, Savitzky et al., Microsc. Microanal. 27 (2021) 712, take it): with
//  a_i = max(I_i, 0)^p for the measured spots and b_j = (|F_j|² exp(−s_j² / 2σ²))^p for the predicted ones (py4DSTEM's
//  intensity, s its excitation error, |s| ≤ 3σ, p = `intensityPower`), each measured spot paired with the nearest predicted
//  spot within the solver's tolerance,
//      C = Σ_paired a_i b_j / sqrt(Σ_all measured a_i² · Σ_predicted in range b_j²),
//  so a measured spot without a prediction and a predicted spot without a measurement both lower C. The predicted spots in
//  range lie within 2 % of the measured |q| range. Search: the tilt t = (t_x, t_y) about lab x and y (the crystal turned by
//  |t| about the in-plane axis t / |t|, `CrystalOrientation.tilted`'s frame) on a grid of `step` out to `reach`, then two
//  grids of step/5 and step/25 (11 x 11) around the best. The scale and in-plane turn stay those of the start.
//
//  Look-alike tilts: a zero-layer pattern's intensities can have more symmetry than the crystal shows down that axis (the
//  half turn about the beam, a mirror of the pattern that is no mirror of the crystal), so two tilts can agree equally. Every
//  refined orientation within `ambiguityMargin` of the best and more than 1° from it is returned as an alternative.
//  No Python function does this step alone; the pins are planted patterns at known tilts (docs/parity.md P44).
//

import Foundation
import simd
import DiffractoCore

/// One orientation after the tilt refinement.
public struct RefinedOrientation: Sendable, Equatable {
    public let orientation: CrystalOrientation
    /// The tilt applied to the start, radians about lab x and y (`CrystalOrientation.tilted`'s frame for small tilts).
    public let tilt: SIMD2<Double>
    /// The agreement C (0 to 1) after and before the refinement.
    public let agreement: Double
    public let startAgreement: Double
    /// The scale of the start (kept).
    public let scale: Double
}

/// The best refined orientation of a solution and the ones that agree as well.
public struct TiltRefinementResult: Sendable, Equatable {
    public let best: RefinedOrientation
    /// Different orientations (more than 1° from `best` and from each other) whose agreement is within the margin of the best:
    /// the intensities cannot tell them apart.
    public let alternatives: [RefinedOrientation]
    /// Every start refined (each candidate and its half turn about the beam), best first.
    public let refined: [RefinedOrientation]
}

extension ZoneAxisSolver {

    public struct TiltSettings: Sendable, Equatable {
        /// The largest tilt searched, radians. Default 4°.
        public var reach: Double
        /// The coarse grid step, radians. Default 0.25°.
        public var step: Double
        /// p of a = I^p and b = (|F|² e^(−s²/2σ²))^p. Default 0.5 (amplitudes).
        public var intensityPower: Double
        /// Refined orientations within this much of the best agreement are alternatives. Default 1e-3.
        public var ambiguityMargin: Double

        public init(reach: Double = 4 * .pi / 180, step: Double = 0.25 * .pi / 180, intensityPower: Double = 0.5,
                    ambiguityMargin: Double = 1e-3) {
            self.reach = reach
            self.step = step
            self.intensityPower = intensityPower
            self.ambiguityMargin = ambiguityMargin
        }

        /// Every value finite and in range (reach ≥ 0, step > 0 with at most 1000 steps, power > 0, margin ≥ 0).
        public var isValid: Bool {
            reach.isFinite && step.isFinite && reach >= 0 && step > 0 && reach / step <= 1000 && intensityPower.isFinite
                && intensityPower > 0 && ambiguityMargin.isFinite && ambiguityMargin >= 0
        }
    }

    /// The agreement C of `orientation` at `scale` with the measured spots (Å⁻¹) and their intensities (any unit; one per
    /// spot): 0 to 1, −1 when no predicted spot falls in the measured range.
    public func intensityAgreement(of orientation: CrystalOrientation, scale: Double, spots: [SIMD2<Double>],
                                   intensities: [Double], intensityPower: Double = 0.5) -> Double {
        guard spots.count == intensities.count, !spots.isEmpty else { return -1 }
        let a = intensities.map { pow(max($0, 0), intensityPower) }
        let qn = spots.map { simd_length($0) }
        return agreement(orientation.matrix, scale: scale, spots: spots, qn: qn, a: a, power: intensityPower)
    }

    /// `start` with its tilt refined on the intensities.
    public func refineTilt(_ start: CrystalOrientation, scale: Double, spots: [SIMD2<Double>], intensities: [Double],
                           settings: TiltSettings = .init()) -> RefinedOrientation? {
        guard settings.isValid, spots.count == intensities.count, !spots.isEmpty, scale.isFinite, scale > 0 else { return nil }
        let a = intensities.map { pow(max($0, 0), settings.intensityPower) }
        let qn = spots.map { simd_length($0) }
        let m0 = start.matrix
        func c(_ t: SIMD2<Double>) -> Double {
            agreement(Self.tilting(m0, by: t), scale: scale, spots: spots, qn: qn, a: a, power: settings.intensityPower)
        }
        let n = Int((settings.reach / settings.step).rounded())
        var best = (value: -2.0, tilt: SIMD2<Double>(0, 0))
        for i in -n...n {
            for k in -n...n where (Double(i * i + k * k)).squareRoot() <= Double(n) {
                let t = SIMD2(Double(i), Double(k)) * settings.step
                let v = c(t)
                if v > best.value { best = (v, t) }
            }
        }
        for fine in [settings.step / 5, settings.step / 25] {
            let centre = best.tilt
            for i in -5...5 {
                for k in -5...5 {
                    let t = centre + SIMD2(Double(i), Double(k)) * fine
                    let v = c(t)
                    if v > best.value { best = (v, t) }
                }
            }
        }
        guard let o = try? CrystalOrientation(matrix: Self.tilting(m0, by: best.tilt)) else { return nil }
        return RefinedOrientation(orientation: o, tilt: best.tilt, agreement: best.value, startAgreement: c(.zero), scale: scale)
    }

    /// How far (radians) from `orientation` the agreement stays within `settings.ambiguityMargin` of its own, over the coarse
    /// tilt grid about it (`reach`, `step`): 0 when every other tilt agrees worse, up to `reach` when the intensities barely
    /// fix the tilt (few spots).
    public func tiltSpread(of orientation: CrystalOrientation, scale: Double, spots: [SIMD2<Double>], intensities: [Double],
                           settings: TiltSettings = .init()) -> Double {
        guard settings.isValid, spots.count == intensities.count, !spots.isEmpty, scale.isFinite, scale > 0 else { return .nan }
        let a = intensities.map { pow(max($0, 0), settings.intensityPower) }
        let qn = spots.map { simd_length($0) }
        let m0 = orientation.matrix
        let c0 = agreement(m0, scale: scale, spots: spots, qn: qn, a: a, power: settings.intensityPower)
        let n = Int((settings.reach / settings.step).rounded())
        var spread = 0.0
        for i in -n...n {
            for k in -n...n where (Double(i * i + k * k)).squareRoot() <= Double(n) {
                let t = SIMD2(Double(i), Double(k)) * settings.step
                guard simd_length(t) > spread else { continue }
                if agreement(Self.tilting(m0, by: t), scale: scale, spots: spots, qn: qn, a: a, power: settings.intensityPower)
                    >= c0 - settings.ambiguityMargin { spread = simd_length(t) }
            }
        }
        return spread
    }

    /// The candidates of `solution` (each and its half turn about the beam) refined on the intensities; nil when there is
    /// none to refine or the input does not fit.
    public func refineTilts(of solution: ZoneAxisSolution, spots: [SIMD2<Double>], intensities: [Double],
                            settings: TiltSettings = .init()) -> TiltRefinementResult? {
        var refined: [RefinedOrientation] = []
        for candidate in solution.candidates {
            for start in [candidate.orientation, candidate.orientation.rotatedInPlane(by: .pi)] {
                if let r = refineTilt(start, scale: candidate.scale, spots: spots, intensities: intensities, settings: settings) {
                    refined.append(r)
                }
            }
        }
        guard !refined.isEmpty else { return nil }
        refined.sort { $0.agreement > $1.agreement }
        let best = refined[0]
        let distinct = 1 * Double.pi / 180
        var alternatives: [RefinedOrientation] = []
        for r in refined.dropFirst() where r.agreement >= best.agreement - settings.ambiguityMargin {
            let others = [best] + alternatives
            if others.allSatisfy({ angle(between: $0.orientation, and: r.orientation, includingHalfTurn: false) > distinct }) {
                alternatives.append(r)
            }
        }
        return TiltRefinementResult(best: best, alternatives: alternatives, refined: refined)
    }

    /// The crystal turned by |t| about the lab axis (t_x, t_y, 0) / |t|: the columns become M Rᵀ.
    static func tilting(_ m: simd_double3x3, by t: SIMD2<Double>) -> simd_double3x3 {
        let angle = simd_length(t)
        guard angle > 1e-12 else { return m }
        let r = simd_double3x3(simd_quatd(angle: angle, axis: SIMD3(t.x / angle, t.y / angle, 0)))
        return m * r.transpose
    }

    func agreement(_ m: simd_double3x3, scale: Double, spots: [SIMD2<Double>], qn: [Double], a: [Double], power: Double)
        -> Double {
        let sigma = settings.pattern.sigmaExcitationError
        let reach = sigma * settings.pattern.excitationMultiple
        let qMin = (qn.min() ?? 0) * 0.98, qMax = (qn.max() ?? 0) * 1.02
        let mt = m.transpose
        var positions: [SIMD2<Double>] = []
        var b: [Double] = []
        var excited = 0
        for j in g.indices {
            let p = mt * g[j]
            let s = (2 * p.z - wavelength * gLength[j] * gLength[j]) / (2 - 2 * wavelength * p.z)
            guard s.isFinite, abs(s) <= reach else { continue }
            excited += 1
            let xy = scale * SIMD2(p.x, p.y)
            let r = simd_length(xy)
            guard r >= qMin, r <= qMax else { continue }
            positions.append(xy)
            b.append(pow(intensity[j] * exp(-s * s / (2 * sigma * sigma)), power))
        }
        guard excited >= 2, !positions.isEmpty else { return -1 }
        var numerator = 0.0
        for (i, q) in spots.enumerated() {
            let radius = max(settings.tolerance * qn[i], settings.floor)
            var nearest = Double.infinity, jBest = 0
            for (j, p) in positions.enumerated() {
                let d = simd_distance(q, p)
                if d < nearest { nearest = d; jBest = j }
            }
            if nearest <= radius { numerator += a[i] * b[jBest] }
        }
        let denominator = (a.reduce(0.0) { $0 + $1 * $1 } * b.reduce(0.0) { $0 + $1 * $1 }).squareRoot()
        return denominator > 0 ? numerator / denominator : -1
    }
}
