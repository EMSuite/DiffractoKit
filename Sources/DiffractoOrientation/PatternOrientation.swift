//
//  PatternOrientation.swift
//  Role: the full orientation of one measured spot pattern of a known structure (DiffractoKit S6, Electron Microscopy
//        Studio's Index). Three steps, each a method the kit already holds:
//          1. the template match (`ACOMMatcher.matches`, py4DSTEM's polar-correlation method: Rauch & Dupuy, Arch. Metall.
//             Mater. 50 (2005) 87; Ophus et al., Microsc. Microanal. 28 (2022) 390) gives its best `templateMatches`
//             orientations at distinct beam directions;
//          2. the zone-axis solver (`ZoneAxisSolver.solve`, indexing by lattice-vector pairs) gives its best
//             `solverCandidates`, with hkl, scale and ties;
//          3. every one of those starts, and its half turn about the beam, gets its in-plane turn and scale fitted to the spot
//             positions (`ZoneAxisSolver.fit`), its tilt refined on the spot intensities (`ZoneAxisSolver.refineTilt`, S5b) and
//             its positions fitted again. The best intensity agreement C wins.
//        Look-alikes: a pattern can be its own mirror image where the crystal has no such mirror (the zero layer of a
//        low-symmetry zone axis). The mirror is a half turn about an in-plane axis, a proper rotation, so the winner turned
//        about each of the pattern's candidate mirror lines (along and across its innermost spots) is fitted and refined too.
//        Every refined orientation within `ambiguityMargin` of the best C and more than 1° from it and from each other is an
//        alternative: the pattern cannot tell them apart. The spread of the tilt is how far from the best tilt the agreement
//        stays within `spreadMargin` (0.005): a pattern of few spots fixes its tilt poorly. An answer with an alternative or
//        a spread over 1° is flagged.
//
//  No Python function does this; py4DSTEM's `match_single_pattern` is the template step alone (a 2° grid). The pins are
//  py4DSTEM-planted patterns with py4DSTEM's own answer (docs/parity.md P45).
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoGeometry

/// Why a ``PatternOrientation`` could not be built.
public enum PatternOrientationError: Error, Sendable, Equatable, CustomStringConvertible {
    case solver(ZoneAxisSolverError)
    case templates(ACOMError)
    case invalidSettings(String)

    public var description: String {
        switch self {
        case .solver(let e): "zone-axis solver: \(e)"
        case .templates(let e): "template library: \(e)"
        case .invalidSettings(let s): "invalid settings: \(s)"
        }
    }
}

/// One orientation of a pattern after the three steps.
public struct IndexedOrientation: Sendable, Equatable {
    /// Which step proposed the start.
    public enum Source: Sendable, Equatable {
        /// The n-th best template match (0 = the best).
        case template(Int)
        /// The n-th candidate of the zone-axis solver.
        case solver(Int)
        /// A look-alike of the winner: the winner turned half about an in-plane axis of the pattern.
        case lookAlike
    }

    /// The full orientation (py4DSTEM's convention, ``CrystalOrientation``).
    public let orientation: CrystalOrientation
    /// The measured |q| over the predicted |g| (1 = the calibration is exact).
    public let scale: Double
    /// The intensity agreement C (0 to 1, `ZoneAxisSolver.intensityAgreement`) that ranks the orientations.
    public let agreement: Double
    /// The positions' fit at this orientation: zone axis, F_β, matched spots, rms distance and each spot's hkl.
    public let positions: ZoneAxisCandidate
    /// How far (radians) from this tilt the agreement stays within the spread margin, on the tilt grid about it.
    public let tiltSpread: Double
    public let source: Source
}

/// The answer for one pattern.
public struct PatternOrientationResult: Sendable, Equatable {
    /// The orientation of the best intensity agreement.
    public let best: IndexedOrientation
    /// Orientations more than 1° from `best` and from each other whose agreement is within the margin of the best: the pattern
    /// cannot tell them apart. Empty when the answer is unique.
    public let alternatives: [IndexedOrientation]
    /// The zone-axis solver's own answer (verdict, ties, refusal) for the same spots.
    public let solution: ZoneAxisSolution

    /// The answer is flagged: an alternative exists (the pattern cannot tell orientations apart), or the tilt is spread wider
    /// than the limit (the intensities fix it poorly). On the planted sets no answer more than 2° off was unflagged.
    public func isFlagged(tiltSpreadLimit: Double = PatternOrientation.Settings.defaultTiltSpreadLimit) -> Bool {
        !alternatives.isEmpty || !(best.tiltSpread <= tiltSpreadLimit)
    }
}

/// Indexes one spot pattern of a known structure to its full orientation. Build it once per structure and voltage (the
/// template library is the costly part).
public struct PatternOrientation: Sendable {

    public struct Settings: Sendable, Equatable {
        /// Template matches refined (py4DSTEM's `num_matches_return`). Default 5.
        public var templateMatches: Int
        /// Zone-axis solver candidates refined. Default 5.
        public var solverCandidates: Int
        /// The tilt search of each start (`ZoneAxisSolver.TiltSettings`; its `ambiguityMargin` is the margin here too).
        public var tilt: ZoneAxisSolver.TiltSettings
        /// The template library: py4DSTEM's excitation error, σ = 0.02 Å⁻¹ (width σ√2), `corr_kernel_size` 0.08 and
        /// `power_intensity` 0.25 by default; the wavelength is set from the voltage.
        public var templates: ACOMTemplateSettings
        /// A look-alike start is refined when its position F_β is within this of the winner's. Default 0.05.
        public var lookAlikeMargin: Double
        /// The tilt spread counts the tilts whose agreement is within this of the orientation's own. Default 0.005.
        public var spreadMargin: Double
        /// `PatternOrientationResult.isFlagged`'s default limit on the tilt spread: 1°.
        public static let defaultTiltSpreadLimit = 1 * Double.pi / 180

        public init(templateMatches: Int = 5, solverCandidates: Int = 5, tilt: ZoneAxisSolver.TiltSettings = .init(),
                    templates: ACOMTemplateSettings = ACOMTemplateSettings(sgWidth: 0.02 * 2.0.squareRoot()),
                    lookAlikeMargin: Double = 0.05, spreadMargin: Double = 0.005) {
            self.templateMatches = templateMatches
            self.solverCandidates = solverCandidates
            self.tilt = tilt
            self.templates = templates
            self.lookAlikeMargin = lookAlikeMargin
            self.spreadMargin = spreadMargin
        }
    }

    public let solver: ZoneAxisSolver
    public let matcher: ACOMMatcher
    public let settings: Settings

    /// - Parameters:
    ///   - kMax: the reflections' reach, Å⁻¹, of the template library (py4DSTEM's `k_max`); the solver reaches 0.1 Å⁻¹
    ///     further. Default 1.5.
    ///   - solverSettings: the zone-axis solver's settings (tolerance, scale window, excitation error).
    public init(structure: CrystalStructure, volts: Double, kMax: Double = 1.5,
                solverSettings: ZoneAxisSolver.Settings = .init(), settings: Settings = .init())
        throws(PatternOrientationError) {
        guard settings.templateMatches >= 0, settings.solverCandidates >= 0,
              settings.templateMatches + settings.solverCandidates > 0 else {
            throw .invalidSettings("templateMatches and solverCandidates must be non-negative, one of them positive")
        }
        guard settings.tilt.isValid, settings.lookAlikeMargin.isFinite, settings.lookAlikeMargin >= 0,
              settings.spreadMargin.isFinite, settings.spreadMargin >= 0 else {
            throw .invalidSettings("tilt settings, lookAlikeMargin and spreadMargin must be finite and in range")
        }
        // the solver's reflections reach past the longest spot at the largest scale it fits
        do { solver = try ZoneAxisSolver(structure: structure, volts: volts, kMax: kMax + 0.1, settings: solverSettings) }
        catch { throw .solver(error) }
        var templates = settings.templates
        templates.excitation = .py4DSTEM(wavelengthAngstrom: ElectronWavelength.angstroms(volts: volts))
        let library: ACOMTemplateLibrary
        do { library = try ACOMTemplateLibrary.generate(crystal: structure, kMax: kMax, settings: templates) }
        catch { throw .templates(error) }
        guard let m = ACOMMatcher(library: library, refinesInPlaneAngle: true) else {
            throw .invalidSettings("the template library's azimuthal length has no FFT plan")
        }
        matcher = m
        self.settings = settings
    }

    /// The full orientation of `spots` (Å⁻¹ in the pattern plane, direct beam removed) with their `intensities` (one per spot,
    /// any unit). Nil when no start explains two or more spots, or the input does not fit (counts differ, fewer than two
    /// finite spots).
    public func index(spots: [SIMD2<Double>], intensities: [Double]) -> PatternOrientationResult? {
        guard spots.count == intensities.count else { return nil }
        let kept = spots.indices.filter {
            spots[$0].x.isFinite && spots[$0].y.isFinite && simd_length(spots[$0]) > 0 && intensities[$0].isFinite
        }
        let q = kept.map { spots[$0] }, intensity = kept.map { max(intensities[$0], 0) }
        guard q.count >= 2 else { return nil }

        var starts: [(CrystalOrientation, Double, IndexedOrientation.Source)] = []
        // 1. templates, in Å⁻¹ as pixels of 1e-3 Å⁻¹ about the origin
        let pixel = 1e-3
        let peaks = zip(q, intensity).map { ACOMPeak(x: Float($0.0.x / pixel), y: Float($0.0.y / pixel), intensity: Float($0.1)) }
        if settings.templateMatches > 0 {
            for (n, r) in matcher.matches(peaks: peaks, originX: 0, originY: 0, invAngstromPerPixel: pixel,
                                          count: settings.templateMatches).enumerated() {
                let m = matcher.library.orientationMatrix(templateIndex: r.templateIndex, inPlaneAngle: r.inPlaneAngle,
                                                          mirrored: r.mirrored)
                guard let o = try? CrystalOrientation(matrix: Self.orthonormalized(m)) else { continue }
                // the template's in-plane angle is an azimuthal bin (2.8° at 128); the position fit only pulls in turns of
                // about 0.6° (the 1 % tolerance on the innermost spots), so the start is the best of turns ±1.5° in 0.5°
                let turns = stride(from: -1.5, through: 1.5, by: 0.5).map { $0 * Double.pi / 180 }
                let fitted = turns.compactMap { solver.fit(o.rotatedInPlane(by: $0), spots: q) }
                    .max { $0.agreement < $1.agreement }
                starts.append((fitted?.orientation ?? o, fitted?.scale ?? 1, .template(n)))
            }
        }
        // 2. the zone-axis solver
        let solution = solver.solve(q, top: max(settings.solverCandidates, 1))
        for (n, c) in solution.candidates.prefix(settings.solverCandidates).enumerated() {
            starts.append((c.orientation, c.scale, .solver(n)))
        }
        // 3. each start and its half turn: positions, tilt, positions
        var refined = starts.flatMap { refine($0.0, scale: $0.1, source: $0.2, q: q, intensity: intensity) }
        guard let first = refined.max(by: { $0.agreement < $1.agreement }) else { return nil }
        for start in lookAlikeStarts(of: first, q: q) {
            refined += refine(start.orientation, scale: start.scale, source: .lookAlike, q: q, intensity: intensity)
        }
        // best first; equal agreements keep their order (templates, solver, look-alikes)
        refined = refined.enumerated().sorted { a, b in
            a.element.agreement != b.element.agreement ? a.element.agreement > b.element.agreement : a.offset < b.offset
        }.map(\.element)
        let best = withSpread(refined[0], q: q, intensity: intensity)
        let distinct = 1 * Double.pi / 180
        var alternatives: [IndexedOrientation] = []
        for r in refined.dropFirst() where r.agreement >= best.agreement - settings.tilt.ambiguityMargin {
            let others = [best] + alternatives
            if others.allSatisfy({ solver.angle(between: $0.orientation, and: r.orientation, includingHalfTurn: false) > distinct }) {
                alternatives.append(withSpread(r, q: q, intensity: intensity))
            }
        }
        return PatternOrientationResult(best: best, alternatives: alternatives, solution: solution)
    }

    // MARK: - steps

    /// One start and its half turn about the beam, each: positions fitted, tilt refined, positions fitted again (kept when C
    /// does not drop). Both are returned: a start and its half turn can agree equally (a look-alike), and only the
    /// alternatives can show that.
    private func refine(_ start: CrystalOrientation, scale: Double, source: IndexedOrientation.Source,
                        q: [SIMD2<Double>], intensity: [Double]) -> [IndexedOrientation] {
        var out: [IndexedOrientation] = []
        for s in [start, start.rotatedInPlane(by: .pi)] {
            let fitted = solver.fit(s, scale: scale, spots: q)
            let o1 = fitted?.orientation ?? s, s1 = fitted?.scale ?? scale
            guard let tilted = solver.refineTilt(o1, scale: s1, spots: q, intensities: intensity, settings: settings.tilt)
            else { continue }
            var o = tilted.orientation, sc = s1, c = tilted.agreement
            if let again = solver.fit(tilted.orientation, scale: s1, spots: q) {
                let c3 = solver.intensityAgreement(of: again.orientation, scale: again.scale, spots: q, intensities: intensity,
                                                   intensityPower: settings.tilt.intensityPower)
                if c3 >= c { o = again.orientation; sc = again.scale; c = c3 }
            }
            guard let positions = solver.fit(o, scale: sc, spots: q) else { continue }
            out.append(IndexedOrientation(orientation: o, scale: sc, agreement: c, positions: positions, tiltSpread: .nan,
                                          source: source))
        }
        return out
    }

    private func withSpread(_ r: IndexedOrientation, q: [SIMD2<Double>], intensity: [Double]) -> IndexedOrientation {
        var tilt = settings.tilt
        tilt.ambiguityMargin = settings.spreadMargin
        return IndexedOrientation(orientation: r.orientation, scale: r.scale, agreement: r.agreement, positions: r.positions,
                           tiltSpread: solver.tiltSpread(of: r.orientation, scale: r.scale, spots: q, intensities: intensity,
                                                         settings: tilt),
                           source: r.source)
    }

    /// The winner turned half about each candidate mirror line of the pattern (along and across its 8 innermost spots), kept
    /// where the positions still fit within `lookAlikeMargin` of the winner's F_β.
    private func lookAlikeStarts(of winner: IndexedOrientation, q: [SIMD2<Double>]) -> [ZoneAxisCandidate] {
        let order = q.indices.sorted { simd_length(q[$0]) < simd_length(q[$1]) }.prefix(8)
        var lines: [Double] = []
        let halfDegree = 0.5 * Double.pi / 180
        for i in order {
            let a = atan2(q[i].y, q[i].x)
            for x in [a, a + .pi / 2] {
                let line = x - .pi * (x / .pi).rounded(.down)       // in [0, π)
                if lines.allSatisfy({ abs(remainder(line - $0, .pi)) > halfDegree }) { lines.append(line) }
            }
        }
        var out: [ZoneAxisCandidate] = []
        for line in lines {
            // the half turn about the lab axis (cos α, sin α, 0): a mirror of the pattern across that line, z to −z
            let c = cos(2 * line), s = sin(2 * line)
            let turn = simd_double3x3(rows: [SIMD3(c, s, 0), SIMD3(s, -c, 0), SIMD3(0, 0, -1)])
            guard let start = try? CrystalOrientation(matrix: winner.orientation.matrix * turn.transpose),
                  let f = solver.fit(start, scale: winner.scale, spots: q),
                  f.agreement >= winner.positions.agreement - settings.lookAlikeMargin else { continue }
            out.append(f)
        }
        return out
    }

    /// The library's matrix (built from a Float in-plane angle) made orthonormal again (Gram-Schmidt on the columns, the third the cross product).
    static func orthonormalized(_ m: simd_double3x3) -> simd_double3x3 {
        let x = simd_normalize(m.columns.0)
        let y = simd_normalize(m.columns.1 - simd_dot(m.columns.1, x) * x)
        return simd_double3x3(columns: (x, y, simd_cross(x, y)))
    }
}
