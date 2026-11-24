//
//  PhaseFinder.swift
//  Role: candidate phases ranked by how well each explains one measured pattern (DiffractoKit S11, #166). Ranked evidence
//        with ties, never a verdict.
//
//  Lineage: search-match on d-spacings (Hanawalt, Rinn and Frevel, Ind. Eng. Chem. Anal. Ed. 10 (1938) 457) as done for
//  electron rings in ProcessDiffraction (Lábár, Ultramicroscopy 103 (2005) 237) and CrysTBox (Klinger and Jäger, J. Appl.
//  Cryst. 48 (2015) 2012), and multi-phase matching for spots as in py4DSTEM. No package in the pinned stack ranks phases
//  from a list, so there is no Python reference: every number per candidate is one the kit already computes for one phase
//  (the kinematic ring list of `kinematicReflections`; the zone-axis solver's scale window, tolerance and agreement), and the
//  order is the rule written below.
//
//  RINGS. For each candidate the scale s (measured |q| over the candidate's |g|) is fitted inside the solver's window: every
//  s = q_i / g_j inside it is tried, the one explaining the most rings with the least rms misfit kept, then refined by least
//  squares over its pairs (kept only if it explains as many). A measured ring is explained when a candidate ring lies within
//  max(tolerance · q, floor) of it. A candidate ring is strong when its Σ|F|² is at least `strongFraction` of the
//  candidate's strongest ring within the pattern's reach (the rings a pattern of this extent could show, not the whole list); the strong rings that fall inside the measured range (widened by the tolerance) are expected, and
//  seen when a measured ring lies within the tolerance. Precision = explained / measured, recall = seen / expected (1 when
//  none is expected), agreement = their F1, as in the zone-axis solver's agreement with β = 1.
//  SPOTS. Each candidate is solved with `ZoneAxisSolver` (its own window, tolerance and agreement); the numbers are its best
//  candidate's.
//  ORDER. Agreement, highest first; then (rings) the rings explained and the rms misfit or distance (spots), then the position in the input.
//  TIES. Rings: the best candidate not yet placed leads a tie, and every later candidate that one more ring of evidence
//  (one more measured ring explained, or one more strong ring seen) would bring up to the leader's agreement joins it; a
//  tie is listed together at its leader's place. One ring is the grain of the evidence: a pattern loses rings (a fifth in
//  the planted set) and kinematic intensities only roughly say which rings are strong, so a lead of less than one ring is
//  not a lead. Spots: the same explained count and an agreement within the solver's tie margin of the tie's first. With a
//  free scale two phases with nearly equal lattices (Al and Au, 0.7 % apart) tie: positions cannot tell them apart,
//  chemistry can.
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoGeometry

public enum PhaseFinder {

    public struct RingSettings: Sendable, Equatable {
        /// A ring is explained within this fraction of its |q| ... Default 0.01 (the zone-axis solver's).
        public var tolerance: Double
        /// ... or within this many Å⁻¹, whichever is larger. Default 0.004 (the solver's).
        public var floor: Double
        /// The scale may be off by up to this fraction. Default 0.03 (the solver's); 0 holds the calibration fixed.
        public var scaleWindow: Double
        /// A candidate ring is strong at this fraction of the candidate's strongest Σ|F|². Default 0.2.
        public var strongFraction: Double

        public init(tolerance: Double = 0.01, floor: Double = 0.004, scaleWindow: Double = 0.03, strongFraction: Double = 0.2) {
            self.tolerance = tolerance
            self.floor = floor
            self.scaleWindow = scaleWindow
            self.strongFraction = strongFraction
        }

        func validated() throws(PhaseFinderError) {
            guard tolerance.isFinite, tolerance >= 0, floor.isFinite, floor >= 0, tolerance > 0 || floor > 0 else {
                throw .invalidSettings("tolerance and floor must be finite, at least 0, and not both 0")
            }
            guard scaleWindow.isFinite, scaleWindow >= 0, scaleWindow < 0.5 else {
                throw .invalidSettings("scaleWindow must be at least 0 and below 0.5")
            }
            guard strongFraction.isFinite, strongFraction > 0, strongFraction <= 1 else {
                throw .invalidSettings("strongFraction must be above 0 and at most 1")
            }
        }
    }

    /// One candidate's evidence on a ring pattern.
    public struct RingScore: Sendable, Equatable {
        /// The candidate's position in the input.
        public let candidate: Int
        /// Measured rings a candidate ring explains, of `measured`.
        public let explained: Int
        public let measured: Int
        /// The candidate's strong rings inside the measured range, and how many of them were seen.
        public let strongExpected: Int
        public let strongSeen: Int
        public var strongMissing: Int { strongExpected - strongSeen }
        /// Measured |q| over the candidate's |g| (1 = the calibration as given).
        public let scale: Double
        /// rms of (q − s g) / q over the explained rings; 0 when none is explained.
        public let misfit: Double
        /// F1 of explained / measured and seen / expected, 0 to 1.
        public let agreement: Double
        /// Candidates with the same tie number are one tie, listed together; numbered from 0 in rank order.
        public let tie: Int
        /// Why the candidate could not be scored (its ring list was refused); it then ranks last with zero counts.
        public let refusal: String?
    }

    /// One candidate's evidence on a spot pattern.
    public struct SpotScore: Sendable, Equatable {
        public let candidate: Int
        /// Measured spots the best orientation explains, of `measured`.
        public let explained: Int
        public let measured: Int
        /// The solver's agreement (F_β of recall and precision), 0 to 1.
        public let agreement: Double
        /// The solver's scale and rms distance of the explained spots (Å⁻¹).
        public let scale: Double
        public let rmsDistance: Double
        /// The best orientation; nil when the solver found none.
        public let best: ZoneAxisCandidate?
        public let tie: Int
        public let refusal: String?
    }

    /// The candidates ranked on measured ring positions `q` (Å⁻¹, no 2π, any order). Throws for an empty or non-positive
    /// list and for settings out of range; a candidate whose ring list is refused is scored last with its refusal.
    public static func rank(rings q: [Double], candidates: [CrystalStructure], settings: RingSettings = .init())
        throws(PhaseFinderError) -> [RingScore] {
        try settings.validated()
        guard !q.isEmpty, q.allSatisfy({ $0.isFinite && $0 > 0 }) else {
            throw .invalidInput("the ring positions must be one or more positive numbers of Å⁻¹")
        }
        let measured = q.sorted()
        // Far enough for the outermost ring at the smallest scale: (q + window(q)) / (1 − w), with a margin.
        let last = measured.last!
        let kMax = max(last * (1 + settings.scaleWindow + settings.tolerance) + settings.floor + 0.01,
                       (last + max(settings.tolerance * last, settings.floor)) / (1 - settings.scaleWindow) * (1 + 1e-6) + 0.01)
        var raw: [(candidate: Int, counts: [Int], scale: Double, misfit: Double, agreement: Double, refusal: String?)] = []
        for (index, structure) in candidates.enumerated() {
            let families: [KinematicFamily]
            do {
                families = try structure.kinematicReflections(kMax: kMax, retainingAll: false).families
            } catch {
                raw.append((index, [0, 0, 0], 1, 0, 0, String(describing: error)))
                continue
            }
            let s = ringScore(measured: measured, families: families, settings: settings)
            raw.append((index, [s.explained, s.seen, s.expected], s.scale, s.misfit, s.agreement, nil))
        }
        let order = raw.indices.sorted { a, b in
            let x = raw[a], y = raw[b]
            if (x.refusal == nil) != (y.refusal == nil) { return x.refusal == nil }
            if x.agreement != y.agreement { return x.agreement > y.agreement }
            if x.counts[0] != y.counts[0] { return x.counts[0] > y.counts[0] }
            if x.misfit != y.misfit { return x.misfit < y.misfit }
            return x.candidate < y.candidate
        }
        // Ties: the first candidate not yet placed leads; every later one that one more ring of evidence (one more
        // measured ring explained, or one more strong ring seen) would bring up to the leader's agreement joins it.
        var placed = Set<Int>()
        var scores: [RingScore] = []
        var tie = -1
        for leader in order where !placed.contains(leader) {
            tie += 1
            let l = raw[leader]
            let members = l.refusal != nil ? [leader] : order.filter { i in
                guard !placed.contains(i), raw[i].refusal == nil else { return false }
                return i == leader || oneRingMore(raw[i].counts, measured: measured.count) >= l.agreement - 1e-12
            }
            for i in members {
                placed.insert(i)
                let r = raw[i]
                scores.append(RingScore(candidate: r.candidate, explained: r.counts[0], measured: measured.count,
                                        strongExpected: r.counts[2], strongSeen: r.counts[1], scale: r.scale, misfit: r.misfit,
                                        agreement: r.agreement, tie: tie, refusal: r.refusal))
            }
        }
        return scores
    }

    /// The candidates ranked on measured spot positions (Å⁻¹ in the pattern plane, direct beam removed) at `volts`, each
    /// solved with `ZoneAxisSolver(structure:volts:kMax:settings:)`. A candidate the solver cannot be built for is scored
    /// last with its refusal.
    public static func rank(spots: [SIMD2<Double>], candidates: [CrystalStructure], volts: Double, kMax: Double = 2,
                            solverSettings: ZoneAxisSolver.Settings = .init()) throws(PhaseFinderError) -> [SpotScore] {
        guard spots.count >= 2, spots.allSatisfy({ $0.x.isFinite && $0.y.isFinite }) else {
            throw .invalidInput("a spot pattern needs two or more finite spots")
        }
        var raw: [(candidate: Int, explained: Int, agreement: Double, scale: Double, rms: Double, best: ZoneAxisCandidate?,
                   refusal: String?)] = []
        for (index, structure) in candidates.enumerated() {
            do {
                let solver = try ZoneAxisSolver(structure: structure, volts: volts, kMax: kMax, settings: solverSettings)
                let best = solver.solve(spots).candidates.first
                raw.append((index, best?.matched ?? 0, best?.agreement ?? 0, best?.scale ?? 1, best?.rmsDistance ?? 0, best, nil))
            } catch {
                raw.append((index, 0, 0, 1, 0, nil, String(describing: error)))
            }
        }
        let order = raw.indices.sorted { a, b in
            let x = raw[a], y = raw[b]
            if (x.refusal == nil) != (y.refusal == nil) { return x.refusal == nil }
            if x.agreement != y.agreement { return x.agreement > y.agreement }
            if x.explained != y.explained { return x.explained > y.explained }
            if x.rms != y.rms { return x.rms < y.rms }
            return x.candidate < y.candidate
        }
        var scores: [SpotScore] = []
        var tie = -1
        var lead: (explained: Int, agreement: Double)? = nil
        for i in order {
            let r = raw[i]
            if r.refusal != nil || lead == nil || r.explained != lead!.explained
                || lead!.agreement - r.agreement > solverSettings.tieMargin {
                tie += 1
                lead = r.refusal == nil ? (r.explained, r.agreement) : nil
            }
            scores.append(SpotScore(candidate: r.candidate, explained: r.explained, measured: spots.count, agreement: r.agreement,
                                    scale: r.scale, rmsDistance: r.rms, best: r.best, tie: tie, refusal: r.refusal))
        }
        return scores
    }

    // MARK: - Rings

    /// F1 of precision e / m and recall seen / expected (1 when nothing is expected).
    static func f1(explained e: Int, seen: Int, expected: Int, measured m: Int) -> Double {
        let p = Double(e) / Double(m), r = expected > 0 ? Double(seen) / Double(expected) : 1
        return p + r > 0 ? 2 * p * r / (p + r) : 0
    }

    /// The best agreement one more ring of evidence would give counts [explained, seen, expected].
    static func oneRingMore(_ c: [Int], measured m: Int) -> Double {
        max(f1(explained: min(c[0] + 1, m), seen: c[1], expected: c[2], measured: m),
            f1(explained: c[0], seen: min(c[1] + 1, c[2]), expected: c[2], measured: m))
    }

    struct RingCount {
        var explained = 0, seen = 0, expected = 0
        var scale = 1.0, misfit = 0.0, agreement = 0.0
    }

    /// The score of one candidate's `families` on the sorted `measured` rings (internal for the tests).
    static func ringScore(measured: [Double], families: [KinematicFamily], settings: RingSettings) -> RingCount {
        let g = families.map(\.family.g).filter { $0 > 0 }
        guard !g.isEmpty else { return RingCount() }
        func window(_ q: Double) -> Double { max(settings.tolerance * q, settings.floor) }
        /// Explained count, sum of squared relative residuals, and the pairs (q, g).
        func explain(_ s: Double) -> (count: Int, sse: Double, pairs: [(Double, Double)]) {
            var count = 0, sse = 0.0
            var pairs: [(Double, Double)] = []
            for q in measured {
                var best = Double.infinity, bestG = 0.0
                for x in g { let d = abs(q - s * x); if d < best { best = d; bestG = x } }
                if best <= window(q) {
                    count += 1
                    sse += (best / q) * (best / q)
                    pairs.append((q, bestG))
                }
            }
            return (count, sse, pairs)
        }
        func better(_ a: (count: Int, sse: Double), _ sa: Double, than b: (count: Int, sse: Double), _ sb: Double) -> Bool {
            if a.count != b.count { return a.count > b.count }
            let ra = a.count > 0 ? a.sse / Double(a.count) : 0, rb = b.count > 0 ? b.sse / Double(b.count) : 0
            if ra != rb { return ra < rb }
            return abs(sa - 1) < abs(sb - 1)
        }
        var trials: [Double] = [1]
        if settings.scaleWindow > 0 {
            for q in measured { for x in g { let s = q / x; if abs(s - 1) <= settings.scaleWindow { trials.append(s) } } }
        }
        var scale = 1.0
        var bestFit = explain(1)
        for s in trials.dropFirst() {
            let fit = explain(s)
            if better((fit.count, fit.sse), s, than: (bestFit.count, bestFit.sse), scale) { bestFit = fit; scale = s }
        }
        if settings.scaleWindow > 0, !bestFit.pairs.isEmpty {
            // Least squares q ≈ s g over the pairs, clamped to the window.
            let num = bestFit.pairs.reduce(0) { $0 + $1.0 * $1.1 }, den = bestFit.pairs.reduce(0) { $0 + $1.1 * $1.1 }
            let refined = min(max(num / den, 1 - settings.scaleWindow), 1 + settings.scaleWindow)
            let fit = explain(refined)
            if fit.count >= bestFit.count && (fit.count > bestFit.count || fit.sse <= bestFit.sse) { bestFit = fit; scale = refined }
        }
        var count = RingCount()
        count.explained = bestFit.count
        count.scale = scale
        count.misfit = bestFit.count > 0 ? (bestFit.sse / Double(bestFit.count)).squareRoot() : 0
        let strongest = families.map(\.intensity).max() ?? 0
        let low = measured.first! - window(measured.first!), high = measured.last! + window(measured.last!)
        for f in families where f.intensity >= settings.strongFraction * strongest {
            let p = scale * f.family.g
            guard p >= low, p <= high else { continue }
            count.expected += 1
            if measured.contains(where: { abs($0 - p) <= window($0) }) { count.seen += 1 }
        }
        count.agreement = f1(explained: count.explained, seen: count.seen, expected: count.expected, measured: measured.count)
        return count
    }
}

public enum PhaseFinderError: Error, Equatable, Sendable, CustomStringConvertible {
    case invalidInput(String)
    case invalidSettings(String)

    public var description: String {
        switch self {
        case .invalidInput(let s): return "The pattern cannot be ranked: \(s)."
        case .invalidSettings(let s): return "The phase finder settings are not valid: \(s)."
        }
    }
}

extension PhaseFinderError: LocalizedError {
    public var errorDescription: String? { description }
}
