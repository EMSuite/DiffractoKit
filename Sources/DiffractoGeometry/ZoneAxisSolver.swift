//
//  ZoneAxisSolver.swift
//  Role: the zone axis and full orientation of ONE measured spot pattern (DiffractoKit S5): from the spot positions
//        (Å⁻¹, the direct beam removed) and a structure, the orientations that explain them, ranked by how well each one's
//        kinematic pattern (`SpotPatternSimulator`, S2) agrees with the measured spots, with the scale, the in-plane turn,
//        every spot's hkl, and a verdict: solved, a tie, or refused.
//
//  Method (indexing by lattice-vector pairs, the textbook route of Williams & Carter, Transmission Electron Microscopy,
//  2nd ed. 2009, ch. 18, and of the pair libraries of diffsims 0.7.0 `VectorLibraryGenerator` and pyxem's vector
//  indexation):
//    1. Seeds: the shortest measured vectors, one per direction (a Friedel mate or a higher order of a row adds nothing).
//       Fewer than two directions is one row of spots, which fixes no zone axis: refused.
//    2. Every ordered pair of seeds (q_i, q_k) at more than 14.5° to each other picks reflection pairs (g_a, g_b) whose
//       lengths match within the scale window plus the tolerance and whose 3-D angle and length ratio match the measured
//       ones; g_a is one member of its symmetry orbit. Readings of the measured pair in 3-D: both spots in the zero layer
//       (z = 0), both on the Ewald sphere (z = λg²/2), or q_i in the zero layer and q_k above or below it (its z from the
//       angle between g_a and g_b). When the best of these explains less than 95 % of the pattern, a second pass reads
//       q_k's height from its length (g_a and g_b near perpendicular) and both spots in the same or opposite first-order
//       layers (the common scale from the angle): a zero layer of one row or none. Each reflection pair fixes a full
//       orientation (TRIAD).
//    3. A candidate within `snapAngle` of a zone axis [u v w] (|u|, |v|, |w| ≤ `maxIndex`) is turned onto it.
//    4. Each candidate predicts its pattern (excitation error ≤ the reach); the scale and in-plane turn are fitted by the
//       closed-form similarity fit on the matched spots (twice), then scored: recall R = the share of measured spots within
//       the tolerance of a predicted one; precision P = the share of the predicted intensity (|F|² times the excitation
//       envelope) inside the measured |q| range, over predicted spots at least as strong as the weakest one matched (the
//       pattern's detection limit), that a measured spot sits on (1 for clicked spots, a subset by choice); agreement
//       F_β = (1 + β²) P R / (β² P + R) with β = `recallWeight`.
//    5. Ranked by agreement, then matched spots, then the matched spots' share of the strongest predicted intensity, then
//       the rms distance. Orientations within 0.2° up to the crystal's proper rotations are one; one that is the half
//       turn about the beam of a better one (the same zero-layer spots for a centrosymmetric pattern) is merged into it and
//       flags it when it scores within the tie margin.
//    6. Verdict over the best 20 distinct orientations: refused when the best agreement is below `minimumAgreement`; a
//       tie when another scores within `tieMargin` of the best; solved otherwise.
//  The defaults were set on the planted set (P43): 8 phases, 1070 patterns, none wrong.
//  No Python function does this step alone (py4DSTEM's match_single_pattern is a template match on a 2° grid, S6); the pins
//  are planted patterns with a known orientation (docs/parity.md P43).
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering

/// Why the solver returned no orientation.
public enum ZoneAxisRefusal: Sendable, Equatable, CustomStringConvertible {
    /// Fewer than two spots.
    case tooFewSpots(Int)
    /// The spots lie on one row through the direct beam: any zone axis that contains that row explains them.
    case oneRow
    /// No reflection pair of the structure matches the measured lengths and angle.
    case noCandidate
    /// The best orientation explains too little of the pattern (its agreement, 0 to 1).
    case lowAgreement(Double)

    public var description: String {
        switch self {
        case .tooFewSpots(let n): return "A zone axis needs at least two spots in two directions; \(n) were given."
        case .oneRow: return "The spots lie on one row through the direct beam, which fixes no zone axis."
        case .noCandidate: return "No pair of the structure's reflections matches the measured lengths and angle."
        case .lowAgreement(let f): return String(format: "The best orientation explains too little of the pattern (agreement %.2f).", f)
        }
    }
}

public enum ZoneAxisSolverError: Error, Equatable, Sendable, CustomStringConvertible {
    /// The structure's reflections were refused (see `SpotPatternError`).
    case pattern(SpotPatternError)
    /// A setting is out of its range; the text names it.
    case invalidSettings(String)

    public var description: String {
        switch self {
        case .pattern(let e): return e.description
        case .invalidSettings(let s): return "The zone-axis solver settings are not valid: \(s)."
        }
    }
}

extension ZoneAxisSolverError: LocalizedError {
    public var errorDescription: String? { description }
}

/// One orientation that explains a measured pattern.
public struct ZoneAxisCandidate: Sendable, Equatable {
    /// The full orientation (py4DSTEM's convention, `CrystalOrientation`).
    public let orientation: CrystalOrientation
    /// The zone axis nearest the beam (|u|, |v|, |w| ≤ the solver's `maxIndex`) and the angle to it, radians.
    public let zoneAxis: ZoneAxis
    public let zoneAxisAngle: Double
    /// The measured |q| over the predicted |g|: the calibration error of the pattern times the lattice's (1 = exact).
    public let scale: Double
    /// F_β of recall and precision (β = the solver's `recallWeight`), 0 to 1.
    public let agreement: Double
    /// The share of measured spots explained.
    public let recall: Double
    /// The share of predicted intensity in the measured range that a measured spot sits on (1 for clicked spots).
    public let precision: Double
    /// The number of measured spots explained.
    public let matched: Int
    /// The rms distance of the explained spots from their predicted positions, Å⁻¹.
    public let rmsDistance: Double
    /// Each measured spot's hkl, in the order given; nil for a spot this orientation does not explain.
    public let assignments: [MillerIndex?]
    /// The half turn about the beam scores within the tie margin and is not a symmetry of the crystal: the pattern cannot
    /// tell the two apart (a centrosymmetric zero layer). The other is `orientation.rotatedInPlane(by: .pi)`.
    public let halfTurnAmbiguous: Bool
}

/// The solver's answer for one pattern.
public struct ZoneAxisSolution: Sendable, Equatable {
    public enum Verdict: Sendable, Equatable {
        /// One orientation stands out.
        case solved
        /// The first `count` candidates (2 or more) score within the tie margin: the pattern cannot tell them apart.
        case tie(count: Int)
        case refused(ZoneAxisRefusal)
    }
    public let verdict: Verdict
    /// The distinct orientations found, best first (at most the `top` asked for); may be non-empty when refused for a low
    /// agreement.
    public let candidates: [ZoneAxisCandidate]
    /// The best candidate when solved, nil otherwise.
    public var solved: ZoneAxisCandidate? {
        if case .solved = verdict { return candidates.first }
        return nil
    }
}

/// The zone axis and full orientation of one spot pattern of a known structure. Build it once per structure and voltage.
public struct ZoneAxisSolver: Sendable {

    public struct Settings: Sendable, Equatable {
        /// A spot is explained within this fraction of its |q| ... Default 0.01 (1 %).
        public var tolerance: Double
        /// ... or within this many Å⁻¹, whichever is larger. Default 0.004.
        public var floor: Double
        /// The pattern's scale (calibration times lattice) may be off by up to this fraction. Default 0.03.
        public var scaleWindow: Double
        /// The excitation envelope and the reach of the predicted spots (σ = 0.02 Å⁻¹ to 3σ by default).
        public var pattern: SpotPatternSimulator.Settings
        /// How many of the shortest measured directions seed reflection pairs. Default 8.
        public var seeds: Int
        /// Zone axes up to this index are snapped to and reported. Default 3.
        public var maxIndex: Int
        /// A candidate within this angle (radians) of a zone axis is turned onto it. Default 3°.
        public var snapAngle: Double
        /// The largest lab z of a higher-layer seed, Å⁻¹. Default 0.15.
        public var higherLayerReach: Double
        /// A different orientation within this much of the best agreement makes a tie. Default 0.05.
        public var tieMargin: Double
        /// A best agreement below this is refused. Default 0.7.
        public var minimumAgreement: Double
        /// β of the agreement F_β = (1 + β²) P R / (β² P + R): recall weighs β times precision. Default 1.5 (a missing
        /// weak spot is commoner than a false one).
        public var recallWeight: Double

        public init(tolerance: Double = 0.01, floor: Double = 0.004, scaleWindow: Double = 0.03,
                    pattern: SpotPatternSimulator.Settings = .init(), seeds: Int = 8, maxIndex: Int = 3,
                    snapAngle: Double = 3 * .pi / 180, higherLayerReach: Double = 0.15, tieMargin: Double = 0.05,
                    minimumAgreement: Double = 0.7, recallWeight: Double = 1.5) {
            self.tolerance = tolerance
            self.floor = floor
            self.scaleWindow = scaleWindow
            self.pattern = pattern
            self.seeds = seeds
            self.maxIndex = maxIndex
            self.snapAngle = snapAngle
            self.higherLayerReach = higherLayerReach
            self.tieMargin = tieMargin
            self.minimumAgreement = minimumAgreement
            self.recallWeight = recallWeight
        }

        func validated() throws(ZoneAxisSolverError) {
            func positive(_ v: Double, _ name: String) throws(ZoneAxisSolverError) {
                guard v.isFinite, v > 0 else { throw .invalidSettings("\(name) must be a finite number above 0") }
            }
            try positive(tolerance, "tolerance")
            try positive(floor, "floor")
            try positive(higherLayerReach, "higherLayerReach")
            try positive(recallWeight, "recallWeight")
            guard scaleWindow.isFinite, scaleWindow >= 0, scaleWindow < 0.5 else {
                throw .invalidSettings("scaleWindow must be at least 0 and below 0.5")
            }
            guard snapAngle.isFinite, snapAngle >= 0 else { throw .invalidSettings("snapAngle must be at least 0") }
            guard seeds >= 2 else { throw .invalidSettings("seeds must be at least 2") }
            guard (1...12).contains(maxIndex) else { throw .invalidSettings("maxIndex must be 1 to 12") }
            guard tieMargin.isFinite, tieMargin >= 0 else { throw .invalidSettings("tieMargin must be at least 0") }
            guard minimumAgreement.isFinite, (0...1).contains(minimumAgreement) else {
                throw .invalidSettings("minimumAgreement must be 0 to 1")
            }
            do { try pattern.validated() } catch { throw .pattern(error) }
        }
    }

    public let structure: CrystalStructure
    public let volts: Double
    public let settings: Settings
    /// λ, Å.
    public let wavelength: Double

    let hkl: [MillerIndex]
    let g: [SIMD3<Double>]
    let gLength: [Double]
    let intensity: [Double]
    /// One member of each orbit under the proper rotations (the lexicographically largest hkl).
    private let canonical: [Bool]
    /// The crystal's proper rotations in its Cartesian frame.
    let rotations: [simd_double3x3]

    /// The solver for `structure` at `volts` (200 kV is 200_000), with the reflections up to `kMax` Å⁻¹ (they must reach
    /// past the longest measured spot).
    public init(structure: CrystalStructure, volts: Double, kMax: Double = 2, settings: Settings = .init())
        throws(ZoneAxisSolverError) {
        try settings.validated()
        let simulator: SpotPatternSimulator
        do { simulator = try SpotPatternSimulator(structure: structure, kMax: kMax) } catch { throw .pattern(error) }
        let lambda = ElectronWavelength.angstroms(volts: volts)
        guard lambda.isFinite, lambda > 0 else { throw .pattern(.invalidVoltage(volts)) }
        self.structure = structure
        self.volts = volts
        self.settings = settings
        self.wavelength = lambda
        let r = simulator.reflections.map(\.structureFactor)
        let hkls = r.map(\.hkl)
        hkl = hkls
        g = r.map(\.g)
        gLength = r.map(\.gLength)
        intensity = r.map(\.intensity)

        // W acts on fractional coordinates (x → W x); on the Cartesian frame it is C = Aᵀ W A⁻ᵀ (Aᵀ: columns a₁, a₂, a₃),
        // and on a reflection h → h W (the orbit of h).
        // A proper rotation of the cell's frame; a cell that does not carry W is skipped (`DoubleTiltHolder`, the one copy).
        let carried = DoubleTiltHolder.cartesianProperRotations(of: structure)
        var ops = carried.map(\.cartesian)
        var integer = carried.map(\.integer)
        if ops.isEmpty { ops = [matrix_identity_double3x3]; integer = [.identity] }
        rotations = ops
        canonical = hkls.map { h in
            var best = (h.h, h.k, h.l)
            for w in integer {
                let m = h.h &* w.rows[0] &+ h.k &* w.rows[1] &+ h.l &* w.rows[2]
                if (m.x, m.y, m.z) > best { best = (m.x, m.y, m.z) }
            }
            return best == (h.h, h.k, h.l)
        }
    }

    /// The orientations that explain `spots` (Å⁻¹ in the pattern plane, direct beam removed, any order), best first, at most
    /// `top`. `clicked`: the spots are a subset a person picked (all in the zero layer), so a predicted spot without a
    /// measured one does not count against an orientation.
    public func solve(_ spots: [SIMD2<Double>], clicked: Bool = false, top: Int = 5) -> ZoneAxisSolution {
        let q = spots.filter { $0.x.isFinite && $0.y.isFinite && simd_length($0) > 0 }
        guard q.count >= 2 else { return ZoneAxisSolution(verdict: .refused(.tooFewSpots(q.count)), candidates: []) }
        let qn = q.map { simd_length($0) }
        let order = q.indices.sorted { qn[$0] < qn[$1] }
        var seeds: [Int] = []
        for i in order {
            let u = q[i] / qn[i]
            if seeds.allSatisfy({ abs(u.x * q[$0].y - u.y * q[$0].x) / qn[$0] > 0.1 }) { seeds.append(i) }
            if seeds.count >= settings.seeds { break }
        }
        guard seeds.count >= 2 else { return ZoneAxisSolution(verdict: .refused(.oneRow), candidates: []) }

        func pass(_ modes: [Mode]) -> [Scored] {
            var out: [Scored] = []
            for i in seeds { for k in seeds where k != i {
                for c in pairCandidates(q[i], q[k], firstIsEarlier: i < k, modes: modes) {
                    if let s = score(snapped(c.orientation) ?? c.orientation, q, qn, scale: c.scale, clicked: clicked) {
                        out.append(s)
                    }
                }
            } }
            return out
        }
        // Clicked spots are zero-layer spots by choice. A full pattern whose best reading explains little may have no two
        // seeds in the zero layer: the second pass reads pairs with both spots in higher layers.
        var scored = pass(clicked ? [.zero] : Mode.first)
        if !clicked, (scored.map(\.agreement).max() ?? 0) < Self.secondPassBelow { scored += pass(Mode.second) }
        guard !scored.isEmpty else { return ZoneAxisSolution(verdict: .refused(.noCandidate), candidates: []) }
        scored.sort { a, b in
            let fa = (a.agreement * 1e6).rounded(), fb = (b.agreement * 1e6).rounded()
            if fa != fb { return fa > fb }
            if a.matched != b.matched { return a.matched > b.matched }
            let sa = (a.support * 1e6).rounded(), sb = (b.support * 1e6).rounded()
            if sa != sb { return sa > sb }
            return a.rms < b.rms
        }

        // Merge orientations equal up to the crystal's proper rotations; the half turn about the beam joins its better twin.
        let turn = simd_double3x3(diagonal: SIMD3(-1, -1, 1))
        var merged: [(Scored, Bool)] = []
        for s in scored {
            var duplicate = false
            for j in merged.indices {
                let m = merged[j].0
                if misorientation(s.matrix, m.matrix) < Self.sameAngle { duplicate = true; break }
                if misorientation(s.matrix * turn, m.matrix) < Self.sameAngle {
                    if s.agreement >= m.agreement - settings.tieMargin { merged[j].1 = true }
                    duplicate = true; break
                }
            }
            if !duplicate {
                merged.append((s, false))
                if merged.count >= max(top, Self.verdictDepth) { break }
            }
        }
        let out = merged.compactMap { s, twin in candidate(s, halfTurnAmbiguous: twin) }
        guard let best = out.first?.agreement else {
            return ZoneAxisSolution(verdict: .refused(.noCandidate), candidates: [])
        }
        if best < settings.minimumAgreement {
            return ZoneAxisSolution(verdict: .refused(.lowAgreement(best)), candidates: Array(out.prefix(max(top, 1))))
        }
        let tied = out.prefix { $0.agreement >= best - settings.tieMargin }.count
        let shown = Array(out.prefix(max(top, tied, 1)))
        return ZoneAxisSolution(verdict: tied > 1 ? .tie(count: tied) : .solved, candidates: shown)
    }

    /// `orientation` at `scale` with its in-plane turn and scale fitted to `spots` (Å⁻¹, as `solve` takes them): the scoring
    /// step of `solve` for an orientation found elsewhere (a template match, a refined tilt). Nil when fewer than two spots
    /// sit on a predicted spot. The scale stays within the settings' window of 1 (plus 0.005); `halfTurnAmbiguous` is false
    /// (only `solve` compares the half turn).
    public func fit(_ orientation: CrystalOrientation, scale: Double = 1, spots: [SIMD2<Double>], clicked: Bool = false)
        -> ZoneAxisCandidate? {
        let q = spots.filter { $0.x.isFinite && $0.y.isFinite && simd_length($0) > 0 }
        guard q.count >= 2, scale.isFinite, scale > 0,
              let s = score(orientation.matrix, q, q.map { simd_length($0) }, scale: scale, clicked: clicked) else { return nil }
        return candidate(s, halfTurnAmbiguous: false)
    }

    private func candidate(_ s: Scored, halfTurnAmbiguous: Bool) -> ZoneAxisCandidate? {
        guard let o = try? CrystalOrientation(matrix: s.matrix) else { return nil }
        let zone = o.nearestZoneAxis(cell: structure.cell, maxIndex: settings.maxIndex)
        return ZoneAxisCandidate(orientation: o, zoneAxis: zone.map { ZoneAxis($0.uvw.x, $0.uvw.y, $0.uvw.z) } ?? ZoneAxis(0, 0, 1),
                                 zoneAxisAngle: zone?.angle ?? .nan, scale: s.scale, agreement: s.agreement,
                                 recall: s.recall, precision: s.precision, matched: s.matched, rmsDistance: s.rms,
                                 assignments: s.assignments.map { $0.map { hkl[$0] } }, halfTurnAmbiguous: halfTurnAmbiguous)
    }

    /// Orientations within this angle (radians) are one (0.2°).
    static let sameAngle = 0.2 * Double.pi / 180
    /// The distinct orientations ranked before the verdict, so that a tie is seen past the `top` asked for.
    static let verdictDepth = 20
    /// The second pass of pair readings runs when the best agreement of the first is below this.
    static let secondPassBelow = 0.95

    // MARK: - candidates

    /// How a measured pair is read in 3-D (z: the lab height of each spot, from the reflections' lengths and angle).
    private enum Mode {
        /// Both spots in the zero layer (z = 0).
        case zero
        /// Both spots on the Ewald sphere (z = λg²/2).
        case ewald
        /// q_i in the zero layer, q_k in a higher one; its height from the 3-D angle (sign ±1).
        case layer(Double)
        /// q_i in the zero layer, q_k in a higher one; its height from its length.
        case length(Double)
        /// Both in the same first-order layer.
        case sameLayer(Double)
        /// Both in opposite first-order layers.
        case oppositeLayer(Double)

        static let first: [Mode] = [.zero, .ewald, .layer(1), .layer(-1)]
        static let second: [Mode] = [.length(1), .length(-1), .sameLayer(1), .sameLayer(-1), .oppositeLayer(1),
                                     .oppositeLayer(-1)]
    }

    /// The orientations one ordered pair of measured vectors fixes under `modes` (zero and ewald only once per pair).
    private func pairCandidates(_ qi: SIMD2<Double>, _ qk: SIMD2<Double>, firstIsEarlier: Bool, modes: [Mode])
        -> [(orientation: simd_double3x3, scale: Double)] {
        let s = settings
        let ni = simd_length(qi), nk = simd_length(qk)
        let sine = (qi.x * qk.y - qi.y * qk.x) / (ni * nk)
        guard abs(sine) >= 0.25 else { return [] }
        let window = s.scaleWindow + s.tolerance
        let ia = g.indices.filter { canonical[$0] && abs(gLength[$0] / ni - 1) <= window }
        let ib = g.indices.filter { abs(gLength[$0] / nk - 1) <= window }
        guard !ia.isEmpty, !ib.isEmpty else { return [] }
        let angleTolerance = 2 * max(s.tolerance, s.floor / min(ni, nk))
        let dot2 = simd_dot(qi, qk)
        let t = dot2 / ni
        var out: [(orientation: simd_double3x3, scale: Double)] = []
        for a in ia {
            for b in ib {
                let cosine = simd_dot(g[a], g[b]) / (gLength[a] * gLength[b])
                let gAngle = acos(min(1, max(-1, cosine)))
                for mode in modes {
                    var zi = 0.0, zk = 0.0
                    switch mode {
                    case .zero:
                        guard firstIsEarlier else { continue }
                    case .ewald:
                        guard firstIsEarlier else { continue }
                        zi = 0.5 * wavelength * gLength[a] * gLength[a]; zk = 0.5 * wavelength * gLength[b] * gLength[b]
                    case .layer(let sign):
                        // q_i in the zero layer, q_k at height z: g_a·g_b = q_i·q_k, so |(q_k, z)| = t / cos.
                        guard abs(cosine) > 0.05, t / cosine > 0 else { continue }
                        let z2 = (t / cosine) * (t / cosine) - nk * nk
                        guard z2 >= 0 else { continue }
                        zk = sign * z2.squareRoot()
                    case .length(let sign):
                        // q_i in the zero layer sets the scale; q_k's height from the length of g_b (g_a and g_b near
                        // perpendicular, where the angle cannot give it).
                        guard abs(cosine) <= 0.2 else { continue }
                        let z2 = (ni * gLength[b] / gLength[a]) * (ni * gLength[b] / gLength[a]) - nk * nk
                        guard z2 > 0 else { continue }
                        zk = sign * z2.squareRoot()
                    case .sameLayer(let sign), .oppositeLayer(let sign):
                        // Both spots in first-order layers, the same one (z_k = z_i) or opposite ones (z_k = −z_i): the
                        // common scale from the 3-D angle, z from the length of g_a.
                        let opposite: Double
                        if case .sameLayer = mode { opposite = 1 } else { opposite = -1 }
                        let ga2 = gLength[a] * gLength[a]
                        let s2 = (dot2 - opposite * ni * ni) / (cosine * gLength[a] * gLength[b] - opposite * ga2)
                        guard s2.isFinite, s2 > 0, abs(s2.squareRoot() - 1) <= window else { continue }
                        let z2 = s2 * ga2 - ni * ni
                        guard z2 > 1e-6 else { continue }
                        zi = sign * z2.squareRoot(); zk = opposite * zi
                    }
                    switch mode {
                    case .layer, .length:
                        guard abs(zk) <= s.higherLayerReach, abs(zk) > 1e-3 else { continue }
                    default: break
                    }
                    let li = (ni * ni + zi * zi).squareRoot(), lk = (nk * nk + zk * zk).squareRoot()
                    let measuredAngle = acos(min(1, max(-1, (dot2 + zi * zk) / (li * lk))))
                    let ratio = (gLength[a] / gLength[b]) / (li / lk)
                    guard abs(gAngle - measuredAngle) <= angleTolerance, abs(ratio - 1) <= 2 * s.tolerance else { continue }
                    var wc = simd_cross(g[a], g[b])
                    let nw = simd_length(wc)
                    guard nw >= 1e-9 else { continue }
                    wc /= nw
                    let u1 = g[a] / gLength[a]
                    let v1 = simd_normalize(SIMD3(qi.x, qi.y, zi))
                    let v2 = SIMD3(qk.x, qk.y, zk)
                    let wl = simd_normalize(simd_cross(v1, v2))
                    let crystal = simd_double3x3(columns: (u1, wc, simd_cross(u1, wc)))
                    let lab = simd_double3x3(columns: (v1, wl, simd_cross(v1, wl)))
                    // lab = R crystal maps the crystal triad onto the lab one; the orientation's columns are the lab axes
                    // in the crystal frame, Rᵀ = crystal labᵀ.
                    let m = crystal * lab.transpose
                    let scale = s.scaleWindow > 0 ? 0.5 * (ni / gLength[a] + nk / gLength[b]) : 1
                    out.append((m, scale))
                }
            }
        }
        return out
    }

    /// The orientation turned onto its nearest zone axis by the smallest crystal rotation; nil when that axis is farther
    /// than `snapAngle`.
    private func snapped(_ m: simd_double3x3) -> simd_double3x3? {
        guard let o = try? CrystalOrientation(matrix: m),
              let zone = o.nearestZoneAxis(cell: structure.cell, maxIndex: settings.maxIndex),
              zone.angle <= settings.snapAngle,
              let d = ZoneAxis(zone.uvw.x, zone.uvw.y, zone.uvw.z).cartesian(in: structure.cell) else { return nil }
        let n = m.columns.2
        let axis = simd_cross(n, d)
        let sine = simd_length(axis)
        guard sine >= 1e-15 else { return m }
        let r = simd_double3x3(simd_quatd(angle: atan2(sine, simd_dot(n, d)), axis: axis / sine))   // n → d
        return r * m
    }

    // MARK: - scoring

    private struct Scored {
        var matrix: simd_double3x3
        var scale: Double
        var agreement: Double
        var recall: Double
        var precision: Double
        var matched: Int
        var support: Double
        var rms: Double
        var assignments: [Int?]
    }

    private func score(_ m: simd_double3x3, _ q: [SIMD2<Double>], _ qn: [Double], scale scale0: Double, clicked: Bool)
        -> Scored? {
        let s = settings
        let reach = s.pattern.sigmaExcitationError * s.pattern.excitationMultiple
        let minusTwoSigmaSquared = -2 * s.pattern.sigmaExcitationError * s.pattern.sigmaExcitationError
        let mt = m.transpose
        var predicted: [SIMD2<Double>] = []
        var index: [Int] = []
        var weight: [Double] = []
        for j in g.indices {
            let p = mt * g[j]
            let sg = (2 * p.z - wavelength * gLength[j] * gLength[j]) / (2 - 2 * wavelength * p.z)
            guard sg.isFinite, abs(sg) <= reach else { continue }
            predicted.append(SIMD2(p.x, p.y))
            index.append(j)
            weight.append(intensity[j] * exp(sg * sg / minusTwoSigmaSquared))
        }
        guard predicted.count >= 2 else { return nil }
        let radius = qn.map { max(s.tolerance * $0, s.floor) }

        func nearest(_ rotation: Double, _ scale: Double) -> (j: [Int], d: [Double]) {
            let c = cos(rotation), sn = sin(rotation)
            var js = [Int](repeating: 0, count: q.count), ds = [Double](repeating: .infinity, count: q.count)
            for (j, p) in predicted.enumerated() {
                let pp = scale * SIMD2(c * p.x - sn * p.y, sn * p.x + c * p.y)
                for i in q.indices {
                    let d = simd_distance(q[i], pp)
                    if d < ds[i] { ds[i] = d; js[i] = j }
                }
            }
            return (js, ds)
        }

        var scale = scale0, rotation = 0.0
        let low = 1 - s.scaleWindow - 0.005, high = 1 + s.scaleWindow + 0.005
        for _ in 0..<2 {
            let (js, ds) = nearest(rotation, scale)
            var num = 0.0, den = 0.0, aa = 0.0, count = 0
            for i in q.indices where ds[i] <= radius[i] {
                let a = predicted[js[i]], b = q[i]
                num += a.x * b.y - a.y * b.x
                den += a.x * b.x + a.y * b.y
                aa += a.x * a.x + a.y * a.y
                count += 1
            }
            guard count >= 2 else { return nil }
            rotation = atan2(num, den)
            if s.scaleWindow > 0 { scale = min(max((num * num + den * den).squareRoot() / aa, low), high) }
        }
        let (js, ds) = nearest(rotation, scale)
        var matched = 0, sumSquares = 0.0
        var assignments: [Int?] = []
        var support = 0.0
        var weakestMatched = Double.infinity
        let maxWeight = weight.max() ?? 1
        for i in q.indices {
            if ds[i] <= radius[i] {
                matched += 1
                sumSquares += ds[i] * ds[i]
                assignments.append(index[js[i]])
                support += weight[js[i]] / maxWeight
                weakestMatched = min(weakestMatched, weight[js[i]])
            } else {
                assignments.append(nil)
            }
        }
        guard matched >= 2 else { return nil }
        let recall = Double(matched) / Double(q.count)
        var precision = 1.0
        if !clicked {
            let qMin = qn.min()! * (1 - s.tolerance), qMax = qn.max()! * (1 + s.tolerance)
            let c = cos(rotation), sn = sin(rotation)
            var inside = 0.0, hit = 0.0
            for (j, p) in predicted.enumerated() {
                let pp = scale * SIMD2(c * p.x - sn * p.y, sn * p.x + c * p.y)
                let r = simd_length(pp)
                // only predicted spots at least as strong as the weakest one matched: the detection limit of the pattern
                guard r >= qMin, r <= qMax, weight[j] >= weakestMatched else { continue }
                inside += weight[j]
                let tolerance = max(s.tolerance * r, s.floor)
                if q.contains(where: { simd_distance($0, pp) <= tolerance }) { hit += weight[j] }
            }
            guard inside > 0 else { return nil }
            precision = hit / inside
        }
        let b2 = s.recallWeight * s.recallWeight
        let agreement = precision + recall == 0 ? 0 : (1 + b2) * precision * recall / (b2 * precision + recall)
        // the in-plane turn by `rotation` about lab z: the columns become M Rzᵀ
        let c = cos(rotation), sn = sin(rotation)
        let rz = simd_double3x3(rows: [SIMD3(c, -sn, 0), SIMD3(sn, c, 0), SIMD3(0, 0, 1)])
        return Scored(matrix: m * rz.transpose, scale: scale, agreement: agreement, recall: recall, precision: precision,
                      matched: matched, support: support / Double(matched), rms: (sumSquares / Double(matched)).squareRoot(),
                      assignments: assignments)
    }

    /// The smallest rotation angle (radians) between two orientation matrices over the crystal's proper rotations.
    func misorientation(_ a: simd_double3x3, _ b: simd_double3x3) -> Double {
        var best = Double.pi
        let at = a.transpose
        for r in rotations {
            let t = at * r * b
            let trace = t.columns.0.x + t.columns.1.y + t.columns.2.z
            best = min(best, acos(min(1, max(-1, (trace - 1) / 2))))
        }
        return best
    }

    /// The smallest angle (radians) between two orientations of this crystal over its proper rotations and the half turn
    /// about the beam, which a centrosymmetric zero layer cannot tell apart.
    public func angle(between a: CrystalOrientation, and b: CrystalOrientation, includingHalfTurn: Bool = true) -> Double {
        let turned = misorientation(a.matrix * simd_double3x3(diagonal: SIMD3(-1, -1, 1)), b.matrix)
        return includingHalfTurn ? min(misorientation(a.matrix, b.matrix), turned) : misorientation(a.matrix, b.matrix)
    }
}
