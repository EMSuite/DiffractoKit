//
//  TwoBeamThickness.swift
//  Role: the specimen thickness from the intensity profile across the diffracted disc of a two-beam CBED pattern, with the
//        extinction distance, and the straight line through the fringe minima as the textbook cross-check.
//
//  Method (Kelly, Jostsons, Blake and Napier, phys. stat. sol. (a) 31, 771 (1975); Allen, Philos. Mag. A 43, 325 (1981)):
//    the two-beam intensity of the diffracted beam at deviation parameter s is
//        I(s) = sin²(π t s_eff) / (ξ s_eff)²,   s_eff² = (s − c)² + 1/ξ²,
//    with t the thickness along the beam, ξ the extinction distance and c the exact Bragg position on the profile; its
//    minima lie at s_n = √(n²/t² − 1/ξ²) for integers n > t/ξ, so (s_n/n)² = 1/t² − (1/ξ²)(1/n²): the fringe line.
//    The fit is the whole profile, as CrysTBox's default: scale·I(s) + background by least squares over (t, ξ, c), the
//    scale and the background solved linearly at every step. The extinction distance stays inside 0.7 to 1.15 of the
//    two-beam value from the phase's structure factor: fixed at that value the other beams bias the thickness by up to 11 %,
//    left free the fit lands on the wrong fringe order (issue #170, the cloud experiment
//    of 23 Nov 2026). Steps:
//      1. c₀ = the profile's centroid above its minimum.
//      2. A grid: thickness from the range's lower end by a ratio of 1.002, 25 extinction distances spaced geometrically
//         over the window; at each point the least-squares scale and background (a scale ≤ 0 is no fit); per thickness the
//         best extinction distance.
//      3. The grid's local minima in thickness, best first; up to four within four times the best's sum of squares are
//         refined by Levenberg-Marquardt on (t, ξ, c) with central-difference derivatives, ξ clamped to the window.
//      4. The best refined fit is the answer. A refined candidate more than 5 % away in thickness whose rms is within 5 % of
//         the best's is returned beside it: the profile cannot tell the two apart.
//      5. The minima the fit predicts inside the profile are counted on each side of c; fewer than two on a side is refused.
//      6. The fringe line: each predicted minimum is matched to the nearest measured minimum (a parabola through three
//         samples) within 0.35 of the gap to the next order, on both sides; the mean distance from c of each matched pair
//         and its order give the line by least squares.
//  No Python package does this; the pins are the closed form and a simulated disc (abTEM 1.0.10), tools/pins/two-beam-thickness
//  (fit.py there is this file's NumPy twin).
//

import Foundation
import DiffractoCore
import DiffractoScattering

/// An intensity profile across the diffracted disc of a two-beam CBED pattern, along g through the disc's centre.
public struct TwoBeamProfile: Sendable, Equatable {
    /// The deviation parameter s of each sample, Å⁻¹, strictly increasing.
    public let deviations: [Double]
    /// The intensity of each sample (any unit; a constant background is fitted).
    public let intensities: [Double]

    public init(deviations: [Double], intensities: [Double]) {
        self.deviations = deviations
        self.intensities = intensities
    }

    /// A profile measured on a pattern, with no calibration: the distance between the centres of the direct and the
    /// diffracted disc is twice the Bragg angle of the reflection (λ|g|), so a sample at `offsets[i]` from the diffracted
    /// disc's centre along g (positive away from the direct disc) has s = (x / D) λ |g|². `offsets` and `discSeparation`
    /// share any unit (pixels); `gLength` is |g| in Å⁻¹ and `wavelength` λ in Å.
    public init(offsets: [Double], discSeparation: Double, gLength: Double, wavelength: Double, intensities: [Double]) {
        let factor = wavelength * gLength * gLength / discSeparation
        self.deviations = offsets.map { $0 * factor }
        self.intensities = intensities
    }
}

/// The straight line through the fringe minima (Allen 1981): (s_n/n)² = 1/t² − (1/ξ²)(1/n²).
public struct TwoBeamFringeLine: Sendable, Equatable {
    /// The thickness from the intercept, Å.
    public let thickness: Double
    /// The extinction distance from the slope, Å; nil when the slope is not negative.
    public let extinctionDistance: Double?
    /// The orders n of the minima the line goes through, and their deviations s_n (Å⁻¹, measured from the Bragg position).
    public let orders: [Int]
    public let deviations: [Double]
}

/// One thickness that explains a profile.
public struct TwoBeamThicknessCandidate: Sendable, Equatable {
    /// Thickness along the beam, Å.
    public let thickness: Double
    /// The fitted extinction distance, Å.
    public let extinctionDistance: Double
    /// The exact Bragg position on the profile, Å⁻¹.
    public let centre: Double
    /// The rms of the fit's residual over the largest |intensity| of the profile.
    public let residual: Double
}

/// The thickness of one two-beam CBED profile.
public struct TwoBeamThicknessFit: Sendable, Equatable {
    /// Thickness along the beam, Å.
    public let thickness: Double
    /// The fitted extinction distance, Å (inside the window around `twoBeamExtinctionDistance`).
    public let extinctionDistance: Double
    /// The two-beam extinction distance the window is centred on, Å.
    public let twoBeamExtinctionDistance: Double
    /// The exact Bragg position on the profile (where s_eff is smallest), Å⁻¹.
    public let centre: Double
    /// The fitted profile is scale · I(s) + background.
    public let scale: Double
    public let background: Double
    /// The rms of the fit's residual over the largest |intensity| of the profile.
    public let residual: Double
    /// The minima the fit predicts inside the profile, below and above the Bragg position.
    public let minimaBelow: Int
    public let minimaAbove: Int
    /// The fringe line through the measured minima nearest the predicted ones; nil when fewer than two orders have a
    /// measured minimum on both sides.
    public let fringeLine: TwoBeamFringeLine?
    /// A second thickness more than 5 % away that fits within 5 % of this one's residual: the profile cannot tell them
    /// apart, and both are the answer.
    public let secondCandidate: TwoBeamThicknessCandidate?

    /// Two thicknesses fit equally.
    public var isAmbiguous: Bool { secondCandidate != nil }
    /// The fitted extinction distance sits on an edge of the window (within 1e-6 of it): the two-beam picture may not hold.
    public var extinctionDistanceAtWindowEdge: Bool {
        let lo = TwoBeamThickness.extinctionWindow.lowerBound * twoBeamExtinctionDistance
        let hi = TwoBeamThickness.extinctionWindow.upperBound * twoBeamExtinctionDistance
        return extinctionDistance <= lo * (1 + 1e-6) || extinctionDistance >= hi * (1 - 1e-6)
    }
    /// The best fit as a candidate.
    public var candidate: TwoBeamThicknessCandidate {
        TwoBeamThicknessCandidate(thickness: thickness, extinctionDistance: extinctionDistance, centre: centre, residual: residual)
    }
}

/// Why no thickness was returned.
public enum TwoBeamThicknessError: Error, Equatable, Sendable, CustomStringConvertible {
    /// The deviations and the intensities differ in count.
    case mismatchedLengths(deviations: Int, intensities: Int)
    /// Fewer samples than `TwoBeamThickness.minimumSamples`.
    case tooFewSamples(Int)
    /// A deviation or an intensity at this index is not a finite number.
    case nonFiniteSample(Int)
    /// The deviation at this index is not larger than the one before.
    case deviationsNotIncreasing(Int)
    /// Every intensity is the same: there are no fringes.
    case flatProfile
    /// The extinction distance is not a finite positive number of Å.
    case invalidExtinctionDistance(Double)
    /// The thickness range is not finite and positive with its lower end below its upper.
    case invalidThicknessRange
    /// The voltage is not a finite positive number of volts.
    case invalidVoltage(Double)
    /// The phase has no structure factor for this reflection (an element without a form factor), or it is zero.
    case noStructureFactor(MillerIndex)
    /// The fit predicts fewer than two minima on a side of the Bragg position inside the profile: too thin a specimen or
    /// too short a profile for the fringes to fix the thickness.
    case tooFewMinima(below: Int, above: Int)
    /// The fringe line needs at least two minima with positive orders and matching counts.
    case invalidFringeMinima
    /// The fit gave no finite answer.
    case noFit

    public var description: String {
        switch self {
        case .mismatchedLengths(let d, let i):
            return "The profile has \(d) deviations and \(i) intensities; they must match."
        case .tooFewSamples(let n):
            return "A disc profile needs at least \(TwoBeamThickness.minimumSamples) samples; \(n) were given."
        case .nonFiniteSample(let i): return "Sample \(i) of the profile is not a finite number."
        case .deviationsNotIncreasing(let i): return "The deviations must increase along the profile; sample \(i) does not."
        case .flatProfile: return "The profile is flat: there are no fringes to measure."
        case .invalidExtinctionDistance(let x): return "The extinction distance must be a finite number of Å above 0; \(x) is not."
        case .invalidThicknessRange: return "The thickness range must be finite and positive, its lower end below its upper."
        case .invalidVoltage(let v): return "The voltage must be a finite number of volts above 0; \(v) is not."
        case .noStructureFactor(let hkl):
            return "The phase gives no structure factor for \(hkl), so its extinction distance is unknown."
        case .tooFewMinima(let below, let above):
            return "The profile shows \(below) fringe minima below and \(above) above the Bragg position; the thickness needs "
                + "at least two on each side. Take a longer profile or a thicker area."
        case .invalidFringeMinima: return "The fringe line needs at least two minima, each with a positive order."
        case .noFit: return "The two-beam fit gave no finite answer for this profile."
        }
    }
}

extension TwoBeamThicknessError: LocalizedError {
    public var errorDescription: String? { description }
}

/// Thickness from a two-beam CBED disc profile.
public enum TwoBeamThickness {

    /// The window of the fitted extinction distance, as fractions of the two-beam value.
    public static let extinctionWindow: ClosedRange<Double> = 0.7...1.15
    /// The default thickness range searched, Å (10 to 500 nm).
    public static let defaultThicknessRange: ClosedRange<Double> = 100...5000
    /// The fewest samples a profile may have.
    public static let minimumSamples = 16

    static let extinctionSteps = 25
    static let thicknessRatio = 1.002
    static let secondMargin = 0.05
    static let secondApart = 0.05
    static let refinedCandidates = 4

    // MARK: The closed form

    /// The two-beam intensity of the diffracted beam, I(s) = sin²(π t s_eff) / (ξ s_eff)², s_eff² = s² + 1/ξ², for a
    /// deviation s (Å⁻¹), a thickness t and an extinction distance ξ (Å). The direct beam carries 1 − I.
    public static func intensity(deviation s: Double, thickness t: Double, extinctionDistance xi: Double) -> Double {
        let se = (s * s + 1 / (xi * xi)).squareRoot()
        let sine = sin(Double.pi * t * se)
        let den = xi * se
        return sine * sine / (den * den)
    }

    /// The deviation of the minimum of order n, s_n = √(n²/t² − 1/ξ²), Å⁻¹; nil when n ≤ t/ξ (no minimum of that order).
    public static func minimumDeviation(order n: Int, thickness t: Double, extinctionDistance xi: Double) -> Double? {
        let nd = Double(n)
        let v = nd * nd / (t * t) - 1 / (xi * xi)
        guard v > 0 else { return nil }
        return v.squareRoot()
    }

    /// The fringe line by least squares through minima of the given orders at the given deviations (Å⁻¹, measured from the
    /// Bragg position): thickness from the intercept, extinction distance from the slope.
    public static func fringeLine(orders: [Int], deviations: [Double]) throws(TwoBeamThicknessError) -> TwoBeamFringeLine {
        guard orders.count == deviations.count, orders.count >= 2, orders.allSatisfy({ $0 > 0 }),
              deviations.allSatisfy({ $0.isFinite }) else { throw .invalidFringeMinima }
        let count = Double(orders.count)
        var xs: [Double] = []
        var ys: [Double] = []
        for k in 0..<orders.count {
            let n = Double(orders[k])
            let q = deviations[k] / n
            xs.append(1 / (n * n))
            ys.append(q * q)
        }
        let xm = xs.reduce(0, +) / count
        let ym = ys.reduce(0, +) / count
        var sxx = 0.0
        var sxy = 0.0
        for k in 0..<xs.count {
            sxx += (xs[k] - xm) * (xs[k] - xm)
            sxy += (xs[k] - xm) * (ys[k] - ym)
        }
        guard sxx > 0 else { throw .invalidFringeMinima }
        let slope = sxy / sxx
        let intercept = ym - slope * xm
        guard intercept > 0, intercept.isFinite else { throw .invalidFringeMinima }
        let xi: Double? = slope < 0 ? 1 / (-slope).squareRoot() : nil
        return TwoBeamFringeLine(thickness: 1 / intercept.squareRoot(), extinctionDistance: xi, orders: orders,
                                 deviations: deviations)
    }

    // MARK: The extinction distance of a phase

    /// The two-beam extinction distance ξ_g = π V cos θ_B / (λ γ |F_g|) of a reflection (Å), with the kit's electron
    /// structure factor F (Å⁻², which carries the 1/V), the relativistic wavelength λ and mass factor γ at `volts`.
    public static func extinctionDistance(of hkl: MillerIndex, in structure: CrystalStructure,
                                          volts: Double) throws(TwoBeamThicknessError) -> Double {
        guard volts.isFinite, volts > 0 else { throw .invalidVoltage(volts) }
        guard let f = StructureFactors.compute(hkl, in: structure), f.magnitude > 0, f.magnitude.isFinite else {
            throw .noStructureFactor(hkl)
        }
        let lambda = ElectronWavelength.angstroms(volts: volts)
        let gamma = 1 + volts / ElectronWavelength.electronRestEnergyEV
        let sinBragg = lambda * f.gLength / 2
        guard sinBragg < 1 else { throw .noStructureFactor(hkl) }
        let cosBragg = (1 - sinBragg * sinBragg).squareRoot()
        let xi = Double.pi * cosBragg / (lambda * gamma * f.magnitude)
        guard xi.isFinite, xi > 0 else { throw .noStructureFactor(hkl) }
        return xi
    }

    // MARK: The fit

    /// The thickness of `profile` for a reflection of a phase at a voltage: the two-beam extinction distance from the
    /// phase's structure factor, then `fit(_:extinctionDistance:thicknessRange:)`.
    public static func fit(_ profile: TwoBeamProfile, reflection hkl: MillerIndex, in structure: CrystalStructure,
                           volts: Double,
                           thicknessRange: ClosedRange<Double> = defaultThicknessRange) throws(TwoBeamThicknessError)
        -> TwoBeamThicknessFit {
        let xi = try extinctionDistance(of: hkl, in: structure, volts: volts)
        return try fit(profile, extinctionDistance: xi, thicknessRange: thicknessRange)
    }

    /// The thickness of `profile`, its extinction distance fitted inside `extinctionWindow` times `extinctionDistance`
    /// (the two-beam value, Å), searched over `thicknessRange` (Å). Refused when the fit predicts fewer than two minima on
    /// a side of the Bragg position.
    public static func fit(_ profile: TwoBeamProfile, extinctionDistance xi0: Double,
                           thicknessRange: ClosedRange<Double> = defaultThicknessRange) throws(TwoBeamThicknessError)
        -> TwoBeamThicknessFit {
        let s = profile.deviations
        let y = profile.intensities
        try validate(s, y)
        guard xi0.isFinite, xi0 > 0 else { throw .invalidExtinctionDistance(xi0) }
        guard thicknessRange.lowerBound.isFinite, thicknessRange.upperBound.isFinite, thicknessRange.lowerBound > 0,
              thicknessRange.upperBound > thicknessRange.lowerBound else { throw .invalidThicknessRange }

        let model = Model(s: s, y: y)
        let c0 = model.centroid()
        let grid = Grid(model: model, xi0: xi0, range: thicknessRange, centre: c0)
        let minima = grid.localMinima()
        guard let first = minima.first else { throw .noFit }

        let lower = [thicknessRange.lowerBound * 0.5, xi0 * extinctionWindow.lowerBound, s[0]]
        let upper = [thicknessRange.upperBound * 2, xi0 * extinctionWindow.upperBound, s[s.count - 1]]
        var refined: [(ssr: Double, p: [Double])] = []
        for i in minima.prefix(refinedCandidates) {
            if i != first, grid.ssr[i] > 4 * grid.ssr[first] { break }
            let start = [grid.thickness[i], grid.bestXi[i], c0]
            let r = model.levenbergMarquardt(start: start, lower: lower, upper: upper)
            if r.ssr.isFinite { refined.append(r) }
        }
        refined.sort { $0.ssr < $1.ssr }
        guard let best = refined.first else { throw .noFit }
        let t = best.p[0]
        let xi = best.p[1]
        let c = best.p[2]
        guard t.isFinite, xi.isFinite, c.isFinite else { throw .noFit }

        let peak = y.map { abs($0) }.max() ?? 1
        let norm = (Double(s.count)).squareRoot() * peak
        var second: TwoBeamThicknessCandidate?
        for r in refined.dropFirst() {
            let apart = abs(r.p[0] / t - 1) > secondApart
            let equal = r.ssr.squareRoot() <= (1 + secondMargin) * best.ssr.squareRoot()
            if apart, equal {
                second = TwoBeamThicknessCandidate(thickness: r.p[0], extinctionDistance: r.p[1], centre: r.p[2],
                                                   residual: r.ssr.squareRoot() / norm)
                break
            }
        }

        let counts = predictedMinima(t: t, xi: xi, c: c, first: s[0], last: s[s.count - 1])
        guard counts.below >= 2, counts.above >= 2 else { throw .tooFewMinima(below: counts.below, above: counts.above) }

        let linear = model.linear(model.shape(t: t, xi: xi, c: c))
        let line = fringeLine(model: model, t: t, xi: xi, c: c)
        return TwoBeamThicknessFit(thickness: t, extinctionDistance: xi, twoBeamExtinctionDistance: xi0, centre: c,
                                   scale: linear.scale, background: linear.background, residual: best.ssr.squareRoot() / norm,
                                   minimaBelow: counts.below, minimaAbove: counts.above, fringeLine: line,
                                   secondCandidate: second)
    }

    static func validate(_ s: [Double], _ y: [Double]) throws(TwoBeamThicknessError) {
        guard s.count == y.count else { throw .mismatchedLengths(deviations: s.count, intensities: y.count) }
        guard s.count >= minimumSamples else { throw .tooFewSamples(s.count) }
        for i in 0..<s.count where !(s[i].isFinite && y[i].isFinite) { throw .nonFiniteSample(i) }
        for i in 1..<s.count where !(s[i] > s[i - 1]) { throw .deviationsNotIncreasing(i) }
        let lo = y.min() ?? 0
        let hi = y.max() ?? 0
        guard hi > lo else { throw .flatProfile }
    }

    /// The predicted minima inside [first, last], below and above c.
    static func predictedMinima(t: Double, xi: Double, c: Double, first: Double,
                                last: Double) -> (below: Int, above: Int) {
        var below = 0
        var above = 0
        var n = Int((t / xi).rounded(.down)) + 1
        while let sn = minimumDeviation(order: n, thickness: t, extinctionDistance: xi) {
            let upIn = c + sn <= last
            let downIn = c - sn >= first
            if !upIn, !downIn { break }
            if upIn { above += 1 }
            if downIn { below += 1 }
            n += 1
        }
        return (below, above)
    }

    /// Step 6: the measured minima nearest the predicted ones, both sides, then the line.
    static func fringeLine(model: Model, t: Double, xi: Double, c: Double) -> TwoBeamFringeLine? {
        let measured = model.measuredMinima()
        guard !measured.isEmpty else { return nil }
        let first = model.s[0]
        let last = model.s[model.s.count - 1]
        var orders: [Int] = []
        var deviations: [Double] = []
        var n = Int((t / xi).rounded(.down)) + 1
        while let sp = minimumDeviation(order: n, thickness: t, extinctionDistance: xi),
              let sn = minimumDeviation(order: n + 1, thickness: t, extinctionDistance: xi) {
            if c + sp > last, c - sp < first { break }
            let tolerance = 0.35 * (sn - sp)
            var pair: [Double] = []
            for sign in [-1.0, 1.0] {
                let target = c + sign * sp
                var bestIndex = 0
                for k in 1..<measured.count where abs(measured[k] - target) < abs(measured[bestIndex] - target) {
                    bestIndex = k
                }
                if abs(measured[bestIndex] - target) < tolerance { pair.append(sign * (measured[bestIndex] - c)) }
            }
            if pair.count == 2 {
                orders.append(n)
                deviations.append(0.5 * (pair[0] + pair[1]))
            }
            n += 1
        }
        guard orders.count >= 2 else { return nil }
        return try? fringeLine(orders: orders, deviations: deviations)
    }

    // MARK: The model

    /// The profile and the least-squares pieces of the fit.
    struct Model {
        let s: [Double]
        let y: [Double]
        let yMean: Double

        init(s: [Double], y: [Double]) {
            self.s = s
            self.y = y
            self.yMean = y.reduce(0, +) / Double(y.count)
        }

        /// Σ s (y − min y) / Σ (y − min y).
        func centroid() -> Double {
            let lo = y.min() ?? 0
            var num = 0.0
            var den = 0.0
            for i in 0..<s.count {
                let w = y[i] - lo
                num += s[i] * w
                den += w
            }
            return num / den
        }

        func shape(t: Double, xi: Double, c: Double) -> [Double] {
            var p = [Double](repeating: 0, count: s.count)
            let inv = 1 / (xi * xi)
            for i in 0..<s.count {
                let d = s[i] - c
                let se = (d * d + inv).squareRoot()
                let sine = sin(Double.pi * t * se)
                let den = xi * se
                p[i] = sine * sine / (den * den)
            }
            return p
        }

        /// The least-squares scale and background of the shape p (centred); a scale ≤ 0 is no fit: the constant mean.
        func linear(_ p: [Double]) -> (scale: Double, background: Double) {
            let pm = p.reduce(0, +) / Double(p.count)
            var sxx = 0.0
            var sxy = 0.0
            for i in 0..<p.count {
                let dp = p[i] - pm
                sxx += dp * dp
                sxy += dp * (y[i] - yMean)
            }
            let a = sxx > 0 ? sxy / sxx : 0
            guard a > 0 else { return (0, yMean) }
            return (a, yMean - a * pm)
        }

        func residuals(_ q: [Double]) -> [Double] {
            let p = shape(t: q[0], xi: q[1], c: q[2])
            let fit = linear(p)
            var r = [Double](repeating: 0, count: p.count)
            for i in 0..<p.count { r[i] = fit.scale * p[i] + fit.background - y[i] }
            return r
        }

        func sumOfSquares(_ r: [Double]) -> Double {
            var total = 0.0
            for v in r { total += v * v }
            return total
        }

        /// Local minima of the samples, a parabola through each and its neighbours (uniform spacing assumed locally).
        func measuredMinima() -> [Double] {
            var out: [Double] = []
            guard y.count >= 3 else { return out }
            for i in 1..<(y.count - 1) where y[i] < y[i - 1] && y[i] <= y[i + 1] {
                let den = y[i - 1] - 2 * y[i] + y[i + 1]
                let ds = 0.5 * (s[i + 1] - s[i - 1])
                let shift = den > 0 ? 0.5 * (y[i - 1] - y[i + 1]) / den * ds : 0
                out.append(s[i] + shift)
            }
            return out
        }

        /// Levenberg-Marquardt on (t, ξ, c), scaled by the start's t and ξ and the profile's span, each step clamped to
        /// [lower, upper]; central differences of relative step 1e-7. Stops when a step moves every parameter by at most
        /// 1e-12 of its scale, when no step lowers the sum of squares, or after 200 iterations.
        func levenbergMarquardt(start: [Double], lower: [Double], upper: [Double]) -> (ssr: Double, p: [Double]) {
            var x = start
            let scale = [start[0], start[1], s[s.count - 1] - s[0]]
            var r = residuals(x)
            var f = sumOfSquares(r)
            var lambda = 1e-3
            for _ in 0..<200 {
                var jac = [[Double]](repeating: [0, 0, 0], count: s.count)
                for k in 0..<3 {
                    let h = 1e-7 * scale[k]
                    var xp = x
                    var xm = x
                    xp[k] += h
                    xm[k] -= h
                    let rp = residuals(xp)
                    let rm = residuals(xm)
                    for i in 0..<s.count { jac[i][k] = (rp[i] - rm[i]) / (2 * h) * scale[k] }
                }
                var jtj = [[Double]](repeating: [0, 0, 0], count: 3)
                var jtr = [0.0, 0.0, 0.0]
                for i in 0..<s.count {
                    for a in 0..<3 {
                        jtr[a] += jac[i][a] * r[i]
                        for b in 0..<3 { jtj[a][b] += jac[i][a] * jac[i][b] }
                    }
                }
                var improved = false
                var done = false
                while lambda < 1e12 {
                    var m = jtj
                    for a in 0..<3 { m[a][a] += lambda * jtj[a][a] }
                    guard let step = solve3(m, [-jtr[0], -jtr[1], -jtr[2]]) else {
                        lambda *= 10
                        continue
                    }
                    var xn = x
                    for k in 0..<3 { xn[k] = min(max(x[k] + step[k] * scale[k], lower[k]), upper[k]) }
                    let rn = residuals(xn)
                    let fn = sumOfSquares(rn)
                    if fn < f {
                        done = true
                        for k in 0..<3 where abs(xn[k] - x[k]) > 1e-12 * abs(scale[k]) { done = false }
                        x = xn
                        r = rn
                        f = fn
                        lambda = max(lambda / 10, 1e-12)
                        improved = true
                        break
                    }
                    lambda *= 10
                }
                if !improved || done { break }
            }
            return (f, x)
        }
    }

    /// Step 2: per thickness of the grid, the best sum of squares over the extinction distances, at c₀.
    struct Grid {
        let thickness: [Double]
        let ssr: [Double]
        let bestXi: [Double]

        init(model: Model, xi0: Double, range: ClosedRange<Double>, centre c0: Double) {
            let lo = TwoBeamThickness.extinctionWindow.lowerBound
            let hi = TwoBeamThickness.extinctionWindow.upperBound
            let count = TwoBeamThickness.extinctionSteps
            let ratio = TwoBeamThickness.thicknessRatio
            var xis: [Double] = []
            for k in 0..<count {
                let fraction = Double(k) / Double(count - 1)
                xis.append(xi0 * lo * pow(hi / lo, fraction))
            }
            let steps = Int((log(range.upperBound / range.lowerBound) / log(ratio)).rounded(.down)) + 1
            var ts: [Double] = []
            for i in 0..<steps { ts.append(range.lowerBound * pow(ratio, Double(i))) }
            let n = model.s.count
            // Per extinction distance, s_eff and the denominator (ξ s_eff)² of every sample.
            var seTable: [[Double]] = []
            var denTable: [[Double]] = []
            for xi in xis {
                var se = [Double](repeating: 0, count: n)
                var den = [Double](repeating: 0, count: n)
                for i in 0..<n {
                    let d = model.s[i] - c0
                    se[i] = (d * d + 1 / (xi * xi)).squareRoot()
                    let q = xi * se[i]
                    den[i] = q * q
                }
                seTable.append(se)
                denTable.append(den)
            }
            var best = [Double](repeating: .infinity, count: steps)
            var bestXi = [Double](repeating: xi0, count: steps)
            var p = [Double](repeating: 0, count: n)
            for j in 0..<steps {
                let t = ts[j]
                for k in 0..<xis.count {
                    let se = seTable[k]
                    let den = denTable[k]
                    for i in 0..<n {
                        let sine = sin(Double.pi * t * se[i])
                        p[i] = sine * sine / den[i]
                    }
                    let fit = model.linear(p)
                    var total = 0.0
                    for i in 0..<n {
                        let r = fit.scale * p[i] + fit.background - model.y[i]
                        total += r * r
                    }
                    if total < best[j] {
                        best[j] = total
                        bestXi[j] = xis[k]
                    }
                }
            }
            self.thickness = ts
            self.ssr = best
            self.bestXi = bestXi
        }

        /// Interior local minima of the sum of squares in thickness, best first.
        func localMinima() -> [Int] {
            guard ssr.count >= 3 else { return [] }
            var out: [Int] = []
            for i in 1..<(ssr.count - 1) where ssr[i] < ssr[i - 1] && ssr[i] <= ssr[i + 1] { out.append(i) }
            return out.sorted { ssr[$0] < ssr[$1] }
        }
    }
}

/// The 3 × 3 system m x = b by Gaussian elimination with partial pivoting; nil when singular.
private func solve3(_ m: [[Double]], _ b: [Double]) -> [Double]? {
    var a = m
    var v = b
    for col in 0..<3 {
        var pivot = col
        for row in (col + 1)..<3 where abs(a[row][col]) > abs(a[pivot][col]) { pivot = row }
        guard abs(a[pivot][col]) > 0, a[pivot][col].isFinite else { return nil }
        if pivot != col {
            a.swapAt(pivot, col)
            v.swapAt(pivot, col)
        }
        for row in (col + 1)..<3 {
            let factor = a[row][col] / a[col][col]
            for k in col..<3 { a[row][k] -= factor * a[col][k] }
            v[row] -= factor * v[col]
        }
    }
    var x = [0.0, 0.0, 0.0]
    for row in stride(from: 2, through: 0, by: -1) {
        var total = v[row]
        for k in (row + 1)..<3 { total -= a[row][k] * x[k] }
        x[row] = total / a[row][row]
    }
    return x
}
