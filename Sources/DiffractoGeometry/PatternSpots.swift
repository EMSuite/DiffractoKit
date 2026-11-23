//
//  PatternSpots.swift
//  Role: every spot of one diffraction pattern (DiffractoKit S4): py4DSTEM's get_maxima_2D (preprocess/utils.py at f050d207:
//        Gaussian smoothing, the 8-neighbour maxima, the edge boundary, brightest first with NumPy's tie order, the filter
//        cascade of filter_2D_maxima, the parabolic sub-pixel step and the bilinear intensity), after the subtraction of the
//        pattern's radial median, then each spot refined with the held `PatternSpot.refine`.
//
//  `maxima(in:settings:)` is get_maxima_2D on the pixels as they are and is the pinned step (P41). `find(in:centre:settings:)`
//  is the route an application calls. The Gaussian is SciPy's `ndimage.gaussian_filter` (truncate 4, mode "reflect", axis 0
//  then axis 1, the symmetric sum of `correlate1d`), written here rather than taken from ImagoFilters so that DiffractoGeometry
//  keeps its single ImagoCore dependency and the kit's ImagoKit floor.
//

import Foundation
import simd

/// One spot of a pattern. Rows and columns are pixel indices of `DiffractionPattern` (row = y, column = x).
public struct FoundSpot: Sendable, Equatable {
    /// The spot's centre: the refined one when `isRefined`, else the maximum's.
    public let row: Double
    public let column: Double
    /// get_maxima_2D's position (the parabolic step on the smoothed image).
    public let maximumRow: Double
    public let maximumColumn: Double
    /// get_maxima_2D's intensity: the smoothed (and, in `find`, background-free) image interpolated at the maximum.
    public let intensity: Double
    /// Whether `PatternSpot.refine` gave the centre (false when refinement is off or it refused the window; then the centre
    /// is the maximum's).
    public let isRefined: Bool

    public init(row: Double, column: Double, maximumRow: Double, maximumColumn: Double, intensity: Double, isRefined: Bool) {
        self.row = row; self.column = column
        self.maximumRow = maximumRow; self.maximumColumn = maximumColumn
        self.intensity = intensity; self.isRefined = isRefined
    }
}

/// Finds the spots of one diffraction pattern.
public enum PatternSpots {

    /// The settings. The maxima's are py4DSTEM's `find_Bragg_disks` values for its call of get_maxima_2D (sigma_cc 2,
    /// minRelativeIntensity 0.005 of the brightest, edgeBoundary 20, maxNumPeaks 70) except the spacing: py4DSTEM's 60 px is a
    /// 4D-STEM detector's, and a selected-area pattern's neighbouring spots lie closer, so 10 px.
    public struct Settings: Sendable, Equatable {
        /// σ of the Gaussian smoothing, pixels (get_maxima_2D `sigma`; 0 for none). Default 2.
        public var sigma: Double
        /// Maxima below this intensity are dropped (`minAbsoluteIntensity`). Default 0.
        public var minimumAbsoluteIntensity: Double
        /// Maxima below this fraction of the `relativeToPeak`-th brightest are dropped (`minRelativeIntensity`). Default 0.005.
        public var minimumRelativeIntensity: Double
        /// Which maximum the relative cut refers to, 0 the brightest (`relativeToPeak`). Default 0.
        public var relativeToPeak: Int
        /// Of two maxima closer than this (pixels) the dimmer is dropped (`minSpacing`). Default 10.
        public var minimumSpacing: Double
        /// Maxima within this many pixels of the edge are dropped; below 1 counts as 1 (`edgeBoundary`). Default 20.
        public var edgeBoundary: Int
        /// At most this many maxima, brightest first (`maxNumPeaks`). Default 70.
        public var maximumCount: Int
        /// `find` subtracts the median of each one-pixel ring about the centre and clips at 0 before the maxima. Default true.
        /// The innermost rings are the direct beam itself, so the subtraction also takes the beam's core away.
        public var subtractsRadialBackground: Bool
        /// The half width of `PatternSpot.refine`'s window, on the background-free image; 0 for no refinement. Default 3.
        public var refinementHalfWindow: Int
        /// `find` drops spots whose maximum lies closer than this to the centre (pixels): 0 keeps the direct beam. Default 0.
        public var beamExclusionRadius: Double

        public init(sigma: Double = 2, minimumAbsoluteIntensity: Double = 0, minimumRelativeIntensity: Double = 0.005,
                    relativeToPeak: Int = 0, minimumSpacing: Double = 10, edgeBoundary: Int = 20, maximumCount: Int = 70,
                    subtractsRadialBackground: Bool = true, refinementHalfWindow: Int = 3, beamExclusionRadius: Double = 0) {
            self.sigma = sigma
            self.minimumAbsoluteIntensity = minimumAbsoluteIntensity
            self.minimumRelativeIntensity = minimumRelativeIntensity
            self.relativeToPeak = relativeToPeak
            self.minimumSpacing = minimumSpacing
            self.edgeBoundary = edgeBoundary
            self.maximumCount = maximumCount
            self.subtractsRadialBackground = subtractsRadialBackground
            self.refinementHalfWindow = refinementHalfWindow
            self.beamExclusionRadius = beamExclusionRadius
        }
    }

    /// Why a pattern or the settings were refused.
    public enum Failure: Error, Sendable, Equatable, CustomStringConvertible {
        case invalidSettings(String)
        /// No pixels, or a pixel count that is not width × height.
        case emptyPattern
        /// A pixel is NaN or infinite.
        case nonFinitePixel
        /// The centre is not finite or lies outside the pattern.
        case centreOutside

        public var description: String {
            switch self {
            case .invalidSettings(let s): "invalid spot-finding settings: \(s)"
            case .emptyPattern: "the pattern has no pixels or its size does not match"
            case .nonFinitePixel: "the pattern has a pixel that is not a finite number"
            case .centreOutside: "the centre is not inside the pattern"
            }
        }
    }

    /// py4DSTEM's get_maxima_2D(subpixel="poly") on the pattern's pixels as they are, with `settings`' sigma and filters
    /// (the background, refinement and beam settings are not used). Brightest first; equal intensities in NumPy's order of
    /// its structured sort (larger row, then larger column, first). `row`/`column` are the maximum's, `isRefined` false.
    public static func maxima(in pattern: DiffractionPattern, settings: Settings = Settings()) throws(Failure) -> [FoundSpot] {
        try validate(settings)
        let values = try pixels(of: pattern)
        return getMaxima(values, rows: pattern.height, columns: pattern.width, settings).map {
            FoundSpot(row: $0.row, column: $0.column, maximumRow: $0.row, maximumColumn: $0.column, intensity: $0.intensity,
                      isRefined: false)
        }
    }

    /// Every spot of a pattern: the radial median about `centre` subtracted (clipped at 0), get_maxima_2D's maxima, spots
    /// nearer the centre than `beamExclusionRadius` dropped, then each refined by `PatternSpot.refine` on the background-free
    /// image. Brightest first.
    /// - Parameter centre: the pattern centre (the unscattered beam), as (row, column) pixels.
    public static func find(in pattern: DiffractionPattern, centre: SIMD2<Double>,
                            settings: Settings = Settings()) throws(Failure) -> [FoundSpot] {
        try validate(settings)
        var values = try pixels(of: pattern)
        let rows = pattern.height, columns = pattern.width
        guard centre.x.isFinite, centre.y.isFinite, centre.x >= 0, centre.x <= Double(rows - 1), centre.y >= 0,
              centre.y <= Double(columns - 1) else { throw .centreOutside }
        if settings.subtractsRadialBackground { values = subtractingRadialMedian(values, rows: rows, columns: columns, centre: centre) }
        var found = getMaxima(values, rows: rows, columns: columns, settings)
        if settings.beamExclusionRadius > 0 {
            found.removeAll { hypot($0.row - centre.x, $0.column - centre.y) < settings.beamExclusionRadius }
        }
        guard settings.refinementHalfWindow > 0 else {
            return found.map { FoundSpot(row: $0.row, column: $0.column, maximumRow: $0.row, maximumColumn: $0.column,
                                         intensity: $0.intensity, isRefined: false) }
        }
        let clean = DiffractionPattern(width: columns, height: rows, pixels: values.map { Float($0) })
        return found.map { m in
            if let r = PatternSpot.refine(pattern: clean, row: m.row, column: m.column, half: settings.refinementHalfWindow) {
                return FoundSpot(row: r.row, column: r.column, maximumRow: m.row, maximumColumn: m.column, intensity: m.intensity,
                                 isRefined: true)
            }
            return FoundSpot(row: m.row, column: m.column, maximumRow: m.row, maximumColumn: m.column, intensity: m.intensity,
                             isRefined: false)
        }
    }

    // MARK: get_maxima_2D

    typealias Maximum = (row: Double, column: Double, intensity: Double)

    static func getMaxima(_ input: [Double], rows: Int, columns: Int, _ s: Settings) -> [Maximum] {
        let ar = s.sigma > 0 ? gaussianFilter(input, rows: rows, columns: columns, sigma: s.sigma) : input
        let eb = max(1, s.edgeBoundary)
        var found: [Maximum] = []
        if rows - eb > eb, columns - eb > eb {
            let w = columns
            for r in eb..<(rows - eb) {
                for c in eb..<(columns - eb) {
                    let i = r * w + c, v = ar[i]
                    // NumPy's rolls: >= toward the following row or column, > toward the preceding one.
                    if v >= ar[i + w], v > ar[i - w], v >= ar[i + 1], v > ar[i - 1],
                       v >= ar[i + w + 1], v > ar[i + w - 1], v >= ar[i - w + 1], v > ar[i - w - 1] {
                        found.append((Double(r), Double(c), v))
                    }
                }
            }
        }
        // np.sort(order="intensity")[::-1]: ascending by intensity, then row, then column; reversed.
        found.sort { a, b in
            if a.intensity != b.intensity { return a.intensity > b.intensity }
            if a.row != b.row { return a.row > b.row }
            return a.column > b.column
        }
        guard !found.isEmpty else { return [] }
        // filter_2D_maxima
        if s.minimumAbsoluteIntensity > 0 { found.removeAll { $0.intensity < s.minimumAbsoluteIntensity } }
        if s.minimumRelativeIntensity > 0, found.count > s.relativeToPeak {
            let reference = found[s.relativeToPeak].intensity
            found.removeAll { $0.intensity / reference < s.minimumRelativeIntensity }
        }
        if s.minimumSpacing > 0 {
            let limit = s.minimumSpacing * s.minimumSpacing
            var deleted = [Bool](repeating: false, count: found.count)
            for i in found.indices where !deleted[i] {
                for j in (i + 1)..<found.count where !deleted[j] {
                    let dr = found[j].row - found[i].row, dc = found[j].column - found[i].column
                    if dr * dr + dc * dc < limit { deleted[j] = true }
                }
            }
            found = found.indices.filter { !deleted[$0] }.map { found[$0] }
        }
        if found.count > s.maximumCount { found = Array(found.prefix(s.maximumCount)) }
        // The parabola through each axis' three pixels; its denominator is positive at a strict maximum.
        let w = columns
        for k in found.indices {
            let r = Int(found[k].row), c = Int(found[k].column), i = r * w + c
            let i0 = ar[i], up = ar[i + w], down = ar[i - w], right = ar[i + 1], left = ar[i - 1]
            let dr = (up - down) / (4 * i0 - 2 * up - 2 * down)
            let dc = (right - left) / (4 * i0 - 2 * right - 2 * left)
            found[k].row += dr.isFinite ? dr : 0
            found[k].column += dc.isFinite ? dc : 0
            found[k].intensity = bilinear(ar, columns: columns, row: found[k].row, column: found[k].column)
        }
        return found
    }

    /// py4DSTEM's linear_interpolation_2D: the floor and ceiling neighbours.
    static func bilinear(_ a: [Double], columns: Int, row: Double, column: Double) -> Double {
        let r0 = Int(row.rounded(.down)), r1 = Int(row.rounded(.up))
        let c0 = Int(column.rounded(.down)), c1 = Int(column.rounded(.up))
        let dr = row - Double(r0), dc = column - Double(c0)
        return (1 - dr) * (1 - dc) * a[r0 * columns + c0] + (1 - dr) * dc * a[r0 * columns + c1]
            + dr * (1 - dc) * a[r1 * columns + c0] + dr * dc * a[r1 * columns + c1]
    }

    // MARK: SciPy's gaussian_filter (order 0, truncate 4, mode "reflect")

    static func gaussianFilter(_ a: [Double], rows: Int, columns: Int, sigma: Double) -> [Double] {
        guard sigma > 1e-15 else { return a }
        let radius = Int(4 * sigma + 0.5)
        var weights = (0...2 * radius).map { exp(-0.5 / (sigma * sigma) * Double(($0 - radius) * ($0 - radius))) }
        let total = weights.reduce(0, +)
        weights = weights.map { $0 / total }
        var out = a
        // Axis 0 (down each column), then axis 1 (along each row).
        var line = [Double](repeating: 0, count: rows)
        for c in 0..<columns {
            for r in 0..<rows { line[r] = out[r * columns + c] }
            let f = correlateSymmetric(line, weights: weights, radius: radius)
            for r in 0..<rows { out[r * columns + c] = f[r] }
        }
        line = [Double](repeating: 0, count: columns)
        for r in 0..<rows {
            for c in 0..<columns { line[c] = out[r * columns + c] }
            let f = correlateSymmetric(line, weights: weights, radius: radius)
            for c in 0..<columns { out[r * columns + c] = f[c] }
        }
        return out
    }

    /// SciPy's correlate1d for a symmetric kernel: w[r]·x[i] + Σⱼ w[r + j]·(x[i + j] + x[i − j]) for j = r down to 1, the line extended by
    /// reflection (d c b a | a b c d | d c b a), repeated for a kernel longer than the line.
    static func correlateSymmetric(_ x: [Double], weights w: [Double], radius: Int) -> [Double] {
        let n = x.count
        func at(_ i: Int) -> Double {
            if n == 1 { return x[0] }
            let period = 2 * n
            var m = i % period
            if m < 0 { m += period }
            return x[m < n ? m : period - 1 - m]
        }
        var out = [Double](repeating: 0, count: n)
        for i in 0..<n {
            var t = x[i] * w[radius]
            // SciPy's order: the outermost pair first, inward (NI_Correlate1D), so that sums equal SciPy's to the bit and ties
            // between neighbours break as they do in py4DSTEM.
            if radius > 0 { for j in stride(from: radius, through: 1, by: -1) { t += (at(i + j) + at(i - j)) * w[radius + j] } }
            out[i] = t
        }
        return out
    }

    // MARK: The radial background

    /// The image less the median of its one-pixel ring about `centre` (ring k holds the pixels at floor(distance) = k; NumPy's
    /// median, the mean of the two middle values for an even count), clipped at 0.
    static func subtractingRadialMedian(_ a: [Double], rows: Int, columns: Int, centre: SIMD2<Double>) -> [Double] {
        var ring = [Int](repeating: 0, count: a.count)
        var maximum = 0
        for r in 0..<rows {
            for c in 0..<columns {
                let k = Int(hypot(Double(r) - centre.x, Double(c) - centre.y))
                ring[r * columns + c] = k
                maximum = max(maximum, k)
            }
        }
        var members = [[Double]](repeating: [], count: maximum + 1)
        for i in a.indices { members[ring[i]].append(a[i]) }
        let median: [Double] = members.map { m in
            guard !m.isEmpty else { return 0 }
            let s = m.sorted(), h = s.count / 2
            return s.count % 2 == 1 ? s[h] : (s[h - 1] + s[h]) / 2
        }
        return a.indices.map { max(a[$0] - median[ring[$0]], 0) }
    }

    // MARK: Checks

    static func pixels(of pattern: DiffractionPattern) throws(Failure) -> [Double] {
        guard pattern.width > 0, pattern.height > 0, pattern.pixels.count == pattern.width * pattern.height else {
            throw .emptyPattern
        }
        guard pattern.pixels.allSatisfy({ $0.isFinite }) else { throw .nonFinitePixel }
        return pattern.pixels.map(Double.init)
    }

    static func validate(_ s: Settings) throws(Failure) {
        func bad(_ text: String) -> Failure { .invalidSettings(text) }
        guard s.sigma.isFinite, s.sigma >= 0, s.sigma <= 1000 else { throw bad("sigma must be in 0...1000") }
        guard s.minimumAbsoluteIntensity.isFinite, s.minimumRelativeIntensity.isFinite, s.minimumRelativeIntensity >= 0 else {
            throw bad("the intensity cuts must be finite, the relative one not negative")
        }
        guard s.relativeToPeak >= 0 else { throw bad("relativeToPeak must not be negative") }
        guard s.minimumSpacing.isFinite, s.minimumSpacing >= 0 else { throw bad("minimumSpacing must be finite and not negative") }
        guard s.edgeBoundary >= 0, s.edgeBoundary <= 1 << 20 else { throw bad("edgeBoundary must be in 0...2^20") }
        guard s.maximumCount >= 1 else { throw bad("maximumCount must be at least 1") }
        guard s.refinementHalfWindow >= 0, s.refinementHalfWindow < 10_000 else { throw bad("refinementHalfWindow must be in 0..<10000") }
        guard s.beamExclusionRadius.isFinite, s.beamExclusionRadius >= 0 else {
            throw bad("beamExclusionRadius must be finite and not negative")
        }
    }
}
