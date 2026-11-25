//
//  PowderPattern.swift
//  Role: the line pattern of a powder: X-ray peaks as pymatgen's `XRDCalculator.get_pattern` computes them (pymatgen
//        2026.9.24, pymatgen-core 2026.10.16), ported step for step so the peaks, their 2θ, intensities and hkl labels
//        agree with pymatgen's (pin P49).
//
//        1. Every reciprocal-lattice point g = h b₁ + k b₂ + l b₃ (no 2π) with |g|² ≤ r_max² + 1e-8 (pymatgen's
//           `find_points_in_spheres` tolerance), r = 2 sin(2θ/2) / λ at the ends of the 2θ range; the lower end filters
//           |g| ≥ r_min only when it is not 0. One Friedel half-space is kept ((h, k, l) lexicographically positive),
//           ordered by |g|, then by descending h, k, l.
//        2. F(g) = Σ_elements f_e(s) exp(−B_e s²) Σ_sites occ exp(2πi h·r), s = |g| / 2, f from ``XRayScatteringFactors``;
//           the elements in the order they first appear in the sites, as pymatgen groups them.
//        3. I = 2 |F|² LP, LP = (1 + cos² 2θ) / (sin² θ cos θ), θ = asin(λ|g|/2) (2: the −g half).
//        4. Reflections whose round(2θ / 1e-5) (half to even) is equal are one peak, at the 2θ and d of the first; the
//           hkl labels are pymatgen's `get_unique_families` of every member and its negative (four indices h k i l on a
//           lattice pymatgen calls hexagonal: two angles within 5° of 90, one within 5° of 60 or 120, the two lengths at
//           the right angles within 0.01 Å).
//        5. Peaks at or below 1e-3 % of the strongest are dropped; intensities are scaled so the strongest is 100.
//
//        A line pattern for planning and comparison: no instrument profile, no texture, no absorption, no anomalous
//        dispersion; thermal motion only as typed Debye-Waller factors.
//        DEPARTURE from pymatgen: a range that holds no reflection gives an empty pattern, where pymatgen raises
//        (max() of an empty sequence). As in pymatgen, the scale is relative to the range: a range that holds only
//        reflections the structure forbids shows their rounding residue scaled to 100 (`absoluteIntensity` tells).
//

import Foundation
import simd
import DiffractoCore

/// The peaks of a powder (Debye-Scherrer) pattern.
public struct PowderPattern: Sendable, Equatable {

    /// What a peak's `position` is.
    public enum Axis: Sendable, Equatable, Hashable {
        /// 2θ in degrees (an X-ray pattern).
        case twoThetaDegrees
        /// q = |g| = 1/d in Å⁻¹, without 2π (an electron ring pattern).
        case scatteringVector
    }

    /// One group of reflections in a peak, named as pymatgen names it.
    public struct Family: Sendable, Equatable, Hashable, CustomStringConvertible {
        /// (h, k, l), or (h, k, i, l) with i = −h − k on a hexagonal lattice: the largest member of the group
        /// (lexicographic), the group being every member whose |indices| are a permutation of each other.
        public let indices: [Int]
        /// The number of members of the group (hkl and −hkl both counted).
        public let multiplicity: Int

        public init(indices: [Int], multiplicity: Int) {
            self.indices = indices
            self.multiplicity = multiplicity
        }

        public var description: String { indices.map(String.init).joined(separator: " ") }
    }

    /// One peak: reflections at one position.
    public struct Peak: Sendable, Equatable {
        /// The groups of reflections at this position, in the order their members first appear (pymatgen's order up to
        /// round-off ties in |g|).
        public let families: [Family]
        /// d of the peak's first reflection (Å).
        public let d: Double
        /// 2θ in degrees, or q in Å⁻¹ (``PowderPattern/axis``).
        public let position: Double
        /// Relative intensity: the strongest peak of the pattern is 100. For electron rings the scaled Σ|F|² per ring
        /// (multiplicity × |F|²); see `profileIntensity` for the height in a measured profile.
        public let intensity: Double
        /// The intensity before scaling. X-rays: 2 Σ|F|² LP over the merged reflections (electrons², LP included), what
        /// pymatgen returns with `scaled=False`. Electron rings (`electron(structure:wavelengthAngstroms:)`): Σ|F|² over the
        /// ring's members (multiplicity × |F|², Å⁻⁴), diffsims' ring intensity. Comparable between patterns computed the same way.
        public let absoluteIntensity: Double
        /// Electron rings only: the ring's height in the azimuthally averaged profile of a polycrystalline film, Σ|F|² · d²
        /// scaled so the strongest ring is 100 (kinematic, no Debye-Waller factor). The integrated ring carries m|F|²d and the
        /// profile divides it by the circumference 2π/d (Vainshtein, Structure Analysis by Electron Diffraction, 1964, ch. 2),
        /// so this is the number to set against a measured profile: Au 220 is 26 against 111's 100 here, 68 in `intensity`.
        /// Nil for an X-ray pattern, whose `intensity` already carries its Lorentz-polarisation factor.
        public var profileIntensity: Double? = nil

        /// Σ multiplicity over the families.
        public var multiplicity: Int { families.reduce(0) { $0 + $1.multiplicity } }
    }

    /// What `position` measures.
    public let axis: Axis
    /// λ in Å.
    public let wavelengthAngstroms: Double
    /// The range the peaks were computed in (2θ in degrees for X-rays).
    public let range: ClosedRange<Double>
    /// The peaks, in ascending position.
    public let peaks: [Peak]

    public init(axis: Axis, wavelengthAngstroms: Double, range: ClosedRange<Double>, peaks: [Peak]) {
        self.axis = axis
        self.wavelengthAngstroms = wavelengthAngstroms
        self.range = range
        self.peaks = peaks
    }

    /// The strongest peak (intensity 100), nil when the pattern is empty.
    public var strongest: Peak? { peaks.max { $0.intensity < $1.intensity } }
}

/// Why a powder pattern was not computed.
public enum PowderPatternError: Error, Equatable, Sendable, CustomStringConvertible {
    /// Atomic numbers the scattering table has no entry for.
    case unsupportedElements([Int])
    /// The structure fails `CrystalStructure.validationIssues`.
    case invalidStructure([String])
    /// λ must be a finite number of Å above 0.
    case invalidWavelength(Double)
    /// The range must be finite, start at 0 or above, end above its start and at most at 180°.
    case invalidRange(lower: Double, upper: Double)
    /// A Debye-Waller factor must be a finite number (Å²); keyed by atomic number.
    case invalidDebyeWallerFactor(z: Int, value: Double)
    /// The cell at this wavelength holds more reflections than `ReflectionList.maximumReflections`.
    case tooManyReflections(estimated: Double)
    /// An intensity came out non-finite at this reflection.
    case nonFiniteIntensity(MillerIndex)
    /// The electron ring list was refused (`kinematicReflections`); the text is its reason.
    case invalidRingList(String)

    public var description: String {
        switch self {
        case .unsupportedElements(let z):
            return "The X-ray scattering table has no entry for Z = \(z.map(String.init).joined(separator: ", ")), so this "
                + "crystal's powder pattern would be wrong."
        case .invalidStructure(let issues):
            return "The structure is not valid, so its powder pattern would be meaningless: \(issues.joined(separator: "; "))."
        case .invalidWavelength(let w):
            return "The wavelength must be a number of Å above 0; \(w) is not."
        case .invalidRange(let lo, let hi):
            return "The 2θ range must run from 0° or more up to at most 180°, its end above its start; \(lo) to \(hi) does not."
        case .invalidDebyeWallerFactor(let z, let b):
            return "The Debye-Waller factor of Z = \(z) must be a finite number of Å²; \(b) is not."
        case .tooManyReflections(let n):
            return String(format: "This cell at this wavelength holds about %.3g reflections, more than a pattern is computed "
                          + "for (%ld); narrow the 2θ range or use a longer wavelength.", n, ReflectionList.maximumReflections)
        case .nonFiniteIntensity(let hkl):
            return "The intensity of \(hkl) is not a finite number."
        case .invalidRingList(let reason):
            return "The ring pattern cannot be computed: \(reason)"
        }
    }
}

/// `localizedDescription` (what an app shows) is the same text as `description`.
extension PowderPatternError: LocalizedError {
    public var errorDescription: String? { description }
}

extension PowderPattern {

    /// pymatgen's peak-merging tolerance, `TWO_THETA_TOL` (degrees).
    public static let twoThetaMergeTolerance = 1e-5
    /// pymatgen's floor, `SCALED_INTENSITY_TOL`: peaks at or below this percentage of the strongest are dropped.
    public static let intensityFloorPercent = 1e-3

    /// The X-ray powder pattern of `structure`, as pymatgen's `XRDCalculator(wavelength, debye_waller_factors)
    /// .get_pattern(structure, two_theta_range)` computes it. Defaults are pymatgen's: Cu Kα, 0 to 90°, no thermal factor.
    /// `debyeWallerFactors`: B in Å² by atomic number, applied as exp(−B s²).
    public static func xray(structure: CrystalStructure, radiation: XRayRadiation = .cuKAlpha,
                            twoThetaRange: ClosedRange<Double> = 0...90,
                            debyeWallerFactors: [Int: Double] = [:],
                            validation: CrystalStructure.ValidationPolicy = .default)
        throws(PowderPatternError) -> PowderPattern {
        let lambda = radiation.wavelengthAngstroms
        guard lambda.isFinite, lambda > 0 else { throw .invalidWavelength(lambda) }
        let lo = twoThetaRange.lowerBound, hi = twoThetaRange.upperBound
        guard lo.isFinite, hi.isFinite, lo >= 0, hi > lo, hi <= 180 else { throw .invalidRange(lower: lo, upper: hi) }
        for (z, b) in debyeWallerFactors.sorted(by: { $0.key < $1.key }) where !b.isFinite {
            throw .invalidDebyeWallerFactor(z: z, value: b)
        }
        let structureSites = structure.sites
        let missing = Array(Set(structureSites.map(\.z).filter { !XRayScatteringFactors.isSupported(z: $0) })).sorted()
        guard missing.isEmpty else { throw .unsupportedElements(missing) }
        let issues = structure.validationIssues(using: validation)
        guard issues.isEmpty else { throw .invalidStructure(issues.map(\.message)) }

        // 1. the reciprocal points in the shell
        let degToRad = Double.pi / 180
        let minR = 2 * sin((lo / 2) * degToRad) / lambda
        let maxR = 2 * sin((hi / 2) * degToRad) / lambda
        let maxR2 = maxR * maxR + 1e-8
        let cell = structure.cell
        let estimate = 4.0 / 3.0 * Double.pi * pow(maxR2.squareRoot(), 3) * cell.volume
        guard estimate.isFinite, estimate <= Double(ReflectionList.maximumReflections) else {
            throw .tooManyReflections(estimated: estimate)
        }
        let b = cell.reciprocalLattice
        let reach = maxR2.squareRoot()
        let tile = cell.realLattice.map { Int((reach * length($0)).rounded(.up)) }
        var points: [(hkl: MillerIndex, g: Double)] = []
        for h in 0...tile[0] {
            for k in (h == 0 ? 0 : -tile[1])...tile[1] {
                for l in -tile[2]...tile[2] {
                    // the Friedel half-space: (h, k, l) lexicographically positive
                    guard h > 0 || k > 0 || (k == 0 && l > 0) else { continue }
                    let g = Double(h) * b[0] + Double(k) * b[1] + Double(l) * b[2]
                    let g2 = dot(g, g)
                    guard g2 <= maxR2 else { continue }
                    let gLen = g2.squareRoot()
                    if minR != 0, gLen < minR { continue }
                    points.append((MillerIndex(h, k, l), gLen))
                }
            }
        }
        points.sort { x, y in
            if x.g != y.g { return x.g < y.g }
            return (x.hkl.h, x.hkl.k, x.hkl.l) > (y.hkl.h, y.hkl.k, y.hkl.l)
        }

        // 2. the sites grouped by element, in the order of first appearance
        var order: [Int] = []
        var groups: [Int: [(position: SIMD3<Double>, occupancy: Double)]] = [:]
        for site in structureSites {
            if groups[site.z] == nil { order.append(site.z); groups[site.z] = [] }
            groups[site.z]!.append((site.fractional, site.occupancy))
        }
        let isHexagonal = pymatgenIsHexagonal(cell)

        // 3 and 4. intensities, merged by 2θ bin
        struct Bin { var intensity: Double; var members: [[Int]]; let twoTheta: Double; let d: Double }
        var bins: [Int: Bin] = [:]
        for p in points {
            let theta = asin(min(max(lambda * p.g / 2, -1), 1))
            let half = p.g / 2
            let s2 = half * half
            var re = 0.0, im = 0.0
            for z in order {
                // pymatgen's form factor, Z − 41.78214 s² Σᵢ aᵢ exp(−bᵢ s²): the kit's one copy, in the same order of operations
                // (the element is known to be in the table; a non-finite s² stays not a number, as it was).
                var fs = XRayScatteringFactors.f(z: z, sSquared: s2) ?? .nan
                if let bw = debyeWallerFactors[z], bw != 0 { fs *= exp(-bw * s2) }
                var sre = 0.0, sim = 0.0
                for site in groups[z]! {
                    let phase = 2 * Double.pi * (Double(p.hkl.h) * site.position.x + Double(p.hkl.k) * site.position.y
                                                 + Double(p.hkl.l) * site.position.z)
                    sre += cos(phase) * site.occupancy
                    sim += sin(phase) * site.occupancy
                }
                re += fs * sre
                im += fs * sim
            }
            let i = re * re + im * im
            let cos2t = cos(2 * theta), sint = sin(theta), cost = cos(theta)
            let lorentz = (1 + cos2t * cos2t) / (sint * sint * cost)
            let intensity = 2 * i * lorentz
            guard intensity.isFinite else { throw .nonFiniteIntensity(p.hkl) }
            let twoTheta = (2 * theta) * (180 / Double.pi)
            let key = Int((twoTheta / twoThetaMergeTolerance).rounded(.toNearestOrEven))
            var hkl = [p.hkl.h, p.hkl.k, p.hkl.l]
            var neg = hkl.map { -$0 }
            if isHexagonal {
                hkl = [hkl[0], hkl[1], -hkl[0] - hkl[1], hkl[2]]
                neg = [neg[0], neg[1], -neg[0] - neg[1], neg[2]]
            }
            if bins[key] != nil {
                bins[key]!.intensity += intensity
                bins[key]!.members.append(hkl)
                bins[key]!.members.append(neg)
            } else {
                bins[key] = Bin(intensity: intensity, members: [hkl, neg], twoTheta: twoTheta, d: 1 / p.g)
            }
        }

        // 5. the floor and the scale
        guard let maxIntensity = bins.values.map(\.intensity).max(), maxIntensity > 0 else {
            return PowderPattern(axis: .twoThetaDegrees, wavelengthAngstroms: lambda, range: twoThetaRange, peaks: [])
        }
        var kept: [Bin] = []
        for key in bins.keys.sorted() {
            let bin = bins[key]!
            if bin.intensity / maxIntensity * 100 > intensityFloorPercent { kept.append(bin) }
        }
        let factor = kept.map(\.intensity).max() ?? maxIntensity
        let peaks = kept.map { bin in
            Peak(families: uniqueFamilies(bin.members), d: bin.d, position: bin.twoTheta,
                 intensity: bin.intensity / (factor / 100), absoluteIntensity: bin.intensity)
        }
        return PowderPattern(axis: .twoThetaDegrees, wavelengthAngstroms: lambda, range: twoThetaRange, peaks: peaks)
    }

    /// pymatgen's `Lattice.is_hexagonal()` (hex_angle_tol 5°, hex_length_tol 0.01 Å, absolute).
    static func pymatgenIsHexagonal(_ cell: UnitCell) -> Bool {
        let lengths = [cell.a, cell.b, cell.c]
        let angles = [cell.alphaDeg, cell.betaDeg, cell.gammaDeg]
        let right = (0..<3).filter { abs(angles[$0] - 90) <= 5 }
        let hex = (0..<3).filter { abs(angles[$0] - 60) <= 5 || abs(angles[$0] - 120) <= 5 }
        return right.count == 2 && hex.count == 1 && abs(lengths[right[0]] - lengths[right[1]]) <= 0.01
    }

    /// pymatgen's `get_unique_families`: members grouped when their |indices| are permutations of each other, each group
    /// named by its largest member (lexicographic), in the order the groups first appear.
    static func uniqueFamilies(_ members: [[Int]]) -> [Family] {
        var signatures: [[Int]] = []
        var groups: [[[Int]]] = []
        for m in members {
            let signature = m.map { abs($0) }.sorted()
            if let i = signatures.firstIndex(of: signature) {
                groups[i].append(m)
            } else {
                signatures.append(signature)
                groups.append([m])
            }
        }
        return groups.map { group in
            let largest = group.max { $0.lexicographicallyPrecedes($1) }!
            return Family(indices: largest, multiplicity: group.count)
        }
    }
}
