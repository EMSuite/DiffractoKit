//
//  PowderPattern+Rings.swift
//  Role: the electron powder pattern (rings) from the held kinematic ring list, a drawn profile of any powder pattern, and
//        which peaks of several shown phases overlap and which are diagnostic (DiffractoKit S10, #165).
//
//  Rings: one peak per ring of `CrystalStructure.kinematicReflections` (DiffractoCore's families of the reflections above
//  the |F| cut, Σ|F|² over the members, that is multiplicity × |F|²; pin P4), at q = |g| in Å⁻¹, scaled so the strongest is
//  100. Nothing is new in that number but the scaling; the families' orbits are the peak's labels, each with its members.
//  Against a measured, azimuthally averaged profile a ring's height goes as Σ|F|² · d² (the integrated ring carries m|F|²d,
//  the profile divides by the circumference; Vainshtein, Structure Analysis by Electron Diffraction): `profileIntensity`,
//  also scaled to 100 (#223). The overlap marks rank and floor by it where a peak has one.
//  Profile: a sum of unit-area Gaussians, one per peak, weighted by the peak's relative intensity, with a typed full width
//  at half maximum. A drawing aid, no instrument function and no texture; tested as its closed form.
//  Marks: a rule on the peak lists of the shown phases. A peak is overlapped when another shown phase (same axis) has a
//  peak of at least `floor` (of 100) within the resolution of it; a phase's diagnostic peak is its strongest peak that is
//  not overlapped. "Strength" is `profileIntensity` where a peak has one (electron rings), else `intensity`. No Python package computes this; it is checked exact on planted lists and against the design's
//  numbers for Al + Mg2Si + Si with Cu Kα at 0.1° (38.51°, 40.16°, 28.47°).
//

import Foundation
import DiffractoCore

extension PowderPattern {

    /// The most points `profile(width:range:step:)` draws: 10 million.
    public static let maximumProfilePoints = 10_000_000

    /// The electron ring pattern of `structure`: the rings of `kinematicReflections(kMax:)` at q = |g| (Å⁻¹, no 2π), their
    /// Σ|F|² (multiplicity × |F|², no Debye-Waller factor, as diffsims and py4DSTEM) as `absoluteIntensity` (Å⁻⁴) and scaled so
    /// the strongest ring is 100 as `intensity`, and their height in an azimuthally averaged profile, Σ|F|² · d² scaled the
    /// same way, as `profileIntensity` (the one to compare with a measured profile). `wavelengthAngstroms` is recorded with the
    /// pattern (it does not move a ring in q). The ring list's refusals come as `invalidRingList` with its reason. Labels
    /// are the rings' orbits, (h k i l) where the ring list names a hexagonal cell's rings in four indices.
    public static func electron(structure: CrystalStructure, wavelengthAngstroms: Double, kMax: Double = 2,
                                cut: Double = StructureFactors.defaultCut,
                                validation: CrystalStructure.ValidationPolicy = .default)
        throws(PowderPatternError) -> PowderPattern {
        guard wavelengthAngstroms.isFinite, wavelengthAngstroms > 0 else { throw .invalidWavelength(wavelengthAngstroms) }
        let list: KinematicReflectionList
        do {
            list = try structure.kinematicReflections(kMax: kMax, cut: cut, validation: validation, retainingAll: false)
        } catch {
            throw .invalidRingList(error.description)
        }
        let strongest = list.families.map(\.intensity).max() ?? 0
        // The profile's height of a ring: Σ|F|² · d² (the integrated m|F|²d over the circumference 2π/d).
        let profile = list.families.map { $0.intensity * $0.family.d * $0.family.d }
        let strongestProfile = profile.max() ?? 0
        let rotations = structure.laueGroup.rotations
        let fourIndices = ReflectionList.usesMillerBravais(laueClass: structure.laueClass, cell: structure.cell)
        let peaks: [Peak] = list.families.indices.map { i in
            let ring = list.families[i]
            return Peak(families: orbitFamilies(ring.family, rotations: rotations, fourIndices: fourIndices), d: ring.family.d,
                        position: ring.family.g, intensity: strongest > 0 ? ring.intensity / strongest * 100 : 0,
                        absoluteIntensity: ring.intensity,
                        profileIntensity: strongestProfile > 0 ? profile[i] / strongestProfile * 100 : 0)
        }
        return PowderPattern(axis: .scatteringVector, wavelengthAngstroms: wavelengthAngstroms, range: 0...kMax, peaks: peaks)
    }

    /// The ring's orbits as labels, each with the number of its members (a member joins the first orbit an operator maps to
    /// it).
    static func orbitFamilies(_ ring: ReflectionFamily, rotations: [IntegerMatrix], fourIndices: Bool = false) -> [Family] {
        var counts = Array(repeating: 0, count: ring.orbits.count)
        for m in ring.members {
            let i = ring.orbits.firstIndex { o in
                rotations.contains { o.applying($0) == m }
            } ?? 0
            if !counts.isEmpty { counts[i] += 1 }
        }
        return ring.orbits.indices.map { i in
            let o = ring.orbits[i]
            return Family(indices: fourIndices ? [o.h, o.k, -o.h - o.k, o.l] : [o.h, o.k, o.l], multiplicity: counts[i])
        }
    }

    /// The drawn profile: at x = lower + i · step for i = 0 … ⌊(upper − lower) / step⌋, Σ over the peaks of
    /// I · exp(−(x − p)² / 2σ²) / (σ √(2π)), σ = `width` / (2 √(2 ln 2)) (`width` is the full width at half maximum, in the
    /// pattern's unit). Empty for a width or step that is not a finite positive number, a range that is not finite, or more
    /// than `maximumProfilePoints` points.
    public func profile(width: Double, range: ClosedRange<Double>, step: Double) -> [Double] {
        guard width.isFinite, width > 0, step.isFinite, step > 0, range.lowerBound.isFinite, range.upperBound.isFinite else {
            return []
        }
        let sigma = width / (2 * (2 * log(2.0)).squareRoot())
        let norm = 1 / (sigma * (2 * Double.pi).squareRoot())
        let n = ((range.upperBound - range.lowerBound) / step + 1e-9).rounded(.down)
        guard n.isFinite, n >= 0, n < Double(Self.maximumProfilePoints) else { return [] }
        let count = Int(n) + 1
        return (0 ..< count).map { i in
            let x = range.lowerBound + Double(i) * step
            var y = 0.0
            for p in peaks {
                let u = (x - p.position) / sigma
                y += p.intensity * norm * exp(-0.5 * u * u)
            }
            return y
        }
    }
}

extension PowderPattern.Peak {
    /// What the marks rank and floor by: `profileIntensity` for an electron ring, else `intensity` (both of 100).
    var strength: Double { profileIntensity ?? intensity }
}

/// How near two peaks must be to be one in a drawing.
public enum PowderResolution: Sendable, Equatable {
    /// Within this many units of the axis (degrees of 2θ, or Å⁻¹).
    case absolute(Double)
    /// Within this fraction of the peak's position (0.01: 1 % of q, the tolerance Index uses for rings).
    case fraction(Double)

    func window(at position: Double) -> Double {
        switch self {
        case .absolute(let w): return w
        case .fraction(let f): return f * abs(position)
        }
    }
}

/// A peak of another shown pattern.
public struct PowderPeakReference: Sendable, Equatable, Hashable {
    public let pattern: Int
    public let peak: Int

    public init(pattern: Int, peak: Int) {
        self.pattern = pattern
        self.peak = peak
    }
}

/// What a peak is among the shown phases.
public struct PowderPeakMark: Sendable, Equatable {
    /// The peaks of other patterns within the resolution, of at least the floor, in pattern then peak order. Empty when
    /// the peak stands alone.
    public let partners: [PowderPeakReference]
    /// The strongest peak of its pattern that is not overlapped.
    public let isDiagnostic: Bool
    public var isOverlapped: Bool { !partners.isEmpty }

    public init(partners: [PowderPeakReference], isDiagnostic: Bool) {
        self.partners = partners
        self.isDiagnostic = isDiagnostic
    }
}

extension PowderPattern {

    /// The marks of every peak of every pattern, in the patterns' and peaks' order. A peak's strength, for the floor and for
    /// the diagnostic peak, is its `profileIntensity` where it has one (electron rings: what a measured profile shows), else
    /// its `intensity`. Patterns on different axes never
    /// overlap. A pattern whose peaks are all overlapped has no diagnostic peak; of equal strongest peaks the first is
    /// diagnostic. A resolution or floor that is not a finite number of at least 0 marks nothing as overlapped.
    public static func marks(_ patterns: [PowderPattern], resolution: PowderResolution, floor: Double = 1) -> [[PowderPeakMark]] {
        let valid: Bool = {
            switch resolution {
            case .absolute(let w), .fraction(let w): return w.isFinite && w >= 0 && floor.isFinite && floor >= 0
            }
        }()
        return patterns.indices.map { i in
            let pattern = patterns[i]
            let partners: [[PowderPeakReference]] = pattern.peaks.map { peak in
                guard valid else { return [] }
                let w = resolution.window(at: peak.position)
                var found: [PowderPeakReference] = []
                for j in patterns.indices where j != i && patterns[j].axis == pattern.axis {
                    for (k, other) in patterns[j].peaks.enumerated()
                    where other.strength >= floor && abs(other.position - peak.position) <= w {
                        found.append(PowderPeakReference(pattern: j, peak: k))
                    }
                }
                return found
            }
            var diagnostic: Int? = nil
            for k in pattern.peaks.indices where partners[k].isEmpty {
                if diagnostic == nil || pattern.peaks[k].strength > pattern.peaks[diagnostic!].strength { diagnostic = k }
            }
            return pattern.peaks.indices.map { PowderPeakMark(partners: partners[$0], isDiagnostic: $0 == diagnostic) }
        }
    }
}
