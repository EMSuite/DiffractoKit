//
//  VariantOrientations.swift
//  Role: the orientations of every variant of a precipitate in a matrix, and which of them look alike in one view
//        (DiffractoKit S7, #162).
//
//  The user states the relationship as a paper does, two parallel pairs; `OrientationRelationship.rotation` turns them into
//  the precipitate-to-matrix rotation R (first pair exact, second fixing the turn about it). The variants are the images
//  g R over the matrix group, two being one when they differ by a precipitate operator on the right (Cayron, Acta Cryst.
//  A62 (2006) 21; `OrientationRelationship.variants`, pinned against orix 0.15.0 `Misorientation.equivalent`, P36 and P52).
//  Each variant's orientation is Rᵢᵀ M, M the matrix orientation (columns lab x, y, z in the matrix's Cartesian frame), and
//  its spots are `SpotPatternSimulator`'s, so they carry the spot pattern's pins (P42).
//
//  Which variants look alike is the rule of the held `PhaseVectorMatcher.indistinguishableVariants`: two variants are one
//  group when every spot of either has a spot of the other within the tolerance, through chains of such pairs. Here it runs
//  on the drawn spot patterns (above a cut relative to each variant's strongest spot), so the grouping is what the screen
//  shows. The held rule keeps a flat slab |g_z| ≤ slab; the drawn pattern keeps the Ewald sphere's excitation error, which is
//  not the same for g and its image under a twofold axis in the foil plane. On β″ in Al down [001] at 200 kV that splits
//  two of the held rule's four pairs (their higher-order spots differ); in the flat-sphere limit the two rules agree
//  (`VariantOrientationsTests`).
//

import Foundation
import simd
import DiffractoCore
import DiffractoGeometry

/// One variant of a precipitate in a matrix, in the matrix's pattern frame.
public struct VariantOrientation: Sendable, Equatable {
    /// The variant's place in `OrientationRelationship.variants` (the matrix group's operator order).
    public let index: Int
    /// Precipitate Cartesian to matrix Cartesian for this variant, Rᵢ = gᵢ R.
    public let relation: simd_double3x3
    /// The precipitate's orientation in the lab: Rᵢᵀ times the matrix orientation.
    public let orientation: CrystalOrientation
    /// The matrix vector the variant's first pair lies along (gᵢ applied to the stated matrix vector, same kind), for the
    /// label "b ∥ [100] Al". nil when the image is not a lattice vector of the matrix (an operator outside the matrix's
    /// lattice, which a structure in its standard setting does not have).
    public let firstPairInMatrix: LatticeVector?
    /// Variants of one group predict the same spots in this view; groups are numbered from 0 in the order of their first
    /// variant.
    public let group: Int
    /// The variant's drawn spot pattern.
    public let pattern: SpotPattern

    public init(index: Int, relation: simd_double3x3, orientation: CrystalOrientation, firstPairInMatrix: LatticeVector?,
                group: Int, pattern: SpotPattern) {
        self.index = index
        self.relation = relation
        self.orientation = orientation
        self.firstPairInMatrix = firstPairInMatrix
        self.group = group
        self.pattern = pattern
    }
}

/// The rule that says which variants look alike in one view.
public struct VariantGroupingSettings: Sendable, Equatable {
    /// Two spots nearer than this are one spot, Å⁻¹. Default 0.02, the held rule's.
    public var toleranceInvAngstrom: Double
    /// A spot counts when its intensity is at least this fraction of its variant's strongest spot. Default 0.01, the
    /// design's display cut; 0 counts every spot of the pattern.
    public var relativeIntensityCut: Double
    /// Spots nearer the centre than this are the direct beam and tell nothing apart, Å⁻¹. Default 0.
    public var directBeamRadiusInvAngstrom: Double
    /// Spots farther from the centre are not seen, Å⁻¹. nil (default) keeps every spot of the pattern.
    public var reachInvAngstrom: Double?

    public init(toleranceInvAngstrom: Double = 0.02, relativeIntensityCut: Double = 0.01,
                directBeamRadiusInvAngstrom: Double = 0, reachInvAngstrom: Double? = nil) {
        self.toleranceInvAngstrom = toleranceInvAngstrom
        self.relativeIntensityCut = relativeIntensityCut
        self.directBeamRadiusInvAngstrom = directBeamRadiusInvAngstrom
        self.reachInvAngstrom = reachInvAngstrom
    }

    func validated() throws(VariantOrientationError) {
        guard toleranceInvAngstrom.isFinite, toleranceInvAngstrom >= 0 else {
            throw .invalidGrouping("toleranceInvAngstrom must be a finite number of at least 0")
        }
        guard relativeIntensityCut.isFinite, relativeIntensityCut >= 0, relativeIntensityCut <= 1 else {
            throw .invalidGrouping("relativeIntensityCut must be between 0 and 1")
        }
        guard directBeamRadiusInvAngstrom.isFinite, directBeamRadiusInvAngstrom >= 0 else {
            throw .invalidGrouping("directBeamRadiusInvAngstrom must be a finite number of at least 0")
        }
        if let reach = reachInvAngstrom, !(reach.isFinite && reach >= 0) {
            throw .invalidGrouping("reachInvAngstrom must be a finite number of at least 0")
        }
    }
}

public enum VariantOrientationError: Error, Equatable, Sendable, CustomStringConvertible {
    /// A pair has a vector of no length, or the two vectors of one side are parallel.
    case invalidRelationship
    /// The phase (precipitate or matrix) has no orientation symmetry in its own cell (`PhaseModel.orientationSymmetry(of:)`).
    case noOrientationSymmetry(String)
    /// The spot pattern was refused (an element without a form factor, the voltage, a setting).
    case spotPattern(SpotPatternError)
    /// A grouping setting is out of its range; the text names it.
    case invalidGrouping(String)

    public var description: String {
        switch self {
        case .invalidRelationship:
            return "The orientation relationship needs two pairs of vectors with length, the two vectors of each crystal not parallel."
        case .noOrientationSymmetry(let phase):
            return "The \(phase) has no orientation symmetry in its own cell; restate it in a standard setting."
        case .spotPattern(let e): return e.description
        case .invalidGrouping(let s): return "The variant grouping settings are not valid: \(s)."
        }
    }
}

extension VariantOrientationError: LocalizedError {
    public var errorDescription: String? { description }
}

extension OrientationRelationship {

    /// Every variant of the relationship `first`, `second` (precipitate ∥ matrix, first exact) at `matrixOrientation`, with its
    /// orientation, its label, its drawn pattern and its group in this view.
    ///
    /// - Parameters:
    ///   - matrixOrientation: the matrix's orientation in the pattern (from Index or a drag).
    ///   - volts: the accelerating voltage (200 kV is 200_000).
    ///   - kMax: the precipitate's reflections out to this |g|, Å⁻¹.
    ///   - settings: the spot pattern's excitation envelope and cuts.
    ///   - grouping: which spots count when variants are compared.
    public static func variantOrientations(first: OrientationRelationship, second: OrientationRelationship,
                                           precipitate: CrystalStructure, matrix: CrystalStructure,
                                           matrixOrientation: CrystalOrientation, volts: Double, kMax: Double = 2,
                                           settings: SpotPatternSimulator.Settings = .init(),
                                           grouping: VariantGroupingSettings = .init())
        throws(VariantOrientationError) -> [VariantOrientation] {
        try grouping.validated()
        let relation = try checkedRotation(first, second, precipitate: precipitate, matrix: matrix)
        guard let ms = PhaseModel.orientationSymmetry(of: matrix) else { throw .noOrientationSymmetry("matrix") }
        guard let ps = PhaseModel.orientationSymmetry(of: precipitate) else { throw .noOrientationSymmetry("precipitate") }
        let relations = variants(of: relation, matrix: ms, precipitate: ps)
        let simulator: SpotPatternSimulator
        do { simulator = try SpotPatternSimulator(structure: precipitate, kMax: kMax) } catch { throw .spotPattern(error) }
        var orientations: [CrystalOrientation] = []
        var patterns: [SpotPattern] = []
        for r in relations {
            let o: CrystalOrientation
            do { o = try CrystalOrientation(matrix: r.transpose * matrixOrientation.matrix) } catch {
                throw .invalidRelationship
            }
            orientations.append(o)
            do { patterns.append(try simulator.pattern(orientation: o, volts: volts, settings: settings)) } catch {
                throw .spotPattern(error)
            }
        }
        var groupOf = Array(repeating: 0, count: relations.count)
        for (n, members) in groups(of: patterns, grouping: grouping).enumerated() { for i in members { groupOf[i] = n } }
        let stated = first.matrix.cartesian(in: matrix)
        return relations.indices.map { i in
            let operatorInMatrix = relations[i] * relation.transpose
            return VariantOrientation(index: i, relation: relations[i], orientation: orientations[i],
                                      firstPairInMatrix: latticeVector(operatorInMatrix * stated, like: first.matrix, in: matrix),
                                      group: groupOf[i], pattern: patterns[i])
        }
    }

    /// The variants grouped by the spots they draw: the held rule of `PhaseVectorMatcher.indistinguishableVariants` on the
    /// spots of each pattern that `grouping` counts. Each group ascending, groups ordered by their first variant; every index
    /// of `patterns` is in one. Settings out of range count every spot (validate them with `variantOrientations`).
    public static func groups(of patterns: [SpotPattern], grouping: VariantGroupingSettings = .init()) -> [[Int]] {
        let spots: [[SIMD2<Double>]] = patterns.map { pattern in
            let strongest = pattern.spots.map(\.intensity).max() ?? 0
            let floor = strongest * grouping.relativeIntensityCut
            return pattern.spots.compactMap { spot -> SIMD2<Double>? in
                let r = simd_length(spot.position)
                guard spot.intensity >= floor, r > grouping.directBeamRadiusInvAngstrom,
                      r <= (grouping.reachInvAngstrom ?? .infinity) else { return nil }
                return spot.position
            }
        }
        return VariantSpotGrouping.groups(spots, toleranceInvAngstrom: grouping.toleranceInvAngstrom)
    }

    /// How far the second pair is from parallel once the first is made exact, radians: the angle between the precipitate's
    /// second vector carried into the matrix and the matrix's second vector. The planted β″ relation ([100] ∥ [320] Al,
    /// [001] ∥ [1 -3 0] Al) gives 105.3° − arccos(−3/√130) = 0.045°.
    public static func secondPairMisfit(first: OrientationRelationship, second: OrientationRelationship,
                                        precipitate: CrystalStructure, matrix: CrystalStructure)
        throws(VariantOrientationError) -> Double {
        let r = try checkedRotation(first, second, precipitate: precipitate, matrix: matrix)
        let a = r * second.candidate.cartesian(in: precipitate), b = second.matrix.cartesian(in: matrix)
        return atan2(simd_length(simd_cross(a, b)), simd_dot(a, b))
    }

    static func checkedRotation(_ first: OrientationRelationship, _ second: OrientationRelationship,
                                precipitate: CrystalStructure, matrix: CrystalStructure)
        throws(VariantOrientationError) -> simd_double3x3 {
        let r: simd_double3x3
        do { r = try rotation(first, second, precipitate: precipitate, matrix: matrix) } catch { throw .invalidRelationship }
        guard Misorientation.isRotation(r) else { throw .invalidRelationship }
        return r
    }

    /// `v` (Cartesian, in `crystal`'s frame) as integer indices of the same kind as `kind`, or nil when an index is more than
    /// 1e-6 from an integer.
    static func latticeVector(_ v: SIMD3<Double>, like kind: LatticeVector, in crystal: CrystalStructure) -> LatticeVector? {
        let real = crystal.cell.realLattice
        let indices: SIMD3<Double>
        switch kind {
        case .plane:
            // g = h a* + k b* + l c*, and a*·a = 1: h = g·a.
            indices = SIMD3(simd_dot(v, real[0]), simd_dot(v, real[1]), simd_dot(v, real[2]))
        case .direction:
            // v = u a + v b + w c: solve with the rows a, b, c as columns.
            let columns = simd_double3x3(columns: (real[0], real[1], real[2]))
            indices = columns.inverse * v
        }
        let rounded = indices.rounded(.toNearestOrAwayFromZero)
        guard indices.x.isFinite, indices.y.isFinite, indices.z.isFinite,
              simd_reduce_max(simd_abs(indices - rounded)) <= 1e-6 else { return nil }
        let ints = SIMD3<Int>(Int(rounded.x), Int(rounded.y), Int(rounded.z))
        switch kind {
        case .plane: return .plane(ints)
        case .direction: return .direction(ints)
        }
    }
}

/// The grouping rule shared by `PhaseVectorMatcher.indistinguishableVariants` and `OrientationRelationship.groups(of:grouping:)`.
enum VariantSpotGrouping {
    /// Union-find over the pairs whose spot sets cover each other within `toleranceInvAngstrom`.
    static func groups(_ spots: [[SIMD2<Double>]], toleranceInvAngstrom: Double) -> [[Int]] {
        func covers(_ a: [SIMD2<Double>], _ b: [SIMD2<Double>]) -> Bool {
            a.allSatisfy { s in b.contains { simd_length($0 - s) <= toleranceInvAngstrom } }
        }
        var parent = Array(spots.indices)
        func root(_ i: Int) -> Int { var r = i; while parent[r] != r { r = parent[r] }; return r }
        for i in spots.indices {
            for j in spots.indices where j > i {
                if covers(spots[i], spots[j]) && covers(spots[j], spots[i]) {
                    let a = root(i), b = root(j)
                    if a != b { parent[max(a, b)] = min(a, b) }
                }
            }
        }
        var groups: [Int: [Int]] = [:]
        for i in spots.indices { groups[root(i), default: []].append(i) }
        return groups.keys.sorted().map { groups[$0]!.sorted() }
    }
}
