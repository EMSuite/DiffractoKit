//
//  VariantOrientationsTests.swift
//  DiffractoKit S7 (#162): the variant orientations of a relationship stated as two pairs, their labels, the second pair's
//  misfit, and the grouping on drawn spot patterns against the held `indistinguishableVariants`. Planted truth: β″ (Andersen
//  et al. 1998), θ′ and T1 in aluminium (PrecipitateFixtures).
//

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoGeometry
import DiffractoOrientation
import DiffractoKitTestSupport

@Suite struct VariantOrientationsTests {
    static let al = CrystalPresets.aluminum
    static let beta = CrystalPresets.betaDoublePrime

    static func pair(_ p: LatticeVector, _ m: LatticeVector) -> OrientationRelationship {
        OrientationRelationship(candidate: p, matrix: m)
    }

    /// β″: [100] ∥ [3 2 0] Al exact, [001] ∥ [1 -3 0] Al.
    static let betaFirst = pair(.direction(SIMD3(1, 0, 0)), .direction(SIMD3(3, 2, 0)))
    static let betaSecond = pair(.direction(SIMD3(0, 0, 1)), .direction(SIMD3(1, -3, 0)))
    /// θ′: (001) ∥ (001) Al, [100] ∥ [100] Al.
    static let thetaFirst = pair(.plane(SIMD3(0, 0, 1)), .plane(SIMD3(0, 0, 1)))
    static let thetaSecond = pair(.direction(SIMD3(1, 0, 0)), .direction(SIMD3(1, 0, 0)))
    /// T1: (0001) ∥ (111) Al, [2 1 0] (<10-10>) ∥ [1 -1 0] Al.
    static let t1First = pair(.plane(SIMD3(0, 0, 1)), .plane(SIMD3(1, 1, 1)))
    static let t1Second = pair(.direction(SIMD3(2, 1, 0)), .direction(SIMD3(1, -1, 0)))

    static func alOrientation(_ uvw: SIMD3<Int> = SIMD3(0, 0, 1), degrees: Double = 12) throws -> CrystalOrientation {
        try CrystalOrientation(zoneAxis: uvw, inPlaneAngle: degrees * .pi / 180, cell: al.cell)
    }

    static func gap(_ a: simd_double3x3, _ b: simd_double3x3) -> Double {
        var g = 0.0
        for c in 0 ..< 3 { for r in 0 ..< 3 { g = max(g, abs(a[c][r] - b[c][r])) } }
        return g
    }

    static func indices(_ v: LatticeVector?) -> SIMD3<Int>? {
        switch v {
        case .plane(let i)?, .direction(let i)?: return i
        case nil: return nil
        }
    }

    /// Up to sign, as a sorted triple of absolute values and as a canonical vector (first non-zero entry positive).
    static func canonical(_ i: SIMD3<Int>) -> SIMD3<Int> {
        let first = [i.x, i.y, i.z].first { $0 != 0 } ?? 1
        return first < 0 ? SIMD3(-i.x, -i.y, -i.z) : i
    }

    /// β″: twelve variants, each orientation the composed rotation Rᵢᵀ M, the relations those of `variants(of:)`, the labels
    /// the twelve images of [3 2 0], and the second pair 0.045° off.
    /// Mutation: the orientation taken as Rᵢ M (not transposed) fails the composed-rotation check.
    @Test func betaDoublePrimeHasTwelveVariantsEachTheComposedRotation() throws {
        let m = try Self.alOrientation()
        let vs = try OrientationRelationship.variantOrientations(first: Self.betaFirst, second: Self.betaSecond, precipitate: Self.beta,
                                                                 matrix: Self.al, matrixOrientation: m, volts: 200_000)
        #expect(vs.count == 12)
        let r = try OrientationRelationship.rotation(Self.betaFirst, Self.betaSecond, precipitate: Self.beta, matrix: Self.al)
        let held = OrientationRelationship.variants(of: r, matrix: CrystalOrientationSymmetry(laueClass: .mThreeBarM),
                                                    precipitate: CrystalOrientationSymmetry(laueClass: .twoOverM))
        #expect(vs.map(\.index) == Array(0 ..< 12))
        for v in vs {
            #expect(Self.gap(v.relation, held[v.index]) < 1e-12)
            #expect(Self.gap(v.orientation.matrix, v.relation.transpose * m.matrix) < 1e-12)
            #expect(v.pattern.orientation == v.orientation && !v.pattern.spots.isEmpty)
            let label = try #require(Self.indices(v.firstPairInMatrix))
            #expect([abs(label.x), abs(label.y), abs(label.z)].sorted() == [0, 2, 3], "\(label)")
            // The label is where the variant puts a of β″.
            let a = simd_normalize(v.relation * Self.beta.cell.realLattice[0])
            let l = simd_normalize(SIMD3(Double(label.x), Double(label.y), Double(label.z)))
            #expect(simd_dot(a, l) > 1 - 1e-12)
        }
        #expect(Set(vs.compactMap { Self.indices($0.firstPairInMatrix) }).count == 12)
        let misfit = try OrientationRelationship.secondPairMisfit(first: Self.betaFirst, second: Self.betaSecond,
                                                                  precipitate: Self.beta, matrix: Self.al)
        #expect(abs(misfit - (105.3 * .pi / 180 - acos(-3 / 130.0.squareRoot()))) < 1e-12)
        #expect(abs(misfit * 180 / .pi - 0.0449) < 1e-4)
    }

    /// θ′: three variants, one per cube axis, their (001) on (100), (010) and (001) of Al; the second pair exact.
    @Test func thetaPrimeHasThreeVariantsOnePerCubeAxis() throws {
        let theta = PrecipitateFixtures.thetaPrime
        #expect(theta.laueClass == .fourOverMmm)
        let vs = try OrientationRelationship.variantOrientations(first: Self.thetaFirst, second: Self.thetaSecond, precipitate: theta,
                                                                 matrix: Self.al, matrixOrientation: try Self.alOrientation(), volts: 200_000)
        #expect(vs.count == 3)
        let labels = Set(vs.compactMap { Self.indices($0.firstPairInMatrix).map(Self.canonical) })
        #expect(labels == [SIMD3(1, 0, 0), SIMD3(0, 1, 0), SIMD3(0, 0, 1)])
        #expect(vs.allSatisfy { if case .plane? = $0.firstPairInMatrix { true } else { false } })
        let misfit = try OrientationRelationship.secondPairMisfit(first: Self.thetaFirst, second: Self.thetaSecond,
                                                                  precipitate: theta, matrix: Self.al)
        #expect(misfit < 1e-12)
    }

    /// T1: four variants, their basal plane on the four {111} of Al.
    @Test func t1HasFourVariantsOnePerOctahedralPlane() throws {
        let t1 = PrecipitateFixtures.t1
        #expect(t1.laueClass == .sixOverMmm)
        let m = try Self.alOrientation(SIMD3(1, 1, 0))
        let vs = try OrientationRelationship.variantOrientations(first: Self.t1First, second: Self.t1Second, precipitate: t1,
                                                                 matrix: Self.al, matrixOrientation: m, volts: 200_000)
        #expect(vs.count == 4)
        for v in vs { #expect(Self.gap(v.orientation.matrix, v.relation.transpose * m.matrix) < 1e-12) }
        let labels = Set(vs.compactMap { Self.indices($0.firstPairInMatrix).map(Self.canonical) })
        #expect(labels == [SIMD3(1, 1, 1), SIMD3(1, 1, -1), SIMD3(1, -1, 1), SIMD3(1, -1, -1)])
        let misfit = try OrientationRelationship.secondPairMisfit(first: Self.t1First, second: Self.t1Second, precipitate: t1, matrix: Self.al)
        #expect(misfit < 1e-9)
    }

    /// The held three-grain case (IndistinguishableVariantsTests: β″ in Al down [001], 12° in plane, reach 0.75 Å⁻¹): with a
    /// flat Ewald sphere (10 MV) and the slab as the envelope, the grouping on drawn patterns is the held grouping exactly.
    /// At 200 kV the sphere's curvature can only split groups: every group is inside one held group.
    /// Mutation: comparing the patterns one way only merges variants whose spots are a subset of another's.
    @Test func theGroupingOnDrawnPatternsIsTheHeldGrouping() throws {
        let reach = IndistinguishableVariantsTests.reach
        let first = Self.pair(.direction(SIMD3(0, 1, 0)), .direction(SIMD3(0, 0, 1)))
        let second = Self.pair(.direction(SIMD3(0, 0, 1)), .direction(SIMD3(3, 1, 0)))
        let m = try Self.alOrientation()
        let candidate = try PhaseOrientationCandidate(name: "beta''", crystal: Self.beta, kMax: reach)
        let slab = candidate.excitationSlabInvAngstrom
        let flat = SpotPatternSimulator.Settings(sigmaExcitationError: slab / 3, excitationMultiple: 3, minimumIntensity: 0)
        let seen = VariantGroupingSettings(toleranceInvAngstrom: 0.02, relativeIntensityCut: 0, directBeamRadiusInvAngstrom: 0.15,
                                           reachInvAngstrom: reach)
        let vs = try OrientationRelationship.variantOrientations(first: first, second: second, precipitate: Self.beta, matrix: Self.al,
                                                                 matrixOrientation: m, volts: 10_000_000, kMax: reach, settings: flat,
                                                                 grouping: seen)
        let held = PhaseVectorMatcher.indistinguishableVariants(vs.map(\.relation), matrixOrientation: m.matrix,
                                                                reflections: candidate.reflections, excitationSlabInvAngstrom: slab,
                                                                reachInvAngstrom: reach)
        #expect(held.filter { $0.count == 1 }.count == 4 && held.filter { $0.count == 2 }.count == 4, "\(held)")
        let drawn = OrientationRelationship.groups(of: vs.map(\.pattern), grouping: seen)
        #expect(drawn == held)
        for v in vs { #expect(drawn[v.group].contains(v.index)) }

        let at200 = try OrientationRelationship.variantOrientations(first: first, second: second, precipitate: Self.beta, matrix: Self.al,
                                                                    matrixOrientation: m, volts: 200_000, kMax: reach, settings: flat,
                                                                    grouping: seen)
        let split = OrientationRelationship.groups(of: at200.map(\.pattern), grouping: seen)
        #expect(split.flatMap { $0 }.sorted() == Array(0 ..< 12))
        for g in split { #expect(held.contains { Set(g).isSubset(of: $0) }, "\(g) crosses the held groups \(held)") }
    }

    /// With nothing seen every variant is one group (the held rule's blind case); every variant is in exactly one group
    /// whatever the cut.
    @Test func whatCountsDecidesTheGroups() throws {
        let vs = try OrientationRelationship.variantOrientations(first: Self.betaFirst, second: Self.betaSecond, precipitate: Self.beta,
                                                                 matrix: Self.al, matrixOrientation: try Self.alOrientation(), volts: 200_000,
                                                                 kMax: 1)
        let patterns = vs.map(\.pattern)
        #expect(OrientationRelationship.groups(of: patterns, grouping: .init(reachInvAngstrom: 0)) == [Array(0 ..< 12)])
        for cut in [0, 0.01, 0.05, 1] {
            let groups = OrientationRelationship.groups(of: patterns, grouping: .init(relativeIntensityCut: cut))
            #expect(groups.flatMap { $0 }.sorted() == Array(0 ..< 12) && groups == groups.sorted { $0[0] < $1[0] })
        }
        #expect(OrientationRelationship.groups(of: []).isEmpty)
    }

    @Test func refusals() throws {
        let m = try Self.alOrientation()
        let parallel = Self.pair(.direction(SIMD3(2, 0, 0)), .direction(SIMD3(6, 4, 0)))
        #expect(throws: VariantOrientationError.invalidRelationship) {
            try OrientationRelationship.variantOrientations(first: Self.betaFirst, second: parallel, precipitate: Self.beta, matrix: Self.al,
                                                            matrixOrientation: m, volts: 200_000)
        }
        #expect(throws: VariantOrientationError.invalidRelationship) {
            try OrientationRelationship.secondPairMisfit(first: Self.pair(.direction(.zero), .direction(SIMD3(1, 0, 0))),
                                                         second: Self.betaSecond, precipitate: Self.beta, matrix: Self.al)
        }
        for bad in [VariantGroupingSettings(toleranceInvAngstrom: -1), VariantGroupingSettings(relativeIntensityCut: 2),
                    VariantGroupingSettings(directBeamRadiusInvAngstrom: .nan), VariantGroupingSettings(reachInvAngstrom: -0.1)] {
            #expect(throws: VariantOrientationError.self) {
                try OrientationRelationship.variantOrientations(first: Self.betaFirst, second: Self.betaSecond, precipitate: Self.beta,
                                                                matrix: Self.al, matrixOrientation: m, volts: 200_000, grouping: bad)
            }
        }
        #expect(throws: VariantOrientationError.spotPattern(.invalidVoltage(-1))) {
            try OrientationRelationship.variantOrientations(first: Self.betaFirst, second: Self.betaSecond, precipitate: Self.beta,
                                                            matrix: Self.al, matrixOrientation: m, volts: -1)
        }
    }
}
