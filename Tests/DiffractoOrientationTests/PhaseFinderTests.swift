//
//  PhaseFinderTests.swift
//  DiffractoKit S11 (#166): candidate phases ranked on one pattern. No Python reference: planted truth. Nine candidates (Al,
//  Au, Si, Fe, Mg, Mg2Si, β″, θ′, T1); ring lists of each phase with 0.3 % noise, a scale error up to 2 % and a fifth of
//  the rings missing, six per phase. The true phase is first or inside the first tie every time, and Al and Au are one tie
//  while the scale is free.
//

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoGeometry
import DiffractoOrientation
import DiffractoKitTestSupport

@Suite struct PhaseFinderTests {

    /// Mg2Si antifluorite (a 6.351 Å), planted; θ′ and T1 are PrecipitateFixtures'.
    static var mg2Si: CrystalStructure {
        let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        var sites = fcc.map { AtomSite(z: 14, fractional: $0) }
        for p in fcc { for o in [SIMD3(0.25, 0.25, 0.25), SIMD3(0.75, 0.75, 0.75)] {
            let f = p + o
            sites.append(AtomSite(z: 12, fractional: f - f.rounded(.down)))
        } }
        return CrystalPresets.structure(a: 6.351, b: 6.351, c: 6.351, sites: sites)
    }

    static let names = ["Al", "Au", "Si", "Fe", "Mg", "Mg2Si", "beta''", "theta'", "T1"]
    static let candidates: [CrystalStructure] = [CrystalPresets.aluminum, CrystalPresets.gold, CrystalPresets.silicon,
                                                 CrystalPresets.iron, CrystalPresets.magnesium, mg2Si,
                                                 CrystalPresets.betaDoublePrime, PrecipitateFixtures.thetaPrime,
                                                 PrecipitateFixtures.t1]

    struct SplitMix {
        var state: UInt64
        mutating func next() -> UInt64 {
            state &+= 0x9E37_79B9_7F4A_7C15
            var z = state
            z = (z ^ (z >> 30)) &* 0xBF58_476D_1CE4_E5B9
            z = (z ^ (z >> 27)) &* 0x94D0_49BB_1331_11EB
            return z ^ (z >> 31)
        }
        mutating func uniform() -> Double { Double(next() >> 11) / Double(1 << 53) }
        mutating func normal() -> Double {
            let u = max(uniform(), 1e-300), v = uniform()
            return (-2 * log(u)).squareRoot() * cos(2 * .pi * v)
        }
    }

    /// A planted ring list of `structure`: its rings of at least 5 % of the strongest Σ|F|² out to 1.2 Å⁻¹ (the first 12), each
    /// kept with probability 0.8 (at least 5), times 1 + a uniform scale error and 1 + 0.3 % normal noise.
    static func plantedRings(_ structure: CrystalStructure, scaleError: Double, rng: inout SplitMix) throws -> [Double] {
        let families = try structure.kinematicReflections(kMax: 1.6, retainingAll: false).families
        let strongest = families.map(\.intensity).max() ?? 0
        let visible = families.filter { $0.intensity >= 0.05 * strongest && $0.family.g <= 1.2 }.map(\.family.g).prefix(12)
        var kept = visible.filter { _ in rng.uniform() >= 0.2 }
        if kept.count < min(5, visible.count) { kept = Array(visible.prefix(5)) }
        let scale = 1 + scaleError * (2 * rng.uniform() - 1)
        return kept.map { $0 * scale * (1 + 0.003 * rng.normal()) }.sorted()
    }

    /// Mutation: a tie rule of equal counts only (no one-ring grain) splits Al and Au on a lost ring and names the wrong one.
    @Test func theTruePhaseIsFirstOrInTheFirstTie() throws {
        var rng = SplitMix(state: 20_261_124)
        var trials = 0, alone = 0
        for (t, truth) in Self.candidates.enumerated() {
            for _ in 0 ..< 6 {
                let q = try Self.plantedRings(truth, scaleError: 0.02, rng: &rng)
                let ranked = try PhaseFinder.rank(rings: q, candidates: Self.candidates)
                #expect(ranked.count == 9 && Set(ranked.map(\.candidate)) == Set(0 ..< 9))
                let first = ranked.filter { $0.tie == 0 }.map(\.candidate)
                #expect(first.contains(t), "\(Self.names[t]) not in the first tie \(first.map { Self.names[$0] }): \(ranked)")
                if t <= 1 { #expect(first.contains(0) && first.contains(1), "Al and Au split on \(Self.names[t]): \(ranked)") }
                if first == [t] { alone += 1 }
                trials += 1
            }
        }
        // Ties are the exception, not the answer: most patterns name one phase.
        #expect(trials == 54 && alone >= 27, "\(alone) of \(trials) alone")
    }

    /// Each number is the counting rule at the reported scale, recounted by brute force; agreement is their F1; the list is in
    /// the stated order (a tie together at its leader's place, the leaders by agreement).
    @Test func theNumbersAreTheRuleAtTheReportedScale() throws {
        var rng = SplitMix(state: 7)
        let q = try Self.plantedRings(CrystalPresets.silicon, scaleError: 0.02, rng: &rng)
        let ranked = try PhaseFinder.rank(rings: q, candidates: Self.candidates)
        // The ring lists reach as far as the finder's: the last ring at the widest scale and tolerance.
        let last = q.last!
        let kMax = max(last * (1 + 0.03 + 0.01) + 0.004 + 0.01, (last + max(0.01 * last, 0.004)) / (1 - 0.03) * (1 + 1e-6) + 0.01)
        for s in ranked {
            let families = try Self.candidates[s.candidate].kinematicReflections(kMax: kMax, retainingAll: false).families
            let window = { (x: Double) in max(0.01 * x, 0.004) }
            let explained = q.filter { x in families.contains { abs(x - s.scale * $0.family.g) <= window(x) } }.count
            #expect(explained == s.explained && s.measured == q.count)
            let strongest = try #require(families.map(\.intensity).max())
            let lo = q.first! - window(q.first!), hi = q.last! + window(q.last!)
            let strong = families.filter { $0.intensity >= 0.2 * strongest && s.scale * $0.family.g >= lo && s.scale * $0.family.g <= hi }
            #expect(strong.count == s.strongExpected, "\(Self.names[s.candidate])")
            let seen = strong.filter { f in q.contains { abs($0 - s.scale * f.family.g) <= window($0) } }.count
            #expect(seen == s.strongSeen && s.strongMissing == s.strongExpected - s.strongSeen)
            #expect(abs(s.scale - 1) <= 0.03 + 1e-12)
            let p = Double(s.explained) / Double(s.measured), r = s.strongExpected > 0 ? Double(s.strongSeen) / Double(s.strongExpected) : 1
            #expect(abs(s.agreement - (p + r > 0 ? 2 * p * r / (p + r) : 0)) < 1e-12)
        }
        // Ties are contiguous and numbered in order; each tie's first has the highest agreement of the candidates not
        // yet placed.
        #expect(ranked.map(\.tie) == ranked.map(\.tie).sorted())
        for (i, s) in ranked.enumerated() where i == 0 || ranked[i - 1].tie != s.tie {
            #expect(ranked[i...].allSatisfy { $0.agreement <= s.agreement })
        }
    }

    /// With the scale held (window 0) the scale is 1; a refused candidate ranks last with its reason; bad input throws.
    @Test func fixedScaleRefusalsAndInput() throws {
        let al = try CrystalPresets.aluminum.kinematicReflections(kMax: 1.2, retainingAll: false).families.map(\.family.g)
        let ranked = try PhaseFinder.rank(rings: Array(al.prefix(6)), candidates: Self.candidates,
                                          settings: .init(scaleWindow: 0))
        #expect(ranked.allSatisfy { $0.scale == 1 })
        #expect(ranked[0].candidate == 0 && ranked[0].explained == 6 && ranked[0].misfit < 1e-12)
        let unknown = CrystalPresets.structure(a: 4, b: 4, c: 4, sites: [AtomSite(z: 150, fractional: [0, 0, 0])])
        let withRefused = try PhaseFinder.rank(rings: Array(al.prefix(6)), candidates: [unknown, CrystalPresets.aluminum])
        #expect(withRefused.map(\.candidate) == [1, 0] && withRefused[1].refusal != nil && withRefused[0].refusal == nil)
        #expect(throws: PhaseFinderError.self) { try PhaseFinder.rank(rings: [], candidates: Self.candidates) }
        #expect(throws: PhaseFinderError.self) { try PhaseFinder.rank(rings: [0.5, -1], candidates: Self.candidates) }
        #expect(throws: PhaseFinderError.self) {
            try PhaseFinder.rank(rings: [0.5], candidates: Self.candidates, settings: .init(scaleWindow: 0.6))
        }
    }

    /// Spots: every number is the zone-axis solver's best candidate for that structure; an Al pattern has Al first or tied.
    @Test func spotScoresAreTheSolversNumbers() throws {
        let volts = 200_000.0
        let o = try CrystalOrientation(zoneAxis: SIMD3(0, 1, 1), inPlaneAngle: 0.3, cell: CrystalPresets.aluminum.cell)
        let pattern = try SpotPattern.simulate(CrystalPresets.aluminum, orientation: o, volts: volts, kMax: 1.2)
        let strongest = pattern.spots.map(\.intensity).max()!
        let spots = pattern.spots.filter { $0.intensity >= 1e-3 * strongest }.map { $0.position * 1.01 }
        let few = [CrystalPresets.aluminum, CrystalPresets.gold, CrystalPresets.silicon, CrystalPresets.iron]
        let ranked = try PhaseFinder.rank(spots: spots, candidates: few, volts: volts)
        for s in ranked {
            let best = try ZoneAxisSolver(structure: few[s.candidate], volts: volts).solve(spots).candidates.first
            #expect(s.best == best && s.explained == (best?.matched ?? 0) && s.agreement == (best?.agreement ?? 0))
            #expect(s.measured == spots.count)
        }
        #expect(ranked.filter { $0.tie == ranked[0].tie }.contains { $0.candidate == 0 }, "\(ranked.map { ($0.candidate, $0.agreement) })")
        #expect(throws: PhaseFinderError.self) { try PhaseFinder.rank(spots: [SIMD2(0.1, 0)], candidates: few, volts: volts) }
    }
}
