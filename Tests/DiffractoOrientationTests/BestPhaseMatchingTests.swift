//
//  BestPhaseMatchingTests.swift
//  DiffractoKit#121: the best phase per position across several structures, with ASTAR's phase reliability. Planted truth: a
//  two-phase polycrystal of aluminium (fcc) and magnesium (hcp) grains, every pattern the zero-order spots of its grain's zone at
//  its in-plane angle, written out here and not taken from the libraries under test; each phase matched with its own library.
//

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct BestPhaseMatchingTests {
    static let q = 0.012, reach = 0.75
    static let origin: Float = 63.5

    struct Grain { let phase: Int; let zone: SIMD3<Double>; let inPlaneDeg: Double }

    static func reflections(_ crystal: CrystalStructure) throws -> [ACOMReflection] {
        try crystal.kinematicReflections(kMax: reach, retainingAll: false).kept
            .map { ACOMReflection(g: $0.reflection.g, intensity: $0.structureFactor.intensity) }
    }

    /// The zero-order pattern of a zone (Cartesian, the crystal's own frame) at an in-plane angle, in pixels with the direct beam.
    static func pattern(_ zone: SIMD3<Double>, inPlaneDeg: Double, _ reflections: [ACOMReflection]) -> [ACOMPeak] {
        let n = simd_normalize(zone)
        let basis = ZoneAxisFrame.detectorBasis(zoneAxis: n)
        let c = cos(inPlaneDeg * .pi / 180), s = sin(inPlaneDeg * .pi / 180)
        var peaks = [ACOMPeak(x: origin, y: origin, intensity: 1)]
        for r in reflections where abs(simd_dot(r.g, n)) < 1e-6 && simd_length(r.g) < reach {
            let x = simd_dot(r.g, basis.columns.0), y = simd_dot(r.g, basis.columns.1)
            peaks.append(ACOMPeak(x: origin + Float((c * x - s * y) / q), y: origin + Float((s * x + c * y) / q), intensity: 1))
        }
        return peaks
    }

    /// Six grains on a 6 x 4 scan, four positions each: Al on [001], [011] and [111], Mg on [0001], [2-1-10] and [01-10]
    /// (Cartesian (0, 0, 1), (1, 0, 0) and (0, 1, 0)), each at its own in-plane angle. Every position's phase is right, its beam
    /// direction lies within one template step of the planted zone (after the phase's symmetry), the runner-up is the other
    /// phase and the reliability is the published formula of the two scores.
    /// Mutations: the runner-up taken from the winning phase (reliability from the winner's own second score); the higher score
    /// losing (`<` for `>`, every phase wrong).
    @Test func aPlantedTwoPhasePolycrystalIsMappedFromItsTwoStructures() throws {
        let crystals = [CrystalPresets.aluminum, CrystalPresets.magnesium]
        let reflections = try crystals.map(Self.reflections)
        let matchers = try crystals.map { crystal in
            try #require(ACOMMatcher(library: try ACOMTemplateLibrary.generate(crystal: crystal, kMax: Self.reach)))
        }
        let grains = [
            Grain(phase: 0, zone: SIMD3(0, 0, 1), inPlaneDeg: 12), Grain(phase: 1, zone: SIMD3(0, 0, 1), inPlaneDeg: 20),
            Grain(phase: 0, zone: SIMD3(0, 1, 1), inPlaneDeg: 40), Grain(phase: 1, zone: SIMD3(1, 0, 0), inPlaneDeg: 5),
            Grain(phase: 0, zone: SIMD3(1, 1, 1), inPlaneDeg: 70), Grain(phase: 1, zone: SIMD3(0, 1, 0), inPlaneDeg: 33),
        ]
        let width = 6, height = 4
        // Grain g fills columns 2(g % 3)..+1 of rows 2(g / 3)..+1.
        func grain(_ i: Int) -> Grain { grains[((i / width) / 2) * 3 + (i % width) / 2] }
        let peaks = (0..<width * height).map { i -> [ACOMPeak] in
            let g = grain(i)
            return Self.pattern(g.zone, inPlaneDeg: g.inPlaneDeg, reflections[g.phase])
        }
        let scans = matchers.map { matcher in
            var scan = ACOMScanResults(width: width, height: height, templateCount: matcher.library.count)
            scan.results = peaks.map { matcher.match(peaks: $0, originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q) }
            return scan
        }
        let libraries = matchers.map(\.library)
        #expect(BestPhaseMatching.areComparable(libraries))
        let map = try #require(BestPhaseMatching.bestPhases(scans, libraries: libraries))
        #expect(map.phaseCounts(phases: 2) == [12, 12])
        for (i, r) in map.results.enumerated() {
            let g = grain(i)
            #expect(r.phaseIndex == g.phase, "position \(i): phase \(r.phaseIndex), planted \(g.phase)")
            #expect(r.orientation.phaseID == r.phaseIndex)
            #expect(r.runnerUpPhaseIndex == 1 - g.phase, "position \(i): runner-up \(r.runnerUpPhaseIndex)")
            #expect(r.runnerUpScore == scans[1 - g.phase].results[i].score)
            #expect(abs(r.reliability - 100 * (1 - r.runnerUpScore / r.score)) < 1e-4)
            #expect(r.reliability > 0, "position \(i): no margin over the other phase")
            // The beam direction within one template step (the library's sampling) of the planted zone, up to the phase's symmetry.
            let library = matchers[r.phaseIndex].library
            let beam = simd_normalize(library.zoneAxes[r.orientation.templateIndex])
            let planted = simd_normalize(g.zone)
            let symmetry = library.orientationSymmetry
            let angle = symmetry.operators.map { acos(min(1, abs(simd_dot($0 * planted, beam)))) * 180 / .pi }.min() ?? 180
            #expect(angle < 3, "position \(i): beam \(angle) degrees from the planted zone")
        }
        let reliabilities = map.results.map(\.reliability).sorted()
        print("BEST PHASE planted: reliability min \(reliabilities.first ?? 0), median \(reliabilities[reliabilities.count / 2])")
    }

    /// The rule on hand-made scores: the higher score wins, the runner-up is the best other phase, the published formula; a tie
    /// goes to the lower index with reliability 0; one phase alone has reliability 100; no template, a zero or a NaN score does
    /// not compete; no phase at all is `none`; scans of different sizes are refused.
    @Test func theRuleOnHandMadeScores() throws {
        func m(_ score: Float, template: Int = 0) -> OrientationResult {
            var r = OrientationResult.empty
            r.templateIndex = template; r.score = score; r.secondScore = score / 2
            return r
        }
        let three = BestPhaseMatching.best(of: [m(0.5), m(0.8), m(0.6)])
        #expect(three.phaseIndex == 1 && three.runnerUpPhaseIndex == 2 && three.runnerUpScore == 0.6)
        #expect(abs(three.reliability - 25) < 1e-4 && three.score == 0.8 && three.orientation.phaseID == 1)
        let later = BestPhaseMatching.best(of: [m(0.6), m(0.8)])
        #expect(later.phaseIndex == 1 && later.runnerUpPhaseIndex == 0, "the first leader becomes the runner-up")
        let tie = BestPhaseMatching.best(of: [m(0.7), m(0.7)])
        #expect(tie.phaseIndex == 0 && tie.runnerUpPhaseIndex == 1 && tie.reliability == 0)
        let alone = BestPhaseMatching.best(of: [m(0.4), m(0.9, template: -1), m(0), m(.nan)])
        #expect(alone.phaseIndex == 0 && alone.runnerUpPhaseIndex == -1 && alone.reliability == 100)
        #expect(BestPhaseMatching.best(of: [m(0), m(0.3, template: -1)]) == .none)
        #expect(BestPhaseMatching.best(of: []) == .none)
        let al = try ACOMTemplateLibrary.generate(crystal: CrystalPresets.aluminum, kMax: Self.reach)
        #expect(BestPhaseMatching.bestPhases([], libraries: []) == nil)
        let a = ACOMScanResults(width: 2, height: 1, templateCount: 1), b = ACOMScanResults(width: 1, height: 2, templateCount: 1)
        #expect(BestPhaseMatching.bestPhases([a, b], libraries: [al, al]) == nil, "scans of different sizes")
        #expect(BestPhaseMatching.bestPhases([a, a], libraries: [al]) == nil, "a library missing")
        #expect(BestPhaseMatching.bestPhases([a, a], libraries: [al, al]) != nil)
        // Another reach draws the pattern on another grid: the scores do not compare, refused.
        let short = try ACOMTemplateLibrary.generate(crystal: CrystalPresets.magnesium, kMax: 0.6)
        #expect(!BestPhaseMatching.areComparable([al, short]))
        #expect(BestPhaseMatching.bestPhases([a, a], libraries: [al, short]) == nil)
        #expect(BestPhaseMatching.publishedReliableAbove == 15)
    }
}
