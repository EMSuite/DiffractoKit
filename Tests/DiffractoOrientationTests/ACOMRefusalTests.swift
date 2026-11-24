//
//  ACOMRefusalTests.swift
//  Role: DiffractoKit#123. Values the ACOM route must refuse rather than carry: a NaN excitation error at the zero of
//        py4DSTEM's denominator 2 − 2λ g·n, a measured peak that is not finite, and an orientation the symmetry reduction
//        refuses. Each test fails on the code before the fix (a trap shows as a crashed test process).
//

import Testing
import simd
import DiffractoCore
@testable import DiffractoOrientation

@Suite struct ACOMRefusalTests {

    /// A simple cubic set of reflections (a = 4 Å, indices to 2), every intensity 1: enough to build a small m-3m library.
    private static func library() throws -> ACOMTemplateLibrary {
        var reflections: [ACOMReflection] = []
        for h in -2...2 { for k in -2...2 { for l in -2...2 where h != 0 || k != 0 || l != 0 {
            reflections.append(ACOMReflection(g: SIMD3(Double(h), Double(k), Double(l)) / 4, intensity: 1))
        } } }
        return try ACOMTemplateLibrary.generate(reflections: reflections, kMax: 1.0, laueClass: .mThreeBarM,
                                                settings: ACOMTemplateSettings(zoneAxisCount: 30))
    }

    private static let peaks: [ACOMPeak] = [
        ACOMPeak(x: 110, y: 100, intensity: 1), ACOMPeak(x: 90, y: 100, intensity: 1),
        ACOMPeak(x: 100, y: 110, intensity: 1), ACOMPeak(x: 100, y: 90, intensity: 1),
        ACOMPeak(x: 110, y: 110, intensity: 0.5), ACOMPeak(x: 90, y: 90, intensity: 0.5),
    ]

    /// py4DSTEM's s = (2 g·n − λ|g|²) / (2 − 2λ g·n) is 0/0 for g·n = 1/λ and |g|² = 2/λ²: at λ = 1 Å, n = z and g = (1, 0, 1)
    /// both terms are exactly 0. py4DSTEM's `np.abs(sg) < kernel` drops that reflection; the kit's `abs(sg) > sgMax` kept it
    /// with a NaN weight. The ordinary reflection beside it is kept.
    @Test func aNaNExcitationErrorIsNotProjected() {
        let settings = ACOMTemplateSettings(excitation: .py4DSTEM(wavelengthAngstrom: 1))
        let reflections = [ACOMReflection(g: SIMD3(1, 0, 1), intensity: 1), ACOMReflection(g: SIMD3(0.2, 0, 0), intensity: 1)]
        let spots = ACOMTemplateLibrary.project(reflections: reflections, zoneAxis: SIMD3(0, 0, 1), settings: settings)
        #expect(spots.count == 1)
        #expect(spots.allSatisfy { $0.weight.isFinite && $0.r.isFinite && $0.azimuth.isFinite })
        #expect(abs((spots.first?.r ?? 0) - 0.2) <= 1e-15)
    }

    /// A NaN origin or pixel size gave every peak a NaN radius. A non-finite peak alone is left out (NonFinitePeakTests).
    @Test func aNonFiniteOriginOrScaleIsRefused() throws {
        let matcher = try #require(ACOMMatcher(library: try Self.library()))
        #expect(matcher.match(peaks: Self.peaks, originX: .nan, originY: 100, invAngstromPerPixel: 0.025) == .empty)
        #expect(matcher.match(peaks: Self.peaks, originX: 100, originY: .infinity, invAngstromPerPixel: 0.025) == .empty)
        #expect(matcher.match(peaks: Self.peaks, originX: 100, originY: 100, invAngstromPerPixel: .nan) == .empty)
        #expect(matcher.match(peaks: Self.peaks, originX: 100, originY: 100, invAngstromPerPixel: 0) == .empty)
        #expect(matcher.experimentalPolarImage(peaks: Self.peaks, originX: 100, originY: 100, invAngstromPerPixel: -1) == nil)
        #expect(matcher.templateScores(peaks: Self.peaks, originX: .nan, originY: 100, invAngstromPerPixel: 0.025).isEmpty)
    }

    /// A finite peak far off the detector gave a radius past Int.max, and `Int(_:)` trapped in the polar image. It lies
    /// outside every ring, so it is left out and the rest of the pattern is matched as without it.
    @Test func aPeakFarOffTheDetectorIsLeftOut() throws {
        let matcher = try #require(ACOMMatcher(library: try Self.library()))
        let far = Self.peaks + [ACOMPeak(x: 1e30, y: 100, intensity: 1)]
        let result = matcher.match(peaks: far, originX: 100, originY: 100, invAngstromPerPixel: 0.025)
        #expect(result == matcher.match(peaks: Self.peaks, originX: 100, originY: 100, invAngstromPerPixel: 0.025))
        #expect(result.templateIndex >= 0)
    }

    /// When the symmetry reduction refuses the matrix (here a non-finite angle), the result is `.empty`; before, it was the
    /// unreduced matrix with disorientation 0.
    @Test func anOrientationTheReductionRefusesIsEmpty() throws {
        let matcher = try #require(ACOMMatcher(library: try Self.library()))
        let refused = matcher.makeResult(bestTemplate: 0, bestBin: 0, bestScore: 0.9, secondScore: 0.1, mirrored: false,
                                         binOffset: .nan)
        #expect(refused == .empty)
        let kept = matcher.makeResult(bestTemplate: 0, bestBin: 3, bestScore: 0.9, secondScore: 0.1, mirrored: false)
        #expect(kept.templateIndex == 0)
        #expect(kept.euler.phi1.isFinite && kept.euler.Phi.isFinite && kept.euler.phi2.isFinite)
    }
}
