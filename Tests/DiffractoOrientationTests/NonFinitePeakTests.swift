//
//  NonFinitePeakTests.swift
//  Role: the orientation sites of DiffractoKit #117 that trapped (a crashed test process on 0.18.1): one peak that is NaN,
//        infinite or far outside the polar image reached `Int(_:)` in the polar binning inside `matchAll`'s worker pool; a
//        grid whose width times height overflows `Int`. A peak that is not finite is now left out, as
//        `PhaseVectorMatcher` leaves it out, and the rest of the pattern is matched; an overflowing grid is refused.
//

import Foundation
import Testing
import simd
import DiffractoCore
@testable import DiffractoOrientation

private func cubicReflections() -> [ACOMReflection] {
    var out: [ACOMReflection] = []
    for h in -3...3 { for k in -3...3 { for l in -3...3 where (h, k, l) != (0, 0, 0) {
        out.append(ACOMReflection(g: SIMD3(Double(h), Double(k), Double(l)) / 4, intensity: 1 / Double(h * h + k * k + l * l)))
    } } }
    return out
}

struct NonFinitePeakTests {

    static func matcher() throws -> ACOMMatcher {
        let settings = ACOMTemplateSettings(nRadial: 32, nAzimuthal: 64, zoneAxisCount: 24,
                                            excitation: .appCurvature(wavelengthAngstrom: 0.025))
        let lib = try ACOMTemplateLibrary.generate(reflections: cubicReflections(), kMax: 1.5, laueClass: .mThreeBarM,
                                                   settings: settings)
        return try #require(ACOMMatcher(library: lib))
    }

    static let finite = [ACOMPeak(x: 0.25, y: 0, intensity: 1), ACOMPeak(x: 0, y: 0.25, intensity: 0.5),
                         ACOMPeak(x: -0.25, y: 0.25, intensity: 0.3), ACOMPeak(x: 0.5, y: 0.5, intensity: 0.2)]

    @Test(arguments: [ACOMPeak(x: .nan, y: 0, intensity: 1), ACOMPeak(x: 0, y: .infinity, intensity: 1),
                      ACOMPeak(x: 0.1, y: 0.1, intensity: .nan), ACOMPeak(x: 0.1, y: 0.1, intensity: .infinity)])
    func aPeakThatIsNotFiniteIsLeftOut(_ bad: ACOMPeak) throws {
        let m = try Self.matcher()
        let clean = try #require(m.experimentalPolarImage(peaks: Self.finite, originX: 0, originY: 0))
        let withBad = try #require(m.experimentalPolarImage(peaks: Self.finite + [bad], originX: 0, originY: 0))
        #expect(withBad == clean)
        // and the scan route, through the worker pool, matches the position
        let scan = m.matchAll(peaks: [Self.finite + [bad], Self.finite], width: 2, height: 1, originX: 0, originY: 0)
        let results = try #require(scan)
        #expect(results.results[0] == results.results[1])
    }

    /// A finite radius far outside the image (1e30 Å⁻¹) deposits nothing rather than trapping in `Int(1e30 / scale)`.
    @Test func aPeakFarOutsideTheImageDepositsNothing() throws {
        let m = try Self.matcher()
        let geo = m.library.geometry
        let far = ACOMTemplateLibrary.buildPolar(spots: [(r: 1e30, azimuth: 0, weight: 1), (r: 1e300, azimuth: 0, weight: 1)],
                                                 geometry: geo, azimuthalBlurBins: 1.5, radialKernel: 0.01)
        #expect(far.allSatisfy { $0 == 0 })
        let nan = ACOMTemplateLibrary.buildPolar(spots: [(r: .nan, azimuth: 0, weight: 1), (r: 0.3, azimuth: .nan, weight: 1),
                                                         (r: 0.3, azimuth: 0, weight: .infinity)],
                                                 geometry: geo, azimuthalBlurBins: 1.5, radialKernel: 0)
        #expect(nan.allSatisfy { $0 == 0 })
    }

    @Test func aGridThatOverflowsIsRefused() throws {
        #expect(!gridHolds(0, width: Int.max, height: 2))
        #expect(gridHolds(6, width: 3, height: 2))
        #expect(OrientationMapAnalysis.eulerText(results: [], width: Int.max, height: 2, x: 0, y: 0) == nil)
        let m = try Self.matcher()
        #expect(m.matchAll(peaks: [], width: 1 << 62, height: 4, originX: 0, originY: 0) == nil)
    }
}
