//
//  NonFiniteInputTests.swift
//  Role: the geometry sites of DiffractoKit #117 that trapped (shown as a crashed test process on 0.18.1): a pattern
//        whose width times height overflows `Int`, a profile bin width so small its bin count is infinite or gigabytes, and
//        a fit with extreme but finite axes whose coefficients overflow to infinity before `Int(_:)`. Each now answers nil
//        (or an empty measurement), as the functions document for other unusable input.
//

import Foundation
import Testing
@testable import DiffractoGeometry

struct NonFiniteInputTests {

    static let huge = DiffractionPattern(width: 1 << 62, height: 8, pixels: [Float](repeating: 1, count: 16))

    @Test func aPatternWhoseSizeOverflowsIsRefused() {
        #expect(PatternCentre.find(pattern: Self.huge) == nil)
        #expect(PatternRadialProfile.about(pattern: Self.huge, centerRow: 1, centerColumn: 1) == nil)
        #expect(PatternSpot.refine(pattern: Self.huge, row: 1, column: 1) == nil)
        #expect(FriedelCentre.find(pattern: Self.huge) == nil)
        #expect(FriedelCentre.find(pattern: Self.huge, mask: [Bool](repeating: true, count: 16)) == nil)
    }

    @Test(arguments: [1e-300, 1e-12])
    func aBinWidthThatNeedsTooManyBinsIsRefused(_ width: Double) {
        let p = PatternMeasurementTests.pattern(n: 32) { _, _ in 1 }
        #expect(PatternRadialProfile.about(pattern: p, centerRow: 16, centerColumn: 16, binWidth: width) == nil)
        let profile = PatternRadialProfile.about(pattern: p, centerRow: 16, centerColumn: 16, binWidth: 0.5)
        #expect(profile?.mean.count == 30)
    }

    static func extremeFit(a: Double, b: Double) -> EllipseCalibrationFit {
        EllipseCalibrationFit(centerQX: 32, centerQY: 32, a: a, b: b, theta: 0.3, normalizedResidual: 0, conicResidual: 0,
                              sampleCount: 100, occupiedAngularBins: 36, model: .conic, profile: nil, profileFallbackReason: nil)
    }

    @Test func extremeAxesOverflowNoBin() {
        let fit = Self.extremeFit(a: 1e-200, b: 1e200)
        let xs = (0..<100).map { Double($0 % 10) + 27 }, ys = (0..<100).map { Double($0 / 10) + 27 }
        #expect(EllipseCalibration.halfMaximumRingSigma(values: [Double](repeating: 1, count: 100), xs: xs, ys: ys, seed: fit,
                                                        innerRadius: 10, outerRadius: 30) == nil)
        #expect(EllipseCalibration.halfMaximumRingSigma(values: [], xs: [], ys: [], seed: Self.extremeFit(a: 20, b: 20),
                                                        innerRadius: 10, outerRadius: .infinity) == nil)
        let ring = PatternMeasurementTests.pattern(n: 64) { y, x in
            let r = ((y - 32) * (y - 32) + (x - 32) * (x - 32)).squareRoot()
            return 10 + 100 * exp(-(r - 20) * (r - 20) / 4)
        }
        let evidence = EllipseCalibration.ringEvidence(pattern: ring, fit: fit, centerQX: 32, centerQY: 32,
                                                       innerRadius: 10, outerRadius: 30)
        #expect(evidence.occupiedBins >= 0)   // measured without a trap
    }
}
