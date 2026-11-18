//
//  EllipseCrystallineRingTests.swift
//  `fitBestAvailable` on a crystalline ring standard: the scan-mean pattern of py4DSTEM's calibrationData_bullseyeProbe.h5
//  (polycrystalline Au, 8,400 patterns; Fixtures/bullseye-mean.f32, 250 x 250 Float32, row-major little-endian, sha256
//  c7cc328d...), centre and annulus as 4DSTEM Studio passes them after its centre-of-mass origin (4DSTEM-Studio#460).
//

import Foundation
import Testing
import DiffractoGeometry

@Suite struct EllipseCrystallineRingTests {
    static let pattern: DiffractionPattern = {
        let url = URL(fileURLWithPath: #filePath).deletingLastPathComponent().appendingPathComponent("Fixtures/bullseye-mean.f32")
        let data = try! Data(contentsOf: url)
        let pixels = data.withUnsafeBytes { raw in
            stride(from: 0, to: raw.count, by: 4).map {
                Float(bitPattern: UInt32(littleEndian: raw.loadUnaligned(fromByteOffset: $0, as: UInt32.self)))
            }
        }
        return DiffractionPattern(width: 250, height: 250, pixels: pixels)
    }()

    /// The amorphous-ring profile settles on b/a 0.68 (38-50 px) and 0.47 (40-48 px) on these sharp Au rings; py4DSTEM's
    /// `fit_ellipse_1D` (0.14.19) gives a 44.77 / b 44.55 and a 44.39 / b 44.17. The fit keeps the 1D ring fit and says why.
    /// Not a parity pin (py4DSTEM's 1D fit ran on a float64 mean): within 0.05 px. Mutation: drop the departure check ->
    /// a 50.18 / b 34.21 (asymmetric Gaussian), red.
    @Test(arguments: [(38.0, 50.0, 44.77, 44.55), (40.0, 48.0, 44.39, 44.17)])
    func aCrystallineMeanPatternKeepsTheRingFit(inner: Double, outer: Double, a: Double, b: Double) throws {
        #expect(Self.pattern.pixels.count == 250 * 250)
        let fit = try EllipseCalibration.fitBestAvailable(pattern: Self.pattern, centerQX: 127.45553588867188,
                                                          centerQY: 127.04931640625, innerRadius: inner, outerRadius: outer)
        #expect(fit.model == .conic)
        #expect(abs(fit.a - a) < 0.05 && abs(fit.b - b) < 0.05, "a \(fit.a) b \(fit.b)")
        #expect(fit.profileFallbackReason?.contains("departs from the 1D ring fit") == true)
    }
}
