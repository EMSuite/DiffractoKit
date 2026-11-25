import Foundation
import Testing
@testable import DiffractoGeometry

/// The Au {111} ring of pyxem's public cross-grating pattern (`au_xgrating_20cm.tif`, annulus 37 to 46 px) read 0.10073 nm⁻¹/px in the
/// app where py4DSTEM's `fit_ellipse_amorphous_ring` (f050d207) gives 0.10053 (EMSuite/Electron-Microscopy-Studio#1002, DiffractoKit#242).
/// CAUSE: py4DSTEM's `leastsq` leaves the central Gaussian's amplitude I0 free, and on that ring it settles at I0 = -1.18 (sigma0 13 px):
/// a background rising outward across the annulus. The kit held I0 >= 0, so the ring's asymmetric Gaussian took up the slope and the
/// radius came out 0.07 px small (42.18 px against 42.25). What would have refuted it: the bounded and the free fit agreeing on that
/// ring, or the start centre moving py4DSTEM's answer (it does not: 0.10052 to 0.10054 from the image centre, the ring's centre and the
/// pattern's centre).
///
/// The planted pattern is py4DSTEM's fitted profile of that ring, rounded: I0 -1.2, sigma0 13, I1 0.057, sigma1 = sigma2 1.8, background
/// 0.083, centre (127.84, 127.38), a 43, b 41.5, theta 1 rad, with a flat central beam of 60 inside 20 px (as the real pattern has one; it
/// lies outside the annulus and only sets the intensity bound). The pinned numbers are py4DSTEM f050d207's `fit_ellipse_amorphous_ring`
/// on this pattern from the start centre (128, 128); the kit held at I0 >= 0 gives R 42.178 px (a NumPy transcription of the fit).
/// Tolerance: py4DSTEM ships no test of this function, so `assert_allclose`'s default for float64, rtol 1e-7 (docs/parity.md); the
/// transcription of the kit's fit lands 1e-10 px from py4DSTEM's R.
struct RingFitCentralGaussianSignTests {
    static let planted = (i0: -1.2, i1: 0.057, sigma0: 13.0, sigma1: 1.8, sigma2: 1.8, background: 0.083,
                          row: 127.84, column: 127.38, a: 43.0, b: 41.5, theta: 1.0)
    // py4DSTEM f050d207, fit_ellipse_amorphous_ring(pattern, (128, 128), (37, 46)): x0, y0, a, b, (a + b) / 2.
    static let reference = (row: 127.8399999963596, column: 127.37999999485223, a: 42.99999996995962,
                            b: 41.49999997813876, radius: 42.24999997404919)

    /// py4DSTEM's `double_sided_gaussian` in (a, b, theta), as `fitAmorphousRing` evaluates it; x is the row, y the column.
    static func pattern(n: Int = 256) -> DiffractionPattern {
        let p = planted
        let s = sin(p.theta), c = cos(p.theta), a2 = p.a * p.a, b2 = p.b * p.b
        let radius = (p.a + p.b) / 2, r2 = radius * radius
        let A = (s * s / b2 + c * c / a2) * r2, B = 2 * (b2 - a2) * s * c / (a2 * b2) * r2, C = (c * c / b2 + s * s / a2) * r2
        var pixels = [Float](repeating: 0, count: n * n)
        for row in 0 ..< n { for column in 0 ..< n {
            let dx = Double(row) - p.row, dy = Double(column) - p.column
            let e2 = max(0, A * dx * dx + B * dx * dy + C * dy * dy)
            let delta = e2.squareRoot() - radius
            let sigma = delta < 0 ? p.sigma1 : p.sigma2
            if hypot(dx, dy) < 20 { pixels[row * n + column] = 60; continue }
            let value = p.i0 * exp(-e2 / (2 * p.sigma0 * p.sigma0)) + p.i1 * exp(-delta * delta / (2 * sigma * sigma)) + p.background
            pixels[row * n + column] = Float(value)
        } }
        return DiffractionPattern(width: n, height: n, pixels: pixels)
    }

    static func close(_ got: Double, _ want: Double, rtol: Double = 1e-7) -> Bool { abs(got - want) <= rtol * abs(want) }

    @Test func aBackgroundRisingOutwardIsANegativeCentralGaussianAsInPy4DSTEM() throws {
        let ring = try RingCalibration.fitRing(pattern: Self.pattern(), centerQX: 128, centerQY: 128, innerRadius: 37, outerRadius: 46)
        let fit = try #require(ring.fit)
        let radius = try #require(ring.radius)
        #expect(fit.model == .asymmetricGaussian)
        #expect(ring.isUsable, Comment(rawValue: "refusals \(ring.refusals)"))
        #expect(Self.close(radius, Self.reference.radius), Comment(rawValue: "R \(radius), py4DSTEM \(Self.reference.radius)"))
        #expect(Self.close(max(fit.a, fit.b), Self.reference.a), Comment(rawValue: "a \(max(fit.a, fit.b))"))
        #expect(Self.close(min(fit.a, fit.b), Self.reference.b), Comment(rawValue: "b \(min(fit.a, fit.b))"))
        #expect(Self.close(fit.centerQX, Self.reference.row), Comment(rawValue: "row \(fit.centerQX)"))
        #expect(Self.close(fit.centerQY, Self.reference.column), Comment(rawValue: "column \(fit.centerQY)"))
        let central = try #require(fit.profile?.centralIntensity)
        #expect(central < 0, Comment(rawValue: "I0 \(central)"))
    }
}
