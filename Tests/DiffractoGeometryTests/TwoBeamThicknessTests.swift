//
//  TwoBeamThicknessTests.swift
//  Role: the thickness from a two-beam CBED disc profile (TwoBeamThickness) against the closed form (Kelly et al. 1975,
//        Allen 1981): tools/pins/two-beam-thickness/pins.py writes the intensity, the minima and the fringe line in NumPy.
//        No Python package does this measurement, so the closed form is the reference: the intensity and the fit on
//        noiseless profiles (40 to 400 nm, extinction distance 50 to 150 nm) at rtol 1e-6, the fringe line through exact
//        minima at rtol 1e-7 (assert_allclose's float64 default). The simulated disc is in the parity target
//        (TwoBeamThicknessDiscTests).
//

import DiffractoKitTestSupport
import Foundation
import Testing
import DiffractoCore
import DiffractoGeometry

struct TwoBeamThicknessTests {
    struct Case: Decodable {
        let thickness: Double
        let extinctionDistance: Double
        let deviation: [Double]
        let intensity: [Double]
        let orders: [Int]
        let minima: [Double]
        let profileHalfWidth: Double
        let profileSamples: Int
        let profileCentre: Double
        let fringeLineThickness: Double
        let fringeLineExtinctionDistance: Double
    }
    struct Pins: Decodable {
        let cases: [Case]
        let siExtinctionDistance: Double
        let siPendelloesungPeriod: Double
    }

    static let pins: Pins = {
        let url = FixtureDirs.fixtures(of: "DiffractoGeometryTests").appendingPathComponent("two-beam-thickness-pins.json")
        let decoder = JSONDecoder()
        decoder.keyDecodingStrategy = .convertFromSnakeCase
        // swiftlint:disable:next force_try
        return try! decoder.decode(Pins.self, from: Data(contentsOf: url))
    }()

    static func close(_ a: Double, _ b: Double, rtol: Double) -> Bool { abs(a - b) <= rtol * abs(b) }

    /// The noiseless profile of a case: scale 2.5, background 0.1, the Bragg position off the sampling's centre.
    static func profile(_ c: Case) -> TwoBeamProfile {
        var s: [Double] = []
        var y: [Double] = []
        for i in 0..<c.profileSamples {
            let x = -c.profileHalfWidth + 2 * c.profileHalfWidth * Double(i) / Double(c.profileSamples - 1)
            s.append(x)
            let v = TwoBeamThickness.intensity(deviation: x - c.profileCentre, thickness: c.thickness,
                                               extinctionDistance: c.extinctionDistance)
            y.append(2.5 * v + 0.1)
        }
        return TwoBeamProfile(deviations: s, intensities: y)
    }

    @Test func theFixtureCoversTheRange() {
        #expect(Self.pins.cases.count == 7)
        #expect(Self.pins.cases.map(\.thickness).min() == 400 && Self.pins.cases.map(\.thickness).max() == 4000)
        #expect(Self.pins.cases.map(\.extinctionDistance).min() == 500 && Self.pins.cases.map(\.extinctionDistance).max() == 1500)
    }

    @Test func theIntensityIsTheClosedForm() {
        for c in Self.pins.cases {
            for k in 0..<c.deviation.count {
                let got = TwoBeamThickness.intensity(deviation: c.deviation[k], thickness: c.thickness,
                                                     extinctionDistance: c.extinctionDistance)
                #expect(Self.close(got, c.intensity[k], rtol: 1e-6), "t \(c.thickness), s \(c.deviation[k]): \(got)")
            }
        }
    }

    @Test func theMinimaAreTheClosedForm() throws {
        for c in Self.pins.cases {
            for k in 0..<c.orders.count {
                let got = try #require(TwoBeamThickness.minimumDeviation(order: c.orders[k], thickness: c.thickness,
                                                                         extinctionDistance: c.extinctionDistance))
                #expect(Self.close(got, c.minima[k], rtol: 1e-7))
                // the intensity is zero there
                let at = TwoBeamThickness.intensity(deviation: got, thickness: c.thickness, extinctionDistance: c.extinctionDistance)
                #expect(at < 1e-20)
            }
            // no minimum of order t/ξ or below
            #expect(TwoBeamThickness.minimumDeviation(order: c.orders[0] - 1, thickness: c.thickness,
                                                      extinctionDistance: c.extinctionDistance) == nil)
        }
    }

    /// Allen's line through exact minima gives back the thickness and the extinction distance (rtol 1e-7), and agrees
    /// with numpy.polyfit through the same points.
    @Test func theFringeLineThroughExactMinima() throws {
        for c in Self.pins.cases {
            let line = try TwoBeamThickness.fringeLine(orders: c.orders, deviations: c.minima)
            #expect(Self.close(line.thickness, c.thickness, rtol: 1e-7), "\(line.thickness) vs \(c.thickness)")
            let xi = try #require(line.extinctionDistance)
            #expect(Self.close(xi, c.extinctionDistance, rtol: 1e-7), "\(xi) vs \(c.extinctionDistance)")
            #expect(Self.close(line.thickness, c.fringeLineThickness, rtol: 1e-7))
            #expect(Self.close(xi, c.fringeLineExtinctionDistance, rtol: 1e-7))
        }
        #expect(throws: TwoBeamThicknessError.invalidFringeMinima) {
            try TwoBeamThickness.fringeLine(orders: [3], deviations: [0.001])
        }
        #expect(throws: TwoBeamThicknessError.invalidFringeMinima) {
            try TwoBeamThickness.fringeLine(orders: [0, 1], deviations: [0.001, 0.002])
        }
    }

    /// The whole-profile fit on noiseless profiles gives back the thickness and the extinction distance within rtol 1e-6,
    /// with the two-beam value anywhere that keeps the true one inside the window, and finds no second candidate.
    @Test(arguments: [0.75, 1.1])
    func theFitGivesBackANoiselessProfile(offset: Double) throws {
        for c in Self.pins.cases {
            let fit = try TwoBeamThickness.fit(Self.profile(c), extinctionDistance: c.extinctionDistance / offset)
            let label = "t \(c.thickness), xi \(c.extinctionDistance), window at \(offset)"
            #expect(Self.close(fit.thickness, c.thickness, rtol: 1e-6), "\(label): t \(fit.thickness)")
            #expect(Self.close(fit.extinctionDistance, c.extinctionDistance, rtol: 1e-6), "\(label): xi \(fit.extinctionDistance)")
            #expect(abs(fit.centre - c.profileCentre) < 1e-7 * c.profileHalfWidth, "\(label): centre \(fit.centre)")
            #expect(Self.close(fit.scale, 2.5, rtol: 1e-6) && abs(fit.background - 0.1) < 1e-6, "\(label)")
            #expect(fit.residual < 1e-6, "\(label): residual \(fit.residual)")
            #expect(fit.secondCandidate == nil && !fit.isAmbiguous, "\(label): \(String(describing: fit.secondCandidate))")
            #expect(fit.minimaBelow >= 4 && fit.minimaAbove >= 4, "\(label)")
            #expect(fit.twoBeamExtinctionDistance == c.extinctionDistance / offset)
            #expect(!fit.extinctionDistanceAtWindowEdge)
            // the fringe line through the sampled minima (parabolas through three samples): close, not exact
            let line = try #require(fit.fringeLine, "\(label)")
            #expect(line.orders.count >= 3 && abs(line.thickness / c.thickness - 1) < 1e-3, "\(label): line \(line.thickness)")
        }
    }

    /// The calibration-free deviation: s = (x / D) λ g².
    @Test func theDeviationFromADiscLine() {
        let lambda = ElectronWavelength.angstroms(volts: 200e3)
        let g = 0.5208
        let p = TwoBeamProfile(offsets: [-10, 0, 25], discSeparation: 100, gLength: g, wavelength: lambda, intensities: [1, 2, 3])
        let factor = lambda * g * g / 100
        #expect(p.deviations == [-10 * factor, 0, 25 * factor])
        #expect(p.intensities == [1, 2, 3])
    }

    /// Si 220 at 200 kV: the two-beam extinction distance from the kit's structure factor is the NumPy pin's (rtol 1e-7)
    /// and within 5 % of the Pendellösung period abTEM 1.0.10's multislice shows at the Bragg position (95.6 nm).
    @Test func theExtinctionDistanceOfSilicon220() throws {
        let si = CrystalPresets.diamond(a: 5.43, z: 14)
        let xi = try TwoBeamThickness.extinctionDistance(of: MillerIndex(2, 2, 0), in: si, volts: 200e3)
        #expect(Self.close(xi, Self.pins.siExtinctionDistance, rtol: 1e-7), "\(xi)")
        #expect(abs(xi / Self.pins.siPendelloesungPeriod - 1) < 0.05, "\(xi)")
        // a forbidden reflection has no extinction distance
        #expect(throws: TwoBeamThicknessError.noStructureFactor(MillerIndex(2, 0, 0))) {
            try TwoBeamThickness.extinctionDistance(of: MillerIndex(2, 0, 0), in: si, volts: 200e3)
        }
        #expect(throws: TwoBeamThicknessError.invalidVoltage(0)) {
            try TwoBeamThickness.extinctionDistance(of: MillerIndex(2, 2, 0), in: si, volts: 0)
        }
    }

    /// A profile with fewer than two minima per side is refused by name, never answered.
    @Test func tooFewMinimaIsRefused() throws {
        // 40 nm, extinction distance 95 nm, out to the first minimum and a little past it
        let t = 400.0
        let xi = 950.0
        let s1 = try #require(TwoBeamThickness.minimumDeviation(order: 1, thickness: t, extinctionDistance: xi))
        let s2 = try #require(TwoBeamThickness.minimumDeviation(order: 2, thickness: t, extinctionDistance: xi))
        let half = 0.5 * (s1 + s2)
        var s: [Double] = []
        var y: [Double] = []
        for i in 0..<201 {
            let x = -half + 2 * half * Double(i) / 200
            s.append(x)
            y.append(TwoBeamThickness.intensity(deviation: x, thickness: t, extinctionDistance: xi))
        }
        #expect(throws: TwoBeamThicknessError.tooFewMinima(below: 1, above: 1)) {
            try TwoBeamThickness.fit(TwoBeamProfile(deviations: s, intensities: y), extinctionDistance: xi)
        }
    }

    @Test func badInputIsRefusedByName() {
        let s = (0..<20).map { Double($0) * 1e-4 }
        let y = s.map { 1 + sin(3000 * $0) }
        #expect(throws: TwoBeamThicknessError.mismatchedLengths(deviations: 20, intensities: 19)) {
            try TwoBeamThickness.fit(TwoBeamProfile(deviations: s, intensities: Array(y.dropLast())), extinctionDistance: 900)
        }
        #expect(throws: TwoBeamThicknessError.tooFewSamples(10)) {
            try TwoBeamThickness.fit(TwoBeamProfile(deviations: Array(s.prefix(10)), intensities: Array(y.prefix(10))),
                                     extinctionDistance: 900)
        }
        var bad = y
        bad[5] = .nan
        #expect(throws: TwoBeamThicknessError.nonFiniteSample(5)) {
            try TwoBeamThickness.fit(TwoBeamProfile(deviations: s, intensities: bad), extinctionDistance: 900)
        }
        var back = s
        back[7] = back[6]
        #expect(throws: TwoBeamThicknessError.deviationsNotIncreasing(7)) {
            try TwoBeamThickness.fit(TwoBeamProfile(deviations: back, intensities: y), extinctionDistance: 900)
        }
        #expect(throws: TwoBeamThicknessError.flatProfile) {
            try TwoBeamThickness.fit(TwoBeamProfile(deviations: s, intensities: [Double](repeating: 2, count: 20)),
                                     extinctionDistance: 900)
        }
        #expect(throws: TwoBeamThicknessError.invalidExtinctionDistance(-1)) {
            try TwoBeamThickness.fit(TwoBeamProfile(deviations: s, intensities: y), extinctionDistance: -1)
        }
        #expect(throws: TwoBeamThicknessError.invalidThicknessRange) {
            try TwoBeamThickness.fit(TwoBeamProfile(deviations: s, intensities: y), extinctionDistance: 900,
                                     thicknessRange: 500...500)
        }
        #expect(TwoBeamThicknessError.tooFewMinima(below: 1, above: 0).localizedDescription.contains("at least two on each side"))
    }
}
