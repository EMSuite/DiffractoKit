//
//  TwoBeamThicknessDiscTests.swift
//  Role: the thickness from a simulated two-beam CBED disc (docs/parity.md P56), read from
//        Tests/DiffractoGeometryTests/Fixtures/two-beam-thickness-disc.json (tools/pins/two-beam-thickness/disc.py: Si 220 at
//        200 kV, abTEM 1.0.10 potential slices, a plane-wave multislice per point of the disc along g, static lattice).
//        No Python package measures thickness from a disc, so the pin is the simulation's known thickness.
//  Done when (plan row MD2): 108 mrad off [001] (a clean two-beam orientation) the thickness is within 3 % at 60, 100, 150,
//  200, 250 and 300 nm with no second candidate, the extinction distance fitted from the kit's Si 220 value; 42 mrad off
//  [001] (a strong line through the disc) at 250 nm two thicknesses fit equally and both are returned. The NumPy twin of
//  the fit (tools/pins/two-beam-thickness/fit.py) gives +0.20 to +2.24 % on the clean cases and the fringe line +0.07 to
//  +1.79 %; at 42 mrad, 250 nm it gives 388.8 and 349.6 nm (the true 249.8 nm fits 27 % worse: the two-beam picture fails
//  there, which is why the answer is both and never one).
//

import DiffractoKitTestSupport
import Foundation
import Testing
import DiffractoCore
import DiffractoGeometry

private struct DiscCase: Decodable {
    let name: String
    let orientation: String
    let thickness: Double
    let deviation: [Double]
    let intensity: [Double]
}

private struct DiscFixture: Decodable {
    let volts: Double
    let g: Double
    let wavelength: Double
    let cases: [DiscCase]
}

struct TwoBeamThicknessDiscTests {
    private static let fixture: DiscFixture = {
        let url = FixtureDirs.fixtures(of: "DiffractoGeometryTests").appendingPathComponent("two-beam-thickness-disc.json")
        // swiftlint:disable:next force_try
        return try! JSONDecoder().decode(DiscFixture.self, from: Data(contentsOf: url))
    }()

    private static let silicon = CrystalPresets.diamond(a: 5.43, z: 14)

    @Test func theFixtureIsTheSimulation() {
        let f = Self.fixture
        #expect(f.cases.count == 7 && f.volts == 200e3)
        #expect(f.cases.allSatisfy({ $0.deviation.count == 261 && $0.intensity.count == 261 }))
        #expect(abs(f.wavelength / ElectronWavelength.angstroms(volts: 200e3) - 1) < 1e-7)
        #expect(abs(f.g * 5.43 / 8.0.squareRoot() - 1) < 1e-12)
    }

    @Test func aCleanDiscGivesTheThicknessWithin3Percent() throws {
        let clean = Self.fixture.cases.filter { $0.orientation == "clean" }
        #expect(clean.count == 6)
        for c in clean {
            let profile = TwoBeamProfile(deviations: c.deviation, intensities: c.intensity)
            let fit = try TwoBeamThickness.fit(profile, reflection: MillerIndex(2, 2, 0), in: Self.silicon, volts: 200e3)
            let error = fit.thickness / c.thickness - 1
            #expect(abs(error) < 0.03, "\(c.name): \(fit.thickness) Å against \(c.thickness) Å (\(100 * error) %)")
            #expect(!fit.isAmbiguous, "\(c.name): second candidate \(String(describing: fit.secondCandidate))")
            #expect(fit.minimaBelow >= 2 && fit.minimaAbove >= 2)
            let line = try #require(fit.fringeLine, "\(c.name): no fringe line")
            #expect(abs(line.thickness / c.thickness - 1) < 0.03, "\(c.name): fringe line \(line.thickness) Å")
        }
    }

    @Test func aStrongLineThroughTheDiscReturnsTwoThicknesses() throws {
        let c = try #require(Self.fixture.cases.first(where: { $0.orientation == "line" }))
        let profile = TwoBeamProfile(deviations: c.deviation, intensities: c.intensity)
        let fit = try TwoBeamThickness.fit(profile, reflection: MillerIndex(2, 2, 0), in: Self.silicon, volts: 200e3)
        let second = try #require(fit.secondCandidate, "\(c.name): one answer, \(fit.thickness) Å")
        #expect(fit.isAmbiguous)
        #expect(abs(second.thickness / fit.thickness - 1) > 0.05)
        #expect(second.residual <= 1.05 * fit.residual)
    }

    /// The same disc measured in pixels: offsets from the diffracted disc's centre and the disc separation give the same
    /// deviations, so the same thickness, with no calibration.
    @Test func theDiscLineInPixelsGivesTheSameThickness() throws {
        let f = Self.fixture
        let c = try #require(f.cases.first(where: { $0.name == "clean-150nm" }))
        let separation = 412.0  // pixels between the disc centres: 2θ_B = λ g
        let perPixel = f.wavelength * f.g / separation  // radians per pixel
        let offsets = c.deviation.map { $0 / f.g / perPixel }
        let profile = TwoBeamProfile(offsets: offsets, discSeparation: separation, gLength: f.g, wavelength: f.wavelength,
                                     intensities: c.intensity)
        for k in 0..<offsets.count { #expect(abs(profile.deviations[k] - c.deviation[k]) <= 1e-12 * abs(c.deviation[k]) + 1e-18) }
        let fit = try TwoBeamThickness.fit(profile, reflection: MillerIndex(2, 2, 0), in: Self.silicon, volts: f.volts)
        #expect(abs(fit.thickness / c.thickness - 1) < 0.03)
    }
}
