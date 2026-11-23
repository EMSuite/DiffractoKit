//
//  SpotPatternTests.swift
//  Role: `SpotPattern` (S2) on closed forms: aluminium down [001] (the zero layer's spots, their excitation error −λg²/2,
//        the envelope, the forbidden reflections and the direct beam absent), an in-plane turn moves every spot, R = L tan 2θ
//        and the detector position, the simulator reused, and the refusals. The py4DSTEM, diffsims and 60-digit pins are in
//        DiffractoKitParityTests/Geometry/SpotPatternPinsTests.swift.
//

import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoGeometry
import DiffractoScattering

@Suite struct SpotPatternTests {

    private let aluminium: CrystalStructure = {
        let fcc = [SIMD3(0.0, 0, 0), SIMD3(0, 0.5, 0.5), SIMD3(0.5, 0, 0.5), SIMD3(0.5, 0.5, 0)]
        return CrystalStructure(cell: .cubic(a: 4.0495), sites: fcc.map { AtomSite(z: 13, fractional: $0) }, operations: [.identity])
    }()

    @Test func aluminiumDown001ShowsTheZeroLayer() throws {
        let o = try CrystalOrientation(zoneAxis: SIMD3(0, 0, 1), cell: aluminium.cell)
        let settings = SpotPatternSimulator.Settings(sigmaExcitationError: 0.02, excitationMultiple: 3, minimumIntensity: 1e-5)
        let p = try SpotPattern.simulate(aluminium, orientation: o, volts: 200_000, kMax: 1.2, settings: settings)
        let lambda = ElectronWavelength.angstroms(volts: 200_000)
        #expect(p.wavelength == lambda)
        let hkls = Set(p.spots.map(\.hkl))
        #expect(!hkls.contains(MillerIndex(0, 0, 0)))
        for want in [MillerIndex(2, 0, 0), MillerIndex(0, 2, 0), MillerIndex(-2, 0, 0), MillerIndex(2, 2, 0), MillerIndex(4, 0, 0)] {
            #expect(hkls.contains(want), "\(want)")
        }
        // Mixed parity is forbidden in fcc (|F| = 0, below the cut); l ≠ 0 lies far off the sphere (s about l/a).
        #expect(!hkls.contains(MillerIndex(1, 0, 0)) && !hkls.contains(MillerIndex(1, 1, 0)))
        #expect(p.spots.allSatisfy { $0.hkl.l == 0 })
        for spot in p.spots {
            let g = Double(spot.hkl.h * spot.hkl.h + spot.hkl.k * spot.hkl.k) / (4.0495 * 4.0495)   // |g|² in the zero layer
            #expect(abs(simd_length_squared(spot.position) - g) <= 1e-14)
            #expect(abs(spot.d - 1 / g.squareRoot()) <= 1e-14 * spot.d)
            // g_z = 0: s = −λ|g|²/2 exactly.
            #expect(abs(spot.excitationError + lambda * g / 2) <= 1e-15 * abs(spot.excitationError))
            let envelope = exp(-spot.excitationError * spot.excitationError / (2 * 0.02 * 0.02))
            #expect(abs(spot.intensity - spot.structureFactorIntensity * envelope) <= 1e-15 * spot.intensity)
            #expect(abs(spot.excitationError) <= 0.06 && spot.intensity > 1e-5)
        }
        // Friedel mates sit opposite each other.
        let s200 = try #require(p.spots.first { $0.hkl == MillerIndex(2, 0, 0) })
        let m200 = try #require(p.spots.first { $0.hkl == MillerIndex(-2, 0, 0) })
        #expect(simd_length(s200.position + m200.position) <= 1e-15)
    }

    @Test func anInPlaneTurnTurnsEverySpot() throws {
        let a = try CrystalOrientation(zoneAxis: SIMD3(1, 1, 2), inPlaneAngle: 0.2, cell: aluminium.cell)
        let b = a.rotatedInPlane(by: 0.5)
        let simulator = try SpotPatternSimulator(structure: aluminium, kMax: 1.5)
        let pa = try simulator.pattern(orientation: a, volts: 300_000), pb = try simulator.pattern(orientation: b, volts: 300_000)
        #expect(pa.spots.map(\.hkl) == pb.spots.map(\.hkl))
        #expect(pa.spots.count >= 8)
        let c = cos(0.5), s = sin(0.5)
        for (x, y) in zip(pa.spots, pb.spots) {
            // m3z on the right turns the pattern counter-clockwise: (x, y) → (c x − s y, s x + c y).
            let turned = SIMD2(c * x.position.x - s * x.position.y, s * x.position.x + c * x.position.y)
            #expect(simd_length(turned - y.position) <= 1e-14)
            #expect(abs(x.excitationError - y.excitationError) <= 1e-16)
        }
        // The simulator and the one-shot call give the same pattern.
        #expect(try SpotPattern.simulate(aluminium, orientation: a, volts: 300_000, kMax: 1.5) == pa)
    }

    @Test func radiusIsLTanTwoThetaAndPointsAlongTheSpot() throws {
        let lambda = ElectronWavelength.angstroms(volts: 200_000)
        for d in [0.5, 1.2, 2.338, 4.0495] {
            let r = SpotPattern.radius(d: d, cameraLength: 800, wavelength: lambda)
            #expect(abs(r - 800 * tan(2 * asin(lambda / (2 * d)))) <= 1e-13 * r)
            // Small angles: R ≈ Lλ/d.
            #expect(abs(r / (800 * lambda / d) - 1) <= 1e-3)
        }
        #expect(SpotPattern.radius(d: 0, cameraLength: 800, wavelength: lambda).isNaN)
        #expect(SpotPattern.radius(d: 2, cameraLength: -1, wavelength: lambda).isNaN)
        #expect(SpotPattern.radius(d: .nan, cameraLength: 800, wavelength: lambda).isNaN)
        // λ/2d at or past 1/√2 (2θ ≥ 90°): no flat detector sees it.
        #expect(SpotPattern.radius(d: 0.7, cameraLength: 800, wavelength: 1.0).isNaN)

        let o = try CrystalOrientation(zoneAxis: SIMD3(0, 1, 1), inPlaneAngle: 1.0, cell: aluminium.cell)
        let p = try SpotPattern.simulate(aluminium, orientation: o, volts: 200_000)
        for spot in p.spots {
            let at = try #require(p.detectorPosition(of: spot, cameraLength: 1000))
            #expect(abs(simd_length(at) - p.radius(of: spot, cameraLength: 1000)) <= 1e-12 * simd_length(at))
            #expect(abs(simd_dot(simd_normalize(at), simd_normalize(spot.position)) - 1) <= 1e-14)
        }
    }

    @Test func refusals() throws {
        let o = try CrystalOrientation(zoneAxis: SIMD3(0, 0, 1), cell: aluminium.cell)
        let simulator = try SpotPatternSimulator(structure: aluminium)
        for v in [0.0, -200_000, .nan, .infinity] {
            do {
                _ = try simulator.pattern(orientation: o, volts: v)
                Issue.record("\(v) V was accepted")
            } catch {
                guard case .invalidVoltage = error else { Issue.record("\(v) V: \(error)"); continue }
            }
        }
        for bad in [SpotPatternSimulator.Settings(sigmaExcitationError: 0), .init(sigmaExcitationError: .nan),
                    .init(excitationMultiple: -1), .init(minimumIntensity: -1e-5), .init(minimumIntensity: .infinity)] {
            #expect(throws: SpotPatternError.self) { try simulator.pattern(orientation: o, volts: 200_000, settings: bad) }
        }
        // No form factor for Z = 120: the kinematic list refuses it, before any F.
        let unknown = CrystalStructure(cell: .cubic(a: 4), sites: [AtomSite(z: 120, fractional: .zero)], operations: [.identity])
        #expect(throws: SpotPatternError.kinematic(.unsupportedElements([120]))) { try SpotPatternSimulator(structure: unknown) }
        #expect(throws: SpotPatternError.self) { try SpotPatternSimulator(structure: aluminium, kMax: 0) }
    }
}
