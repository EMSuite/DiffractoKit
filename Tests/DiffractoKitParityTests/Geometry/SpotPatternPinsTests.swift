//
//  SpotPatternPinsTests.swift
//  Role: `SpotPattern` (S2) against its references, read from Tests/DiffractoGeometryTests/Fixtures:
//    spot-pattern-pins.json (tools/pins/spot-pattern/pins-py4dstem.py, py4DSTEM f050d207): 7 phases at 80, 200 and 300 kV,
//      9 orientations each (zone axes turned, the same tilted, random): every spot's h k l, qx, qy, intensity and excitation
//      error from `generate_diffraction_pattern` and `excitation_errors`;
//    spot-pattern-diffsims-pins.json (tools/pins/spot-pattern/pins-diffsims.py, diffsims 0.7.0): the second check, 7 phases
//      x 8 random rotations at 200 kV;
//    spot-radius-pins.json (tools/pins/spot-pattern/radius-pins.py, mpmath at 60 digits): R = L tan 2θ, 540 cases.
//  Tolerances (docs/parity.md P42): positions rtol 1e-7 (assert_allclose's float64 default; py4DSTEM ships no tests) with
//  a floor of 1e-12 Å⁻¹; intensities 5e-5 relative (py4DSTEM keeps F in complex64); excitation errors 1e-8 Å⁻¹ absolute
//  (py4DSTEM's wavelength uses 6- and 7-digit constants, 1e-7 relative from CODATA's: 5e-9 Å⁻¹ at |g| 1.5).
//  diffsims: positions 1e-7, |F|² 1e-6 relative, s 1e-5 Å⁻¹ (plan row S2), its two documented deviations left out.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoGeometry
import DiffractoScattering

private struct PinCell: Decodable { let a, b, c, alpha, beta, gamma: Double }
private struct PinSite: Decodable { let frac: [Double]; let occupancy: Double; let z: Int }
private struct PinPattern: Decodable {
    let label: String
    let volts: Double
    let wavelength: Double
    let matrix: [Double]
    /// h, k, l, qx, qy, intensity, excitation error
    let spots: [[Double]]
}
private struct PinPhase: Decodable { let name: String; let cell: PinCell; let sites: [PinSite]; let patterns: [PinPattern] }
private struct PinSettings: Decodable {
    let kMax, sigmaExcitationError, tolExcitationErrorMult, tolIntensity, tolStructureFactor: Double
}
private struct Py4DSTEMPins: Decodable { let versions: [String: String]; let settings: PinSettings; let phases: [PinPhase] }

private struct DiffsimsPattern: Decodable {
    let matrix: [Double]
    let kabschRms: Double
    /// h, k, l, x, y, z, intensity, |F|²
    let spots: [[Double]]
}
private struct DiffsimsPhase: Decodable { let name: String; let patterns: [DiffsimsPattern] }
private struct DiffsimsPins: Decodable {
    let versions: [String: String]
    let reciprocalRadius, maxExcitationError, volts, wavelength: Double
    let phases: [DiffsimsPhase]
}
private struct RadiusPins: Decodable { let digits: Int; let cases: [[Double]] }

private enum Pins {
    static func load<T: Decodable>(_ name: String) -> T {
        let url = FixtureDirs.fixtures(of: "DiffractoGeometryTests").appendingPathComponent(name)
        let decoder = JSONDecoder()
        decoder.keyDecodingStrategy = .convertFromSnakeCase
        do { return try decoder.decode(T.self, from: Data(contentsOf: url)) } catch { fatalError("\(name): \(error)") }
    }
    static let py4dstem: Py4DSTEMPins = load("spot-pattern-pins.json")
    static let diffsims: DiffsimsPins = load("spot-pattern-diffsims-pins.json")
    static let radius: RadiusPins = load("spot-radius-pins.json")
}

private func matrix(_ r: [Double]) -> simd_double3x3 {
    simd_double3x3(rows: [SIMD3(r[0], r[1], r[2]), SIMD3(r[3], r[4], r[5]), SIMD3(r[6], r[7], r[8])])
}

private func structure(_ p: PinPhase) -> CrystalStructure {
    let c = p.cell
    return CrystalStructure(cell: UnitCell(a: c.a, b: c.b, c: c.c, alphaDeg: c.alpha, betaDeg: c.beta, gammaDeg: c.gamma),
                            sites: p.sites.map { AtomSite(z: $0.z, fractional: SIMD3($0.frac[0], $0.frac[1], $0.frac[2]),
                                                          occupancy: $0.occupancy) },
                            operations: [.identity])
}

private func key(_ r: [Double]) -> MillerIndex { MillerIndex(Int(r[0]), Int(r[1]), Int(r[2])) }

@Suite("SpotPattern against py4DSTEM f050d207, diffsims 0.7.0 and 60-digit R = L tan 2θ")
struct SpotPatternPinsTests {

    @Test func pinsCarryTheirVersionsAndSizes() {
        #expect(Pins.py4dstem.versions["py4DSTEM"] == "0.14.19")
        #expect(Pins.py4dstem.phases.count == 7)
        #expect(Pins.py4dstem.phases.allSatisfy { $0.patterns.count == 27 })
        #expect(Pins.diffsims.versions["diffsims"] == "0.7.0")
        #expect(Pins.diffsims.phases.count == 7)
        #expect(Pins.radius.cases.count == 540 && Pins.radius.digits == 60)
    }

    /// The same spots as `generate_diffraction_pattern` (000 left out), each at the same place with the same intensity and
    /// excitation error. Every py4DSTEM spot is the kit's; a spot the kit alone has must sit within the rounding of a cut:
    /// |s| at σ·multiple, the intensity at the floor, |F| at the structure-factor cut, or |g| at k_max.
    @Test func spotsEqualPy4DSTEMsPatterns() throws {
        let pins = Pins.py4dstem
        let s = pins.settings
        let settings = SpotPatternSimulator.Settings(sigmaExcitationError: s.sigmaExcitationError,
                                                     excitationMultiple: s.tolExcitationErrorMult,
                                                     minimumIntensity: s.tolIntensity)
        let sgMax = s.sigmaExcitationError * s.tolExcitationErrorMult
        var worst = (position: 0.0, intensity: 0.0, excitation: 0.0), count = 0, edge = 0
        for phase in pins.phases {
            let simulator = try SpotPatternSimulator(structure: structure(phase), kMax: s.kMax,
                                                     structureFactorCut: s.tolStructureFactor)
            for p in phase.patterns {
                let o = try CrystalOrientation(matrix: matrix(p.matrix))
                let pattern = try simulator.pattern(orientation: o, volts: p.volts, settings: settings)
                #expect(abs(pattern.wavelength - p.wavelength) <= 2e-7 * p.wavelength, "\(phase.name): λ")
                var ours: [MillerIndex: SimulatedSpot] = [:]
                for spot in pattern.spots { ours[spot.hkl] = spot }
                let label = "\(phase.name) \(p.label) \(p.volts)"
                func nearACut(_ spot: SimulatedSpot) -> Bool {
                    abs(abs(spot.excitationError) - sgMax) <= 1e-7
                        || abs(spot.intensity - s.tolIntensity) <= 1e-4 * s.tolIntensity
                        || abs(spot.structureFactorIntensity.squareRoot() - s.tolStructureFactor) <= 1e-4 * s.tolStructureFactor
                        || abs(1 / spot.d - s.kMax) <= 1e-12
                }
                var theirs = Set<MillerIndex>()
                for r in p.spots {
                    let hkl = key(r)
                    theirs.insert(hkl)
                    guard let spot = ours[hkl] else {
                        Issue.record("\(label): \(hkl) is py4DSTEM's spot only (s \(r[6]), intensity \(r[5]))")
                        continue
                    }
                    count += 1
                    let dx = abs(spot.position.x - r[3]), dy = abs(spot.position.y - r[4])
                    worst.position = max(worst.position, dx, dy)
                    #expect(dx <= 1e-7 * abs(r[3]) + 1e-12 && dy <= 1e-7 * abs(r[4]) + 1e-12,
                            "\(label) \(hkl): (\(spot.position.x), \(spot.position.y)) against (\(r[3]), \(r[4]))")
                    let di = abs(spot.intensity - r[5]) / r[5]
                    worst.intensity = max(worst.intensity, di)
                    #expect(di <= 5e-5, "\(label) \(hkl): intensity \(spot.intensity) against \(r[5])")
                    let ds = abs(spot.excitationError - r[6])
                    worst.excitation = max(worst.excitation, ds)
                    #expect(ds <= 1e-8, "\(label) \(hkl): s \(spot.excitationError) against \(r[6])")
                }
                for spot in pattern.spots where !theirs.contains(spot.hkl) {
                    edge += 1
                    #expect(nearACut(spot), "\(label): \(spot.hkl) is the kit's spot only (s \(spot.excitationError))")
                }
            }
        }
        #expect(count > 5000)
        print("S2WORST py4DSTEM position \(worst.position) intensity \(worst.intensity) s \(worst.excitation) over \(count) spots; kit-only at a cut \(edge)")
    }

    /// diffsims: the same positions in its own lab frame, the same |F|² (×V², diffsims' F has no 1/V), and the kit's s equal
    /// to py4DSTEM's form evaluated on diffsims' own (x, y, z), which is −s of diffsims' exact form to 1.1e-5 Å⁻¹ at |s| ≤ 0.03
    /// (deviation 2, measured in the generator's notes). Reflections with |F| at or below the kit's cut are left out
    /// (deviation 1: diffsims keeps them).
    @Test func spotsAgreeWithDiffsims() throws {
        let pins = Pins.diffsims
        let phases = Dictionary(uniqueKeysWithValues: Pins.py4dstem.phases.map { ($0.name, $0) })
        let lambda = ElectronWavelength.angstroms(volts: pins.volts)
        #expect(abs(lambda / pins.wavelength - 1) <= 1e-9)
        var worst = (position: 0.0, f2: 0.0, s: 0.0, form: 0.0), count = 0, leftOut = 0
        for d in pins.phases {
            let phase = try #require(phases[d.name])
            let crystal = structure(phase)
            let simulator = try SpotPatternSimulator(structure: crystal, kMax: pins.reciprocalRadius, structureFactorCut: 1e-4)
            let settings = SpotPatternSimulator.Settings(sigmaExcitationError: pins.maxExcitationError / 3, excitationMultiple: 3,
                                                         minimumIntensity: 0)
            let volume2 = crystal.cell.volume * crystal.cell.volume
            for p in d.patterns {
                #expect(p.kabschRms <= 1e-9)
                let pattern = try simulator.pattern(orientation: try CrystalOrientation(matrix: matrix(p.matrix)), volts: pins.volts,
                                                    settings: settings)
                var ours: [MillerIndex: SimulatedSpot] = [:]
                for spot in pattern.spots { ours[spot.hkl] = spot }
                for r in p.spots {
                    let xyz = SIMD3(r[3], r[4], r[5])
                    // py4DSTEM's form on diffsims' geometry; within 1e-5 of the kit's cut it may fall on either side.
                    let sPy = (2 * xyz.z - lambda * simd_length_squared(xyz)) / (2 - 2 * lambda * xyz.z)
                    guard let spot = ours[key(r)] else {
                        leftOut += 1
                        let f = StructureFactorProbe.magnitude(key(r), in: crystal)
                        // py4DSTEM's form may pass 0.03 where diffsims' exact one does not (the forms differ by up to 1.13e-5).
                        #expect(f <= 1e-4 * (1 + 1e-6)
                                    || (abs(sPy) > pins.maxExcitationError && abs(sPy) - pins.maxExcitationError <= 1.2e-5),
                                "\(d.name): diffsims' \(key(r)) is missing (|F| \(f), s \(sPy))")
                        continue
                    }
                    count += 1
                    let dp = max(abs(spot.position.x - xyz.x), abs(spot.position.y - xyz.y))
                    worst.position = max(worst.position, dp)
                    #expect(dp <= 1e-7 * simd_length(SIMD2(xyz.x, xyz.y)) + 1e-12, "\(d.name) \(spot.hkl): position")
                    let df = abs(spot.structureFactorIntensity * volume2 / r[7] - 1)
                    worst.f2 = max(worst.f2, df)
                    #expect(df <= 1e-6, "\(d.name) \(spot.hkl): |F|² \(spot.structureFactorIntensity * volume2) against \(r[7])")
                    let ds = abs(spot.excitationError - sPy)
                    worst.s = max(worst.s, ds)
                    #expect(ds <= 1e-5, "\(d.name) \(spot.hkl): s \(spot.excitationError) against \(sPy)")
                    let r0 = 1 / pins.wavelength
                    let exact = r0 - (r0 * r0 - xyz.x * xyz.x - xyz.y * xyz.y).squareRoot() - xyz.z
                    worst.form = max(worst.form, abs(spot.excitationError + exact))
                    // The independent check of s: diffsims' exact distance to the sphere, opposite in sign (deviation 2; measured
                    // 1.127e-5 at |s| ≤ 0.03).
                    #expect(abs(spot.excitationError + exact) <= 1.2e-5, "\(d.name) \(spot.hkl): s \(spot.excitationError) against −\(exact)")
                }
            }
        }
        #expect(count > 1000)
        print("S2WORST diffsims position \(worst.position) |F|² \(worst.f2) s \(worst.s) over \(count) spots; left out \(leftOut); the two forms of s differ by \(worst.form)")
    }

    /// R = L tan 2θ, sin θ = λ/2d, against 60-digit arithmetic: 3 voltages × 5 camera lengths × 36 spacings.
    @Test func radiusIsLTanTwoTheta() {
        var worst = 0.0
        for c in Pins.radius.cases {
            let lambda = ElectronWavelength.angstroms(volts: c[0])
            #expect(abs(lambda - c[3]) <= 1e-7 * c[3])
            let r = SpotPattern.radius(d: c[2], cameraLength: c[1], wavelength: lambda)
            worst = max(worst, abs(r - c[4]) / c[4])
            #expect(abs(r - c[4]) <= 1e-7 * c[4], "\(c)")
        }
        print("S2WORST radius \(worst)")
    }
}

/// |F| of one reflection, for the left-out check.
private enum StructureFactorProbe {
    static func magnitude(_ hkl: MillerIndex, in s: CrystalStructure) -> Double {
        StructureFactors.compute(hkl, in: s)?.magnitude ?? .nan
    }
}
