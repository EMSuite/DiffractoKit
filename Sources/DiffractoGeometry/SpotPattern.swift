//
//  SpotPattern.swift
//  Role: the kinematic spot pattern of a crystal structure at any orientation (DiffractoKit S2): py4DSTEM f050d207's
//        `Crystal.generate_diffraction_pattern` without precession or foil normal, every spot with hkl, position, d,
//        excitation error and intensity, and its distance on a film or detector at a camera length, R = L tan 2θ.
//
//  The reflections are the kit's kinematic list (`CrystalStructure.kinematicReflections`: every hkl with 0 < |g| ≤ k_max and
//  |F| above the cut, F from py4DSTEM's scattering-factor table). For an orientation M (columns lab x, y, z in the crystal's
//  Cartesian frame, `CrystalOrientation`) each g is seen at Mᵀg; its excitation error is py4DSTEM's
//      s = (2 g_z − λ |g|²) / (2 − 2 λ g_z)        (Crystal.excitation_errors, k₀ = (0, 0, −1/λ), foil normal −z),
//  and a spot is kept when |s| ≤ σ · multiple and |F|² exp(−s² / 2σ²) > the intensity floor. The direct beam (000) is not a
//  spot of the list. py4DSTEM's pins: P42 in docs/parity.md.
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering

/// One spot of a simulated pattern.
public struct SimulatedSpot: Sendable, Equatable {
    public let hkl: MillerIndex
    /// Where the spot appears: (qx, qy) = the lab x and y of Mᵀg, Å⁻¹ (no 2π). py4DSTEM's qx, qy.
    public let position: SIMD2<Double>
    /// d = 1/|g|, Å.
    public let d: Double
    /// s = (2 g_z − λ|g|²) / (2 − 2λ g_z), Å⁻¹: positive inside the Ewald sphere (py4DSTEM's sign).
    public let excitationError: Double
    /// |F|² exp(−s² / 2σ²), Å⁻⁴ (py4DSTEM's intensity).
    public let intensity: Double
    /// |F|², Å⁻⁴, before the excitation envelope.
    public let structureFactorIntensity: Double

    public init(hkl: MillerIndex, position: SIMD2<Double>, d: Double, excitationError: Double, intensity: Double,
                structureFactorIntensity: Double) {
        self.hkl = hkl
        self.position = position
        self.d = d
        self.excitationError = excitationError
        self.intensity = intensity
        self.structureFactorIntensity = structureFactorIntensity
    }
}

public enum SpotPatternError: Error, Equatable, Sendable, CustomStringConvertible {
    /// The structure's kinematic list was refused (an element without a form factor, an invalid structure, k_max).
    case kinematic(KinematicError)
    /// The voltage is not a finite number of volts above 0.
    case invalidVoltage(Double)
    /// A setting is out of its range; the text names it.
    case invalidSettings(String)

    public var description: String {
        switch self {
        case .kinematic(let e): return e.description
        case .invalidVoltage(let v): return "The accelerating voltage must be a finite number of volts above 0; \(v) is not."
        case .invalidSettings(let s): return "The spot pattern settings are not valid: \(s)."
        }
    }
}

extension SpotPatternError: LocalizedError {
    public var errorDescription: String? { description }
}

/// A structure's reflections, ready for patterns at any orientation and voltage. Build it once per structure.
public struct SpotPatternSimulator: Sendable {

    /// The excitation envelope and the cuts of `pattern(orientation:volts:settings:)`.
    public struct Settings: Sendable, Equatable {
        /// σ of the envelope exp(−s² / 2σ²), Å⁻¹. Default 0.02 (py4DSTEM's own default is 0.01; 0.02 is the design's choice
        /// for a selected-area pattern of a thin foil, which shows more of the higher-order Laue zones).
        public var sigmaExcitationError: Double
        /// Spots with |s| ≤ σ · multiple are kept. Default 3 (py4DSTEM's).
        public var excitationMultiple: Double
        /// Spots with an intensity above this are kept, Å⁻⁴. Default 1e-5 (py4DSTEM's `tol_intensity`).
        public var minimumIntensity: Double

        public init(sigmaExcitationError: Double = 0.02, excitationMultiple: Double = 3, minimumIntensity: Double = 1e-5) {
            self.sigmaExcitationError = sigmaExcitationError
            self.excitationMultiple = excitationMultiple
            self.minimumIntensity = minimumIntensity
        }

        func validated() throws(SpotPatternError) {
            guard sigmaExcitationError.isFinite, sigmaExcitationError > 0 else {
                throw .invalidSettings("sigmaExcitationError must be a finite number above 0")
            }
            guard excitationMultiple.isFinite, excitationMultiple > 0 else {
                throw .invalidSettings("excitationMultiple must be a finite number above 0")
            }
            guard minimumIntensity.isFinite, minimumIntensity >= 0 else {
                throw .invalidSettings("minimumIntensity must be a finite number of at least 0")
            }
        }
    }

    public let structure: CrystalStructure
    public let kMax: Double
    public let structureFactorCut: Double
    /// The kinematic reflections above the cut, in the kinematic list's order (by |g|, then DiffractoCore's hkl order).
    public let reflections: [KinematicReflection]

    /// The reflections of `structure` with 0 < |g| ≤ `kMax` (Å⁻¹) and |F| > `structureFactorCut` (Å⁻², the kit's and py4DSTEM's
    /// default 1e-4). Throws what the kinematic list refuses.
    public init(structure: CrystalStructure, kMax: Double = 2, structureFactorCut: Double = StructureFactors.defaultCut)
        throws(SpotPatternError) {
        do {
            let list = try structure.kinematicReflections(kMax: kMax, cut: structureFactorCut, retainingAll: false)
            reflections = list.kept
        } catch {
            throw .kinematic(error)
        }
        self.structure = structure
        self.kMax = kMax
        self.structureFactorCut = structureFactorCut
    }

    /// The pattern at `orientation` and `volts` (200 kV is 200_000).
    public func pattern(orientation: CrystalOrientation, volts: Double, settings: Settings = .init())
        throws(SpotPatternError) -> SpotPattern {
        try settings.validated()
        let lambda = ElectronWavelength.angstroms(volts: volts)
        guard lambda.isFinite, lambda > 0 else { throw .invalidVoltage(volts) }
        let sgMax = settings.sigmaExcitationError * settings.excitationMultiple
        let minusTwoSigmaSquared = -2 * settings.sigmaExcitationError * settings.sigmaExcitationError
        var spots: [SimulatedSpot] = []
        for r in reflections {
            let f = r.structureFactor
            let p = orientation.labCoordinates(of: f.g)
            let gSquared = p.x * p.x + p.y * p.y + p.z * p.z
            let s = (2 * p.z - lambda * gSquared) / (2 - 2 * lambda * p.z)
            guard s.isFinite, abs(s) <= sgMax else { continue }
            let intensity = f.intensity * exp(s * s / minusTwoSigmaSquared)
            guard intensity > settings.minimumIntensity else { continue }
            spots.append(SimulatedSpot(hkl: f.hkl, position: SIMD2(p.x, p.y), d: f.d, excitationError: s,
                                       intensity: intensity, structureFactorIntensity: f.intensity))
        }
        return SpotPattern(orientation: orientation, volts: volts, wavelength: lambda, settings: settings, spots: spots)
    }
}

/// The spots of one structure at one orientation and voltage.
public struct SpotPattern: Sendable, Equatable {
    public let orientation: CrystalOrientation
    public let volts: Double
    /// λ, Å (`ElectronWavelength`).
    public let wavelength: Double
    public let settings: SpotPatternSimulator.Settings
    /// The spots, in the order of `SpotPatternSimulator.reflections` (by |g|, then DiffractoCore's hkl order).
    public let spots: [SimulatedSpot]

    /// The pattern of `structure` at `orientation`: `SpotPatternSimulator(structure:kMax:structureFactorCut:)` then
    /// `pattern(orientation:volts:settings:)`. Build the simulator once to simulate many orientations of one structure.
    public static func simulate(_ structure: CrystalStructure, orientation: CrystalOrientation, volts: Double, kMax: Double = 2,
                                settings: SpotPatternSimulator.Settings = .init()) throws(SpotPatternError) -> SpotPattern {
        try SpotPatternSimulator(structure: structure, kMax: kMax).pattern(orientation: orientation, volts: volts,
                                                                            settings: settings)
    }

    /// The distance of a spot of spacing d from the direct beam on a flat detector at camera length L: R = L tan 2θ with
    /// sin θ = λ/2d, in the unit of L (mm in, mm out). NaN where 2θ reaches 90° or an input is not finite and positive
    /// (`RingCalibration.radius(q:cameraLengthPixels:wavelengthAngstroms:)`, P8 and P24).
    public static func radius(d: Double, cameraLength: Double, wavelength: Double) -> Double {
        guard d.isFinite, d > 0 else { return .nan }
        return RingCalibration.radius(q: 1 / d, cameraLengthPixels: cameraLength, wavelengthAngstroms: wavelength)
    }

    /// R = L tan 2θ of `spot` at camera length L (in the unit of L).
    public func radius(of spot: SimulatedSpot, cameraLength: Double) -> Double {
        Self.radius(d: spot.d, cameraLength: cameraLength, wavelength: wavelength)
    }

    /// Where `spot` lies on the detector, in the unit of L: R along the direction of its (qx, qy). nil when R is NaN or the
    /// spot has no direction in the pattern plane.
    public func detectorPosition(of spot: SimulatedSpot, cameraLength: Double) -> SIMD2<Double>? {
        let r = radius(of: spot, cameraLength: cameraLength)
        let q = simd_length(spot.position)
        guard r.isFinite, q > 0, q.isFinite else { return nil }
        return spot.position * (r / q)
    }
}
