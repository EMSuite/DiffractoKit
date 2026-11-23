//
//  XRayScatteringFactors.swift
//  Role: X-ray atomic scattering factors as pymatgen's XRDCalculator computes them,
//        f(s) = Z − 41.78214 s² Σᵢ aᵢ exp(−bᵢ s²), s = sin θ / λ = |g| / 2 (Å⁻¹), with pymatgen's four (a, b) pairs per
//        element (`XRayScatteringTable.swift`, written by tools/pins/powder-xrd/make-table.py). The arithmetic is
//        pymatgen's, in the same order. No anomalous dispersion, as in pymatgen.
//

import Foundation

/// X-ray atomic scattering factors (pymatgen's `XRDCalculator` table and formula).
public enum XRayScatteringFactors {

    /// The largest atomic number the table covers (Cf).
    public static let maximumAtomicNumber = 98

    public static func isSupported(z: Int) -> Bool { z >= 1 && z <= maximumAtomicNumber }

    /// f(Z, s²) in electrons, s² = (sin θ / λ)² in Å⁻²; nil for an element outside the table or an s² that is not a
    /// finite number ≥ 0.
    public static func f(z: Int, sSquared: Double) -> Double? {
        guard isSupported(z: z), sSquared.isFinite, sSquared >= 0 else { return nil }
        let p = parameters[z - 1]
        var sum = 0.0
        for i in 0..<4 { sum += p[2 * i] * exp(-p[2 * i + 1] * sSquared) }
        return Double(z) - 41.78214 * sSquared * sum
    }

    /// f(Z, g²) with g² = |g|² in Å⁻² (s = |g| / 2): the signature of a ``FormFactor``, so the X-ray factor can be passed
    /// to ``StructureFactors/compute(_:in:formFactor:)``.
    public static func f(z: Int, gSquared: Double) -> Double? {
        guard gSquared.isFinite, gSquared >= 0 else { return nil }
        let half = gSquared.squareRoot() / 2
        return f(z: z, sSquared: half * half)
    }
}

/// An X-ray wavelength: one of pymatgen's named lines, or any wavelength.
public struct XRayRadiation: Sendable, Hashable, Codable, CustomStringConvertible {
    /// pymatgen's key ("CuKa", "MoKa1" …), nil for a typed wavelength.
    public let name: String?
    /// λ in Å.
    public let wavelengthAngstroms: Double

    /// A typed wavelength in Å.
    public init(wavelengthAngstroms: Double) {
        self.name = nil
        self.wavelengthAngstroms = wavelengthAngstroms
    }

    private init(_ name: String, _ wavelength: Double) {
        self.name = name
        self.wavelengthAngstroms = wavelength
    }

    /// Cu Kα (the Kα1/Kα2 average), 1.54184 Å: pymatgen's default.
    public static let cuKAlpha = XRayRadiation("CuKa", 1.54184)
    /// Mo Kα, 0.71073 Å.
    public static let moKAlpha = XRayRadiation("MoKa", 0.71073)

    /// pymatgen's `WAVELENGTHS`, in its order (Å).
    public static let pymatgenLines: [XRayRadiation] = [
        XRayRadiation("CuKa", 1.54184), XRayRadiation("CuKa2", 1.54439), XRayRadiation("CuKa1", 1.54056),
        XRayRadiation("CuKb1", 1.39222),
        XRayRadiation("MoKa", 0.71073), XRayRadiation("MoKa2", 0.71359), XRayRadiation("MoKa1", 0.70930),
        XRayRadiation("MoKb1", 0.63229),
        XRayRadiation("CrKa", 2.29100), XRayRadiation("CrKa2", 2.29361), XRayRadiation("CrKa1", 2.28970),
        XRayRadiation("CrKb1", 2.08487),
        XRayRadiation("FeKa", 1.93735), XRayRadiation("FeKa2", 1.93998), XRayRadiation("FeKa1", 1.93604),
        XRayRadiation("FeKb1", 1.75661),
        XRayRadiation("CoKa", 1.79026), XRayRadiation("CoKa2", 1.79285), XRayRadiation("CoKa1", 1.78896),
        XRayRadiation("CoKb1", 1.63079),
        XRayRadiation("AgKa", 0.560885), XRayRadiation("AgKa2", 0.563813), XRayRadiation("AgKa1", 0.559421),
        XRayRadiation("AgKb1", 0.497082),
    ]

    /// The named line, by pymatgen's key.
    public static func named(_ name: String) -> XRayRadiation? { pymatgenLines.first { $0.name == name } }

    public var description: String {
        if let name { return "\(name) (\(wavelengthAngstroms) Å)" }
        return "\(wavelengthAngstroms) Å"
    }
}
