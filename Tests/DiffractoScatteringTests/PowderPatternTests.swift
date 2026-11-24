//
//  PowderPatternTests.swift
//  The X-ray powder pattern (PowderPattern.xray) and pymatgen's X-ray form factor (XRayScatteringFactors) on closed forms
//  and refusals. pymatgen's own numbers are pin P49 (DiffractoKitParityTests/Scattering/PowderXRDPinsTests.swift).
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoScattering

struct PowderPatternTests {

    static func fcc(a: Double, z: Int) -> CrystalStructure {
        let sites: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        return CrystalStructure(cell: .cubic(a: a), sites: sites.map { AtomSite(z: z, fractional: $0) }, operations: [])
    }

    static func hcp(a: Double, c: Double, z: Int) -> CrystalStructure {
        let sites: [SIMD3<Double>] = [[1.0 / 3, 2.0 / 3, 0.25], [2.0 / 3, 1.0 / 3, 0.75]]
        return CrystalStructure(cell: .hexagonal(a: a, c: c), sites: sites.map { AtomSite(z: z, fractional: $0) }, operations: [])
    }

    /// f(Z, 0) = Z for every element of the table; f falls with s; the g² form is s = |g| / 2; outside the table nil.
    /// RED UNDER: the 41.78214 constant or the sign; s = |g|.
    @Test func formFactorClosedForms() throws {
        for z in 1...XRayScatteringFactors.maximumAtomicNumber {
            #expect(XRayScatteringFactors.f(z: z, sSquared: 0) == Double(z))
        }
        let near = try #require(XRayScatteringFactors.f(z: 13, sSquared: 0.01))
        let far = try #require(XRayScatteringFactors.f(z: 13, sSquared: 0.25))
        #expect(near < 13 && far < near && far > 0)
        let viaG = try #require(XRayScatteringFactors.f(z: 79, gSquared: 0.36))
        let viaS = try #require(XRayScatteringFactors.f(z: 79, sSquared: 0.09))
        #expect(abs(viaG - viaS) <= 1e-12 * viaS)
        #expect(XRayScatteringFactors.f(z: 0, sSquared: 0) == nil)
        #expect(XRayScatteringFactors.f(z: 99, sSquared: 0) == nil)
        #expect(XRayScatteringFactors.f(z: 13, sSquared: -1) == nil)
        #expect(XRayScatteringFactors.f(z: 13, sSquared: .nan) == nil)
    }

    /// Al with Cu Kα: 111 at 2θ = 2 asin(λ√3 / 2a), d = a/√3, the strongest (100), multiplicity 8; five peaks to 90°
    /// (111, 200, 220, 311, 222); 200 has multiplicity 6.
    /// RED UNDER: θ for 2θ; a half-space without the factor 2 changes nothing here, the multiplicity counting only one sign
    /// does (4 for 8).
    @Test func aluminiumClosedForm() throws {
        let a = 4.0495
        let pattern = try PowderPattern.xray(structure: Self.fcc(a: a, z: 13))
        #expect(pattern.axis == .twoThetaDegrees && pattern.wavelengthAngstroms == 1.54184 && pattern.range == 0...90)
        #expect(pattern.peaks.count == 5)
        let first = try #require(pattern.peaks.first)
        let twoTheta = 2 * asin(1.54184 * 3.0.squareRoot() / (2 * a)) * 180 / Double.pi
        #expect(abs(first.position - twoTheta) <= 1e-10)
        #expect(abs(first.d - a / 3.0.squareRoot()) <= 1e-12)
        #expect(abs(first.intensity - 100) <= 1e-12 && first.multiplicity == 8)
        #expect(first.families == [PowderPattern.Family(indices: [1, 1, 1], multiplicity: 8)])
        #expect(pattern.peaks.map { $0.families.first?.indices ?? [] } == [[1, 1, 1], [2, 0, 0], [2, 2, 0], [3, 1, 1], [2, 2, 2]])
        #expect(pattern.peaks[1].multiplicity == 6)
        #expect(zip(pattern.peaks, pattern.peaks.dropFirst()).allSatisfy { $0.position < $1.position })
    }

    /// A hexagonal lattice is labelled with four indices: Mg's first peak is 1 0 -1 0 with multiplicity 6.
    @Test func hexagonalLabelsHaveFourIndices() throws {
        let pattern = try PowderPattern.xray(structure: Self.hcp(a: 3.2094, c: 5.2108, z: 12))
        let first = try #require(pattern.peaks.first)
        #expect(first.families == [PowderPattern.Family(indices: [1, 0, -1, 0], multiplicity: 6)])
        #expect(pattern.peaks.allSatisfy { peak in
            peak.families.allSatisfy { f in f.indices.count == 4 && f.indices[2] == -f.indices[0] - f.indices[1] }
        })
    }

    /// A range with no reflection gives an empty pattern; a lower bound drops the peaks below it. A range that holds only
    /// absent reflections gives them scaled to 100, as pymatgen does (Al at 0 to 30°: the forbidden 100 at 21.95°, its
    /// absolute intensity rounding residue): the scale is relative to what the range holds.
    @Test func rangeBounds() throws {
        let al = Self.fcc(a: 4.0495, z: 13)
        #expect(try PowderPattern.xray(structure: al, twoThetaRange: 0...20).peaks.isEmpty)
        let absent = try PowderPattern.xray(structure: al, twoThetaRange: 0...30)
        #expect(absent.peaks.count == 1 && abs((absent.peaks.first?.position ?? 0) - 21.949) < 1e-3)
        #expect((absent.peaks.first?.absoluteIntensity ?? 1) < 1e-20)
        let upper = try PowderPattern.xray(structure: al, twoThetaRange: 40...90)
        #expect(upper.peaks.count == 4 && upper.peaks.allSatisfy { $0.position >= 40 })
        #expect(abs((upper.strongest?.intensity ?? 0) - 100) <= 1e-12)
    }

    /// A Debye-Waller factor damps high angles: Au's last peak loses relative intensity with B = 0.62 Å².
    @Test func debyeWallerDampsHighAngles() throws {
        let au = Self.fcc(a: 4.0782, z: 79)
        let cold = try PowderPattern.xray(structure: au)
        let warm = try PowderPattern.xray(structure: au, debyeWallerFactors: [79: 0.62])
        #expect(cold.peaks.count == warm.peaks.count)
        let c = try #require(cold.peaks.last), w = try #require(warm.peaks.last)
        #expect(w.intensity < c.intensity && w.absoluteIntensity < c.absoluteIntensity)
    }

    /// Refusals: an element outside the table, a wavelength that is not a positive number, a range outside 0 to 180° or
    /// empty, a Debye-Waller factor that is not finite.
    @Test func refusals() throws {
        let es = CrystalStructure(cell: .cubic(a: 5), sites: [AtomSite(z: 99, fractional: .zero)], operations: [])
        #expect(throws: PowderPatternError.unsupportedElements([99])) { try PowderPattern.xray(structure: es) }
        let al = Self.fcc(a: 4.0495, z: 13)
        for w in [0.0, -1.0, Double.nan, Double.infinity] {
            do throws(PowderPatternError) {
                _ = try PowderPattern.xray(structure: al, radiation: XRayRadiation(wavelengthAngstroms: w))
                Issue.record("wavelength \(w) accepted")
            } catch {
                guard case .invalidWavelength(let got) = error else { Issue.record("\(error)"); continue }
                #expect(got.isNaN ? w.isNaN : got == w)
            }
        }
        #expect(throws: PowderPatternError.invalidRange(lower: -5, upper: 90)) {
            try PowderPattern.xray(structure: al, twoThetaRange: -5...90)
        }
        #expect(throws: PowderPatternError.invalidRange(lower: 0, upper: 181)) {
            try PowderPattern.xray(structure: al, twoThetaRange: 0...181)
        }
        #expect(throws: PowderPatternError.invalidRange(lower: 30, upper: 30)) {
            try PowderPattern.xray(structure: al, twoThetaRange: 30...30)
        }
        #expect(throws: PowderPatternError.invalidDebyeWallerFactor(z: 13, value: .infinity)) {
            try PowderPattern.xray(structure: al, debyeWallerFactors: [13: .infinity])
        }
    }

    /// pymatgen's lines by name.
    @Test func namedLines() {
        #expect(XRayRadiation.pymatgenLines.count == 24)
        #expect(XRayRadiation.named("CuKa") == .cuKAlpha)
        #expect(XRayRadiation.named("MoKa")?.wavelengthAngstroms == 0.71073)
        #expect(XRayRadiation.named("XxKa") == nil)
    }
}
