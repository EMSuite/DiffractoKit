//
//  PowderRingsTests.swift
//  DiffractoKit S10 (#165): the electron ring pattern from the held ring list, the drawn profile as its closed form, and the
//  overlap and diagnostic marks on planted lists and on Al + Mg2Si + Si with Cu Kα at 0.1° (the design's numbers, from
//  pymatgen 2026.9.24 peak lists, P49).
//

import Testing
import Foundation
import DiffractoCore
import DiffractoScattering

@Suite struct PowderRingsTests {
    static let lambda200kV = 0.025079

    static var mg2Si: CrystalStructure {
        let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        var sites = fcc.map { AtomSite(z: 14, fractional: $0) }
        for p in fcc { for o in [SIMD3(0.25, 0.25, 0.25), SIMD3(0.75, 0.75, 0.75)] {
            let f = p + o
            sites.append(AtomSite(z: 12, fractional: f - f.rounded(.down)))
        } }
        return CrystalPresets.structure(a: 6.351, b: 6.351, c: 6.351, sites: sites)
    }

    /// The ring pattern is the held ring list: one peak per ring, at |g|, with its Σ|F|², members and orbits; strongest 100.
    /// Mutation: the multiplicity taken from the orbit count instead of the members fails Al {111} (8).
    @Test func ringsAreTheHeldRingList() throws {
        let al = CrystalPresets.aluminum
        let rings = try PowderPattern.electron(structure: al, wavelengthAngstroms: Self.lambda200kV, kMax: 2)
        let held = try al.kinematicReflections(kMax: 2).families
        #expect(rings.axis == .scatteringVector && rings.range == 0...2 && rings.wavelengthAngstroms == Self.lambda200kV)
        #expect(rings.peaks.count == held.count && !held.isEmpty)
        for (p, f) in zip(rings.peaks, held) {
            #expect(p.position == f.family.g && p.d == f.family.d && p.absoluteIntensity == f.intensity)
            #expect(p.multiplicity == f.family.multiplicity)
            #expect(p.families.map { MillerIndex($0.indices[0], $0.indices[1], $0.indices[2]) } == f.family.orbits)
        }
        #expect(rings.strongest?.intensity == 100)
        #expect(rings.peaks[0].families == [PowderPattern.Family(indices: [1, 1, 1], multiplicity: 8)])
        // {333} and {511} share |g| in a cubic cell: two orbits, 8 and 24 members.
        let shared = try #require(rings.peaks.first { $0.families.count == 2 })
        #expect(Set(shared.families.map(\.multiplicity)) == [8, 24] && shared.multiplicity == 32)
        // A hexagonal cell's rings are labelled in four indices, i = −(h + k).
        let mg = try PowderPattern.electron(structure: CrystalPresets.magnesium, wavelengthAngstroms: Self.lambda200kV, kMax: 1)
        #expect(!mg.peaks.isEmpty && mg.peaks.allSatisfy { $0.families.allSatisfy { $0.indices.count == 4 && $0.indices[2] == -$0.indices[0] - $0.indices[1] } })
        #expect(throws: PowderPatternError.invalidWavelength(0)) {
            try PowderPattern.electron(structure: al, wavelengthAngstroms: 0)
        }
        let unknown = CrystalPresets.structure(a: 4, b: 4, c: 4, sites: [AtomSite(z: 150, fractional: [0, 0, 0])])
        #expect(throws: PowderPatternError.invalidRingList(KinematicError.unsupportedElements([150]).description)) {
            try PowderPattern.electron(structure: unknown, wavelengthAngstroms: Self.lambda200kV)
        }
    }

    /// A ring's `intensity` is its Σ|F|² (multiplicity × |F|²) and its `profileIntensity` is Σ|F|² · d², each with the strongest
    /// ring at 100 (#223). Au at 200 kV is fcc: every allowed member has F = 4 f_e(|g|²) / V, so the closed form of a ring
    /// with h² + k² + l² = N is m · 16 f_e² / V² · a² / N, from the Lobato factors alone. 111 : 200 : 220 : 311 is
    /// 100 : 46.26 : 25.66 : 26.93 in the profile and 100 : 61.67 : 68.42 : 98.75 in Σ|F|². Mutation: d in place of d² moves
    /// every ratio; marks that floor by `intensity` keep Au 222 (29.94 in Σ|F|², 7.49 in the profile) above a floor of 20.
    @Test func profileIntensitiesAreSumFSquaredTimesDSquared() throws {
        let a = 4.0782, volume = a * a * a
        let au = try PowderPattern.electron(structure: CrystalPresets.gold, wavelengthAngstroms: Self.lambda200kV, kMax: 2)
        func closedForm(_ n: Int, _ m: Int) throws -> (sum: Double, profile: Double) {
            let f = try #require(ElectronScatteringFactors.f(z: 79, gSquared: Double(n) / (a * a)))
            let sum = Double(m) * 16 * f * f / (volume * volume)
            return (sum, sum * a * a / Double(n))
        }
        let rings: [(n: Int, m: Int)] = [(3, 8), (4, 6), (8, 12), (11, 24)]
        let truth = try rings.map { try closedForm($0.n, $0.m) }
        let strongestProfile = try #require(au.peaks.compactMap(\.profileIntensity).max())
        #expect(strongestProfile == 100 && au.peaks[0].profileIntensity == 100, "111 is the strongest ring in the profile")
        for (i, ring) in rings.enumerated() {
            let p = au.peaks[i]
            #expect(p.multiplicity == ring.m)
            #expect(abs(p.d - a / Double(ring.n).squareRoot()) <= 1e-12 * p.d)
            #expect(abs(p.absoluteIntensity - truth[i].sum) <= 1e-9 * truth[i].sum, "Σ|F|² of ring \(i)")
            let profile = try #require(p.profileIntensity)
            let expected = truth[i].profile / truth[0].profile * 100
            #expect(abs(profile - expected) <= 1e-9 * expected, "profile of ring \(i): \(profile) against \(expected)")
            #expect(abs(p.intensity - truth[i].sum / truth[0].sum * 100) <= 1e-9 * p.intensity)
        }
        // The numbers of the issue, to four figures.
        let shown = au.peaks.prefix(4).map { (($0.profileIntensity ?? 0) * 100).rounded() / 100 }
        #expect(shown == [100, 46.26, 25.66, 26.93], "\(shown)")
        #expect(au.peaks.allSatisfy { ($0.profileIntensity ?? -1) >= 0 && ($0.profileIntensity ?? 101) <= 100 })
        // X-ray peaks carry their Lorentz-polarisation factor already: no profile intensity.
        #expect(try PowderPattern.xray(structure: CrystalPresets.aluminum).peaks.allSatisfy { $0.profileIntensity == nil })
        // The marks floor electron rings by the profile: against Al, Au 311 (26.93) stays above a floor of 20 and Au 222 (7.49)
        // falls under it, where Σ|F|² kept 222 (29.94).
        let al = try PowderPattern.electron(structure: CrystalPresets.aluminum, wavelengthAngstroms: Self.lambda200kV, kMax: 1.2)
        let marks = PowderPattern.marks([al, au], resolution: .fraction(0.01), floor: 20)
        #expect(marks[0][3].isOverlapped, "Al 311 against Au 311 at 26.93")
        #expect(!marks[0][4].isOverlapped, "Al 222 against Au 222 at 7.49: under the floor in the profile")
    }

    /// The profile is Σ I exp(−(x − p)² / 2σ²) / (σ√2π) with σ = FWHM / (2√(2 ln 2)), at 1e-12 relative; a lone peak falls to
    /// half at ± FWHM / 2 and its area is its intensity.
    @Test func profileIsTheClosedForm() throws {
        let p = try PowderPattern.xray(structure: CrystalPresets.aluminum)
        let width = 0.1, step = 0.01
        let y = p.profile(width: width, range: 20...90, step: step)
        #expect(y.count == 7001)
        let sigma = width / (2 * (2 * log(2.0)).squareRoot())
        for i in stride(from: 0, to: y.count, by: 37) {
            let x = 20 + Double(i) * step
            let expected = p.peaks.reduce(0.0) { s, k in
                s + k.intensity * exp(-(x - k.position) * (x - k.position) / (2 * sigma * sigma)) / (sigma * (2 * Double.pi).squareRoot())
            }
            #expect(abs(y[i] - expected) <= 1e-12 * max(1, abs(expected)), "at \(x)")
        }
        let first = p.peaks[0]
        let around = p.profile(width: width, range: (first.position - width / 2)...(first.position + width / 2), step: width / 2)
        #expect(around.count == 3)
        #expect(abs(around[0] / around[1] - 0.5) < 1e-9 && abs(around[2] / around[1] - 0.5) < 1e-9)
        let fine = p.profile(width: width, range: (first.position - 1)...(first.position + 1), step: 0.001)
        #expect(abs(fine.reduce(0, +) * 0.001 - first.intensity) < 1e-6 * first.intensity)
        #expect(p.profile(width: 0, range: 0...1, step: 0.1).isEmpty && p.profile(width: 1, range: 0...1, step: -1).isEmpty)
        #expect(p.profile(width: 1, range: 0...1, step: 1e-300).isEmpty && p.profile(width: 1, range: -1e308...1e308, step: 1).isEmpty)
    }

    /// Al + Mg2Si + Si, Cu Kα, 0 to 90°, resolution 0.1°: Al 1 of 5 peaks overlapped (78.31°), diagnostic 111 at 38.51°;
    /// Mg2Si 1 of 11 (78.21°), diagnostic 40.16°; Si none of 6, diagnostic 28.47° (pymatgen's peak lists with the same rule).
    @Test func diagnosticPeaksOfAlMg2SiSi() throws {
        let patterns = try [CrystalPresets.aluminum, Self.mg2Si, CrystalPresets.silicon].map { try PowderPattern.xray(structure: $0) }
        let marks = PowderPattern.marks(patterns, resolution: .absolute(0.1))
        let expected: [(count: Int, overlapped: [Double], diagnostic: Double)] = [
            (5, [78.3069], 38.50570523371996), (11, [78.2092], 40.1598145357157), (6, [], 28.466011929742603),
        ]
        for (i, e) in expected.enumerated() {
            #expect(patterns[i].peaks.count == e.count && marks[i].count == e.count)
            let over = marks[i].indices.filter { marks[i][$0].isOverlapped }.map { patterns[i].peaks[$0].position }
            #expect(over.count == e.overlapped.count && zip(over, e.overlapped).allSatisfy { abs($0 - $1) < 1e-4 }, "\(over)")
            let diagnostic = marks[i].indices.filter { marks[i][$0].isDiagnostic }
            #expect(diagnostic.count == 1)
            #expect(abs(patterns[i].peaks[diagnostic[0]].position - e.diagnostic) < 1e-7)
        }
        // Al's 78.31° and Mg2Si's 78.21° name each other.
        let al = try #require(marks[0].firstIndex { $0.isOverlapped })
        let mg = try #require(marks[1].firstIndex { $0.isOverlapped })
        #expect(marks[0][al].partners == [PowderPeakReference(pattern: 1, peak: mg)])
        #expect(marks[1][mg].partners == [PowderPeakReference(pattern: 0, peak: al)])
    }

    /// Planted: Al and Au rings sit 0.7 % apart. At 1 % every ring of each is overlapped and neither has a diagnostic ring;
    /// at 0.5 % none is, and each one's diagnostic ring is its strongest. A floor above 100 overlaps nothing; patterns on
    /// different axes never overlap.
    @Test func marksOnPlantedLists() throws {
        let al = try PowderPattern.electron(structure: CrystalPresets.aluminum, wavelengthAngstroms: Self.lambda200kV, kMax: 1.2)
        let au = try PowderPattern.electron(structure: CrystalPresets.gold, wavelengthAngstroms: Self.lambda200kV, kMax: 1.2)
        #expect(al.peaks.count == au.peaks.count)
        let wide = PowderPattern.marks([al, au], resolution: .fraction(0.01))
        #expect(wide.allSatisfy { $0.allSatisfy { $0.isOverlapped && !$0.isDiagnostic } })
        #expect(wide[0].enumerated().allSatisfy { $0.element.partners == [PowderPeakReference(pattern: 1, peak: $0.offset)] })
        let narrow = PowderPattern.marks([al, au], resolution: .fraction(0.005))
        for (m, p) in zip(narrow, [al, au]) {
            #expect(m.allSatisfy { !$0.isOverlapped })
            let d = m.indices.filter { m[$0].isDiagnostic }
            #expect(d.count == 1 && p.peaks[d[0]].intensity == 100)
        }
        #expect(PowderPattern.marks([al, au], resolution: .fraction(0.01), floor: 101).allSatisfy { $0.allSatisfy { !$0.isOverlapped } })
        let xray = try PowderPattern.xray(structure: CrystalPresets.aluminum)
        #expect(PowderPattern.marks([al, xray], resolution: .absolute(100)).allSatisfy { $0.allSatisfy { !$0.isOverlapped } })
        #expect(PowderPattern.marks([], resolution: .absolute(0.1)).isEmpty)
    }
}
