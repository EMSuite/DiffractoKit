//
//  PowderXRDPinsTests.swift
//  The kit against pymatgen's own XRDCalculator (pin P45):
//    pymatgen-powder-xrd-pins.json  tools/pins/powder-xrd/pins.py: pymatgen 2026.9.24 (pymatgen-core 2026.10.16)
//  27 cases: the 19 committed structures and planted Fe and Mg2Si with Cu Kα at 0 to 90° (pymatgen's defaults), Mo Kα,
//  a 20 to 60° range, 0 to 180°, Debye-Waller factors. Public API only.
//  Tolerances: the table bit for bit; the peak count, d's peak order and the hkl labels with their multiplicities exact;
//  2θ, d, scaled and unscaled intensity at 1e-7 relative (assert_allclose's default: pymatgen ships no tests; its
//  SCALED_INTENSITY_TOL is a peak-drop threshold, not a comparison tolerance).
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoScattering

private struct XPinSite: Decodable { let z: Int; let frac: [Double]; let occupancy: Double }
private struct XPinStructure: Decodable { let name: String; let cell: [Double]; let sites: [XPinSite] }
private struct XPinTableRow: Decodable { let symbol: String; let z: Int; let coefficients: [[Double]] }
private struct XPinFamily: Decodable { let indices: [Int]; let multiplicity: Int }
private struct XPinPeak: Decodable {
    let twoTheta: Double
    let intensity: Double
    let absoluteIntensity: Double
    let d: Double
    let families: [XPinFamily]
}
private struct XPinCase: Decodable {
    let structure: String
    let radiation: String?
    let wavelength: Double
    let twoThetaRange: [Double]
    let debyeWallerFactors: [String: Double]
    let isHexagonalLattice: Bool
    let peakCount: Int
    let peaks: [XPinPeak]
}
private struct XPins: Decodable {
    let wavelengths: [String: Double]
    let structures: [XPinStructure]
    let table: [XPinTableRow]
    let cases: [XPinCase]
}

private enum PowderFixtures {
    static let url = FixtureDirs.fixtures(of: "DiffractoScatteringTests").appendingPathComponent("pymatgen-powder-xrd-pins.json")
    static let result: Result<XPins, any Error> = Result {
        let d = JSONDecoder()
        d.keyDecodingStrategy = .convertFromSnakeCase
        return try d.decode(XPins.self, from: Data(contentsOf: url))
    }
    static func pins() throws -> XPins { try result.get() }

    static func structure(_ name: String) throws -> CrystalStructure {
        let p = try #require(try pins().structures.first { $0.name == name }, "no pinned structure \(name)")
        let c = p.cell
        let cell = UnitCell(a: c[0], b: c[1], c: c[2], alphaDeg: c[3], betaDeg: c[4], gammaDeg: c[5])
        let sites = p.sites.map { AtomSite(z: $0.z, fractional: SIMD3($0.frac[0], $0.frac[1], $0.frac[2]), occupancy: $0.occupancy) }
        return CrystalStructure(cell: cell, sites: sites, operations: [], name: name)
    }
}

private func relative(_ a: Double, _ b: Double) -> Double { abs(a - b) / max(abs(b), .leastNormalMagnitude) }

private func sortedFamilies(_ f: [[Int]: Int]) -> [String] {
    f.map { "\($0.key):\($0.value)" }.sorted()
}

struct PowderXRDPinsTests {

    /// pymatgen's table: every symbol's four (a, b) pairs equal the kit's row Z bit for bit (D is Z = 1, equal to H).
    /// RED UNDER: a coefficient retyped or rounded; a row shifted by one Z.
    @Test func tableBitForBit() throws {
        let pins = try PowderFixtures.pins()
        #expect(pins.table.count == 99)
        #expect(XRayScatteringFactors.parameters.count == XRayScatteringFactors.maximumAtomicNumber)
        var bad: [String] = []
        for row in pins.table {
            let want = row.coefficients.flatMap { $0 }
            guard row.z >= 1, row.z <= XRayScatteringFactors.parameters.count else { bad.append("\(row.symbol) Z \(row.z)"); continue }
            let got = XRayScatteringFactors.parameters[row.z - 1]
            if want.count != 8 || got.count != 8 { bad.append("\(row.symbol) width \(got.count)"); continue }
            for j in 0..<8 where want[j].bitPattern != got[j].bitPattern { bad.append("\(row.symbol) p\(j) \(got[j]) vs \(want[j])") }
        }
        #expect(bad.isEmpty, "table differs: \(bad.prefix(5))")
    }

    /// pymatgen's named wavelengths, every one equal and in its order.
    @Test func wavelengthsEqualPymatgens() throws {
        let pins = try PowderFixtures.pins()
        #expect(XRayRadiation.pymatgenLines.count == pins.wavelengths.count)
        for line in XRayRadiation.pymatgenLines {
            let name = try #require(line.name)
            #expect(pins.wavelengths[name] == line.wavelengthAngstroms, "\(name)")
        }
        #expect(XRayRadiation.cuKAlpha.wavelengthAngstroms == 1.54184)
    }

    /// Each case: the peak count, then per peak 2θ, d, scaled and unscaled intensity within 1e-7 relative and the hkl
    /// labels with their multiplicities exactly (as a set: pymatgen lists them in the order the members first appear).
    /// RED UNDER: the Lorentz-polarisation factor or the factor 2 dropped; s = |g| for |g| / 2; the merge tolerance or the
    /// floor changed; the half-space or the order of the members changed (the first member's 2θ and d); three indices on a
    /// hexagonal lattice.
    @Test(arguments: 0..<27)
    func caseMatchesPymatgen(_ index: Int) throws {
        let pins = try PowderFixtures.pins()
        #expect(pins.cases.count == 27)
        guard index < pins.cases.count else { return }
        let c = pins.cases[index]
        let label = "\(c.structure) \(c.radiation ?? "\(c.wavelength)") \(c.twoThetaRange) B \(c.debyeWallerFactors)"
        let structure = try PowderFixtures.structure(c.structure)
        var radiation = XRayRadiation(wavelengthAngstroms: c.wavelength)
        if let name = c.radiation { radiation = try #require(XRayRadiation.named(name), "no line \(name)") }
        #expect(radiation.wavelengthAngstroms == c.wavelength)
        var dw: [Int: Double] = [:]
        for (z, b) in c.debyeWallerFactors { dw[try #require(Int(z))] = b }
        let pattern = try PowderPattern.xray(structure: structure, radiation: radiation,
                                             twoThetaRange: c.twoThetaRange[0]...c.twoThetaRange[1], debyeWallerFactors: dw)
        #expect(pattern.axis == .twoThetaDegrees)
        #expect(pattern.peaks.count == c.peakCount, "\(label): \(pattern.peaks.count) peaks, pymatgen \(c.peakCount)")
        guard pattern.peaks.count == c.peaks.count else { return }
        var worst = 0.0, worstAt = ""
        var labelErrors: [String] = []
        for (got, want) in zip(pattern.peaks, c.peaks) {
            for (name, a, b) in [("2θ", got.position, want.twoTheta), ("d", got.d, want.d),
                                 ("I", got.intensity, want.intensity), ("I_abs", got.absoluteIntensity, want.absoluteIntensity)] {
                let e = relative(a, b)
                if e > worst { worst = e; worstAt = "\(name) at \(want.twoTheta)" }
            }
            var g: [[Int]: Int] = [:], w: [[Int]: Int] = [:]
            for f in got.families { g[f.indices] = f.multiplicity }
            for f in want.families { w[f.indices] = f.multiplicity }
            if sortedFamilies(g) != sortedFamilies(w) {
                labelErrors.append("\(want.twoTheta): \(sortedFamilies(g)) vs \(sortedFamilies(w))")
            }
            let width = c.isHexagonalLattice ? 4 : 3
            if got.families.contains(where: { $0.indices.count != width }) { labelErrors.append("\(want.twoTheta): not \(width) indices") }
        }
        #expect(worst <= 1e-7, "\(label): worst relative difference \(worst) (\(worstAt))")
        #expect(labelErrors.isEmpty, "\(label): \(labelErrors.prefix(3))")
    }

    /// The 2.4 design's numbers: Al with Cu Kα has 5 peaks from 10 to 90°, the strongest 111 at 38.51°, multiplicity 8;
    /// Mg2Si's strongest at 40.16°, Si's at 28.47°.
    @Test func designNumbers() throws {
        let al = try PowderPattern.xray(structure: try PowderFixtures.structure("Al"))
        #expect(al.peaks.count == 5)
        let strongest = try #require(al.strongest)
        #expect(abs(strongest.position - 38.51) < 0.005 && abs(strongest.intensity - 100) <= 1e-12)
        #expect(strongest.families == [PowderPattern.Family(indices: [1, 1, 1], multiplicity: 8)])
        let mg2si = try #require(try PowderPattern.xray(structure: try PowderFixtures.structure("Mg2Si")).strongest)
        #expect(abs(mg2si.position - 40.16) < 0.005)
        let si = try #require(try PowderPattern.xray(structure: try PowderFixtures.structure("Si")).strongest)
        #expect(abs(si.position - 28.47) < 0.005)
    }
}
