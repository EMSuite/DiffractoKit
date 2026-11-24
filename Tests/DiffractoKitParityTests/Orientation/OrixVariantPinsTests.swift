import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoOrientation

/// Pins: tools/pins/orix/variant-pins.py (orix 0.15.0, P52). orix `Misorientation(R, symmetry=(Gp, Gm)).equivalent()` for the
/// planted β″, θ′ and T1 relations in Al. Tolerance atol 1e-4 (orix test_misorientation.py) on the matrices; the rotation of
/// the two pairs and the misfit at 1e-9 (the same construction in numpy; measured at round-off); the class count exact.
private struct VariantPins: Decodable {
    struct Case: Decodable {
        let name, precipitateSymmetry: String
        let R: [[Double]]
        let equivalents: [[[Double]]]
        let variantCount: Int
        let secondPairMisfit: Double
    }
    let orix: String
    let cases: [Case]
}

private func mat(_ r: [[Double]]) -> simd_double3x3 { .init(rows: r.map { SIMD3($0[0], $0[1], $0[2]) }) }
private func gap(_ a: simd_double3x3, _ b: simd_double3x3) -> Double {
    var g = 0.0
    for c in 0 ..< 3 { for r in 0 ..< 3 { g = max(g, abs(a[c][r] - b[c][r])) } }
    return g
}

@Suite struct OrixVariantPinsTests {
    static let al = CrystalPresets.aluminum

    static func setup(_ name: String) -> (CrystalStructure, OrientationRelationship, OrientationRelationship)? {
        switch name {
        case "beta-double-prime":
            return (CrystalPresets.betaDoublePrime,
                    OrientationRelationship(candidate: .direction(SIMD3(1, 0, 0)), matrix: .direction(SIMD3(3, 2, 0))),
                    OrientationRelationship(candidate: .direction(SIMD3(0, 0, 1)), matrix: .direction(SIMD3(1, -3, 0))))
        case "theta-prime":
            return (PrecipitateFixtures.thetaPrime,
                    OrientationRelationship(candidate: .plane(SIMD3(0, 0, 1)), matrix: .plane(SIMD3(0, 0, 1))),
                    OrientationRelationship(candidate: .direction(SIMD3(1, 0, 0)), matrix: .direction(SIMD3(1, 0, 0))))
        case "t1":
            return (PrecipitateFixtures.t1,
                    OrientationRelationship(candidate: .plane(SIMD3(0, 0, 1)), matrix: .plane(SIMD3(1, 1, 1))),
                    OrientationRelationship(candidate: .direction(SIMD3(2, 1, 0)), matrix: .direction(SIMD3(1, -1, 0))))
        default: return nil
        }
    }

    /// Mutation: the representatives taken as R g (the operator on the wrong side) leave orix's set Gm R Gp.
    @Test func variantsAreOrixClassesOfTheEquivalents() throws {
        let file = FixtureDirs.target("DiffractoOrientationTests").appendingPathComponent("Fixtures/orix-variant-pins.json")
        let pins = try JSONDecoder().decode(VariantPins.self, from: Data(contentsOf: file))
        #expect(pins.orix == "0.15.0")
        #expect(pins.cases.map(\.name) == ["beta-double-prime", "theta-prime", "t1"])
        let m = try CrystalOrientation(zoneAxis: SIMD3(0, 0, 1), inPlaneAngle: 0.2, cell: Self.al.cell)
        for c in pins.cases {
            let (precipitate, first, second) = try #require(Self.setup(c.name))
            let r = try OrientationRelationship.rotation(first, second, precipitate: precipitate, matrix: Self.al)
            #expect(gap(r, mat(c.R)) < 1e-9, "\(c.name)")
            let vs = try OrientationRelationship.variantOrientations(first: first, second: second, precipitate: precipitate,
                                                                     matrix: Self.al, matrixOrientation: m, volts: 200_000, kMax: 1)
            #expect(vs.count == c.variantCount, "\(c.name)")
            let equivalents = c.equivalents.map(mat)
            for v in vs {
                #expect(equivalents.contains { gap($0, v.relation) < 1e-4 }, "\(c.name) variant \(v.index) is not an orix equivalent")
            }
            let ps = try #require(PhaseModel.orientationSymmetry(of: precipitate))
            for e in equivalents {
                #expect(vs.contains { v in ps.operators.contains { gap(v.relation * $0, e) < 1e-4 } },
                        "\(c.name): an orix equivalent in no variant's class")
            }
            let misfit = try OrientationRelationship.secondPairMisfit(first: first, second: second, precipitate: precipitate, matrix: Self.al)
            #expect(abs(misfit - c.secondPairMisfit) < 1e-9, "\(c.name)")
        }
    }
}
