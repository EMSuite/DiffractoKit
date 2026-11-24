//
//  PairDoubleDiffractionTests.swift
//  DiffractoKit S8 (#163): double diffraction between a precipitate and its matrix on the drawn spot patterns. No Python
//  reference: a closed rule, checked against a brute-force sum of the two spot lists, and its counts for β″, θ′ and T1 on
//  [001], [110] and [112] Al against the ones recorded from py4DSTEM's spot patterns (tools/pins/double-diffraction/pins.py).
//

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoGeometry
import DiffractoKitTestSupport

private struct Pins: Decodable {
    struct Case: Decodable {
        let precipitate: String
        let zone: [Int]
        let matrixSpots, precipitateSpots, positions, nearCut: Int
        let routes: [[Double]]
        let matrix: [Double]
    }
    let py4DSTEM: String
    let volts, kMax, cut: Double
    let cases: [Case]
}

@Suite struct PairDoubleDiffractionTests {
    static let al = CrystalPresets.aluminum

    /// Planted θ′ and T1 cells (lattices published, sites for the Laue class), as tools/pins/variant-pins.py.
    static var thetaPrime: CrystalStructure {
        CrystalPresets.structure(a: 4.04, b: 4.04, c: 5.80, sites: [
            AtomSite(z: 29, fractional: [0, 0, 0]), AtomSite(z: 29, fractional: [0.5, 0.5, 0.5]),
            AtomSite(z: 13, fractional: [0, 0.5, 0.25]), AtomSite(z: 13, fractional: [0.5, 0, 0.25]),
            AtomSite(z: 13, fractional: [0, 0.5, 0.75]), AtomSite(z: 13, fractional: [0.5, 0, 0.75]),
        ])
    }
    static var t1: CrystalStructure {
        CrystalPresets.structure(a: 4.965, b: 4.965, c: 9.345, gammaDeg: 120, sites: [
            AtomSite(z: 29, fractional: [0, 0, 0]),
            AtomSite(z: 13, fractional: [1.0 / 3.0, 2.0 / 3.0, 0.5]), AtomSite(z: 13, fractional: [2.0 / 3.0, 1.0 / 3.0, 0.5]),
            AtomSite(z: 3, fractional: [0, 0, 0.5]),
        ])
    }

    enum Kind { case plane, direction }
    static func cartesian(_ kind: Kind, _ i: SIMD3<Double>, in s: CrystalStructure) -> SIMD3<Double> {
        let rows = kind == .direction ? s.cell.realLattice : s.cell.reciprocalLattice
        return i.x * rows[0] + i.y * rows[1] + i.z * rows[2]
    }
    static func triad(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> simd_double3x3 {
        let e1 = simd_normalize(a), e2 = simd_normalize(b - simd_dot(b, e1) * e1)
        return simd_double3x3(columns: (e1, e2, simd_cross(e1, e2)))
    }

    /// The precipitate and its first-variant relation (precipitate Cartesian to Al Cartesian), the pins' pairs.
    static func setup(_ name: String) throws -> (CrystalStructure, simd_double3x3) {
        let p: CrystalStructure
        let pairs: [(Kind, SIMD3<Double>, Kind, SIMD3<Double>)]
        switch name {
        case "beta''":
            p = CrystalPresets.betaDoublePrime
            pairs = [(.direction, [1, 0, 0], .direction, [3, 2, 0]), (.direction, [0, 0, 1], .direction, [1, -3, 0])]
        case "theta'":
            p = thetaPrime
            pairs = [(.plane, [0, 0, 1], .plane, [0, 0, 1]), (.direction, [1, 0, 0], .direction, [1, 0, 0])]
        case "T1":
            p = t1
            pairs = [(.plane, [0, 0, 1], .plane, [1, 1, 1]), (.direction, [2, 1, 0], .direction, [1, -1, 0])]
        default:
            throw CocoaError(.featureUnsupported)
        }
        let r = triad(cartesian(pairs[0].2, pairs[0].3, in: al), cartesian(pairs[1].2, pairs[1].3, in: al))
            * triad(cartesian(pairs[0].0, pairs[0].1, in: p), cartesian(pairs[1].0, pairs[1].1, in: p)).transpose
        return (p, r)
    }

    static func patterns(_ name: String, zone: SIMD3<Int>) throws -> (SpotPattern, SpotPattern) {
        let (p, r) = try setup(name)
        let m = try CrystalOrientation(zoneAxis: zone, inPlaneAngle: 0, cell: al.cell)
        let a = try SpotPatternSimulator(structure: al, kMax: 1.5).pattern(orientation: m, volts: 200_000)
        let b = try SpotPatternSimulator(structure: p, kMax: 1.5)
            .pattern(orientation: try CrystalOrientation(matrix: r.transpose * m.matrix), volts: 200_000)
        return (a, b)
    }

    struct Key: Hashable { let a: MillerIndex; let b: MillerIndex }

    /// The positions are the brute-force sums: every pair above the cut that misses the direct beam and every kinematic spot
    /// is a route of exactly one position, at its sum; each position's routes are ordered by weight, positions by their
    /// strongest route. Mutation: weights taken from absolute intensities (not relative to each pattern's strongest) change
    /// which routes pass the cut.
    @Test func positionsAreTheBruteForceSums() throws {
        for name in ["beta''", "theta'", "T1"] {
            for zone in [SIMD3(0, 0, 1), SIMD3(1, 1, 0), SIMD3(1, 1, 2)] {
                let (a, b) = try Self.patterns(name, zone: zone)
                let added = a.doubleDiffraction(with: b)
                let maxA = try #require(a.spots.map(\.intensity).max()), maxB = try #require(b.spots.map(\.intensity).max())
                let kinematic = (a.spots + b.spots).map(\.position)
                var expected: [Key: (SIMD2<Double>, Double)] = [:]
                for x in a.spots { for y in b.spots {
                    let w = (x.intensity / maxA) * (y.intensity / maxB)
                    let p = x.position + y.position
                    guard w >= 0.1, simd_length(p) > 1e-6, !kinematic.contains(where: { simd_distance($0, p) <= 1e-6 }) else { continue }
                    expected[Key(a: x.hkl, b: y.hkl)] = (p, w)
                } }
                var seen = Set<Key>()
                for spot in added {
                    #expect(spot.routes.map(\.weight) == spot.routes.map(\.weight).sorted(by: >))
                    #expect(spot.routes.contains { r in expected[Key(a: r.first.hkl, b: r.second.hkl)]?.0 == spot.position })
                    for r in spot.routes {
                        #expect(r.first.crystal == 0 && r.second.crystal == 1)
                        let key = Key(a: r.first.hkl, b: r.second.hkl)
                        let e = try #require(expected[key], "\(name) \(zone): a route the brute force does not have")
                        #expect(r.weight == e.1 && simd_distance(e.0, spot.position) <= 1e-6)
                        #expect(seen.insert(key).inserted, "a route on two positions")
                    }
                }
                #expect(seen == Set(expected.keys), "\(name) \(zone)")
                let top = added.map { $0.routes[0].weight }
                #expect(top == top.sorted(by: >))
                for i in added.indices { for j in added.indices where j > i {
                    #expect(simd_distance(added[i].position, added[j].position) > 1e-6)
                } }
            }
        }
    }

    /// The counts and routes for β″, θ′ and T1 on [001], [110] and [112] Al equal the recorded ones (py4DSTEM spot
    /// patterns, P55); no recorded route lies within 1e-3 of the cut, so the sets compare exactly, weights at 1e-4.
    @Test func countsEqualTheRecordedOnes() throws {
        let file = FixtureDirs.target("DiffractoGeometryTests").appendingPathComponent("Fixtures/double-diffraction-pins.json")
        let pins = try JSONDecoder().decode(Pins.self, from: Data(contentsOf: file))
        #expect(pins.py4DSTEM == "0.14.19" && pins.cut == 0.1 && pins.cases.count == 9)
        for c in pins.cases {
            let zone = SIMD3(c.zone[0], c.zone[1], c.zone[2])
            let (a, b) = try Self.patterns(c.precipitate, zone: zone)
            let m = try CrystalOrientation(zoneAxis: zone, inPlaneAngle: 0, cell: Self.al.cell).matrix
            for col in 0 ..< 3 { for row in 0 ..< 3 { #expect(abs(m[col][row] - c.matrix[row * 3 + col]) < 1e-9) } }
            #expect(a.spots.count == c.matrixSpots && b.spots.count == c.precipitateSpots, "\(c.precipitate) \(c.zone)")
            let added = a.doubleDiffraction(with: b)
            #expect(added.count == c.positions, "\(c.precipitate) \(c.zone): \(added.count) positions, recorded \(c.positions)")
            #expect(c.nearCut == 0)
            var recorded: [Key: Double] = [:]
            for r in c.routes {
                recorded[Key(a: MillerIndex(Int(r[0]), Int(r[1]), Int(r[2])), b: MillerIndex(Int(r[3]), Int(r[4]), Int(r[5])))] = r[6]
            }
            let routes = added.flatMap(\.routes)
            #expect(Set(routes.map { Key(a: $0.first.hkl, b: $0.second.hkl) }) == Set(recorded.keys), "\(c.precipitate) \(c.zone)")
            for r in routes {
                if let w = recorded[Key(a: r.first.hkl, b: r.second.hkl)] { #expect(abs(r.weight - w) <= 1e-4 * w) }
            }
        }
    }

    /// A cut of 0 keeps every pair off the kinematic spots; a cut above 1 keeps none; an invalid cut gives nothing.
    @Test func theCut() throws {
        let (a, b) = try Self.patterns("theta'", zone: SIMD3(0, 0, 1))
        let all = a.doubleDiffraction(with: b, routeWeightCut: 0)
        #expect(all.count >= a.doubleDiffraction(with: b).count && !all.isEmpty)
        #expect(all.allSatisfy { $0.routes.allSatisfy { $0.weight >= 0 } })
        #expect(a.doubleDiffraction(with: b, routeWeightCut: 1.01).isEmpty)
        #expect(a.doubleDiffraction(with: b, routeWeightCut: .nan).isEmpty && a.doubleDiffraction(with: b, routeWeightCut: -1).isEmpty)
    }
}
