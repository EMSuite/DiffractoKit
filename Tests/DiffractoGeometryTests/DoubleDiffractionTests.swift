//
//  DoubleDiffractionTests.swift
//  Role: `DoubleDiffraction` (S3) against the extinction rules, a closed form: within one crystal the added spots are
//        exactly the reflections of the zone that the structure forbids and the window excites. Silicon down [110] adds
//        the {002} and {222} family (h + k + l = 4n + 2, all even: 002, 2-22, 006, ...), down [001] and [111] none; fcc
//        (Al) and bcc (Fe) add none down any zone (their allowed sets are closed under addition); Mg down [2-1-10] adds 0001
//        (l odd with h − k = 3n). Two crystals: every spot between them is the sum of its routes and on no kinematic spot.
//

import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoGeometry

@Suite struct DoubleDiffractionTests {

    private static func structure(_ cell: UnitCell, _ z: Int, _ sites: [SIMD3<Double>]) -> CrystalStructure {
        CrystalStructure(cell: cell, sites: sites.map { AtomSite(z: z, fractional: $0) }, operations: [.identity])
    }
    private static let fcc = [SIMD3(0.0, 0, 0), SIMD3(0, 0.5, 0.5), SIMD3(0.5, 0, 0.5), SIMD3(0.5, 0.5, 0)]
    private let silicon = Self.structure(.cubic(a: 5.4309), 14, Self.fcc + Self.fcc.map { $0 + SIMD3(0.25, 0.25, 0.25) })
    private let aluminium = Self.structure(.cubic(a: 4.0495), 13, Self.fcc)
    private let iron = Self.structure(.cubic(a: 2.8665), 26, [SIMD3(0, 0, 0), SIMD3(0.5, 0.5, 0.5)])
    private let magnesium = Self.structure(.hexagonal(a: 3.2094, c: 5.2108), 12, [SIMD3(1.0 / 3, 2.0 / 3, 0.25), SIMD3(2.0 / 3, 1.0 / 3, 0.75)])
    private let volts = 200_000.0

    /// The reflections of zone `uvw` within k_max that `forbidden` names and the window (|s| ≤ 0.06) excites.
    private func expected(_ s: CrystalStructure, _ uvw: SIMD3<Int>, _ o: CrystalOrientation, kMax: Double,
                          forbidden: (MillerIndex) -> Bool) -> Set<MillerIndex> {
        let lambda = ElectronWavelength.angstroms(volts: volts)
        var out = Set<MillerIndex>()
        let n = Int((kMax * max(s.cell.a, s.cell.c)).rounded(.up)) + 1
        for h in -n...n { for k in -n...n { for l in -n...n {
            let hkl = MillerIndex(h, k, l)
            guard !hkl.isZero, h * uvw.x + k * uvw.y + l * uvw.z == 0, forbidden(hkl) else { continue }
            let g = s.cell.gVector(hkl)
            guard simd_length(g) <= kMax else { continue }
            let p = o.labCoordinates(of: g)
            let sg = (2 * p.z - lambda * simd_length_squared(p)) / (2 - 2 * lambda * p.z)
            if abs(sg) <= 0.06 { out.insert(hkl) }
        } } }
        return out
    }

    private func added(_ s: CrystalStructure, _ uvw: SIMD3<Int>, kMax: Double = 2) throws -> (Set<MillerIndex>, CrystalOrientation, DoubleDiffractionPatterns) {
        let o = try CrystalOrientation(zoneAxis: uvw, inPlaneAngle: 0.3, cell: s.cell)
        let r = try DoubleDiffraction.spots(simulator: try SpotPatternSimulator(structure: s, kMax: kMax), orientation: o, volts: volts)
        for spot in r.added {
            #expect(spot.crystal == 0)
            for route in spot.routes {
                let sum = MillerIndex(route.first.hkl.h + route.second.hkl.h, route.first.hkl.k + route.second.hkl.k,
                                      route.first.hkl.l + route.second.hkl.l)
                #expect(sum == spot.hkl, "\(route) does not add to \(String(describing: spot.hkl))")
                #expect(route.weight > 0)
            }
            #expect(zip(spot.routes, spot.routes.dropFirst()).allSatisfy { $0.weight >= $1.weight })
        }
        return (Set(r.added.compactMap(\.hkl)), o, r)
    }

    @Test func siliconDown110AddsThe002And222Family() throws {
        let (got, o, r) = try added(silicon, SIMD3(1, 1, 0))
        // Diamond: all even with h + k + l = 4n + 2 is forbidden (fcc and the second sublattice cancel).
        let want = expected(silicon, SIMD3(1, 1, 0), o, kMax: 2) {
            $0.h % 2 == 0 && $0.k % 2 == 0 && $0.l % 2 == 0 && (($0.h + $0.k + $0.l) % 4 + 4) % 4 == 2
        }
        #expect(got == want)
        #expect(got.contains(MillerIndex(0, 0, 2)) && got.contains(MillerIndex(0, 0, -2)))
        #expect(got.contains(MillerIndex(2, -2, 2)) && got.contains(MillerIndex(-2, 2, -2)))
        // 002 by 1-11 then -111 (both {111}, both spots of the pattern).
        let s002 = try #require(r.added.first { $0.hkl == MillerIndex(0, 0, 2) })
        #expect(s002.routes.contains { Set([$0.first.hkl, $0.second.hkl]) == Set([MillerIndex(1, -1, 1), MillerIndex(-1, 1, 1)]) })
        let kinematic = Set(r.patterns[0].spots.map(\.hkl))
        #expect(got.isDisjoint(with: kinematic))
    }

    @Test(arguments: [SIMD3(0, 0, 1), SIMD3(1, 1, 1)])
    func siliconDown001And111AddsNone(_ uvw: SIMD3<Int>) throws {
        #expect(try added(silicon, uvw).0.isEmpty)
    }

    @Test(arguments: [SIMD3(0, 0, 1), SIMD3(0, 1, 1), SIMD3(1, 1, 1), SIMD3(1, 1, 2), SIMD3(1, 2, 3), SIMD3(0, 1, 2)])
    func fccAndBccAddNone(_ uvw: SIMD3<Int>) throws {
        #expect(try added(aluminium, uvw).0.isEmpty)
        #expect(try added(iron, uvw).0.isEmpty)
    }

    @Test func magnesiumDown2110Adds0001() throws {
        // [2 -1 -1 0] in four indices is [1 0 0] in three (u = 2U + V, v = 2V + U, w = W).
        let (got, o, _) = try added(magnesium, SIMD3(1, 0, 0), kMax: 1.2)
        let want = expected(magnesium, SIMD3(1, 0, 0), o, kMax: 1.2) { ($0.l % 2 != 0) && (($0.h - $0.k) % 3 == 0) }
        #expect(got == want)
        #expect(got.contains(MillerIndex(0, 0, 1)) && got.contains(MillerIndex(0, 0, -1)))
    }

    @Test func spotsBetweenTwoCrystalsAreSumsOfTheirRoutes() throws {
        let simulator = try SpotPatternSimulator(structure: aluminium, kMax: 1.5)
        let a = try CrystalOrientation(zoneAxis: SIMD3(0, 0, 1), cell: aluminium.cell)
        let b = a.rotatedInPlane(by: .pi / 4)
        let crystals = [DiffractingCrystal(simulator: simulator, orientation: a), DiffractingCrystal(simulator: simulator, orientation: b)]
        let r = try DoubleDiffraction.spots(crystals: crystals, volts: volts)
        #expect(r.patterns.count == 2)
        let between = r.added.filter { $0.crystal == nil }
        #expect(!between.isEmpty)
        let kinematic = r.patterns.flatMap { $0.spots.map(\.position) }
        let position: [Int: [MillerIndex: SIMD2<Double>]] = [0: Dictionary(uniqueKeysWithValues: r.patterns[0].spots.map { ($0.hkl, $0.position) }),
                                                              1: Dictionary(uniqueKeysWithValues: r.patterns[1].spots.map { ($0.hkl, $0.position) })]
        for spot in between {
            #expect(spot.hkl == nil)
            // Both orders are routes.
            #expect(spot.routes.contains { $0.first.crystal == 1 } && spot.routes.contains { $0.first.crystal == 0 })
            #expect(kinematic.allSatisfy { simd_distance($0, spot.position) > 1e-7 })
            #expect(abs(spot.excitationError) <= 0.06)
            for route in spot.routes {
                let p1 = try #require(position[route.first.crystal]?[route.first.hkl])
                let p2 = try #require(position[route.second.crystal]?[route.second.hkl])
                #expect(simd_distance(p1 + p2, spot.position) <= 2e-6)
            }
        }
        // One crystal alone is the single-crystal call.
        let one = try DoubleDiffraction.spots(crystals: [crystals[0]], volts: volts)
        #expect(one == (try DoubleDiffraction.spots(simulator: simulator, orientation: a, volts: volts)))
        #expect(one.added.isEmpty)
        #expect(throws: SpotPatternError.self) { try DoubleDiffraction.spots(crystals: [], volts: volts) }
    }
}
