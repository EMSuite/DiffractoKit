//
//  OrixKitConsistencyTests.swift
//  Role: the moved orix port (DiffractoKit#179) and the kit's own orientation types say the same thing:
//    * one rotation value: `Rotation(CrystalOrientation)` is the orientation's matrix (orix's `qu2om` of
//      `orixQuaternion`, pin P40's conversion), and the round trip gives the orientation back;
//    * one point-group table: for each of the 11 Laue classes, the proper rotations of orix's Laue group are the operators
//      of `CrystalOrientationSymmetry` (as sets, within 1e-12). Two settings differ by convention and are compared in
//      the kit's: 2/m with the twofold along y (the kit's b-unique cell; orix's `C2h` has it along z), and -3m in the
//      setting whose twofolds orix's `D3d` uses;
//    * one colour key: orix's `DirectionColorKeyTSL` of the port and the kit's `.orix` IPF key (P26) give the same colour for
//      every Laue class (within the P26 tolerance, 1e-6), the kit's direction taken into orix's frame by the key's own
//      sector frame (so this holds the folding and the colours in orix's frame; the frames themselves are P26's pins);
//    * one cell: `CrystalLattice(UnitCell)` and `CrystalLattice.unitCell` carry the same six parameters.
//

import DiffractoCore
import Foundation
import Testing
import simd
@testable import DiffractoOrientation

@Suite("The orix port and the kit's orientation types agree")
struct OrixKitConsistencyTests {

    /// Row i, column j of a simd matrix (simd stores columns).
    private static func entry(_ m: simd_double3x3, _ i: Int, _ j: Int) -> Double { m[j][i] }

    private static func rows(_ m: simd_double3x3) -> [Double] { (0..<3).flatMap { i in (0..<3).map { j in entry(m, i, j) } } }

    @Test func aRotationIsTheOrientationsMatrixAndRoundTrips() throws {
        var generator = SplitMix64(seed: 179)
        for n in 0..<200 {
            let q = SIMD4<Double>((0..<4).map { _ in Double.random(in: -1...1, using: &generator) })
            guard simd_length(q) > 0.1 else { continue }
            let o = try CrystalOrientation(orixQuaternion: q)
            let r = Rotation(o)
            let gap = zip(r.matrix.flat, Self.rows(o.matrix)).map { abs($0 - $1) }.max() ?? .nan
            #expect(gap <= 1e-12, "case \(n): Rotation's matrix differs from the orientation's by \(gap)")
            let back = try CrystalOrientation(r)
            let angle = acos(min(1, max(-1, ((o.matrix.transpose * back.matrix)[0][0] + (o.matrix.transpose * back.matrix)[1][1]
                                               + (o.matrix.transpose * back.matrix)[2][2] - 1) / 2)))
            #expect(angle <= 1e-7, "case \(n): round trip off by \(angle) rad")
        }
        #expect(throws: CrystalOrientationError.notARotation) { try CrystalOrientation(Rotation(SIMD4(1, 0, 0, 0), improper: true)) }
    }

    /// The proper rotations of orix's group for a Laue class, in the kit's setting.
    private static func orixProperRotations(_ laue: LaueClass) throws -> [simd_double3x3] {
        let group: PointGroup
        switch laue {
        case .twoOverM: group = PointGroup.fromGenerators([PointGroup.c2y, PointGroup.ci], name: "2/m (y)")
        default: group = try #require(PointGroup.named(laue.symbol), "orix has no group \(laue.symbol)")
        }
        return group.elements.filter { !$0.improper }.map { r in
            let m = r.matrix
            return simd_double3x3(rows: [SIMD3(m[0, 0], m[0, 1], m[0, 2]), SIMD3(m[1, 0], m[1, 1], m[1, 2]),
                                         SIMD3(m[2, 0], m[2, 1], m[2, 2])])
        }
    }

    private static func sameSet(_ a: [simd_double3x3], _ b: [simd_double3x3]) -> Bool {
        func close(_ x: simd_double3x3, _ y: simd_double3x3) -> Bool {
            zip(rows(x), rows(y)).allSatisfy { abs($0 - $1) <= 1e-12 }
        }
        return a.count == b.count && a.allSatisfy { x in b.contains { close(x, $0) } } && b.allSatisfy { y in a.contains { close(y, $0) } }
    }

    @Test(arguments: LaueClass.allCases)
    func theLaueGroupsProperRotationsAreTheKitsOperators(_ laue: LaueClass) throws {
        let orix = try Self.orixProperRotations(laue)
        if laue == .threeBarM {
            let matches = ThreeBarMSetting.allCases.filter {
                Self.sameSet(orix, CrystalOrientationSymmetry(laueClass: laue, setting: $0).operators)
            }
            #expect(matches == [.threeBarMOne], "-3m: orix's D3d matches the kit's settings \(matches)")
        } else {
            let kit = CrystalOrientationSymmetry(laueClass: laue).operators
            #expect(Self.sameSet(orix, kit), "\(laue.symbol): orix \(orix.count) rotations, the kit \(kit.count) operators")
        }
    }

    @Test func aLatticeCarriesTheKitsCell() {
        let cell = UnitCell(a: 3.2094, b: 3.2094, c: 5.2105, alphaDeg: 90, betaDeg: 90, gammaDeg: 120)
        let lattice = CrystalLattice(cell)
        #expect(lattice == CrystalLattice(a: 3.2094, b: 3.2094, c: 5.2105, alpha: 90, beta: 90, gamma: 120))
        let back = lattice.unitCell
        #expect([back.a, back.b, back.c, back.alphaDeg, back.betaDeg, back.gammaDeg] == [cell.a, cell.b, cell.c, 90, 90, 120])
    }

    @Test(arguments: LaueClass.allCases)
    func theTwoTSLKeysGiveTheSameColour(_ laue: LaueClass) throws {
        let symmetry = CrystalOrientationSymmetry(laueClass: laue, setting: .threeBarMOne)
        let group = try #require(PointGroup.named(laue.symbol))
        let key = DirectionColorKeyTSL(group)
        let sector = OrixIPFSector.sector(for: laue)
        let frame = laue == .threeBarM ? matrix_identity_double3x3 : sector.frame
        var generator = SplitMix64(seed: 179)
        var worst = 0.0
        for _ in 0..<500 {
            let d = simd_normalize(SIMD3<Double>((0..<3).map { _ in Double.random(in: -1...1, using: &generator) }))
            let kit = symmetry.orixIPFColor(d)
            let port = key.color(of: [frame.transpose * d])[0]
            worst = max(worst, simd_reduce_max(simd_abs(kit - port)))
        }
        #expect(worst <= 1e-6, "\(laue.symbol): the two keys differ by \(worst)")
    }
}

/// A seeded generator, so a failing draw comes back on the next run.
private struct SplitMix64: RandomNumberGenerator {
    var state: UInt64
    init(seed: UInt64) { state = seed }
    mutating func next() -> UInt64 {
        state &+= 0x9E37_79B9_7F4A_7C15
        var z = state
        z = (z ^ (z >> 30)) &* 0xBF58_476D_1CE4_E5B9
        z = (z ^ (z >> 27)) &* 0x94D0_49BB_1331_11EB
        return z ^ (z >> 31)
    }
}
