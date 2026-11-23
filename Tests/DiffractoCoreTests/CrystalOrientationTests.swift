//
//  CrystalOrientationTests.swift
//  Role: `CrystalOrientation` (S1) on closed forms: the zone axis is the third column, turns add, a tilt moves every pattern
//        vector by the lab rotation, the nearest zone axis and its angle, the refusals, the written form and the quaternion
//        round trip. The py4DSTEM and orix pins are in DiffractoKitParityTests/Core/CrystalOrientationPinsTests.swift.
//

import Foundation
import Testing
import simd
import DiffractoCore

@Suite struct CrystalOrientationTests {

    private let cubic = UnitCell.cubic(a: 4.0495)
    private let hexagonal = UnitCell.hexagonal(a: 3.209, c: 5.211)
    private let monoclinic = UnitCell.monoclinic(a: 15.16, b: 4.05, c: 6.74, betaDeg: 105.3)

    private func gap(_ a: simd_double3x3, _ b: simd_double3x3) -> Double {
        var w = 0.0
        for c in 0..<3 { for r in 0..<3 { w = max(w, abs(a[c][r] - b[c][r])) } }
        return w
    }

    @Test(arguments: [SIMD3(0, 0, 1), SIMD3(1, 1, 0), SIMD3(1, 1, 1), SIMD3(-2, 1, 3), SIMD3(0, 1, 0)])
    func theThirdColumnIsTheZoneAxis(_ uvw: SIMD3<Int>) throws {
        for cell in [cubic, hexagonal, monoclinic] {
            let o = try CrystalOrientation(zoneAxis: uvw, inPlaneAngle: 0.4, cell: cell)
            let d = Double(uvw.x) * cell.realLattice[0] + Double(uvw.y) * cell.realLattice[1] + Double(uvw.z) * cell.realLattice[2]
            #expect(simd_length(o.beamDirection - simd_normalize(d)) <= 1e-14)
            // The zone axis is excited: g·[uvw] = 0 puts g in the pattern plane (lab z = 0).
            let g = cell.gVector(MillerIndex(uvw.y, -uvw.x, 0))
            #expect(abs(o.labCoordinates(of: g).z) <= 1e-14)
            let n = try #require(o.nearestZoneAxis(cell: cell))
            #expect(n.uvw == uvw && n.angle <= 1e-7, "\(n)")
        }
    }

    @Test func inPlaneTurnsAddAndTurnThePatternCounterClockwise() throws {
        let a = try CrystalOrientation(zoneAxis: SIMD3(1, 1, 2), inPlaneAngle: 0.3, cell: hexagonal)
        let b = try CrystalOrientation(zoneAxis: SIMD3(1, 1, 2), inPlaneAngle: 1.1, cell: hexagonal)
        #expect(gap(a.rotatedInPlane(by: 0.8).matrix, b.matrix) <= 1e-14)
        let g = hexagonal.gVector(MillerIndex(1, -1, 0))
        let p = a.labCoordinates(of: g), q = b.labCoordinates(of: g)
        let turn = atan2(p.x * q.y - p.y * q.x, p.x * q.x + p.y * q.y)
        #expect(abs(turn - 0.8) <= 1e-13)
        #expect(abs(simd_length(SIMD2(p.x, p.y)) - simd_length(SIMD2(q.x, q.y))) <= 1e-15)
    }

    /// A tilt R acts on the pattern as p → R p, and the tilt away from a zone axis is the angle to it.
    @Test func aTiltMovesThePatternByTheLabRotationAndAwayFromTheZone() throws {
        let o = try CrystalOrientation(zoneAxis: SIMD3(0, 0, 1), inPlaneAngle: 0.2, cell: cubic)
        let angles = SIMD2(0.02, -0.035)
        let t = o.tilted(by: angles)
        let rx = simd_double3x3(simd_quatd(angle: angles.x, axis: SIMD3(1, 0, 0)))
        let ry = simd_double3x3(simd_quatd(angle: angles.y, axis: SIMD3(0, 1, 0)))
        for hkl in [MillerIndex(2, 0, 0), MillerIndex(1, 1, 1), MillerIndex(-3, 1, 1)] {
            let g = cubic.gVector(hkl)
            #expect(simd_length(t.labCoordinates(of: g) - ry * rx * o.labCoordinates(of: g)) <= 1e-14)
        }
        let n = try #require(t.nearestZoneAxis(cell: cubic))
        #expect(n.uvw == SIMD3(0, 0, 1))
        let expected = acos((ry * rx * SIMD3(0, 0, 1)).z)
        #expect(abs(n.angle - expected) <= 1e-12)
        // About x alone: the angle is the tilt.
        let one = try #require(o.tilted(by: SIMD2(0.0349, 0)).nearestZoneAxis(cell: cubic))
        #expect(abs(one.angle - 0.0349) <= 1e-13)
        // A tilt and its inverse (y first undone, then x) return the orientation.
        let back = try CrystalOrientation(matrix: t.matrix * ry * rx)
        #expect(gap(back.matrix, o.matrix) <= 1e-14)
    }

    @Test func theNearestZoneAxisPointsAtTheSourceAndPrefersLowIndices() throws {
        // Down [0 0 -1] the nearest is [0 0 -1] itself, not [0 0 1].
        let down = try CrystalOrientation(zoneAxis: SIMD3(0, 0, -1), cell: cubic)
        #expect(down.nearestZoneAxis(cell: cubic)?.uvw == SIMD3(0, 0, -1))
        // Halfway between [001] and [0 -1 1] (22.5° from each): equal angles go to the smaller largest index, then the smaller
        // index sum, so [001] (sum 1) and not [0 -1 1] (sum 2).
        let between = try CrystalOrientation(zoneAxis: SIMD3(0, 0, 1), cell: cubic).tilted(by: SIMD2(-Double.pi / 8, 0))
        let n = try #require(between.nearestZoneAxis(cell: cubic, maxIndex: 1))
        #expect(abs(n.angle - .pi / 8) <= 1e-12)
        #expect(n.uvw == SIMD3(0, 0, 1))
        // With indices to 3, [0 -1 2] and [0 -1 3] both lie 4.065° away (atan(1/2) + atan(1/3) = 45°, so the two are equal to
        // rounding): the smaller largest index decides, [0 -1 2].
        let fine = try #require(between.nearestZoneAxis(cell: cubic, maxIndex: 3))
        #expect(abs(fine.angle - (atan2(1.0, 2.0) - .pi / 8)) <= 1e-12)
        #expect(abs(fine.uvw.y) == 1 && fine.uvw.z == 2 && fine.uvw.x == 0)
        #expect(between.nearestZoneAxis(cell: cubic, maxIndex: 0) == nil)
    }

    @Test func refusals() {
        #expect(throws: CrystalOrientationError.zeroZoneAxis) { try CrystalOrientation(zoneAxis: .zero, cell: cubic) }
        #expect(throws: CrystalOrientationError.notFinite) {
            try CrystalOrientation(zoneAxis: SIMD3(0, 0, 1), inPlaneAngle: .nan, cell: cubic)
        }
        #expect(throws: CrystalOrientationError.invalidCell) {
            try CrystalOrientation(zoneAxis: SIMD3(0, 0, 1), cell: UnitCell(a: 4, b: 4, c: -4, alphaDeg: 90, betaDeg: 90, gammaDeg: 90))
        }
        #expect(throws: CrystalOrientationError.notARotation) {
            try CrystalOrientation(matrix: simd_double3x3(diagonal: SIMD3(1, 1, -1)))
        }
        #expect(throws: CrystalOrientationError.notARotation) {
            try CrystalOrientation(matrix: simd_double3x3(diagonal: SIMD3(1, 1, 1.001)))
        }
        #expect(throws: CrystalOrientationError.notARotation) {
            try CrystalOrientation(matrix: simd_double3x3(diagonal: SIMD3(1, .nan, 1)))
        }
        #expect(throws: CrystalOrientationError.notFinite) { try CrystalOrientation(orixQuaternion: .zero) }
        #expect(throws: CrystalOrientationError.notFinite) { try CrystalOrientation(orixQuaternion: SIMD4(1, .infinity, 0, 0)) }
    }

    @Test func theWrittenFormCarriesItsVersionAndRoundTrips() throws {
        let o = try CrystalOrientation(zoneAxis: SIMD3(-2, 1, 3), inPlaneAngle: 2.2, cell: monoclinic)
        let data = try JSONEncoder().encode(o)
        let back = try JSONDecoder().decode(CrystalOrientation.self, from: data)
        #expect(back == o)
        #expect(back.hashValue == o.hashValue)
        let object = try #require(try JSONSerialization.jsonObject(with: data) as? [String: Any])
        #expect(object["schemaVersion"] as? Int == 1)
        let wrongVersion = Data(#"{"schemaVersion": 2, "columns": [[1,0,0],[0,1,0],[0,0,1]]}"#.utf8)
        #expect(throws: DecodingError.self) { try JSONDecoder().decode(CrystalOrientation.self, from: wrongVersion) }
        let notRotation = Data(#"{"schemaVersion": 1, "columns": [[1,0,0],[0,1,0],[0,0,2]]}"#.utf8)
        #expect(throws: DecodingError.self) { try JSONDecoder().decode(CrystalOrientation.self, from: notRotation) }
    }

    @Test func theQuaternionRoundTripsAndIsCanonical() throws {
        for (uvw, angle) in [(SIMD3(1, 2, 3), 0.3), (SIMD3(0, 0, 1), 0.0), (SIMD3(1, -1, 0), 3.0)] {
            let o = try CrystalOrientation(zoneAxis: uvw, inPlaneAngle: angle, cell: monoclinic)
            let q = o.orixQuaternion
            #expect(abs(simd_length(q) - 1) <= 1e-15 && q.x >= 0)
            #expect(gap(try CrystalOrientation(orixQuaternion: q).matrix, o.matrix) <= 1e-14)
            #expect(gap(try CrystalOrientation(orixQuaternion: -2 * q).matrix, o.matrix) <= 1e-14)
        }
        // A half turn about x: (0, 1, 0, 0), not (0, −1, 0, 0).
        let half = try CrystalOrientation(matrix: simd_double3x3(diagonal: SIMD3(1, -1, -1)))
        #expect(simd_length(half.orixQuaternion - SIMD4(0, 1, 0, 0)) <= 1e-15)
    }
}
