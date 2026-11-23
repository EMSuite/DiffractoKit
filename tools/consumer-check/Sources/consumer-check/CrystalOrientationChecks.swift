//
//  CrystalOrientationChecks.swift
//  Role: `CrystalOrientation` (DiffractoCore) and its `ZoneAxis` initialiser (DiffractoGeometry) used from outside the package.
//        Closed forms only; the py4DSTEM and orix pins are the package's own tests.
//

import Foundation
import simd
import DiffractoCore
import DiffractoGeometry

func runCrystalOrientationChecks(_ r: inout Report) throws {
    let cell = UnitCell.cubic(a: 4.0495)
    let o = try CrystalOrientation(zoneAxis: ZoneAxis(0, 0, 1), inPlaneAngle: 0, cell: cell)
    r.check("[001] of a cubic cell is the identity", simd_length(o.matrix.columns.2 - SIMD3(0, 0, 1)) < 1e-15
        && simd_length(o.matrix.columns.0 - SIMD3(1, 0, 0)) < 1e-15)
    let tilted = o.tilted(by: SIMD2(0.05, 0))
    let near = try unwrap(tilted.nearestZoneAxisValue(cell: cell))
    r.equal("the nearest zone axis after a tilt", near.axis, ZoneAxis(0, 0, 1))
    r.near("the angle to it is the tilt", near.angle, 0.05, 1e-13)
    let g = cell.gVector(MillerIndex(2, 0, 0))
    r.near("200 turned by 90 degrees lies on lab y", o.rotatedInPlane(by: .pi / 2).labCoordinates(of: g).y, 2 / 4.0495, 1e-15)
    r.near("the quaternion of the identity", o.orixQuaternion.x, 1, 1e-15)
    r.check("a non-rotation is refused", (try? CrystalOrientation(matrix: simd_double3x3(diagonal: SIMD3(1, 1, -1)))) == nil)
    let data = try JSONEncoder().encode(tilted)
    r.equal("the written form round-trips", try JSONDecoder().decode(CrystalOrientation.self, from: data), tilted)
}
