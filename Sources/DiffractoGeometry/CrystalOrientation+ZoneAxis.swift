//
//  CrystalOrientation+ZoneAxis.swift
//  Role: `CrystalOrientation` (DiffractoCore) with DiffractoGeometry's `ZoneAxis` type: the same calls as the integer triple.
//

import simd
import DiffractoCore

extension CrystalOrientation {
    /// ``CrystalOrientation/init(zoneAxis:inPlaneAngle:cell:)`` with a ``ZoneAxis``.
    public init(zoneAxis: ZoneAxis, inPlaneAngle: Double = 0, cell: UnitCell) throws(CrystalOrientationError) {
        try self.init(zoneAxis: SIMD3(zoneAxis.u, zoneAxis.v, zoneAxis.w), inPlaneAngle: inPlaneAngle, cell: cell)
    }

    /// ``CrystalOrientation/nearestZoneAxis(cell:maxIndex:)`` as a ``ZoneAxis``.
    public func nearestZoneAxisValue(cell: UnitCell, maxIndex: Int = 3) -> (axis: ZoneAxis, angle: Double)? {
        nearestZoneAxis(cell: cell, maxIndex: maxIndex).map { (ZoneAxis($0.uvw.x, $0.uvw.y, $0.uvw.z), $0.angle) }
    }
}
