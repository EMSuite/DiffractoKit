//
//  Rotation+CrystalOrientation.swift
//  Role: one rotation value in the kit (DiffractoKit#179). orix's `Rotation` (the quaternion of the moved orix port) and
//        `CrystalOrientation` (DiffractoCore, #138) convert through `CrystalOrientation.orixQuaternion`, the conversion pin
//        P40 holds against orix's `Orientation.from_matrix`: no second matrix or Euler convention.
//

import DiffractoCore

extension Rotation {
    /// orix's rotation of `orientation`: its quaternion `orientation.orixQuaternion` (unit, a ≥ 0), proper.
    public init(_ orientation: CrystalOrientation) {
        self.init(unit: orientation.orixQuaternion)
    }
}

extension CrystalOrientation {
    /// The orientation of a proper orix rotation (its quaternion through `init(orixQuaternion:)`). Throws for an improper
    /// rotation, which no crystal orientation is, and for a quaternion that is not finite.
    public init(_ rotation: Rotation) throws(CrystalOrientationError) {
        guard !rotation.improper else { throw .notARotation }
        try self.init(orixQuaternion: rotation.q)
    }
}
