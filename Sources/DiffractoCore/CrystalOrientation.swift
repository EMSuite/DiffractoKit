//
//  CrystalOrientation.swift
//  Role: the one orientation value of a crystal against the microscope (DiffractoKit S1): a rotation matrix in py4DSTEM's
//        convention, built from a zone axis and an in-plane angle, tilted and turned in the lab frame, the nearest low-index
//        zone axis and the angle to it, and the conversion to and from orix quaternions. The spot pattern, the zone-axis
//        solver, Index and the crystal view all read this one matrix.
//
//  Convention (py4DSTEM f050d207 `Orientation.matrix`, `Crystal.parse_orientation`, `generate_diffraction_pattern`):
//    * the columns of `matrix` are the lab x, y and z axes written in the crystal's Cartesian frame (the `UnitCell` frame:
//      a along x, b in the xy plane);
//    * lab x and y span the pattern; lab z points from the sample to the electron source, so the beam travels along −z and
//      the third column is the beam direction pointing at the source (the zone axis as py4DSTEM writes it);
//    * a reciprocal vector g (crystal Cartesian, Å⁻¹, no 2π) appears in the pattern at (Mᵀg).x, (Mᵀg).y.
//  orix 0.15.0's `Orientation` holds the same matrix (`to_matrix`), so the quaternion (a, b, c, d) of orix is the quaternion of
//  this matrix (Misorientation.swift names the pins that hold it).
//
//  The frame of a zone axis is py4DSTEM's `parse_orientation` default, NOT the kit's older `ZoneAxisFrame.detectorBasis`
//  (a different in-plane reference): an in-plane angle of 0 here is py4DSTEM's pattern of that zone axis.
//

import Foundation
import simd

/// Why a ``CrystalOrientation`` was refused.
public enum CrystalOrientationError: Error, Sendable, Equatable, CustomStringConvertible {
    /// The matrix is not a proper rotation (orthonormal within 1e-8, determinant +1 within 1e-8) or not finite.
    case notARotation
    /// The zone axis is [0 0 0] or has no length in the cell.
    case zeroZoneAxis
    /// The cell has issues (``UnitCell/issues``) and has no Cartesian frame to write a direction in.
    case invalidCell
    /// An angle or a quaternion component is not finite, or the quaternion has no length.
    case notFinite

    public var description: String {
        switch self {
        case .notARotation: "the matrix is not a proper rotation"
        case .zeroZoneAxis: "the zone axis has no length"
        case .invalidCell: "the cell is not valid"
        case .notFinite: "an angle or a quaternion is not finite or has no length"
        }
    }
}

/// The orientation of a crystal in the microscope: a proper rotation whose columns are the lab x, y and z axes in the crystal's
/// Cartesian frame (py4DSTEM's `Orientation.matrix`). Lab z points at the electron source; g appears at (Mᵀg).xy.
public struct CrystalOrientation: Sendable, Hashable, Codable, CustomStringConvertible {
    /// Columns: the lab x, y and z axes in the crystal's Cartesian frame. Always a proper rotation.
    public let matrix: simd_double3x3

    /// The largest departure from orthonormality and from determinant 1 a matrix may have.
    public static let rotationTolerance = 1e-8

    /// An orientation from its matrix (columns lab x, y, z in crystal coordinates).
    /// - Throws: ``CrystalOrientationError/notARotation`` unless the matrix is finite, orthonormal and of determinant +1, each
    ///   within ``rotationTolerance``.
    public init(matrix: simd_double3x3) throws(CrystalOrientationError) {
        guard Self.isRotation(matrix) else { throw .notARotation }
        self.matrix = matrix
    }

    /// The orientation that looks down the zone axis [u v w] of `cell`, turned in the plane by `inPlaneAngle`.
    ///
    /// The zone-axis frame is py4DSTEM's `parse_orientation(zone_axis_lattice=[u v w])`: z = u a₁ + v a₂ + w a₃ normalised;
    /// the x reference is ŷ × z when |z_z| > 1 − 1e-6, else (0, 0, −1); y = z × x, x = y × z. The in-plane angle multiplies on
    /// the right as py4DSTEM's matcher does (`orientation_matrix @ m3z`, m3z = [[c, s, 0], [−s, c, 0], [0, 0, 1]]): the pattern
    /// turns by +`inPlaneAngle` (counter-clockwise from lab x to lab y). The third column is the zone axis.
    /// - Parameters:
    ///   - zoneAxis: [u v w], a real-space direction of the cell (not reduced; [0 0 2] is [0 0 1]).
    ///   - inPlaneAngle: radians.
    ///   - cell: the crystal's cell.
    public init(zoneAxis uvw: SIMD3<Int>, inPlaneAngle: Double = 0, cell: UnitCell) throws(CrystalOrientationError) {
        guard inPlaneAngle.isFinite else { throw .notFinite }
        guard cell.isValid, cell.realLattice.count == 3 else { throw .invalidCell }
        let direction = Double(uvw.x) * cell.realLattice[0] + Double(uvw.y) * cell.realLattice[1]
            + Double(uvw.z) * cell.realLattice[2]
        let length = simd_length(direction)
        guard uvw != .zero, length.isFinite, length > 1e-12 else { throw .zeroZoneAxis }
        let basis = Self.zoneAxisFrame(direction)
        try self.init(matrix: basis * Self.inPlaneMatrix(inPlaneAngle))
    }

    /// An orientation from orix's quaternion (a, b, c, d), a the real part, as `Orientation(q).to_matrix()` gives it. The
    /// quaternion is normalised first; q and −q are the same orientation.
    public init(orixQuaternion q: SIMD4<Double>) throws(CrystalOrientationError) {
        let n = simd_length(q)
        guard q.x.isFinite, q.y.isFinite, q.z.isFinite, q.w.isFinite, n.isFinite, n > 1e-12 else { throw .notFinite }
        let u = q / n
        try self.init(matrix: simd_double3x3(simd_quatd(real: u.x, imag: SIMD3(u.y, u.z, u.w))))
    }

    /// orix's quaternion (a, b, c, d) of this orientation, a the real part, `Orientation.from_matrix(matrix).data`: unit length,
    /// a ≥ 0. At a half turn (a = 0) the sign is the one whose first non-zero imaginary component is positive.
    public var orixQuaternion: SIMD4<Double> {
        let q = simd_quatd(matrix)
        var v = SIMD4(q.real, q.imag.x, q.imag.y, q.imag.z)
        v /= simd_length(v)
        let negate: Bool
        if abs(v.x) < 1e-14 {
            v.x = 0
            negate = [v.y, v.z, v.w].first(where: { abs($0) > 1e-14 }).map { $0 < 0 } ?? false
        } else {
            negate = v.x < 0
        }
        return negate ? -v : v
    }

    /// The beam direction pointing at the source, a unit vector in the crystal's Cartesian frame (the third column).
    public var beamDirection: SIMD3<Double> { matrix.columns.2 }

    /// Where a crystal vector lies in the lab frame, Mᵀv: a reciprocal vector g appears in the pattern at (x, y); z is its
    /// component toward the source.
    public func labCoordinates(of v: SIMD3<Double>) -> SIMD3<Double> { matrix.transpose * v }

    /// The orientation turned by `angle` radians about lab z (the beam): the pattern turns counter-clockwise from lab x toward
    /// lab y. `CrystalOrientation(zoneAxis: z, inPlaneAngle: a, cell: c).rotatedInPlane(by: b)` is the orientation at a + b.
    public func rotatedInPlane(by angle: Double) -> CrystalOrientation {
        CrystalOrientation(unchecked: matrix * Self.inPlaneMatrix(angle))
    }

    /// The orientation after the crystal is tilted about lab x by `angles.x` and then about lab y by `angles.y` (radians, right
    /// hand about each axis; a double-tilt holder whose axes lie along lab x and y). The tilt R acts in the lab frame, so the
    /// new matrix is M Rᵀ with R = R_y(angles.y) R_x(angles.x), and every pattern vector moves as p → R p.
    public func tilted(by angles: SIMD2<Double>) -> CrystalOrientation {
        let rx = simd_double3x3(simd_quatd(angle: angles.x, axis: SIMD3(1, 0, 0)))
        let ry = simd_double3x3(simd_quatd(angle: angles.y, axis: SIMD3(0, 1, 0)))
        return CrystalOrientation(unchecked: matrix * (ry * rx).transpose)
    }

    /// The low-index zone axis nearest the beam and the angle between them (radians): over every [u v w] with
    /// |u|, |v|, |w| ≤ `maxIndex` and no common factor, the one whose direction u a₁ + v a₂ + w a₃ makes the smallest angle
    /// with ``beamDirection``. It points at the source (the angle is below 90°). Equal angles (within 1e-12 rad) go to the
    /// smaller largest index, then the smaller index sum, then the larger (u, v, w) in order. nil when `maxIndex` is not in
    /// 1...64 or the cell is not valid.
    public func nearestZoneAxis(cell: UnitCell, maxIndex: Int = 3) -> (uvw: SIMD3<Int>, angle: Double)? {
        guard maxIndex >= 1, maxIndex <= 64, cell.isValid, cell.realLattice.count == 3 else { return nil }
        let n = beamDirection
        var best: (uvw: SIMD3<Int>, angle: Double, key: (Int, Int))?
        for u in -maxIndex...maxIndex { for v in -maxIndex...maxIndex { for w in -maxIndex...maxIndex {
            guard Self.gcd(Self.gcd(abs(u), abs(v)), abs(w)) == 1 else { continue }
            let d = Double(u) * cell.realLattice[0] + Double(v) * cell.realLattice[1] + Double(w) * cell.realLattice[2]
            let cosine = simd_dot(d, n)
            guard cosine > 0 else { continue }
            let angle = atan2(simd_length(simd_cross(d, n)), cosine)
            let key = (max(abs(u), abs(v), abs(w)), abs(u) + abs(v) + abs(w))
            if let b = best {
                if angle > b.angle + 1e-12 { continue }
                if angle >= b.angle - 1e-12 {
                    if key > b.key { continue }
                    if key == b.key, (u, v, w) < (b.uvw.x, b.uvw.y, b.uvw.z) { continue }
                }
            }
            best = (SIMD3(u, v, w), angle, key)
        } } }
        return best.map { ($0.uvw, $0.angle) }
    }

    public var description: String {
        let c = [matrix.columns.0, matrix.columns.1, matrix.columns.2]
        return "CrystalOrientation(x: \(c[0]), y: \(c[1]), z: \(c[2]))"
    }

    // MARK: Hashable, Codable

    public static func == (a: CrystalOrientation, b: CrystalOrientation) -> Bool {
        a.matrix.columns.0 == b.matrix.columns.0 && a.matrix.columns.1 == b.matrix.columns.1 && a.matrix.columns.2 == b.matrix.columns.2
    }

    public func hash(into hasher: inout Hasher) {
        for c in [matrix.columns.0, matrix.columns.1, matrix.columns.2] { hasher.combine(c.x); hasher.combine(c.y); hasher.combine(c.z) }
    }

    /// The written form's version. Version 1: `{"schemaVersion": 1, "columns": [[x], [y], [z]]}`, the three columns of the
    /// matrix (lab x, y and z in crystal coordinates).
    public static let schemaVersion = 1

    private enum CodingKeys: String, CodingKey { case schemaVersion, columns }

    public init(from decoder: any Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        let version = try c.decode(Int.self, forKey: .schemaVersion)
        guard version == Self.schemaVersion else {
            throw DecodingError.dataCorruptedError(forKey: .schemaVersion, in: c,
                                                   debugDescription: "CrystalOrientation schema version \(version) is not 1")
        }
        let columns = try c.decode([[Double]].self, forKey: .columns)
        guard columns.count == 3, columns.allSatisfy({ $0.count == 3 }) else {
            throw DecodingError.dataCorruptedError(forKey: .columns, in: c, debugDescription: "three columns of three numbers")
        }
        let m = simd_double3x3(columns: (SIMD3(columns[0][0], columns[0][1], columns[0][2]),
                                         SIMD3(columns[1][0], columns[1][1], columns[1][2]),
                                         SIMD3(columns[2][0], columns[2][1], columns[2][2])))
        guard let checked = try? CrystalOrientation(matrix: m) else {
            throw DecodingError.dataCorruptedError(forKey: .columns, in: c, debugDescription: "not a proper rotation")
        }
        self = checked
    }

    public func encode(to encoder: any Encoder) throws {
        var c = encoder.container(keyedBy: CodingKeys.self)
        try c.encode(Self.schemaVersion, forKey: .schemaVersion)
        let cols = [matrix.columns.0, matrix.columns.1, matrix.columns.2].map { [$0.x, $0.y, $0.z] }
        try c.encode(cols, forKey: .columns)
    }

    // MARK: Helpers

    private init(unchecked matrix: simd_double3x3) { self.matrix = matrix }

    /// py4DSTEM's `parse_orientation` frame of a direction (not necessarily unit): columns x, y, z.
    static func zoneAxisFrame(_ direction: SIMD3<Double>) -> simd_double3x3 {
        var x: SIMD3<Double>
        // `lattice_to_cartesian` returns the unit direction, so the |z_z| test is on the normalised z, as here.
        let zn = direction / simd_length(direction)
        if abs(zn.z) > 1 - 1e-6 { x = simd_cross(SIMD3(0, 1, 0), zn) } else { x = SIMD3(0, 0, -1) }
        let y = simd_cross(zn, x)
        x = simd_cross(y, zn)
        return simd_double3x3(columns: (x / simd_length(x), y / simd_length(y), zn))
    }

    /// py4DSTEM's m3z(γ) = [[c, s, 0], [−s, c, 0], [0, 0, 1]] (rows).
    static func inPlaneMatrix(_ angle: Double) -> simd_double3x3 {
        let c = cos(angle), s = sin(angle)
        return simd_double3x3(columns: (SIMD3(c, -s, 0), SIMD3(s, c, 0), SIMD3(0, 0, 1)))
    }

    static func isRotation(_ m: simd_double3x3) -> Bool {
        let cols = [m.columns.0, m.columns.1, m.columns.2]
        guard cols.allSatisfy({ $0.x.isFinite && $0.y.isFinite && $0.z.isFinite }) else { return false }
        let g = m.transpose * m
        for i in 0..<3 { for j in 0..<3 where abs(g[i][j] - (i == j ? 1 : 0)) > rotationTolerance { return false } }
        return abs(simd_determinant(m) - 1) <= rotationTolerance
    }

    static func gcd(_ a: Int, _ b: Int) -> Int { IntegerMath.gcd(a, b) }
}
