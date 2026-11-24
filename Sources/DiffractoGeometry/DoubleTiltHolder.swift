//
//  DoubleTiltHolder.swift
//  Role: the tilts of a double-tilt holder that bring a crystal direction onto the beam, the orientation after a tilt, and the
//        calibration of the holder's α axis in the pattern from two indexed patterns of one crystal (2.4 tilt planner, S12).
//
//  Model (the same planted truth as tools/pins/tilt-planner/planted.py):
//    a sample vector s goes to the pattern frame as p = Rz(φ) Rx(α) Ry(β) s: α about the holder's axis, β about the axis α
//    carries, φ the angle of the α axis in the pattern from lab x toward lab y (it depends on the microscope and the camera
//    length, so it is calibrated). With the kit's orientation matrix M (columns lab x, y, z in the crystal frame) this is
//    Mᵀ = Rz(φ) Rx(α) Ry(β) S for a fixed S. The tilts that bring a crystal direction t (S t in the sample frame at zero
//    tilt) onto the beam are α = asin t_y, β = atan2(−t_x, t_z).
//  Published tools for the same job: ALPHABETA (Cautaerts et al., J. Microsc. 2019, doi 10.1111/jmi.12774), τompas (Xie
//  and Zhang, J. Appl. Cryst. 2020), pycotem (Mompiou and Xie), KSpaceNavigator (Duden et al.), CrysTBox. None exposes a
//  callable function at a pinned version; the kit is held on planted cases (pin P51).
//  Signs and axis directions differ between holders and microscopes; nothing is assumed: tilts are given only with a
//  calibrated axis angle.
//

import Foundation
import simd
import DiffractoCore

/// A double-tilt holder: its limits and, once calibrated, the angle of its α axis in the pattern. Angles in degrees.
public struct DoubleTiltHolder: Sendable, Codable, Equatable {

    /// Holder tilts (degrees).
    public struct Tilts: Sendable, Codable, Equatable, Hashable {
        public var alpha: Double
        public var beta: Double
        public init(alpha: Double, beta: Double) {
            self.alpha = alpha
            self.beta = beta
        }
        public static let zero = Tilts(alpha: 0, beta: 0)
    }

    /// The tilts that bring a direction onto the beam, and whether the holder reaches them.
    public struct Target: Sendable, Equatable {
        public let tilts: Tilts
        /// |α| and |β| within the limits.
        public let reachable: Bool
        /// The angle between the direction and the beam now (degrees).
        public let angleFromBeamDeg: Double
    }

    /// ±limit of α (degrees).
    public var alphaLimitDeg: Double
    /// ±limit of β (degrees).
    public var betaLimitDeg: Double
    /// The angle of the α axis in the pattern, from lab x toward lab y (degrees); nil until calibrated.
    public var axisAngleDeg: Double?

    public init(alphaLimitDeg: Double = 30, betaLimitDeg: Double = 30, axisAngleDeg: Double? = nil) {
        self.alphaLimitDeg = alphaLimitDeg
        self.betaLimitDeg = betaLimitDeg
        self.axisAngleDeg = axisAngleDeg
    }

    /// The tilts that bring `direction` (crystal Cartesian, e.g. u a₁ + v a₂ + w a₃) onto the beam, from `orientation`
    /// measured at `current`. nil while the axis is not calibrated, or for a direction that is not finite and non-zero.
    public func tilts(toBring direction: SIMD3<Double>, orientation: CrystalOrientation, at current: Tilts) -> Target? {
        guard let phi = axisAngleDeg, phi.isFinite, current.alpha.isFinite, current.beta.isFinite else { return nil }
        let n = length(direction)
        guard n.isFinite, n > 0 else { return nil }
        let lab = orientation.matrix.transpose * (direction / n)
        let sample = Self.holder(phi: phi, current).transpose * lab
        let alpha = asin(min(max(sample.y, -1), 1)) * 180 / .pi
        let beta = atan2(-sample.x, sample.z) * 180 / .pi
        let tilts = Tilts(alpha: alpha, beta: beta)
        let reachable = abs(alpha) <= alphaLimitDeg && abs(beta) <= betaLimitDeg
        let angle = atan2((lab.x * lab.x + lab.y * lab.y).squareRoot(), lab.z) * 180 / .pi
        return Target(tilts: tilts, reachable: reachable, angleFromBeamDeg: angle)
    }

    /// The tilts that bring the zone axis [u v w] of `cell` onto the beam.
    public func tilts(toZoneAxis uvw: SIMD3<Int>, cell: UnitCell, orientation: CrystalOrientation, at current: Tilts) -> Target? {
        guard cell.isValid, cell.realLattice.count == 3, uvw != .zero else { return nil }
        let d = Double(uvw.x) * cell.realLattice[0] + Double(uvw.y) * cell.realLattice[1] + Double(uvw.z) * cell.realLattice[2]
        return tilts(toBring: d, orientation: orientation, at: current)
    }

    /// The orientation after the holder moves from `current` to `new`, from `orientation` measured at `current`. nil while
    /// the axis is not calibrated.
    public func orientation(_ orientation: CrystalOrientation, at current: Tilts, movedTo new: Tilts) -> CrystalOrientation? {
        guard let phi = axisAngleDeg, phi.isFinite, current.alpha.isFinite, current.beta.isFinite, new.alpha.isFinite,
              new.beta.isFinite else { return nil }
        // M₁ᵀ = H₁ H₀ᵀ M₀ᵀ, so M₁ = M₀ H₀ H₁ᵀ
        let m = orientation.matrix * Self.holder(phi: phi, current) * Self.holder(phi: phi, new).transpose
        return try? CrystalOrientation(matrix: m)
    }

    // MARK: Calibration

    /// The α axis found from two indexed patterns of one crystal at two known tilts.
    public struct AxisCalibration: Sendable, Equatable {
        /// The angle of the α axis in the pattern, in [0, 360) degrees.
        public let axisAngleDeg: Double
        /// The rotation the stage reports between the two tilts (degrees): the larger, the better the axis is fixed
        /// (planted, 0.2° orientation error per pattern: median 0.48°, 95th percentile 1.4° for 10 to 20°).
        public let tiltChangeDeg: Double
        /// The angle between the measured rotation and the stage's rotation at the found axis (degrees): how well the pair
        /// fits a holder of this model.
        public let misfitDeg: Double
        /// The second pattern fits only turned by a half turn about the beam: its indexing picked the other of the two
        /// orientations one spot pattern cannot tell apart. Use `secondOrientation.rotatedInPlane(by: .pi)`.
        public let secondNeedsHalfTurn: Bool
        /// The misfit of the other half turn (degrees); when it is also below the limit the pair cannot tell them apart.
        public let otherHalfTurnMisfitDeg: Double
    }

    /// Why a calibration was refused.
    public enum CalibrationError: Error, Sendable, Equatable, CustomStringConvertible {
        /// The stage moved less than ``minimumTiltChangeDeg``: the axis is not determined.
        case tiltChangeTooSmall(Double)
        /// The best fit leaves more than ``maximumMisfitDeg`` between the measured and the stage's rotation: the β mechanism
        /// is not a pure rotation, a tilt was read wrong, or a pattern was indexed wrong.
        case poorFit(misfitDeg: Double)
        /// A tilt is not finite.
        case notFinite

        public var description: String {
            switch self {
            case .tiltChangeTooSmall(let d):
                return String(format: "The two patterns are only %.1f° apart; tilt by at least %.0f° (10° or more is better) "
                              + "to find the tilt axis.", d, DoubleTiltHolder.minimumTiltChangeDeg)
            case .poorFit(let m):
                return String(format: "The two patterns are %.1f° from what the stage's tilts predict for any tilt axis; check "
                              + "the tilts read from the files and the indexing of both patterns.", m)
            case .notFinite:
                return "A tilt is not a finite number."
            }
        }
    }

    /// The smallest stage rotation a calibration accepts (degrees).
    public static let minimumTiltChangeDeg = 2.0
    /// The largest misfit a calibration accepts (degrees).
    public static let maximumMisfitDeg = 1.0

    /// The α axis angle from two orientations of one crystal, `a` measured at `tiltsA` and `b` at `tiltsB`: the angle φ that
    /// makes Rz(φ) Q Rz(φ)ᵀ closest to the measured rotation M_bᵀ M_a, Q the stage's rotation between the two tilts (a 0.5°
    /// scan, then a golden-section refinement). Both half turns of `b` about the beam are tried. `symmetry`: proper rotations
    /// of the crystal in its Cartesian frame (``properRotations(of:)``), tried between the two so that two patterns indexed to
    /// different symmetric equivalents still pair; empty takes the identity alone.
    public static func calibrateAxis(_ a: CrystalOrientation, at tiltsA: Tilts, _ b: CrystalOrientation, at tiltsB: Tilts,
                                     symmetry: [simd_double3x3] = []) throws(CalibrationError) -> AxisCalibration {
        guard tiltsA.alpha.isFinite, tiltsA.beta.isFinite, tiltsB.alpha.isFinite, tiltsB.beta.isFinite else { throw .notFinite }
        let q = holder(phi: 0, tiltsB) * holder(phi: 0, tiltsA).transpose
        let change = rotationAngle(q) * 180 / .pi
        guard change >= minimumTiltChangeDeg else { throw .tiltChangeTooSmall(change) }
        let ops = symmetry.isEmpty ? [matrix_identity_double3x3] : symmetry
        var results: [(misfit: Double, phi: Double, half: Bool)] = []
        for half in [false, true] {
            let mb = half ? b.matrix * rz(.pi) : b.matrix
            var best: (misfit: Double, phi: Double)?
            for g in ops {
                let r = mb.transpose * g * a.matrix
                let fit = bestAxis(measured: r, stage: q)
                if best == nil || fit.misfit < best!.misfit { best = fit }
            }
            results.append((best!.misfit, best!.phi, half))
        }
        results.sort { $0.misfit < $1.misfit }
        let chosen = results[0]
        let misfitDeg = chosen.misfit * 180 / .pi
        guard misfitDeg <= maximumMisfitDeg else { throw .poorFit(misfitDeg: misfitDeg) }
        var phiDeg = (chosen.phi * 180 / .pi).truncatingRemainder(dividingBy: 360)
        if phiDeg < 0 { phiDeg += 360 }
        return AxisCalibration(axisAngleDeg: phiDeg, tiltChangeDeg: change, misfitDeg: misfitDeg,
                               secondNeedsHalfTurn: chosen.half, otherHalfTurnMisfitDeg: results[1].misfit * 180 / .pi)
    }

    /// The proper rotations of a structure's Laue group in its Cartesian frame, for ``calibrateAxis(_:at:_:at:symmetry:)``:
    /// W acts on fractional coordinates, so on the Cartesian frame it is A W A⁻¹ (A: columns a₁, a₂, a₃), as the zone-axis
    /// solver takes it. A rotation the cell does not carry (not orthogonal within 1e-6) is skipped.
    public static func properRotations(of structure: CrystalStructure) -> [simd_double3x3] {
        cartesianProperRotations(of: structure).map(\.cartesian)
    }

    /// The same, each with the integer matrix W it comes from (the zone-axis solver keeps both; one copy of the loop, #226).
    static func cartesianProperRotations(of structure: CrystalStructure) -> [(cartesian: simd_double3x3, integer: IntegerMatrix)] {
        let lattice = structure.cell.realLattice
        guard lattice.count == 3 else { return [] }
        let a = simd_double3x3(columns: (lattice[0], lattice[1], lattice[2]))
        let inverse = a.inverse
        var out: [(cartesian: simd_double3x3, integer: IntegerMatrix)] = []
        for w in structure.laueGroup.properRotations {
            let wm = simd_double3x3(rows: w.rows.map { SIMD3<Double>(Double($0.x), Double($0.y), Double($0.z)) })
            let c = a * wm * inverse
            let e = c * c.transpose - matrix_identity_double3x3
            let err = [e.columns.0, e.columns.1, e.columns.2].map { simd_reduce_max(simd_abs($0)) }.max() ?? 1
            if err <= 1e-6 { out.append((c, w)) }
        }
        return out
    }

    // MARK: Geometry

    /// H = Rz(φ) Rx(α) Ry(β), angles in degrees.
    static func holder(phi: Double, _ t: Tilts) -> simd_double3x3 {
        let d = Double.pi / 180
        return rz(phi * d) * rx(t.alpha * d) * ry(t.beta * d)
    }

    static func rx(_ a: Double) -> simd_double3x3 {
        let c = cos(a), s = sin(a)
        return simd_double3x3(rows: [SIMD3(1, 0, 0), SIMD3(0, c, -s), SIMD3(0, s, c)])
    }

    static func ry(_ a: Double) -> simd_double3x3 {
        let c = cos(a), s = sin(a)
        return simd_double3x3(rows: [SIMD3(c, 0, s), SIMD3(0, 1, 0), SIMD3(-s, 0, c)])
    }

    static func rz(_ a: Double) -> simd_double3x3 {
        let c = cos(a), s = sin(a)
        return simd_double3x3(rows: [SIMD3(c, -s, 0), SIMD3(s, c, 0), SIMD3(0, 0, 1)])
    }

    /// The rotation angle of a rotation matrix (radians), from its quaternion: accurate near 0.
    static func rotationAngle(_ r: simd_double3x3) -> Double {
        let q = simd_quatd(r)
        return 2 * atan2(length(q.imag), abs(q.real))
    }

    /// φ (radians) minimising the angle between `measured` and Rz(φ) `stage` Rz(φ)ᵀ, and that angle.
    static func bestAxis(measured: simd_double3x3, stage: simd_double3x3) -> (misfit: Double, phi: Double) {
        func f(_ p: Double) -> Double {
            let z = rz(p)
            return rotationAngle(measured * (z * stage * z.transpose).transpose)
        }
        let step = 0.5 * Double.pi / 180
        var p0 = 0.0, f0 = Double.infinity
        for i in 0..<720 {
            let p = Double(i) * step
            let v = f(p)
            if v < f0 { f0 = v; p0 = p }
        }
        var lo = p0 - step, hi = p0 + step
        let g = (5.0.squareRoot() - 1) / 2
        for _ in 0..<60 {
            let c = hi - g * (hi - lo), d = lo + g * (hi - lo)
            if f(c) < f(d) { hi = d } else { lo = c }
        }
        let p = (lo + hi) / 2
        return (f(p), p)
    }
}
