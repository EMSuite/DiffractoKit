//
//  DoubleTiltHolder.swift
//  Role: the tilts of a double-tilt holder that bring a crystal direction onto the beam, the orientation after a tilt, and the
//        calibration of the holder's α axis in the pattern from two indexed patterns of one crystal (2.4 tilt planner, S12).
//
//  Model (the same planted truth as tools/pins/tilt-planner/planted.py):
//    a sample vector s goes to the pattern frame as p = Rz(φ) Rx(α) Ry(σβ) s: α about the holder's axis, β about the axis
//    α carries, φ the angle of the α axis in the pattern from lab x toward lab y (it depends on the microscope and the
//    camera length, so it is calibrated), σ = ±1 the sense of the β readout (holders and microscopes differ, and the
//    image rotation between stage and camera can add a mirror; a half turn of φ reverses α and β together, nothing but σ
//    reverses β alone, #222). With the kit's orientation matrix M (columns lab x, y, z in the crystal frame) this is
//    Mᵀ = Rz(φ) Rx(α) Ry(σβ) S for a fixed S. The tilts that bring a crystal direction t (S t in the sample frame at zero
//    tilt) onto the beam are α = asin t_y, β = σ atan2(−t_x, t_z).
//  Published tools for the same job: ALPHABETA (Cautaerts et al., J. Microsc. 2019, doi 10.1111/jmi.12774), τompas (Xie
//  and Zhang, J. Appl. Cryst. 2020), pycotem (Mompiou and Xie), KSpaceNavigator (Duden et al.), CrysTBox. None exposes a
//  callable function at a pinned version; the kit is held on planted cases (pin P51).
//  Signs and axis directions differ between holders and microscopes; nothing is assumed: tilts are given only with a
//  calibrated axis angle, and the calibration fits both β senses and refuses a pair that cannot tell them apart.
//

import Foundation
import simd
import DiffractoCore

/// A double-tilt holder: its limits, the sense of its β readout and, once calibrated, the angle of its α axis in the
/// pattern. Angles in degrees.
public struct DoubleTiltHolder: Sendable, Codable, Equatable {

    /// The sense σ of the β readout in the model p = Rz(φ) Rx(α) Ry(σβ) s: `.positive` (+1) or `.negative` (−1, β reads
    /// the other way). The raw value is σ.
    public enum BetaSense: Int, Sendable, Codable, Equatable, Hashable, CaseIterable {
        case positive = 1
        case negative = -1

        /// σ as a number, +1 or −1.
        public var sign: Double { Double(rawValue) }
    }

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
    /// The sense of the β readout; found by ``calibrateAxis(_:at:_:at:symmetry:betaSense:)`` or set by hand.
    public var betaSense: BetaSense

    public init(alphaLimitDeg: Double = 30, betaLimitDeg: Double = 30, axisAngleDeg: Double? = nil,
                betaSense: BetaSense = .positive) {
        self.alphaLimitDeg = alphaLimitDeg
        self.betaLimitDeg = betaLimitDeg
        self.axisAngleDeg = axisAngleDeg
        self.betaSense = betaSense
    }

    private enum CodingKeys: String, CodingKey {
        case alphaLimitDeg, betaLimitDeg, axisAngleDeg, betaSense
    }

    /// A holder saved before the β sense existed decodes with `.positive`, the model it was saved under.
    public init(from decoder: any Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        alphaLimitDeg = try c.decode(Double.self, forKey: .alphaLimitDeg)
        betaLimitDeg = try c.decode(Double.self, forKey: .betaLimitDeg)
        axisAngleDeg = try c.decodeIfPresent(Double.self, forKey: .axisAngleDeg)
        betaSense = try c.decodeIfPresent(BetaSense.self, forKey: .betaSense) ?? .positive
    }

    /// The tilts that bring `direction` (crystal Cartesian, e.g. u a₁ + v a₂ + w a₃) onto the beam, from `orientation`
    /// measured at `current`. nil while the axis is not calibrated, or for a direction that is not finite and non-zero.
    public func tilts(toBring direction: SIMD3<Double>, orientation: CrystalOrientation, at current: Tilts) -> Target? {
        guard let phi = axisAngleDeg, phi.isFinite, current.alpha.isFinite, current.beta.isFinite else { return nil }
        let n = length(direction)
        guard n.isFinite, n > 0 else { return nil }
        let lab = orientation.matrix.transpose * (direction / n)
        let sample = Self.holder(phi: phi, current, sense: betaSense).transpose * lab
        let alpha = asin(min(max(sample.y, -1), 1)) * 180 / .pi
        let beta = betaSense.sign * atan2(-sample.x, sample.z) * 180 / .pi
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
        let m = orientation.matrix * Self.holder(phi: phi, current, sense: betaSense)
            * Self.holder(phi: phi, new, sense: betaSense).transpose
        return try? CrystalOrientation(matrix: m)
    }

    // MARK: Calibration

    /// The α axis and the β sense found from two indexed patterns of one crystal at two known tilts.
    public struct AxisCalibration: Sendable, Equatable {
        /// The angle of the α axis in the pattern, in [0, 360) degrees.
        public let axisAngleDeg: Double
        /// The sense of the β readout the pair fits (the one given, when the calibration was asked for one sense).
        public let betaSense: BetaSense
        /// The best misfit of the other β sense (degrees), above ``DoubleTiltHolder/maximumMisfitDeg`` whenever the
        /// calibration was accepted; nil when the calibration was asked for one sense.
        public let otherSenseMisfitDeg: Double?
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

        /// The holder with this axis and β sense and the given limits.
        public func holder(alphaLimitDeg: Double = 30, betaLimitDeg: Double = 30) -> DoubleTiltHolder {
            DoubleTiltHolder(alphaLimitDeg: alphaLimitDeg, betaLimitDeg: betaLimitDeg, axisAngleDeg: axisAngleDeg,
                             betaSense: betaSense)
        }
    }

    /// Why a calibration was refused.
    public enum CalibrationError: Error, Sendable, Equatable, CustomStringConvertible {
        /// The stage moved less than ``minimumTiltChangeDeg``: the axis is not determined.
        case tiltChangeTooSmall(Double)
        /// The best fit of either β sense leaves more than ``maximumMisfitDeg`` between the measured and the stage's
        /// rotation, while the two rotate by the same angle: the axis of the measured rotation is out of place, so α and β
        /// were swapped, the β mechanism is not a pure rotation about the axis α carries, or a pattern was indexed wrong.
        case poorFit(misfitDeg: Double)
        /// The measured rotation between the two patterns (`measuredDeg`, the closest of the half turns and symmetric
        /// equivalents) and the stage's rotation between the two tilts differ in angle by more than ``maximumMisfitDeg``;
        /// no axis or β sense can fix that: a tilt was read wrong (or in other units) or a
        /// pattern was indexed wrong.
        case rotationMismatch(measuredDeg: Double, stageDeg: Double)
        /// Both β senses fit the pair within ``maximumMisfitDeg`` (`misfitDeg` the better, `otherSenseMisfitDeg` the
        /// other), so the pair cannot tell which way β turns; a pair that changes only one tilt never can. No axis is given.
        case ambiguousSense(misfitDeg: Double, otherSenseMisfitDeg: Double)
        /// A tilt is not finite.
        case notFinite

        public var description: String {
            switch self {
            case .tiltChangeTooSmall(let d):
                return String(format: "The two patterns are only %.1f° apart; tilt by at least %.0f° (10° or more is better) "
                              + "to find the tilt axis.", d, DoubleTiltHolder.minimumTiltChangeDeg)
            case .poorFit(let m):
                return String(format: "The two patterns are %.1f° from what the stage's tilts predict for any tilt axis and "
                              + "either β sense, although they turn by the angle the stage says: the turn is about the "
                              + "wrong axis. Check that α and β were not swapped and the indexing of both patterns; a β "
                              + "that does not turn about the axis α carries does not fit this holder model.", m)
            case .rotationMismatch(let measured, let stage):
                return String(format: "The two patterns are %.1f° apart but the stage's tilts say %.1f°; no tilt axis can "
                              + "fix that. Check the tilts read from the files (and their units) and the indexing of both "
                              + "patterns.", measured, stage)
            case .ambiguousSense(let m, let other):
                return String(format: "Both senses of β fit these two patterns (%.2f° and %.2f°), so they cannot tell "
                              + "which way β turns; a pair that changes only one tilt never can. Record a pattern whose α "
                              + "and β both differ by a few degrees from one of these and calibrate with that pair, or set "
                              + "the β sense by hand.", m, other)
            case .notFinite:
                return "A tilt is not a finite number."
            }
        }
    }

    /// The smallest stage rotation a calibration accepts (degrees).
    public static let minimumTiltChangeDeg = 2.0
    /// The largest misfit a calibration accepts (degrees).
    public static let maximumMisfitDeg = 1.0

    /// The α axis angle and the β sense from two orientations of one crystal, `a` measured at `tiltsA` and `b` at `tiltsB`:
    /// for each β sense σ, the angle φ that makes Rz(φ) Q_σ Rz(φ)ᵀ closest to the measured rotation M_bᵀ M_a, Q_σ the stage's
    /// rotation between the two tilts (a 0.5° scan, then a golden-section refinement); the sense that fits better is taken.
    /// Both half turns of `b` about the beam are tried. `symmetry`: proper rotations of the crystal in its Cartesian frame
    /// (``properRotations(of:)``), tried between the two so that two patterns indexed to different symmetric equivalents
    /// still pair; empty takes the identity alone. `betaSense`: nil fits both senses and refuses with
    /// ``CalibrationError/ambiguousSense(misfitDeg:otherSenseMisfitDeg:)`` when both fit within ``maximumMisfitDeg`` (a pair
    /// that changes only one tilt always does); a sense fits that sense alone (one set by hand).
    public static func calibrateAxis(_ a: CrystalOrientation, at tiltsA: Tilts, _ b: CrystalOrientation, at tiltsB: Tilts,
                                     symmetry: [simd_double3x3] = [],
                                     betaSense: BetaSense? = nil) throws(CalibrationError) -> AxisCalibration {
        guard tiltsA.alpha.isFinite, tiltsA.beta.isFinite, tiltsB.alpha.isFinite, tiltsB.beta.isFinite else { throw .notFinite }
        // The stage turns by the same angle for both senses: Q₋ = X Q₊ X with the mirror X = diag(−1, 1, 1), which keeps
        // Rx(α) and reverses Ry(β).
        let change = rotationAngle(holder(phi: 0, tiltsB) * holder(phi: 0, tiltsA).transpose) * 180 / .pi
        guard change >= minimumTiltChangeDeg else { throw .tiltChangeTooSmall(change) }
        let senses = betaSense.map { [$0] } ?? BetaSense.allCases
        let ops = symmetry.isEmpty ? [matrix_identity_double3x3] : symmetry
        struct Fit { let sense: BetaSense; let misfit: Double; let phi: Double; let half: Bool }
        var perSense: [(best: Fit, otherHalfMisfit: Double)] = []
        // The measured turn closest in angle to the stage's: the misfit of any axis is at least their difference.
        var closestTurnDeg = Double.infinity
        for sense in senses {
            let q = holder(phi: 0, tiltsB, sense: sense) * holder(phi: 0, tiltsA, sense: sense).transpose
            var results: [Fit] = []
            for half in [false, true] {
                let mb = half ? b.matrix * rz(.pi) : b.matrix
                var best: Fit?
                for g in ops {
                    let r = mb.transpose * g * a.matrix
                    let turn = rotationAngle(r) * 180 / .pi
                    if abs(turn - change) < abs(closestTurnDeg - change) { closestTurnDeg = turn }
                    let fit = bestAxis(measured: r, stage: q)
                    if best == nil || fit.misfit < best!.misfit {
                        best = Fit(sense: sense, misfit: fit.misfit, phi: fit.phi, half: half)
                    }
                }
                results.append(best!)
            }
            results.sort { $0.misfit < $1.misfit }
            perSense.append((results[0], results[1].misfit * 180 / .pi))
        }
        perSense.sort { $0.best.misfit < $1.best.misfit }
        let chosen = perSense[0]
        let misfitDeg = chosen.best.misfit * 180 / .pi
        guard misfitDeg <= maximumMisfitDeg else {
            if abs(closestTurnDeg - change) > maximumMisfitDeg {
                throw .rotationMismatch(measuredDeg: closestTurnDeg, stageDeg: change)
            }
            throw .poorFit(misfitDeg: misfitDeg)
        }
        let otherSenseMisfitDeg: Double? = perSense.count > 1 ? perSense[1].best.misfit * 180 / .pi : nil
        if let other = otherSenseMisfitDeg, other <= maximumMisfitDeg {
            throw .ambiguousSense(misfitDeg: misfitDeg, otherSenseMisfitDeg: other)
        }
        var phiDeg = (chosen.best.phi * 180 / .pi).truncatingRemainder(dividingBy: 360)
        if phiDeg < 0 { phiDeg += 360 }
        return AxisCalibration(axisAngleDeg: phiDeg, betaSense: chosen.best.sense, otherSenseMisfitDeg: otherSenseMisfitDeg,
                               tiltChangeDeg: change, misfitDeg: misfitDeg, secondNeedsHalfTurn: chosen.best.half,
                               otherHalfTurnMisfitDeg: chosen.otherHalfMisfit)
    }

    /// The proper rotations of a structure's Laue group in its Cartesian frame, for ``calibrateAxis(_:at:_:at:symmetry:betaSense:)``:
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

    /// H = Rz(φ) Rx(α) Ry(σβ), angles in degrees.
    static func holder(phi: Double, _ t: Tilts, sense: BetaSense = .positive) -> simd_double3x3 {
        let d = Double.pi / 180
        return rz(phi * d) * rx(t.alpha * d) * ry(sense.sign * t.beta * d)
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
