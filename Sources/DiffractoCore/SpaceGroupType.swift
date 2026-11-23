//
//  SpaceGroupType.swift
//  Role: the space-group type of a set of operations, as spglib 2.8.0 `get_spacegroup_type_from_symmetry` names it (IT
//        number, standard Hall number, short international symbol). The operations are compared with the 530 settings
//        spglib lists (`SpaceGroupTable.hallSymbols`, read by `HallSymbol`): same rotations, same centring, and
//        translations that one origin shift carries onto the setting's. The shift is not searched: the congruences
//        (I − R) p ≡ t − t₀ modulo the lattice are solvable exactly when every row the integer elimination of the stacked
//        (I − R) clears to zero leaves an integer on the right (a row echelon form over the integers, on the primitive
//        basis of the centring, where the lattice is Z³).
//

import simd

/// A space-group type and the setting the operations were written in.
public struct SpaceGroupType: Sendable, Equatable, CustomStringConvertible {
    /// The IT number, 1 ... 230.
    public let number: Int
    /// spglib's standard Hall number of the type (1 ... 530), as `get_spacegroup_type_from_symmetry` gives it.
    public let hallNumber: Int
    /// spglib's short international symbol of the standard setting ("P2_1/c", "Fm-3m", "R-3m").
    public let internationalSymbol: String
    /// The Hall number of the setting the operations are written in (one of spglib's 530), up to an origin shift.
    public let settingHallNumber: Int
    /// That setting's Hall symbol.
    public let settingHallSymbol: String

    public var description: String { "\(internationalSymbol) (No. \(number))" }

    /// The type of `operations` (a whole group, centring translations included, in the basis of one of spglib's 530
    /// settings with any origin); nil when they are no such group: a basis none of the settings uses (a reduced or
    /// rotated cell; use the conventional setting's operations), a set that is not closed, or translations that no
    /// origin shift carries onto a setting's within `translationTolerance` (fractional).
    public static func of(_ operations: [SymmetryOperation], translationTolerance: Double = 1e-3) -> SpaceGroupType? {
        guard !operations.isEmpty, translationTolerance.isFinite, translationTolerance >= 0 else { return nil }
        guard operations.allSatisfy({ $0.translation.x.isFinite && $0.translation.y.isFinite && $0.translation.z.isFinite })
        else { return nil }
        let rotations = Set(operations.map(\.rotation))
        guard let centring = Centring(operations: operations, tolerance: translationTolerance) else { return nil }
        // Operations of one rotation differ by centring vectors only, and each rotation carries every centring vector
        // (else the set is no whole group of this lattice).
        let inverse = centring.primitiveBasis.inverse
        var cosets: [IntegerMatrix: [SIMD3<Double>]] = [:]
        for op in operations {
            guard let t0 = cosets[op.rotation]?.first else { cosets[op.rotation] = [op.translation]; continue }
            guard Centring.isLatticeVector(inverse * (op.translation - t0), tolerance: translationTolerance) else { return nil }
            if !cosets[op.rotation, default: []].contains(where: {
                Centring.isLatticeVector($0 - op.translation, tolerance: translationTolerance)
            }) {
                cosets[op.rotation, default: []].append(op.translation)
            }
        }
        guard cosets.values.allSatisfy({ $0.count == 1 + centring.vectors.count }) else { return nil }
        for (index, setting) in Settings.all.enumerated() where setting.rotations == rotations && setting.centring == centring {
            if originShiftExists(from: setting, to: operations, centring: centring, tolerance: translationTolerance) {
                let hallNumber = index + 1
                let number = setting.number
                let standard = SpaceGroupTable.standardHallNumbers[number - 1]
                return SpaceGroupType(number: number, hallNumber: standard,
                                      internationalSymbol: SpaceGroupTable.internationalSymbols[number - 1],
                                      settingHallNumber: hallNumber,
                                      settingHallSymbol: SpaceGroupTable.hallSymbols[hallNumber - 1])
            }
        }
        return nil
    }

    // MARK: - Settings

    /// One of spglib's 530 settings: its rotations, centring, one translation per rotation and its IT number.
    struct Settings {
        let rotations: Set<IntegerMatrix>
        let centring: Centring
        let translationOf: [IntegerMatrix: SIMD3<Double>]
        let number: Int

        static let all: [Settings] = SpaceGroupTable.hallSymbols.enumerated().map { index, symbol in
            // Every symbol is read exactly by HallSymbol (P37); a failure here is a table error.
            let ops = (try? HallSymbol.operations(symbol)) ?? [.identity]
            var translationOf: [IntegerMatrix: SIMD3<Double>] = [:]
            for op in ops where translationOf[op.rotation] == nil { translationOf[op.rotation] = op.translation }
            return Settings(rotations: Set(ops.map(\.rotation)),
                            centring: Centring(operations: ops, tolerance: 1e-9) ?? .primitive,
                            translationOf: translationOf, number: SpaceGroupTable.numberOfHallNumber[index])
        }
    }

    // MARK: - Centring

    /// The centring of the settings, read from the pure translations, with the primitive basis that makes its lattice Z³.
    enum Centring: CaseIterable {
        case primitive, a, b, c, body, face, rhombohedral

        var vectors: [SIMD3<Double>] {
            switch self {
            case .primitive: []
            case .a: [SIMD3(0, 0.5, 0.5)]
            case .b: [SIMD3(0.5, 0, 0.5)]
            case .c: [SIMD3(0.5, 0.5, 0)]
            case .body: [SIMD3(0.5, 0.5, 0.5)]
            case .face: [SIMD3(0, 0.5, 0.5), SIMD3(0.5, 0, 0.5), SIMD3(0.5, 0.5, 0)]
            case .rhombohedral: [SIMD3(2.0 / 3, 1.0 / 3, 1.0 / 3), SIMD3(1.0 / 3, 2.0 / 3, 2.0 / 3)]
            }
        }

        /// The primitive basis vectors (columns) in the setting's fractional coordinates.
        var primitiveBasis: simd_double3x3 {
            switch self {
            case .primitive: matrix_identity_double3x3
            case .a: simd_double3x3(SIMD3(1, 0, 0), SIMD3(0, 0.5, 0.5), SIMD3(0, -0.5, 0.5))
            case .b: simd_double3x3(SIMD3(0.5, 0, 0.5), SIMD3(0, 1, 0), SIMD3(-0.5, 0, 0.5))
            case .c: simd_double3x3(SIMD3(0.5, 0.5, 0), SIMD3(-0.5, 0.5, 0), SIMD3(0, 0, 1))
            case .body: simd_double3x3(SIMD3(-0.5, 0.5, 0.5), SIMD3(0.5, -0.5, 0.5), SIMD3(0.5, 0.5, -0.5))
            case .face: simd_double3x3(SIMD3(0, 0.5, 0.5), SIMD3(0.5, 0, 0.5), SIMD3(0.5, 0.5, 0))
            case .rhombohedral:
                simd_double3x3(SIMD3(2.0 / 3, 1.0 / 3, 1.0 / 3), SIMD3(-1.0 / 3, 1.0 / 3, 1.0 / 3),
                               SIMD3(-1.0 / 3, -2.0 / 3, 1.0 / 3))
            }
        }

        /// The centring whose vectors are exactly the non-zero pure translations of `operations` (mod 1, within
        /// `tolerance`); nil for any other set.
        init?(operations: [SymmetryOperation], tolerance: Double) {
            var pure: [SIMD3<Double>] = []
            for op in operations where op.rotation == .identity {
                let t = op.translationMod1
                if Centring.isLatticeVector(t, tolerance: tolerance) { continue }
                if !pure.contains(where: { Centring.isLatticeVector($0 - t, tolerance: tolerance) }) { pure.append(t) }
            }
            guard let match = Centring.allCases.first(where: { centring in
                centring.vectors.count == pure.count && centring.vectors.allSatisfy { v in
                    pure.contains { Centring.isLatticeVector($0 - v, tolerance: tolerance) }
                }
            }) else { return nil }
            self = match
        }

        static func isLatticeVector(_ v: SIMD3<Double>, tolerance: Double) -> Bool {
            abs(v.x - v.x.rounded()) <= tolerance && abs(v.y - v.y.rounded()) <= tolerance
                && abs(v.z - v.z.rounded()) <= tolerance
        }
    }

    // MARK: - Origin shift

    /// Whether one origin shift p carries `setting`'s translations onto those of `operations`: for each rotation R,
    /// (I − R) p ≡ t − t₀ modulo the centring lattice. On the primitive basis B (lattice Z³) that is
    /// (I − B⁻¹RB) q ≡ B⁻¹(t − t₀) mod Z³; integer row operations on the stacked system bring (I − B⁻¹RB) to row echelon
    /// form, and a solution exists exactly when each zero row leaves an integer right-hand side.
    static func originShiftExists(from setting: Settings, to operations: [SymmetryOperation], centring: Centring,
                                  tolerance: Double) -> Bool {
        let basis = centring.primitiveBasis
        let inverse = basis.inverse
        var rows: [SIMD3<Int>] = []
        var rhs: [Double] = []
        var seen = Set<IntegerMatrix>()
        for op in operations where seen.insert(op.rotation).inserted {
            guard let t0 = setting.translationOf[op.rotation] else { return false }
            let r = simd_double3x3(rows: op.rotation.rows.map { SIMD3<Double>(Double($0.x), Double($0.y), Double($0.z)) })
            let primitive = inverse * r * basis
            let difference = inverse * (op.translation - t0)
            for row in 0..<3 {
                // Row `row` of I − B⁻¹RB; B⁻¹RB is an integer matrix for every symmetry of the centred lattice.
                var integers = SIMD3<Int>(0, 0, 0)
                for column in 0..<3 {
                    let value = (row == column ? 1 : 0) - primitive[column][row]
                    let rounded = value.rounded()
                    guard abs(value - rounded) <= 1e-9 else { return false }
                    integers[column] = Int(rounded)
                }
                rows.append(integers)
                rhs.append(difference[row])
            }
        }
        return integerSystemIsConsistent(rows: rows, rhs: rhs, tolerance: tolerance)
    }

    /// Whether A q ≡ b (mod Z^m) has a real solution q, A integer m×3: integer row operations (unimodular, so the
    /// modulus stays Z^m) bring A to row echelon form; a zero row must leave an integer in b. The tolerance on that
    /// integer grows with the multiplier the elimination applied to b.
    static func integerSystemIsConsistent(rows: [SIMD3<Int>], rhs: [Double], tolerance: Double) -> Bool {
        var a = rows
        var b = rhs
        var weight = [Double](repeating: 1, count: rows.count)  // a bound on the sum of |multipliers| behind each b
        var pivotRow = 0
        for column in 0..<3 where pivotRow < a.count {
            while true {
                // The row from pivotRow on with the smallest non-zero entry in this column becomes the pivot.
                var best: Int?
                for i in pivotRow..<a.count where a[i][column] != 0 {
                    if best == nil || abs(a[i][column]) < abs(a[best!][column]) { best = i }
                }
                guard let pivot = best else { break }
                a.swapAt(pivotRow, pivot)
                b.swapAt(pivotRow, pivot)
                weight.swapAt(pivotRow, pivot)
                var cleared = true
                for i in (pivotRow + 1)..<a.count where a[i][column] != 0 {
                    let q = a[i][column] / a[pivotRow][column]
                    a[i] &-= q &* a[pivotRow]
                    b[i] -= Double(q) * b[pivotRow]
                    weight[i] += Double(abs(q)) * weight[pivotRow]
                    if a[i][column] != 0 { cleared = false }
                }
                if cleared { break }
            }
            if a[pivotRow][column] != 0 { pivotRow += 1 }
        }
        for i in pivotRow..<a.count {
            let residual = abs(b[i] - b[i].rounded())
            if residual > tolerance * weight[i] { return false }
        }
        return true
    }
}

extension CrystalStructure {
    /// The space-group type of this structure's operations, as spglib names it (`SpaceGroupType.of`), read in the
    /// conventional setting when there is one (a file's reduced or primitive cell is in no standard basis), else in the
    /// structure's own; nil when neither is in the basis of one of spglib's 530 settings. Known limit: the conventional
    /// cell is chosen from the lattice metric alone, so a cell whose metric is more symmetric than its atoms (an fcc
    /// cell strained along c and written in its cubic axes; an mm2 structure on a cube) can be in no setting's basis
    /// and reads nil, where spglib, which takes the axes from the group, names it (docs/parity.md). This is what the atoms show;
    /// `spaceGroupNumber` is what the source stated, and the two can differ (a relaxed cell reads lower): show both.
    public var spaceGroupType: SpaceGroupType? {
        let tolerance = max(1e-3, 10 * translationUncertainty)
        if let conventional = conventionalSetting,
           let type = SpaceGroupType.of(conventional.operations,
                                        translationTolerance: max(tolerance, 10 * conventional.translationUncertainty)) {
            return type
        }
        return SpaceGroupType.of(operations, translationTolerance: tolerance)
    }
}
