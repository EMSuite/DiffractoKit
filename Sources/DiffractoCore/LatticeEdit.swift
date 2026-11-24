//
//  LatticeEdit.swift
//  Role: change the lattice parameters of a structure within its crystal system. The parameters that may change, and
//        those that move with them, are read from the structure's own Laue group in its own cell; the fractional sites
//        are kept; the symmetry is read again on the new cell, and an edit that would change the Laue class is refused
//        with its reason. No Python function does this step; the d-spacings of the edited cells are pinned to pymatgen.
//

import Foundation
import simd

/// One of the six cell parameters.
public enum LatticeParameter: String, Sendable, Equatable, Hashable, CaseIterable, CustomStringConvertible {
    case a, b, c, alpha, beta, gamma

    public var description: String {
        switch self {
        case .a, .b, .c: rawValue
        case .alpha: "α"
        case .beta: "β"
        case .gamma: "γ"
        }
    }

    var isLength: Bool { self == .a || self == .b || self == .c }

    /// The axis of a length (0, 1, 2), or the axis opposite an angle (α is between b and c: 0).
    var index: Int {
        switch self {
        case .a, .alpha: 0
        case .b, .beta: 1
        case .c, .gamma: 2
        }
    }

    static func length(_ axis: Int) -> LatticeParameter { [.a, .b, .c][axis] }
    static func angle(opposite axis: Int) -> LatticeParameter { [.alpha, .beta, .gamma][axis] }
}

/// Why a lattice edit was refused. Each case reads as a sentence for a scientist.
public enum LatticeEditError: Error, Sendable, Equatable, CustomStringConvertible {
    /// The parameter is fixed by the symmetry, or moves with another one that is the one to edit.
    case notFree(LatticeParameter, free: [LatticeParameter], system: CrystalSystem?)
    /// A length that is not a positive finite number, or an angle outside (0°, 180°).
    case invalidValue(LatticeParameter, Double)
    /// The six parameters do not make a cell (`UnitCell.issues`).
    case invalidCell(String)
    /// The symmetry read on the new cell has another Laue class.
    case symmetryChanged(from: LaueClass?, to: LaueClass?)
    /// A rotation of the given operations does not keep the new cell's metric.
    case operationsDoNotFit
    /// The structure's symmetry was not read (`symmetryRefusal`, or operations that are not a group): there is no
    /// crystal system to edit within.
    case symmetryNotRead

    public var description: String {
        switch self {
        case .notFree(let p, let free, let system):
            let names = free.map(\.description).joined(separator: ", ")
            let kind = system.map { "this \($0) structure" } ?? "this structure"
            return "\(p) cannot be set on its own in \(kind): the symmetry fixes it or ties it to another parameter. "
                + "The parameters that can be set are \(names.isEmpty ? "none" : names)."
        case .invalidValue(let p, let v):
            return p.isLength
                ? "The length \(p) must be a positive number of ångströms; \(v) is not."
                : "The angle \(p) must lie between 0° and 180°; \(v) does not."
        case .invalidCell(let reason):
            return "These parameters do not make a cell: \(reason)"
        case .symmetryChanged(let from, let to):
            return "With these parameters the cell reads the Laue class \(to?.symbol ?? "of no group"), not "
                + "\(from?.symbol ?? "the one it had"): the edit would change the crystal's symmetry, so it is refused."
        case .operationsDoNotFit:
            return "The structure's symmetry operations do not fit a cell with these parameters, so the edit is refused."
        case .symmetryNotRead:
            return "The symmetry of this structure was not read, so there is no crystal system to keep: the edit is refused."
        }
    }
}

extension CrystalStructure {

    /// The parameters of `cell` that can be set without changing the crystal system, in the order a, b, c, α, β, γ:
    /// one per group of parameters the symmetry ties together (cubic: a; tetragonal and hexagonal: a and c; trigonal
    /// in hexagonal axes: a and c, in rhombohedral axes: a and α; orthorhombic: a, b, c; monoclinic: a, b, c and the
    /// angle at the unique axis; triclinic: all six). Read from the Laue group's rotations in this cell: two axes a
    /// rotation exchanges have one length, and an angle is free when every rotation carries its two axes onto two axes
    /// (up to sign) without forcing it to 90°; each parameter offered is then checked, moved alone, against every
    /// rotation. In a primitive cell of a centred lattice, where rotations carry axes onto sums of axes, this offers
    /// less than the crystal system allows (often one scale of all three lengths); edit the conventional cell for the
    /// rest. A structure whose symmetry was not read (`symmetryRefusal`) offers none.
    public var freeLatticeParameters: [LatticeParameter] {
        guard symmetryRefusal == nil, laueClass != nil else { return [] }
        let ties = LatticeTies(rotations: laueGroup.rotations, cell: cell)
        return LatticeParameter.allCases.filter { ties.isRepresentative($0) }
    }

    /// The structure with some of its free cell parameters set (`freeLatticeParameters`); the parameters tied to them
    /// follow (equal lengths take the value; where the free length is one scale of three unequal lengths, they keep
    /// their ratios), the fixed ones keep their values bit for bit, the fractional sites are kept. The symmetry is read again
    /// on the new cell with the same settings (`replacingCell`). Refused: a parameter that is not free, a value that
    /// is not a length or an angle, parameters that make no cell, and an edit that changes the Laue class.
    public func editingLattice(_ values: [LatticeParameter: Double]) throws(LatticeEditError) -> CrystalStructure {
        guard symmetryRefusal == nil, laueClass != nil else { throw .symmetryNotRead }
        let ties = LatticeTies(rotations: laueGroup.rotations, cell: cell)
        let free = LatticeParameter.allCases.filter { ties.isRepresentative($0) }
        for p in LatticeParameter.allCases {
            guard let v = values[p] else { continue }
            guard free.contains(p) else { throw .notFree(p, free: free, system: crystalSystem) }
            if p.isLength {
                guard v.isFinite, v > 0 else { throw .invalidValue(p, v) }
            } else {
                guard v.isFinite, v > 0, v < 180 else { throw .invalidValue(p, v) }
            }
        }
        var lengths = [cell.a, cell.b, cell.c]
        var angles = [cell.alphaDeg, cell.betaDeg, cell.gammaDeg]
        // The lengths of a group move together: equal lengths take the value, and the unequal lengths of a scale (a
        // primitive cell of a centred lattice) keep their ratios to it.
        let old = lengths
        for axis in 0..<3 {
            let representative = ties.lengthRepresentative(axis)
            guard let v = values[representative] else { continue }
            let reference = old[representative.index]
            lengths[axis] = old[axis] == reference ? v : old[axis] * (v / reference)
        }
        for axis in 0..<3 {
            guard let tie = ties.angleRepresentative(axis), let v = values[tie.parameter] else { continue }
            angles[axis] = tie.sign > 0 ? v : 180 - v
        }
        let edited = UnitCell(a: lengths[0], b: lengths[1], c: lengths[2],
                              alphaDeg: angles[0], betaDeg: angles[1], gammaDeg: angles[2])
        return try replacingCell(edited)
    }

    /// The structure on another cell with the same fractional sites, its symmetry read again: by the same detection
    /// with the same settings (`detectionSettings`) for a detected structure, and for given operations (or the
    /// conventional-setting structure) the same operations, refused when one of their rotations does not keep the new
    /// metric. Refused when the cell is not valid or the Laue class changes.
    public func replacingCell(_ newCell: UnitCell) throws(LatticeEditError) -> CrystalStructure {
        if let issue = newCell.issues.first { throw .invalidCell(issue.description) }
        let rebuilt: CrystalStructure
        if setting == .file, let s = detectionSettings {
            rebuilt = CrystalStructure.detectingSymmetry(
                cell: newCell, sites: sites, fileOperations: fileOperations, spaceGroupNumber: spaceGroupNumber,
                name: name, coarsestHalfStep: s.coarsestHalfStep, metricLengthTolerance: s.metricLengthTolerance,
                metricAngleToleranceDeg: s.metricAngleToleranceDeg, tolerance: s.tolerance,
                conventional: s.conventional, policy: s.policy)
        } else {
            guard LatticeTies.keepsMetric(laueGroup.rotations, of: newCell, as: cell) else { throw .operationsDoNotFit }
            if setting == .conventional {
                rebuilt = CrystalStructure(conventionalCell: newCell, sites: sites, operations: operations,
                                           laueGroup: laueGroup, laueClass: laueClass,
                                           translationUncertainty: translationUncertainty)
            } else {
                rebuilt = CrystalStructure(cell: newCell, sites: sites, operations: operations,
                                           spaceGroupNumber: spaceGroupNumber, name: name)
            }
        }
        guard rebuilt.laueClass == laueClass else { throw .symmetryChanged(from: laueClass, to: rebuilt.laueClass) }
        // The class of the conventional setting must hold too. The conventional cell is read from the metric alone, so a
        // centring that only the atoms show is not seen: a hexagonal crystal in an orthohexagonal C cell reads mmm in
        // both, and an edit of a alone is an edit within mmm (CI, 24 Nov).
        let before = conventionalSetting?.laueClass, after = rebuilt.conventionalSetting?.laueClass
        guard before == after else { throw .symmetryChanged(from: before, to: after) }
        return rebuilt
    }
}

/// Which cell parameters a Laue group ties together, in the cell its rotations are written in (x' = R x, so column i
/// of R is the image of axis i).
struct LatticeTies {
    /// The representative axis of each length (the lowest axis of its group).
    private var lengthRoot: [Int]
    /// For each angle (by the axis opposite it): nil when fixed, else its representative angle and the sign that
    /// relates their cosines (+1 equal angles, -1 supplementary).
    private var angleTie: [(root: Int, sign: Int)?]

    init(rotations: [IntegerMatrix], cell: UnitCell) {
        // The image of axis i under R: (axis j, sign) when column i is ±e_j, nil otherwise.
        func image(_ r: IntegerMatrix, _ i: Int) -> (axis: Int, sign: Int)? {
            let column = [r.rows[0][i], r.rows[1][i], r.rows[2][i]]
            let nonzero = column.indices.filter { column[$0] != 0 }
            guard nonzero.count == 1, abs(column[nonzero[0]]) == 1 else { return nil }
            return (nonzero[0], column[nonzero[0]])
        }

        // Lengths: union of the axes a rotation exchanges.
        var parent = [0, 1, 2]
        func find(_ x: Int) -> Int { var x = x; while parent[x] != x { x = parent[x] }; return x }
        for r in rotations {
            for i in 0..<3 {
                if let target = image(r, i) {
                    let ri = find(i), rj = find(target.axis)
                    if ri != rj { parent[max(ri, rj)] = min(ri, rj) }
                }
            }
        }
        lengthRoot = (0..<3).map { find($0) }

        // Angles: the angle opposite axis k is between the other two axes. Each rotation either carries that pair onto
        // another pair (a tie, with the product of the signs) or not (the angle is fixed). A rotation that carries the
        // pair onto itself with a sign product of -1 forces the cosine to 0: fixed.
        func pair(_ k: Int) -> (Int, Int) { k == 0 ? (1, 2) : (k == 1 ? (0, 2) : (0, 1)) }
        func opposite(_ i: Int, _ j: Int) -> Int { 3 - i - j }
        var fixed = [false, false, false]
        var angleParent = [0, 1, 2]
        var angleSign = [1, 1, 1]   // sign relative to its parent
        func angleFind(_ x: Int) -> (root: Int, sign: Int) {
            var x = x, sign = 1
            while angleParent[x] != x { sign *= angleSign[x]; x = angleParent[x] }
            return (x, sign)
        }
        for r in rotations {
            for k in 0..<3 {
                let (i, j) = pair(k)
                guard let first = image(r, i), let second = image(r, j), first.axis != second.axis else {
                    fixed[k] = true
                    continue
                }
                let target = opposite(first.axis, second.axis)
                let sign = first.sign * second.sign
                let fromK = angleFind(k), fromTarget = angleFind(target)
                let rk = fromK.root, sk = fromK.sign, rt = fromTarget.root, st = fromTarget.sign
                if rk == rt {
                    // Same group: the relation must agree, else the cosine is forced to 0.
                    if sk * st != sign { fixed[k] = true }
                } else {
                    let (low, high) = (min(rk, rt), max(rk, rt))
                    angleParent[high] = low
                    // cos(high) relates to cos(low) by the product of the signs along both paths and the rotation's.
                    angleSign[high] = sk * st * sign
                }
            }
        }
        // A fixed angle fixes every angle it is tied to.
        var groupFixed = [false, false, false]
        for k in 0..<3 where fixed[k] { groupFixed[angleFind(k).root] = true }
        angleTie = (0..<3).map { k -> (root: Int, sign: Int)? in
            let found = angleFind(k)
            return groupFixed[found.root] ? nil : (root: found.root, sign: found.sign)
        }

        // The rule reads each rotation's action on axes; a rotation that carries an axis onto a sum of axes (a primitive
        // cell of a centred lattice) also ties lengths to angles. Each offered parameter is checked: moved alone, with
        // its ties, every rotation must still keep the metric. A length that fails gives way to one scale of all three;
        // an angle that fails is fixed.
        let base = [cell.a, cell.b, cell.c, cell.alphaDeg, cell.betaDeg, cell.gammaDeg]
        let roots = lengthRoot, readTies = angleTie
        func moved(_ p: LatticeParameter) -> UnitCell {
            var q = base
            if p.isLength {
                for axis in 0..<3 where roots[axis] == roots[p.index] { q[axis] = base[axis] * 1.01 }
            } else if let tie = readTies[p.index] {
                let v = base[3 + tie.root] + 1
                for k in 0..<3 {
                    if let t = readTies[k], t.root == tie.root { q[3 + k] = t.sign > 0 ? v : 180 - v }
                }
            }
            return UnitCell(a: q[0], b: q[1], c: q[2], alphaDeg: q[3], betaDeg: q[4], gammaDeg: q[5])
        }
        let lengthsHold = (0..<3).filter { roots[$0] == $0 }
            .allSatisfy { LatticeTies.keepsMetric(rotations, of: moved(.length($0)), as: cell) }
        if !lengthsHold { lengthRoot = [0, 0, 0] }
        var checkedTies = readTies
        for k in 0..<3 {
            guard let tie = readTies[k], tie.root == k else { continue }
            if !LatticeTies.keepsMetric(rotations, of: moved(.angle(opposite: k)), as: cell) {
                for m in 0..<3 where readTies[m]?.root == k { checkedTies[m] = nil }
            }
        }
        angleTie = checkedTies
    }

    /// The largest |RᵀGR − G|ᵢⱼ / √(Gᵢᵢ Gⱼⱼ) over the rotations: 0 for exact symmetries, about the written precision for
    /// a cell typed to a few digits (a primitive bcc cell at 109.4712°), NaN or infinite for no cell.
    static func metricResidual(_ rotations: [IntegerMatrix], of cell: UnitCell) -> Double {
        let g = cell.metricTensor
        var worst = 0.0
        for rotation in rotations {
            let r = rotation.rows.map { SIMD3<Double>(Double($0.x), Double($0.y), Double($0.z)) }
            for i in 0..<3 {
                for j in 0..<3 {
                    var value = 0.0
                    for k in 0..<3 { for l in 0..<3 { value += r[k][i] * g[k][l] * r[l][j] } }
                    let residual = abs(value - g[i][j]) / (g[i][i] * g[j][j]).squareRoot()
                    guard residual.isFinite else { return .infinity }
                    worst = max(worst, residual)
                }
            }
        }
        return worst
    }

    /// Whether the rotations keep `cell`'s metric as well as they keep `reference`'s: within ten times the reference's
    /// residual, and never asking for better than 1e-9. A rotation that is not a symmetry of the cell is off by the
    /// size of the distortion (1e-2 for the 1 % move of the check), far above either.
    static func keepsMetric(_ rotations: [IntegerMatrix], of cell: UnitCell, as reference: UnitCell) -> Bool {
        let allowed = max(1e-9, 10 * metricResidual(rotations, of: reference))
        return metricResidual(rotations, of: cell) <= allowed
    }

    func lengthRepresentative(_ axis: Int) -> LatticeParameter { .length(lengthRoot[axis]) }

    /// The representative angle of the angle opposite `axis` and the sign of the tie; nil when it is fixed.
    func angleRepresentative(_ axis: Int) -> (parameter: LatticeParameter, sign: Int)? {
        angleTie[axis].map { (LatticeParameter.angle(opposite: $0.root), $0.sign) }
    }

    func isRepresentative(_ p: LatticeParameter) -> Bool {
        if p.isLength { return lengthRoot[p.index] == p.index }
        guard let tie = angleTie[p.index] else { return false }
        return tie.root == p.index
    }
}
