//
//  PhaseList.swift
//  Moved from Electron Microscopy Studio's Sources/TEMSEM/CrystalMap (at 7b09bf5) by DiffractoKit#180; the bundles its moved
//  tests read are in Tests/DiffractoOrientationTests/Fixtures/orix-sem, the others stay with the app (docs/parity.md).
//  orix 0.15.0's `Phase` and `PhaseList` (orix.crystal_map._phase, orix.crystal_map._phase_list): the phases of a crystal map
//  by id, ascending, each with its name, space group, point group, lattice and colour; the `not_indexed` phase (id −1, no
//  point group, white) that a map with unindexed points carries. GPL-3.0 (orix is GPL-3.0-or-later), see NOTICE.
//
//  Colours follow orix: a phase given without one takes the next of matplotlib's Tableau colours (`tab:blue`, `tab:orange`,
//  ...) not yet used, in that order; orix continues into the xkcd names after ten phases, the port cycles the ten.
//
//  On the kit's types: a phase's symmetry is the kit's `PointGroup` (orix's, `Orix/PointGroup.swift`) and its lattice the kit's
//  `CrystalLattice` (`crystalLattice`). The phase itself stays orix's `Phase`, a file's header: the kit's `PhaseModel` and
//  `PhaseDefinition` hold a `CrystalStructure` with atoms, which an EBSD map's header does not carry, and the six lattice numbers
//  keep the file's length unit, which `CrystalLattice` does not hold.
//

import Foundation

/// orix's `Phase`: one crystal phase of a map.
public nonisolated struct CrystalPhase: Sendable, Equatable {
    public var name: String
    /// The formula when the file writes one beside the name (`.ang`).
    public var formula: String?
    public var spaceGroup: Int?
    /// The point group as orix names it (`m-3m`, `6/mmm`, `622`, `432`, ...); nil for `not_indexed`.
    public var pointGroupName: String?
    /// a, b, c (the file's unit, `latticeUnit`), alpha, beta, gamma (degrees).
    public var lattice: [Double]?
    /// A matplotlib colour name (`tab:blue`, `white`).
    public var color: String
    /// The unit of a, b, c where the source states it (.ang and .ctf: Å; an EMsoft master pattern: nm); nil where it does not (orix's
    /// .h5 keeps whatever its writer had, and kikuchipy keeps EMsoft's nm), and only then is the unit guessed (macEM#614).
    public var latticeUnit: LatticeUnit?

    public enum LatticeUnit: String, Sendable, Equatable, Codable {
        case angstrom, nanometre
    }

    public init(name: String, formula: String? = nil, spaceGroup: Int? = nil, pointGroup: String? = nil, lattice: [Double]? = nil,
                 color: String = "tab:blue", latticeUnit: LatticeUnit? = nil) {
        self.name = name; self.formula = formula; self.spaceGroup = spaceGroup; self.pointGroupName = pointGroup; self.lattice = lattice
        self.color = color; self.latticeUnit = latticeUnit
    }

    /// a, b, c in ångström, alpha, beta, gamma in degrees: by the stated unit, else by the guess a source with no unit leaves (every
    /// length under 2 is nm: no crystal has a cell edge under 2 Å, and EBSD phases rarely one over 2 nm). nil without six parameters.
    public var latticeInAngstrom: [Double]? {
        guard var p = lattice, p.count == 6 else { return nil }
        let nanometres: Bool
        switch latticeUnit {
        case .angstrom: nanometres = false
        case .nanometre: nanometres = true
        case nil: nanometres = p[0] < 2 && p[1] < 2 && p[2] < 2
        }
        if nanometres { p[0] *= 10; p[1] *= 10; p[2] *= 10 }
        return p
    }

    /// The point group as an orientation symmetry; nil when the phase names none or orix knows no group of that name.
    public var pointGroup: PointGroup? { pointGroupName.flatMap(PointGroup.named) }

    /// The lattice as the kit's `CrystalLattice` (orix's `Phase.structure.lattice`): the six parameters as the file gives them when
    /// all are finite and positive, else orix's default lattice of the point group's crystal system. The file's length unit is kept
    /// (`latticeUnit`, `latticeInAngstrom`).
    public var crystalLattice: CrystalLattice {
        if let l = lattice, l.count == 6, l.allSatisfy({ $0.isFinite && $0 > 0 }) {
            return CrystalLattice(a: l[0], b: l[1], c: l[2], alpha: l[3], beta: l[4], gamma: l[5])
        }
        return CrystalLattice.defaultLattice(system: pointGroup?.system)
    }

    /// orix's `Phase.color_rgb` for the colours a phase list hands out: 0...1 per channel.
    public var colorRGB: SIMD3<Double>? { PhaseList.rgb[color] }

    /// `PhaseList.add_not_indexed`'s phase: no point group, white, orix's default lattice (1, 1, 1, 90°, 90°, 90°).
    public static let notIndexed = CrystalPhase(name: PhaseList.notIndexedName, lattice: [1, 1, 1, 90, 90, 90], color: "white")
}

/// orix's `PhaseList`: phases by id in ascending order.
public nonisolated struct PhaseList: Sendable, Equatable {
    public static let notIndexedName = "not_indexed"
    public static let notIndexedID = -1

    /// matplotlib's `TABLEAU_COLORS` in its order: the colours orix hands out first.
    public static let tableau = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink",
                                  "tab:gray", "tab:olive", "tab:cyan"]
    /// The hex values of `tableau` and `white`, as matplotlib defines them.
    public static let rgb: [String: SIMD3<Double>] = {
        let hex: [String: UInt32] = ["tab:blue": 0x1f77b4, "tab:orange": 0xff7f0e, "tab:green": 0x2ca02c, "tab:red": 0xd62728,
                                     "tab:purple": 0x9467bd, "tab:brown": 0x8c564b, "tab:pink": 0xe377c2, "tab:gray": 0x7f7f7f,
                                     "tab:olive": 0xbcbd22, "tab:cyan": 0x17becf, "white": 0xffffff, "black": 0x000000]
        return hex.mapValues { SIMD3(Double(($0 >> 16) & 0xff) / 255, Double(($0 >> 8) & 0xff) / 255, Double($0 & 0xff) / 255) }
    }()

    public private(set) var ids: [Int]
    public private(set) var phases: [CrystalPhase]

    /// Phases with their ids; sorted by id as orix's `OrderedDict(sorted(...))`. A phase whose colour is already taken gets the
    /// next free Tableau colour, as orix's constructor does.
    public init(_ entries: [(id: Int, phase: CrystalPhase)] = []) {
        var used: [String] = []
        var fixed: [(id: Int, phase: CrystalPhase)] = []
        for (id, phase) in entries {
            var p = phase
            if used.contains(p.color), let free = PhaseList.tableau.first(where: { !used.contains($0) }) { p.color = free }
            used.append(p.color)
            fixed.append((id, p))
        }
        fixed.sort { $0.id < $1.id }
        ids = fixed.map(\.id)
        phases = fixed.map(\.phase)
    }

    /// orix's `PhaseList(names=..., point_groups=..., ids=...)` path: phases without a colour take Tableau colours in order.
    public init(ids: [Int], phases: [CrystalPhase], assignColors: Bool) {
        var list: [(id: Int, phase: CrystalPhase)] = []
        for (i, (id, phase)) in zip(ids, phases).enumerated() {
            var p = phase
            if assignColors { p.color = PhaseList.tableau[i % PhaseList.tableau.count] }
            list.append((id, p))
        }
        self.init(list)
    }

    /// orix's `PhaseList(ids=unique_phase_ids)`: unnamed phases (orix's `Phase()` has name ""), no symmetry, Tableau colours.
    static func unnamed(ids: [Int]) -> PhaseList {
        PhaseList(ids: ids, phases: ids.map { _ in CrystalPhase(name: "") }, assignColors: true)
    }

    public var count: Int { ids.count }
    public var names: [String] { phases.map(\.name) }
    public var colors: [String] { phases.map(\.color) }
    public var pointGroupNames: [String?] { phases.map(\.pointGroupName) }
    public var entries: [(id: Int, phase: CrystalPhase)] { Array(zip(ids, phases)) }

    public subscript(id id: Int) -> CrystalPhase? {
        ids.firstIndex(of: id).map { phases[$0] }
    }

    /// `id_from_name`.
    public func id(named name: String) -> Int? {
        names.firstIndex(of: name).map { ids[$0] }
    }

    public func contains(id: Int) -> Bool { ids.contains(id) }

    /// `add_not_indexed`: the phase at id −1, named `not_indexed`, without point group, white.
    public mutating func addNotIndexed() {
        if let i = ids.firstIndex(of: PhaseList.notIndexedID) { phases[i] = .notIndexed; return }
        ids.insert(PhaseList.notIndexedID, at: 0)
        phases.insert(.notIndexed, at: 0)
    }

    /// Adds or replaces the phase at `id`, keeping ids ascending; a colour already used by another phase moves to the next free
    /// Tableau colour (orix's `add`).
    public mutating func set(_ phase: CrystalPhase, id: Int) {
        var p = phase
        let others = entries.filter { $0.id != id }.map(\.phase.color)
        if others.contains(p.color), let free = PhaseList.tableau.first(where: { !others.contains($0) }) { p.color = free }
        if let i = ids.firstIndex(of: id) { phases[i] = p; return }
        let at = ids.firstIndex { $0 > id } ?? ids.count
        ids.insert(id, at: at)
        phases.insert(p, at: at)
    }

    /// `__delitem__` by id.
    public mutating func remove(id: Int) {
        guard let i = ids.firstIndex(of: id) else { return }
        ids.remove(at: i)
        phases.remove(at: i)
    }

    /// The phases whose ids are in `keep`, in id order (orix's `phases[ids]`).
    func restricted(to keep: Set<Int>) -> PhaseList {
        var out = self
        for id in ids where !keep.contains(id) { out.remove(id: id) }
        return out
    }
}
