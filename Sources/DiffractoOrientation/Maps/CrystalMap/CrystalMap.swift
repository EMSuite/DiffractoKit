//
//  CrystalMap.swift
//  Moved from Electron Microscopy Studio's Sources/TEMSEM/CrystalMap (at 7b09bf5) by DiffractoKit#180; the bundles its moved
//  tests read are in Tests/DiffractoOrientationTests/Fixtures/orix-sem, the others stay with the app (docs/parity.md).
//  orix 0.15.0's `CrystalMap` (orix.crystal_map.crystal_map): one orientation map for every reader and, later, every indexing
//  result. Rotations (unit quaternions, lab to crystal) and phase ids per point, the `PhaseList`, properties by name, x and y
//  with their steps, the scan unit, and orix's `is_in_data` mask: slicing by phase, by a boolean mask or by rows and columns
//  keeps the whole map and hides points, so a slice keeps its place in the original grid (`row`, `col`, `mapData`). The
//  `not_indexed` phase (id −1) is added when any point carries −1. GPL-3.0 (orix is GPL-3.0-or-later), see NOTICE.
//
//  Parity: `tools/parity/sem` bundles `m3/crystal-map-sdss` (shape, steps, phases, slicing, map data against orix, exact) and
//  `m3/ti-neighbour-misorientations` / `m3/ti-neighbour-misorientations-full` (figure F7), docs/parity.md.
//

import Foundation

public nonisolated struct CrystalMap: Sendable, Equatable {
    public enum Refusal: Error, Equatable, Sendable, LocalizedError {
        /// A point count that differs from the length of another per-point array.
        case lengthMismatch(String)
        /// `orientations` on points of more than one phase (orix's ValueError).
        case moreThanOnePhase([String])
        /// A neighbour operation on points that do not fill their rectangle.
        case notARectangle(shape: [Int], points: Int)
        /// A phase list with fewer phases than the data's phase ids (orix's `phases` setter).
        case tooFewPhases
        /// x and y that make no grid (`checkedSlices`, `maximumCells`): why, as a clause.
        case coordinates(String)
        /// A plain-array map given both Euler angles and quaternions, or neither.
        case orientations

        public var errorDescription: String? {
            switch self {
            case .lengthMismatch(let what): "The crystal map's \(what) has another length than its rotations."
            case .moreThanOnePhase(let names): "The points hold the phases \(names.joined(separator: ", ")); this needs one phase."
            case .notARectangle(let shape, let points): "The \(points) points do not fill a \(shape.map(String.init).joined(separator: " × ")) map."
            case .tooFewPhases: "There must be at least as many phases as there are phase ids."
            case .coordinates(let why): "The crystal map's coordinates make no grid: \(why)."
            case .orientations: "A crystal map takes Euler angles or quaternions per point, one of the two."
            }
        }
    }

    /// The most grid cells a map may span, 2^28 (2 GB for one Double per cell), far above any real scan. orix has no bound and
    /// allocates whatever the coordinates ask for; a map past it is refused instead.
    public static let maximumCells = 1 << 28

    // MARK: Stored (every point of the original map; `isInData` hides points)

    /// Unit quaternions (a, b, c, d), one per point.
    public private(set) var rotations: [SIMD4<Double>]
    /// Bunge Euler angles (phi1, Phi, phi2) in radians per point: the file's values where a reader gave them, else from the
    /// rotations (`qu2eu`).
    public private(set) var euler: [SIMD3<Double>]
    public private(set) var allPhaseIDs: [Int]
    public var phases: PhaseList
    /// Every other per-point column (scores, fit, confidence, ...), by its name in the file, one value per point.
    public private(set) var properties: [String: [Double]]
    /// x and y of every point in the scan unit; nil where the map has no such axis.
    public private(set) var allX: [Double]?
    public private(set) var allY: [Double]?
    /// orix's `dx` and `dy` (from every point, in or out of data): `stepSize` of `allX` and `allY`, 0 without the axis. Set once
    /// with the coordinates, which nothing changes after, so an access is not a pass over the map (#624).
    public private(set) var dx: Double
    public private(set) var dy: Double
    public var scanUnit: String
    public private(set) var isInData: [Bool]
    /// The shape of the whole map from the coordinates, rows then columns (orix's `_original_shape`).
    public private(set) var originalShape: [Int]

    // MARK: Construction

    /// orix's `CrystalMap.__init__`. `phaseIDs` nil is one phase (id 0); `x` and `y` both nil is a 1D map with x = 0, 1, 2, ...;
    /// a `phaseList` is matched to the ids in the data as orix matches it (superfluous phases dropped from the top, missing ones
    /// added with a free colour, then the list re-keyed to the data's ids in order); a −1 in the ids adds `not_indexed`.
    public init(rotations: [SIMD4<Double>], euler: [SIMD3<Double>]? = nil, phaseIDs: [Int]? = nil, x: [Double]? = nil, y: [Double]? = nil,
                 phaseList: PhaseList? = nil, properties: [String: [Double]] = [:], scanUnit: String = "px",
                 isInData: [Bool]? = nil) throws {
        let n = rotations.count
        let ids = phaseIDs ?? Array(repeating: 0, count: n)
        guard ids.count == n else { throw Refusal.lengthMismatch("phase ids") }
        if let euler, euler.count != n { throw Refusal.lengthMismatch("Euler angles") }
        if let x, x.count != n { throw Refusal.lengthMismatch("x") }
        if let y, y.count != n { throw Refusal.lengthMismatch("y") }
        if let isInData, isInData.count != n { throw Refusal.lengthMismatch("mask") }
        for (name, values) in properties where values.count != n { throw Refusal.lengthMismatch("property \(name)") }

        self.rotations = rotations.map { Rotation($0).q }
        self.euler = euler ?? rotations.map { RotationConversions.qu2eu(Rotation($0).q) }
        self.allPhaseIDs = ids
        self.allX = (x == nil && y == nil) ? (0..<n).map(Double.init) : x
        self.allY = y
        self.dx = self.allX.map(CrystalMap.stepSize) ?? 0
        self.dy = y.map(CrystalMap.stepSize) ?? 0
        self.properties = properties
        self.scanUnit = scanUnit
        self.isInData = isInData ?? Array(repeating: true, count: n)

        var unique = Array(Set(ids)).sorted()
        let includeNotIndexed = unique.first == PhaseList.notIndexedID
        if includeNotIndexed { unique.removeFirst() }
        if var list = phaseList {
            list.remove(id: PhaseList.notIndexedID)
            var different = list.count - unique.count
            if different > 0 {
                for id in list.ids.reversed() where different > 0 && !unique.contains(id) {
                    list.remove(id: id)
                    different -= 1
                }
            } else if different < 0 {
                var used = list.colors
                var entries: [(id: Int, phase: CrystalPhase)] = []
                for id in unique {
                    if let p = list[id: id] { entries.append((id, p)) } else {
                        let color = PhaseList.tableau.first { !used.contains($0) } ?? "tab:blue"
                        used.append(color)
                        entries.append((id, CrystalPhase(name: "", color: color)))
                    }
                }
                list = PhaseList(entries)
            }
            // Re-keyed to the data's ids in order (orix's `dict(zip(new_ids, values))`).
            self.phases = PhaseList(zip(unique, list.phases).map { (id: $0.0, phase: $0.1) })
        } else {
            self.phases = PhaseList.unnamed(ids: unique)
        }
        if includeNotIndexed { self.phases.addNotIndexed() }
        self.originalShape = []
        // orix's `_original_shape` (`shape(onlyInData: false)`), from coordinates that make a grid of at most `maximumCells` cells.
        let extents = try CrystalMap.checkedSlices([("y", allY ?? []), ("x", allX ?? [])]).map { $0.stop - $0.start }
        var cells = 1
        for extent in extents {
            let product = cells.multipliedReportingOverflow(by: extent)
            guard !product.overflow, product.partialValue <= CrystalMap.maximumCells else {
                throw Refusal.coordinates("they span more than \(CrystalMap.maximumCells) grid cells")
            }
            cells = product.partialValue
        }
        self.originalShape = extents
    }

    /// A reader's map: Bunge Euler angles in radians (lab to crystal, orix's `Rotation.from_euler`) kept as given; `isInData` hides
    /// the points a file keeps outside the map (orix's `is_in_data`).
    public init(euler: [SIMD3<Double>], phaseIDs: [Int]? = nil, x: [Double]? = nil, y: [Double]? = nil, phaseList: PhaseList? = nil,
                 properties: [String: [Double]] = [:], scanUnit: String = "px", isInData: [Bool]? = nil) throws {
        try self.init(rotations: euler.map { RotationConversions.eu2qu($0) }, euler: euler, phaseIDs: phaseIDs, x: x, y: y,
                      phaseList: phaseList, properties: properties, scanUnit: scanUnit, isInData: isInData)
    }

    /// A reader's map from plain arrays, as a file reader that knows no `CrystalMap` hands it over: Bunge Euler angles in radians
    /// (lab to crystal) or unit quaternions, one per point (exactly one of the two); the phase ids, x, y and the named per-point
    /// properties as in `init(rotations:)`; and the file's phase headers, each by its id with its name, formula, space group,
    /// point group (orix's name), lattice (a, b, c in `latticeUnit`, alpha, beta, gamma in degrees), lattice unit
    /// (`CrystalPhase.LatticeUnit`'s raw value, nil or an unknown name where the file states none) and matplotlib colour (nil for
    /// orix's default, the next free Tableau colour). The headers become a `PhaseList` and are matched to the ids as orix matches
    /// them; the map is the one `init(euler:)` or `init(rotations:)` makes from that list.
    public init(euler: [SIMD3<Double>]? = nil, quaternions: [SIMD4<Double>]? = nil, phaseIDs: [Int]? = nil, x: [Double]? = nil,
                y: [Double]? = nil, properties: [String: [Double]] = [:],
                phaseHeaders: [(id: Int, name: String, formula: String?, spaceGroup: Int?, pointGroup: String?, lattice: [Double]?,
                                latticeUnit: String?, color: String?)],
                scanUnit: String = "px", isInData: [Bool]? = nil) throws {
        let entries: [(id: Int, phase: CrystalPhase)] = phaseHeaders.map { h in
            (id: h.id, phase: CrystalPhase(name: h.name, formula: h.formula, spaceGroup: h.spaceGroup, pointGroup: h.pointGroup,
                                           lattice: h.lattice, color: h.color ?? "tab:blue",
                                           latticeUnit: h.latticeUnit.flatMap { CrystalPhase.LatticeUnit(rawValue: $0) }))
        }
        let list = PhaseList(entries)
        switch (euler, quaternions) {
        case (let euler?, nil):
            try self.init(euler: euler, phaseIDs: phaseIDs, x: x, y: y, phaseList: list, properties: properties, scanUnit: scanUnit,
                          isInData: isInData)
        case (nil, let quaternions?):
            try self.init(rotations: quaternions, phaseIDs: phaseIDs, x: x, y: y, phaseList: list, properties: properties,
                          scanUnit: scanUnit, isInData: isInData)
        default:
            throw Refusal.orientations
        }
    }

    /// orix's `create_coordinate_arrays`: x and y of a row-major map of `shape` (rows, columns) with `steps` (dy, dx).
    public static func coordinates(shape: (rows: Int, columns: Int), steps: (dy: Double, dx: Double) = (1, 1)) -> (x: [Double], y: [Double]) {
        var x: [Double] = [], y: [Double] = []
        x.reserveCapacity(shape.rows * shape.columns); y.reserveCapacity(shape.rows * shape.columns)
        for r in 0..<shape.rows {
            for c in 0..<shape.columns {
                x.append(Double(c) * steps.dx)
                y.append(Double(r) * steps.dy)
            }
        }
        return (x, y)
    }

    // MARK: Points in data

    /// Every point of the map, in data or not (the readers' row count).
    public var pointCount: Int { allPhaseIDs.count }

    /// The indices of the points in data (orix's `id`).
    public var ids: [Int] { isInData.indices.filter { isInData[$0] } }
    /// orix's `size`.
    public var size: Int { isInData.reduce(0) { $0 + ($1 ? 1 : 0) } }
    public var phaseIDs: [Int] { inData(allPhaseIDs) }
    public var rotationsInData: [SIMD4<Double>] { inData(rotations) }
    public var eulerInData: [SIMD3<Double>] { inData(euler) }
    public func property(_ name: String) -> [Double]? { properties[name].map { inData($0) } }
    public var propertyNames: [String] { properties.keys.sorted() }

    /// orix's `x`: nil when there are none or all are equal (a step of 0: the coordinates are finite, `checkedSlices`, so two
    /// distinct values are a step above 0).
    public var x: [Double]? {
        guard let allX, dx != 0 else { return nil }
        return inData(allX)
    }

    public var y: [Double]? {
        guard let allY, dy != 0 else { return nil }
        return inData(allY)
    }

    /// orix's `shape`: rows then columns of the points in data, from their coordinates.
    public var shape: [Int] { shape(onlyInData: true) }

    /// orix's `is_indexed` per point in data.
    public var isIndexed: [Bool] { phaseIDs.map { $0 != PhaseList.notIndexedID } }
    public var allIndexed: Bool { !phaseIDs.contains(PhaseList.notIndexedID) }

    /// orix's `phases_in_data`.
    public var phasesInData: PhaseList { phases.restricted(to: Set(phaseIDs)) }

    /// orix's `row` and `col`: each point's row and column in the original grid, less the smallest of the points in data.
    public var row: [Int] { gridIndices(rows: true) }
    public var col: [Int] { gridIndices(rows: false) }

    /// orix's `orientations`: the rotations in data with the symmetry of their one phase (C1 where it names none).
    public func orientations() throws -> Orientations {
        let inPhases = phasesInData
        guard inPhases.count == 1 else { throw Refusal.moreThanOnePhase(inPhases.names) }
        return Orientations(quaternions: rotationsInData, symmetry: inPhases.phases[0].pointGroup ?? PointGroup.c1)
    }

    // MARK: Slicing (orix's `__getitem__`: a copy with a narrower mask)

    /// The points of the phases named, or of every indexed phase for "indexed".
    public func points(ofPhases names: [String]) -> CrystalMap {
        var keep = Array(repeating: false, count: size)
        let pids = phaseIDs
        for name in names {
            for (id, phase) in phases.entries where phase.name == name {
                for i in pids.indices where pids[i] == id { keep[i] = true }
            }
            if name.lowercased() == "indexed" {
                for i in pids.indices where pids[i] != PhaseList.notIndexedID { keep[i] = true }
            }
        }
        return masked(keep)
    }

    public subscript(phase name: String) -> CrystalMap { points(ofPhases: [name]) }

    /// The points in data where `mask` (one entry per point in data) is true.
    public func masked(_ mask: [Bool]) -> CrystalMap {
        precondition(mask.count == size, "a mask has one entry per point in data")
        var copy = self
        var next = Array(repeating: false, count: isInData.count)
        for (k, i) in ids.enumerated() where mask[k] { next[i] = true }
        copy.isInData = next
        return copy
    }

    /// orix's slicing by rows and columns of the current data extent (`xmap[rows, columns]`); a 1D map slices its one axis by
    /// `rows`. Ranges are clamped to the extent as numpy clamps slices. Unlike orix, which re-shows hidden points inside the
    /// window (it overwrites the mask there), a point hidden before stays hidden: a phase slice cut to rows keeps its phase.
    public func sliced(rows: Range<Int>? = nil, columns: Range<Int>? = nil) -> CrystalMap {
        let current = shape
        let extents = slices(onlyInData: true)
        let base = slices(onlyInData: false)
        let original = originalShape
        var copy = self
        var next = isInData
        func clamp(_ r: Range<Int>?, _ n: Int) -> Range<Int> {
            guard let r else { return 0..<n }
            let lo = min(max(r.lowerBound, 0), n), hi = min(max(r.upperBound, lo), n)
            return lo..<hi
        }
        if original.count == 2, current.count == 2, extents.count == 2, base.count == 2 {
            let rr = clamp(rows, current[0]), cc = clamp(columns, current[1])
            for i in next.indices where next[i] {
                let r = i / original[1] + base[0].start - extents[0].start
                let c = i % original[1] + base[1].start - extents[1].start
                next[i] = rr.contains(r) && cc.contains(c)
            }
        } else if current.count == 1, extents.count == 1, base.count == 1 {
            let rr = clamp(rows, current[0])
            for i in next.indices where next[i] { next[i] = rr.contains(i + base[0].start - extents[0].start) }
        }
        copy.isInData = next
        return copy
    }

    // MARK: Map data

    /// orix's `get_map_data(item)` for a per-point array of the points in data: the data extent's rows x columns (row-major),
    /// points not in data set to `fill` (NaN as orix's default).
    public func mapData(_ values: [Double], fill: Double = .nan) -> [Double] {
        precondition(values.count == size, "map data takes one value per point in data")
        var full = Array(repeating: fill, count: isInData.count)
        for (k, i) in ids.enumerated() { full[i] = values[k] }
        return cropToExtent(full)
    }

    /// `get_map_data(name)` of a property or of "phase_id" (as Double).
    public func mapData(named name: String, fill: Double = .nan) -> [Double]? {
        if name == "phase_id" { return mapData(phaseIDs.map(Double.init), fill: fill) }
        return property(name).map { mapData($0, fill: fill) }
    }

    // MARK: Helpers

    private func inData<T>(_ values: [T]) -> [T] {
        var out: [T] = []
        out.reserveCapacity(values.count)
        for i in values.indices where isInData[i] { out.append(values[i]) }
        return out
    }

    /// The coordinate arrays orix's shape uses, y first then x, of the points in data or of all points.
    private func coordinateArrays(onlyInData: Bool) -> [[Double]?] {
        if onlyInData { return [y, x] }
        return [allY, allX]
    }

    /// orix `_data_slices_from_coordinates`: per axis with coordinates and a non-zero step, `round(min / step) ..< round(max /
    /// step + 1)` (numpy's round half to even).
    public func slices(onlyInData: Bool) -> [(start: Int, stop: Int)] {
        CrystalMap.dataSlices(coordinates: coordinateArrays(onlyInData: onlyInData).map { $0 ?? [] }, steps: [dy, dx])
    }

    /// orix `_data_slices_from_coordinates` on coordinate arrays and their steps (an empty array or a zero step is skipped).
    public static func dataSlices(coordinates: [[Double]], steps: [Double]) -> [(start: Int, stop: Int)] {
        var out: [(start: Int, stop: Int)] = []
        for (coords, step) in zip(coordinates, steps) {
            guard let lo = coords.min(), let hi = coords.max(), step != 0 else { continue }
            out.append((Int((lo / step).rounded(.toNearestOrEven)), Int((hi / step + 1).rounded(.toNearestOrEven))))
        }
        return out
    }

    /// `dataSlices` of named coordinate arrays with their own steps (`stepSize`), refused where its arithmetic cannot run: a value
    /// that is not finite, or a slice bound past 2^53 steps (a step far finer than the spread, as one jittered coordinate makes).
    /// Coordinates that pass give `dataSlices`'s slices unchanged.
    public static func checkedSlices(_ axes: [(name: String, values: [Double])]) throws -> [(start: Int, stop: Int)] {
        for axis in axes where !axis.values.allSatisfy(\.isFinite) {
            throw Refusal.coordinates("\(axis.name) holds a value that is not a finite number")
        }
        let steps = axes.map { stepSize($0.values) }
        for (axis, step) in zip(axes, steps) {
            guard let lo = axis.values.min(), let hi = axis.values.max(), step != 0 else { continue }
            guard abs((lo / step).rounded(.toNearestOrEven)) < 0x1p53, abs((hi / step + 1).rounded(.toNearestOrEven)) < 0x1p53 else {
                throw Refusal.coordinates("\(axis.name) spans more than 2^53 times its smallest step")
            }
        }
        return dataSlices(coordinates: axes.map { $0.values }, steps: steps)
    }

    private func shape(onlyInData: Bool) -> [Int] { slices(onlyInData: onlyInData).map { $0.stop - $0.start } }

    /// orix `_step_size_from_coordinates`: `unique[1] - unique[0]` of the sorted unique values, 0 for fewer than two. The two
    /// smallest distinct values are found in one pass, without the set and the sort (#624); a NaN is passed over (a map's
    /// coordinates are finite, `checkedSlices`), and −0 and 0 are one value, as in the set.
    public static func stepSize(_ coordinates: [Double]) -> Double {
        var lowest = 0.0, second = 0.0, distinct = 0
        for v in coordinates where !v.isNaN {
            switch distinct {
            case 0:
                lowest = v; distinct = 1
            case 1:
                if v < lowest { second = lowest; lowest = v; distinct = 2 } else if v > lowest { second = v; distinct = 2 }
            default:
                if v < lowest { second = lowest; lowest = v } else if v > lowest, v < second { second = v }
            }
        }
        return distinct == 2 ? second - lowest : 0
    }

    private func gridIndices(rows: Bool) -> [Int] {
        var grid = originalShape
        if grid.count == 1 { grid = allX == nil ? [grid[0], 1] : [1, grid[0]] }
        guard grid.count == 2, grid[1] > 0 else { return [] }
        let values = ids.map { rows ? $0 / grid[1] : $0 % grid[1] }
        guard let low = values.min() else { return [] }
        return values.map { $0 - low }
    }

    /// The whole map's row-major array cut to the current data extent (orix's `reshaped_array[slices]`, offset by where the
    /// whole map starts). A cell past the last point is NaN, with one axis as with two: a scan stopped mid-row, or a one-row map
    /// whose coordinates leave a gap, spans more cells than it has points, and the array keeps the grid's length.
    private func cropToExtent(_ full: [Double]) -> [Double] {
        let original = originalShape
        let extents = slices(onlyInData: true)
        let base = slices(onlyInData: false)
        guard original.count == 2, extents.count == 2, base.count == 2 else {
            guard let e = extents.first, let b = base.first else { return full }
            let start = max(e.start - b.start, 0), stop = max(e.stop - b.start, start)
            return (start..<stop).map { $0 < full.count ? full[$0] : .nan }
        }
        var out: [Double] = []
        out.reserveCapacity((extents[0].stop - extents[0].start) * (extents[1].stop - extents[1].start))
        for r in (extents[0].start - base[0].start)..<(extents[0].stop - base[0].start) {
            for c in (extents[1].start - base[1].start)..<(extents[1].stop - base[1].start) {
                let i = r * original[1] + c
                out.append(i < full.count ? full[i] : .nan)
            }
        }
        return out
    }
}
