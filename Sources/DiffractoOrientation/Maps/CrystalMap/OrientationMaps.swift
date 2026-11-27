//
//  OrientationMaps.swift
//  Moved from Electron Microscopy Studio's Sources/TEMSEM/CrystalMap (at 7b09bf5) by DiffractoKit#180; the bundles its moved
//  tests read are in Tests/DiffractoOrientationTests/Fixtures/orix-sem, the others stay with the app (docs/parity.md).
//  The numbers of the pictures the app's Maps room draws of a `CrystalMap`, as orix 0.15.0's `CrystalMap.plot` makes them:
//  per-point colours (the IPF colours of each phase with `IPFColorKeyTSL` of its point group along a
//  sample direction, or each point's phase colour) laid into the map's rows and columns by `CrystalMap.get_map_data`
//  (`CrystalMap.mapData`), a point not in data or not indexed black in an IPF map as orix leaves it (its colour array starts at
//  zero). The layers, the rasters, the hover's words, the colour key's labels and the unit cell stay with the app's room
//  (DiffractoKit#180). GPL-3.0 (orix is GPL-3.0-or-later), see NOTICE.
//
//  Parity: `tools/parity/sem` bundle `m3/maps-room-sdss` (the IPF-X, IPF-Y, IPF-Z and phase map images of orix-data's two-phase
//  steel, atol 1e-6 for the IPF colours, exact for the phase colours), docs/parity.md.
//

import Foundation

/// Per-point colours of a crystal map and the map's grid, as orix's `CrystalMap.plot` lays them out.
public nonisolated enum OrientationMaps {
    /// The map's picture size: rows then columns of the points in data (orix's `shape`); a 1D map is one row.
    public static func gridShape(_ map: CrystalMap) -> (rows: Int, columns: Int) {
        let shape = map.shape
        switch shape.count {
        case 2: return (shape[0], shape[1])
        case 1: return (1, shape[0])
        default: return (0, 0)
        }
    }

    /// orix's IPF colours of every point in data, in data order: each phase with `IPFColorKeyTSL(point_group, direction)`, a point
    /// not indexed or of a phase without a point group black (0, 0, 0), as orix's plot leaves it.
    public static func ipfColors(_ map: CrystalMap, direction: SIMD3<Double>) -> [SIMD3<Double>] {
        let ids = map.phaseIDs
        let rotations = map.rotationsInData
        var rgb = [SIMD3<Double>](repeating: .zero, count: ids.count)
        for (id, phase) in map.phases.entries where id != PhaseList.notIndexedID {
            guard let group = phase.pointGroup else { continue }
            let rows = ids.indices.filter { ids[$0] == id }
            guard !rows.isEmpty else { continue }
            let colours = IPFColorKeyTSL(group, direction: direction).color(of: rows.map { rotations[$0] })
            for (r, c) in zip(rows, colours) { rgb[r] = c }
        }
        return rgb
    }

    /// Each point in data in its phase's colour (orix's `Phase.color_rgb`; the not-indexed phase white), black for a colour
    /// matplotlib's Tableau list does not name.
    public static func phaseColors(_ map: CrystalMap) -> [SIMD3<Double>] {
        var byID: [Int: SIMD3<Double>] = [:]
        for (id, phase) in map.phases.entries { byID[id] = phase.colorRGB ?? .zero }
        return map.phaseIDs.map { byID[$0] ?? .zero }
    }

    /// Per-point colours (in data order) laid into the map's rows and columns (`CrystalMap.mapData` per channel), black where no
    /// point is in data.
    public static func image(_ colours: [SIMD3<Double>], of map: CrystalMap) -> [SIMD3<Double>] {
        let r = map.mapData(colours.map(\.x), fill: 0), g = map.mapData(colours.map(\.y), fill: 0), b = map.mapData(colours.map(\.z), fill: 0)
        return r.indices.map { SIMD3(r[$0], g[$0], b[$0]) }
    }

    /// The index into the map's per-point arrays (`euler`, `allPhaseIDs`, `properties`) of every grid cell, row-major, −1 where no
    /// point is in data. A map whose points stop short of its grid (a scan stopped mid-row) has no point in the cells past its last
    /// one: `mapData` leaves them NaN, and they read −1 too.
    public static func pointIndexGrid(_ map: CrystalMap) -> [Int] {
        map.mapData(map.ids.map(Double.init), fill: -1).map { $0.isFinite ? Int($0) : -1 }
    }
}
