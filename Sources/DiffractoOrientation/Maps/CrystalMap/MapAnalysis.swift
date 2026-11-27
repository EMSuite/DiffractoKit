//
//  MapAnalysis.swift
//  Moved from Electron Microscopy Studio's Sources/TEMSEM/CrystalMap (at 7b09bf5) by DiffractoKit#180; the bundles its moved
//  tests read are in Tests/DiffractoOrientationTests/Fixtures/orix-sem, the others stay with the app (docs/parity.md).
//  What the app's Maps room computes from a crystal map beyond its colours (the app's brief P5-04): the grains, their
//  boundaries, KAM and GOS (DefDAP 1.3.1's functions, `GrainAnalysis`), the GND density (`NyeTensor`, Pantleon 2008 as DefDAP
//  computes it), the grain table (size, equivalent diameter, mean orientation, GOS, neighbours) with the ASTM E112 grain size
//  number from the mean lineal intercept, and each grain's Schmid factors and slip traces (`SlipTraceAnalysis`). Every number
//  is one of those ported functions on the map's own arrays, per grid cell, row-major, NaN where no value is. The room's layers,
//  rasters, grain colours, boundary and outline segments, histograms and the texture pane (`MapTexture`) stay with the app
//  (DiffractoKit#180).
//
//  The frame: the map's Euler angles as the reader gave them (orix's convention, lab to crystal), so the Schmid load and the
//  slip traces are in the frame the IPF maps colour; DefDAP's EDAX frame change is not applied (D-SEM-101). The quaternion
//  DefDAP builds from the Bunge angles is orix's `eu2qu` term for term, so a grain's mean orientation reads back through
//  `RotationConversions.qu2eu`.
//
//  Parity: the grain map, KAM and GOS of orix-data's two-phase steel against `m5/grains-sdss` (the app's `MapsAnalysisTests`, which
//  read the map with the app's `.ang` reader), docs/parity.md.
//

import Foundation

/// The parameters of a map's grain analysis (the app's Maps room holds them in its inspector's Grains and Misorientation sections).
public nonisolated struct MapAnalysisSettings: Sendable, Equatable, Hashable {
    /// The misorientation (°) above which two neighbours are on different grains.
    public var boundaryDegrees: Double = 5
    /// Grains of fewer points are left out (their points read −2, DefDAP's mark).
    public var minimumSize: Int = 10
    /// KAM drops a neighbour more misoriented than this (°).
    public var kamMaxDegrees: Double = 5

    public init(boundaryDegrees: Double = 5, minimumSize: Int = 10, kamMaxDegrees: Double = 5) {
        self.boundaryDegrees = boundaryDegrees
        self.minimumSize = minimumSize
        self.kamMaxDegrees = kamMaxDegrees
    }
}

/// One row of the grain table.
public nonisolated struct GrainRow: Sendable, Equatable, Identifiable {
    /// The grain's number on the grain map (1, 2, ...).
    public let id: Int
    /// The crystal map's phase id, nil for a grain DefDAP could not give one phase.
    public let phaseID: Int?
    public let phaseName: String
    public let points: Int
    /// Points × step², in the scan unit squared.
    public let area: Double
    /// The diameter of the circle of the grain's area, in the scan unit.
    public let diameter: Double
    /// The mean orientation (DefDAP's `calc_average_ori`) as Bunge angles in degrees.
    public let euler: SIMD3<Double>
    /// The mean orientation as a quaternion (lab to crystal).
    public let orientation: SIMD4<Double>
    /// The grain orientation spread, degrees.
    public let gosDegrees: Double
    /// Grains sharing a boundary with this one.
    public let neighbours: Int
    /// The centroid: column, row.
    public let centroid: SIMD2<Double>
}

/// A grain's slip systems under one load: per plane group its trace and the largest Schmid factor of its directions.
public nonisolated struct GrainSlip: Sendable, Equatable {
    public struct Plane: Sendable, Equatable {
        /// "(111)".
        public let label: String
        /// Radians from the map's +y axis, counter-clockwise; NaN for a plane parallel to the map.
        public let traceAngle: Double
        /// Radians between the plane and the map plane.
        public let inclination: Double
        /// The largest Schmid factor of the plane's systems.
        public let schmid: Double
        /// DefDAP's trace colour of the group.
        public let colour: String
    }

    public let grain: Int
    public let planes: [Plane]
    public var maximumSchmid: Double { planes.map(\.schmid).max() ?? .nan }

    /// Each grain's `maximumSchmid` by grain number, the first entry of a grain listed twice: a table looks a grain up in it per cell
    /// (#271), where a search of the list per cell is O(rows × grains).
    public static func maximumSchmidByGrain(_ slip: [GrainSlip]) -> [Int: Double] {
        var byGrain: [Int: Double] = [:]
        for entry in slip where byGrain[entry.grain] == nil { byGrain[entry.grain] = entry.maximumSchmid }
        return byGrain
    }
}

public nonisolated struct MapAnalysis: Sendable {
    public enum Refusal: Error, Equatable, Sendable, LocalizedError {
        /// The points in data do not fill the map's rectangle in row order (a sliced map).
        case notAGrid
        case grains(String)

        public var errorDescription: String? {
            switch self {
            case .notAGrid: "Grains need the whole map: this one is a slice of its scan."
            case .grains(let text): text
            }
        }
    }

    public let settings: MapAnalysisSettings
    public let rows: Int
    public let columns: Int
    public let grainMap: GrainMap
    public let grains: GrainAnalysis
    /// Per grid cell, the GND density in m⁻² (NaN where not indexed); nil when the first phase has no crystal structure DefDAP knows.
    public let gnd: [Double]?
    /// The grain table, in grain number order.
    public let table: [GrainRow]
    /// The crystal map's phase id of each of DefDAP's phase numbers (index = number − 1).
    public let phaseIDs: [Int]
    /// The mean lineal intercept in the scan unit; nil without a boundary.
    public let meanLinealIntercept: Double?
    /// The scan unit of the lengths.
    public let scanUnit: String
    /// Why `gnd` is nil although the phase has a structure DefDAP knows (a step in no length unit); nil otherwise.
    public let gndRefusal: String?

    /// The analysis of `map` with `settings` (a whole map: the points in data fill rows × columns in row order). `isCancelled` is asked
    /// through the grain analysis (`GrainAnalysis.run`) and before the GND and the table; it throws `CancellationError` once it says so
    /// (the Maps room's Stop).
    public static func run(_ map: CrystalMap, settings: MapAnalysisSettings = MapAnalysisSettings(),
                            isCancelled: () -> Bool = { false }) throws -> MapAnalysis {
        let (rows, columns) = OrientationMaps.gridShape(map)
        guard rows * columns == map.size, OrientationMaps.pointIndexGrid(map) == map.ids else { throw Refusal.notAGrid }
        let grainMap: GrainMap
        do { grainMap = try GrainMap(map: map) } catch { throw Refusal.grains(error.localizedDescription) }
        let analysis: GrainAnalysis
        do {
            analysis = try GrainAnalysis.run(grainMap, boundaryDegrees: settings.boundaryDegrees, minimumSize: settings.minimumSize,
                                             kam: .symmetric(thresholdDegrees: settings.kamMaxDegrees), isCancelled: isCancelled)
        } catch let cancelled as CancellationError {
            throw cancelled
        } catch {
            throw Refusal.grains(error.localizedDescription)
        }
        if isCancelled() { throw CancellationError() }
        let indexed = map.phases.entries.filter { $0.id != PhaseList.notIndexedID }
        let phaseIDs = indexed.map(\.id)

        // GND: DefDAP takes the primary phase's crystal structure for the whole map; the step in micrometres.
        var gnd: [Double]?
        var gndRefusal: String?
        if let first = indexed.first, let group = first.phase.pointGroup, let structure = DefDAPCrystalStructure(pointGroup: group) {
            // A step in no length unit refuses the GND rather than guessing µm (macEM#611).
            if let micrometres = Self.micrometres(grainMap.step, unit: map.scanUnit) {
                let result = NyeTensor.compute(quaternions: grainMap.quaternions, rows: rows, columns: columns, stepMicrons: micrometres,
                                               structure: structure)
                gnd = zip(result.gnd, grainMap.phase).map { $1 == 0 ? .nan : $0 }
            } else {
                gndRefusal = Self.gndNeedsALength(unit: map.scanUnit)
            }
        }
        if isCancelled() { throw CancellationError() }

        // The table.
        let segmentation = analysis.segmentation
        let gos = analysis.statistics.gos
        var neighbourCount = [Int](repeating: 0, count: segmentation.grains.count)
        for pair in analysis.neighbours {
            if pair.a < neighbourCount.count { neighbourCount[pair.a] += 1 }
            if pair.b < neighbourCount.count { neighbourCount[pair.b] += 1 }
        }
        let areas = segmentation.areas, diameters = segmentation.equivalentCircleDiameters
        var table: [GrainRow] = []
        table.reserveCapacity(segmentation.grains.count)
        for (k, grain) in segmentation.grains.enumerated() {
            let phaseID = grain.phaseIndex.map { phaseIDs[$0] }
            let name = phaseID.flatMap { id in map.phases[id: id]?.name } ?? "mixed"
            var sum = SIMD2<Double>(0, 0)
            for p in grain.points { sum += SIMD2(Double(p % columns), Double(p / columns)) }
            let q = analysis.statistics.averageOrientations[k]
            table.append(GrainRow(id: grain.id + 1, phaseID: phaseID, phaseName: name, points: grain.size, area: areas[k],
                                  diameter: diameters[k], euler: RotationConversions.qu2eu(q) * (180 / Double.pi), orientation: q,
                                  gosDegrees: gos[k] * 180 / Double.pi, neighbours: neighbourCount[k],
                                  centroid: sum / Double(max(grain.size, 1))))
        }
        return MapAnalysis(settings: settings, rows: rows, columns: columns, grainMap: grainMap, grains: analysis, gnd: gnd, table: table,
                           phaseIDs: phaseIDs,
                           meanLinealIntercept: GrainSegmentation.meanLinealIntercept(boundaries: analysis.boundaries, step: grainMap.step, stepY: grainMap.stepY),
                           scanUnit: map.scanUnit, gndRefusal: gndRefusal)
    }

    /// The same off the caller's actor; a cancelled task stops at its next row or grain and throws `CancellationError`.
    @concurrent public static func runAsync(_ map: CrystalMap, settings: MapAnalysisSettings) async throws -> MapAnalysis {
        try run(map, settings: settings, isCancelled: { Task.isCancelled })
    }

    /// A length in `unit` in micrometres (the readers write "um"; nm and mm convert); nil for a unit that is no length ("px", or one
    /// the app does not know), which before macEM#611 was taken as micrometres.
    public static func micrometres(_ length: Double, unit: String) -> Double? {
        switch unit.lowercased() {
        case "nm": length / 1000
        case "mm": length * 1000
        case "um", "µm", "μm", "micron", "microns": length
        default: nil
        }
    }

    /// Why the map has no GND density although its phase has a structure: its step is in no length unit.
    static func gndNeedsALength(unit: String) -> String {
        "No GND density: the map's step is in \(unit.isEmpty ? "no unit" : "\"\(unit)\"") rather than a length, so a density per m² cannot be given. The grains, KAM and GOS do not need one."
    }

    // MARK: Numbers

    /// The grains' mean equivalent diameter by number and weighted by area, in the scan unit.
    public var meanDiameter: (byNumber: Double, byArea: Double)? {
        guard !table.isEmpty else { return nil }
        let n = Double(table.count)
        let byNumber = table.reduce(0) { $0 + $1.diameter } / n
        let totalArea = table.reduce(0) { $0 + $1.area }
        let byArea = totalArea > 0 ? table.reduce(0) { $0 + $1.diameter * $1.area } / totalArea : .nan
        return (byNumber, byArea)
    }

    /// ASTM E112's grain size number from the mean lineal intercept; nil without a boundary or with a scan unit that is no length.
    public var astmGrainSizeNumber: Double? {
        guard let l = meanLinealIntercept, let um = Self.micrometres(l, unit: scanUnit) else { return nil }
        return GrainSegmentation.astmGrainSizeNumber(meanLinealInterceptMillimetres: um / 1000)
    }

    /// The grain number (1, 2, ...) at a grid cell; nil for a point of no grain (too small, or not indexed).
    public func grain(atRow row: Int, column: Int) -> Int? {
        guard row >= 0, row < rows, column >= 0, column < columns else { return nil }
        let g = grains.segmentation.grainMap[row * columns + column]
        return g > 0 ? g : nil
    }

    // MARK: Slip

    /// DefDAP's slip-system file for a phase: hexagonal → `hexagonal_withca`; cubic with space group 225 → `cubic_fcc`, 229 →
    /// `cubic_bcc`; a cubic phase whose file names no space group (EDAX's .ang) by its name: "austenite" or "fcc" → `cubic_fcc`,
    /// "ferrite", "martensite" or "bcc" → `cubic_bcc` (D-SEM-101); nil otherwise, and the inspector's Slip section picks.
    public static func defaultSlipFile(for phase: CrystalPhase) -> SlipSystemFile? {
        guard let group = phase.pointGroup, let structure = DefDAPCrystalStructure(pointGroup: group) else { return nil }
        switch structure {
        case .hexagonal: return .hexagonalWithCA
        case .cubic:
            switch phase.spaceGroup {
            case 225: return .cubicFCC
            case 229: return .cubicBCC
            case .some: return nil
            case .none:
                let name = phase.name.lowercased()
                if ["austenite", "fcc"].contains(where: { name.contains($0) }) { return .cubicFCC }
                if ["ferrite", "martensite", "bcc"].contains(where: { name.contains($0) }) { return .cubicBCC }
                return nil
            }
        }
    }

    /// The slip systems of `file` for a phase (DefDAP's, without systems whose direction leaves their plane); nil when the phase's
    /// crystal structure does not fit the file.
    public static func slipSystems(_ file: SlipSystemFile, phase: CrystalPhase) -> SlipSystemSet? {
        guard let group = phase.pointGroup, let structure = DefDAPCrystalStructure(pointGroup: group) else { return nil }
        let crystal: SlipCrystal
        switch structure {
        case .cubic:
            guard !file.isHexagonal else { return nil }
            crystal = .cubic
        case .hexagonal:
            guard file.isHexagonal else { return nil }
            let l = phase.lattice ?? []
            crystal = .hexagonal(cOverA: l.count >= 3 && l[0] > 0 ? l[2] / l[0] : 1.633)
        }
        return (try? SlipSystemSet.load(file, crystal: crystal))?.consistent
    }

    /// Every grain's slip under `load` (a sample direction), for the phases `files` gives systems.
    public func slip(load: SIMD3<Double>, files: [Int: SlipSystemFile], phases: PhaseList) -> [GrainSlip] {
        var systems: [Int: SlipSystemSet] = [:]
        for (id, file) in files {
            if let phase = phases[id: id], let set = Self.slipSystems(file, phase: phase) { systems[id] = set }
        }
        return table.compactMap { row in
            guard let id = row.phaseID, let set = systems[id] else { return nil }
            let schmid = SlipTraceAnalysis.schmidFactors(orientation: row.orientation, load: load, systems: set)
            let traces = SlipTraceAnalysis.traces(orientation: row.orientation, systems: set)
            let planes = set.groups.indices.map { k in
                GrainSlip.Plane(label: set.groups[k].first?.planeLabel ?? "", traceAngle: traces[k].angle, inclination: traces[k].inclination,
                                schmid: schmid[k].max() ?? .nan,
                                // DefDAP's `add_traces`: a trace past the colour list takes the list's last colour.
                                colour: k < set.traceColours.count ? set.traceColours[k] : (set.traceColours.last ?? "black"))
            }
            return GrainSlip(grain: row.id, planes: planes)
        }
    }
}
