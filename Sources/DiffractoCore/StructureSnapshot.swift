//
//  StructureSnapshot.swift
//  Role: a structure as a value a document or a library can store (Codable, with a schema version) and rebuild
//        exactly, and a content fingerprint that says whether two structures hold the same cell, sites and symmetry
//        reading. No app type: the library id, the provenance line and the revision stay in the app.
//

import CryptoKit
import Foundation
import simd

extension CrystalStructure {

    /// A structure as stored data: the cell, the sites, the name, the stated IT number and what its symmetry was read
    /// from, with the fingerprint of the structure it was made from. `structure()` rebuilds the structure: a snapshot of
    /// a structure from `detectingSymmetry` (a CIF, a preset, a Materials Project cell) runs the same detection with the
    /// same settings, so the rebuilt structure equals the original as long as the detection of this kit version reads
    /// the same; `structure().fingerprint == fingerprint` says that it does. A structure with given operations, and the
    /// conventional-setting structure of `ConventionalSetting`, are rebuilt from their stored operations.
    ///
    /// Every number is stored as it is in the structure, except that a negative zero is stored as zero (JSON has no
    /// negative zero that survives every reader). JSON cannot hold NaN or infinity: `JSONEncoder` refuses a snapshot
    /// with a non-finite value (a cell or site value, which `validationIssues` reports as not usable, or a position
    /// tolerance or written half-step of the detection that is not finite).
    public struct Snapshot: Sendable, Equatable, Codable {
        /// The schema this kit writes. A snapshot of a newer schema is refused on decoding with a sentence.
        public static let currentSchemaVersion = 1

        public let schemaVersion: Int
        public let name: String?
        public let cell: Cell
        public let sites: [Site]
        /// The IT number the source stated (`CrystalStructure.spaceGroupNumber`).
        public let spaceGroupNumber: Int?
        public let symmetry: Symmetry
        /// `CrystalStructure.fingerprint` of the structure the snapshot was made from.
        public let fingerprint: String

        /// The six cell parameters (Å and degrees), as `UnitCell` holds them.
        public struct Cell: Sendable, Equatable, Codable {
            public let a: Double
            public let b: Double
            public let c: Double
            public let alphaDeg: Double
            public let betaDeg: Double
            public let gammaDeg: Double
        }

        /// One `AtomSite`.
        public struct Site: Sendable, Equatable, Codable {
            public let atomicNumber: Int
            public let fractional: SIMD3<Double>
            public let occupancy: Double
            public let label: String?
            /// U_iso, Å² (`AtomSite.isotropicDisplacement`).
            public let isotropicDisplacement: Double?
        }

        /// One symmetry operation: the rotation's rows and the translation.
        public struct Operation: Sendable, Equatable, Codable {
            /// Row by row, nine integers.
            public let rotation: [Int]
            public let translation: SIMD3<Double>

            private enum CodingKeys: String, CodingKey { case rotation, translation }

            public init(from decoder: any Decoder) throws {
                let container = try decoder.container(keyedBy: CodingKeys.self)
                rotation = try container.decode([Int].self, forKey: .rotation)
                translation = try container.decode(SIMD3<Double>.self, forKey: .translation)
                guard rotation.count == 9 else {
                    throw DecodingError.dataCorruptedError(
                        forKey: .rotation, in: container,
                        debugDescription: "A symmetry operation's rotation holds nine integers; this one holds \(rotation.count).")
                }
            }
        }

        /// What the symmetry of the structure is read from.
        public enum Symmetry: Sendable, Equatable, Codable {
            /// `init(cell:sites:operations:)`: the operations as given.
            case given(operations: [Operation])
            /// `detectingSymmetry`: the file's operation loop (nil when there was none) and the detection's settings.
            case detected(fileOperations: [Operation]?, settings: Settings)
            /// `ConventionalSetting.structure`: the operations detected in the conventional cell and their uncertainty.
            case conventional(operations: [Operation], translationUncertainty: Double)

            private enum CodingKeys: String, CodingKey {
                case source, operations, fileOperations, settings, translationUncertainty
            }

            public init(from decoder: any Decoder) throws {
                let container = try decoder.container(keyedBy: CodingKeys.self)
                let source = try container.decode(String.self, forKey: .source)
                switch source {
                case "given":
                    self = .given(operations: try container.decode([Operation].self, forKey: .operations))
                case "detected":
                    self = .detected(fileOperations: try container.decodeIfPresent([Operation].self, forKey: .fileOperations),
                                     settings: try container.decode(Settings.self, forKey: .settings))
                case "conventional":
                    self = .conventional(
                        operations: try container.decode([Operation].self, forKey: .operations),
                        translationUncertainty: try container.decode(Double.self, forKey: .translationUncertainty))
                default:
                    throw DecodingError.dataCorruptedError(
                        forKey: .source, in: container,
                        debugDescription: "The symmetry source \"\(source)\" is not one this version knows (given, detected, conventional).")
                }
            }

            public func encode(to encoder: any Encoder) throws {
                var container = encoder.container(keyedBy: CodingKeys.self)
                switch self {
                case .given(let operations):
                    try container.encode("given", forKey: .source)
                    try container.encode(operations, forKey: .operations)
                case .detected(let fileOperations, let settings):
                    try container.encode("detected", forKey: .source)
                    try container.encodeIfPresent(fileOperations, forKey: .fileOperations)
                    try container.encode(settings, forKey: .settings)
                case .conventional(let operations, let translationUncertainty):
                    try container.encode("conventional", forKey: .source)
                    try container.encode(operations, forKey: .operations)
                    try container.encode(translationUncertainty, forKey: .translationUncertainty)
                }
            }
        }

        /// `CrystalStructure.DetectionSettings` as stored data.
        public struct Settings: Sendable, Equatable, Codable {
            public let coarsestHalfStep: Double
            public let metricLengthTolerance: Double
            public let metricAngleToleranceDeg: Double
            /// "fractional" or "cartesian".
            public let toleranceKind: String
            public let tolerance: Double
            public let conventional: Bool
            public let occupancyTolerance: Double
            public let translationGrid: Double
            public let gridSnapFactor: Double
            public let precisionLimitFactor: Double

            private enum CodingKeys: String, CodingKey {
                case coarsestHalfStep, metricLengthTolerance, metricAngleToleranceDeg, toleranceKind, tolerance, conventional
                case occupancyTolerance, translationGrid, gridSnapFactor, precisionLimitFactor
            }

            public init(from decoder: any Decoder) throws {
                let container = try decoder.container(keyedBy: CodingKeys.self)
                coarsestHalfStep = try container.decode(Double.self, forKey: .coarsestHalfStep)
                metricLengthTolerance = try container.decode(Double.self, forKey: .metricLengthTolerance)
                metricAngleToleranceDeg = try container.decode(Double.self, forKey: .metricAngleToleranceDeg)
                toleranceKind = try container.decode(String.self, forKey: .toleranceKind)
                guard toleranceKind == "fractional" || toleranceKind == "cartesian" else {
                    throw DecodingError.dataCorruptedError(
                        forKey: .toleranceKind, in: container,
                        debugDescription: "The position tolerance kind \"\(toleranceKind)\" is not one this version knows (fractional, cartesian).")
                }
                tolerance = try container.decode(Double.self, forKey: .tolerance)
                conventional = try container.decode(Bool.self, forKey: .conventional)
                occupancyTolerance = try container.decode(Double.self, forKey: .occupancyTolerance)
                translationGrid = try container.decode(Double.self, forKey: .translationGrid)
                gridSnapFactor = try container.decode(Double.self, forKey: .gridSnapFactor)
                precisionLimitFactor = try container.decode(Double.self, forKey: .precisionLimitFactor)
            }
        }

        private enum CodingKeys: String, CodingKey {
            case schemaVersion, name, cell, sites, spaceGroupNumber, symmetry, fingerprint
        }

        /// The snapshot of `structure`.
        public init(_ structure: CrystalStructure) {
            schemaVersion = Self.currentSchemaVersion
            name = structure.name
            let cell = structure.cell
            self.cell = Cell(a: positiveZero(cell.a), b: positiveZero(cell.b), c: positiveZero(cell.c),
                             alphaDeg: positiveZero(cell.alphaDeg), betaDeg: positiveZero(cell.betaDeg),
                             gammaDeg: positiveZero(cell.gammaDeg))
            sites = structure.sites.map {
                Site(atomicNumber: $0.z, fractional: positiveZero($0.fractional), occupancy: positiveZero($0.occupancy),
                     label: $0.label, isotropicDisplacement: $0.isotropicDisplacement.map { positiveZero($0) })
            }
            spaceGroupNumber = structure.spaceGroupNumber
            if structure.setting == .conventional {
                symmetry = .conventional(operations: structure.operations.map(Operation.init),
                                         translationUncertainty: positiveZero(structure.translationUncertainty))
            } else if let settings = structure.detectionSettings {
                symmetry = .detected(fileOperations: structure.fileOperations?.map(Operation.init),
                                     settings: Settings(settings))
            } else {
                symmetry = .given(operations: structure.operations.map(Operation.init))
            }
            fingerprint = structure.fingerprint
        }

        public init(from decoder: any Decoder) throws {
            let container = try decoder.container(keyedBy: CodingKeys.self)
            let version = try container.decode(Int.self, forKey: .schemaVersion)
            guard version >= 1, version <= Self.currentSchemaVersion else {
                throw DecodingError.dataCorruptedError(
                    forKey: .schemaVersion, in: container,
                    debugDescription: "This structure was saved with schema \(version); this version of DiffractoKit reads "
                        + "schemas 1 to \(Self.currentSchemaVersion). A newer version of the app wrote it.")
            }
            schemaVersion = version
            name = try container.decodeIfPresent(String.self, forKey: .name)
            cell = try container.decode(Cell.self, forKey: .cell)
            sites = try container.decode([Site].self, forKey: .sites)
            spaceGroupNumber = try container.decodeIfPresent(Int.self, forKey: .spaceGroupNumber)
            symmetry = try container.decode(Symmetry.self, forKey: .symmetry)
            fingerprint = try container.decode(String.self, forKey: .fingerprint)
        }

        public func encode(to encoder: any Encoder) throws {
            var container = encoder.container(keyedBy: CodingKeys.self)
            try container.encode(schemaVersion, forKey: .schemaVersion)
            try container.encodeIfPresent(name, forKey: .name)
            try container.encode(cell, forKey: .cell)
            try container.encode(sites, forKey: .sites)
            try container.encodeIfPresent(spaceGroupNumber, forKey: .spaceGroupNumber)
            try container.encode(symmetry, forKey: .symmetry)
            try container.encode(fingerprint, forKey: .fingerprint)
        }

        /// The structure again: the same detection with the same settings, or the stored operations (see the type's
        /// note). Compare its `fingerprint` with `fingerprint` to learn whether this kit version reads the symmetry
        /// as the one that made the snapshot.
        public func structure() -> CrystalStructure {
            let unitCell = UnitCell(a: cell.a, b: cell.b, c: cell.c,
                                    alphaDeg: cell.alphaDeg, betaDeg: cell.betaDeg, gammaDeg: cell.gammaDeg)
            let atomSites = sites.map {
                AtomSite(z: $0.atomicNumber, fractional: $0.fractional, occupancy: $0.occupancy, label: $0.label,
                         isotropicDisplacement: $0.isotropicDisplacement)
            }
            switch symmetry {
            case .given(let operations):
                return CrystalStructure(cell: unitCell, sites: atomSites, operations: operations.map(\.operation),
                                        spaceGroupNumber: spaceGroupNumber, name: name)
            case .detected(let fileOperations, let settings):
                let detection = settings.detectionSettings
                return CrystalStructure.detectingSymmetry(
                    cell: unitCell, sites: atomSites, fileOperations: fileOperations?.map(\.operation),
                    spaceGroupNumber: spaceGroupNumber, name: name, coarsestHalfStep: detection.coarsestHalfStep,
                    metricLengthTolerance: detection.metricLengthTolerance,
                    metricAngleToleranceDeg: detection.metricAngleToleranceDeg, tolerance: detection.tolerance,
                    conventional: detection.conventional, policy: detection.policy)
            case .conventional(let operations, let translationUncertainty):
                let ops = operations.map(\.operation)
                let group = PointSymmetry.laueGroup(of: ops)
                return CrystalStructure(conventionalCell: unitCell, sites: atomSites, operations: ops, laueGroup: group,
                                        laueClass: group.laueClass, translationUncertainty: translationUncertainty)
            }
        }
    }

    /// SHA-256 (hex) of the structure's content: the six cell parameters, every site's atomic number, fractional
    /// position, occupancy and displacement, the setting, and the symmetry reading (the operations, translations mod 1
    /// on a 1e-6 grid, so that a last-bit difference in a refined translation almost never reads as a different reading). Equal for
    /// equal content whatever the order of the sites or of the operations; any change of a cell parameter or of a site's
    /// element, position, occupancy or displacement changes it. Not in it: the name, the site labels, the stated IT number, the file's own operation loop and the settings
    /// the detection was given (two reads that find the same operations have the same fingerprint). A negative zero
    /// counts as zero.
    public var fingerprint: String {
        var bytes = Data("DiffractoKit structure fingerprint 1".utf8)
        bytes.append(0)
        for value in [cell.a, cell.b, cell.c, cell.alphaDeg, cell.betaDeg, cell.gammaDeg] { appendDouble(value, to: &bytes) }

        let siteRecords: [Data] = sites.map { site in
            var record = Data()
            appendInt(site.z, to: &record)
            for value in [site.fractional.x, site.fractional.y, site.fractional.z, site.occupancy] {
                appendDouble(value, to: &record)
            }
            if let u = site.isotropicDisplacement {
                record.append(1)
                appendDouble(u, to: &record)
            } else {
                record.append(0)
            }
            return record
        }
        appendInt(siteRecords.count, to: &bytes)
        for record in siteRecords.sorted(by: { $0.lexicographicallyPrecedes($1) }) { bytes.append(record) }

        bytes.append(setting == .conventional ? 1 : 0)
        let operationRecords = Set(operations.map { op -> Data in
            var record = Data()
            for row in op.rotation.rows { for entry in [row.x, row.y, row.z] { appendInt(entry, to: &record) } }
            for t in [op.translation.x, op.translation.y, op.translation.z] {
                appendInt(gridTranslation(t), to: &record)
            }
            return record
        })
        appendInt(operationRecords.count, to: &bytes)
        for record in operationRecords.sorted(by: { $0.lexicographicallyPrecedes($1) }) { bytes.append(record) }

        return SHA256.hash(data: bytes).map { String(format: "%02x", $0) }.joined()
    }
}

extension CrystalStructure.Snapshot.Operation {
    init(_ op: SymmetryOperation) {
        rotation = op.rotation.rows.flatMap { [$0.x, $0.y, $0.z] }
        translation = positiveZero(op.translation)
    }

    /// Nine entries: `init` and decoding make sure of it.
    var operation: SymmetryOperation {
        SymmetryOperation(
            rotation: IntegerMatrix(SIMD3(rotation[0], rotation[1], rotation[2]), SIMD3(rotation[3], rotation[4], rotation[5]),
                                    SIMD3(rotation[6], rotation[7], rotation[8])),
            translation: translation)
    }
}

extension CrystalStructure.Snapshot.Settings {
    init(_ s: CrystalStructure.DetectionSettings) {
        coarsestHalfStep = positiveZero(s.coarsestHalfStep)
        metricLengthTolerance = positiveZero(s.metricLengthTolerance)
        metricAngleToleranceDeg = positiveZero(s.metricAngleToleranceDeg)
        switch s.tolerance {
        case .fractional(let t):
            toleranceKind = "fractional"
            tolerance = positiveZero(t)
        case .cartesian(let t):
            toleranceKind = "cartesian"
            tolerance = positiveZero(t)
        }
        conventional = s.conventional
        occupancyTolerance = positiveZero(s.policy.occupancyTolerance)
        translationGrid = positiveZero(s.policy.translationGrid)
        gridSnapFactor = positiveZero(s.policy.gridSnapFactor)
        precisionLimitFactor = positiveZero(s.policy.precisionLimitFactor)
    }

    /// The kind is "fractional" or "cartesian": `init(from:)` and `init(_:)` make sure of it.
    var detectionSettings: CrystalStructure.DetectionSettings {
        CrystalStructure.DetectionSettings(
            coarsestHalfStep: coarsestHalfStep, metricLengthTolerance: metricLengthTolerance,
            metricAngleToleranceDeg: metricAngleToleranceDeg,
            tolerance: toleranceKind == "cartesian" ? .cartesian(tolerance) : .fractional(tolerance),
            conventional: conventional,
            policy: SymmetryDetection.Policy(occupancyTolerance: occupancyTolerance, translationGrid: translationGrid,
                                             gridSnapFactor: gridSnapFactor, precisionLimitFactor: precisionLimitFactor))
    }
}

/// Zero for a negative zero, the value otherwise.
private func positiveZero(_ value: Double) -> Double { value == 0 ? 0 : value }

private func positiveZero(_ v: SIMD3<Double>) -> SIMD3<Double> {
    SIMD3(positiveZero(v.x), positiveZero(v.y), positiveZero(v.z))
}

/// The bits of a value, little-endian; a negative zero as zero, every NaN as one NaN.
private func appendDouble(_ value: Double, to data: inout Data) {
    let canonical = value.isNaN ? Double.nan : positiveZero(value)
    withUnsafeBytes(of: canonical.bitPattern.littleEndian) { data.append(contentsOf: $0) }
}

private func appendInt(_ value: Int, to data: inout Data) {
    withUnsafeBytes(of: Int64(value).littleEndian) { data.append(contentsOf: $0) }
}

/// A translation mod 1 on the 1e-6 grid, 0 ..< 1_000_000; a non-finite translation as -1.
private func gridTranslation(_ t: Double) -> Int {
    guard t.isFinite else { return -1 }
    let wrapped = t - t.rounded(.down)
    let step = (wrapped * 1e6).rounded()
    guard step.isFinite else { return -1 }
    return Int(step) % 1_000_000
}
