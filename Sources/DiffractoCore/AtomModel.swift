//
//  AtomModel.swift
//  Role: an atom model built elsewhere (ASE, atomsk, VESTA, a molecular-dynamics run) read from extended XYZ or POSCAR:
//        the atoms as atomic numbers and Cartesian positions in Å, the cell's three rows, periodicity, the number of
//        frames. Plain data with no symmetry: a crystal with symmetry is `CrystalStructure` (CIF). The two readers are
//        `AtomModel+XYZ.swift` and `AtomModel+POSCAR.swift`; the byte scanner they share is `AtomModel+Text.swift`.
//
//  Reference: `ase.io.read` of ASE 3.29.0 (`ase/io/extxyz.py` `read_xyz`, `ase/io/vasp.py` `read_vasp`), pinned by
//  tools/pins/atom-models/pins.py (DiffractoKit #169, plan row MD1). Where the kit departs from ASE it refuses, by name
//  (`AtomModelError`), what ASE would read with help from another file or turn into a model no simulation can use: a
//  POSCAR without its element line, the dummy element X, a number written as nan or inf, a cell of zero volume.
//

import Foundation

/// Atoms read from an extended XYZ or POSCAR file: what `ase.io.read` of ASE 3.29.0 returns for it.
public struct AtomModel: Sendable, Equatable {

    /// The file formats the kit reads atom models from.
    public enum Format: String, Sendable, Equatable, CaseIterable {
        /// Extended XYZ, and plain XYZ as its special case (ASE's `extxyz`; a `.xyz` file).
        case extendedXYZ
        /// VASP 5 POSCAR or CONTCAR (ASE's `vasp`).
        case poscar
    }

    /// The format the model was read from.
    public let format: Format
    /// Atomic number of each atom, 1 ... 118, in file order.
    public let atomicNumbers: [Int]
    /// Cartesian position of each atom in Å, in file order. Positions are as the file states them (XYZ, Cartesian
    /// POSCAR at scale 1) or converted as ASE converts them (Direct coordinates times the cell; a scale factor).
    public let positions: [SIMD3<Double>]
    /// The cell's three vectors a, b, c in Å (ASE's `atoms.cell` rows). A vector the file does not give is zero: all
    /// three for a plain XYZ file, so `hasCell` is false.
    public let cell: [SIMD3<Double>]
    /// Periodicity along a, b and c (ASE's `atoms.pbc`): true for each vector of a POSCAR; for XYZ, the file's `pbc`,
    /// else true along all three when it gives a `Lattice`, else false.
    public let periodic: [Bool]
    /// How many frames the file holds. Of several (a trajectory), the last is the one read, as `ase.io.read` does.
    public let frameCount: Int
    /// The comment line: POSCAR's first line, the XYZ frame's second line, without surrounding whitespace.
    public let comment: String

    /// The number of atoms.
    public var count: Int { atomicNumbers.count }

    /// Whether the file gives at least one cell vector.
    public var hasCell: Bool { cell.contains { $0 != SIMD3<Double>(0, 0, 0) } }

    /// Atoms per element, by atomic number.
    public var countsByAtomicNumber: [Int: Int] {
        var counts: [Int: Int] = [:]
        for z in atomicNumbers { counts[z, default: 0] += 1 }
        return counts
    }

    package init(format: Format, atomicNumbers: [Int], positions: [SIMD3<Double>], cell: [SIMD3<Double>],
                 periodic: [Bool], frameCount: Int, comment: String) {
        self.format = format
        self.atomicNumbers = atomicNumbers
        self.positions = positions
        self.cell = cell
        self.periodic = periodic
        self.frameCount = frameCount
        self.comment = comment
    }
}

/// Why an atom-model file is refused. Every case names the line (counted from 1) or the text a person can find in
/// their file.
public enum AtomModelError: Error, Equatable, Sendable, LocalizedError {
    /// The file has more than `AtomModel.maximumBytes` bytes.
    case tooLarge(bytes: Int)
    /// The bytes hold a NUL byte, so the file is not text (a binary file, or UTF-16).
    case notText
    /// The file's name says neither XYZ nor POSCAR.
    case unknownFormat(fileName: String)
    /// The file holds no atoms: it is empty, or its count is 0.
    case noAtoms
    /// A line does not hold what the format puts there. `expected` says what, in a few words.
    case malformedLine(line: Int, expected: String)
    /// A number is not one (`1.0d0`, `?`), or is nan, inf or beyond the range of a double.
    case badNumber(text: String, line: Int)
    /// An element symbol (or atomic number) names no element 1 ... 118. The dummy `X` is refused here too.
    case unknownElement(symbol: String, line: Int)
    /// The file ends before the atoms it announces: `expected` atoms, `found` lines.
    case missingAtoms(expected: Int, found: Int)
    /// A POSCAR without the element line (VASP 4): ASE takes the elements from a POTCAR or OUTCAR beside it; the kit
    /// reads no other file. Add the element line above the counts (VASP 5).
    case noElementLine
    /// A POSCAR's counts line names more species than the element line.
    case countsWithoutElements(line: Int, elements: Int, counts: Int)
    /// A POSCAR's scale line holds neither one number nor three positive numbers.
    case badScale(line: Int)
    /// The cell has zero volume (a POSCAR), so Direct coordinates and the scale mean nothing.
    case degenerateCell
    /// An extended XYZ `Lattice` that is not nine numbers.
    case badLattice(line: Int)
    /// An extended XYZ `pbc` that is not one or three of T, F, True, False (or numbers).
    case badPeriodicity(line: Int)
    /// An extended XYZ `Properties` value the reader cannot follow: `reason` says why.
    case badProperties(line: Int, reason: String)
    /// The comment line of an extended XYZ frame cannot be read as key=value pairs (it starts with `=`).
    case badCommentLine(line: Int)

    /// The case's name without its values.
    public var name: String {
        switch self {
        case .tooLarge: return "tooLarge"
        case .notText: return "notText"
        case .unknownFormat: return "unknownFormat"
        case .noAtoms: return "noAtoms"
        case .malformedLine: return "malformedLine"
        case .badNumber: return "badNumber"
        case .unknownElement: return "unknownElement"
        case .missingAtoms: return "missingAtoms"
        case .noElementLine: return "noElementLine"
        case .countsWithoutElements: return "countsWithoutElements"
        case .badScale: return "badScale"
        case .degenerateCell: return "degenerateCell"
        case .badLattice: return "badLattice"
        case .badPeriodicity: return "badPeriodicity"
        case .badProperties: return "badProperties"
        case .badCommentLine: return "badCommentLine"
        }
    }

    /// One sentence for the person reading it.
    public var errorDescription: String? {
        switch self {
        case .tooLarge(let bytes):
            return "This file is \(bytes) bytes; an atom model larger than \(AtomModel.maximumBytes) bytes is not read."
        case .notText:
            return "This file is not text (it holds a NUL byte), so it cannot be an XYZ or POSCAR file."
        case .unknownFormat(let fileName):
            return "\"\(fileName)\" is neither an XYZ file (.xyz) nor a POSCAR (POSCAR, CONTCAR or .vasp)."
        case .noAtoms:
            return "The file holds no atoms."
        case .malformedLine(let line, let expected):
            return "Line \(line) should hold \(expected)."
        case .badNumber(let text, let line):
            return "Line \(line): \"\(text)\" is not a finite number."
        case .unknownElement(let symbol, let line):
            return "Line \(line): \"\(symbol)\" is not an element."
        case .missingAtoms(let expected, let found):
            return "The file announces \(expected) atoms but ends after \(found)."
        case .noElementLine:
            return "This POSCAR has no element line (VASP 4 format). Add the element symbols on a line above the counts."
        case .countsWithoutElements(let line, let elements, let counts):
            return "Line \(line) gives \(counts) counts for \(elements) elements."
        case .badScale(let line):
            return "Line \(line) should hold one scale factor, or three positive ones."
        case .degenerateCell:
            return "The cell has zero volume."
        case .badLattice(let line):
            return "Line \(line): Lattice should be nine numbers, the vectors a, b and c."
        case .badPeriodicity(let line):
            return "Line \(line): pbc should be one or three of T and F."
        case .badProperties(let line, let reason):
            return "Line \(line): the Properties value cannot be read (\(reason))."
        case .badCommentLine(let line):
            return "Line \(line) cannot be read as key=value pairs."
        }
    }
}

extension AtomModel {

    /// The most bytes an atom-model file may have: 512 MiB (about ten million atoms of XYZ). One byte more is refused.
    public static let maximumBytes = 512 * 1024 * 1024

    /// The format a file's name says, as ASE's `filetype` does for these two: `.xyz` and `.extxyz` are XYZ; `.vasp`,
    /// `.poscar` and a name that holds POSCAR or CONTCAR are POSCAR. Nil for any other name.
    public static func format(forFileName fileName: String) -> Format? {
        let base = URL(fileURLWithPath: fileName).lastPathComponent
        let ext = URL(fileURLWithPath: base).pathExtension.lowercased()
        switch ext {
        case "xyz", "extxyz": return .extendedXYZ
        case "vasp", "poscar", "contcar": return .poscar
        default: break
        }
        let upper = base.uppercased()
        if upper.contains("POSCAR") || upper.contains("CONTCAR") { return .poscar }
        return nil
    }

    /// Reads an atom model from a file's bytes, the format from the file's name (`format(forFileName:)`). Throws
    /// `AtomModelError`.
    public static func read(bytes: Data, fileName: String) throws -> AtomModel {
        guard let format = format(forFileName: fileName) else { throw AtomModelError.unknownFormat(fileName: fileName) }
        return try read(bytes: bytes, format: format)
    }

    /// Reads an atom model from a file's bytes in the given format. Throws `AtomModelError`: the size and NUL checks
    /// first, then what the reader finds.
    public static func read(bytes: Data, format: Format) throws -> AtomModel {
        guard bytes.count <= maximumBytes else { throw AtomModelError.tooLarge(bytes: bytes.count) }
        let text = AtomModelText(bytes: bytes)
        guard !text.holdsNUL else { throw AtomModelError.notText }
        switch format {
        case .extendedXYZ: return try readXYZ(text)
        case .poscar: return try readPOSCAR(text)
        }
    }

    /// Reads an atom model from text (UTF-8), in the given format.
    public static func read(_ text: String, format: Format) throws -> AtomModel {
        try read(bytes: Data(text.utf8), format: format)
    }
}
