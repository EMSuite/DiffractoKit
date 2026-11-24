//
//  CIFExport.swift
//  Role: a crystal structure written as CIF (DDL1 tags, as VESTA, pymatgen and CrystalMaker read and write them), with the
//        symmetry the structure holds: its space-group number and Hermann–Mauguin symbol, the operation list, one site per
//        orbit, occupancy and U_iso. Numbers are Swift's shortest exact form (`Double.description`), so a value read back is
//        the value written, to the bit. Before it returns, the writer reads its own text back with `CIFImport`; if the
//        expanded sites differ from the structure's in any bit, it writes P1 with every site instead and says so in the
//        file and in the outcome. A wrong asymmetric unit cannot leave the kit.
//        Which operations: the file's own loop when the structure came from one and every listed operation is a symmetry
//        (`fileOperations`, `fileOperationsAreSymmetries`), else the structure's operations. The file's loop is preferred
//        because it built the sites: expanding the same base sites with the same operations repeats the import's arithmetic.
//        Pins: C1 (pymatgen-core 2026.10.16 `CifParser` reads the text to the structure's cell and sites) and C2 (the kit's
//        own round trip, equal to the bit), suite design "Files, figures and exports".
//

import Foundation

/// Why a structure could not be written as CIF.
public enum CIFExportError: Error, Equatable, Sendable, CustomStringConvertible, LocalizedError {
    /// A site's atomic number names no element.
    case unknownElement(z: Int)
    /// A cell length or angle, a coordinate, an occupancy or a displacement is NaN or infinite.
    case nonFiniteValue(String)
    /// Even the P1 text, every site listed, does not read back to the structure: `reason` says what differs. A structure
    /// built in code that `CIFImport` would refuse (two atoms closer than its limit), or one with a coordinate outside
    /// [0, 1) (the reader wraps it), ends here.
    case readBackDiffers(String)

    public var description: String {
        switch self {
        case let .unknownElement(z): "No element has the atomic number \(z), so the site cannot be written as CIF."
        case let .nonFiniteValue(what): "\(what) is not a finite number, so the structure cannot be written as CIF."
        case let .readBackDiffers(reason): "The CIF written does not read back to the structure: \(reason)."
        }
    }

    public var errorDescription: String? { description }
}

public enum CIFExport {
    /// What the file carries besides the structure.
    public struct Options: Sendable, Equatable {
        /// The `data_` block's name; nil takes the structure's name, else its formula. Spaces and characters a CIF name
        /// cannot hold become `_`.
        public var blockName: String?
        /// Free lines written as `#` comments at the head of the file: the methods record's paragraph, "lattice edited".
        public var comments: [String]
        /// Where the structure came from (a file name, a Materials Project id with its date and licence): `#` lines under
        /// "Source:". A source is a comment, never a tag another program would read as a database code.
        public var sourceLines: [String]
        /// The program that wrote the file, as `_audit_creation_method`.
        public var creator: String

        public init(blockName: String? = nil, comments: [String] = [], sourceLines: [String] = [], creator: String = "DiffractoKit") {
            self.blockName = blockName
            self.comments = comments
            self.sourceLines = sourceLines
            self.creator = creator
        }
    }

    /// How the symmetry went into the file.
    public enum Symmetry: Sendable, Equatable {
        /// The operation list and one site per orbit: `operations` listed, `sites` written of the structure's `expanded`.
        case kept(operations: Int, sites: Int, expanded: Int)
        /// Every site in P1: the structure holds no operation but the identity (`reason` nil), or the symmetric text did
        /// not read back to the structure (`reason` says what differed; the file says it too).
        case p1(reason: String?)
    }

    public struct Outcome: Sendable, Equatable {
        /// The CIF text, ending in a newline.
        public let text: String
        public let symmetry: Symmetry
    }

    /// The structure as CIF. Throws only when a value cannot be written or the P1 text itself does not read back.
    public static func text(of structure: CrystalStructure, options: Options = Options()) throws(CIFExportError) -> Outcome {
        try checkFinite(structure)
        let labels = try siteLabels(structure.sites)
        if let written = writtenOperations(of: structure) {
            let operations = written.operations
            let orbits = orbitRepresentatives(structure.sites, operations: operations)
            let text = try cifText(structure, options: options, operations: operations, operationsFromFile: written.fromFile,
                                   sites: orbits.representatives, labels: labels, p1Reason: nil)
            // Every image reads back with its orbit's label.
            let expected = orbits.owner.map { labels[$0] }
            switch readBackDifference(text, structure: structure, labels: expected) {
            case .none:
                return Outcome(text: text, symmetry: .kept(operations: operations.count, sites: orbits.representatives.count,
                                                           expanded: structure.sites.count))
            case let .some(reason):
                let p1 = try cifText(structure, options: options, operations: [.identity], operationsFromFile: false,
                                     sites: Array(structure.sites.indices), labels: labels, p1Reason: reason)
                if let p1Difference = readBackDifference(p1, structure: structure, labels: labels) { throw .readBackDiffers(p1Difference) }
                return Outcome(text: p1, symmetry: .p1(reason: reason))
            }
        }
        let p1 = try cifText(structure, options: options, operations: [.identity], operationsFromFile: false,
                             sites: Array(structure.sites.indices), labels: labels, p1Reason: nil)
        if let difference = readBackDifference(p1, structure: structure, labels: labels) { throw .readBackDiffers(difference) }
        return Outcome(text: p1, symmetry: .p1(reason: nil))
    }

    // MARK: - Symmetry

    /// The operations the file lists, identity first, and whether they are the file's loop; nil when the structure holds the
    /// identity alone.
    static func writtenOperations(of structure: CrystalStructure) -> (operations: [SymmetryOperation], fromFile: Bool)? {
        let fromFile = structure.fileOperationsAreSymmetries == true
            && (structure.fileOperations ?? []).contains(where: { !$0.isIdentity })
        let source = fromFile ? structure.fileOperations ?? [] : structure.operations
        var distinct: [SymmetryOperation] = [.identity]
        for op in source where !distinct.contains(where: { $0.isEquivalent(to: op) }) { distinct.append(op) }
        return distinct.count > 1 ? (distinct, fromFile) : nil
    }

    /// One site of each orbit, the first in the structure's order: a site joins an earlier one's orbit when it has the same
    /// element, occupancy and displacement and sits within 1e-6 (fractional, per axis, modulo 1) of one of its images. Its
    /// label is not compared: an image takes its representative's label, as in any CIF of an asymmetric unit, so the
    /// Si2 ... Si8 of a P1 file read back as Si1.
    /// The read-back decides whether the choice was right; this only has to be right for the structures it can be.
    /// `owner[j]` is the representative of site j's orbit.
    static func orbitRepresentatives(_ sites: [AtomSite], operations: [SymmetryOperation]) -> (representatives: [Int], owner: [Int]) {
        var owner = Array(sites.indices)
        var covered = [Bool](repeating: false, count: sites.count)
        var representatives: [Int] = []
        for (i, site) in sites.enumerated() where !covered[i] {
            representatives.append(i)
            covered[i] = true
            let images = operations.map { wrapFractional($0.apply(site.fractional)) }
            for j in (i + 1)..<max(i + 1, sites.count) where !covered[j] {
                let other = sites[j]
                guard other.z == site.z, other.occupancy == site.occupancy,
                      other.isotropicDisplacement == site.isotropicDisplacement else { continue }
                if images.contains(where: { image in
                    (0..<3).allSatisfy { fractionalAxisDelta(image[$0], other.fractional[$0]) <= 1e-6 }
                }) {
                    covered[j] = true
                    owner[j] = i
                }
            }
        }
        return (representatives, owner)
    }

    /// What differs between the structure and `text` read back by `CIFImport` (cell, and the expanded sites as a multiset:
    /// element, coordinates, occupancy, displacement and label, each equal); nil when nothing does.
    static func readBackDifference(_ text: String, structure: CrystalStructure, labels: [String]) -> String? {
        let read: CrystalStructure
        do { read = try CIFImport.read(text).structure } catch { return "the text does not read (\(error))" }
        let a = structure.cell, b = read.cell
        guard a.a == b.a, a.b == b.b, a.c == b.c, a.alphaDeg == b.alphaDeg, a.betaDeg == b.betaDeg, a.gammaDeg == b.gammaDeg else {
            return "the cell read back is not the structure's"
        }
        guard read.sites.count == structure.sites.count else {
            return "the operations expand to \(read.sites.count) sites, the structure has \(structure.sites.count)"
        }
        let written = zip(structure.sites, labels).map { SiteKey($0, label: $1) }.sorted()
        let back = read.sites.map { SiteKey($0, label: $0.label ?? "") }.sorted()
        guard written == back else { return "the expanded sites are not the structure's to the bit" }
        return nil
    }

    /// A site as the read-back compares it, ordered for a multiset comparison.
    struct SiteKey: Comparable {
        let z: Int
        let position: [Double]
        let occupancy: Double
        let displacement: Double
        let label: String

        init(_ site: AtomSite, label: String) {
            z = site.z
            position = [site.fractional.x, site.fractional.y, site.fractional.z]
            occupancy = site.occupancy
            displacement = site.isotropicDisplacement ?? -.infinity
            self.label = label
        }

        static func < (l: SiteKey, r: SiteKey) -> Bool {
            if l.z != r.z { return l.z < r.z }
            if l.position != r.position { return l.position.lexicographicallyPrecedes(r.position) }
            if l.occupancy != r.occupancy { return l.occupancy < r.occupancy }
            if l.displacement != r.displacement { return l.displacement < r.displacement }
            return l.label < r.label
        }
    }

    // MARK: - Text

    private static func checkFinite(_ structure: CrystalStructure) throws(CIFExportError) {
        let c = structure.cell
        for (name, value) in [("a", c.a), ("b", c.b), ("c", c.c), ("alpha", c.alphaDeg), ("beta", c.betaDeg), ("gamma", c.gammaDeg)]
            where !value.isFinite {
            throw .nonFiniteValue("The cell's \(name)")
        }
        for (i, site) in structure.sites.enumerated() {
            let values = [site.fractional.x, site.fractional.y, site.fractional.z, site.occupancy, site.isotropicDisplacement ?? 0]
            if !values.allSatisfy(\.isFinite) { throw .nonFiniteValue("A value of site \(i + 1)") }
            if Element.symbol(forAtomicNumber: site.z) == nil { throw .unknownElement(z: site.z) }
        }
    }

    /// Each site's label as written: its own when it is one CIF token, else its element and its number in the file
    /// (`Al3`). Labels are not made unique: the reader does not need them to be, and a file's own labels are kept.
    static func siteLabels(_ sites: [AtomSite]) throws(CIFExportError) -> [String] {
        var labels: [String] = []
        for (i, site) in sites.enumerated() {
            guard let symbol = Element.symbol(forAtomicNumber: site.z) else { throw .unknownElement(z: site.z) }
            if let own = site.label, isBareToken(own) { labels.append(own) } else { labels.append("\(symbol)\(i + 1)") }
        }
        return labels
    }

    /// Whether `text` can stand unquoted as a CIF value: no whitespace or quote, not starting with a character or word CIF
    /// reserves, and not `?` or `.`.
    static func isBareToken(_ text: String) -> Bool {
        guard let first = text.first, text != "?", text != ".",
              !"_$'\";[]".contains(first),
              text.allSatisfy({ character in
                  (character.asciiValue ?? 0) > 32 && character != "'" && character != "\"" && character != "#"
              })
        else { return false }
        let lower = text.lowercased()
        return !["data_", "loop_", "save_", "global_", "stop_"].contains { lower.hasPrefix($0) }
    }

    /// `text` as one CIF value: bare when it can be, else in single quotes, else in double quotes; nil when it holds both
    /// quotes or a line break.
    static func value(_ text: String) -> String? {
        if isBareToken(text) { return text }
        guard !text.contains(where: \.isNewline) else { return nil }
        if !text.contains("'") { return "'\(text)'" }
        if !text.contains("\"") { return "\"\(text)\"" }
        return nil
    }

    /// The `data_` block's name: one token of letters, digits and `-_.()`.
    static func blockName(_ structure: CrystalStructure, options: Options) -> String {
        let raw = options.blockName ?? structure.name ?? structure.formula
        let allowed = Set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_.()")
        let name = String(raw.map { allowed.contains($0) ? $0 : "_" })
        return name.isEmpty ? "structure" : name
    }

    /// A number as Swift prints it shortest (`Double(text)` gives the same double), always positional: `0.00005`, never
    /// `5e-05`, because the reader judges a coordinate's precision by its decimals and reads no exponent in an operation.
    static func number(_ value: Double) -> String {
        let shortest = "\(value)"
        guard let e = shortest.firstIndex(where: { $0 == "e" || $0 == "E" }), let exponent = Int(shortest[shortest.index(after: e)...])
        else { return shortest }
        var mantissa = String(shortest[..<e])
        let negative = mantissa.hasPrefix("-")
        if negative { mantissa.removeFirst() }
        let dot = mantissa.firstIndex(of: ".")
        let integerDigits = dot.map { mantissa.distance(from: mantissa.startIndex, to: $0) } ?? mantissa.count
        var digits = mantissa.filter(\.isNumber)
        let point = integerDigits + exponent
        var text: String
        if point <= 0 {
            text = "0." + String(repeating: "0", count: -point) + digits
        } else if point >= digits.count {
            text = digits + String(repeating: "0", count: point - digits.count) + ".0"
        } else {
            digits.insert(".", at: digits.index(digits.startIndex, offsetBy: point))
            text = digits
        }
        // Trailing zeros after the point say nothing ("0.000050" is 0.00005); keep one digit after it.
        while text.hasSuffix("0"), let p = text.firstIndex(of: "."), text.distance(from: p, to: text.endIndex) > 2 { text.removeLast() }
        return (negative ? "-" : "") + text
    }

    /// An operation as `x,y,z` text whose translation reads back to the bit: a fraction n/d when `Double(n) / Double(d)` is
    /// the translation exactly (the reader divides the same way), else the shortest decimal.
    static func xyz(_ op: SymmetryOperation) -> String {
        let names = ["x", "y", "z"]
        var parts: [String] = []
        for axis in 0..<3 {
            var s = ""
            let row = op.rotation.rows[axis]
            for (i, coefficient) in [row.x, row.y, row.z].enumerated() where coefficient != 0 {
                let magnitude = abs(coefficient) == 1 ? "" : "\(abs(coefficient))"
                s += (coefficient < 0 ? "-" : (s.isEmpty ? "" : "+")) + magnitude + names[i]
            }
            let t = op.translation[axis]
            if t != 0 {
                s += (t < 0 ? "-" : (s.isEmpty ? "" : "+")) + fraction(abs(t))
            }
            parts.append(s.isEmpty ? "0" : s)
        }
        return parts.joined(separator: ",")
    }

    private static func fraction(_ t: Double) -> String {
        for denominator in [2, 3, 4, 6, 8, 12, 24] {
            let numerator = (t * Double(denominator)).rounded()
            if numerator >= 1, numerator < 1000, numerator / Double(denominator) == t { return "\(Int(numerator))/\(denominator)" }
        }
        return number(t)
    }

    private static func cifText(_ structure: CrystalStructure, options: Options, operations: [SymmetryOperation],
                                operationsFromFile: Bool, sites: [Int], labels: [String],
                                p1Reason: String?) throws(CIFExportError) -> String {
        var lines: [String] = []
        for comment in options.comments { lines.append(contentsOf: commentLines(comment)) }
        if !options.sourceLines.isEmpty {
            lines.append("# Source:")
            for source in options.sourceLines { lines.append(contentsOf: commentLines("  " + source)) }
        }
        if let p1Reason {
            lines.append(contentsOf: commentLines("Written in P1 with every site: the symmetric text did not read back to the structure (\(p1Reason))."))
        }
        lines.append("data_\(blockName(structure, options: options))")
        if let creator = value(options.creator) { lines.append("_audit_creation_method \(creator)") }
        if let name = structure.name, let quoted = value(name) { lines.append("_chemical_name_common \(quoted)") }
        if let formula = value(structure.formula) { lines.append("_chemical_formula_sum \(formula)") }
        let c = structure.cell
        lines.append("_cell_length_a \(number(c.a))")
        lines.append("_cell_length_b \(number(c.b))")
        lines.append("_cell_length_c \(number(c.c))")
        lines.append("_cell_angle_alpha \(number(c.alphaDeg))")
        lines.append("_cell_angle_beta \(number(c.betaDeg))")
        lines.append("_cell_angle_gamma \(number(c.gammaDeg))")
        if operations.count == 1 {
            lines.append("_space_group_name_H-M_alt 'P 1'")
            lines.append("_space_group_IT_number 1")
        } else {
            let type = SpaceGroupType.of(operations)
            // The symbol names the standard setting: written only when the operations are in it.
            if let type, type.settingHallNumber == type.hallNumber, let symbol = value(type.internationalSymbol) {
                lines.append("_space_group_name_H-M_alt \(symbol)")
            }
            // The file's loop goes with the number the file stated; detected operations with the number of their own group
            // (a P1 file found to be Fm-3m stated 1, which the reader would hold against 192 operations).
            let number = operationsFromFile ? structure.spaceGroupNumber ?? type?.number : type?.number
            if let number { lines.append("_space_group_IT_number \(number)") }
        }
        lines.append("loop_")
        lines.append("_space_group_symop_operation_xyz")
        for op in operations { lines.append("'\(xyz(op))'") }
        let withDisplacement = sites.contains { structure.sites[$0].isotropicDisplacement != nil }
        lines.append("loop_")
        lines.append(contentsOf: ["_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y",
                                  "_atom_site_fract_z", "_atom_site_occupancy"])
        if withDisplacement { lines.append("_atom_site_U_iso_or_equiv") }
        for i in sites {
            let site = structure.sites[i]
            guard let symbol = Element.symbol(forAtomicNumber: site.z) else { throw .unknownElement(z: site.z) }
            var row = [labels[i], symbol, number(site.fractional.x), number(site.fractional.y), number(site.fractional.z),
                       number(site.occupancy)]
            if withDisplacement { row.append(site.isotropicDisplacement.map(number) ?? "?") }
            lines.append(row.joined(separator: " "))
        }
        return lines.joined(separator: "\n") + "\n"
    }

    /// `text` as `#` lines, one per line of it.
    private static func commentLines(_ text: String) -> [String] {
        text.split(omittingEmptySubsequences: false, whereSeparator: \.isNewline).map { "# " + $0 }
    }
}
