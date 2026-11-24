//
//  AtomModel+XYZ.swift
//  Role: the extended XYZ reader (plain XYZ is its special case), as ASE 3.29.0 `read_xyz` with its default index -1
//        (ase/io/extxyz.py 569) reads a file: the frames are counted (a blank header line ends the file), the last is
//        read. Its comment line is key=value pairs (`key_val_str_to_dict`, line 67: quotes, braces and brackets
//        delimit a value, a backslash escapes, a key without a value is true); `Lattice` is nine numbers, the vectors
//        a, b, c in turn; `pbc` is one or three of T/F (also True, false, ...); `Properties` names the columns
//        (`parse_properties`, line 338; default `species:S:1:pos:R:3`). Species are capitalised (`MG` is Mg), a `Z`
//        column wins over species, a `Lattice` without `pbc` is periodic along all three. A frame followed by `VEC1`
//        ... `VEC3` lines takes its cell from them and its comment is not parsed, as `_read_xyz_frame` (line 390) does.
//        Every R and I column is read as a number, as ASE converts every column.
//

import Foundation

extension AtomModel {

    /// Where a frame starts and what its header announced.
    private struct XYZFrame {
        let position: Int
        let lineNumber: Int
        let atoms: Int
        let vectors: Int
    }

    /// One property of `Properties`: its name, type letter (R, I, S or L), first column and number of columns.
    private struct XYZColumn {
        let name: String
        let type: Character
        let first: Int
        let count: Int
    }

    static func readXYZ(_ source: AtomModelText) throws -> AtomModel {
        var text = source

        // Count the frames, as read_xyz scans the file: header, comment, the atoms' lines, then VEC lines.
        var last: XYZFrame?
        var frameCount = 0
        while true {
            let start = text.position
            let startLine = text.lineNumber
            guard let header = text.next(), !text.isBlank(header) else { break }
            let headerFields = text.fields(header, limit: 2)
            guard headerFields.count == 1, let n = text.integer(headerFields[0]), n >= 0 else {
                throw AtomModelError.malformedLine(line: text.lineNumber, expected: "the number of atoms")
            }
            _ = text.next()
            var skipped = 0
            while skipped < n, text.next() != nil { skipped += 1 }
            var vectors = 0
            while true {
                let before = text.position
                let beforeLine = text.lineNumber
                guard let line = text.next() else { break }
                if startsWithVEC(text, line) {
                    vectors += 1
                    if vectors > 3 {
                        throw AtomModelError.malformedLine(line: text.lineNumber, expected: "at most three VEC lines")
                    }
                } else {
                    text.position = before
                    text.lineNumber = beforeLine
                    break
                }
            }
            last = XYZFrame(position: start, lineNumber: startLine, atoms: n, vectors: vectors)
            frameCount += 1
        }
        guard let frame = last, frame.atoms > 0 else { throw AtomModelError.noAtoms }

        // Read the last frame.
        text.position = frame.position
        text.lineNumber = frame.lineNumber
        _ = text.next()
        var comment = ""
        if let commentLine = text.next() {
            comment = text.string(text.trimmed(commentLine))
        }
        let commentLineNumber = text.lineNumber

        var info: [String: String] = [:]
        if frame.vectors == 0, !comment.isEmpty {
            guard let pairs = keyValues(comment) else { throw AtomModelError.badCommentLine(line: commentLineNumber) }
            info = pairs
        }

        var cell: [SIMD3<Double>] = [SIMD3<Double>(0, 0, 0), SIMD3<Double>(0, 0, 0), SIMD3<Double>(0, 0, 0)]
        var periodic = [false, false, false]
        if let lattice = info["Lattice"] {
            let values = tokens(lattice).compactMap(decimal)
            guard values.count == 9, tokens(lattice).count == 9 else {
                throw AtomModelError.badLattice(line: commentLineNumber)
            }
            for j in 0..<3 {
                cell[j] = SIMD3<Double>(values[3 * j], values[3 * j + 1], values[3 * j + 2])
            }
            periodic = [true, true, true]
        }
        if let pbc = info["pbc"] {
            guard let flags = periodicity(pbc) else { throw AtomModelError.badPeriodicity(line: commentLineNumber) }
            periodic = flags
        }

        // The columns.
        let columns = try properties(info["Properties"] ?? "species:S:1:pos:R:3", line: commentLineNumber)
        let width = columns.reduce(0) { $0 + $1.count }
        var positionColumn: XYZColumn?
        var numberColumn: XYZColumn?
        var speciesColumn: XYZColumn?
        for column in columns {
            switch column.name {
            case "pos", "positions": positionColumn = column
            case "Z", "numbers": numberColumn = column
            case "species", "symbols": speciesColumn = column
            default: break
            }
        }
        guard let pos = positionColumn, pos.type == "R", pos.count == 3 else {
            throw AtomModelError.badProperties(line: commentLineNumber, reason: "no pos column of three reals")
        }
        if let z = numberColumn, z.type != "I" || z.count != 1 {
            throw AtomModelError.badProperties(line: commentLineNumber, reason: "Z is not one integer column")
        }
        if numberColumn == nil {
            guard let s = speciesColumn, s.type == "S", s.count == 1 else {
                throw AtomModelError.badProperties(line: commentLineNumber, reason: "no species or Z column")
            }
        }
        var realColumns: [Int] = []
        var integerColumns: [Int] = []
        for column in columns {
            for c in column.first..<(column.first + column.count) {
                if column.type == "R" { realColumns.append(c) }
                if column.type == "I" { integerColumns.append(c) }
            }
        }

        var atomicNumbers: [Int] = []
        var positions: [SIMD3<Double>] = []
        let capacity = min(frame.atoms, text.bytes.count / 4 + 1)
        atomicNumbers.reserveCapacity(capacity)
        positions.reserveCapacity(capacity)
        var symbolCache: [String: Int] = [:]
        for k in 0..<frame.atoms {
            guard let line = text.next() else { throw AtomModelError.missingAtoms(expected: frame.atoms, found: k) }
            let lineNumber = text.lineNumber
            let fields = text.fields(line, limit: width)
            guard fields.count == width else {
                throw AtomModelError.malformedLine(line: lineNumber, expected: "\(width) columns (Properties)")
            }
            var values: [Double] = [0, 0, 0]
            for c in realColumns {
                guard let v = text.double(fields[c]) else {
                    throw AtomModelError.badNumber(text: text.string(fields[c]), line: lineNumber)
                }
                let offset = c - pos.first
                if offset >= 0, offset < 3 { values[offset] = v }
            }
            for c in integerColumns where text.integer(fields[c]) == nil {
                throw AtomModelError.badNumber(text: text.string(fields[c]), line: lineNumber)
            }
            positions.append(SIMD3<Double>(values[0], values[1], values[2]))

            if let zColumn = numberColumn {
                let field = fields[zColumn.first]
                guard let z = text.integer(field), z >= 1, z <= Element.maximumAtomicNumber else {
                    throw AtomModelError.unknownElement(symbol: text.string(field), line: lineNumber)
                }
                atomicNumbers.append(z)
            } else if let sColumn = speciesColumn {
                let raw = text.string(fields[sColumn.first])
                if let z = symbolCache[raw] {
                    atomicNumbers.append(z)
                } else {
                    guard let z = Element.atomicNumber(forSymbol: capitalised(raw)) else {
                        throw AtomModelError.unknownElement(symbol: raw, line: lineNumber)
                    }
                    symbolCache[raw] = z
                    atomicNumbers.append(z)
                }
            }
        }

        // VEC lines after the atoms give the cell, one vector each, and make that direction periodic.
        for ln in 0..<frame.vectors {
            guard let line = text.next() else {
                throw AtomModelError.malformedLine(line: text.lineNumber + 1, expected: "VEC\(ln + 1)")
            }
            let fields = text.fields(line)
            let label = fields.isEmpty ? "" : text.string(fields[0])
            guard fields.count == 4, label == "VEC\(ln + 1)" else {
                throw AtomModelError.malformedLine(line: text.lineNumber, expected: "VEC\(ln + 1) and three numbers")
            }
            var v: [Double] = [0, 0, 0]
            for k in 0..<3 {
                guard let x = text.double(fields[k + 1]) else {
                    throw AtomModelError.badNumber(text: text.string(fields[k + 1]), line: text.lineNumber)
                }
                v[k] = x
            }
            cell[ln] = SIMD3<Double>(v[0], v[1], v[2])
            periodic[ln] = true
        }

        return AtomModel(format: .extendedXYZ, atomicNumbers: atomicNumbers, positions: positions, cell: cell,
                         periodic: periodic, frameCount: frameCount, comment: comment)
    }

    // MARK: - The comment line

    /// Whether a line, without its leading whitespace, starts with `VEC`.
    private static func startsWithVEC(_ text: AtomModelText, _ line: Range<Int>) -> Bool {
        let t = text.trimmed(line)
        guard t.count >= 3 else { return false }
        let b = text.bytes
        return b[t.lowerBound] == 86 && b[t.lowerBound + 1] == 69 && b[t.lowerBound + 2] == 67
    }

    /// Python's `str.isspace()` for one character.
    static func isPythonSpace(_ s: Unicode.Scalar) -> Bool {
        if s.value >= 0x1C && s.value <= 0x1F { return true }
        return s.properties.isWhitespace
    }

    /// The comment line as key=value pairs, as `key_val_str_to_dict` splits it (values are the raw strings; a key
    /// without a value is "T"). Nil where ASE fails (a line that starts with `=`).
    static func keyValues(_ line: String) -> [String: String]? {
        let delimiters: [Unicode.Scalar: Unicode.Scalar] = ["'": "'", "\"": "\"", "{": "}", "[": "]"]
        let scalars = Array(line.unicodeScalars)
        var lower = 0
        var upper = scalars.count
        while lower < upper, isPythonSpace(scalars[lower]) { lower += 1 }
        while upper > lower, isPythonSpace(scalars[upper - 1]) { upper -= 1 }

        var pairs: [[String]] = [[""]]
        var closing: Unicode.Scalar?
        var escaped = false
        func add(_ ch: Unicode.Scalar) {
            let p = pairs.count - 1
            let e = pairs[p].count - 1
            pairs[p][e].unicodeScalars.append(ch)
        }
        for ch in scalars[lower..<upper] {
            if escaped {
                add(ch)
                escaped = false
            } else if ch == "\\" {
                escaped = true
            } else if let close = closing {
                if ch == close { closing = nil } else { add(ch) }
            } else if let close = delimiters[ch] {
                closing = close
            } else if isPythonSpace(ch) {
                if pairs.count == 1, pairs[0].count == 1, pairs[0][0].isEmpty { continue }
                let lastPair = pairs[pairs.count - 1]
                if lastPair[lastPair.count - 1].isEmpty { continue }
                pairs.append([""])
            } else if ch == "=" {
                let lastPair = pairs[pairs.count - 1]
                if lastPair.count == 1, lastPair[0].isEmpty { pairs.removeLast() }
                if pairs.isEmpty { return nil }
                pairs[pairs.count - 1].append("")
            } else {
                add(ch)
            }
        }
        var result: [String: String] = [:]
        for pair in pairs where !pair.isEmpty {
            if pair.count == 1 {
                result[pair[0]] = "T"
            } else {
                result[pair[0]] = pair[1...].joined(separator: "=")
            }
        }
        return result
    }

    /// The value split as `re.findall(r'[^\s,]+', value)`.
    static func tokens(_ value: String) -> [Substring] {
        value.split(whereSeparator: { ch in
            if ch == "," { return true }
            return ch.unicodeScalars.allSatisfy { isPythonSpace($0) }
        })
    }

    /// A finite decimal number: digits, signs, a point and an exponent only.
    static func decimal(_ token: Substring) -> Double? {
        guard !token.isEmpty else { return nil }
        for u in token.utf8 {
            let allowed = (u >= 48 && u <= 57) || u == 43 || u == 45 || u == 46 || u == 69 || u == 101
            if !allowed { return nil }
        }
        guard let v = Double(token), v.isFinite else { return nil }
        return v
    }

    /// `pbc` as three flags: one or three of T, F, true, false, True, False, TRUE, FALSE, or of numbers (non-zero is
    /// true), as `key_val_str_to_dict` and `Atoms(pbc=...)` read it. Nil for anything else.
    static func periodicity(_ value: String) -> [Bool]? {
        let parts = tokens(value)
        guard parts.count == 1 || parts.count == 3 else { return nil }
        var flags: [Bool] = []
        let numbers = parts.compactMap(decimal)
        if numbers.count == parts.count {
            flags = numbers.map { $0 != 0 }
        } else {
            let words: [String: Bool] = ["T": true, "F": false, "true": true, "false": false, "True": true,
                                         "False": false, "TRUE": true, "FALSE": false]
            for part in parts {
                guard let flag = words[String(part)] else { return nil }
                flags.append(flag)
            }
        }
        return flags.count == 1 ? [flags[0], flags[0], flags[0]] : flags
    }

    /// The columns of a `Properties` value (`parse_properties`): name:type:count triplets.
    private static func properties(_ value: String, line: Int) throws -> [XYZColumn] {
        let parts = value.split(separator: ":", omittingEmptySubsequences: false)
        var columns: [XYZColumn] = []
        var names: Set<String> = []
        var next = 0
        for t in 0..<(parts.count / 3) {
            let name = String(parts[3 * t])
            let typeText = String(parts[3 * t + 1])
            guard let count = Int(parts[3 * t + 2]), count >= 1, count <= 1000 else {
                throw AtomModelError.badProperties(line: line, reason: "the column count of \(name) is not a whole number")
            }
            guard typeText == "R" || typeText == "I" || typeText == "S" || typeText == "L" else {
                throw AtomModelError.badProperties(line: line, reason: "unknown type \(typeText) of \(name)")
            }
            guard !names.contains(name) else {
                throw AtomModelError.badProperties(line: line, reason: "\(name) is named twice")
            }
            names.insert(name)
            columns.append(XYZColumn(name: name, type: Character(typeText), first: next, count: count))
            next += count
            guard next <= 10_000 else {
                throw AtomModelError.badProperties(line: line, reason: "more than 10000 columns")
            }
        }
        return columns
    }

    /// Python's `str.capitalize()` for a symbol: the first letter upper case, the rest lower case.
    static func capitalised(_ symbol: String) -> String {
        guard let first = symbol.first else { return symbol }
        return first.uppercased() + symbol.dropFirst().lowercased()
    }
}
