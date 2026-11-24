//
//  AtomModel+Text.swift
//  Role: the byte scanner the XYZ and POSCAR readers share: lines as Python's universal newlines split them (\n, \r\n
//        and a lone \r), fields as `str.split()` splits them (ASCII whitespace), and numbers as Python's `float()` and
//        `int()` read them, except that a number must be finite and is written in decimal (no `nan`, `inf`, hex or
//        `1_000`). It works on the file's bytes, not on a String, so a model of a million atoms reads in well under a
//        second.
//

import Foundation

struct AtomModelText {
    /// The file's bytes.
    let bytes: [UInt8]
    /// The first byte after the current line.
    var position = 0
    /// The number of the line `next()` returned last, counted from 1 (0 before the first).
    var lineNumber = 0

    init(bytes data: Data) {
        bytes = [UInt8](data)
    }

    /// Whether the bytes hold a NUL byte.
    var holdsNUL: Bool { bytes.contains(0) }

    /// Whether every line has been read.
    var atEnd: Bool { position >= bytes.count }

    /// The next line's byte range, without its line break; nil at the end of the file. A last line without a line
    /// break is a line.
    mutating func next() -> Range<Int>? {
        let end = bytes.count
        guard position < end else { return nil }
        let start = position
        var i = start
        while i < end {
            let b = bytes[i]
            if b == 10 || b == 13 { break }
            i += 1
        }
        let line = start..<i
        if i < end {
            if bytes[i] == 13, i + 1 < end, bytes[i + 1] == 10 { i += 2 } else { i += 1 }
        }
        position = i
        lineNumber += 1
        return line
    }

    /// Python's `str.isspace()` for one ASCII byte: space, \t, \n, \v, \f, \r and the separators 0x1C ... 0x1F.
    static func isSpace(_ b: UInt8) -> Bool {
        b == 32 || (b >= 9 && b <= 13) || (b >= 28 && b <= 31)
    }

    /// The whitespace-separated fields of a line (`str.split()`), at most `limit` of them (the rest are not looked at).
    func fields(_ line: Range<Int>, limit: Int = Int.max) -> [Range<Int>] {
        var result: [Range<Int>] = []
        var i = line.lowerBound
        let end = line.upperBound
        while i < end, result.count < limit {
            while i < end, AtomModelText.isSpace(bytes[i]) { i += 1 }
            if i >= end { break }
            let start = i
            while i < end, !AtomModelText.isSpace(bytes[i]) { i += 1 }
            result.append(start..<i)
        }
        return result
    }

    /// Whether a line holds only whitespace (`line.strip() == ''`).
    func isBlank(_ line: Range<Int>) -> Bool {
        for i in line where !AtomModelText.isSpace(bytes[i]) { return false }
        return true
    }

    /// The line without surrounding whitespace.
    func trimmed(_ line: Range<Int>) -> Range<Int> {
        var lower = line.lowerBound
        var upper = line.upperBound
        while lower < upper, AtomModelText.isSpace(bytes[lower]) { lower += 1 }
        while upper > lower, AtomModelText.isSpace(bytes[upper - 1]) { upper -= 1 }
        return lower..<upper
    }

    /// The bytes of a range as text (UTF-8; a byte that is not UTF-8 becomes U+FFFD).
    func string(_ range: Range<Int>) -> String {
        String(decoding: bytes[range], as: UTF8.self)
    }

    /// A finite decimal number (Python's `float()` on the field), or nil. Digits, one sign, a point and an exponent
    /// only, so `nan`, `inf`, `0x1p3` and `1.0d0` are not numbers; `1e999` is not finite and is not one either.
    func double(_ field: Range<Int>) -> Double? {
        guard !field.isEmpty else { return nil }
        for i in field {
            let b = bytes[i]
            let allowed = (b >= 48 && b <= 57) || b == 43 || b == 45 || b == 46 || b == 69 || b == 101
            if !allowed { return nil }
        }
        guard let value = Double(String(decoding: bytes[field], as: UTF8.self)), value.isFinite else { return nil }
        return value
    }

    /// An integer (Python's `int()` on the field: an optional sign, then digits), or nil; nil beyond ±10^15.
    func integer(_ field: Range<Int>) -> Int? {
        var i = field.lowerBound
        let end = field.upperBound
        guard i < end else { return nil }
        var negative = false
        if bytes[i] == 43 || bytes[i] == 45 {
            negative = bytes[i] == 45
            i += 1
        }
        guard i < end else { return nil }
        var value = 0
        while i < end {
            let b = bytes[i]
            guard b >= 48, b <= 57 else { return nil }
            value = value * 10 + Int(b - 48)
            if value > 1_000_000_000_000_000 { return nil }
            i += 1
        }
        return negative ? -value : value
    }
}
