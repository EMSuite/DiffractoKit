//
//  AtomModel+POSCAR.swift
//  Role: the VASP 5 POSCAR / CONTCAR reader, line for line as ASE 3.29.0 `read_vasp_configuration` (ase/io/vasp.py 167)
//        reads it: the comment line; one scale factor (a negative one is the cell's volume) or three positive ones; the
//        three cell vectors; the element line (a POTCAR label such as `Na_pv/6a2f546d` gives `Na`); the counts (a `!`
//        starts a comment); `Selective dynamics` (its flags are ignored); Direct or Cartesian (first letter C or K);
//        one position per atom. What follows the positions (velocities, an MD block) is not read.
//
//  Departures, each a refusal by name: no element line (VASP 4: ASE takes the elements from a POTCAR or OUTCAR beside
//  the file) is `noElementLine`; a cell of zero volume is `degenerateCell`.
//

import Foundation

extension AtomModel {

    static func readPOSCAR(_ source: AtomModelText) throws -> AtomModel {
        var text = source
        guard let titleLine = text.next() else { throw AtomModelError.noAtoms }
        let comment = text.string(text.trimmed(titleLine))

        // The scale: up to three numbers, read until one is not a number (ASE's parse_poscar_scaling_factor).
        guard let scaleLine = text.next() else {
            throw AtomModelError.malformedLine(line: 2, expected: "the scale factor")
        }
        let scaleLineNumber = text.lineNumber
        var scale: [Double] = []
        for field in text.fields(scaleLine, limit: 3) {
            guard let value = text.double(field) else { break }
            scale.append(value)
        }
        guard scale.count == 1 || scale.count == 3 else { throw AtomModelError.badScale(line: scaleLineNumber) }
        if scale.count == 3, scale[0] < 0 || scale[1] < 0 || scale[2] < 0 {
            throw AtomModelError.badScale(line: scaleLineNumber)
        }

        // The cell, as written.
        var written: [SIMD3<Double>] = []
        for _ in 0..<3 {
            written.append(try vector(&text, expected: "a cell vector (three numbers)"))
        }

        // A negative scale is the volume: the factor that gives the cell that volume (np.cbrt(-scale / det(cell))).
        let factor: SIMD3<Double>
        if scale.count == 1 {
            var s = scale[0]
            if s < 0 {
                let volume = determinant(written)
                guard volume != 0 else { throw AtomModelError.degenerateCell }
                s = cbrt(-s / volume)
                guard s.isFinite else { throw AtomModelError.degenerateCell }
            }
            factor = SIMD3<Double>(s, s, s)
        } else {
            factor = SIMD3<Double>(scale[0], scale[1], scale[2])
        }
        let cell = written.map { $0 * factor }
        let lengthA = (cell[0] * cell[0]).sum().squareRoot()
        let lengthB = (cell[1] * cell[1]).sum().squareRoot()
        let lengthC = (cell[2] * cell[2]).sum().squareRoot()
        let lengths = lengthA * lengthB * lengthC
        let volume = determinant(cell)
        guard volume.isFinite, abs(volume) > 1e-10 * lengths else { throw AtomModelError.degenerateCell }

        // The element line (VASP 5). A line whose first field is an integer is the counts line of VASP 4.
        guard let elementLine = text.next() else {
            throw AtomModelError.malformedLine(line: text.lineNumber + 1, expected: "the element symbols")
        }
        let elementLineNumber = text.lineNumber
        let elementFields = text.fields(elementLine)
        guard let firstElement = elementFields.first else {
            throw AtomModelError.malformedLine(line: elementLineNumber, expected: "the element symbols")
        }
        if text.integer(firstElement) != nil { throw AtomModelError.noElementLine }
        let symbols: [String] = elementFields.map { field in
            let label = text.string(field)
            let potcar = label.split(separator: "/", maxSplits: 1, omittingEmptySubsequences: false).first ?? ""
            let element = potcar.split(separator: "_", maxSplits: 1, omittingEmptySubsequences: false).first ?? ""
            return String(element)
        }

        // The counts; a field holding `!` and all after it are a comment.
        guard let countsLine = text.next() else {
            throw AtomModelError.malformedLine(line: elementLineNumber + 1, expected: "the number of atoms of each element")
        }
        let countsLineNumber = text.lineNumber
        var counts: [Int] = []
        for field in text.fields(countsLine) {
            if text.bytes[field].contains(33) { break }
            guard let n = text.integer(field), n >= 0 else {
                throw AtomModelError.malformedLine(line: countsLineNumber, expected: "the number of atoms of each element")
            }
            counts.append(n)
        }
        guard !counts.isEmpty else {
            throw AtomModelError.malformedLine(line: countsLineNumber, expected: "the number of atoms of each element")
        }
        guard counts.count <= symbols.count else {
            throw AtomModelError.countsWithoutElements(line: countsLineNumber, elements: symbols.count, counts: counts.count)
        }
        var total = 0
        for n in counts {
            total += n
            guard total <= 1_000_000_000 else {
                throw AtomModelError.malformedLine(line: countsLineNumber, expected: "the number of atoms of each element")
            }
        }
        guard total > 0 else { throw AtomModelError.noAtoms }
        var speciesNumbers: [Int] = []
        for index in counts.indices {
            guard let z = Element.atomicNumber(forSymbol: symbols[index]) else {
                throw AtomModelError.unknownElement(symbol: symbols[index], line: elementLineNumber)
            }
            speciesNumbers.append(z)
        }

        // Selective dynamics (first letter S), then Direct or Cartesian (first letter C or K is Cartesian).
        var modeLine = try nonBlankTrimmedLine(&text, expected: "Direct or Cartesian")
        if text.bytes[modeLine.lowerBound] == 83 || text.bytes[modeLine.lowerBound] == 115 {
            modeLine = try nonBlankTrimmedLine(&text, expected: "Direct or Cartesian")
        }
        let mode = text.bytes[modeLine.lowerBound]
        let cartesian = mode == 67 || mode == 99 || mode == 75 || mode == 107

        // One position per atom; each atom's element follows the counts. Nothing is reserved beyond what the file's
        // bytes can hold, so an announced count far beyond the file is refused without allocating for it.
        let capacity = min(total, text.bytes.count / 6 + 1)
        var atomicNumbers: [Int] = []
        var positions: [SIMD3<Double>] = []
        atomicNumbers.reserveCapacity(capacity)
        positions.reserveCapacity(capacity)
        var species = 0
        var left = counts[0]
        for k in 0..<total {
            guard !text.atEnd else { throw AtomModelError.missingAtoms(expected: total, found: k) }
            while left == 0 {
                species += 1
                left = counts[species]
            }
            atomicNumbers.append(speciesNumbers[species])
            left -= 1
            let p = try vector(&text, expected: "a position (three numbers)")
            if cartesian {
                positions.append(p * factor)
            } else {
                let alongA: SIMD3<Double> = p.x * cell[0]
                let alongB: SIMD3<Double> = p.y * cell[1]
                let alongC: SIMD3<Double> = p.z * cell[2]
                positions.append(alongA + alongB + alongC)
            }
        }
        return AtomModel(format: .poscar, atomicNumbers: atomicNumbers, positions: positions, cell: cell,
                         periodic: [true, true, true], frameCount: 1, comment: comment)
    }

    /// The first three numbers of the next line.
    private static func vector(_ text: inout AtomModelText, expected: String) throws -> SIMD3<Double> {
        guard let line = text.next() else {
            throw AtomModelError.malformedLine(line: text.lineNumber + 1, expected: expected)
        }
        let fields = text.fields(line, limit: 3)
        guard fields.count == 3 else { throw AtomModelError.malformedLine(line: text.lineNumber, expected: expected) }
        var values: [Double] = [0, 0, 0]
        for k in 0..<3 {
            guard let v = text.double(fields[k]) else {
                throw AtomModelError.badNumber(text: text.string(fields[k]), line: text.lineNumber)
            }
            values[k] = v
        }
        return SIMD3<Double>(values[0], values[1], values[2])
    }

    /// The next line without surrounding whitespace; refused when it is missing or blank.
    private static func nonBlankTrimmedLine(_ text: inout AtomModelText, expected: String) throws -> Range<Int> {
        guard let line = text.next() else {
            throw AtomModelError.malformedLine(line: text.lineNumber + 1, expected: expected)
        }
        let trimmed = text.trimmed(line)
        guard !trimmed.isEmpty else { throw AtomModelError.malformedLine(line: text.lineNumber, expected: expected) }
        return trimmed
    }

    /// a · (b × c) of three row vectors.
    static func determinant(_ rows: [SIMD3<Double>]) -> Double {
        let a = rows[0], b = rows[1], c = rows[2]
        let cross = SIMD3<Double>(b.y * c.z - b.z * c.y, b.z * c.x - b.x * c.z, b.x * c.y - b.y * c.x)
        return (a * cross).sum()
    }
}
