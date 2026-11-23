//
//  CIFCorpusFuzzTests.swift
//  Role: the reader never crashes (DiffractoKit #117, plan row K4). The corpus is every CIF the package commits (the
//        structure fixtures of all test targets and the edge cases of `cif-edge-cases/`); each is read whole, then mutated
//        by a seeded generator (numbers swapped for nan, inf, 1e999, 1e-320, ? and .; lines dropped, doubled, swapped and
//        cut; keywords and stray bytes inserted). A read may throw anything; a trap shows as a crashed test process.
//        Deterministic: the seed is fixed, so a failure reproduces.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import DiffractoCore

struct CIFCorpusFuzzTests {

    /// SplitMix64: a fixed, portable sequence.
    struct Generator: RandomNumberGenerator {
        var state: UInt64
        mutating func next() -> UInt64 {
            state &+= 0x9E37_79B9_7F4A_7C15
            var z = state
            z = (z ^ (z >> 30)) &* 0xBF58_476D_1CE4_E5B9
            z = (z ^ (z >> 27)) &* 0x94D0_49BB_1331_11EB
            return z ^ (z >> 31)
        }
    }

    /// Every committed `.cif` below Tests/, sorted by path so the run is the same everywhere.
    static func corpus() throws -> [(name: String, bytes: Data)] {
        let root = FixtureDirs.tests
        guard let walker = FileManager.default.enumerator(at: root, includingPropertiesForKeys: nil) else { return [] }
        var files: [URL] = []
        for case let url as URL in walker where url.pathExtension.lowercased() == "cif" { files.append(url) }
        return try files.sorted { $0.path < $1.path }.map { ($0.lastPathComponent, try Data(contentsOf: $0)) }
    }

    static let replacements = ["nan", "-nan", "inf", "-inf", "1e999", "-1e999", "1e308", "-1e308", "1e-320", "0", "-0", "?",
                               ".", "1e15", "0.5(", "(3)", "''", "1/0"]
    static let insertions = ["loop_", "data_", "data_x", "_cell_length_a", "_atom_site_fract_x", "'", "\"", ";", "#",
                             "_symmetry_equiv_pos_as_xyz", "'x,y,z,'", "global_", "save_", "\u{FEFF}", "\t", "_"]

    static func mutate(_ text: String, using rng: inout Generator) -> String {
        var lines = text.components(separatedBy: "\n")
        guard !lines.isEmpty else { return text }
        let edits = Int.random(in: 1...4, using: &rng)
        for _ in 0..<edits {
            let i = Int.random(in: 0..<max(1, lines.count), using: &rng)
            switch Int.random(in: 0..<8, using: &rng) {
            case 0, 1, 2:
                // swap one number of the line for a hostile one
                var tokens = lines[i].split(separator: " ", omittingEmptySubsequences: false).map(String.init)
                let numeric = tokens.indices.filter { Double(tokens[$0].split(separator: "(").first.map(String.init) ?? "") != nil }
                if let k = numeric.randomElement(using: &rng) {
                    tokens[k] = replacements.randomElement(using: &rng)!
                    lines[i] = tokens.joined(separator: " ")
                }
            case 3: lines.remove(at: i); if lines.isEmpty { lines = [""] }
            case 4: lines.insert(lines[i], at: i)
            case 5: lines.swapAt(i, Int.random(in: 0..<lines.count, using: &rng))
            case 6: lines.insert(insertions.randomElement(using: &rng)!, at: i)
            default:
                // cut the file in the middle of this line
                let cut = Int.random(in: 0...lines[i].count, using: &rng)
                lines[i] = String(lines[i].prefix(cut))
                lines.removeSubrange((i + 1)...)
            }
            if lines.isEmpty { lines = [""] }
        }
        return lines.joined(separator: "\n")
    }

    /// A read that returns or throws is a pass; what is read is used the way an app uses it.
    static func exercise(_ text: String) {
        guard let doc = try? CIFImport.read(text) else { return }
        _ = doc.structure.formula
        _ = doc.structure.validationIssues
        _ = doc.warnings.map(\.description)
    }

    @Test func theCorpusIsReadOrRefusedNeverCrashes() throws {
        let corpus = try Self.corpus()
        // 26 structure fixtures and 9 edge cases on 0.19: the walk must find them, or the test proves nothing.
        #expect(corpus.count >= 35, "\(corpus.map(\.name))")
        for file in corpus {
            if let loaded = try? CrystalLoading.load(bytes: file.bytes, fileName: file.name) {
                _ = loaded.structure.formula
            }
        }
    }

    @Test func mutatedCorpusFilesNeverCrashTheReader() throws {
        var rng = Generator(state: 0x4B34_0117)   // "K4", #117
        var reads = 0
        for file in try Self.corpus() {
            guard let text = CrystalLoading.text(of: file.bytes) else { continue }
            for _ in 0..<40 {
                Self.exercise(Self.mutate(text, using: &rng))
                reads += 1
            }
        }
        #expect(reads >= 35 * 40)
    }

    /// Bytes, not text: flipped and inserted bytes go through the loader's encoding step too.
    @Test func mutatedBytesNeverCrashTheLoader() throws {
        var rng = Generator(state: 0xB17E5)
        for file in try Self.corpus() where !file.bytes.isEmpty {
            for _ in 0..<10 {
                var bytes = [UInt8](file.bytes)
                for _ in 0..<Int.random(in: 1...8, using: &rng) {
                    let i = Int.random(in: 0..<bytes.count, using: &rng)
                    if Bool.random(using: &rng) { bytes[i] = UInt8.random(in: 0...255, using: &rng) }
                    else { bytes.insert(UInt8.random(in: 0x80...0xFF, using: &rng), at: i) }
                }
                if let loaded = try? CrystalLoading.load(bytes: Data(bytes), fileName: file.name) {
                    _ = loaded.structure.formula
                }
            }
        }
    }

    /// Symmetry-operation text on its own: hostile numbers and a raised magnitude bound never trap.
    @Test func hostileOperationTextIsRefusedNotTrapped() {
        var rng = Generator(state: 0x5E7)
        let parts = ["x", "-y", "z", "1/2", "+", "-", "1e300", "100000000000000000000x", "nan", "inf", "0/0", "1/0", ",", ".5",
                     "3x", "x+y+z+x+y+z"]
        for _ in 0..<2000 {
            let text = (0..<Int.random(in: 1...9, using: &rng)).map { _ in parts.randomElement(using: &rng)! }.joined()
            _ = try? SymmetryOperation(xyz: text, maximumMagnitude: 1e300)
            _ = try? SymmetryOperation(xyz: text)
        }
        // The case of #117: a coefficient past Int's range under a raised bound; and coefficients whose wrapping
        // determinant is Int.min (2^63), whose abs(_:) trapped.
        #expect(throws: SymmetryOperationError.self) {
            try SymmetryOperation(xyz: "100000000000000000000x,y,z", maximumMagnitude: 1e300)
        }
        #expect(throws: SymmetryOperationError.self) {
            try SymmetryOperation(xyz: "2097152x,2097152y,2097152z", maximumMagnitude: 1e300)
        }
    }
}
