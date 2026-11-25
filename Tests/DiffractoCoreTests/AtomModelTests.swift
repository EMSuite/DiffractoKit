//
//  AtomModelTests.swift
//  Role: the atom-model readers, extended XYZ and POSCAR (DiffractoKit #169, plan row MD1). Every file of the corpus
//        (ASE 3.29.0's own test inputs, exports by ASE's writers, the layouts of VESTA and atomsk) against what
//        `ase.io.read` returns (`tools/pins/atom-models/pins.py` -> `atom-models-pins.json`): elements and counts exact,
//        numbers the file states exact, converted ones at rtol 1e-7 (assert_allclose's float64 default, as ASE's own
//        reader tests use). Then the refusals by name, a fuzz run that never crashes, and a million atoms in under 2 s.
//

import DiffractoKitTestSupport
import Foundation
import Testing
@testable import DiffractoCore

struct AtomModelTests {

    struct Pins: Decodable {
        struct Case: Decodable {
            let file: String
            let format: String
            let refusal: String?
            let numbers: [Int]?
            let positions: [[Double]]?
            let cell: [[Double]]?
            let pbc: [Bool]?
            let frames: Int?
            let positionsExact: Bool?
            let cellExact: Bool?
        }
        let ase: String
        let rtol: Double
        let cases: [Case]
    }

    static let corpusDirectory = Fixtures.directory.appendingPathComponent("atom-models")

    static func pins() throws -> Pins {
        let url = Fixtures.directory.appendingPathComponent("atom-models-pins.json")
        return try JSONDecoder().decode(Pins.self, from: Data(contentsOf: url))
    }

    static func bytes(_ name: String) throws -> Data {
        try Data(contentsOf: corpusDirectory.appendingPathComponent(name))
    }

    /// assert_allclose(actual, desired, rtol) with atol 0.
    static func close(_ actual: Double, _ desired: Double, rtol: Double) -> Bool {
        abs(actual - desired) <= rtol * abs(desired)
    }

    static func vector(_ v: SIMD3<Double>) -> [Double] { [v.x, v.y, v.z] }

    // MARK: - The corpus against ASE

    @Test func thePinsComeFromTheNamedASE() throws {
        let pins = try Self.pins()
        #expect(pins.ase == "3.29.0")
        #expect(pins.rtol == 1e-7)
        #expect(pins.cases.count >= 25)
        #expect(pins.cases.contains { $0.format == "xyz" } && pins.cases.contains { $0.format == "poscar" })
    }

    @Test func everyCorpusFileEqualsASE() throws {
        let pins = try Self.pins()
        for pin in pins.cases {
            let data = try Self.bytes(pin.file)
            if let refusal = pin.refusal {
                do {
                    _ = try AtomModel.read(bytes: data, fileName: pin.file)
                    Issue.record("\(pin.file) should be refused as \(refusal)")
                } catch let error as AtomModelError {
                    #expect(error.name == refusal, "\(pin.file)")
                }
                continue
            }
            let model = try AtomModel.read(bytes: data, fileName: pin.file)
            #expect(model.format == (pin.format == "xyz" ? .extendedXYZ : .poscar), "\(pin.file)")
            #expect(model.atomicNumbers == pin.numbers, "\(pin.file)")
            #expect(model.periodic == pin.pbc, "\(pin.file)")
            #expect(model.frameCount == pin.frames, "\(pin.file)")
            let positions = pin.positions ?? []
            let cell = pin.cell ?? []
            try #require(model.positions.count == positions.count, "\(pin.file)")
            try #require(cell.count == 3, "\(pin.file)")
            for (i, expected) in positions.enumerated() {
                let actual = Self.vector(model.positions[i])
                if pin.positionsExact == true {
                    #expect(actual == expected, "\(pin.file) atom \(i)")
                } else {
                    for k in 0..<3 {
                        #expect(Self.close(actual[k], expected[k], rtol: pins.rtol) || abs(expected[k]) < 1e-12
                                && abs(actual[k]) < 1e-12, "\(pin.file) atom \(i) axis \(k): \(actual[k]) vs \(expected[k])")
                    }
                }
            }
            for row in 0..<3 {
                let actual = Self.vector(model.cell[row])
                if pin.cellExact == true {
                    #expect(actual == cell[row], "\(pin.file) cell row \(row)")
                } else {
                    for k in 0..<3 {
                        #expect(Self.close(actual[k], cell[row][k], rtol: pins.rtol) || abs(cell[row][k]) < 1e-12
                                && abs(actual[k]) < 1e-12, "\(pin.file) cell \(row) \(k): \(actual[k]) vs \(cell[row][k])")
                    }
                }
            }
        }
    }

    @Test func theCasesThatSayWhatTheFormatPromises() throws {
        // A plain XYZ file has no cell; of several frames the last is read and the count is reported.
        let plain = try AtomModel.read(bytes: Self.bytes("export-si-plain.xyz"), fileName: "export-si-plain.xyz")
        #expect(!plain.hasCell)
        #expect(plain.periodic == [false, false, false])
        let trajectory = try AtomModel.read(bytes: Self.bytes("export-cu-three-frames.xyz"),
                                            fileName: "export-cu-three-frames.xyz")
        #expect(trajectory.frameCount == 3)
        #expect(trajectory.count == 12)
        #expect(trajectory.countsByAtomicNumber == [29: 12])
        // Upper-case species are capitalised (MG is Mg).
        let upper = try AtomModel.read(bytes: Self.bytes("ase-uppercase-blank-line.xyz"), fileName: "x.xyz")
        #expect(upper.atomicNumbers.first == 12)
        // POSCAR: a POTCAR label on the element line, a negative scale is the volume, selective dynamics ignored.
        let hashes = try AtomModel.read(bytes: Self.bytes("ase-poscar-hashes.vasp"), fileName: "POSCAR")
        #expect(hashes.atomicNumbers == [11])
        let negative = try AtomModel.read(bytes: Self.bytes("ase-poscar-negative-scale.vasp"), fileName: "CONTCAR")
        #expect(abs(abs(AtomModel.determinant(negative.cell)) - 120) < 1e-9)
        let selective = try AtomModel.read(bytes: Self.bytes("export-al-vacancy-selective.vasp"), fileName: "a.vasp")
        #expect(selective.count == 107)
        #expect(selective.hasCell && selective.periodic == [true, true, true])
    }

    // MARK: - Names and formats

    @Test func theFormatComesFromTheName() {
        #expect(AtomModel.format(forFileName: "model.xyz") == .extendedXYZ)
        #expect(AtomModel.format(forFileName: "/tmp/run.EXTXYZ") == .extendedXYZ)
        #expect(AtomModel.format(forFileName: "POSCAR") == .poscar)
        #expect(AtomModel.format(forFileName: "CONTCAR_relaxed") == .poscar)
        #expect(AtomModel.format(forFileName: "si.vasp") == .poscar)
        #expect(AtomModel.format(forFileName: "si.cif") == nil)
        #expect(throws: AtomModelError.unknownFormat(fileName: "si.cif")) {
            try AtomModel.read(bytes: Data("1\n\nSi 0 0 0\n".utf8), fileName: "si.cif")
        }
    }

    // MARK: - Refusals by name

    static func refusal(_ text: String, _ format: AtomModel.Format) -> String? {
        do {
            _ = try AtomModel.read(text, format: format)
            return nil
        } catch let error as AtomModelError {
            #expect(error.errorDescription?.isEmpty == false)
            return error.name
        } catch {
            return "\(error)"
        }
    }

    static let poscarHead = "Si\n1.0\n5.43 0 0\n0 5.43 0\n0 0 5.43\n"

    @Test func poscarRefusals() {
        // VASP 4: the counts follow the cell directly.
        #expect(Self.refusal(Self.poscarHead + "2\nDirect\n0 0 0\n0.25 0.25 0.25\n", .poscar) == "noElementLine")
        #expect(Self.refusal(Self.poscarHead + "Si\n1 1\nDirect\n0 0 0\n0.25 0.25 0.25\n", .poscar)
                == "countsWithoutElements")
        #expect(Self.refusal(Self.poscarHead + "Xq\n1\nDirect\n0 0 0\n", .poscar) == "unknownElement")
        #expect(Self.refusal(Self.poscarHead + "Si\n2\nDirect\n0 0 0\n", .poscar) == "missingAtoms")
        #expect(Self.refusal(Self.poscarHead + "Si\n1\nDirect\nnan 0 0\n", .poscar) == "badNumber")
        #expect(Self.refusal(Self.poscarHead + "Si\n1\nDirect\n0 0\n", .poscar) == "malformedLine")
        #expect(Self.refusal(Self.poscarHead + "Si\n1\n\n0 0 0\n", .poscar) == "malformedLine")
        #expect(Self.refusal(Self.poscarHead + "Si\n0\nDirect\n", .poscar) == "noAtoms")
        #expect(Self.refusal("Si\n1.0 2.0\n5.43 0 0\n0 5.43 0\n0 0 5.43\nSi\n1\nDirect\n0 0 0\n", .poscar) == "badScale")
        #expect(Self.refusal("Si\n1 -1 1\n5.43 0 0\n0 5.43 0\n0 0 5.43\nSi\n1\nDirect\n0 0 0\n", .poscar) == "badScale")
        #expect(Self.refusal("Si\n1.0\n5.43 0 0\n5.43 0 0\n0 0 5.43\nSi\n1\nDirect\n0 0 0\n", .poscar)
                == "degenerateCell")
        #expect(Self.refusal("", .poscar) == "noAtoms")
    }

    @Test func poscarDetails() throws {
        // A `!` comment on the counts line, Cartesian by its first letter K, the comment kept.
        let model = try AtomModel.read(Self.poscarHead + "Si O\n1 1 ! two atoms\nkartesian\n0 0 0\n1 1 1\n", format: .poscar)
        #expect(model.atomicNumbers == [14, 8])
        #expect(model.positions[1] == SIMD3<Double>(1, 1, 1))
        #expect(model.comment == "Si")
        // More element symbols than counts: the extra symbols are ignored, as in ASE.
        let extra = try AtomModel.read(Self.poscarHead + "Si O\n1\nDirect\n0.5 0.5 0.5\n", format: .poscar)
        #expect(extra.atomicNumbers == [14])
        #expect(extra.positions[0] == SIMD3<Double>(2.715, 2.715, 2.715))
        // Windows line breaks and old Mac ones read the same.
        let crlf = try AtomModel.read((Self.poscarHead + "Si\n1\nDirect\n0 0 0\n").replacingOccurrences(of: "\n", with: "\r\n"),
                                      format: .poscar)
        let cr = try AtomModel.read((Self.poscarHead + "Si\n1\nDirect\n0 0 0\n").replacingOccurrences(of: "\n", with: "\r"),
                                    format: .poscar)
        #expect(crlf == cr)
    }

    @Test func xyzRefusals() {
        #expect(Self.refusal("", .extendedXYZ) == "noAtoms")
        #expect(Self.refusal("two\n\nSi 0 0 0\n", .extendedXYZ) == "malformedLine")
        #expect(Self.refusal("0\n\n", .extendedXYZ) == "noAtoms")
        #expect(Self.refusal("2\n\nSi 0 0 0\n", .extendedXYZ) == "missingAtoms")
        #expect(Self.refusal("1\n\nQq 0 0 0\n", .extendedXYZ) == "unknownElement")
        #expect(Self.refusal("1\n\nX 0 0 0\n", .extendedXYZ) == "unknownElement")
        #expect(Self.refusal("1\n\nSi 0 0 inf\n", .extendedXYZ) == "badNumber")
        #expect(Self.refusal("1\n\nSi 0 0 1e999\n", .extendedXYZ) == "badNumber")
        #expect(Self.refusal("1\n\nSi 0 0\n", .extendedXYZ) == "malformedLine")
        #expect(Self.refusal("1\nLattice=\"1 0 0 0 1 0 0 0\"\nSi 0 0 0\n", .extendedXYZ) == "badLattice")
        #expect(Self.refusal("1\npbc=\"T T\"\nSi 0 0 0\n", .extendedXYZ) == "badPeriodicity")
        #expect(Self.refusal("1\nProperties=species:S:1\nSi\n", .extendedXYZ) == "badProperties")
        #expect(Self.refusal("1\nProperties=species:Q:1:pos:R:3\nSi 0 0 0\n", .extendedXYZ) == "badProperties")
        #expect(Self.refusal("1\nProperties=pos:R:3\n0 0 0\n", .extendedXYZ) == "badProperties")
        #expect(Self.refusal("1\nProperties=species:S:1:pos:R:3:tag:I:1\nSi 0 0 0 1.5\n", .extendedXYZ) == "badNumber")
        #expect(Self.refusal("1\n=x\nSi 0 0 0\n", .extendedXYZ) == "badCommentLine")
        #expect(Self.refusal("1\n\nSi 0 0 0\nVEC1 1 0 0\nVEC2 0 1 0\nVEC3 0 0 1\nVEC4 1 1 1\n", .extendedXYZ)
                == "malformedLine")
        #expect(Self.refusal("1\n\0\nSi 0 0 0\n", .extendedXYZ) == "notText")
    }

    @Test func xyzDetails() throws {
        // A Z column wins over species; `pbc` as numbers; `key = value` with spaces; a comment kept as written.
        let z = try AtomModel.read("1\nProperties=species:S:1:pos:R:3:Z:I:1 pbc = \"1 0 1\" note='a b'\nC 1 2 3 14\n",
                                   format: .extendedXYZ)
        #expect(z.atomicNumbers == [14])
        #expect(z.periodic == [true, false, true])
        #expect(z.positions == [SIMD3<Double>(1, 2, 3)])
        // A Lattice without pbc is periodic along all three; Lattice values may be separated by commas.
        let lattice = try AtomModel.read("1\nLattice=\"2,0,0 0,3,0 0,0,4\"\nAl 0 0 0\n", format: .extendedXYZ)
        #expect(lattice.cell == [SIMD3<Double>(2, 0, 0), SIMD3<Double>(0, 3, 0), SIMD3<Double>(0, 0, 4)])
        #expect(lattice.periodic == [true, true, true])
        // Extra columns beyond Properties are ignored; extra spaces too.
        let extra = try AtomModel.read("1\n\n  Fe   1.5  2.5   3.5   9 9 9\n", format: .extendedXYZ)
        #expect(extra.atomicNumbers == [26])
        #expect(extra.positions == [SIMD3<Double>(1.5, 2.5, 3.5)])
        // A blank header line ends the file: what follows is not read.
        let stop = try AtomModel.read("1\n\nFe 0 0 0\n\n1\n\nCu 0 0 0\n", format: .extendedXYZ)
        #expect(stop.frameCount == 1 && stop.atomicNumbers == [26])
    }

    @Test func theCommentLineIsSplitAsASESplitsIt() {
        #expect(AtomModel.keyValues("a=1 b=\"x y\" c") == ["a": "1", "b": "x y", "c": "T"])
        #expect(AtomModel.keyValues("key = value") == ["key": "value"])
        #expect(AtomModel.keyValues("k={1 2} l=[3,4] m='5'") == ["k": "1 2", "l": "3,4", "m": "5"])
        #expect(AtomModel.keyValues("a=b=c") == ["a": "b=c"])
        #expect(AtomModel.keyValues("a=x\\ y") == ["a": "x y"])
        #expect(AtomModel.keyValues("=a") == nil)
        #expect(AtomModel.keyValues("Cubic bulk silicon") == ["Cubic": "T", "bulk": "T", "silicon": "T"])
        #expect(AtomModel.periodicity("True false T") == [true, false, true])
        #expect(AtomModel.periodicity("F") == [false, false, false])
        #expect(AtomModel.periodicity("yes") == nil)
    }

    // MARK: - Never crashes

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

    static let hostile = ["nan", "inf", "-inf", "1e999", "1e-320", "0", "-1", "?", ".", "99999999999999999999", "Lattice=",
                          "pbc=", "Properties=species:S:1:pos:R:3:x:I:9999", "VEC1", "VEC9 1 2 3", "\"", "'", "{", "[", "\\",
                          "Selective", "Cartesian", "Direct", "X", "Na_pv/1", "!", "=", "1000000000", "-5"]

    static func mutate(_ text: String, using rng: inout Generator) -> String {
        var lines = text.components(separatedBy: "\n")
        for _ in 0..<Int.random(in: 1...4, using: &rng) {
            if lines.isEmpty { lines = [""] }
            let i = Int.random(in: 0..<lines.count, using: &rng)
            switch Int.random(in: 0..<7, using: &rng) {
            case 0, 1:
                var tokens = lines[i].split(separator: " ").map(String.init)
                if tokens.isEmpty { tokens = [""] }
                tokens[Int.random(in: 0..<tokens.count, using: &rng)] = hostile.randomElement(using: &rng)!
                lines[i] = tokens.joined(separator: " ")
            case 2: lines.remove(at: i)
            case 3: lines.insert(lines[i], at: i)
            case 4: lines.swapAt(i, Int.random(in: 0..<lines.count, using: &rng))
            case 5: lines.insert(hostile.randomElement(using: &rng)!, at: i)
            default:
                lines[i] = String(lines[i].prefix(Int.random(in: 0...lines[i].count, using: &rng)))
                lines.removeSubrange((i + 1)...)
            }
        }
        return lines.joined(separator: "\n")
    }

    @Test func mutatedCorpusFilesNeverCrashTheReader() throws {
        var rng = Generator(state: 0x169)
        let names = try FileManager.default.contentsOfDirectory(atPath: Self.corpusDirectory.path).sorted()
        #expect(names.count >= 25)
        var reads = 0
        for name in names {
            let text = String(decoding: try Self.bytes(name), as: UTF8.self)
            for format in AtomModel.Format.allCases {
                for _ in 0..<40 {
                    if let model = try? AtomModel.read(Self.mutate(text, using: &rng), format: format) {
                        #expect(model.atomicNumbers.count == model.positions.count)
                        #expect(model.cell.count == 3 && model.periodic.count == 3)
                        #expect(model.positions.allSatisfy { $0.x.isFinite && $0.y.isFinite && $0.z.isFinite })
                        _ = model.countsByAtomicNumber
                    }
                    reads += 1
                }
            }
        }
        #expect(reads >= 25 * 2 * 40)
    }

    @Test func randomBytesNeverCrashTheReader() {
        var rng = Generator(state: 0xB17E5_169)
        let alphabet = Array("0123456789 .-+eE\n\r\tSiAlOVECLatice=\"'{}[]\\!:RISpbcTF_/".utf8)
        for _ in 0..<2000 {
            let n = Int.random(in: 0..<200, using: &rng)
            let bytes = (0..<n).map { _ in alphabet.randomElement(using: &rng)! }
            for format in AtomModel.Format.allCases {
                _ = try? AtomModel.read(bytes: Data(bytes), format: format)
            }
        }
    }

    /// A header that announces far more atoms than the file holds is refused without reserving memory for them.
    @Test func aHugeAnnouncedCountIsRefusedCheaply() {
        #expect(Self.refusal("999999999999999\n\nSi 0 0 0\n", .extendedXYZ) == "missingAtoms")
        #expect(Self.refusal(Self.poscarHead + "Si\n999999999\nDirect\n0 0 0\n", .poscar) == "missingAtoms")
    }

    // MARK: - Speed

    @Test func aMillionAtomsReadInUnderTwoSeconds() throws {
        let n = 1_000_000
        let rows = ["Si 1.35772500 2.71545000 0.00000000\n", "O -0.76323900 4.07317500 12.12345678\n",
                    "Al 10.86180000 0.00000000 5.43090000\n", "Sr 3.90500000 3.90500000 3.90500000\n"]
        var bytes: [UInt8] = []
        bytes.reserveCapacity(n * 40 + 200)
        bytes += Array("\(n)\nLattice=\"100 0 0 0 100 0 0 0 100\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n".utf8)
        let encoded = rows.map { Array($0.utf8) }
        for i in 0..<n { bytes += encoded[i % 4] }
        let data = Data(bytes)
        let clock = ContinuousClock()
        var model: AtomModel?
        let elapsed = try clock.measure { model = try AtomModel.read(bytes: data, format: .extendedXYZ) }
        #expect(model?.count == n)
        #expect(model?.countsByAtomicNumber == [14: n / 4, 8: n / 4, 13: n / 4, 38: n / 4])
        #expect(model?.positions[1] == SIMD3<Double>(-0.763239, 4.073175, 12.12345678))
        // A wall-clock bound measures the shared CI runner's load, so it is checked only in the opt-in benchmark run.
        if ProcessInfo.processInfo.environment["DIFFRACTOKIT_BENCHMARK"] != nil {
            #expect(elapsed < .seconds(2), "\(elapsed)")
        }
    }
}
