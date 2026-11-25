//
//  CIFExportTests.swift
//  A crystal structure written as CIF (`CIFExport`, DK1). Pin C2: every structure of the corpus (P2's 19, six with an
//  edited lattice, three Materials Project imports) written and read back by `CIFImport` has the structure's cell and its
//  expanded sites, as a multiset, equal to the bit, or, where the file says it kept the symmetry within pymatgen's reading
//  tolerance, within 1e-4 modulo 1 and with the structure's content fingerprint. The test reads back itself, so a writer
//  whose own read-back check were wrong would still be caught. Every structure of the corpus keeps its symmetry; the
//  seven hexagonal ones with sites at 1/3 and 2/3 do so within the tolerance (#220). The planted structures whose
//  operations do not expand their sites fall back to P1 and say so. Pin C1 (pymatgen reads the same text) is
//  `CIFExportPymatgenTests` in the parity target.
//  Red before: `CIFExport` did not exist; the seven hexagonal structures were written in P1 before #220.
//

import DiffractoKitTestSupport
import Foundation
import Testing
@testable import DiffractoCore

struct CIFExportTests {
    /// A site as compared: element, coordinates, occupancy, displacement, each with `==` on its double.
    private struct Key: Comparable {
        let z: Int, x: [Double], occupancy: Double, u: Double
        init(_ s: AtomSite) {
            z = s.z; x = [s.fractional.x, s.fractional.y, s.fractional.z]; occupancy = s.occupancy; u = s.isotropicDisplacement ?? -1
        }
        static func < (l: Key, r: Key) -> Bool {
            if l.z != r.z { return l.z < r.z }
            if l.x != r.x { return l.x.lexicographicallyPrecedes(r.x) }
            if l.occupancy != r.occupancy { return l.occupancy < r.occupancy }
            return l.u < r.u
        }
    }

    private func expectSameStructure(_ read: CrystalStructure, _ original: CrystalStructure, _ name: String) {
        let a = original.cell, b = read.cell
        #expect([a.a, a.b, a.c, a.alphaDeg, a.betaDeg, a.gammaDeg] == [b.a, b.b, b.c, b.alphaDeg, b.betaDeg, b.gammaDeg], "\(name) cell")
        #expect(read.sites.map(Key.init).sorted() == original.sites.map(Key.init).sorted(), "\(name) sites")
    }

    /// The same cell to the bit, and the same sites as a multiset within `tolerance` modulo 1 (element, occupancy and
    /// displacement equal).
    private func expectSameStructure(_ read: CrystalStructure, _ original: CrystalStructure, _ name: String, within tolerance: Double) {
        let a = original.cell, b = read.cell
        #expect([a.a, a.b, a.c, a.alphaDeg, a.betaDeg, a.gammaDeg] == [b.a, b.b, b.c, b.alphaDeg, b.betaDeg, b.gammaDeg], "\(name) cell")
        #expect(read.sites.count == original.sites.count, "\(name) site count")
        var unmatched = Array(original.sites.indices)
        for site in read.sites {
            let found = unmatched.firstIndex { i in
                let mine = original.sites[i]
                return mine.z == site.z && mine.occupancy == site.occupancy && mine.isotropicDisplacement == site.isotropicDisplacement
                    && (0..<3).allSatisfy { fractionalAxisDelta(mine.fractional[$0], site.fractional[$0]) <= tolerance }
            }
            if let found { unmatched.remove(at: found) } else { Issue.record("\(name): \(site.fractional) is not a site") }
        }
    }

    static let tolerantNote = "# Symmetry kept within pymatgen's reading tolerance (1e-4)"

    /// The hexagonal structures with sites at 1/3 and 2/3, written in P1 before #220, and their IT numbers.
    static let hexagonal: [(name: String, number: Int)] = [
        ("magnesium_P63mmc", 194), ("ws2_2H", 194), ("cdi2_P-3m1", 164), ("ideal_hcp_P1", 194), ("magnesium_P63mmc_c5.3", 194),
        ("cdi2_P-3m1_a4.3", 164), ("input-magnesiumJSON.json", 194),
    ]

    @Test func everyStructureOfTheCorpusReadsBackAndKeepsItsSymmetry() throws {
        let corpus = try CIFExportCorpus.structures()
        #expect(corpus.count == 19 + 6 + 3)
        var p1: [String] = []
        var tolerant: [String] = []
        for (name, structure) in corpus {
            let outcome = try CIFExport.text(of: structure, options: .init(comments: ["corpus \(name)"]))
            let read = try CIFImport.read(outcome.text).structure
            if case .p1 = outcome.symmetry { p1.append(name) }
            if outcome.text.contains(Self.tolerantNote) {
                tolerant.append(name)
                expectSameStructure(read, structure, name, within: 1e-4)
            } else {
                expectSameStructure(read, structure, name)
            }
        }
        #expect(p1.isEmpty, "written in P1: \(p1)")
        // Only the hexagonal structures need the tolerance; every other text reads back to the bit.
        #expect(Set(tolerant) == Set(Self.hexagonal.map(\.name)), "kept within the tolerance: \(tolerant)")
    }

    /// #220: the seven keep their space group, and read back to the same cell, sites within pymatgen's 1e-4 and the same
    /// content fingerprint, so a phase saved as CIF and added again is the same phase.
    @Test func hexagonalStructuresKeepTheirSymmetryAndTheirFingerprint() throws {
        let corpus = try CIFExportCorpus.structures()
        for (name, number) in Self.hexagonal {
            guard let structure = corpus.first(where: { $0.name == name })?.structure else { Issue.record("\(name) missing"); continue }
            let outcome = try CIFExport.text(of: structure)
            guard case let .kept(operations, sites, expanded) = outcome.symmetry else {
                Issue.record("\(name) in P1: \(outcome.symmetry)"); continue
            }
            #expect(operations > 1 && sites < expanded && expanded == structure.sites.count, "\(name): \(outcome.symmetry)")
            #expect(outcome.text.contains("_space_group_IT_number \(number)\n"), "\(name)")
            #expect(outcome.text.contains(Self.tolerantNote), "\(name)")
            #expect(!outcome.text.contains("Written in P1"), "\(name)")
            let read = try CIFImport.read(outcome.text).structure
            expectSameStructure(read, structure, name, within: 1e-4)
            #expect(read.fingerprint == structure.fingerprint, "\(name): fingerprint")
            #expect(read.operations.count == structure.operations.count, "\(name): operations")
        }
    }

    /// pymatgen's `frac_tolerance` rule: a coordinate whose ratio to 1/3 or 2/3 is within 1e-4 of 1 is that fraction.
    @Test func coordinatesNearAThirdAreWrittenAsPymatgenReadsThem() {
        #expect(CIFExport.exactThird(0.333333) == 1.0 / 3)
        #expect(CIFExport.exactThird(0.33334) == 1.0 / 3)
        #expect(CIFExport.exactThird(0.6666667) == 2.0 / 3)
        #expect(CIFExport.exactThird(0.3332) == 0.3332)
        #expect(CIFExport.exactThird(0.25) == 0.25)
        #expect(CIFExport.exactThird(0) == 0)
        #expect(CIFExport.number(CIFExport.exactThird(0.333333)) == "0.3333333333333333")
    }

    @Test func aFileLoopIsWrittenAsItsAsymmetricUnit() throws {
        let rocksalt = try CIFImport.read(try Fixtures.text("rocksalt_Fm-3m_192ops.cif")).structure
        let outcome = try CIFExport.text(of: rocksalt)
        guard case let .kept(operations, sites, expanded) = outcome.symmetry else {
            Issue.record("rock salt in P1: \(outcome.symmetry)"); return
        }
        #expect(sites == 2)
        #expect(expanded == rocksalt.sites.count)
        #expect(operations == rocksalt.fileOperations?.count)
        #expect(outcome.text.contains("_space_group_IT_number 225"))
        #expect(outcome.text.contains("loop_\n_space_group_symop_operation_xyz\n'x,y,z'\n"))
        #expect(outcome.text.hasSuffix("\n"))
    }

    /// Gold with one site moved off its exact image by `offset` along x.
    private func plantedGold(offset: Double) throws -> CrystalStructure {
        let fcc = try CIFImport.read(try Fixtures.text("gold_centring_ops.cif")).structure
        var sites = fcc.sites
        let moved = sites[1].fractional + SIMD3(offset, 0, 0)
        sites[1] = AtomSite(z: sites[1].z, fractional: moved, occupancy: sites[1].occupancy, label: sites[1].label,
                            isotropicDisplacement: sites[1].isotropicDisplacement)
        return CrystalStructure(cell: fcc.cell, sites: sites, operations: fcc.operations, spaceGroupNumber: 225, name: "gold")
    }

    /// A site 1e-9 off its exact image: not to the bit, but within pymatgen's reading tolerance, so the symmetry is kept and
    /// the file says how; the read-back is the structure within 1e-4 and has its fingerprint.
    @Test func aSiteJustOffItsExactImageKeepsTheSymmetryWithinTheTolerance() throws {
        let planted = try plantedGold(offset: 1e-9)
        let outcome = try CIFExport.text(of: planted)
        guard case let .kept(_, sites, expanded) = outcome.symmetry else { Issue.record("P1: \(outcome.symmetry)"); return }
        #expect(sites == 1 && expanded == 4)
        #expect(outcome.text.contains(Self.tolerantNote))
        let read = try CIFImport.read(outcome.text).structure
        expectSameStructure(read, planted, "planted gold", within: 1e-4)
        #expect(read.fingerprint == planted.fingerprint)
    }

    /// A site 1e-3 off its exact image, beyond pymatgen's tolerance: its own orbit, whose images the read-back cannot
    /// place, and the file is P1 with every site and the reason.
    @Test func aSiteOffItsImageBeyondTheToleranceFallsBackToP1() throws {
        let planted = try plantedGold(offset: 1e-3)
        let outcome = try CIFExport.text(of: planted)
        guard case let .p1(reason) = outcome.symmetry else { Issue.record("kept: \(outcome.symmetry)"); return }
        #expect(reason != nil)
        #expect(outcome.text.contains("# Written in P1 with every site"))
        #expect(!outcome.text.contains(Self.tolerantNote))
        #expect(outcome.text.contains("_space_group_IT_number 1\n"))
        expectSameStructure(try CIFImport.read(outcome.text).structure, planted, "planted gold")
    }

    /// Operations given that the sites do not have (an inversion for a lone atom off the centre): the representative
    /// expands to two sites, the structure has one, and the file is P1.
    @Test func operationsTheSitesDoNotHaveFallBackToP1() throws {
        let cell = UnitCell(a: 4, b: 5, c: 6, alphaDeg: 90, betaDeg: 90, gammaDeg: 90)
        let lone = [AtomSite(z: 14, fractional: SIMD3(0.1, 0.2, 0.3))]
        let inversion = try SymmetryOperation(xyz: "-x,-y,-z")
        let planted = CrystalStructure(cell: cell, sites: lone, operations: [.identity, inversion])
        let outcome = try CIFExport.text(of: planted)
        guard case let .p1(reason) = outcome.symmetry else { Issue.record("kept: \(outcome.symmetry)"); return }
        #expect(reason?.contains("2 sites") == true, "\(String(describing: reason))")
        expectSameStructure(try CIFImport.read(outcome.text).structure, planted, "lone silicon")
    }

    @Test func aStructureOfTheIdentityAloneIsP1WithoutAReason() throws {
        let cell = UnitCell(a: 3, b: 3, c: 3, alphaDeg: 90, betaDeg: 90, gammaDeg: 90)
        let structure = CrystalStructure(cell: cell, sites: [AtomSite(z: 79, fractional: SIMD3(0.1, 0.2, 0.3))], operations: [.identity])
        let outcome = try CIFExport.text(of: structure)
        #expect(outcome.symmetry == .p1(reason: nil))
        #expect(!outcome.text.contains("Written in P1"))
    }

    /// Numbers are Swift's shortest form (0.1 is "0.1", never 0.10000000000000001); a translation is the fraction that
    /// divides to it exactly; labels, block name, comments, source and creator are written as given or made legal.
    @Test func theTextCarriesItsNumbersLabelsAndNotes() throws {
        let cell = UnitCell(a: 3.1, b: 3.1, c: 5.2, alphaDeg: 90, betaDeg: 90, gammaDeg: 120)
        let sites = [AtomSite(z: 12, fractional: SIMD3(0.1, 0.25, 0.7), occupancy: 0.5, label: nil, isotropicDisplacement: 0.012),
                     AtomSite(z: 12, fractional: SIMD3(0.4, 0.6, 0.9), occupancy: 1, label: "my site")]
        let structure = CrystalStructure(cell: cell, sites: sites, operations: [.identity], name: "Mg test")
        let options = CIFExport.Options(blockName: "Mg test #1", comments: ["lattice edited", "two\nlines"],
                                        sourceLines: ["mp-153, 2026-11-24, CC BY 4.0"], creator: "Electron Microscopy Studio 2.4.0")
        let text = try CIFExport.text(of: structure, options: options).text
        #expect(text.hasPrefix("# lattice edited\n# two\n# lines\n# Source:\n#   mp-153, 2026-11-24, CC BY 4.0\ndata_Mg_test__1\n"))
        #expect(text.contains("_audit_creation_method 'Electron Microscopy Studio 2.4.0'"))
        #expect(text.contains("_chemical_name_common 'Mg test'"))
        #expect(text.contains("_cell_length_a 3.1\n"))
        #expect(text.contains("_cell_angle_gamma 120.0\n"))
        #expect(text.contains("Mg1 Mg 0.1 0.25 0.7 0.5 0.012\n"))
        #expect(text.contains("Mg2 Mg 0.4 0.6 0.9 1.0 ?\n"), "a label with a space is replaced, a missing U is ?")
        let third = SymmetryOperation(rotation: .identity, translation: SIMD3(1.0 / 3, 2.0 / 3, 0.5))
        #expect(CIFExport.xyz(third) == "x+1/3,y+2/3,z+1/2")
        let odd = SymmetryOperation(rotation: .identity, translation: SIMD3(0.123, 0, -0.25))
        #expect(CIFExport.xyz(odd) == "x+0.123,y,z-1/4")
        #expect(try SymmetryOperation(xyz: CIFExport.xyz(third)).translation == third.translation)
        // Positional, never an exponent: the reader takes a coordinate's precision from its decimals.
        for value in [0.00005, -1.25e-7, 1e-17, 123456.0, 1.5e20, 0.1, 0] {
            let text = CIFExport.number(value)
            #expect(!text.contains("e"), "\(text)")
            #expect(Double(text) == value, "\(text)")
        }
        #expect(CIFExport.number(0.00005) == "0.00005")
        #expect(CIFExport.isBareToken("O#1") == false)
    }

    @Test func valuesThatCannotBeWrittenAreRefusedByName() {
        let cell = UnitCell(a: 3, b: 3, c: 3, alphaDeg: 90, betaDeg: 90, gammaDeg: 90)
        let unknown = CrystalStructure(cell: cell, sites: [AtomSite(z: 200, fractional: SIMD3(0, 0, 0))], operations: [.identity])
        #expect(throws: CIFExportError.unknownElement(z: 200)) { try CIFExport.text(of: unknown) }
        let nan = CrystalStructure(cell: cell, sites: [AtomSite(z: 13, fractional: SIMD3(.nan, 0, 0))], operations: [.identity])
        #expect(throws: CIFExportError.nonFiniteValue("A value of site 1")) { try CIFExport.text(of: nan) }
    }
}
