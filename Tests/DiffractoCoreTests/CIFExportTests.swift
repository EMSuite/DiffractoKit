//
//  CIFExportTests.swift
//  A crystal structure written as CIF (`CIFExport`, DK1). Pin C2: every structure of the corpus (P2's 19, six with an
//  edited lattice, three Materials Project imports) written and read back by `CIFImport` has the structure's cell and its
//  expanded sites, as a multiset, equal to the bit. The test reads back itself, so a writer whose own read-back check were
//  wrong would still be caught. Symmetry is kept for the files whose loop built their sites; the two planted structures
//  whose operations do not expand their sites exactly fall back to P1 and say so. Pin C1 (pymatgen reads the same text)
//  is `CIFExportPymatgenTests` in the parity target.
//  Red before: `CIFExport` did not exist.
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

    @Test func everyStructureOfTheCorpusReadsBackToTheBit() throws {
        let corpus = try CIFExportCorpus.structures()
        #expect(corpus.count == 19 + 6 + 3)
        var p1: [String] = []
        for (name, structure) in corpus {
            let outcome = try CIFExport.text(of: structure, options: .init(comments: ["corpus \(name)"]))
            let read = try CIFImport.read(outcome.text).structure
            expectSameStructure(read, structure, name)
            if case .p1 = outcome.symmetry { p1.append(name) }
        }
        // The files whose own loop built their sites keep their symmetry (the same base sites, the same operations), and so
        // does silicon's P1 file, whose detected operations map its dyadic sites exactly.
        for name in ["rocksalt_Fm-3m_192ops", "pyrite_Pa-3", "beta_double_prime", "rocksalt_Fm-3m_a5.7", "silicon_Fd-3m",
                     "silicon_Fd-3m_a5.5"] {
            #expect(!p1.contains(name), "\(name) fell back to P1")
        }
        print("CIFExport corpus written in P1: \(p1)")
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

    /// A structure built with operations that do not reproduce its sites to the bit: the orbit is found within 1e-6, the
    /// read-back expands the representative exactly, the bits differ, and the file is P1 with every site and the reason.
    @Test func aSiteOffItsExactImageFallsBackToP1() throws {
        let fcc = try CIFImport.read(try Fixtures.text("gold_centring_ops.cif")).structure
        var sites = fcc.sites
        let moved = sites[1].fractional + SIMD3(1e-9, 0, 0)
        sites[1] = AtomSite(z: sites[1].z, fractional: moved, occupancy: sites[1].occupancy, label: sites[1].label,
                            isotropicDisplacement: sites[1].isotropicDisplacement)
        let planted = CrystalStructure(cell: fcc.cell, sites: sites, operations: fcc.operations, spaceGroupNumber: 225, name: "gold")
        let outcome = try CIFExport.text(of: planted)
        guard case let .p1(reason) = outcome.symmetry else { Issue.record("kept: \(outcome.symmetry)"); return }
        #expect(reason?.contains("to the bit") == true, "\(String(describing: reason))")
        #expect(outcome.text.contains("# Written in P1 with every site"))
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
