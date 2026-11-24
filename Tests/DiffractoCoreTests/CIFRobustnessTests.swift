//
//  CIFRobustnessTests.swift
//  Role: the CIFs people actually have (DiffractoKit #117, plan row K4): unknown and inapplicable values, deuterium and
//        tritium, dummy atoms, Latin-1 text and a byte-order mark, DDLm names, a publication block before the structure,
//        a tag without its value, displacement parameters; mmCIF and mCIF refused by name; a number that is not finite
//        refused with a typed error; and the public paths that trapped on such a number. Each refusal or read here
//        failed (or trapped, which shows as a crashed test process) on 0.18.1.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore

struct CIFRobustnessTests {

    static let edgeCases = FixtureDirs.fixtures(of: "DiffractoCoreTests").appendingPathComponent("cif-edge-cases")

    static func bytes(_ name: String) throws -> Data {
        try Data(contentsOf: edgeCases.appendingPathComponent(name + ".cif"))
    }

    static func read(_ name: String) throws -> CIFDocument {
        try CrystalLoading.load(bytes: try bytes(name), fileName: name + ".cif").document
    }

    private func refusal(_ cif: String) -> CIFImportError? {
        do { _ = try CIFImport.read(cif, fileBaseName: "t"); Issue.record("admitted"); return nil }
        catch { return error }
    }

    private static let gold = CIFImportTests.goldP1

    // MARK: Numbers that are not finite

    /// `Double("nan")`, `Double("inf")` and `Double("1e999")` all parse; on 0.18.1 the NaN reached the symmetry detection's
    /// cell list and trapped in `Int(NaN)`. Each is now refused, naming the tag and the text.
    @Test(arguments: ["nan", "NaN", "inf", "-inf", "infinity", "1e999", "-1e999"])
    func aNonFiniteCoordinateIsRefused(_ text: String) {
        let cif = Self.gold.replacingOccurrences(of: "Au2 Au 0.5 0.5 0", with: "Au2 Au \(text) 0.5 0")
        #expect(refusal(cif) == .nonFiniteNumber(tag: "_atom_site_fract_x (Au2)", value: text))
    }

    @Test func aNonFiniteCellOrOccupancyIsRefused() {
        #expect(refusal(Self.gold.replacingOccurrences(of: "_cell_length_b 4.0782", with: "_cell_length_b inf"))
                == .nonFiniteNumber(tag: "_cell_length_b", value: "inf"))
        #expect(refusal(Self.gold.replacingOccurrences(of: "_cell_angle_gamma 90", with: "_cell_angle_gamma nan(2)"))
                == .nonFiniteNumber(tag: "_cell_angle_gamma", value: "nan(2)"))
        let occupancy = Self.gold.replacingOccurrences(of: "_atom_site_fract_z\n", with: "_atom_site_fract_z\n_atom_site_occupancy\n")
            .replacingOccurrences(of: "Au1 Au 0 0 0", with: "Au1 Au 0 0 0 nan")
            .replacingOccurrences(of: "Au2 Au 0.5 0.5 0", with: "Au2 Au 0.5 0.5 0 1")
            .replacingOccurrences(of: "Au3 Au 0.5 0 0.5", with: "Au3 Au 0.5 0 0.5 1")
            .replacingOccurrences(of: "Au4 Au 0 0.5 0.5", with: "Au4 Au 0 0.5 0.5 1")
        #expect(refusal(occupancy) == .nonFiniteNumber(tag: "_atom_site_occupancy (Au1)", value: "nan"))
    }

    @Test func theNaNFixtureIsRefusedThroughTheLoader() throws {
        #expect(throws: CIFImportError.nonFiniteNumber(tag: "_atom_site_fract_x (Au1)", value: "nan")) {
            try Self.read("nan_coordinate")
        }
    }

    /// The public routes past the importer: a site that is not finite gives the identity and a refusal, never a trap.
    @Test func theDetectionRefusesASiteThatIsNotFinite() {
        let sites = [AtomSite(z: 79, fractional: SIMD3(0, 0, 0)), AtomSite(z: 79, fractional: SIMD3(.nan, 0.5, 0), label: "Au2")]
        let structure = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0782), sites: sites)
        #expect(structure.operations == [.identity])
        #expect(structure.symmetryRefusal?.contains("Au2") == true)
        #expect(!structure.isUsable)
        #expect(structure.validationIssues.map(\.code) == ["invalid_fractional_coordinate"])
        let flat = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0782), sites: sites, conventional: false)
        #expect(flat.operations == [.identity])
        #expect(!flat.isUsable)
        let infinite = [AtomSite(z: 79, fractional: SIMD3(0, .infinity, 0))]
        let detection = SymmetryDetection.detect(sites: infinite, cell: .cubic(a: 4.0782))
        #expect(detection.operations == [.identity])
        #expect(detection.refusal != nil)
    }

    @Test func aMaterialsProjectSiteThatIsNotFiniteIsMalformed() {
        func document(abc: [Double], occupancy: Double = 1, a: Double = 4.05) -> MaterialsProjectDocument {
            MaterialsProjectDocument(
                materialID: "mp-134",
                structure: MaterialsProjectStructure(
                    lattice: MaterialsProjectLattice(a: a, b: 4.05, c: 4.05, alpha: 90, beta: 90, gamma: 90),
                    sites: [MaterialsProjectSite(label: "Al", species: [MaterialsProjectSpecies(element: "Al", occu: occupancy)],
                                                 abc: abc)]),
                symmetry: MaterialsProjectSymmetry(crystalSystem: "Cubic", symbol: "Pm-3m", number: 221))
        }
        for doc in [document(abc: [.nan, 0, 0]), document(abc: [0, .infinity, 0]), document(abc: [0, 0, 0], occupancy: .nan),
                    document(abc: [0, 0, 0], a: .nan)] {
            do {
                _ = try MaterialsProjectImport.structure(from: doc)
                Issue.record("admitted")
            } catch let caught as MaterialsProjectImport.Failure {
                guard case .malformed = caught else { Issue.record("\(caught)"); continue }
            } catch { Issue.record("\(error)") }
        }
    }

    /// An occupancy of 1e30 is no structure, but its formula must not trap in `Int(1e30)`.
    @Test func theFormulaOfAHugeOccupancyDoesNotTrap() {
        let s = CrystalStructure(cell: .cubic(a: 4), sites: [AtomSite(z: 13, fractional: .zero, occupancy: 1e30)],
                                 operations: [.identity])
        #expect(s.formula == "Al1e+30")
        let nan = CrystalStructure(cell: .cubic(a: 4), sites: [AtomSite(z: 13, fractional: .zero, occupancy: .nan)],
                                   operations: [.identity])
        #expect(nan.formula == "Alnan")
    }

    // MARK: Unknown and inapplicable values

    /// `?` and `.` occupancies are CIF's default 1; an unknown type symbol falls back to the label; unknown U is nil.
    @Test func unknownValuesTakeTheirDefaults() throws {
        let doc = try Self.read("unknown_values")
        let s = doc.structure
        #expect(s.sites.count == 4)
        #expect(s.formula == "Au4")
        #expect(s.sites.allSatisfy { $0.occupancy == 1 })
        #expect(doc.declaredSpaceGroupSymbol == nil)
        #expect(doc.declaredSpaceGroupNumber == nil)
        let byLabel = Dictionary(uniqueKeysWithValues: s.sites.map { ($0.label ?? "", $0) })
        #expect(byLabel["Au1"]?.isotropicDisplacement == nil)
        #expect(byLabel["Au2"]?.isotropicDisplacement == nil)
        #expect(byLabel["Au3"]?.isotropicDisplacement == 0.0081)
        #expect(byLabel["Au4"]?.isotropicDisplacement == 0.0081)
    }

    @Test func anUnknownCoordinateOrCellParameterIsRefusedByName() {
        #expect(refusal(Self.gold.replacingOccurrences(of: "Au2 Au 0.5 0.5 0", with: "Au2 Au 0.5 ? 0"))
                == .unknownValue(tag: "_atom_site_fract_y (Au2)"))
        #expect(refusal(Self.gold.replacingOccurrences(of: "_cell_length_c 4.0782", with: "_cell_length_c ."))
                == .unknownValue(tag: "_cell_length_c"))
    }

    // MARK: Isotopes, dummy atoms, displacement

    @Test func deuteriumAndTritiumAreHydrogenAndBIsReadAsU() throws {
        let s = try Self.read("isotopes_P1").structure
        #expect(s.formula == "H2 O1")
        #expect(s.sites.map(\.z).sorted() == [1, 1, 8])
        let u = try #require(s.sites.first { $0.z == 8 }?.isotropicDisplacement)
        #expect(abs(u - 1.2 / (8 * Double.pi * Double.pi)) < 1e-15)
        #expect(Element.normalizedSymbol(fromCIF: "D1") == "H")
        #expect(Element.normalizedSymbol(fromCIF: "T") == "H")
        #expect(Element.normalizedSymbol(fromCIF: "Dy3+") == "Dy")
        #expect(Element.normalizedSymbol(fromCIF: "Ta") == "Ta")
    }

    @Test func uWinsOverBWhenBothAreWritten() throws {
        let cif = Self.gold.replacingOccurrences(of: "_atom_site_fract_z\n", with: "_atom_site_fract_z\n_atom_site_B_iso_or_equiv\n_atom_site_U_iso_or_equiv\n")
            .replacingOccurrences(of: "Au1 Au 0 0 0", with: "Au1 Au 0 0 0 0.6 0.0076(3)")
            .replacingOccurrences(of: "Au2 Au 0.5 0.5 0", with: "Au2 Au 0.5 0.5 0 0.6 ?")
            .replacingOccurrences(of: "Au3 Au 0.5 0 0.5", with: "Au3 Au 0.5 0 0.5 ? ?")
            .replacingOccurrences(of: "Au4 Au 0 0.5 0.5", with: "Au4 Au 0 0.5 0.5 . 0.0076")
        let sites = try CIFImport.read(cif).structure.sites
        let u = Dictionary(uniqueKeysWithValues: sites.map { ($0.label ?? "", $0.isotropicDisplacement) })
        #expect(u["Au1"] == 0.0076)
        #expect(u["Au2"] == 0.6 / (8 * Double.pi * Double.pi))
        #expect(u["Au3"].map { $0 == nil } == true)
        #expect(u["Au4"] == 0.0076)
    }

    /// The displacement rides with every image of the site through the expansion.
    @Test func theDisplacementFollowsTheSymmetryExpansion() throws {
        let cif = """
        data_Al
        _cell_length_a 4.0495
        _cell_length_b 4.0495
        _cell_length_c 4.0495
        _cell_angle_alpha 90
        _cell_angle_beta 90
        _cell_angle_gamma 90
        loop_
        _symmetry_equiv_pos_as_xyz
        'x,y,z'
        'x,y+1/2,z+1/2'
        'x+1/2,y,z+1/2'
        'x+1/2,y+1/2,z'
        loop_
        _atom_site_label
        _atom_site_fract_x
        _atom_site_fract_y
        _atom_site_fract_z
        _atom_site_U_iso_or_equiv
        Al1 0 0 0 0.0101

        """
        let s = try CIFImport.read(cif).structure
        #expect(s.sites.count == 4)
        #expect(s.sites.allSatisfy { $0.isotropicDisplacement == 0.0101 })
    }

    /// A no-break space between tokens (0xA0 in a Latin-1 file) separates them; on 0.18.1 it hung the tokenizer.
    @Test func aNoBreakSpaceSeparatesTokens() throws {
        var bytes = [UInt8](Data(Self.gold.utf8))
        let line = [UInt8]("Au2 Au 0.5 0.5 0".utf8)
        let at = try #require((0...(bytes.count - line.count)).first { Array(bytes[$0..<($0 + line.count)]) == line })
        bytes[at + 3] = 0xA0                          // "Au2\u{A0}Au 0.5 0.5 0"
        #expect(try CrystalLoading.load(bytes: Data(bytes), fileName: "nbsp.cif").structure.formula == "Au4")
        #expect(try CIFImport.read(Self.gold.replacingOccurrences(of: "Au2 Au", with: "Au2\u{A0}Au")).structure.sites.count == 4)
    }

    /// A displacement that is no number never refuses the file: nothing depends on it.
    @Test func anUnreadableDisplacementIsUnknown() throws {
        let cif = Self.gold.replacingOccurrences(of: "_atom_site_fract_z\n", with: "_atom_site_fract_z\n_atom_site_U_iso_or_equiv\n")
            .replacingOccurrences(of: "Au1 Au 0 0 0", with: "Au1 Au 0 0 0 nan")
            .replacingOccurrences(of: "Au2 Au 0.5 0.5 0", with: "Au2 Au 0.5 0.5 0 abc")
            .replacingOccurrences(of: "Au3 Au 0.5 0 0.5", with: "Au3 Au 0.5 0 0.5 0.01")
            .replacingOccurrences(of: "Au4 Au 0 0.5 0.5", with: "Au4 Au 0 0.5 0.5 0.01")
        let sites = try CIFImport.read(cif).structure.sites
        #expect(sites.filter { $0.isotropicDisplacement == nil }.count == 2)
    }

    @Test func dummyAtomsAreLeftOutWithAWarning() throws {
        let doc = try Self.read("dummy_atoms")
        #expect(doc.structure.formula == "Au4")
        #expect(doc.warnings == [.dummyAtomsSkipped(labels: ["Q1", "Lp1", "X1"])])
        #expect(doc.warnings[0].description.contains("Q1, Lp1, X1"))
        // Only dummy atoms: nothing to scatter.
        let onlyDummy = Self.gold.replacingOccurrences(of: "Au", with: "Xx")
        #expect(refusal(onlyDummy) == .noAtomSites)
        // An unknown symbol that is no dummy is still refused by name.
        #expect(refusal(Self.gold.replacingOccurrences(of: "Au1 Au 0 0 0", with: "Zz1 Zz 0 0 0"))
                == .unknownElementSymbol("Zz", atomSiteLabel: "Zz1"))
    }

    // MARK: Encodings

    @Test func latin1AndAByteOrderMarkAreRead() throws {
        let latin1 = try Self.read("latin1_name")
        #expect(latin1.displayName == "Or natif (échantillon)")
        #expect(latin1.structure.formula == "Au4")
        let bom = try Self.read("utf8_bom")
        #expect(bom.dataBlockName == "bom")
        #expect(bom.structure.formula == "Au4")
        // The BOM as text, too (a caller that decoded the file itself).
        #expect(try CIFImport.read("\u{FEFF}" + Self.gold).dataBlockName == "gold")
        // The fingerprint stays the hash of the bytes as given.
        let raw = try Self.bytes("utf8_bom")
        #expect(try CrystalLoading.load(bytes: raw, fileName: "b.cif").sha256 == CrystalLoading.fingerprint(of: raw))
    }

    @Test func bytesWithANulAreNotText() {
        #expect(CrystalLoading.text(of: Data([0x64, 0x00, 0x61])) == nil)
        // UTF-16 holds NUL bytes for ASCII text.
        #expect(CrystalLoading.text(of: Data("data_x".utf16.flatMap { [UInt8($0 & 0xFF), UInt8($0 >> 8)] })) == nil)
        #expect(CrystalLoading.text(of: Data([0x80, 0xFF, 0x41])) != nil)
        #expect(CrystalLoading.text(of: Data([0xEF, 0xBB, 0xBF, 0x41])) == "A")
    }

    // MARK: Blocks, names and keywords

    /// A publication block before the structure: its space group does not leak into the structure's block (on 0.18.1 its
    /// 'F m -3 m' refused this P1 file), a tag without a value does not swallow the next block or loop, a quoted value that
    /// starts with an underscore is a value, and DDLm (CIF 2) names read as their DDL1 aliases.
    @Test func theStructureBlockIsReadAlone() throws {
        let doc = try Self.read("two_blocks_ddlm")
        #expect(doc.dataBlockName == "gold")
        #expect(doc.displayName == "_gold_")
        #expect(doc.declaredSpaceGroupSymbol == nil)
        #expect(doc.structure.formula == "Au4")
        #expect(doc.structure.laueClass == .mThreeBarM)
    }

    @Test func aTagWithoutAValueDoesNotSwallowTheLoop() throws {
        let cif = Self.gold.replacingOccurrences(of: "loop_\n", with: "_exptl_crystal_colour\nloop_\n")
        #expect(try CIFImport.read(cif).structure.sites.count == 4)
        let atEnd = Self.gold + "_exptl_crystal_colour\n"
        #expect(try CIFImport.read(atEnd).structure.sites.count == 4)
    }

    @Test func macromolecularAndMagneticCIFsAreRefusedByName() throws {
        #expect(throws: CIFImportError.macromolecularCIF) { try Self.read("macromolecular") }
        #expect(throws: CIFImportError.magneticCIF) { try Self.read("magnetic") }
        #expect(CIFImportError.macromolecularCIF.description.contains("mmCIF"))
        #expect(CIFImportError.magneticCIF.description.contains("mCIF"))
        // A DDLm small-molecule CIF with `_atom_site.label` is no mmCIF.
        #expect(try Self.read("two_blocks_ddlm").structure.sites.count == 4)
        // Nor is one that lists Cartesian coordinates beside its fractional ones (`_atom_site.Cartn_x` is core DDLm too, #225).
        let cartesian = Self.gold.replacingOccurrences(of: "_atom_site_", with: "_atom_site.")
            .replacingOccurrences(of: "_atom_site.fract_z\n", with: "_atom_site.fract_z\n_atom_site.Cartn_x\n")
            .replacingOccurrences(of: "Au1 Au 0 0 0", with: "Au1 Au 0 0 0 0")
            .replacingOccurrences(of: "Au2 Au 0.5 0.5 0", with: "Au2 Au 0.5 0.5 0 2.0391")
            .replacingOccurrences(of: "Au3 Au 0.5 0 0.5", with: "Au3 Au 0.5 0 0.5 2.0391")
            .replacingOccurrences(of: "Au4 Au 0 0.5 0.5", with: "Au4 Au 0 0.5 0.5 0")
        #expect(cartesian.contains("_atom_site.Cartn_x"))
        #expect(try CIFImport.read(cartesian).structure.sites.count == 4)
    }

    @Test func theNewErrorsShowTheirText() {
        let errors: [CIFImportError] = [.nonFiniteNumber(tag: "_cell_length_a", value: "nan"), .unknownValue(tag: "_cell_length_a"),
                                        .macromolecularCIF, .magneticCIF]
        for error in errors {
            #expect(!error.description.isEmpty)
            #expect((error as Error).localizedDescription == error.description)
        }
    }
}
