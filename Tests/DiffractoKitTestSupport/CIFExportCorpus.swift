//
//  CIFExportCorpus.swift
//  Role: the structures pins C1 and C2 write as CIF (suite design "Files, figures and exports"): P2's 19 (the CIF fixtures
//        and β″ with its P1 expansion), six of them with an edited lattice, and the three Materials Project imports of the
//        mp-import pins. Built the way the app builds them: `CIFImport.read`, `editingLattice`, `MaterialsProjectImport`.
//

import DiffractoCore
import Foundation

package enum CIFExportCorpus {
    /// The CIF fixtures of P2, by pin name (β″ is built from its generator's text).
    package static let cifFiles: [(name: String, file: String)] = [
        ("gold_P1", "gold_P1.cif"), ("gold_centring_ops", "gold_centring_ops.cif"), ("silicon_Fd-3m", "silicon_Fd-3m.cif"),
        ("magnesium_P63mmc", "magnesium_P63mmc.cif"), ("rocksalt_Fm-3m_192ops", "rocksalt_Fm-3m_192ops.cif"),
        ("alpha_quartz_P3221", "alpha_quartz_P3221.cif"), ("calcite_R-3c_hex", "calcite_R-3c_hex.cif"),
        ("fluorapatite_P63m", "fluorapatite_P63m.cif"), ("pyrite_Pa-3", "pyrite_Pa-3.cif"), ("ws2_2H", "ws2_2H.cif"),
        ("cdi2_P-3m1", "cdi2_P-3m1.cif"), ("tetragonal_a_equals_c", "tetragonal_a_equals_c.cif"),
        ("ideal_hcp_P1", "ideal_hcp_P1.cif"), ("magnesium_P63mmc_24ops", "magnesium_P63mmc_24ops.cif"),
        ("tetragonal_a_equals_c_4ops", "tetragonal_a_equals_c_4ops.cif"),
        ("alpha_quartz_P3221_as_seeded", "alpha_quartz_P3221_as_seeded.cif"), ("pyrite_Pa-3_as_seeded", "pyrite_Pa-3_as_seeded.cif"),
    ]

    /// The six lattice edits: a free parameter of each, moved by a few per cent.
    package static let edits: [(name: String, base: String, parameter: LatticeParameter, value: Double)] = [
        ("silicon_Fd-3m_a5.5", "silicon_Fd-3m", .a, 5.5), ("magnesium_P63mmc_c5.3", "magnesium_P63mmc", .c, 5.3),
        ("rocksalt_Fm-3m_a5.7", "rocksalt_Fm-3m_192ops", .a, 5.7), ("alpha_quartz_P3221_c5.5", "alpha_quartz_P3221", .c, 5.5),
        ("calcite_R-3c_hex_a5.0", "calcite_R-3c_hex", .a, 5.0), ("cdi2_P-3m1_a4.3", "cdi2_P-3m1", .a, 4.3),
    ]

    /// The three served Materials Project inputs of the mp-import pins.
    package static let materialsProject = ["input-aluminium-v3.json", "input-magnesiumJSON.json", "input-siliconJSON.json"]

    /// Every structure of the corpus, by name, in the order above (P2, the edits, the imports).
    package static func structures() throws -> [(name: String, structure: CrystalStructure)] {
        var out: [(name: String, structure: CrystalStructure)] = []
        let beta = try CIFImport.read(Fixtures.betaDoublePrimeCIF, fileBaseName: "beta_double_prime").structure
        out.append(("beta_double_prime", beta))
        out.append(("beta_double_prime_P1", CrystalStructure.detectingSymmetry(cell: beta.cell, sites: beta.sites)))
        for (name, file) in cifFiles {
            out.append((name, try CIFImport.read(try Fixtures.text(file), fileBaseName: name).structure))
        }
        for edit in edits {
            guard let base = out.first(where: { $0.name == edit.base })?.structure else { continue }
            out.append((edit.name, try base.editingLattice([edit.parameter: edit.value])))
        }
        for input in materialsProject {
            let data = try Data(contentsOf: Fixtures.directory.appendingPathComponent("mp-import/" + input))
            let document = try MaterialsProjectImport.firstDocument(in: MaterialsProjectImport.decode(data))
            out.append((input, try MaterialsProjectImport.structure(from: document).structure))
        }
        return out
    }
}
