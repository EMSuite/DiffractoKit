//
//  CIFBlocksTests.swift
//  Role: a CIF of several structures (DiffractoKit #157, plan row K5): its blocks are listed with their names, and one
//        is read by name, with nothing of the other blocks in it. Readers are exact.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import DiffractoCore

struct CIFBlocksTests {

    static func text() throws -> String { try Fixtures.text("cif-edge-cases/three_structures.cif") }

    @Test func theBlocksAreListedWithTheirNames() throws {
        let blocks = try CIFImport.blocks(in: try Self.text())
        #expect(blocks.map(\.name) == ["global", "gold", "Magnesium", "silicon"])
        #expect(blocks.map(\.definesStructure) == [false, true, true, true])
        #expect(blocks.map(\.chemicalName) == [nil, "gold", "Mg", nil])
        #expect(blocks.map(\.siteCount) == [0, 4, 1, 8])
        // The same through the loader, behind its size and encoding refusals.
        let bytes = Data(try Self.text().utf8)
        #expect(try CrystalLoading.blocks(bytes: bytes) == blocks)
        #expect(throws: CrystalLoading.Refusal.notUTF8) { try CrystalLoading.blocks(bytes: Data([0x64, 0x00])) }
    }

    @Test func oneStructureIsReadByName() throws {
        let text = try Self.text()
        let gold = try CIFImport.read(text, block: "gold")
        #expect(gold.dataBlockName == "gold")
        #expect(gold.structure.formula == "Au4")
        #expect(gold.structure.laueClass == .mThreeBarM)
        #expect(gold.declaredSpaceGroupNumber == nil)

        // CIF names are compared without regard to case; the block's own tags only (no other block's number).
        let magnesium = try CIFImport.read(text, block: "magnesium")
        #expect(magnesium.dataBlockName == "Magnesium")
        #expect(magnesium.structure.formula == "Mg2")
        #expect(magnesium.declaredSpaceGroupNumber == 194)
        #expect(magnesium.hadSymmetryOperations)
        #expect(magnesium.structure.laueClass == .sixOverMmm)

        let silicon = try CrystalLoading.load(bytes: Data(text.utf8), fileName: "three.cif", block: "silicon")
        #expect(silicon.structure.formula == "Si8")
        #expect(silicon.document.displayName == "silicon")
        #expect(silicon.structure.cell.a == 5.4309)
    }

    @Test func withoutANameTheFileOfSeveralStructuresIsRefusedWithTheirNames() throws {
        let text = try Self.text()
        #expect(throws: CIFImportError.multipleStructures(blockNames: ["gold", "Magnesium", "silicon"])) {
            try CIFImport.read(text)
        }
        let message = CIFImportError.multipleStructures(blockNames: ["gold", "Magnesium"]).description
        #expect(message.contains("choose one"))
    }

    @Test func aNameTheFileDoesNotHoldIsRefusedWithTheOnesItDoes() throws {
        let text = try Self.text()
        #expect(throws: CIFImportError.blockNotFound("copper", available: ["gold", "Magnesium", "silicon"])) {
            try CIFImport.read(text, block: "copper")
        }
        #expect(CIFImportError.blockNotFound("copper", available: ["gold"]).description.contains("gold"))
        // The publication block holds no structure.
        #expect(throws: CIFImportError.noAtomSites) { try CIFImport.read(text, block: "global") }
    }

    /// A file of one structure reads the same by its name as without one.
    @Test func aSingleStructureReadsTheSameByName() throws {
        let text = try Fixtures.text("gold_P1.cif")
        let blocks = try CIFImport.blocks(in: text)
        let name = try #require(blocks.first(where: \.definesStructure)?.name)
        #expect(try CIFImport.read(text, block: name) == CIFImport.read(text))
    }
}
