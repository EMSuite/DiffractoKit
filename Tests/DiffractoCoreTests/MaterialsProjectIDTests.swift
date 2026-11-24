//
//  MaterialsProjectIDTests.swift
//  Role: Materials Project material ids in both forms (DiffractoKit #203). The main API serves every id in emmet-core
//        0.88.3's alphabetical form (`AlphaID`: base 26, "a" is 0, padded to eight letters); people and the keyless
//        route know the legacy `mp-<digits>`. The pairs below are the ones the service answered on 2026-11-24 (the
//        keyed recordings of #151): asked for mp-134 it served mp-aaaaaafe.
//

import Foundation
import Testing
import DiffractoCore

struct MaterialsProjectIDTests {

    typealias Route = MaterialsProjectImport

    @Test(arguments: [("mp-aaaaaafe", 134), ("mp-aaaaaaft", 149), ("mp-aaacjkpl", 1_079_297), ("mp-aaaaaaaa", 0),
                      ("mp-ft", 149), ("mp-134", 134), ("mp-0134", 134)])
    func bothFormsGiveTheSameNumber(id: String, number: Int) {
        #expect(Route.materialNumber(id) == number)
    }

    @Test func theDisplayFormIsLegacyUpToTheCutPoint() {
        #expect(Route.displayMaterialID("mp-aaaaaafe") == "mp-134")
        #expect(Route.displayMaterialID("mp-aaacjkpl") == "mp-1079297")
        #expect(Route.displayMaterialID("mp-134") == "mp-134")
        #expect(Route.displayMaterialID("mp-3347529") == "mp-3347529")
        #expect(Route.displayMaterialID("mp-3347530") == "mp-aaahilze")
        #expect(Route.displayMaterialID("mp-aaahilze") == "mp-aaahilze")
        #expect(Route.legacyIDLimit == 3_347_529)
    }

    @Test func anythingElseIsRefused() {
        for bad in ["", "mp-", "mp-12a", "MP-12", "12", "mp- 12", "mp-AAAAAAFE", "mp-aaa1", "mp-ä", "../info",
                    "mp-aaaaaaaaaaaaa", "mp-9999999999999"] {
            #expect(Route.materialNumber(bad) == nil, "\(bad)")
            #expect(Route.displayMaterialID(bad) == nil, "\(bad)")
        }
        #expect(throws: Route.Failure.invalidMaterialID("mp-AAAAAAFE")) {
            try Route.request(materialID: "mp-AAAAAAFE", apiKey: "KEY")
        }
    }

    @Test func requestsAskForTheFormPeopleKnow() throws {
        let keyed = try Route.request(materialID: "mp-aaaaaafe", apiKey: "KEY")
        #expect(keyed.url?.query?.contains("material_ids=mp-134&") == true)
        let keyless = try MaterialsProjectOptimade.entryRequest(materialID: "mp-aaaaaafe")
        #expect(keyless.url?.absoluteString == "https://optimade.materialsproject.org/v1/structures/mp-134")
    }

    @Test func servedIDsAreReadInTheFormPeopleKnow() throws {
        let summary = Data(#"{"data":[{"material_id":"mp-aaaaaafe","formula_pretty":"Al"}],"meta":{}}"#.utf8)
        #expect(try Route.firstDocument(in: Route.decode(summary)).materialID == "mp-134")
        let unknown = Data(#"{"data":[{"material_id":"task-12","formula_pretty":"Al"}],"meta":{}}"#.utf8)
        #expect(try Route.firstDocument(in: Route.decode(unknown)).materialID == "task-12")
        let search = Data(#"{"data":[{"material_id":"mp-aaacjkpl","energy_above_hull":0.1},{"material_id":"mp-aaaaaaft","energy_above_hull":0.1}],"meta":{"total_doc":2}}"#.utf8)
        guard case .rows(let rows, _) = try Route.decodeSearch(search) else {
            Issue.record("expected rows")
            return
        }
        #expect(rows.map(\.materialID) == ["mp-149", "mp-1079297"])
    }
}
