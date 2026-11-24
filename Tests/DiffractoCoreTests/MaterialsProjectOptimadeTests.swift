//
//  MaterialsProjectOptimadeTests.swift
//  Role: the keyless Materials Project route (DiffractoKit #150, plan row K3a). Requests are built with no network;
//        answers are decoded from inline bodies shaped like the service's and from the recorded answers in
//        `Fixtures/mp-optimade/` (recorded by `tools/materials-project/record.py` on GitHub's Linux runner, with a
//        manifest of URL, date and SHA-256); an empty service is told apart from no results.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import DiffractoCore

struct MaterialsProjectOptimadeTests {

    typealias Route = MaterialsProjectOptimade

    static func query(of request: URLRequest) -> [String: String] {
        let components = URLComponents(url: request.url!, resolvingAgainstBaseURL: false)!
        var items: [String: String] = [:]
        for item in components.queryItems ?? [] { items[item.name] = item.value }
        return items
    }

    @Test func searchRequestsAreBuiltWithNoNetwork() throws {
        let request = try Route.searchRequest(.elements(["Si", "Mg"], only: true))
        #expect(request.httpMethod == "GET")
        #expect(request.url?.scheme == "https")
        #expect(request.url?.host == "optimade.materialsproject.org")
        #expect(request.url?.path == "/v1/structures")
        #expect(request.value(forHTTPHeaderField: "X-API-KEY") == nil)
        let items = Self.query(of: request)
        #expect(items["filter"] == "elements HAS ALL \"Mg\",\"Si\" AND nelements=2")
        #expect(items["response_fields"] == "chemical_formula_reduced,nsites,_mp_stability")
        #expect(items["page_limit"] == "500")
        // Every reserved character of the filter is percent-encoded, "=" included.
        let raw = try #require(request.url?.absoluteString)
        #expect(raw.contains("filter=elements%20HAS%20ALL%20%22Mg%22%2C%22Si%22%20AND%20nelements%3D2"))

        #expect(try Route.filter(for: .elements(["Al"], only: false)) == "elements HAS ALL \"Al\"")
        #expect(try Route.filter(for: .formula("CuAl2")) == "chemical_formula_reduced=\"Al2Cu\"")
        #expect(try Route.filter(for: .formula("Al4 Cu2")) == "chemical_formula_reduced=\"Al2Cu\"")
        #expect(try Route.filter(for: .formula("O3Al2")) == "chemical_formula_reduced=\"Al2O3\"")
        #expect(try Route.filter(for: .formula("Al")) == "chemical_formula_reduced=\"Al\"")

        let entry = try Route.entryRequest(materialID: "mp-134")
        #expect(entry.url?.absoluteString == "https://optimade.materialsproject.org/v1/structures/mp-134")
    }

    @Test func aRequestThatCannotBeBuiltIsRefused() {
        #expect(throws: Route.RequestError.noElements) { try Route.searchRequest(.elements([], only: false)) }
        #expect(throws: Route.RequestError.unknownElement("Qq")) { try Route.searchRequest(.elements(["Al", "Qq"], only: false)) }
        #expect(throws: Route.RequestError.unknownElement("al")) { try Route.searchRequest(.elements(["al"], only: false)) }
        #expect(throws: Route.RequestError.unreadableFormula("Al(OH)3")) { try Route.searchRequest(.formula("Al(OH)3")) }
        #expect(throws: Route.RequestError.unknownElement("Xy")) { try Route.searchRequest(.formula("Xy2")) }
        #expect(throws: Route.RequestError.invalidMaterialID("mp-12a")) { try Route.entryRequest(materialID: "mp-12a") }
        #expect(throws: Route.RequestError.invalidMaterialID("../info")) { try Route.entryRequest(materialID: "../info") }
    }

    static func searchBody(rows: String, dataAvailable: Int, more: Bool) -> Data {
        Data("""
        {"data":[\(rows)],"included":[],"links":{"next":null},
         "meta":{"api_version":"1.2.0","more_data_available":\(more),"data_returned":0,"data_available":\(dataAvailable)}}
        """.utf8)
    }

    static func row(_ id: String, _ formula: String, _ hull: [String: Double]) -> String {
        let stability = hull.map { "\"\($0.key)\":{\"thermo_id\":\"\(id)_X\",\"energy_above_hull\":\($0.value)}" }
            .joined(separator: ",")
        return """
        {"id":"\(id)","type":"structures","attributes":{"nsites":1,"chemical_formula_reduced":"\(formula)","_mp_stability":{\(stability)}}}
        """
    }

    /// No rows from a service that holds structures is "no results"; no rows from a service that holds none is the
    /// service's own failure (July 2026), told apart so the app can say "try again later".
    @Test func anEmptyServiceIsToldApartFromNoResults() throws {
        #expect(try Route.decodeSearch(Self.searchBody(rows: "", dataAvailable: 154_387, more: false)) == .noResults)
        #expect(try Route.decodeSearch(Self.searchBody(rows: "", dataAvailable: 0, more: false)) == .serviceEmpty)
    }

    @Test func rowsAreSortedByTheLowestEnergyAboveTheHull() throws {
        let rows = [
            Self.row("mp-3", "Al", ["r2scan": 0.2, "gga_gga+u": 0.05]),
            Self.row("mp-2", "Al", ["gga_gga+u": 0.0]),
            Self.row("mp-1", "Al", [:]),
            Self.row("mp-10", "Al", ["gga_gga+u": 0.0]),
        ].joined(separator: ",")
        let answer = try Route.decodeSearch(Self.searchBody(rows: rows, dataAvailable: 154_387, more: true))
        guard case .rows(let decoded, let more) = answer else {
            Issue.record("expected rows, got \(answer)")
            return
        }
        #expect(more)
        // Equal energies by the id's number, as the keyed search orders them (#226): mp-2 before mp-10.
        #expect(decoded.map(\.materialID) == ["mp-2", "mp-10", "mp-3", "mp-1"])
        #expect(decoded[2].lowestEnergyAboveHull == 0.05)
        #expect(decoded[3].lowestEnergyAboveHull == nil)
    }

    @Test func aServiceErrorAndABodyThatIsNotAnAnswerAreRefused() {
        let error = Data(#"{"errors":[{"status":"400","title":"Bad Request","detail":"filter could not be parsed"}],"meta":{}}"#.utf8)
        #expect(throws: Route.AnswerError.serviceError(status: "400", message: "filter could not be parsed")) {
            try Route.decodeSearch(error)
        }
        #expect(throws: Route.AnswerError.self) { try Route.decodeSearch(Data("<html>".utf8)) }
        #expect(throws: Route.AnswerError.self) { try Route.decodeEntry(Data(#"{"data":null,"meta":{}}"#.utf8)) }
    }

    /// mp-134 as the service serves it: the one-atom primitive fcc cell, rounded to 1e-6 Å. It converts to the
    /// conventional Fm-3m cell of four atoms.
    static let aluminiumEntry = Data("""
    {"data":{"id":"mp-134","type":"structures","attributes":{"last_modified":"2023-03-11T22:01:11.232000",
     "elements":["Al"],"nelements":1,"chemical_formula_reduced":"Al",
     "lattice_vectors":[[2.473329,0.0,1.427977],[0.824443,2.331877,1.427977],[0.0,0.0,2.855955]],
     "space_group_it_number":null,"cartesian_site_positions":[[0.0,0.0,0.0]],"nsites":1,
     "species":[{"name":"Al","chemical_symbols":["Al"],"concentration":[1.0],"mass":null}],
     "species_at_sites":["Al"]}},"meta":{"data_available":154387,"more_data_available":false},"links":{"next":null}}
    """.utf8)

    @Test func anEntryDecodesAndConvertsToItsConventionalCell() throws {
        let entry = try Route.decodeEntry(Self.aluminiumEntry)
        #expect(entry.materialID == "mp-134")
        #expect(entry.lastModified == "2023-03-11T22:01:11.232000")
        #expect(entry.latticeVectors[0] == SIMD3(2.473329, 0.0, 1.427977))
        #expect(entry.speciesAtSites == ["Al"])
        let imported = try Route.structure(from: entry)
        #expect(imported.cell == .conventional(centring: "F"))
        #expect(imported.structure.sites.count == 4)
        #expect(imported.structure.laueClass == .mThreeBarM)
        #expect(imported.structure.name == "mp-134")
        let a = imported.structure.cell.a
        #expect(abs(a - 4.039) < 1e-3)
        let cell = imported.structure.cell
        #expect(abs(cell.b - a) < 1e-9 && abs(cell.c - a) < 1e-9)
        #expect([cell.alphaDeg, cell.betaDeg, cell.gammaDeg].allSatisfy { abs($0 - 90) < 1e-9 })
    }

    /// A functional whose stability is null is dropped; the row and the search still decode.
    @Test func aNullStabilityIsDropped() throws {
        let row = #"{"id":"mp-7","type":"structures","attributes":{"nsites":1,"chemical_formula_reduced":"Al","_mp_stability":{"r2scan":null,"gga_gga+u":{"energy_above_hull":0.1}}}}"#
        let answer = try Route.decodeSearch(Self.searchBody(rows: row, dataAvailable: 154_387, more: false))
        guard case .rows(let rows, _) = answer else {
            Issue.record("expected rows, got \(answer)")
            return
        }
        #expect(rows.first?.energyAboveHull == ["gga_gga+u": 0.1])
    }

    static func chiralEntry(_ rows: String) -> Data {
        Data("""
        {"data":{"id":"mp-999","type":"structures","attributes":{
         "lattice_vectors":\(rows),
         "cartesian_site_positions":[[0.0,0.0,0.0],[1.3,0.4,0.2],[0.3,1.7,0.5],[0.6,0.5,2.1]],
         "species_at_sites":["Al","Si","Mg","Fe"],
         "species":[{"name":"Al","chemical_symbols":["Al"],"concentration":[1.0]},
                    {"name":"Si","chemical_symbols":["Si"],"concentration":[1.0]},
                    {"name":"Mg","chemical_symbols":["Mg"],"concentration":[1.0]},
                    {"name":"Fe","chemical_symbols":["Fe"],"concentration":[1.0]}]}},"meta":{}}
        """.utf8)
    }

    /// A left-handed basis (a and b swapped) is the same crystal: the structure is the right-handed one, not its mirror
    /// image (four different atoms at general positions in P1 have no mirror).
    @Test func aLeftHandedBasisIsNotMirrored() throws {
        let right = try Route.structure(from: Route.decodeEntry(Self.chiralEntry("[[5.1,0.0,0.0],[1.2,6.3,0.0],[0.4,0.9,7.2]]")))
        let left = try Route.structure(from: Route.decodeEntry(Self.chiralEntry("[[1.2,6.3,0.0],[5.1,0.0,0.0],[0.4,0.9,7.2]]")))
        #expect(left == right)
        #expect(right.structure.sites.count == 4)
    }

    @Test func anEntryThatCannotBeUsedIsRefused() throws {
        let entry = try Route.decodeEntry(Self.aluminiumEntry)
        func edited(_ species: [String], _ symbols: [String]) -> MaterialsProjectOptimade.Entry {
            MaterialsProjectOptimadeTestEntry.make(entry, speciesAtSites: species, symbols: symbols)
        }
        #expect(throws: MaterialsProjectImport.Failure.self) { try Route.structure(from: edited(["Zz"], ["Al"])) }
        #expect(throws: MaterialsProjectImport.Failure.self) { try Route.structure(from: edited(["Al"], ["Qq"])) }
        #expect(throws: MaterialsProjectImport.Failure.self) { try Route.structure(from: edited([], ["Al"])) }
    }
}

/// Builds edited entries from a decoded one, through the decoder, so the tests use the public route only.
enum MaterialsProjectOptimadeTestEntry {
    static func make(_ entry: MaterialsProjectOptimade.Entry, speciesAtSites: [String], symbols: [String]) -> MaterialsProjectOptimade.Entry {
        let positions = speciesAtSites.isEmpty ? "[]" : "[[0.0,0.0,0.0]]"
        let sites = speciesAtSites.map { "\"\($0)\"" }.joined(separator: ",")
        let chemical = symbols.map { "\"\($0)\"" }.joined(separator: ",")
        let body = """
        {"data":{"id":"\(entry.materialID)","type":"structures","attributes":{
         "lattice_vectors":[[2.473329,0.0,1.427977],[0.824443,2.331877,1.427977],[0.0,0.0,2.855955]],
         "cartesian_site_positions":\(positions),"species_at_sites":[\(sites)],
         "species":[{"name":"Al","chemical_symbols":[\(chemical)],"concentration":[1.0]}]}},"meta":{}}
        """
        return try! MaterialsProjectOptimade.decodeEntry(Data(body.utf8))
    }
}
