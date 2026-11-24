//
//  MaterialsProjectSearch.swift
//  Role: the keyed Materials Project search (DiffractoKit #151, plan row K3b): `GET /materials/summary/` by elements,
//        chemical system or formula with the user's own key, sorted on the server by the energy above the hull, the
//        decoding of its rows, and `GET /heartbeat` for the database version a phase records. The requests are built
//        here and sent by the app; the kit holds no key and does no network access. The chosen row's structure is
//        fetched with `MaterialsProjectImport.request(materialID:apiKey:)` and converted by `structure(from:)`.
//
//  The parameter names are mp-api 0.47.2's (`SummaryRester.search`: `elements`, `chemsys`, `formula`, `deprecated`,
//  `_fields`, `_sort_fields`, `_limit`, and `batch_id_not_eq=gnome_r2scan_statics`, its `include_gnome=False`), and the
//  row fields emmet-core 0.88.3's `SummaryDoc`. GNoME entries are left out because their licence (CC BY-NC) is not the
//  rest of the database's (CC BY 4.0).
//

import Foundation

extension MaterialsProjectImport {

    /// The main API's base URL (HTTPS, the only host this route names).
    public static let baseURL = URL(string: "https://api.materialsproject.org")!
    /// The rows a search asks for at most, as for the keyless route: more matches than this ask to narrow the search.
    public static let maximumRows = 500
    /// The small fields of a search row.
    public static let searchFields = "material_id,formula_pretty,nsites,energy_above_hull,is_stable,theoretical,deprecated,symmetry"
    /// The batch a search leaves out: GNoME's structures, licensed CC BY-NC.
    public static let excludedBatch = "gnome_r2scan_statics"

    /// `GET /materials/summary/` for a query with the user's key: `elements=Mg,Si` (entries holding all of them),
    /// `chemsys=Mg-Si` (`only`: those elements and no other), or `formula=Al2Cu` (reduced, as the keyless route sends
    /// it); deprecated entries and GNoME's left out, the most stable first, at most `maximumRows` rows.
    public static func searchRequest(_ query: MaterialsProjectOptimade.Query,
                                     apiKey: String) throws(MaterialsProjectOptimade.RequestError) -> URLRequest {
        var items: [(String, String)]
        switch query {
        case .elements(let raw, let only):
            let symbols = try MaterialsProjectOptimade.elementSymbols(raw)
            items = only ? [("chemsys", symbols.joined(separator: "-"))] : [("elements", symbols.joined(separator: ","))]
        case .formula(let text):
            items = [("formula", try MaterialsProjectOptimade.reducedFormula(text))]
        }
        items += [("deprecated", "false"), ("batch_id_not_eq", excludedBatch), ("_fields", searchFields),
                  ("_sort_fields", "energy_above_hull"), ("_limit", String(maximumRows))]
        var request = get(path: "materials/summary/", items: items)
        request.setValue(apiKey, forHTTPHeaderField: "X-API-KEY")
        return request
    }

    /// `GET /heartbeat`: the service's state and its database version. It needs no key and is sent without one.
    public static func heartbeatRequest() -> URLRequest {
        get(path: "heartbeat", items: [])
    }

    private static func get(path: String, items: [(String, String)]) -> URLRequest {
        // Spelled out, not appended: the summary path keeps its trailing slash, as the service routes it.
        var components = URLComponents(string: baseURL.absoluteString + "/" + path)!
        if !items.isEmpty { components.queryItems = items.map { URLQueryItem(name: $0.0, value: $0.1) } }
        var request = URLRequest(url: components.url!)
        request.httpMethod = "GET"
        request.setValue("application/json", forHTTPHeaderField: "Accept")
        return request
    }

    /// One search row.
    public struct SearchRow: Sendable, Equatable {
        public let materialID: String
        public let formula: String?
        public let siteCount: Int?
        /// Energy above the hull (eV/atom) of the database's preferred functional mix.
        public let energyAboveHull: Double?
        public let isStable: Bool?
        /// No experimental (ICSD) entry is known for the structure: calculated only.
        public let theoretical: Bool?
        public let deprecated: Bool?
        /// The space group the Materials Project states (symbol, number, crystal system), found at its own tolerance.
        public let symmetry: MaterialsProjectSymmetry?

        public init(materialID: String, formula: String?, siteCount: Int?, energyAboveHull: Double?, isStable: Bool?,
                    theoretical: Bool?, deprecated: Bool?, symmetry: MaterialsProjectSymmetry?) {
            self.materialID = materialID
            self.formula = formula
            self.siteCount = siteCount
            self.energyAboveHull = energyAboveHull
            self.isStable = isStable
            self.theoretical = theoretical
            self.deprecated = deprecated
            self.symmetry = symmetry
        }
    }

    /// What a search answered.
    public enum SearchAnswer: Sendable, Equatable {
        /// Matching rows, the most stable first (lowest energy above the hull, a row without one last, then the id);
        /// `more`: the server holds more matches than were returned (narrow the search).
        case rows([SearchRow], more: Bool)
        /// No entry matches.
        case noResults
    }

    /// Decodes a search answer. The server sorts; the rows are sorted again here so that equal energies come in the
    /// order of their ids. An `errors` entry throws `Failure.apiError`; a body that is not the summary envelope throws
    /// `Failure.malformed`.
    public static func decodeSearch(_ data: Data) throws(Failure) -> SearchAnswer {
        let envelope: SearchEnvelope
        do {
            envelope = try JSONDecoder().decode(SearchEnvelope.self, from: data)
        } catch {
            throw .malformed("the search answer is not the documented summary response.")
        }
        if let first = envelope.errors?.first { throw .apiError(code: first.code, message: first.message) }
        // A refusal (`{"detail": "..."}`, the service's 401, 403 and 422 bodies) has no data array: not "no results".
        guard let data = envelope.data else {
            throw .malformed("the search answer is not the documented summary response.")
        }
        let rows = data.map { item in
            SearchRow(materialID: item.materialID, formula: item.formulaPretty, siteCount: item.nsites,
                      energyAboveHull: item.energyAboveHull, isStable: item.isStable, theoretical: item.theoretical,
                      deprecated: item.deprecated, symmetry: item.symmetry)
        }
        guard !rows.isEmpty else { return .noResults }
        let sorted = rows.sorted { lhs, rhs in
            let l = lhs.energyAboveHull ?? .infinity, r = rhs.energyAboveHull ?? .infinity
            return l != r ? l < r : lhs.materialID < rhs.materialID
        }
        return .rows(sorted, more: (envelope.meta?.totalDoc ?? sorted.count) > sorted.count)
    }

    /// The database version a heartbeat answer states ("2026_09_28"), which a phase from this route records.
    public static func databaseVersion(from data: Data) throws(Failure) -> String {
        guard let heartbeat = try? JSONDecoder().decode(Heartbeat.self, from: data),
              let version = heartbeat.dbVersion, !version.isEmpty else {
            throw .malformed("the heartbeat answer states no database version.")
        }
        return version
    }

    private struct SearchEnvelope: Decodable {
        struct Item: Decodable {
            let materialID: String
            let formulaPretty: String?
            let nsites: Int?
            let energyAboveHull: Double?
            let isStable: Bool?
            let theoretical: Bool?
            let deprecated: Bool?
            let symmetry: MaterialsProjectSymmetry?

            enum CodingKeys: String, CodingKey {
                case materialID = "material_id"
                case formulaPretty = "formula_pretty"
                case nsites
                case energyAboveHull = "energy_above_hull"
                case isStable = "is_stable"
                case theoretical, deprecated, symmetry
            }
        }
        struct Meta: Decodable {
            let totalDoc: Int?
            enum CodingKeys: String, CodingKey { case totalDoc = "total_doc" }
        }
        let data: [Item]?
        let errors: [MaterialsProjectError]?
        let meta: Meta?
    }

    private struct Heartbeat: Decodable {
        let dbVersion: String?
        enum CodingKeys: String, CodingKey { case dbVersion = "db_version" }
    }
}
