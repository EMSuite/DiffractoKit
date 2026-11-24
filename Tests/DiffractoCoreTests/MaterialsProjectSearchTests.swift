//
//  MaterialsProjectSearchTests.swift
//  Role: the keyed Materials Project search (DiffractoKit #151, plan row K3b). Requests are built with no network;
//        answers are decoded from inline bodies shaped like emmet-core 0.88.3's `SummaryDoc` and, once recorded with
//        the owner's key by `tools/materials-project/record-keyed.py`, from the bodies in `Fixtures/mp-keyed/`, which
//        must decode to the rows its manifest states and hold no key. No test touches the network.
//

import CryptoKit
import DiffractoKitTestSupport
import Foundation
import Testing
import DiffractoCore

struct MaterialsProjectSearchTests {

    typealias Route = MaterialsProjectImport

    static func query(of request: URLRequest) -> [String: String] {
        let components = URLComponents(url: request.url!, resolvingAgainstBaseURL: false)!
        var items: [String: String] = [:]
        for item in components.queryItems ?? [] { items[item.name] = item.value }
        return items
    }

    @Test func searchRequestsAreBuiltWithNoNetwork() throws {
        let request = try Route.searchRequest(.elements(["Si", "Mg"], only: true), apiKey: "KEY")
        #expect(request.httpMethod == "GET")
        #expect(request.url?.scheme == "https")
        #expect(request.url?.host == "api.materialsproject.org")
        #expect(request.url?.path == "/materials/summary")
        #expect(request.value(forHTTPHeaderField: "X-API-KEY") == "KEY")
        #expect(Self.query(of: request) == [
            "chemsys": "Mg-Si", "deprecated": "false", "batch_id_not_eq": "gnome_r2scan_statics",
            "_fields": "material_id,formula_pretty,nsites,energy_above_hull,is_stable,theoretical,deprecated,symmetry",
            "_sort_fields": "energy_above_hull", "_limit": "500",
        ])
        // The key travels in the header only.
        #expect(request.url?.absoluteString.contains("KEY") == false)

        let elements = try Route.searchRequest(.elements(["Si", "Al", "Si"], only: false), apiKey: "KEY")
        #expect(Self.query(of: elements)["elements"] == "Al,Si")
        #expect(Self.query(of: elements)["chemsys"] == nil)
        let formula = try Route.searchRequest(.formula("Al4 Cu2"), apiKey: "KEY")
        #expect(Self.query(of: formula)["formula"] == "Al2Cu")
        #expect(formula.url?.absoluteString.hasPrefix("https://api.materialsproject.org/materials/summary/?formula=Al2Cu&") == true)
        let byID = try Route.request(materialID: "mp-134", apiKey: "KEY")
        #expect(byID.url?.absoluteString.hasPrefix("https://api.materialsproject.org/materials/summary/?") == true)

        let heartbeat = Route.heartbeatRequest()
        #expect(heartbeat.url?.absoluteString == "https://api.materialsproject.org/heartbeat")
        #expect(heartbeat.value(forHTTPHeaderField: "X-API-KEY") == nil)
    }

    @Test func aSearchThatCannotBeBuiltIsRefused() {
        typealias Refusal = MaterialsProjectOptimade.RequestError
        #expect(throws: Refusal.noElements) { try Route.searchRequest(.elements([], only: true), apiKey: "KEY") }
        #expect(throws: Refusal.unknownElement("al")) { try Route.searchRequest(.elements(["al"], only: false), apiKey: "KEY") }
        #expect(throws: Refusal.unreadableFormula("Al(OH)3")) { try Route.searchRequest(.formula("Al(OH)3"), apiKey: "KEY") }
    }

    static func row(_ id: String, eah: String, nsites: Int = 1) -> String {
        """
        {"material_id":"\(id)","formula_pretty":"Al","nsites":\(nsites),"energy_above_hull":\(eah),"is_stable":\(eah == "0.0"),
         "theoretical":false,"deprecated":false,
         "symmetry":{"crystal_system":"Cubic","symbol":"Fm-3m","number":225,"point_group":"m-3m","symprec":0.1,
                     "angle_tolerance":5.0,"version":"2.5.0"}}
        """
    }

    static func body(_ rows: [String], total: Int) -> Data {
        Data("""
        {"data":[\(rows.joined(separator: ","))],"errors":[],
         "meta":{"api_version":"0.47.2","time_stamp":"2026-11-24T00:00:00","total_doc":\(total),"max_limit":1000}}
        """.utf8)
    }

    @Test func rowsDecodeMostStableFirst() throws {
        let answer = try Route.decodeSearch(Self.body(
            [Self.row("mp-3", eah: "0.05"), Self.row("mp-10", eah: "0.0"), Self.row("mp-1", eah: "null"),
             Self.row("mp-2", eah: "0.0")], total: 4))
        guard case .rows(let rows, let more) = answer else {
            Issue.record("expected rows, got \(answer)")
            return
        }
        #expect(!more)
        // Equal energies come in the order of the ids' numbers, not of their text.
        #expect(rows.map(\.materialID) == ["mp-2", "mp-10", "mp-3", "mp-1"])
        #expect(rows[1] == Route.SearchRow(
            materialID: "mp-10", formula: "Al", siteCount: 1, energyAboveHull: 0, isStable: true, theoretical: false,
            deprecated: false, symmetry: MaterialsProjectSymmetry(crystalSystem: "Cubic", symbol: "Fm-3m", number: 225,
                                                                  pointGroup: "m-3m")))
        #expect(rows[2].energyAboveHull == 0.05)
        #expect(rows[3].energyAboveHull == nil)
        let unknownEnergy = try #require(rows.first { $0.materialID == "mp-1" })
        #expect(try Route.decodeSearch(Self.body([Self.row("mp-1", eah: "null")], total: 501)) ==
                .rows([unknownEnergy], more: true))
        #expect(try Route.decodeSearch(Self.body([], total: 0)) == .noResults)
    }

    @Test func anErrorAndABodyThatIsNotAnAnswerAreRefused() {
        let error = Data(#"{"data":[],"errors":[{"code":403,"message":"invalid api key"}],"meta":{}}"#.utf8)
        #expect(throws: Route.Failure.apiError(code: 403, message: "invalid api key")) { try Route.decodeSearch(error) }
        let notAnAnswer = Route.Failure.malformed("the search answer is not the documented summary response.")
        #expect(throws: notAnAnswer) { try Route.decodeSearch(Data(#"{"detail":"Not authenticated"}"#.utf8)) }
        #expect(throws: notAnAnswer) { try Route.decodeSearch(Data("<html>".utf8)) }
    }

    @Test func theHeartbeatGivesTheDatabaseVersion() throws {
        let body = Data(#"{"status":"OK","message":"","version":"0.47.2","db_version":"2026_09_28","access_controlled_batch_ids":["gnome_r2scan_statics"]}"#.utf8)
        #expect(try Route.databaseVersion(from: body) == "2026_09_28")
        let none = Route.Failure.malformed("the heartbeat answer states no database version.")
        #expect(throws: none) { try Route.databaseVersion(from: Data(#"{"status":"OK"}"#.utf8)) }
        #expect(throws: none) { try Route.databaseVersion(from: Data(#"{"db_version":""}"#.utf8)) }
    }

    // MARK: Recorded answers

    static let directory = FixtureDirs.fixtures(of: "DiffractoCoreTests").appendingPathComponent("mp-keyed")
    static let recorded = FileManager.default.fileExists(atPath: directory.appendingPathComponent("manifest.json").path)

    struct Manifest: Decodable {
        struct File: Decodable {
            let file: String
            let url: String
            let sha256: String
            let bytes: Int
        }
        struct Row: Decodable {
            let materialID: String
            let formula: String?
            let siteCount: Int?
            let energyAboveHull: Double?
            let isStable: Bool?
            let theoretical: Bool?
            let deprecated: Bool?
            let symmetry: Symmetry?
        }
        struct Symmetry: Decodable {
            let crystalSystem: String?
            let symbol: String?
            let number: Int?
            let pointGroup: String?
            let hall: String?
        }
        struct Search: Decodable {
            let file: String
            let more: Bool
            let rows: [Row]
        }
        let databaseVersion: String
        let files: [File]
        let searches: [Search]
        let summaries: [String]
    }

    /// Recorded summaries the keyed route refuses, as designed: a centred primitive cell whose stated system does not
    /// fit is refused, not guessed. mp-623511 (Sn) is stated I4/mmm at the Materials Project's symprec 0.1, but its
    /// conventional cell has a = 3.453 and b = 3.551 (2.8 % apart), so it is body-centred orthorhombic at the kit's
    /// tolerance; the keyless route, which states no symmetry, reads it as oI (DiffractoKit #218).
    static let refusedSummaries: Set<String> = ["mp-623511"]

    static func manifest() throws -> Manifest {
        try JSONDecoder().decode(Manifest.self, from: Data(contentsOf: directory.appendingPathComponent("manifest.json")))
    }

    @Test(.enabled(if: MaterialsProjectSearchTests.recorded, "the keyed answers are recorded on the owner's Mac (record-keyed.py)"))
    func recordedAnswersDecodeExactly() throws {
        let manifest = try Self.manifest()
        for file in manifest.files {
            let data = try Data(contentsOf: Self.directory.appendingPathComponent(file.file))
            #expect(data.count == file.bytes, "\(file.file)")
            #expect(SHA256.hash(data: data).map { String(format: "%02x", $0) }.joined() == file.sha256, "\(file.file)")
        }
        let heartbeat = try Data(contentsOf: Self.directory.appendingPathComponent("heartbeat.json"))
        #expect(try Route.databaseVersion(from: heartbeat) == manifest.databaseVersion)
        #expect(!manifest.searches.isEmpty)
        for search in manifest.searches {
            let answer = try Route.decodeSearch(Data(contentsOf: Self.directory.appendingPathComponent(search.file)))
            let expected: Route.SearchAnswer = search.rows.isEmpty ? .noResults : .rows(search.rows.map { row in
                // The manifest keeps the ids as served (alphabetical); the rows give the form people know.
                Route.SearchRow(
                    materialID: Route.displayMaterialID(row.materialID) ?? row.materialID, formula: row.formula, siteCount: row.siteCount,
                    energyAboveHull: row.energyAboveHull, isStable: row.isStable, theoretical: row.theoretical,
                    deprecated: row.deprecated,
                    symmetry: row.symmetry.map {
                        MaterialsProjectSymmetry(crystalSystem: $0.crystalSystem, symbol: $0.symbol, number: $0.number,
                                                 pointGroup: $0.pointGroup, hall: $0.hall)
                    })
            }, more: search.more)
            #expect(answer == expected, "\(search.file)")
        }
        for id in manifest.summaries {
            let response = try Route.decode(Data(contentsOf: Self.directory.appendingPathComponent("summary-\(id).json")))
            let document = try Route.firstDocument(in: response)
            #expect(document.materialID == id)
            if Self.refusedSummaries.contains(id) {
                // The case and the id, not the sentence: a rewording of the message keeps the behaviour pinned.
                do {
                    _ = try Route.structure(from: document)
                    Issue.record("\(id) was read; the keyed route refuses it (#218)")
                } catch {
                    guard case .nonStandardCell(let text) = error else {
                        Issue.record("\(id): expected nonStandardCell, got \(error)")
                        continue
                    }
                    #expect(text.hasPrefix("\(id):"), "\(text)")
                }
            } else {
                _ = try Route.structure(from: document)
            }
        }
        #expect(Self.refusedSummaries.isSubset(of: manifest.summaries))
    }

    /// No recorded Materials Project answer holds a key: no header name, and no 32-character run of letters and digits
    /// (the form of a Materials Project key) standing on its own. Runs over every recording, keyed or not.
    @Test func noRecordedAnswerHoldsAKey() throws {
        let fixtures = FixtureDirs.fixtures(of: "DiffractoCoreTests")
        let token = try NSRegularExpression(pattern: "(?<![A-Za-z0-9])[A-Za-z0-9]{32}(?![A-Za-z0-9])")
        var searched = 0
        for folder in ["mp-keyed", "mp-optimade", "mp-import"] {
            let directory = fixtures.appendingPathComponent(folder)
            guard let names = try? FileManager.default.contentsOfDirectory(atPath: directory.path) else { continue }
            for name in names {
                let text = String(decoding: try Data(contentsOf: directory.appendingPathComponent(name)), as: UTF8.self)
                let lowered = text.lowercased()
                #expect(!lowered.contains("x-api-key") && !lowered.contains("api_key") && !lowered.contains("apikey"),
                        "\(folder)/\(name)")
                let range = NSRange(text.startIndex..., in: text)
                #expect(token.firstMatch(in: text, range: range) == nil, "\(folder)/\(name)")
                searched += 1
            }
        }
        #expect(searched > 0)
    }
}
