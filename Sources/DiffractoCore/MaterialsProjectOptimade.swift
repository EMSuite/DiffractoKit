//
//  MaterialsProjectOptimade.swift
//  Role: the keyless Materials Project route (OPTIMADE, https://optimade.materialsproject.org/v1): the search and entry
//        requests (built here, sent by the app), the decoding of their answers, and the conversion of an entry to a
//        structure in its conventional cell. No key. The server rounds lattice vectors and Cartesian positions to 1e-6 Å
//        and states no space group, so the symmetry is read from the atoms at a metric tolerance derived from that
//        rounding (`metricLengthTolerance`). Checked against recorded answers of the service and against pymatgen.
//

import Foundation
import simd

/// The keyless Materials Project route: OPTIMADE requests, answers and structures.
public enum MaterialsProjectOptimade {

    /// The service's base URL (HTTPS, the only host this route names).
    public static let baseURL = URL(string: "https://optimade.materialsproject.org/v1")!
    /// The rows a search asks for at most: more matches than this ask the user to narrow the search.
    public static let maximumRows = 500
    /// The three small fields of a search row.
    public static let searchFields = "chemical_formula_reduced,nsites,_mp_stability"
    /// The relative metric tolerance of the conversion. The server rounds every lattice vector component to 1e-6 Å, so
    /// an entry of the metric G = L Lᵀ moves by up to about 2 × 0.5e-6 × |a| Å² against its scale |a|², a relative
    /// 1e-6 / |a|; with |a| ≥ 1 Å that is at most 1e-6. The tolerance is ten times that. It is still far below any
    /// real departure from a lattice's symmetry that a calculated cell shows.
    public static let metricLengthTolerance = 1e-5

    /// What a search asks for.
    public enum Query: Sendable, Equatable {
        /// Entries holding all these elements; `only`: and no other element.
        case elements([String], only: Bool)
        /// Entries of this formula, in any order of its elements and any multiple of its counts ("CuAl2", "Al4Cu2").
        case formula(String)
    }

    /// Why a request could not be built.
    public enum RequestError: Error, Sendable, Equatable, CustomStringConvertible {
        case noElements
        case unknownElement(String)
        case unreadableFormula(String)
        case invalidMaterialID(String)

        public var description: String {
            switch self {
            case .noElements: "Name at least one element to search for."
            case .unknownElement(let symbol): "\"\(symbol)\" is not an element symbol."
            case .unreadableFormula(let text): "\"\(text)\" is not a formula this search can read (element symbols with whole counts, such as Al2Cu)."
            case .invalidMaterialID(let id): "\"\(id)\" is not a Materials Project id (expected \"mp-\" and digits)."
            }
        }
    }

    // MARK: Requests

    /// The OPTIMADE filter of a query: `elements HAS ALL "Mg","Si"` (with `AND nelements=2` for `only`), or
    /// `chemical_formula_reduced="Al2Cu"` (the elements in alphabetical order, the counts divided by their greatest
    /// common divisor, a count of 1 left out: the form the server stores).
    public static func filter(for query: Query) throws(RequestError) -> String {
        switch query {
        case .elements(let raw, let only):
            var symbols: [String] = []
            for symbol in raw {
                let trimmed = symbol.trimmingCharacters(in: .whitespaces)
                guard Element.atomicNumber(forSymbol: trimmed) != nil,
                      trimmed == Element.normalizedSymbol(fromCIF: trimmed) else { throw .unknownElement(symbol) }
                if !symbols.contains(trimmed) { symbols.append(trimmed) }
            }
            guard !symbols.isEmpty else { throw .noElements }
            let list = symbols.sorted().map { "\"\($0)\"" }.joined(separator: ",")
            return "elements HAS ALL \(list)" + (only ? " AND nelements=\(symbols.count)" : "")
        case .formula(let text):
            return "chemical_formula_reduced=\"\(try reducedFormula(text))\""
        }
    }

    /// The formula as OPTIMADE's `chemical_formula_reduced` holds it.
    public static func reducedFormula(_ text: String) throws(RequestError) -> String {
        let trimmed = text.replacingOccurrences(of: " ", with: "")
        guard let counts = MaterialsProjectImport.parseFormula(trimmed), !counts.isEmpty else {
            throw .unreadableFormula(text)
        }
        for symbol in counts.keys where Element.atomicNumber(forSymbol: symbol) == nil { throw .unknownElement(symbol) }
        func gcd(_ a: Int, _ b: Int) -> Int { b == 0 ? a : gcd(b, a % b) }
        let divisor = counts.values.reduce(0, gcd)
        return counts.keys.sorted().map { symbol in
            let n = counts[symbol]! / max(divisor, 1)
            return n == 1 ? symbol : "\(symbol)\(n)"
        }.joined()
    }

    /// GET /structures with the query's filter, the three small fields and at most `maximumRows` rows.
    public static func searchRequest(_ query: Query) throws(RequestError) -> URLRequest {
        let items = [("filter", try filter(for: query)), ("response_fields", searchFields),
                     ("page_limit", String(maximumRows))]
        return get(path: "structures", items: items)
    }

    /// GET /structures/<id>: one entry with its lattice and sites.
    public static func entryRequest(materialID: String) throws(RequestError) -> URLRequest {
        guard MaterialsProjectImport.isValidMaterialID(materialID) else { throw .invalidMaterialID(materialID) }
        return get(path: "structures/\(materialID)", items: [])
    }

    private static func get(path: String, items: [(String, String)]) -> URLRequest {
        var components = URLComponents(url: baseURL.appendingPathComponent(path), resolvingAgainstBaseURL: false)!
        // Every character but the unreserved ones is percent-encoded: a filter holds spaces, quotes, commas and "=".
        let unreserved = CharacterSet(charactersIn: "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-._~")
        if !items.isEmpty {
            components.percentEncodedQuery = items.map { name, value in
                "\(name)=\(value.addingPercentEncoding(withAllowedCharacters: unreserved)!)"
            }.joined(separator: "&")
        }
        var request = URLRequest(url: components.url!)
        request.httpMethod = "GET"
        request.setValue("application/vnd.api+json, application/json", forHTTPHeaderField: "Accept")
        return request
    }

    // MARK: Answers

    /// One search row.
    public struct Row: Sendable, Equatable {
        public let materialID: String
        public let formula: String?
        public let siteCount: Int?
        /// Energy above the hull (eV/atom) by functional ("r2scan", "gga_gga+u", "gga_gga+u_r2scan"; the keys vary by entry).
        public let energyAboveHull: [String: Double]

        /// The lowest energy above the hull over the functionals given; nil when none is.
        public var lowestEnergyAboveHull: Double? { energyAboveHull.values.filter(\.isFinite).min() }
    }

    /// What a search answered.
    public enum SearchAnswer: Sendable, Equatable {
        /// Matching rows, the most stable first (lowest energy above the hull, then the id); `more`: the server holds more
        /// matches than were returned (narrow the search).
        case rows([Row], more: Bool)
        /// The service holds structures, and none matches.
        case noResults
        /// The service answered with no structures at all (`data_available` 0): its data is missing, as for two days
        /// in July 2026. Not "no results": the search should be tried again later.
        case serviceEmpty
    }

    /// Why an answer could not be used.
    public enum AnswerError: Error, Sendable, Equatable, CustomStringConvertible {
        /// The body is not the OPTIMADE envelope.
        case notAnAnswer(String)
        /// The service sent errors (JSON:API `errors`: status, title, detail).
        case serviceError(status: String?, message: String)

        public var description: String {
            switch self {
            case .notAnAnswer(let reason): "The Materials Project answer could not be read: \(reason)"
            case .serviceError(let status, let message):
                "The Materials Project service answered with an error\(status.map { " (\($0))" } ?? ""): \(message)"
            }
        }
    }

    /// Decodes a search answer. Rows are sorted by stability here, since the server cannot sort on `_mp_stability`.
    public static func decodeSearch(_ data: Data) throws(AnswerError) -> SearchAnswer {
        let envelope: SearchEnvelope = try decodeEnvelope(data)
        let rows = (envelope.data ?? []).map { item in
            Row(materialID: item.id, formula: item.attributes?.chemicalFormulaReduced, siteCount: item.attributes?.nsites,
                energyAboveHull: (item.attributes?.mpStability ?? [:]).compactMapValues { $0?.energyAboveHull })
        }
        if rows.isEmpty {
            return envelope.meta?.dataAvailable == 0 ? .serviceEmpty : .noResults
        }
        let sorted = rows.sorted { lhs, rhs in
            let l = lhs.lowestEnergyAboveHull ?? .infinity, r = rhs.lowestEnergyAboveHull ?? .infinity
            return l != r ? l < r : lhs.materialID < rhs.materialID
        }
        return .rows(sorted, more: envelope.meta?.moreDataAvailable ?? false)
    }

    /// One entry: the lattice vectors and Cartesian positions as served (Å, rounded to 1e-6), the species.
    public struct Entry: Sendable, Equatable {
        public let materialID: String
        public let formula: String?
        public let lastModified: String?
        public let latticeVectors: [SIMD3<Double>]
        public let cartesianPositions: [SIMD3<Double>]
        public let speciesAtSites: [String]
        public let species: [Species]
    }

    /// An OPTIMADE species: a name used by `speciesAtSites`, its chemical symbols and their concentrations.
    public struct Species: Sendable, Equatable {
        public let name: String
        public let chemicalSymbols: [String]
        public let concentration: [Double]
    }

    /// Decodes an entry answer (`GET /structures/<id>`).
    public static func decodeEntry(_ data: Data) throws(AnswerError) -> Entry {
        let envelope: EntryEnvelope = try decodeEnvelope(data)
        guard let item = envelope.data else { throw .notAnAnswer("it holds no entry.") }
        let a = item.attributes
        func vectors(_ rows: [[Double]]?, _ what: String) throws(AnswerError) -> [SIMD3<Double>] {
            guard let rows, rows.allSatisfy({ $0.count == 3 }) else { throw .notAnAnswer("\(item.id) has no \(what).") }
            return rows.map { SIMD3($0[0], $0[1], $0[2]) }
        }
        let lattice = try vectors(a?.latticeVectors, "lattice vectors")
        guard lattice.count == 3 else { throw .notAnAnswer("\(item.id) has \(lattice.count) lattice vectors, not 3.") }
        return Entry(
            materialID: item.id, formula: a?.chemicalFormulaReduced, lastModified: a?.lastModified,
            latticeVectors: lattice, cartesianPositions: try vectors(a?.cartesianSitePositions, "site positions"),
            speciesAtSites: a?.speciesAtSites ?? [],
            species: (a?.species ?? []).map {
                Species(name: $0.name, chemicalSymbols: $0.chemicalSymbols ?? [], concentration: $0.concentration ?? [])
            })
    }

    // MARK: Structure

    /// The entry as a structure in the conventional cell of its lattice, named by its id. The fractional positions are
    /// the Cartesian ones in the served lattice; the lattice and its conventional cell are read at
    /// `metricLengthTolerance`, the cell's parameters idealised to its family, the sites carried into it with the
    /// lattice centring; the symmetry is then detected there. Refused (`MaterialsProjectImport.Failure.malformed`):
    /// a number that is not finite, a site without a species, an unknown element, a lattice that cannot be reduced.
    public static func structure(from entry: Entry) throws(MaterialsProjectImport.Failure) -> MaterialsProjectImport.ImportedStructure {
        let id = entry.materialID
        var rows = entry.latticeVectors
        guard rows.count == 3, rows.allSatisfy({ $0.x.isFinite && $0.y.isFinite && $0.z.isFinite }) else {
            throw .malformed("\(id)'s lattice holds a number that is not finite.")
        }
        guard !entry.cartesianPositions.isEmpty, entry.cartesianPositions.count == entry.speciesAtSites.count else {
            throw .malformed("\(id) lists \(entry.cartesianPositions.count) positions and \(entry.speciesAtSites.count) species.")
        }
        let servedVolume = simd_double3x3(columns: (rows[0], rows[1], rows[2])).determinant
        guard servedVolume.isFinite, abs(servedVolume) > 0 else { throw .malformed("\(id)'s lattice has no volume.") }
        // A left-handed basis is the same lattice with a and b swapped (the fractional x and y swap with them); a cell
        // rebuilt from lengths and angles is always right-handed, so without the swap a chiral structure would mirror.
        if servedVolume < 0 { rows.swapAt(0, 1) }
        let basis = simd_double3x3(columns: (rows[0], rows[1], rows[2]))
        let toFractional = basis.inverse
        let speciesByName = Dictionary(entry.species.map { ($0.name, $0) }, uniquingKeysWith: { first, _ in first })
        var sites: [AtomSite] = []
        for (index, (position, name)) in zip(entry.cartesianPositions, entry.speciesAtSites).enumerated() {
            guard position.x.isFinite, position.y.isFinite, position.z.isFinite else {
                throw .malformed("Site \(index) of \(id) has a position that is not finite.")
            }
            let fractional = toFractional * position
            guard let species = speciesByName[name], !species.chemicalSymbols.isEmpty,
                  species.chemicalSymbols.count == species.concentration.count else {
                throw .malformed("Site \(index) of \(id) names the species \"\(name)\", which the entry does not describe.")
            }
            for (symbol, concentration) in zip(species.chemicalSymbols, species.concentration) {
                // OPTIMADE writes a vacancy share as the symbol "vacancy"; it is no atom.
                if symbol == "vacancy" { continue }
                guard concentration.isFinite else {
                    throw .malformed("Site \(index) of \(id) has a concentration that is not finite.")
                }
                guard let z = Element.atomicNumber(forSymbol: symbol) else {
                    throw .malformed("Site \(index) of \(id) has an unrecognized element symbol \"\(symbol)\".")
                }
                sites.append(AtomSite(z: z, fractional: fractional, occupancy: concentration))
            }
        }
        guard !sites.isEmpty else { throw .malformed("\(id) lists no atoms.") }

        func angle(_ u: SIMD3<Double>, _ v: SIMD3<Double>) -> Double {
            acos(min(1, max(-1, dot(u, v) / (length(u) * length(v))))) * 180 / .pi
        }
        let served = UnitCell(a: length(rows[0]), b: length(rows[1]), c: length(rows[2]),
                              alphaDeg: angle(rows[1], rows[2]), betaDeg: angle(rows[0], rows[2]),
                              gammaDeg: angle(rows[0], rows[1]))
        let read = CrystalStructure.detectingSymmetry(cell: served, sites: sites, name: id,
                                                      metricLengthTolerance: metricLengthTolerance)
        guard let conventional = read.conventionalSetting, conventional.reductionConverged else {
            throw .malformed("\(id)'s lattice could not be reduced to find its conventional cell.")
        }
        let structure = CrystalStructure.detectingSymmetry(cell: conventional.cell, sites: conventional.sites, name: id,
                                                           metricLengthTolerance: metricLengthTolerance)
        // A P lattice that was only reduced or relabelled is no conversion from a centred primitive cell.
        let centring = conventional.bravais.centring
        let outcome: CellStandardisation.Outcome = conventional.transformation == .identity || centring == "P"
            ? .asServed : .conventional(centring: centring)
        return MaterialsProjectImport.ImportedStructure(structure: structure, cell: outcome)
    }

    // MARK: Envelopes

    private static func decodeEnvelope<T: Decodable & HasErrors>(_ data: Data) throws(AnswerError) -> T {
        let envelope: T
        do {
            envelope = try JSONDecoder().decode(T.self, from: data)
        } catch {
            throw .notAnAnswer("\(error)")
        }
        if let first = envelope.errors?.first {
            throw .serviceError(status: first.status, message: first.detail ?? first.title ?? "no detail given")
        }
        return envelope
    }

    private protocol HasErrors { var errors: [ErrorObject]? { get } }

    private struct ErrorObject: Decodable {
        let status: String?
        let title: String?
        let detail: String?
    }

    private struct Meta: Decodable {
        let dataAvailable: Int?
        let moreDataAvailable: Bool?
        enum CodingKeys: String, CodingKey {
            case dataAvailable = "data_available", moreDataAvailable = "more_data_available"
        }
    }

    private struct Stability: Decodable {
        let energyAboveHull: Double?
        enum CodingKeys: String, CodingKey { case energyAboveHull = "energy_above_hull" }
    }

    private struct SpeciesObject: Decodable {
        let name: String
        let chemicalSymbols: [String]?
        let concentration: [Double]?
        enum CodingKeys: String, CodingKey {
            case name, chemicalSymbols = "chemical_symbols", concentration
        }
    }

    private struct Attributes: Decodable {
        let chemicalFormulaReduced: String?
        let nsites: Int?
        // A functional with a null entry is dropped, not a reason to refuse the whole search.
        let mpStability: [String: Stability?]?
        let lastModified: String?
        let latticeVectors: [[Double]]?
        let cartesianSitePositions: [[Double]]?
        let speciesAtSites: [String]?
        let species: [SpeciesObject]?
        enum CodingKeys: String, CodingKey {
            case chemicalFormulaReduced = "chemical_formula_reduced", nsites, mpStability = "_mp_stability"
            case lastModified = "last_modified", latticeVectors = "lattice_vectors"
            case cartesianSitePositions = "cartesian_site_positions", speciesAtSites = "species_at_sites", species
        }
    }

    private struct Item: Decodable {
        let id: String
        let attributes: Attributes?
    }

    private struct SearchEnvelope: Decodable, HasErrors {
        let data: [Item]?
        let meta: Meta?
        let errors: [ErrorObject]?
    }

    private struct EntryEnvelope: Decodable, HasErrors {
        let data: Item?
        let meta: Meta?
        let errors: [ErrorObject]?
    }
}
