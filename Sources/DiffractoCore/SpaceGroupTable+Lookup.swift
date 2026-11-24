//
//  SpaceGroupTable+Lookup.swift
//  Role: the setting a CIF with no operation list is expanded in, read from the generated `SpaceGroupTable` in
//        pymatgen-core 2026.10.16 `CifParser.get_symops` order: a Hermann–Mauguin symbol pymatgen reads, else the IT
//        number, else a short or full international symbol of spglib's that pymatgen does not read.
//

extension SpaceGroupTable {
    struct Setting: Equatable {
        let hallSymbol: String
        /// The file's symbol the setting was read from, as written; nil when it was read from the IT number.
        let symbol: String?
    }

    /// The setting of the first of `symbols` pymatgen reads, else of `number`, else of the first of `symbols` spglib
    /// names; nil when none is known, or when a symbol pymatgen would read first is one of `leftOutSymbols`.
    static func setting(symbols: [String], number: Int?) -> Setting? {
        for symbol in symbols {
            let key = sanitized(symbol)
            if let hall = symbolSettings[key] { return setting(hall, symbol: symbol) }
            // pymatgen reads this symbol (before the number) in a setting the kit does not hold: refused, never the
            // number's setting in its place.
            if leftOutSymbols.contains(key) { return nil }
        }
        if let number, (1...numberSettings.count).contains(number) {
            return setting(numberSettings[number - 1], symbol: nil)
        }
        for symbol in symbols {
            if let hall = spglibSymbolSettings[sanitized(symbol)] { return setting(hall, symbol: symbol) }
        }
        return nil
    }

    /// pymatgen's key: white space and underscores removed ("P 6_3/m m c" and "P63/mmc" are one symbol); case kept.
    static func sanitized(_ symbol: String) -> String {
        String(symbol.filter { !$0.isWhitespace && $0 != "_" })
    }

    private static func setting(_ hallNumber: Int, symbol: String?) -> Setting {
        Setting(hallSymbol: hallSymbols[hallNumber - 1], symbol: symbol)
    }
}
