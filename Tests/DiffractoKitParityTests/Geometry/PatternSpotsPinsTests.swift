//
//  PatternSpotsPinsTests.swift
//  Role: `PatternSpots` (S4) against py4DSTEM f050d207's get_maxima_2D, read from
//        Tests/DiffractoGeometryTests/Fixtures/pattern-spots-pins.json (tools/pins/pattern-spots/pins.py): six planted spot
//        patterns (256 x 256 counts) and four settings sets. The maxima and their order are exact; positions and intensities
//        at numpy.testing.assert_allclose's float64 default (rtol 1e-7), because py4DSTEM ships no test for the function.
//        Then the pass mark of the plan row: on the same images `find` misses no more true spots and reports no more false
//        ones than py4DSTEM's get_maxima_2D with the same settings.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoGeometry

private struct SpotPins: Decodable {
    struct Settings: Decodable {
        let sigma: Double, minAbsoluteIntensity: Double, minRelativeIntensity: Double, relativeToPeak: Int
        let minSpacing: Double, edgeBoundary: Int, maxNumPeaks: Int
    }
    struct Image: Decodable {
        let label: String
        let centre: [Double]
        /// The counts, uint16 little-endian, row-major, base64. (Named without digits: `convertFromSnakeCase` does not map
        /// "u16le_b64" onto a Swift name the same way on every Foundation.)
        let countsBase64: String
        /// row, column, height over the background, h, k, l
        let truth: [[Double]]
        /// per settings name: row, column, intensity
        let maxima: [String: [[Double]]]
    }
    let versions: [String: String]
    let size: Int
    let settings: [String: Settings]
    let images: [Image]

    static let shared: SpotPins = {
        let url = FixtureDirs.fixtures(of: "DiffractoGeometryTests").appendingPathComponent("pattern-spots-pins.json")
        let decoder = JSONDecoder()
        decoder.keyDecodingStrategy = .convertFromSnakeCase
        do { return try decoder.decode(SpotPins.self, from: Data(contentsOf: url)) } catch { fatalError("pattern-spots-pins.json: \(error)") }
    }()
}

private func pattern(_ image: SpotPins.Image, size: Int) throws -> DiffractionPattern {
    let data = try #require(Data(base64Encoded: image.countsBase64))
    #expect(data.count == size * size * 2)
    var pixels = [Float](repeating: 0, count: size * size)
    data.withUnsafeBytes { raw in
        for i in 0..<(size * size) { pixels[i] = Float(UInt16(littleEndian: raw.loadUnaligned(fromByteOffset: 2 * i, as: UInt16.self))) }
    }
    return DiffractionPattern(width: size, height: size, pixels: pixels)
}

private func kitSettings(_ s: SpotPins.Settings) -> PatternSpots.Settings {
    PatternSpots.Settings(sigma: s.sigma, minimumAbsoluteIntensity: s.minAbsoluteIntensity, minimumRelativeIntensity: s.minRelativeIntensity,
                          relativeToPeak: s.relativeToPeak, minimumSpacing: s.minSpacing, edgeBoundary: s.edgeBoundary,
                          maximumCount: s.maxNumPeaks)
}

/// Misses (true spots of height ≥ 60 counts, the direct beam left out, with no found spot within 3 px) and false spots
/// (found spots farther than 8 px from the centre with no true spot within 3 px).
private func rates(_ found: [SIMD2<Double>], _ image: SpotPins.Image) -> (miss: Int, visible: Int, wrong: Int, found: Int) {
    let centre = SIMD2(image.centre[0], image.centre[1])
    let truth = image.truth.map { SIMD2($0[0], $0[1]) }
    let visible = image.truth.filter { $0[2] >= 60 && !($0[3] == 0 && $0[4] == 0 && $0[5] == 0) }.map { SIMD2($0[0], $0[1]) }
    let off = found.filter { simd_distance($0, centre) > 8 }
    let miss = visible.filter { v in !off.contains { simd_distance($0, v) <= 3 } }.count
    let wrong = off.filter { f in !truth.contains { simd_distance($0, f) <= 3 } }.count
    return (miss, visible.count, wrong, off.count)
}

@Suite("PatternSpots against py4DSTEM f050d207 get_maxima_2D")
struct PatternSpotsPinsTests {

    @Test func pinsCarryTheirVersionsAndSizes() {
        #expect(SpotPins.shared.versions["py4DSTEM"] == "0.14.19")
        #expect(SpotPins.shared.images.count == 6)
        #expect(SpotPins.shared.settings.count == 4)
    }

    /// The same maxima in the same order, positions and intensities at rtol 1e-7.
    @Test(arguments: ["py4dstem_defaults", "kit_defaults", "strict", "unsmoothed"])
    func maximaEqualGetMaxima2D(_ name: String) throws {
        let pins = SpotPins.shared
        let settings = kitSettings(try #require(pins.settings[name]))
        var worst = 0.0
        for image in pins.images {
            let want = try #require(image.maxima[name])
            let got = try PatternSpots.maxima(in: try pattern(image, size: pins.size), settings: settings)
            #expect(got.count == want.count, "\(image.label) \(name): \(got.count) maxima against \(want.count)")
            for (g, w) in zip(got, want) {
                let d = [abs(g.row - w[0]) - 1e-7 * abs(w[0]), abs(g.column - w[1]) - 1e-7 * abs(w[1]),
                         abs(g.intensity - w[2]) - 1e-7 * abs(w[2])]
                worst = max(worst, abs(g.row - w[0]), abs(g.column - w[1]))
                #expect(d.allSatisfy { $0 <= 1e-12 }, "\(image.label) \(name): (\(g.row), \(g.column), \(g.intensity)) against \(w)")
            }
        }
        print("S4WORST \(name) position \(worst)")
    }

    /// The pass mark: `find` (radial median subtracted, the kit's defaults, refined) misses no more and reports no more false
    /// spots than py4DSTEM's get_maxima_2D with the same maxima settings on the raw image. Measured when written: py4DSTEM
    /// 0 of 165 missed, 155 of 414 false; the kit 0 of 165, 153 of 414 (the false spots are noise maxima above py4DSTEM's
    /// relative cut of 0.005; docs/parity.md P41).
    @Test func findIsNoWorseThanPy4DSTEM() throws {
        let pins = SpotPins.shared
        var py = (miss: 0, visible: 0, wrong: 0, found: 0), kit = py
        for image in pins.images {
            let p = try pattern(image, size: pins.size)
            let theirs = rates(try #require(image.maxima["kit_defaults"]).map { SIMD2($0[0], $0[1]) }, image)
            let spots = try PatternSpots.find(in: p, centre: SIMD2(image.centre[0], image.centre[1]))
            let ours = rates(spots.map { SIMD2($0.row, $0.column) }, image)
            #expect(ours.miss <= theirs.miss, "\(image.label): misses \(ours.miss) against \(theirs.miss)")
            #expect(ours.wrong <= theirs.wrong, "\(image.label): false spots \(ours.wrong) against \(theirs.wrong)")
            py = (py.miss + theirs.miss, py.visible + theirs.visible, py.wrong + theirs.wrong, py.found + theirs.found)
            kit = (kit.miss + ours.miss, kit.visible + ours.visible, kit.wrong + ours.wrong, kit.found + ours.found)
        }
        print("S4RATES py4DSTEM miss \(py.miss)/\(py.visible) false \(py.wrong)/\(py.found); kit miss \(kit.miss)/\(kit.visible) false \(kit.wrong)/\(kit.found)")
    }
}
