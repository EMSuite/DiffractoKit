//
//  PatternSpotsTests.swift
//  Role: `PatternSpots` (S4) on planted closed forms and refusals. The py4DSTEM pins are in
//        DiffractoKitParityTests/Geometry/PatternSpotsPinsTests.swift.
//

import Foundation
import Testing
import simd
import DiffractoGeometry

@Suite struct PatternSpotsTests {

    /// A 96 × 96 pattern: a flat background, a beam at `centre` and Gaussian spots (σ 1.5 px) at the given offsets.
    private func planted(centre: SIMD2<Double>, spots: [(SIMD2<Double>, Double)], background: Float = 20) -> DiffractionPattern {
        let n = 96
        var pixels = [Float](repeating: background, count: n * n)
        let all = [(SIMD2<Double>(0, 0), 4000.0)] + spots
        for r in 0..<n { for c in 0..<n {
            var v = 0.0
            for (o, h) in all {
                let d = SIMD2(Double(r), Double(c)) - (centre + o)
                v += h * exp(-simd_length_squared(d) / (2 * 1.5 * 1.5))
            }
            pixels[r * n + c] += Float(v)
        } }
        return DiffractionPattern(width: n, height: n, pixels: pixels)
    }

    @Test func plantedSpotsAreFoundAtTheirCentres() throws {
        let centre = SIMD2(47.3, 48.6)
        let offsets: [(SIMD2<Double>, Double)] = [(SIMD2(20.2, 0.4), 900), (SIMD2(-20.2, -0.4), 800), (SIMD2(0.3, 21.1), 700),
                                                  (SIMD2(-0.3, -21.1), 600)]
        let found = try PatternSpots.find(in: planted(centre: centre, spots: offsets), centre: centre,
                                          settings: .init(edgeBoundary: 5, beamExclusionRadius: 4))
        #expect(found.count == 4)
        for (o, _) in offsets {
            let want = centre + o
            let hit = try #require(found.min { simd_distance(SIMD2($0.row, $0.column), want) < simd_distance(SIMD2($1.row, $1.column), want) })
            // The 7-px window of the default refinement cuts the spot's tail asymmetrically: measured 0.02 to 0.05 px.
            #expect(simd_distance(SIMD2(hit.row, hit.column), want) <= 0.1, "\(hit) against \(want)")
            #expect(hit.isRefined)
        }
        // Brightest first.
        #expect(zip(found, found.dropFirst()).allSatisfy { $0.intensity >= $1.intensity })
        // The radial median takes the beam's core away (its innermost rings are the beam). Without the subtraction the beam is
        // the first spot.
        let withBeam = try PatternSpots.find(in: planted(centre: centre, spots: offsets), centre: centre,
                                             settings: .init(edgeBoundary: 5, subtractsRadialBackground: false))
        #expect(withBeam.count == 5)
        #expect(simd_distance(SIMD2(withBeam[0].row, withBeam[0].column), centre) <= 0.1)
        // Without refinement the maximum's position is the spot's.
        let plain = try PatternSpots.find(in: planted(centre: centre, spots: offsets), centre: centre,
                                          settings: .init(edgeBoundary: 5, refinementHalfWindow: 0))
        #expect(plain.allSatisfy { !$0.isRefined && $0.row == $0.maximumRow && $0.column == $0.maximumColumn })
    }

    @Test func theSpacingKeepsTheBrighterOfTwoCloseSpots() throws {
        let centre = SIMD2(47.0, 47.0)
        let p = planted(centre: centre, spots: [(SIMD2(20, 0), 900), (SIMD2(20, 7), 500)])
        let wide = try PatternSpots.find(in: p, centre: centre, settings: .init(minimumSpacing: 10, edgeBoundary: 5, beamExclusionRadius: 4))
        #expect(wide.count == 1 && abs(wide[0].column - 47) < 0.5)
        let narrow = try PatternSpots.find(in: p, centre: centre, settings: .init(minimumSpacing: 5, edgeBoundary: 5, beamExclusionRadius: 4))
        #expect(narrow.count == 2)
    }

    @Test func refusals() {
        let p = DiffractionPattern(width: 8, height: 8, pixels: [Float](repeating: 1, count: 64))
        #expect(throws: PatternSpots.Failure.emptyPattern) { try PatternSpots.maxima(in: DiffractionPattern(width: 8, height: 8, pixels: [1])) }
        var nan = p.pixels; nan[9] = .nan
        #expect(throws: PatternSpots.Failure.nonFinitePixel) { try PatternSpots.maxima(in: DiffractionPattern(width: 8, height: 8, pixels: nan)) }
        #expect(throws: PatternSpots.Failure.centreOutside) { try PatternSpots.find(in: p, centre: SIMD2(9, 2)) }
        #expect(throws: PatternSpots.Failure.centreOutside) { try PatternSpots.find(in: p, centre: SIMD2(.nan, 2)) }
        for bad in [PatternSpots.Settings(sigma: -1), .init(minimumRelativeIntensity: .nan), .init(relativeToPeak: -1),
                    .init(minimumSpacing: .infinity), .init(edgeBoundary: -2), .init(maximumCount: 0),
                    .init(refinementHalfWindow: -1), .init(beamExclusionRadius: -1)] {
            #expect(throws: (any Error).self) { try PatternSpots.maxima(in: p, settings: bad) }
        }
        // A flat pattern has no strict maximum (edge 1, so the 6 × 6 inside is searched), and a pattern inside its edge
        // boundary has none either.
        #expect((try? PatternSpots.maxima(in: p, settings: .init(sigma: 0, minimumRelativeIntensity: 0, edgeBoundary: 1)))?.isEmpty == true)
        var peak = p.pixels; peak[3 * 8 + 4] = 5
        let one = try? PatternSpots.maxima(in: DiffractionPattern(width: 8, height: 8, pixels: peak),
                                           settings: .init(sigma: 0, minimumRelativeIntensity: 0, edgeBoundary: 1))
        #expect(one?.count == 1 && one?.first?.row == 3 && one?.first?.column == 4)
        #expect((try? PatternSpots.maxima(in: DiffractionPattern(width: 8, height: 8, pixels: peak)))?.isEmpty == true)
    }
}
