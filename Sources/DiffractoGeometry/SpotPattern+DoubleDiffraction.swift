//
//  SpotPattern+DoubleDiffraction.swift
//  Role: double diffraction between two crystals in one beam, a precipitate in its matrix (DiffractoKit S8, #163), on the
//        two drawn spot patterns.
//
//  Rule (Gjønnes and Moodie, Acta Cryst. 19 (1965) 65; Williams and Carter, Transmission Electron Microscopy, 2nd ed. (2009),
//  on double diffraction at precipitates): a beam diffracted by a spot g₁ of one crystal is diffracted again by a spot g₂ of
//  the other and leaves along g₁ + g₂. Every sum of a drawn spot of this pattern and a drawn spot of the other is a
//  position, unless it lies on a kinematic spot of either pattern or on the direct beam (within `DoubleDiffraction.coincidence`,
//  in the pattern plane). Each route carries a weight, the product of the two spots' intensities relative to their own
//  pattern's strongest, used only to order the routes and cut them (default 0.1); it is never an intensity, which kinematic
//  theory does not give. Sums that land on one position are one position with all their routes. No Python package
//  computes this; it is checked against a brute-force sum of the two spot lists.
//  The sum is taken in the pattern plane: each crystal's spots are already the ones its own excitation error lets through.
//  `DoubleDiffraction.spots(crystals:volts:settings:)` is the stricter three-dimensional rule (the sum excited too) with
//  routes within each crystal; this is the one drawn over a precipitate's variants.
//

import Foundation
import simd
import DiffractoCore

/// A position double diffraction adds between two crystals.
public struct PairDoubleDiffractionSpot: Sendable, Equatable {
    /// (qx, qy) of g₁ + g₂, Å⁻¹: the sum of the first pair found in the order of the two spot lists; every route's sum
    /// lies within `DoubleDiffraction.coincidence` of it.
    public let position: SIMD2<Double>
    /// Every route that lands here: `first` a spot of the pattern called on (crystal 0), `second` a spot of the other
    /// (crystal 1); the largest weight first.
    public let routes: [DoubleDiffractionRoute]
}

extension SpotPattern {

    /// The design's route weight cut: 0.1 of the product of the two relative intensities.
    public static let defaultRouteWeightCut = 0.1

    /// The positions double diffraction adds between this pattern and `other` (both in one beam, one frame), with their
    /// routes, the strongest route first. Empty for a cut that is not a finite number of at least 0, or a pattern without
    /// a spot of positive intensity.
    public func doubleDiffraction(with other: SpotPattern, routeWeightCut: Double = defaultRouteWeightCut)
        -> [PairDoubleDiffractionSpot] {
        guard routeWeightCut.isFinite, routeWeightCut >= 0 else { return [] }
        let maxA = spots.map(\.intensity).max() ?? 0, maxB = other.spots.map(\.intensity).max() ?? 0
        guard maxA > 0, maxB > 0 else { return [] }
        let eps = DoubleDiffraction.coincidence
        func cell(_ p: SIMD2<Double>) -> SIMD2<Int> {
            SIMD2(Int((p.x / eps).rounded(.down)), Int((p.y / eps).rounded(.down)))
        }
        func near<T>(_ grid: [SIMD2<Int>: [T]], _ p: SIMD2<Double>, _ position: (T) -> SIMD2<Double>) -> T? {
            let c = cell(p)
            for dx in -1...1 { for dy in -1...1 {
                if let hit = grid[c &+ SIMD2(dx, dy)]?.first(where: { simd_distance(position($0), p) <= eps }) { return hit }
            } }
            return nil
        }
        var kinematic: [SIMD2<Int>: [SIMD2<Double>]] = [:]
        for s in spots + other.spots { kinematic[cell(s.position), default: []].append(s.position) }
        var positions: [SIMD2<Double>] = []
        var routes: [[DoubleDiffractionRoute]] = []
        var grid: [SIMD2<Int>: [Int]] = [:]
        for a in spots {
            let wa = a.intensity / maxA
            for b in other.spots {
                let w = wa * (b.intensity / maxB)
                guard w >= routeWeightCut else { continue }
                let p = a.position + b.position
                guard simd_length(p) > eps, near(kinematic, p, { $0 }) == nil else { continue }
                let route = DoubleDiffractionRoute(first: .init(crystal: 0, hkl: a.hkl), second: .init(crystal: 1, hkl: b.hkl),
                                                   weight: w)
                if let k = near(grid, p, { positions[$0] }) {
                    routes[k].append(route)
                } else {
                    grid[cell(p), default: []].append(positions.count)
                    positions.append(p)
                    routes.append([route])
                }
            }
        }
        func order(_ x: DoubleDiffractionRoute, _ y: DoubleDiffractionRoute) -> Bool {
            if x.weight != y.weight { return x.weight > y.weight }
            if x.first.hkl != y.first.hkl { return x.first.hkl < y.first.hkl }
            return x.second.hkl < y.second.hkl
        }
        return positions.indices.map { PairDoubleDiffractionSpot(position: positions[$0], routes: routes[$0].sorted(by: order)) }
            .sorted { x, y in
                if x.routes[0].weight != y.routes[0].weight { return x.routes[0].weight > y.routes[0].weight }
                if x.position.x != y.position.x { return x.position.x < y.position.x }
                return x.position.y < y.position.y
            }
    }
}
