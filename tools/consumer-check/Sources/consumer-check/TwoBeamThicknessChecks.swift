//
//  TwoBeamThicknessChecks.swift
//  Role: the thickness from a two-beam CBED disc (`TwoBeamThickness`), used from outside the package: a noiseless profile
//        of the closed form at 150 nm and an extinction distance of 90 nm gives both back, a disc line in pixels gives the
//        same deviations, Si 220's extinction distance at 200 kV is near 96.7 nm, and a short profile is refused by name.
//        Called from main.swift.
//

import Foundation
import DiffractoCore
import DiffractoGeometry

func runTwoBeamThicknessChecks(_ r: inout Report) throws {
    let t = 1500.0
    let xi = 900.0
    let edge = try unwrap(TwoBeamThickness.minimumDeviation(order: 6, thickness: t, extinctionDistance: xi))
    var s: [Double] = []
    var y: [Double] = []
    for i in 0..<201 {
        let x = -edge + 2 * edge * Double(i) / 200
        s.append(x)
        y.append(3 * TwoBeamThickness.intensity(deviation: x, thickness: t, extinctionDistance: xi) + 0.2)
    }
    let fit = try TwoBeamThickness.fit(TwoBeamProfile(deviations: s, intensities: y), extinctionDistance: 1000)
    r.near("thickness of a noiseless profile", fit.thickness, t, 1e-6 * t)
    r.near("extinction distance of a noiseless profile", fit.extinctionDistance, xi, 1e-6 * xi)
    r.equal("one answer", fit.isAmbiguous, false)
    r.equal("the fringe line is drawn", fit.fringeLine != nil, true)
    let line = try TwoBeamThickness.fringeLine(orders: [3, 4, 5], deviations: [3, 4, 5].map {
        TwoBeamThickness.minimumDeviation(order: $0, thickness: t, extinctionDistance: xi) ?? 0
    })
    r.near("fringe line through exact minima", line.thickness, t, 1e-7 * t)
    let pixels = TwoBeamProfile(offsets: [10], discSeparation: 200, gLength: 0.5, wavelength: 0.025, intensities: [1])
    r.near("deviation from a disc line in pixels", pixels.deviations[0], 10.0 / 200 * 0.025 * 0.25, 1e-18)
    let si = CrystalPresets.diamond(a: 5.431, z: 14)
    let xiSi = try TwoBeamThickness.extinctionDistance(of: MillerIndex(2, 2, 0), in: si, volts: 200e3)
    r.near("Si 220 extinction distance at 200 kV", xiSi, 967.0, 0.5)
    // 40 nm out to between its first and second minima: one minimum a side
    let s1 = try unwrap(TwoBeamThickness.minimumDeviation(order: 1, thickness: 400, extinctionDistance: 950))
    let s2 = try unwrap(TwoBeamThickness.minimumDeviation(order: 2, thickness: 400, extinctionDistance: 950))
    let half = 0.5 * (s1 + s2)
    let thin = (0..<201).map { -half + 2 * half * Double($0) / 200 }
    let thinY = thin.map { TwoBeamThickness.intensity(deviation: $0, thickness: 400, extinctionDistance: 950) }
    var refused = false
    do {
        _ = try TwoBeamThickness.fit(TwoBeamProfile(deviations: thin, intensities: thinY), extinctionDistance: 950)
    } catch TwoBeamThicknessError.tooFewMinima {
        refused = true
    } catch {}
    r.equal("too few minima refused by name", refused, true)
}
