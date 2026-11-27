//
//  CrystalMapPlainArrayTests.swift
//  The public surface the map code gained with its move from Electron Microscopy Studio (DiffractoKit#180), through a plain
//  `import` as a consumer sees it: a crystal map built from a reader's plain arrays (Euler angles or quaternions, phase ids, x, y,
//  named properties and the phase headers as plain values, so a file reader never names `CrystalMap`) is the map the `PhaseList`
//  initialisers build; a phase's lattice as the kit's `CrystalLattice`; the slip analysis of one orientation as the kit's
//  `Rotation`. Red before: none of `init(euler:quaternions:...:phaseHeaders:...)`, `CrystalPhase.crystalLattice` or the `Rotation`
//  forms of `SlipTraceAnalysis` exists.
//

import DiffractoOrientation
import Foundation
import Testing

@Suite struct CrystalMapPlainArrayTests {
    /// A 2 × 3 map of two phases and one point not indexed, on a 0.5 µm grid.
    private static let euler: [SIMD3<Double>] = (0..<6).map { SIMD3(0.1 * Double($0), 0.2, 0.3) }
    private static let ids = [1, 1, 2, 2, -1, 1]
    private static let grid = CrystalMap.coordinates(shape: (rows: 2, columns: 3), steps: (dy: 0.5, dx: 0.5))
    private static let iq: [String: [Double]] = ["iq": [10, 11, 12, 13, 14, 15]]
    private static let headers: [(id: Int, name: String, formula: String?, spaceGroup: Int?, pointGroup: String?, lattice: [Double]?,
                                  latticeUnit: String?, color: String?)] = [
        (id: 1, name: "austenite", formula: "Fe", spaceGroup: 225, pointGroup: "m-3m", lattice: [3.6, 3.6, 3.6, 90, 90, 90],
         latticeUnit: "angstrom", color: nil),
        (id: 2, name: "ferrite", formula: nil, spaceGroup: 229, pointGroup: "m-3m", lattice: [2.87, 2.87, 2.87, 90, 90, 90],
         latticeUnit: nil, color: nil),
    ]

    /// The same headers as `CrystalPhase`s in a `PhaseList`, as the app's readers build them.
    private static var phaseList: PhaseList {
        PhaseList([
            (id: 1, phase: CrystalPhase(name: "austenite", formula: "Fe", spaceGroup: 225, pointGroup: "m-3m",
                                        lattice: [3.6, 3.6, 3.6, 90, 90, 90], latticeUnit: .angstrom)),
            (id: 2, phase: CrystalPhase(name: "ferrite", spaceGroup: 229, pointGroup: "m-3m", lattice: [2.87, 2.87, 2.87, 90, 90, 90])),
        ])
    }

    @Test func eulerAnglesAndPhaseHeadersMakeThePhaseListsMap() throws {
        let plain = try CrystalMap(euler: Self.euler, phaseIDs: Self.ids, x: Self.grid.x, y: Self.grid.y, properties: Self.iq,
                                   phaseHeaders: Self.headers, scanUnit: "um")
        let listed = try CrystalMap(euler: Self.euler, phaseIDs: Self.ids, x: Self.grid.x, y: Self.grid.y, phaseList: Self.phaseList,
                                    properties: Self.iq, scanUnit: "um")
        #expect(plain == listed)
        #expect(plain.phases.ids == [-1, 1, 2])
        #expect(plain.phases.names == ["not_indexed", "austenite", "ferrite"])
        // A header without a colour takes orix's default, then the next free Tableau colour.
        #expect(plain.phases.colors == ["white", "tab:blue", "tab:orange"])
        #expect(plain.phases[id: 1]?.latticeUnit == .angstrom)
        #expect(plain.phases[id: 2]?.latticeUnit == nil)
        #expect(plain.shape == [2, 3])
        #expect(plain.mapData(named: "iq") == [10, 11, 12, 13, 14, 15])
    }

    @Test func quaternionsMakeTheRotationsMap() throws {
        let q = Self.euler.map { RotationConversions.eu2qu($0) }
        let plain = try CrystalMap(quaternions: q, phaseIDs: Self.ids, x: Self.grid.x, y: Self.grid.y, phaseHeaders: Self.headers)
        let listed = try CrystalMap(rotations: q, phaseIDs: Self.ids, x: Self.grid.x, y: Self.grid.y, phaseList: Self.phaseList)
        #expect(plain == listed)
    }

    @Test func bothOrNeitherOrientationArrayIsRefused() {
        let q = Self.euler.map { RotationConversions.eu2qu($0) }
        #expect(throws: CrystalMap.Refusal.orientations) {
            try CrystalMap(euler: Self.euler, quaternions: q, phaseHeaders: Self.headers)
        }
        #expect(throws: CrystalMap.Refusal.orientations) {
            try CrystalMap(phaseHeaders: Self.headers)
        }
    }

    @Test func aPhasesLatticeIsTheKitsCrystalLattice() {
        let ti = CrystalPhase(name: "Ti", pointGroup: "622", lattice: [2.951, 2.951, 4.684, 90, 90, 120])
        #expect(ti.crystalLattice == CrystalLattice(a: 2.951, b: 2.951, c: 4.684, alpha: 90, beta: 90, gamma: 120))
        let none = CrystalPhase(name: "x", pointGroup: "6/mmm")
        #expect(none.crystalLattice == CrystalLattice.defaultLattice(system: PointGroup.named("6/mmm")?.system))
    }

    /// The cube orientation under a load along [001]: (111)[01-1] has |cos φ cos λ| = 1/√6, and (111) cuts the map plane 45° from +y.
    @Test func theSlipAnalysisTakesTheKitsRotation() throws {
        let fcc = try SlipSystemSet.load(.cubicFCC, crystal: .cubic)
        let sf = SlipTraceAnalysis.schmidFactors(orientation: Rotation.identity, load: SIMD3(0, 0, 1), systems: fcc)
        #expect(abs(sf[0][0] - 1 / (3.0.squareRoot() * 2.0.squareRoot())) <= 1e-15)
        #expect(SlipTraceAnalysis.maximumSchmidFactor(orientation: Rotation.identity, load: SIMD3(0, 0, 1), systems: fcc)
            == sf.joined().max())
        let t = SlipTraceAnalysis.traces(orientation: Rotation.identity, systems: fcc)
        #expect(abs(t[0].angle - Double.pi / 4) <= 1e-15)
    }
}
