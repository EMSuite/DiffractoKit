//
//  PrecipitateFixtures.swift
//  Role: the planted precipitates of the variant tests (DiffractoKit S7, #162) and their relationships in aluminium, the
//        same numbers as tools/pins/orix/variant-pins.py. Planted cells: the lattices are the published ones, the sites are
//        chosen only to give the published Laue class (the variants depend on the cell and the class, not on the atoms).
//

import Foundation
import DiffractoCore

package enum PrecipitateFixtures {
    /// θ′ (Al2Cu) cell, a 4.04, c 5.80 Å (Silcock, Heal and Hardy 1953): Cu body-centred, Al at (0 ½ ¼) and its images; 4/mmm.
    package static var thetaPrime: CrystalStructure {
        CrystalPresets.structure(a: 4.04, b: 4.04, c: 5.80, sites: [
            AtomSite(z: 29, fractional: [0, 0, 0]), AtomSite(z: 29, fractional: [0.5, 0.5, 0.5]),
            AtomSite(z: 13, fractional: [0, 0.5, 0.25]), AtomSite(z: 13, fractional: [0.5, 0, 0.25]),
            AtomSite(z: 13, fractional: [0, 0.5, 0.75]), AtomSite(z: 13, fractional: [0.5, 0, 0.75]),
        ], name: "theta' (planted)")
    }

    /// T1 (Al2CuLi) cell, a 4.965, c 9.345 Å (Huang and Ardell 1987), sites of the AlB2 type plus Li on the axis; 6/mmm.
    package static var t1: CrystalStructure {
        CrystalPresets.structure(a: 4.965, b: 4.965, c: 9.345, gammaDeg: 120, sites: [
            AtomSite(z: 29, fractional: [0, 0, 0]),
            AtomSite(z: 13, fractional: [1.0 / 3.0, 2.0 / 3.0, 0.5]), AtomSite(z: 13, fractional: [2.0 / 3.0, 1.0 / 3.0, 0.5]),
            AtomSite(z: 3, fractional: [0, 0, 0.5]),
        ], name: "T1 (planted)")
    }
}
