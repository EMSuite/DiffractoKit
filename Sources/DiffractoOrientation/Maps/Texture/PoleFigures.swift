//
//  PoleFigures.swift
//  Moved from Electron Microscopy Studio's Sources/TEMSEM/Texture (at 7b09bf5) by DiffractoKit#180; the bundles its moved
//  tests read are in Tests/DiffractoOrientationTests/Fixtures/orix-sem, the others stay with the app (docs/parity.md).
//  Pole figures and inverse pole figure densities of a set of orientations, as orix 0.15.0 draws them (GPL-3.0-or-later,
//  see NOTICE): a {hkl} or {hkil} family made whole by `Miller.symmetrise(unique=True)` under the phase's point group, each
//  pole turned into the sample frame by the inverse orientation (`(~O).outer(m)`, orix's pole figure), the sample direction
//  turned into the crystal frame by the orientation (`O * v`, the inverse pole figure), and the densities of
//  `PoleDensityFunction`. The rasters the app's Maps room draws of them (`TextureRaster`) stay with the app (DiffractoKit#180).
//
//  Parity: `tools/parity/sem` bundles `m5/texture-ti`, `m5/texture-cubic`, `m5/texture-ti-full` (figure F13).
//

import Foundation

public nonisolated enum PoleFigures {
    /// `Miller.symmetrise(unique=True)` of each vector under `group`, then `.unit`: every element applied (an improper one
    /// negates), the images made unique as `Vector3d.unique` makes them (rounded to 10 decimals, first appearances kept),
    /// then unit. The families' poles in the order the vectors were given; within a family orix sorts, which no density
    /// depends on.
    public static func symmetrise(_ vectors: [SIMD3<Double>], group: PointGroup) -> [SIMD3<Double>] {
        vectors.flatMap { v in Vector3.unique(group.elements.map { $0 * v }).map(Vector3.unit) }
    }

    /// The poles of a family given as Miller indices in a lattice (three indices hkl, or four hkil), symmetrised by `group`.
    public static func poles(family: [Double], lattice: CrystalLattice, group: PointGroup) throws -> [SIMD3<Double>] {
        let m: Miller
        switch family.count {
        case 3: m = Miller(hkl: [SIMD3(family[0], family[1], family[2])], lattice: lattice)
        case 4: m = try Miller(hkil: [SIMD4(family[0], family[1], family[2], family[3])], lattice: lattice)
        default: throw Refusal.indices(family.count)
        }
        return symmetrise(m.xyz, group: group)
    }

    public enum Refusal: Swift.Error, Equatable, Sendable, CustomStringConvertible {
        case indices(Int)

        public var description: String {
            switch self {
            case .indices(let n): "A pole family takes three (hkl) or four (hkil) indices, not \(n)."
            }
        }
    }

    /// `(~O).outer(poles)`: per orientation, each pole in the sample frame; orientation-major.
    public static func sampleDirections(_ orientations: Orientations, poles: [SIMD3<Double>]) -> [SIMD3<Double>] {
        var out: [SIMD3<Double>] = []
        out.reserveCapacity(orientations.count * poles.count)
        for q in orientations.quaternions {
            let inverse = Rotation(q).inverse
            for p in poles { out.append(inverse * p) }
        }
        return out
    }

    /// `O * v`: the sample direction in each crystal frame (the inverse pole figure's vectors).
    public static func crystalDirections(_ orientations: Orientations, sampleDirection: SIMD3<Double> = Vector3.z) -> [SIMD3<Double>] {
        orientations.quaternions.map { Rotation($0) * sampleDirection }
    }

    /// The pole figure of one family: orix's `pole_density_function` of the sample directions on a hemisphere. The directions are
    /// binned as they are made, not held (N × |poles| vectors, #621): the same counts as `PoleDensityFunction.counts` of
    /// `sampleDirections`, in the same order, each bin a sum of ones.
    public static func poleFigure(_ orientations: Orientations, poles: [SIMD3<Double>], resolution: Double = 1, sigma: Double = 5,
                                   hemisphere: EqualAreaGrid.Hemisphere = .upper) -> PoleDensityMap {
        let grid = EqualAreaGrid(resolution: resolution, hemisphere: hemisphere)
        var counts = [Double](repeating: 0, count: grid.binCount)
        for q in orientations.quaternions {
            let inverse = Rotation(q).inverse
            for p in poles {
                if let k = grid.bin(of: inverse * p) { counts[k] += 1 }
            }
        }
        return PoleDensityFunction.density(counts: counts, resolution: resolution, sigma: sigma, hemisphere: hemisphere)
    }

    /// The inverse pole figure density of one sample direction in `laue`'s fundamental sector.
    public static func inversePoleFigure(_ orientations: Orientations, sampleDirection: SIMD3<Double> = Vector3.z, laue: PointGroup,
                                          resolution: Double = 1, sigma: Double = 5) -> PoleDensityMap {
        PoleDensityFunction.compute(crystalDirections(orientations, sampleDirection: sampleDirection), resolution: resolution,
                                    sigma: sigma, symmetry: laue)
    }
}
