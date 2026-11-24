//
//  Bonds.swift
//  Role: the bonds a crystal view draws, by VESTA's rule as pymatgen holds it (CutOffDictNN.from_preset("vesta_2019"),
//        pymatgen-core 2026.10.16): two atoms are bonded when their distance is strictly below the cut-off of their
//        element pair; a pair the table does not list has no bond, so metals and most intermetallics show none. The
//        view draws what this returns and computes no distance itself.
//

import simd

/// One bond seen from one of its atoms: `site` is bonded to the atom of `neighbour` translated by `image` lattice
/// vectors, at `sites[neighbour].fractional + image` (fractional, as the coordinates are given: an image is relative
/// to the stored coordinate, not to a wrapped one, as pymatgen's `jimage`).
public struct Bond: Sendable, Hashable {
    /// Index of the atom in `CrystalStructure.sites`.
    public let site: Int
    /// Index of the bonded atom in `CrystalStructure.sites`.
    public let neighbour: Int
    /// The lattice translation of the bonded atom.
    public let image: SIMD3<Int>
    /// Centre-to-centre distance, Å.
    public let length: Double

    public init(site: Int, neighbour: Int, image: SIMD3<Int>, length: Double) {
        self.site = site
        self.neighbour = neighbour
        self.image = image
        self.length = length
    }

    /// The same bond seen from the other atom.
    public var reversed: Bond { Bond(site: neighbour, neighbour: site, image: SIMD3<Int>.zero &- image, length: length) }
}

/// VESTA's bond cut-offs per element pair (Å), from the pinned pymatgen file.
public enum BondCutoffs {
    /// The cut-off of the pair of atomic numbers `z1` and `z2` (either order), or nil when the table lists no bond for
    /// the pair. 905 pairs: the file's 912 less its seven deuterium and "NH" rows (a site here holds only an atomic
    /// number, and the kit reads D as H).
    public static func vesta(_ z1: Int, _ z2: Int) -> Double? {
        table[key(z1, z2)]
    }

    private static func key(_ z1: Int, _ z2: Int) -> Int { min(z1, z2) &* 1024 &+ max(z1, z2) }

    private static let table: [Int: Double] = {
        var table: [Int: Double] = [:]
        table.reserveCapacity(VestaBondCutoffs.pairs.count)
        for (z1, z2, cutoff) in VestaBondCutoffs.pairs { table[key(z1, z2)] = cutoff }
        return table
    }()
}

extension CrystalStructure {
    /// Every bond of every atom by VESTA's cut-offs (`BondCutoffs.vesta`): each bond appears twice, once from each
    /// atom (`reversed`), so the bonds of one atom are its coordination. Sorted by site, then neighbour, then image.
    /// Two atoms of one element within 1e-8 Å are one atom and never bonded (pymatgen's coincidence rule, which also
    /// covers an atom and itself); two elements at one position (a mixed site written as two sites) are a zero-length
    /// bond from both atoms when their pair is listed (pymatgen lists it from the first site of the pair only, an
    /// artefact of its neighbour search: docs/parity.md). Occupancy is not read: a partly occupied
    /// site bonds by its element as VESTA draws it (pymatgen, which keys by species string such as "Si:0.5", gives
    /// such a site none). Empty for a cell whose volume is not finite and positive.
    public var bonds: [Bond] { bonds(cutoff: BondCutoffs.vesta) }

    /// Every bond by `cutoff` (Å per pair of atomic numbers, nil for no bond), each from both atoms as `bonds`.
    public func bonds(cutoff: (Int, Int) -> Double?) -> [Bond] {
        let lattice = cell.realLattice
        guard !sites.isEmpty, cell.volume.isFinite, cell.volume > 0,
              lattice.allSatisfy({ $0.x.isFinite && $0.y.isFinite && $0.z.isFinite })
        else { return [] }
        // The spacing of the lattice planes normal to each axis: a translation by n along axis k moves an atom by at
        // least |n + Δf_k| · height_k, so the images within a cut-off r span |Δf_k + n| ≤ r / height_k.
        let heights = [
            cell.volume / length(cross(lattice[1], lattice[2])),
            cell.volume / length(cross(lattice[2], lattice[0])),
            cell.volume / length(cross(lattice[0], lattice[1])),
        ]
        var result: [Bond] = []
        for i in sites.indices {
            for j in sites.indices {
                guard let r = cutoff(sites[i].z, sites[j].z), r.isFinite, r > 0 else { continue }
                let delta = sites[j].fractional - sites[i].fractional
                guard delta.x.isFinite, delta.y.isFinite, delta.z.isFinite else { continue }
                let bounds = (0..<3).map { axis -> (Double, Double) in
                    let reach = r / heights[axis]
                    return ((-reach - delta[axis]).rounded(.down), (reach - delta[axis]).rounded(.up))
                }
                // A coordinate millions of cells away is no structure to draw; it would only overflow the range.
                guard bounds.allSatisfy({ $0.0 > -1e6 && $0.1 < 1e6 }) else { continue }
                let ranges = bounds.map { Int($0.0)...Int($0.1) }
                for nx in ranges[0] {
                    for ny in ranges[1] {
                        for nz in ranges[2] {
                            let image = SIMD3<Int>(nx, ny, nz)
                            let f = delta + SIMD3<Double>(Double(nx), Double(ny), Double(nz))
                            let cartesian = f.x * lattice[0] + f.y * lattice[1] + f.z * lattice[2]
                            let distance = length(cartesian)
                            if distance <= 1e-8, sites[i].z == sites[j].z { continue }
                            if distance < r {
                                result.append(Bond(site: i, neighbour: j, image: image, length: distance))
                            }
                        }
                    }
                }
            }
        }
        return result
    }
}
