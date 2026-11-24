# Status: DiffractoKit 0.21.0 <!-- x-release-please-version -->

DiffractoKit 0.21.0. <!-- x-release-please-version -->
It builds on ImagoKit `"0.2.0"..<"1.0.0"`.

Live:

- Cells in the seven systems, Niggli reduction, the conventional cell.
- CIF import with symmetry detection and a completeness check by IT number, by symbol or by Hall symbol (a CIF with no
  operation loop is expanded from its Hall symbol, Hermann–Mauguin symbol or IT number, in pymatgen's setting, with a warning; an operation list that is not its space group is refused, and one that disagrees with the
  Hall symbol warns); space-group-free symmetry from the atoms; reflections and rings.
- The space-group type the atoms show (`CrystalStructure.spaceGroupType`: IT number, spglib's short symbol and standard
  Hall number, the setting written in), equal to spglib on all 530 settings; the number a source states is kept beside it.
  A cell whose metric is more symmetric than its atoms (Amm2 written on a cube) is not named yet.
- The CIFs people have: `?` and `.` values, deuterium and tritium, dummy atoms (left out with a warning), U and B iso,
  Latin-1 and a byte-order mark, DDLm names, a publication block before the structure; mmCIF and mCIF refused by name,
  a number that is not finite refused with a typed error, and a seeded fuzz of the committed CIFs that never crashes (#117).
- A cell that is singular, or whose volume, occupancy or metric is not a number, is refused by `UnitCell`,
  `StructureFactors.compute` and the symmetry's metric test instead of computing on an identity inverse or NaN (#147).
- A structure as stored data: `CrystalStructure.Snapshot` (Codable, with a schema version) rebuilds the same
  structure, and `CrystalStructure.fingerprint` changes with any cell parameter, any
  site's element, position, occupancy or displacement, or the symmetry operations (#148).
- A Materials Project structure with no key: `MaterialsProjectOptimade` searches the keyless OPTIMADE service by
  elements or formula (most stable first; an empty service is told apart from no results) and reads an entry to its
  conventional cell, checked against pymatgen on 14 recorded entries, one per Bravais lattice (#150).
- A Materials Project search with the user's own key: `MaterialsProjectImport.searchRequest` asks the main API by
  elements, chemical system or formula, most stable first, deprecated and GNoME entries left out; `decodeSearch` reads
  the rows (stability, calculated-only and deprecated flags, the stated space group) and `databaseVersion(from:)` the
  database version. The kit holds no key; the answers recorded with the owner's key are still to come (#151).
- Lattice edits within the crystal system: `CrystalStructure.editingLattice` sets the free parameters (cubic a;
  tetragonal and hexagonal a, c; trigonal a, c in hexagonal axes and a, α in rhombohedral axes; orthorhombic a, b, c;
  monoclinic a, b, c, β; triclinic all six), keeps the fractional sites and refuses an edit that changes the Laue class
  (#149).
- A CIF of several structures: `CIFImport.blocks(in:)` lists its `data_` blocks and `read(_:block:)` reads one by name (#157).
- Kinematic structure factors (py4DSTEM's scattering-factor table, bit for bit).
- Powder X-ray patterns (`PowderPattern.xray`): pymatgen's `XRDCalculator` peaks, 2θ, intensities and hkl labels for any
  structure, Cu Kα and pymatgen's other lines or any wavelength, a 2θ range and Debye-Waller factors; its X-ray scattering
  table bit for bit (`XRayScatteringFactors`).
- Electron powder rings (`PowderPattern.electron`) from the held ring list with relative intensities, a drawn profile
  of any powder pattern (`profile(width:range:step:)`), and the overlapped and diagnostic peaks of several shown phases
  (`PowderPattern.marks`).
- Bonds for a crystal view by VESTA's cut-offs, as pymatgen's `CutOffDictNN.from_preset("vesta_2019")` gives them
  (`CrystalStructure.bonds`): metals show none; the table is generated from pymatgen's file.
- Pattern geometry: ring and ellipse fits, calibration, zone axes, indexing; the ring-sequence assignment of a known crystal,
  Δq from several rings, and the single-pattern beam centre (`PatternCentre.find` refuses a blocked beam; `FriedelCentre`
  finds the centre of a pattern with one), radial profile and spot refinement.
- The spots of one pattern (`PatternSpots`): py4DSTEM's `get_maxima_2D` exactly, then a radial background, the beam left
  out and each spot refined.
- The spot pattern of a structure at any orientation and voltage (`SpotPattern`): py4DSTEM's kinematic pattern, every spot
  with hkl, position, d, excitation error and intensity, and its distance R = L tan 2θ at a camera length.
- The zone axis and full orientation of one measured spot pattern (`ZoneAxisSolver`): lattice-vector pairs scored on every
  spot against the kinematic pattern, the scale and in-plane turn fitted, every spot's hkl, and a verdict (solved, a tie,
  or refused for one row of spots or a poor fit); the half turn about the beam is flagged where the crystal lacks it.
- The tilt planner's kit half: the stereographic projection (orix's), the zone-axis map of an orientation (`ZoneAxisMap`:
  poles to index 3 and the traces of the lowest plane families with their Bragg angles) and a double-tilt holder
  (`DoubleTiltHolder`): the α and β that bring a zone axis onto the beam, whether the holder reaches them, and the tilt
  axis calibrated from two indexed patterns at two tilts, which also settles the half turn of one spot pattern.
- The spots double diffraction adds to one crystal's pattern or between crystals (a precipitate in its matrix), each with
  its routes (`DoubleDiffraction`).
- Double diffraction between a precipitate and its matrix on the two drawn patterns, each position with its routes,
  weighted by the relative intensities and cut at 0.1 (`SpotPattern.doubleDiffraction(with:routeWeightCut:)`).
- The tilt of an indexed pattern off its zone axis, refined on the spot intensities (`ZoneAxisSolver.refineTilts`), with
  the look-alike tilts that agree as well returned beside the best.
- The full orientation of one spot pattern (`PatternOrientation.index`): template matches and zone-axis solutions each
  refined on positions and intensities, the best agreement with the look-alikes that agree as well, and a flag when the
  pattern cannot decide (an alternative, or a tilt the intensities fix poorly).
- The thickness from the profile across a two-beam CBED disc (`TwoBeamThickness`): the whole profile fitted by the
  two-beam intensity with the extinction distance inside a window about the phase's structure-factor value, the fringe
  line as the cross-check, and both thicknesses when two fit equally; profiles with fewer than two minima a side are
  refused. Checked on closed forms and an abTEM-simulated Si 220 disc (within 3 %), not yet on a measured pattern.
- The ACOM matcher (CPU reference, Metal option, optional in-plane refinement) and vector-matched phase maps.
- The best phase per position across several structures, with ASTAR's phase reliability (`BestPhaseMatching`).
- One orientation value, `CrystalOrientation` (DiffractoCore): py4DSTEM's matrix from a zone axis and an in-plane angle,
  tilts and turns in the lab frame, the nearest zone axis and the angle to it, orix quaternions in and out.
- Orientation symmetry and IPF keys, misorientation between two groups, `Orientation` from Euler angles with `dot`,
  `dotOuter` and `distanceMatrix`, orientation-relationship variants (beta'' in Al: 12), orientation-map planes.
- Variant orientations of a relationship stated as two pairs (`OrientationRelationship.variantOrientations`): each
  variant's orientation, label, drawn spot pattern and group of look-alikes in the view, and the second pair's misfit
  (beta'', theta' and T1 in Al: 12, 3 and 4).
- The phase finder (`PhaseFinder`): candidate phases ranked on one ring or spot pattern by the numbers ring indexing and
  the zone-axis solver already compute, with ties where the pattern cannot tell phases apart; ranked evidence, never a
  verdict.

Since 0.10.0: 0.11.0 (2026-10-30) dropped the unused public API; 0.12.0 (2026-11-01) added the items above, made the ACOM
Metal path 1.9 to 3.2 times faster and the ring profile fit stop runs that cannot improve, and retired the apps-parity dumps.
Details and the migration notes: `CHANGELOG.md`.

## In progress, not wired: phase orientations from the structure alone (#81)

The aim: a crystal structure in, and the scan's phases and their orientations out, with no zone axis typed. Five steps are
merged as kit functions; none is used by `PhaseVectorMatcher.map`. State of 2026-11-11.

| Step | PR | What it is |
|---|---|---|
| 1 | #82 | The matrix orientation of every position, its spots predicted and removed (`matrixOrientations`, `removingMatrixSpots`) |
| 2 | #84 | One orientation per grain (`regularizingMatrixOrientations`) and the guarded free search of the residual (`searchOrientations`) |
| 3 | #85 | Relationships learned from the strong positions, refined, and every position scored against their variants (`learnOrientationRelationships`, `classifyAgainstVariants`); grain overlaps in the removal |
| 4 | #86 | The two passes as one call, scored on the whole pattern, firm rule by default, neighbour rule as a switch (`classifyInTwoPasses`) |
| 5 | #88 | The rational relationships near a learned one, listed and never chosen (`rationalRelationships`) |
| guard | #95 | A displaced matrix pattern is never a phase (inside `classifyInTwoPasses`) |
| input | #96 | The detector's distortion measured on the matrix's spots, for a scan with no ring (`matrixSpotDistortion`) |
| guard | #103 | A neighbouring grain's spots are no phase's vectors: a grain boundary is not called a precipitate (inside `classifyInTwoPasses`) |
| report | #107 | Which variants predict the same spots within the reach and cannot be told apart (`indistinguishableVariants`) |
| one call | #99 | Steps 1 to 4 chained with the values they share, returning a `PhaseMap` (`mapFromStructures`); the call an application makes |

Where it stands against the bar (typed zone axes through `map`):

- Thronsen et al. dataset A, stride 3: 28 771 of 29 241 positions (98.39 %) with 103 matrix positions called a phase,
  against 28 858 (98.69 %) and 104. With the neighbour rule 98.51 % and 112. The positions that decide it show one vector of
  their own: the precipitate edges.
- The synthetic three-grain demo: relationship 0.03 degrees from the planted one, all nine precipitates found, nothing
  labelled outside them.
- Planted pairs (`TwoPassOtherMaterialsTests`):

| Pair | Relationship learned | Positions |
|---|---|---|
| alpha-Ti (hcp) + beta-Ti, Burgers | 2.65 degrees off, 12 variants for 6 | all right |
| alpha-Fe (bcc) + cementite, Bagaryatsky | 1.44 degrees off | all right |
| Al + beta'' (monoclinic), Andersen | exact | all right |
| Ni + gamma' Ni3Al, cube on cube | nothing learned | no gamma' found (known gap, pinned by the test) |

It is wired only when Thronsen A reaches the typed result on both numbers and nothing is worse on the demo or the pairs.

Next step: lattice-correspondence proposals in `searchOrientations`. For each position, besides the template proposals,
propose the precipitate orientations that put low-index precipitate lattice vectors onto matrix lattice vectors at that
position's matrix orientation (the free search does not take the matrix orientation today), and score them on the residual.
A phase on the matrix lattice (gamma') leaves only its superlattice spots after the removal and the template match never
proposes its orientation. The case to turn round is Ni + gamma' in `TwoPassOtherMaterialsTests`. Then before and after on
Thronsen A, the demo and the four pairs.

Parked: generating `map`'s entries (zone axis per variant, in-plane pair, slab) from a relationship. The typed Thronsen
entries carry hand-set values (the T1 zone [0 -4 1], a slab of 0.3 for theta'), and reproducing them by rule risks tuning
to that dataset.

Open beside it: the tolerances in Å⁻¹ should follow the disk size and the calibration (#87). Measurements, per step:
`docs/parity.md`.

## Tests

`DiffractoKitParityTests` holds the pinned comparisons against the Python references; the other test targets hold the
fast unit tests; `DiffractoKitTestSupport` holds the helpers both use.

What is checked against what: `docs/parity.md`. Timings against Python: `docs/benchmarks.md`. What is open: the GitHub
issues; known limits: `docs/parity.md`. Where it is going: `ROADMAP.md`.
