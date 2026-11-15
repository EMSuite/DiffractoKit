//
//  BestPhaseMatching.swift
//  Role: the best phase per position across several structures, with a reliability value (DiffractoKit#121, the kit side of
//        mac4DSTEM#252 "Grains and phases (no matrix)"). Every position is matched against each phase's own template library
//        (`ACOMMatcher`, `ACOMMatcher.matchAll`, unchanged); the phase whose best template scores highest wins, and the
//        reliability says how far ahead of the best OTHER phase it is.
//
//  The method is ASTAR's (Rauch & Dupuy, Arch. Metall. Mater. 50 (2005) 87; Rauch et al., Z. Kristallogr. 225 (2010) 103):
//  the correlation index of each phase's best template, and the phase reliability
//
//          R = 100 · (1 − Q2 / Q1),   Q1 the winner's index, Q2 the best index of any other phase.
//
//  ASTAR is closed software; the formula and the reading of it (R of 15 or more is taken as a sure phase) are from its published
//  descriptions, not from its source. The scores are comparable between phases because the kit's correlation is normalised on
//  both sides: the template and the measured polar image are each ring-mean-free and of unit norm (`ACOMTemplateLibrary`,
//  `ACOMMatcher.experimentalPolarImage`), as ASTAR's index divides by both norms, provided the libraries share one polar grid
//  and one set of settings (`areComparable`). A best-of-N over templates also favours the phase sampled more densely: at the
//  default 600 beam directions per library a small fundamental zone (m-3m) is sampled finer than a large one (6/mmm). py4DSTEM's correlation is not normalised on
//  the measured side, and its `CrystalPhase` decides between phases by a non-negative least-squares fit of the pattern to each
//  phase's best simulated pattern instead (`quantify_phase`, its `phase_reliability` a difference of residuals); that method
//  is not built here.
//

import Foundation

/// One position's competition between phases.
public struct BestPhaseResult: Sendable, Equatable {
    /// The winning phase (an index into the phases handed in), or -1 when no phase matched the pattern.
    public var phaseIndex: Int
    /// The winner's orientation result (its `phaseID` is `phaseIndex`); `OrientationResult.empty` when none.
    public var orientation: OrientationResult
    /// The best other phase, or -1 when no other phase scored above 0.
    public var runnerUpPhaseIndex: Int
    /// The best other phase's correlation index, 0 when none.
    public var runnerUpScore: Float
    /// ASTAR's phase reliability, 100 · (1 − runner-up / winner): 100 when no other phase scored, 0 when none matched or two
    /// phases tie. It says how far ahead of the other phases the winner is, not how well it matches: a winner at a score of
    /// 0.02 with no other phase scoring has reliability 100, so read it together with `score`.
    public var reliability: Float

    /// The winner's correlation index, 0 when none.
    public var score: Float { phaseIndex >= 0 ? orientation.score : 0 }

    public init(phaseIndex: Int, orientation: OrientationResult, runnerUpPhaseIndex: Int, runnerUpScore: Float, reliability: Float) {
        self.phaseIndex = phaseIndex; self.orientation = orientation
        self.runnerUpPhaseIndex = runnerUpPhaseIndex; self.runnerUpScore = runnerUpScore; self.reliability = reliability
    }

    /// No phase.
    public static let none = BestPhaseResult(phaseIndex: -1, orientation: .empty, runnerUpPhaseIndex: -1, runnerUpScore: 0,
                                             reliability: 0)
}

/// The best phase at every position of a scan, row-major.
public struct BestPhaseMap: Sendable, Equatable {
    public let width: Int
    public let height: Int
    public var results: [BestPhaseResult]

    public init(width: Int, height: Int, results: [BestPhaseResult]) {
        self.width = width; self.height = height; self.results = results
    }

    /// How many positions each phase won, by phase index.
    public func phaseCounts(phases: Int) -> [Int] {
        var counts = [Int](repeating: 0, count: max(0, phases))
        for r in results where counts.indices.contains(r.phaseIndex) { counts[r.phaseIndex] += 1 }
        return counts
    }
}

public enum BestPhaseMatching {
    /// The reliability above which ASTAR's descriptions take the phase as sure.
    public static let publishedReliableAbove: Float = 15

    /// One position: the phases' own best matches, in phase order (each from that phase's matcher). The highest score wins;
    /// a tie goes to the lower phase index. A phase whose match has no template (`templateIndex < 0`) or a score of 0 does not
    /// compete.
    public static func best(of matches: [OrientationResult]) -> BestPhaseResult {
        var first = -1, second = -1
        for (k, m) in matches.enumerated() where m.templateIndex >= 0 && m.score > 0 && m.score.isFinite {
            if first < 0 || m.score > matches[first].score {
                second = first; first = k
            } else if second < 0 || m.score > matches[second].score {
                second = k
            }
        }
        guard first >= 0 else { return .none }
        var winner = matches[first]
        winner.phaseID = first
        let q1 = winner.score, q2 = second >= 0 ? matches[second].score : 0
        return BestPhaseResult(phaseIndex: first, orientation: winner, runnerUpPhaseIndex: second, runnerUpScore: q2,
                               reliability: max(0, min(100, 100 * (1 - q2 / q1))))
    }

    /// Whether the libraries' scores can be compared: the measured pattern is drawn on each library's own polar grid with its
    /// own settings (``ACOMMatcher/experimentalPolarImage(peaks:originX:originY:invAngstromPerPixel:)``), so two libraries
    /// with another reach, grid or kernel score one pattern against two different images. Same geometry and same settings.
    public static func areComparable(_ libraries: [ACOMTemplateLibrary]) -> Bool {
        guard let first = libraries.first else { return false }
        return libraries.allSatisfy { $0.geometry == first.geometry && $0.settings == first.settings }
    }

    /// A scan: each phase's results over the same positions (`ACOMMatcher.matchAll` once per phase, in phase order), and the
    /// libraries they were matched with, in the same order.
    /// - Returns: nil when there is no phase, the libraries are not as many as the scans or not comparable
    ///   (``areComparable(_:)``), or the scans differ in size or do not fill their grid.
    public static func bestPhases(_ scans: [ACOMScanResults], libraries: [ACOMTemplateLibrary]) -> BestPhaseMap? {
        guard let first = scans.first, libraries.count == scans.count, areComparable(libraries),
              first.results.count == first.width * first.height,
              scans.allSatisfy({ $0.width == first.width && $0.height == first.height && $0.results.count == first.results.count })
        else { return nil }
        let results = first.results.indices.map { i in best(of: scans.map { $0.results[i] }) }
        return BestPhaseMap(width: first.width, height: first.height, results: results)
    }
}
