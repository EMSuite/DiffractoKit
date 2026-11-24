//
//  IntegerMath.swift
//  The kit's one greatest common divisor (#226: six private copies, some taking the magnitude and some not).
//

package enum IntegerMath {
    /// The greatest common divisor of |a| and |b|, never negative; `gcd(0, 0)` is 0 (callers that divide floor the result).
    package static func gcd(_ a: Int, _ b: Int) -> Int {
        var a = abs(a), b = abs(b)
        while b != 0 { (a, b) = (b, a % b) }
        return a
    }
}
