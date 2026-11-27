//
//  OrixCheckedConversion.swift
//  Role: the total conversions the colour keys need (Electron Microscopy Studio's `CheckedConversion`, moved with them,
//        DiffractoKit#179): a NaN colour is drawn as 0 and a size that is not a number is nil, never a trap.
//

enum OrixCheckedConversion {
    /// A fraction as a byte, rounded to the nearest: below 0 is 0, above 1 is 255, NaN is 0.
    @inline(__always)
    static func unitByte(_ value: Double) -> UInt8 {
        value.isNaN ? 0 : UInt8((min(max(value, 0), 1) * 255).rounded())
    }

    /// `unitByte` for a Float.
    @inline(__always)
    static func unitByte(_ value: Float) -> UInt8 { unitByte(Double(value)) }

    /// `Int(value)`, truncation toward zero, when `value` is finite and the result fits an Int; nil otherwise.
    @inline(__always)
    static func int(_ value: Double) -> Int? {
        value.isFinite ? Int(exactly: value.rounded(.towardZero)) : nil
    }
}
