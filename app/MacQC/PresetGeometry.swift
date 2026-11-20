//
//  PresetGeometry.swift
//  ΨLab 15.0
//
//  The geometry parameters of the molecules of the menu that came with 15.0 (H₂⁺, H₃⁺, H₅⁺, N₂, O₂,
//  F₂, Cl₂, NO, NO₂, CO, CO₂, N₂O, H₂CO, CH₄), as the sidebar shows them for H₂, LiH, H₂O, H₂O₂ and
//  NH₃: bond lengths in Å and angles in degrees, in the molecule's symmetry. The defaults are
//  measured from the backend's coordinates (molecule.cpp, copied here), so unchanged parameters
//  run the backend's own geometry; changed ones are sent as Cartesian coordinates (Å) in the
//  backend's atom order.
//

import Foundation

/// A geometry parameter: a bond length (Å) or an angle (degrees).
nonisolated struct PresetParameter: Sendable, Hashable {
    let label: String       // e.g. "R(N–O)", "Angle O–N–O"
    let angular: Bool

    var unit: String { angular ? "°" : "Å" }

    /// The range an input is kept in.
    var range: ClosedRange<Double> { angular ? 30.0...180.0 : 0.3...6.0 }
}

nonisolated enum PresetGeometry {
    /// The backend's coordinates (Å) of the molecules with parameters (atoms in the backend's order).
    private static let reference: [String: [[Double]]] = {
        let b: Double = 0.529177210544  // Å per bohr (molecule.h)
        func pair(_ d: Double) -> [[Double]] { [[0.0, 0.0, -0.5 * d], [0.0, 0.0, 0.5 * d]] }
        var r: [String: [[Double]]] = [:]
        r["H2+"] = pair(1.9972 * b)
        let rc: Double = 1.6500 * b / 3.0.squareRoot()
        var h3: [[Double]] = []
        for k in 0..<3 {
            let t: Double = (90.0 + 120.0 * Double(k)) * Double.pi / 180.0
            h3.append([rc * cos(t), rc * sin(t), 0.0])
        }
        r["H3+"] = h3
        r["H5+"] = [[0.0, 0.0, -0.169235], [0.0, 0.405338, -1.050020], [0.0, -0.405338, -1.050020],
                    [0.380706, 0.0, 1.134638], [-0.380706, 0.0, 1.134638]]
        r["N2"] = pair(1.09768)
        r["O2"] = pair(1.20752)
        r["F2"] = pair(1.41193)
        r["Cl2"] = pair(1.9879)
        r["NO"] = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.15077]]
        let half: Double = 0.5 * 134.1 * Double.pi / 180.0
        let sx: Double = 1.1934 * sin(half), sz: Double = -1.1934 * cos(half)
        r["NO2"] = [[0.0, 0.0, 0.0], [sx, 0.0, sz], [-sx, 0.0, sz]]
        r["CO"] = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.1283]]
        r["CO2"] = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.16], [0.0, 0.0, -1.16]]
        r["N2O"] = [[0.0, 0.0, -1.1273], [0.0, 0.0, 0.0], [0.0, 0.0, 1.1851]]
        r["H2CO"] = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.2033], [0.0, 0.934192, -0.581709], [0.0, -0.934192, -0.581709]]
        let c: Double = 0.62758
        r["CH4"] = [[0.0, 0.0, 0.0], [c, c, c], [-c, -c, c], [-c, c, -c], [c, -c, -c]]
        return r
    }()

    /// The parameters of a molecule of the menu (none for the atoms).
    static func parameters(_ key: String) -> [PresetParameter] {
        func r(_ s: String) -> PresetParameter { PresetParameter(label: s, angular: false) }
        func a(_ s: String) -> PresetParameter { PresetParameter(label: s, angular: true) }
        switch key {
        case "H2+", "N2", "O2", "F2", "Cl2": return [r("R")]
        case "NO": return [r("R(N–O)")]
        case "CO": return [r("R(C–O)")]
        case "H3+": return [r("R(H–H)")]
        case "H5+": return [r("R(H1–H2), R(H1–H3)"), r("R(H2–H3)"), r("R(H4–H5)"), r("R(H1–X)")]
        case "NO2": return [r("R(N–O)"), a("Angle O–N–O")]
        case "CO2": return [r("R(C–O)")]
        case "N2O": return [r("R(N–N)"), r("R(N–O)")]
        case "H2CO": return [r("R(C–O)"), r("R(C–H)"), a("Angle H–C–H")]
        case "CH4": return [r("R(C–H)")]
        default: return []
        }
    }

    /// The symmetry the parameters keep (as the notes of H₂O and H₂O₂).
    static func note(_ key: String) -> String? {
        switch key {
        case "H3+": return "D₃ₕ: an equilateral triangle."
        case "H5+": return "C₂ᵥ: an H₃⁺ unit H1 H2 H3 (H2–H3 at the base) and an H₂ unit H4 H5 perpendicular to it; X is the midpoint of H4–H5, on the axis through H1."
        case "NO2": return "C₂ᵥ: both N–O bonds equal."
        case "CO2": return "D∞ₕ: linear, both C–O bonds equal."
        case "N2O": return "C∞ᵥ: linear N–N–O."
        case "H2CO": return "C₂ᵥ: both C–H bonds equal."
        case "CH4": return "T_d: four equal C–H bonds, tetrahedral angles."
        default: return nil
        }
    }

    /// The parameters of coordinates (Å) in the backend's atom order (symmetric pairs averaged), or
    /// nil if the atoms do not fit the molecule.
    static func measure(_ key: String, _ x: [[Double]]) -> [Double]? {
        let n = parameters(key).count
        guard n > 0, x.count == (reference[key]?.count ?? -1), x.allSatisfy({ $0.count == 3 }) else { return nil }
        func d(_ i: Int, _ j: Int) -> Double {
            var s: Double = 0
            for c in 0..<3 {
                let t: Double = x[i][c] - x[j][c]
                s += t * t
            }
            return s.squareRoot()
        }
        func angle(_ i: Int, _ j: Int, _ k: Int) -> Double {
            let u = (0..<3).map { x[i][$0] - x[j][$0] }, v = (0..<3).map { x[k][$0] - x[j][$0] }
            let c = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (d(i, j) * d(k, j))
            return acos(min(max(c, -1), 1)) * 180 / Double.pi
        }
        switch key {
        case "H2+", "N2", "O2", "F2", "Cl2", "NO", "CO": return [d(0, 1)]
        case "H3+": return [(d(0, 1) + d(1, 2) + d(0, 2)) / 3]
        case "H5+":
            let m = (0..<3).map { 0.5 * (x[3][$0] + x[4][$0]) }
            let dx = (0..<3).map { x[0][$0] - m[$0] }
            return [0.5 * (d(0, 1) + d(0, 2)), d(1, 2), d(3, 4), (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).squareRoot()]
        case "NO2": return [0.5 * (d(0, 1) + d(0, 2)), angle(1, 0, 2)]
        case "CO2": return [0.5 * (d(0, 1) + d(0, 2))]
        case "N2O": return [d(0, 1), d(1, 2)]
        case "H2CO": return [d(0, 1), 0.5 * (d(0, 2) + d(0, 3)), angle(2, 0, 3)]
        case "CH4": return [(d(0, 1) + d(0, 2) + d(0, 3) + d(0, 4)) / 4]
        default: return nil
        }
    }

    /// The default parameters (the backend's geometry).
    static func defaults(_ key: String) -> [Double] {
        guard let x = reference[key] else { return [] }
        return measure(key, x) ?? []
    }

    /// Cartesian coordinates (Å, the backend's atom order) of parameters, or nil if they do not
    /// give a geometry (H₅⁺: R(H1–H2) must exceed half of R(H2–H3)).
    static func coordinates(_ key: String, _ v: [Double]) -> [[Double]]? {
        guard v.count == parameters(key).count, !v.isEmpty, v.allSatisfy({ $0 > 0 && $0.isFinite }) else { return nil }
        let deg: Double = Double.pi / 180.0
        let o: [Double] = [0.0, 0.0, 0.0]
        switch key {
        case "H2+", "N2", "O2", "F2", "Cl2":
            let z: Double = 0.5 * v[0]
            return [[0.0, 0.0, -z], [0.0, 0.0, z]]
        case "NO", "CO":
            return [o, [0.0, 0.0, v[0]]]
        case "H3+":
            let rc: Double = v[0] / 3.0.squareRoot()
            var x: [[Double]] = []
            for k in 0..<3 {
                let t: Double = (90.0 + 120.0 * Double(k)) * deg
                x.append([rc * cos(t), rc * sin(t), 0.0])
            }
            return x
        case "H5+":
            let h: Double = v[0] * v[0] - 0.25 * v[1] * v[1]
            guard h > 1e-6 else { return nil }
            let z: Double = -h.squareRoot(), y: Double = 0.5 * v[1], x: Double = 0.5 * v[2]
            return [o, [0.0, y, z], [0.0, -y, z], [x, 0.0, v[3]], [-x, 0.0, v[3]]]
        case "NO2":
            let half: Double = 0.5 * v[1] * deg
            let x: Double = v[0] * sin(half), z: Double = -v[0] * cos(half)
            return [o, [x, 0.0, z], [-x, 0.0, z]]
        case "CO2":
            return [o, [0.0, 0.0, v[0]], [0.0, 0.0, -v[0]]]
        case "N2O":
            return [[0.0, 0.0, -v[0]], o, [0.0, 0.0, v[1]]]
        case "H2CO":
            let half: Double = 0.5 * v[2] * deg
            let y: Double = v[1] * sin(half), z: Double = -v[1] * cos(half)
            return [o, [0.0, 0.0, v[0]], [0.0, y, z], [0.0, -y, z]]
        case "CH4":
            let c: Double = v[0] / 3.0.squareRoot()
            return [o, [c, c, c], [-c, -c, c], [-c, c, -c], [c, -c, -c]]
        default:
            return nil
        }
    }

    /// The parameters with their units for a text line ("R(N–O) = 1.1934 Å, angle O–N–O = 134.10°").
    static func describe(_ key: String, _ v: [Double]) -> String {
        zip(parameters(key), v).map { p, x in
            p.angular ? String(format: "%@ = %.2f°", p.label, x) : String(format: "%@ = %.4f Å", p.label, x)
        }.joined(separator: ", ")
    }
}

// MARK: - The parameters of the selected molecule of the menu

extension AppModel {
    /// The geometry parameters of the selected molecule of the menu (the defaults until changed).
    var presetParameterValues: [Double] {
        presetValues[presetKey] ?? PresetGeometry.defaults(presetKey)
    }

    func presetParameter(_ i: Int) -> Double {
        let v = presetParameterValues
        return i < v.count ? v[i] : 0
    }

    /// Sets a parameter, kept in its range; a change that gives no geometry is not taken.
    func setPresetParameter(_ i: Int, _ value: Double) {
        let params = PresetGeometry.parameters(presetKey)
        var v = presetParameterValues
        guard i < v.count, i < params.count, value.isFinite else { return }
        v[i] = min(max(value, params[i].range.lowerBound), params[i].range.upperBound)
        guard PresetGeometry.coordinates(presetKey, v) != nil else { return }
        presetValues[presetKey] = v
    }

    /// The parameters differ from the backend's geometry.
    var presetGeometryChanged: Bool {
        let d = PresetGeometry.defaults(presetKey), v = presetParameterValues
        return v.count == d.count && zip(v, d).contains { abs($0 - $1) > 1e-9 }
    }

    func resetPresetGeometry() {
        presetValues[presetKey] = nil
    }

    /// The atoms (Å) sent to the backend when the parameters were changed (nil: the backend's
    /// geometry of the molecule).
    var presetAtoms: [CustomAtom]? {
        guard molecule == .preset, presetGeometryChanged, let p = selectedPreset,
              let x = PresetGeometry.coordinates(presetKey, presetParameterValues), x.count == p.symbols.count else { return nil }
        return zip(p.symbols, x).map { s, c in CustomAtom(symbol: s, Z: Elements.atomicNumber(s), x: c[0], y: c[1], z: c[2]) }
    }

    /// The parameters of the converged optimisation of the displayed result, if it belongs to the
    /// selected molecule of the menu.
    var optimizedPresetParameters: (values: [Double], transitionState: Bool)? {
        guard molecule == .preset, !PresetGeometry.parameters(presetKey).isEmpty, let out = result, out.isGeneral,
              out.input.molecule == presetKey, let opt = out.optimization, opt.converged else { return nil }
        let x = opt.finalGeometry.atoms.map { a in a.position.map { $0 * Units.bohrToAngstrom } }
        guard let v = PresetGeometry.measure(presetKey, x) else { return nil }
        return (v, opt.isTransitionState)
    }

    /// Copying the optimised parameters changes the inputs.
    var canUseOptimizedPresetGeometry: Bool {
        guard let o = optimizedPresetParameters else { return false }
        return zip(o.values, presetParameterValues).contains { abs($0 - $1) > 1e-9 }
    }

    /// Copies the optimised parameters into the inputs (the start of the next job).
    func useOptimizedPresetGeometry() {
        guard let o = optimizedPresetParameters else { return }
        let params = PresetGeometry.parameters(presetKey)
        let v = zip(o.values, params).map { x, p in min(max(x, p.range.lowerBound), p.range.upperBound) }
        if PresetGeometry.coordinates(presetKey, v) != nil { presetValues[presetKey] = v }
    }
}
