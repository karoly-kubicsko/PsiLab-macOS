//
//  Formatting.swift
//  MacQC
//
//  Number formatting helpers (fixed decimals, "E_h" units).
//

import Foundation

nonisolated enum Fmt {
    /// Energy with 8 decimals, e.g. "-2.86066199". Unit not included.
    static func energy(_ x: Double) -> String {
        String(format: "%.8f", x)
    }

    /// Energy with 8 decimals and unit, e.g. "-2.86066199 E_h".
    static func energyWithUnit(_ x: Double) -> String {
        String(format: "%.8f E_h", x)
    }

    /// Matrix / integral entries with 6 decimals (as in the text report).
    static func value6(_ x: Double) -> String {
        String(format: "%.6f", x)
    }

    static func value8(_ x: Double) -> String {
        String(format: "%.8f", x)
    }

    static func scientific(_ x: Double) -> String {
        String(format: "%.3e", x)
    }

    /// Energy with 10 decimals (general molecules: agrees with the backend's report).
    static func energy10(_ x: Double) -> String {
        String(format: "%.10f", x)
    }

    /// Integer count with thousands separators ("11,029,041"), or "5.62e+11" when very large.
    static func count(_ x: Double) -> String {
        if x >= 1.0e10 { return String(format: "%.2e", x) }
        let formatter = NumberFormatter()
        formatter.numberStyle = .decimal
        formatter.groupingSeparator = ","
        formatter.usesGroupingSeparator = true
        formatter.maximumFractionDigits = 0
        return formatter.string(from: NSNumber(value: x.rounded())) ?? String(format: "%.0f", x)
    }

    /// Seconds as "less than a second", "about 25 s", "about 4 min", "about 1.5 h".
    static func duration(_ seconds: Double) -> String {
        if seconds < 1 { return "less than a second" }
        if seconds < 90 { return String(format: "about %.0f s", seconds) }
        if seconds < 5400 { return String(format: "about %.0f min", seconds / 60) }
        return String(format: "about %.1f h", seconds / 3600)
    }

    /// Memory in MB as "43 MB" or "2.8 GB".
    static func memory(_ megabytes: Double) -> String {
        if megabytes < 1024 { return String(format: "%.0f MB", max(megabytes, 1)) }
        return String(format: "%.1f GB", megabytes / 1024)
    }

    /// One entry per element (atoms of the same element share the basis set).
    static func uniqueAtoms(_ atoms: [MacQCAtomBasis]) -> [MacQCAtomBasis] {
        var seen = Set<String>()
        return atoms.filter { seen.insert($0.symbol).inserted }
    }

    /// Chemical formula with subscripts: "H2O" -> "H₂O".
    static func formula(_ name: String) -> String {
        let subscripts: [Character: Character] = ["0": "₀", "1": "₁", "2": "₂", "3": "₃", "4": "₄",
                                                  "5": "₅", "6": "₆", "7": "₇", "8": "₈", "9": "₉"]
        return String(name.map { subscripts[$0] ?? $0 })
    }

    static func bohr(_ x: Double) -> String {
        String(format: "%.4f bohr", x)
    }

    /// Harmonic wavenumber with "i" for imaginary modes, e.g. "1648.52" or "523.11i".
    static func wavenumber(_ m: MacQCNormalMode) -> String {
        m.imaginary ? String(format: "%.2fi", abs(m.wavenumber)) : String(format: "%.2f", m.wavenumber)
    }

    /// A geometric parameter: "R(O1-H2) = 0.958852 Å" or "A(H2-O1-H3) = 104.3424°".
    static func parameter(_ p: MacQCParameter) -> String {
        p.angle ? String(format: "%@ = %.4f°", p.name, p.value)
                : String(format: "%@ = %.6f Å (%.6f bohr)", p.name, p.value * Units.bohrToAngstrom, p.value)
    }

    /// Value of a geometric parameter (bohr or degrees) as "0.958852 Å" or "104.3424°".
    static func parameterValue(_ v: Double, angle: Bool) -> String {
        angle ? String(format: "%.4f°", v) : String(format: "%.6f Å", v * Units.bohrToAngstrom)
    }

    /// Signed change of a geometric parameter (bohr or degrees).
    static func parameterChange(_ d: Double, angle: Bool) -> String {
        angle ? String(format: "%+.4f°", d) : String(format: "%+.6f Å", d * Units.bohrToAngstrom)
    }

    /// The geometry of a calculation in one line: "R(O–H) = 0.9670 Å (1.8274 bohr), R(O–O) = …,
    /// angle O–O–H = 102.32°, dihedral H–O–O–H = 113.70°" (`plain`: "deg" instead of "°").
    static func moleculeGeometry(_ mol: MacQCMoleculeInfo, plain: Bool = false) -> String {
        let deg = plain ? " deg" : "°"
        let choice = MoleculeChoice(rawValue: mol.name)
        var text = "\(choice?.bondLabel ?? "R") = " + lengthBoth(mol.R)
        if let roo = mol.ROO {
            text += ", R(O–O) = " + lengthBoth(roo)
        }
        if let a = mol.angle, let label = choice?.angleLabel {
            text += String(format: ", %@ = %.2f%@", label.replacingOccurrences(of: "Angle", with: "angle"), a, deg)
        }
        if let d = mol.dihedral {
            text += String(format: ", dihedral H–O–O–H = %.2f%@", d, deg)
        }
        return text
    }

    /// The parameters of a geometry as label/value rows (the molecule panels of the Optimisation
    /// and Scan tabs): H₂O₂ R(O–O), both R(O–H), both O–O–H angles and the dihedral angle; NH₃
    /// the three R(N–H), the three H–N–H angles and 360° − ΣH–N–H; H₂O both R(O–H) and the angle;
    /// diatomics R.
    static func geometryRows(_ g: MacQCGeometry) -> [(String, String)] {
        let a = Units.bohrToAngstrom
        if let r1 = g.R1, let r2 = g.R2, let r3 = g.R3 {
            var rows = [("R(N–H)", String(format: "%.4f, %.4f, %.4f Å", r1 * a, r2 * a, r3 * a))]
            if let a1 = g.angle1, let a2 = g.angle2, let a3 = g.angle3 {
                rows.append(("H–N–H", String(format: "%.2f°, %.2f°, %.2f°", a1, a2, a3)))
            }
            if let p = g.pyramidal {
                rows.append(("360° − ΣH–N–H", String(format: "%.2f°", p)))
            }
            return rows
        }
        if let roo = g.ROO, let r1 = g.R1, let r2 = g.R2 {
            var rows = [("R(O–O)", String(format: "%.6f Å", roo * a)),
                        ("R(O–H)", String(format: "%.5f, %.5f Å", r1 * a, r2 * a))]
            if let a1 = g.angle1, let a2 = g.angle2 {
                rows.append(("O–O–H", String(format: "%.3f°, %.3f°", a1, a2)))
            }
            if let d = g.dihedral {
                rows.append(("H–O–O–H", String(format: "%.3f°", d)))
            }
            return rows
        }
        if let r1 = g.R1, let r2 = g.R2, let angle = g.angle {
            return [("R₁(O–H)", String(format: "%.6f Å", r1 * a)),
                    ("R₂(O–H)", String(format: "%.6f Å", r2 * a)),
                    ("Angle", String(format: "%.4f°", angle))]
        }
        if let r = g.R { return [("R", String(format: "%.6f Å", r * a))] }
        return []
    }

    /// The rows of `geometryRows` in one line.
    static func geometryLine(_ g: MacQCGeometry) -> String {
        geometryRows(g).map { "\($0.0) = \($0.1)" }.joined(separator: "; ")
    }

    /// Point group with subscripts: "C2v" -> "C₂ᵥ", "D3h" -> "D₃ₕ", "Cs" -> "Cₛ".
    static func pointGroup(_ name: String) -> String {
        let subscripts: [Character: Character] = ["1": "₁", "2": "₂", "3": "₃", "4": "₄", "5": "₅", "6": "₆",
                                                  "h": "ₕ", "v": "ᵥ", "s": "ₛ", "i": "ᵢ"]
        var out = ""
        for (i, c) in name.enumerated() {
            out.append(i == 0 ? c : (subscripts[c] ?? c))
        }
        return out
    }

    /// The Hessian of a transition-state search: "CalcFC (1 computed), 1 negative eigenvalue".
    static func tsHessian(_ opt: MacQCOptimization) -> String {
        var text = opt.tsHessian == "calcall" ? "CalcAll" : "CalcFC"
        text += " (\(opt.hessianEvaluations ?? 0) computed)"
        if let n = opt.negativeEigenvalues {
            text += ", \(n) negative eigenvalue" + (n == 1 ? "" : "s")
        }
        return text
    }

    /// A value of a scanned coordinate: "1.4556 Å" or "113.70°".
    static func scanValue(_ v: Double, angular: Bool) -> String {
        angular ? String(format: "%.2f°", v) : String(format: "%.4f Å", v)
    }

    /// Length given in bohr, shown as "0.7743 Å (1.4632 bohr)".
    static func lengthBoth(_ bohr: Double, decimals: Int = 4) -> String {
        let digits = Int32(decimals)  // "%.*f" takes a C int
        return String(format: "%.*f Å (%.*f bohr)", digits, bohr * Units.bohrToAngstrom, digits, bohr)
    }
}

nonisolated enum Units {
    /// 1 bohr in angstrom (CODATA 2022).
    static let bohrToAngstrom = 0.529177210544
    /// J/mol per hartree (CODATA 2018, as the backend's thermochemistry).
    static let joulePerMolPerHartree = 2_625_499.639_479_9
    /// J per thermochemical calorie.
    static let joulePerCalorie = 4.184
    /// kcal/mol per hartree.
    static let kcalPerMolPerHartree = joulePerMolPerHartree / 4184
}
