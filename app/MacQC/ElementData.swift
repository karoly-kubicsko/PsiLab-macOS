//
//  ElementData.swift
//  MacQC
//
//  ΨLab 9.0: the elements H … Og (Z = 1 … 118) for custom molecules: symbols, covalent radii
//  (bonds of the drawings; the same table as the backend: Pyykkö & Atsumi for Z ≤ 36, Cordero et
//  al. 2008 for Z = 37 … 96, 1.50 Å above) and CPK colours (Jmol) for the drawings.
//

import SwiftUI

nonisolated enum Elements {
    static let symbols: [String] = [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
        "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
        "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
        "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
        "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
        "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
        "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
        "Fl", "Mc", "Lv", "Ts", "Og"]

    /// Covalent radii in Å, index Z − 1.
    static let covalentRadii: [Double] = [
        0.32, 0.46, 1.33, 1.02, 0.85, 0.75, 0.71, 0.63, 0.64, 0.67, 1.55, 1.39, 1.26, 1.16, 1.11, 1.03,
        0.99, 0.96, 1.96, 1.71, 1.48, 1.36, 1.34, 1.22, 1.19, 1.16, 1.11, 1.10, 1.12, 1.18, 1.24, 1.21,
        1.21, 1.16, 1.14, 1.17, 2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44,
        1.42, 1.39, 1.39, 1.38, 1.39, 1.40, 2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96,
        1.94, 1.92, 1.92, 1.89, 1.90, 1.87, 1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32,
        1.45, 1.46, 1.48, 1.40, 1.50, 1.50, 2.60, 2.21, 2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.69,
        1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,
        1.50, 1.50, 1.50, 1.50, 1.50, 1.50]

    static func symbol(_ Z: Int) -> String {
        Z >= 1 && Z <= symbols.count ? symbols[Z - 1] : "?"
    }

    /// Atomic number of a symbol (any case), 0 if unknown.
    static func atomicNumber(_ symbol: String) -> Int {
        let s = symbol.lowercased()
        return (symbols.firstIndex { $0.lowercased() == s } ?? -1) + 1
    }

    /// Covalent radius in bohr.
    static func covalentRadiusBohr(_ Z: Int) -> Double {
        let r = Z >= 1 && Z <= covalentRadii.count ? covalentRadii[Z - 1] : 1.5
        return r / Units.bohrToAngstrom
    }

    /// ΨLab 12.0: van der Waals radius in Å for the space-filling and ball-and-stick views of the
    /// 3D viewer: Bondi (1964), with Mantina et al. (2009) for the main-group elements Bondi left
    /// out; 2.0 Å for the elements of neither.
    static func vanDerWaalsRadius(_ Z: Int) -> Double {
        switch Z {
        case 1: return 1.20
        case 2: return 1.40
        case 3: return 1.82
        case 4: return 1.53
        case 5: return 1.92
        case 6: return 1.70
        case 7: return 1.55
        case 8: return 1.52
        case 9: return 1.47
        case 10: return 1.54
        case 11: return 2.27
        case 12: return 1.73
        case 13: return 1.84
        case 14: return 2.10
        case 15: return 1.80
        case 16: return 1.80
        case 17: return 1.75
        case 18: return 1.88
        case 19: return 2.75
        case 20: return 2.31
        case 28: return 1.63
        case 29: return 1.40
        case 30: return 1.39
        case 31: return 1.87
        case 32: return 2.11
        case 33: return 1.85
        case 34: return 1.90
        case 35: return 1.85
        case 36: return 2.02
        case 37: return 3.03
        case 38: return 2.49
        case 46: return 1.63
        case 47: return 1.72
        case 48: return 1.58
        case 49: return 1.93
        case 50: return 2.17
        case 51: return 2.06
        case 52: return 2.06
        case 53: return 1.98
        case 54: return 2.16
        case 55: return 3.43
        case 56: return 2.68
        case 78: return 1.75
        case 79: return 1.66
        case 80: return 1.55
        case 81: return 1.96
        case 82: return 2.02
        case 83: return 2.07
        case 84: return 1.97
        case 85: return 2.02
        case 86: return 2.20
        case 87: return 3.48
        case 88: return 2.83
        case 92: return 1.86
        default: return 2.00
        }
    }

    /// CPK colour (Jmol), as red, green, blue in 0 … 1; grey-blue for elements not listed.
    static func cpk(_ Z: Int) -> (Double, Double, Double) {
        switch Z {
        case 1: return (0.93, 0.93, 0.93)
        case 2: return (0.75, 0.95, 1.0)
        case 3: return (0.75, 0.45, 1.0)
        case 4: return (0.76, 1.0, 0.0)
        case 5: return (1.0, 0.71, 0.71)
        case 6: return (0.45, 0.45, 0.45)
        case 7: return (0.19, 0.31, 0.97)
        case 8: return (0.92, 0.16, 0.12)
        case 9: return (0.56, 0.88, 0.31)
        case 10: return (0.70, 0.89, 0.96)
        case 11: return (0.67, 0.36, 0.95)
        case 12: return (0.54, 1.0, 0.0)
        case 13: return (0.75, 0.65, 0.65)
        case 14: return (0.94, 0.78, 0.63)
        case 15: return (1.0, 0.50, 0.0)
        case 16: return (0.95, 0.85, 0.19)
        case 17: return (0.12, 0.94, 0.12)
        case 18: return (0.50, 0.82, 0.89)
        case 19: return (0.56, 0.25, 0.83)
        case 20: return (0.24, 1.0, 0.0)
        case 26: return (0.88, 0.40, 0.20)
        case 29: return (0.78, 0.50, 0.20)
        case 30: return (0.49, 0.50, 0.69)
        case 31: return (0.76, 0.56, 0.56)
        case 32: return (0.40, 0.56, 0.56)
        case 33: return (0.74, 0.50, 0.89)
        case 34: return (1.0, 0.63, 0.0)
        case 35: return (0.65, 0.16, 0.16)
        case 36: return (0.36, 0.72, 0.82)
        case 53: return (0.58, 0.0, 0.58)
        case 54: return (0.26, 0.62, 0.69)
        case 79: return (1.0, 0.82, 0.14)
        default: return (0.55, 0.58, 0.68)
        }
    }
}
