//
//  CustomMolecules.swift
//  MacQC
//
//  ΨLab 9.0: custom molecules. A structure file (.xyz in Å, the comment line is the name, which may
//  contain LaTeX; or a Gaussian input .gjf/.com with charge and multiplicity) is read by the
//  backend (macqc_parse_structure) and kept by the app: the molecule menu lists every kept
//  molecule by its name, like H₂ or H₂O. An optimised geometry can be kept as well; it carries its
//  level (method/basis set, e.g. CCSD/aug-cc-pCVTZ) and its energy. The list is stored as JSON in
//  the app's Application Support folder (Molecules.json).
//

import Foundation

/// An atom of a custom molecule; coordinates in Å, as in the structure file.
nonisolated struct CustomAtom: Codable, Sendable, Equatable, Hashable {
    var symbol: String
    var Z: Int
    var x: Double
    var y: Double
    var z: Double

    /// The atom of a backend result (coordinates in bohr).
    init(_ atom: MacQCAtom) {
        symbol = atom.symbol
        Z = atom.Z
        let p = (0..<3).map { $0 < atom.position.count ? atom.position[$0] * Units.bohrToAngstrom : 0 }
        x = p[0]
        y = p[1]
        z = p[2]
    }

    init(symbol: String, Z: Int, x: Double, y: Double, z: Double) {
        self.symbol = symbol
        self.Z = Z
        self.x = x
        self.y = y
        self.z = z
    }
}

/// A custom molecule kept by the app.
nonisolated struct SavedMolecule: Codable, Sendable, Equatable, Identifiable, Hashable {
    var id: UUID = UUID()
    /// The name: the comment line of the .xyz file or the title of the Gaussian input (LaTeX).
    var title: String
    var atoms: [CustomAtom]          // Å
    var charge: Int = 0
    var multiplicity: Int = 1
    /// Where the coordinates come from: the file name, or the calculation that optimised them.
    var source: String
    /// nil: the coordinates of the file; else the level of the optimisation, e.g. "CCSD/aug-cc-pCVTZ".
    var level: String? = nil
    var transitionState = false
    var energy: Double? = nil        // E_h at the optimised geometry
    var added = Date()
    var warnings: [String] = []

    var electrons: Int { atoms.reduce(0) { $0 + $1.Z } - charge }

    /// Molecular formula in Hill order: "CO2", "CH4", "H2O" (with the charge: "NH4+", "SO4 2-").
    var formula: String {
        var counts: [(String, Int)] = []
        for a in atoms {
            if let k = counts.firstIndex(where: { $0.0 == a.symbol }) {
                counts[k].1 += 1
            } else {
                counts.append((a.symbol, 1))
            }
        }
        let carbon = counts.contains { $0.0 == "C" }
        func rank(_ s: String) -> Int { carbon ? (s == "C" ? 0 : (s == "H" ? 1 : 2)) : 2 }
        counts.sort { rank($0.0) != rank($1.0) ? rank($0.0) < rank($1.0) : $0.0 < $1.0 }
        var f = counts.map { $0.0 + ($0.1 > 1 ? "\($0.1)" : "") }.joined()
        if charge == 1 { f += "+" } else if charge == -1 { f += "-" } else if charge != 0 {
            f += " \(abs(charge))" + (charge > 0 ? "+" : "-")
        }
        return f
    }

    /// The name as LaTeX source (the formula with subscripts when the file gave no name).
    var nameSource: String {
        let t = title.trimmingCharacters(in: .whitespaces)
        if !t.isEmpty { return t }
        var s = ""
        var digits = ""
        for c in formula {
            if c.isNumber {
                digits.append(c)
            } else {
                if !digits.isEmpty { s += "$_{\(digits)}$"; digits = "" }
                s.append(c)
            }
        }
        if !digits.isEmpty { s += "$_{\(digits)}$" }
        return s
    }

    /// The level as shown: "CCSD/aug-cc-pCVTZ" or "TS · CCSD/aug-cc-pCVTZ".
    var levelText: String? {
        level.map { (transitionState ? "TS · " : "") + $0 }
    }

    /// Plain text for menus: "CO₂", "CO₂ · CCSD/aug-cc-pCVTZ".
    var menuTitle: String {
        LaTeX.plain(nameSource) + (levelText.map { " · " + $0 } ?? "")
    }

    /// Labels with numbers: "C1", "O2", …
    var atomLabels: [String] {
        atoms.enumerated().map { "\($0.element.symbol)\($0.offset + 1)" }
    }

    /// The molecule as an XYZ file (Å); the comment line is the name, followed by the level and
    /// the energy of an optimised geometry.
    func xyzText() -> String {
        var comment = title
        if let l = levelText {
            comment += (comment.isEmpty ? "" : " ") + "(optimised " + l
            if let e = energy { comment += String(format: ", E = %.10f E_h", e) }
            comment += ")"
        }
        var lines = ["\(atoms.count)", comment]
        for a in atoms {
            let symbol = a.symbol.padding(toLength: max(2, a.symbol.count), withPad: " ", startingAt: 0)
            lines.append(symbol + String(format: " %16.10f %16.10f %16.10f", a.x, a.y, a.z))
        }
        return lines.joined(separator: "\n") + "\n"
    }

    /// Same molecule and geometry (for not keeping a file twice).
    func sameStructure(as other: SavedMolecule) -> Bool {
        guard title == other.title, charge == other.charge, level == other.level,
              atoms.count == other.atoms.count else { return false }
        for (a, b) in zip(atoms, other.atoms) {
            if a.Z != b.Z || abs(a.x - b.x) > 1e-8 || abs(a.y - b.y) > 1e-8 || abs(a.z - b.z) > 1e-8 { return false }
        }
        return true
    }
}

/// The structure file as read by the backend.
nonisolated struct ParsedStructure: Decodable, Sendable {
    var error: String?
    var format: String?
    var title: String?
    var charge: Int?
    var multiplicity: Int?
    var hasCharge: Bool?
    var warnings: [String]
    var atoms: [CustomAtom]
}

/// The kept molecules on disk (Application Support/ΨLab/Molecules.json).
nonisolated enum MoleculeStore {
    static var fileURL: URL? {
        guard let base = FileManager.default.urls(for: .applicationSupportDirectory, in: .userDomainMask).first else {
            return nil
        }
        return base.appendingPathComponent("ΨLab", isDirectory: true).appendingPathComponent("Molecules.json")
    }

    static func load() -> [SavedMolecule] {
        guard let url = fileURL, let data = try? Data(contentsOf: url) else { return [] }
        let decoder = JSONDecoder()
        decoder.dateDecodingStrategy = .iso8601
        return (try? decoder.decode([SavedMolecule].self, from: data)) ?? []
    }

    /// Writes the list; returns an error message or nil.
    @discardableResult
    static func save(_ molecules: [SavedMolecule]) -> String? {
        guard let url = fileURL else { return "The Application Support folder is not available." }
        do {
            try FileManager.default.createDirectory(at: url.deletingLastPathComponent(), withIntermediateDirectories: true)
            let encoder = JSONEncoder()
            encoder.dateEncodingStrategy = .iso8601
            encoder.outputFormatting = [.prettyPrinted, .sortedKeys]
            try encoder.encode(molecules).write(to: url, options: .atomic)
            return nil
        } catch {
            return error.localizedDescription
        }
    }
}
