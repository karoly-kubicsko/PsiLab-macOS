//
//  AppModelFiles.swift
//  MacQC
//
//  ΨLab 13.0: output files that can be opened again, and geometries as Gaussian input.
//
//  Every output file ends (before its termination banner) with the basis set in full, the
//  molecular orbitals of the Orbitals tab and the archive: the complete results as JSON between
//  the lines "@@PSILAB-ARCHIVE-BEGIN 1" and "@@PSILAB-ARCHIVE-END" (the backend writes it for the
//  app and for psilab; the text export of HeH⁺ in STO-nG adds it here). File ▸ Open Output File…
//  reads it back, or a JSON export, and shows every tab as after the calculation — the orbitals
//  included, which are computed again from the basis set and the coefficients of the archive.
//

import SwiftUI
import AppKit
import UniformTypeIdentifiers

/// The archive section of an output file.
nonisolated enum OutputArchive {
    static let beginMarker = "@@PSILAB-ARCHIVE-BEGIN"
    static let endMarker = "@@PSILAB-ARCHIVE-END"

    /// The JSON of the (last) archive in the text of an output file, nil if there is none or it
    /// is incomplete: the lines between the markers joined without line breaks.
    static func json(in text: String) -> String? {
        guard let begin = text.range(of: beginMarker, options: .backwards) else { return nil }
        var data = ""
        var first = true
        // Character.isNewline: "\n", "\r" and "\r\n" (one Character in Swift) end a line
        for line in text[begin.upperBound...].split(omittingEmptySubsequences: false, whereSeparator: \.isNewline) {
            if first {   // the rest of the marker line (" 1")
                first = false
                continue
            }
            if line.hasPrefix(endMarker) { return data }
            data += line
        }
        return nil
    }

    /// An archive section for results whose output file has none (HeH⁺ in STO-nG): the results
    /// as JSON (without the output text), 100 characters per line.
    static func text(for out: MacQCOutput) -> String? {
        var copy = out
        copy.output = nil
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        guard let data = try? encoder.encode(copy), var json = String(data: data, encoding: .utf8),
              json.hasPrefix("{") else { return nil }
        // the "archive" object of the backend's archives: format, program and time
        let version = Bundle.main.object(forInfoDictionaryKey: "CFBundleShortVersionString") as? String ?? ""
        let formatter = DateFormatter()
        formatter.locale = Locale(identifier: "en_US_POSIX")
        formatter.dateFormat = "yyyy-MM-dd HH:mm:ss"
        let header = ["format": "1", "program": "\"ΨLab \(version) (app)\"", "written": "\"\(formatter.string(from: Date()))\""]
        let archive = "\"archive\":{" + ["format", "program", "written"].map { "\"\($0)\":\(header[$0]!)" }.joined(separator: ",") + "}"
        let body = String(json.dropFirst())   // without the opening brace
        json = "{" + archive + (body == "}" ? "" : ",") + body
        var lines = [
            "",
            "---- ΨLab archive: the complete results, read by File > Open Output File… of the app ------",
            "    JSON (format 1, written by the app); the lines between the markers are joined without line breaks",
            beginMarker + " 1",
        ]
        var rest = Substring(json)
        while !rest.isEmpty {
            lines.append(String(rest.prefix(100)))
            rest = rest.dropFirst(100)
        }
        lines.append(endMarker)
        return lines.joined(separator: "\n") + "\n"
    }
}

/// An output file (or JSON export) shown instead of a calculation of this session.
struct OpenedOutputFile: Equatable {
    var name: String
    /// "ΨLab 13.0" and the time the file was written, from its archive (empty if unknown).
    var program: String
    var written: String
}

/// The "archive" object of an archive: the program and the time it was written.
private nonisolated struct ArchiveHeader: Decodable {
    struct Info: Decodable {
        var format: Int?
        var program: String?
        var written: String?
    }
    var archive: Info?
}

extension AppModel {
    /// Output files of ΨLab (.out; .log and .txt as well) and JSON exports.
    static var outputFileTypes: [UTType] {
        var types = ["out", "log"].compactMap { UTType(filenameExtension: $0) }
        types.append(contentsOf: [.plainText, .json])
        return types
    }

    /// File ▸ Open Output File…
    func openOutputFile() {
        guard !isRunning && !isScanning else { return }
        let panel = NSOpenPanel()
        panel.title = "Open an Output File"
        panel.message = "An output file of ΨLab 13.0 or later (the app's text export or psilab's .out file) or a JSON export: every tab shows its results, the orbitals included."
        panel.prompt = "Open"
        panel.allowedContentTypes = AppModel.outputFileTypes
        panel.allowsMultipleSelection = false
        panel.canChooseDirectories = false
        guard panel.runModal() == .OK, let url = panel.url else { return }
        loadOutputFile(from: url)
    }

    /// Reads an output file or a JSON export and shows its results (the inputs in the sidebar stay
    /// as they are). `showSummary`: select the Summary tab.
    func loadOutputFile(from url: URL, showSummary: Bool = true) {
        let name = url.lastPathComponent
        let text: String
        do {
            text = try String(contentsOf: url, encoding: .utf8)
        } catch {
            errorMessage = "“\(name)” could not be read: \(error.localizedDescription)"
            return
        }
        let trimmed = text.trimmingCharacters(in: .whitespacesAndNewlines)
        let isJSON = url.pathExtension.lowercased() == "json" || trimmed.hasPrefix("{")
        guard let json = isJSON ? trimmed : OutputArchive.json(in: text) else {
            errorMessage = "“\(name)” has no ΨLab archive section, so its results cannot be shown. Output files of ΨLab 13.0 and later end with one (the text export of the app and psilab’s .out files); JSON exports of any version can be opened as well."
            return
        }
        var out: MacQCOutput
        do {
            out = try JSONDecoder().decode(MacQCOutput.self, from: Data(json.utf8))
        } catch {
            errorMessage = "The results in “\(name)” could not be read: \(error.localizedDescription)"
            return
        }
        let header = try? JSONDecoder().decode(ArchiveHeader.self, from: Data(json.utf8))
        if out.isGeneral && (out.output ?? "").isEmpty && !isJSON {
            out.output = text   // the Output tab shows the file as it is
        }
        result = out
        errorMessage = out.error
        savedMessage = nil
        outputText = out.isGeneral ? (out.output ?? "") : ""
        outputGeneration += 1
        openedFile = OpenedOutputFile(name: name, program: header?.archive?.program ?? "",
                                      written: header?.archive?.written ?? "")
        if out.scan != nil { scanPage = .relaxed }
        if showSummary { selectedTab = .summary }
    }

    /// The text of File ▸ Export Results as Text…: the output file, with an archive section added
    /// where it has none (HeH⁺ in STO-nG), so that the file can be opened again.
    func exportableText() -> String? {
        guard let out = result, let text = fullTextReport() else { return nil }
        if OutputArchive.json(in: text) != nil { return text }
        guard let archive = OutputArchive.text(for: out) else { return text }
        return text + archive
    }

    // MARK: Gaussian input

    /// Saves a molecule as a Gaussian input file (.gjf): see SavedMolecule.gaussianInputText.
    /// `method`/`basis`: ΨLab's names of the level for the route; `orbitals`: the basis functions
    /// ΨLab used (written out where Gaussian has no such basis set), if known.
    func exportGaussianInput(_ m: SavedMolecule, method: String, basis: String, orbitals: MacQCOrbitals?) {
        let panel = NSSavePanel()
        panel.title = "Save Coordinates as Gaussian Input"
        panel.allowedContentTypes = [UTType(filenameExtension: "gjf") ?? .plainText]
        panel.canCreateDirectories = true
        var stem = m.formula.replacingOccurrences(of: " ", with: "")
        if let level = m.level {
            stem += "-" + level.replacingOccurrences(of: "/", with: "-").replacingOccurrences(of: " ", with: "")
                .replacingOccurrences(of: ",", with: "")
        }
        panel.nameFieldStringValue = stem + ".gjf"
        guard panel.runModal() == .OK, let url = panel.url else { return }
        do {
            try m.gaussianInputText(method: method, basis: basis, orbitals: orbitals)
                .write(to: url, atomically: true, encoding: .utf8)
        } catch {
            errorMessage = "Could not write the file: \(error.localizedDescription)"
        }
    }

    /// The level of the optimised geometry of the displayed result (for its Gaussian input) with
    /// its basis functions, else the level of the sidebar.
    var gaussianLevel: (method: String, basis: String, orbitals: MacQCOrbitals?) {
        if let out = result, out.isGeneral, out.optimization != nil {
            return (out.input.method, out.input.basis, out.orbitals)
        }
        return (method.rawValue, currentBasis, orbitalsForBasis(currentBasis))
    }

    /// The basis functions of the displayed result if it used `basis` (the functions of an element
    /// do not depend on the geometry).
    func orbitalsForBasis(_ basis: String) -> MacQCOrbitals? {
        guard let out = result, out.isGeneral, out.input.basis.lowercased() == basis.lowercased() else { return nil }
        return out.orbitals
    }
}

extension SavedMolecule {
    /// Gaussian's keyword for a ΨLab method; nil where Gaussian has none (CCSDT, CCSDTQ, Full CI).
    /// ΨLab correlates all electrons: Gaussian's coupled cluster freezes the core unless asked
    /// for Full.
    static func gaussianMethod(_ method: String) -> String? {
        var m = method.uppercased()
        // 15.0: the reference prefixes of Gaussian (UMP2, ROMP2, ...)
        var prefix = ""
        if ["RHF", "UHF", "ROHF"].contains(m) { return m }
        if m.hasPrefix("RO") && m != "RO" { prefix = "RO"; m.removeFirst(2) }
        else if m.hasPrefix("U") { prefix = "U"; m.removeFirst() }
        switch m {
        case "HF": return prefix.isEmpty ? "HF" : prefix + "HF"
        case "MP2", "MP3", "MP4", "MP5": return prefix + (m == "MP4" ? "MP4(SDTQ,Full)" : m + "(Full)")
        case "CISD": return prefix + "CISD(Full)"
        case "QCISD": return prefix + "QCISD(Full)"
        case "QCISD(T)": return prefix + "QCISD(T,Full)"
        case "CCSD": return prefix + "CCSD(Full)"
        case "CCSD(T)": return prefix + "CCSD(T,Full)"
        default: return nil
        }
    }

    /// Text for Gaussian's title section: ASCII without the characters it must not contain
    /// (@ # ! – _ \, and the hyphen as well) and without LaTeX markup.
    static func gaussianTitle(_ text: String) -> String {
        var out = ""
        for c in text.unicodeScalars {
            switch c {
            case "$", "{", "}", "^", "_", "\\", "@", "#", "!": continue
            case "-", "\u{2013}", "\u{2212}": out.append(" ")
            default:
                if c.isASCII && c.value >= 32 { out.unicodeScalars.append(c) }
            }
        }
        return out.split(separator: " ").joined(separator: " ")
    }

    /// The molecule as a Gaussian input file:
    ///  - comment lines ("!") with the name, the level and energy of an optimised geometry, and
    ///    what the route cannot say;
    ///  - route "#p method/basis": HF, CCSD(Full) or CCSD(T,Full) (all electrons correlated, as in
    ///    ΨLab; HF where Gaussian has no such method); by Gaussian's name the sets it has built in
    ///    (14.0: the 6-31G family, whose default 6D 7F ΨLab uses as well, the 6-311G family and
    ///    Def2SVP, Def2TZVP, Def2TZVPP, Def2QZVP, Def2QZVPP, all with Gaussian's default function
    ///    types; aug-cc-pVXZ); the others (aug-cc-pCVXZ, def2-mTZVP, def2-mTZVPP, def2-TZVPPD,
    ///    def2-QZVPPD) as "Gen 5D 7F" with ΨLab's functions written out (for aug-cc-pCVXZ H and He
    ///    with the aug-cc-pVXZ functions, as ΨLab uses them), or by name per element when the
    ///    functions are not known;
    ///  - title (the name, ASCII without the characters the title section must not contain);
    ///  - charge and multiplicity, Cartesian coordinates in Å, and the basis section if any.
    /// ΨLab reads it back as a custom molecule (comments are skipped).
    func gaussianInputText(method: String, basis: String, orbitals: MacQCOrbitals?) -> String {
        let keyword = SavedMolecule.gaussianMethod(method)
        let name = SavedMolecule.gaussianTitle(title).isEmpty ? SavedMolecule.gaussianTitle(formula)
                                                              : SavedMolecule.gaussianTitle(title)
        var lines: [String] = []
        if let l = levelText {
            var line = "! \(name): " + (transitionState ? "transition state at " : "geometry optimised at ") + l + " with PsiLab"
            if let e = energy { line += String(format: ", E = %.10f E_h", e) }
            lines.append(line)
        } else {
            lines.append("! \(name): coordinates from " + source)
        }
        if keyword == nil {
            lines.append("! Gaussian has no \(MethodChoice.title(of: method)): the route asks for HF")
        } else if keyword != "HF" {
            lines.append("! all electrons correlated, as in PsiLab (Full: no frozen core)")
        }
        // the basis set
        var route = SavedMolecule.gaussianBasisKeyword(basis) ?? basis
        var basisSection: [String] = []
        if BasisSets.isCartesianD(basis) {
            lines.append("! \(basis): Cartesian d functions (6D, Gaussian's default for the 6-31G family), as in PsiLab")
        }
        let heavy = atoms.contains { $0.Z > 36 }
        if heavy, let section = SavedMolecule.gaussianBasisSection(for: atoms.map(\.Z), orbitals: orbitals) {
            // 14.0: elements beyond Kr (the def2 sets with their ECPs, I in 6-311G): the functions and
            // potentials ΨLab used, exactly (Gaussian's built-in def2 sets lack Ce-Lu)
            let ecp = SavedMolecule.gaussianECPSection(for: atoms.map(\.Z), orbitals: orbitals)
            route = ecp.isEmpty ? "Gen 5D 7F" : "GenECP 5D 7F"
            lines.append("! basis set \(basis) as used by PsiLab, written out below"
                         + (ecp.isEmpty ? "" : ", with its effective core potentials"))
            basisSection = ecp.isEmpty ? section : section + [""] + ecp
        } else if heavy && basis.lowercased().hasPrefix("def2") {
            lines.append("! \(basis): PsiLab uses the def2 effective core potentials for the elements beyond Kr;")
            lines.append("! give them with GenECP (basissetexchange.org) if your Gaussian does not apply them")
        } else if SavedMolecule.gaussianBasisKeyword(basis) == nil && !basis.lowercased().hasPrefix("aug-cc-pcv") {
            // 14.0: def2 sets that Gaussian does not have built in (spherical, 5D 7F)
            route = "Gen 5D 7F"
            if let section = SavedMolecule.gaussianBasisSection(for: atoms.map(\.Z), orbitals: orbitals) {
                lines.append("! basis set \(basis) as used by PsiLab, written out below")
                basisSection = section
            } else {
                lines.append("! basis set \(basis) by name: if your Gaussian does not have it,")
                lines.append("! replace the names by its functions (e.g. from basissetexchange.org)")
                var seen: [Int] = []
                for Z in atoms.map(\.Z) where !seen.contains(Z) {
                    seen.append(Z)
                    basisSection += ["\(Elements.symbol(Z)) 0", basis, "****"]
                }
            }
        } else if basis.lowercased().hasPrefix("aug-cc-pcv") {
            route = "Gen 5D 7F"
            let partner = "aug-cc-pV" + basis.dropFirst("aug-cc-pCV".count)
            if let section = SavedMolecule.gaussianBasisSection(for: atoms.map(\.Z), orbitals: orbitals) {
                lines.append("! basis set \(basis) as used by PsiLab (\(partner) for H and He), written out below")
                basisSection = section
            } else {
                lines.append("! basis set \(basis) (\(partner) for H and He), by name: if your Gaussian does not have it,")
                lines.append("! replace the names by its functions (e.g. from basissetexchange.org)")
                var seen: [Int] = []
                for Z in atoms.map(\.Z) where !seen.contains(Z) {
                    seen.append(Z)
                    basisSection += ["\(Elements.symbol(Z)) 0", Z <= 2 ? partner : basis, "****"]
                }
            }
        }
        lines += ["#p \(keyword ?? "HF")/\(route)", "", name.isEmpty ? "molecule" : name, "", "\(charge) \(multiplicity)"]
        for a in atoms {
            let symbol = a.symbol.padding(toLength: max(2, a.symbol.count), withPad: " ", startingAt: 0)
            lines.append(symbol + String(format: " %16.10f %16.10f %16.10f", a.x, a.y, a.z))
        }
        lines.append("")
        if !basisSection.isEmpty {
            lines += basisSection
            lines.append("")
        }
        return lines.joined(separator: "\n") + "\n"
    }

    /// 14.0: Gaussian's keyword of a ΨLab basis set that Gaussian has built in (with the same
    /// default function types as ΨLab: 6D 7F for the 6-31G family, 5D 7F for the others), nil for
    /// aug-cc-pCVXZ, def2-mTZVP, def2-mTZVPP, def2-TZVPPD and def2-QZVPPD.
    static func gaussianBasisKeyword(_ basis: String) -> String? {
        let b = basis.lowercased()
        if b.hasPrefix("6-31") || b.hasPrefix("aug-cc-pv") { return basis }
        switch b {
        case "sto-3g": return "STO-3G"
        case "def2-svp": return "Def2SVP"
        case "def2-tzvp": return "Def2TZVP"
        case "def2-tzvpp": return "Def2TZVPP"
        case "def2-qzvp": return "Def2QZVP"
        case "def2-qzvpp": return "Def2QZVPP"
        default: return nil
        }
    }

    /// Gaussian's general basis input for the elements `Z` from the shells of `orbitals` (one atom
    /// of each element): shell type, number of primitives, exponents and the contraction
    /// coefficients of normalised primitives (ΨLab stores them with the primitive normalisation).
    /// nil if an element is missing.
    static func gaussianBasisSection(for Z: [Int], orbitals: MacQCOrbitals?) -> [String]? {
        guard let orbitals else { return nil }
        let letters = ["S", "P", "D", "F", "G", "H", "I"]
        var section: [String] = []
        var seen: [Int] = []
        for z in Z where !seen.contains(z) {
            seen.append(z)
            guard let atom = orbitals.atoms.firstIndex(where: { $0.Z == z }) else { return nil }
            section.append("\(Elements.symbol(z)) 0")
            for shell in orbitals.shells where shell.atom == atom {
                guard shell.l < letters.count, shell.exponents.count == shell.coefficients.count else { return nil }
                section.append(String(format: "%@ %3ld 1.00", letters[shell.l], shell.exponents.count))
                for (a, c) in zip(shell.exponents, shell.coefficients) {
                    section.append(String(format: "  %20.10E  %20.10E", a, c / primitiveNorm(l: shell.l, exponent: a)))
                }
            }
            section.append("****")
        }
        return section
    }

    /// 14.0: Gaussian's ECP input for the elements `Z` that have an effective core potential in
    /// `orbitals` (one atom of each element): "I 0", "I-ECP L n_core", then the local part and the
    /// parts U_l − U_L, each a title, the number of terms and "2 exponent coefficient" lines.
    static func gaussianECPSection(for Z: [Int], orbitals: MacQCOrbitals?) -> [String] {
        guard let orbitals else { return [] }
        let letters = ["s", "p", "d", "f", "g", "h", "i"]
        var section: [String] = []
        var seen: [Int] = []
        for z in Z where !seen.contains(z) {
            seen.append(z)
            guard let atom = orbitals.atoms.first(where: { $0.Z == z }), let ecp = atom.ecp,
                  ecp.L < letters.count else { continue }
            let symbol = Elements.symbol(z)
            let top = letters[ecp.L]
            section.append("\(symbol) 0")
            section.append("\(symbol)-ECP \(ecp.L) \(ecp.core)")
            var parts: [(String, [[Double]])] = [("\(top) potential", ecp.local)]
            for (l, terms) in ecp.projected.enumerated() where l < letters.count {
                parts.append(("\(letters[l])-\(top) potential", terms))
            }
            for (title, terms) in parts {
                section.append(title)
                // a part without terms (all coefficients zero) as one zero term
                let list = terms.isEmpty ? [[1.0, 0.0]] : terms
                section.append(String(format: "  %ld", list.count))
                for t in list where t.count == 2 {
                    section.append(String(format: "2  %18.10f  %18.10f", t[0], t[1]))
                }
            }
        }
        return section
    }

    /// Normalisation of the primitive x^l exp(-a r^2): (2a/π)^(3/4) (4a)^(l/2) / √((2l−1)!!).
    static func primitiveNorm(l: Int, exponent a: Double) -> Double {
        var doubleFactorial = 1.0
        var k = 2 * l - 1
        while k > 1 {
            doubleFactorial *= Double(k)
            k -= 2
        }
        return pow(2 * a / Double.pi, 0.75) * pow(4 * a, Double(l) / 2) / doubleFactorial.squareRoot()
    }
}
