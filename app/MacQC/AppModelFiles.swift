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
        for line in text[begin.upperBound...].split(separator: "\n", omittingEmptySubsequences: false) {
            if first {   // the rest of the marker line (" 1")
                first = false
                continue
            }
            var l = Substring(line)
            if l.hasSuffix("\r") { l = l.dropLast() }
            if l.hasPrefix(endMarker) { return data }
            data += l
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
        guard let data = try? encoder.encode(copy), let json = String(data: data, encoding: .utf8) else { return nil }
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

    /// Saves a molecule as a Gaussian input file (.gjf): the route for `method`/`basis` (ΨLab's
    /// names), the name and level as title, charge and multiplicity, Cartesian coordinates in Å.
    func exportGaussianInput(_ m: SavedMolecule, method: String, basis: String) {
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
            try m.gaussianInputText(method: method, basis: basis).write(to: url, atomically: true, encoding: .utf8)
        } catch {
            errorMessage = "Could not write the file: \(error.localizedDescription)"
        }
    }

    /// Method and basis set of the optimised geometry of the displayed result (for its Gaussian
    /// input), else those of the sidebar.
    var gaussianLevel: (method: String, basis: String) {
        if let out = result, out.isGeneral, out.optimization != nil {
            return (out.input.method, out.input.basis)
        }
        return (method.rawValue, currentBasis)
    }
}

extension SavedMolecule {
    /// Gaussian's keyword for a ΨLab method (backend name or title); nil where Gaussian has none
    /// (CCSDT, CCSDTQ and Full CI are not methods of Gaussian).
    static func gaussianMethod(_ method: String) -> String? {
        switch method.uppercased() {
        case "RHF", "HF": return "HF"
        case "CCSD": return "CCSD"
        case "CCSD(T)": return "CCSD(T)"
        default: return nil
        }
    }

    /// The molecule as a Gaussian input file: route "#p method/basis" (HF where Gaussian has no
    /// such method; the title says so), title (name; level and energy of an optimised geometry),
    /// charge and multiplicity, Cartesian coordinates in Å, and the closing blank line. ΨLab reads
    /// it back as a custom molecule (File ▸ Open Structure File…).
    func gaussianInputText(method: String, basis: String) -> String {
        let keyword = SavedMolecule.gaussianMethod(method)
        var title = title.trimmingCharacters(in: .whitespaces)
        if title.isEmpty { title = formula }
        if let l = levelText {
            title += " (" + (transitionState ? "transition state at " : "optimised at ") + l
            if let e = energy { title += String(format: ", E = %.10f E_h", e) }
            title += ", ΨLab)"
        }
        if keyword == nil {
            title += " Gaussian has no \(MethodChoice.title(of: method)): the route asks for HF."
        }
        var lines = ["#p \(keyword ?? "HF")/\(basis)", "", title, "", "\(charge) \(multiplicity)"]
        for a in atoms {
            let symbol = a.symbol.padding(toLength: max(2, a.symbol.count), withPad: " ", startingAt: 0)
            lines.append(symbol + String(format: " %16.10f %16.10f %16.10f", a.x, a.y, a.z))
        }
        lines.append("")
        return lines.joined(separator: "\n") + "\n"
    }
}
