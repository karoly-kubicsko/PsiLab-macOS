//
//  AppModelCustom.swift
//  MacQC
//
//  ΨLab 9.0: the custom molecules of the app model. Opening a structure file (File ▸ Open
//  Structure File…, or "Custom molecule…" in the molecule menu) reads it with the backend and
//  keeps the molecule: it appears in the molecule menu by its name. After a converged geometry
//  optimisation (or transition-state search) of any molecule, its coordinates can be kept as
//  well, with the level of the optimisation (method/basis set) and the energy, or saved as an
//  .xyz file.
//

import SwiftUI
import AppKit
import UniformTypeIdentifiers

/// An item of the molecule menu: a built-in molecule, a kept custom molecule, or "Custom
/// molecule…" (opens a structure file).
enum MoleculeSelection: Hashable {
    case builtIn(MoleculeChoice)
    case saved(UUID)
    case openFile
}

extension AppModel {

    // MARK: Selection

    /// The selected custom molecule (molecule == .custom).
    var customMolecule: SavedMolecule? {
        guard let id = customID else { return nil }
        return savedMolecules.first { $0.id == id }
    }

    /// The selection of the molecule menu. Choosing "Custom molecule…" opens the file panel; the
    /// menu keeps the previous molecule until a file has been read.
    var moleculeSelection: MoleculeSelection {
        get {
            if molecule == .custom, let id = customID { return .saved(id) }
            return .builtIn(molecule == .custom ? .h2o : molecule)
        }
        set {
            switch newValue {
            case .builtIn(let m):
                molecule = m
            case .saved(let id):
                selectSaved(id)
            case .openFile:
                // after the menu has closed
                Task { @MainActor in self.openStructureFile() }
            }
        }
    }

    func selectSaved(_ id: UUID) {
        guard savedMolecules.contains(where: { $0.id == id }) else { return }
        customID = id
        molecule = .custom
    }

    /// Number of atoms and their labels ("C1", "O2", …) of the selected molecule.
    var atomCount: Int { molecule == .custom ? (customMolecule?.atoms.count ?? 0) : molecule.atomCount }
    var atomLabels: [String] { molecule == .custom ? (customMolecule?.atomLabels ?? []) : molecule.atomLabels }

    /// The name of the selected molecule as LaTeX source (custom) or plain text.
    var moleculeNameSource: String {
        if molecule == .custom { return customMolecule?.nameSource ?? "Custom molecule" }
        return molecule.title
    }

    // MARK: Reading structure files

    /// .xyz, .gjf, .com, .gau and plain text (recognised from its content).
    static var structureTypes: [UTType] {
        var types = ["xyz", "gjf", "com", "gau"].compactMap { UTType(filenameExtension: $0) }
        types.append(.plainText)
        return types
    }

    /// The file panel for a structure file; the molecule is kept and selected.
    func openStructureFile() {
        guard !isRunning && !isScanning else { return }
        let panel = NSOpenPanel()
        panel.title = "Open a Structure File"
        panel.message = "An XYZ file (coordinates in Å; the comment line is the name of the molecule, LaTeX allowed, e.g. CO$_2$) or a Gaussian input file (.gjf, .com)."
        panel.prompt = "Open"
        panel.allowedContentTypes = AppModel.structureTypes
        panel.allowsMultipleSelection = false
        panel.canChooseDirectories = false
        guard panel.runModal() == .OK, let url = panel.url else { return }
        loadStructure(from: url)
    }

    func loadStructure(from url: URL) {
        let text: String
        do {
            text = try String(contentsOf: url, encoding: .utf8)
        } catch {
            if let latin = try? String(contentsOf: url, encoding: .isoLatin1) {
                text = latin
            } else {
                reportStructureError("\(url.lastPathComponent) could not be read: \(error.localizedDescription)")
                return
            }
        }
        addStructure(text: text, fileName: url.lastPathComponent, select: true)
    }

    /// Reads a structure given as text and keeps the molecule (a molecule already kept with the
    /// same name, charge and coordinates is selected instead of being added twice).
    @discardableResult
    func addStructure(text: String, fileName: String, select: Bool) -> SavedMolecule? {
        let parsed: ParsedStructure
        do {
            parsed = try MacQCEngine.parseStructure(text, fileName: fileName)
        } catch {
            let message = (error as? LocalizedError)?.errorDescription ?? error.localizedDescription
            reportStructureError("\(fileName): \(message)")
            return nil
        }
        var m = SavedMolecule(title: (parsed.title ?? "").trimmingCharacters(in: .whitespaces),
                              atoms: parsed.atoms, source: fileName)
        m.charge = parsed.charge ?? 0
        m.multiplicity = parsed.multiplicity ?? 1
        m.warnings = parsed.warnings
        if let existing = savedMolecules.first(where: { $0.sameStructure(as: m) }) {
            if select { selectSaved(existing.id) }
            return existing
        }
        savedMolecules.append(m)
        persistMolecules()
        if select { selectSaved(m.id) }
        return m
    }

    func reportStructureError(_ message: String) {
        structureError = message
        showStructureError = true
    }

    // MARK: Keeping and editing

    func persistMolecules() {
        if let problem = MoleculeStore.save(savedMolecules) {
            errorMessage = "The molecule list could not be saved: \(problem)"
        }
    }

    /// Changes the selected custom molecule and saves the list.
    func updateCustom(_ change: (inout SavedMolecule) -> Void) {
        guard let id = customID, let k = savedMolecules.firstIndex(where: { $0.id == id }) else { return }
        change(&savedMolecules[k])
        persistMolecules()
    }

    /// The name of the selected custom molecule (LaTeX source).
    var customTitle: String {
        get { customMolecule?.title ?? "" }
        set { updateCustom { $0.title = newValue } }
    }

    /// The charge of the selected custom molecule.
    var customCharge: Int {
        get { customMolecule?.charge ?? 0 }
        set { updateCustom { $0.charge = newValue } }
    }

    /// Removes the selected custom molecule from the menu (H₂O is selected).
    func removeCustomMolecule() {
        guard let id = customID else { return }
        savedMolecules.removeAll { $0.id == id }
        persistMolecules()
        customID = nil
        molecule = .h2o
    }

    /// Saves a molecule as an .xyz file (Å; the comment line is its name and level).
    func exportXYZ(_ m: SavedMolecule) {
        let panel = NSSavePanel()
        panel.title = "Save Coordinates as XYZ"
        panel.allowedContentTypes = [UTType(filenameExtension: "xyz") ?? .plainText]
        panel.canCreateDirectories = true
        var stem = m.formula.replacingOccurrences(of: " ", with: "")
        if let level = m.level {
            stem += "-" + level.replacingOccurrences(of: "/", with: "-").replacingOccurrences(of: " ", with: "")
                .replacingOccurrences(of: ",", with: "")
        }
        panel.nameFieldStringValue = stem + ".xyz"
        guard panel.runModal() == .OK, let url = panel.url else { return }
        do {
            try m.xyzText().write(to: url, atomically: true, encoding: .utf8)
        } catch {
            errorMessage = "Could not write the file: \(error.localizedDescription)"
        }
    }

    // MARK: Results of custom molecules

    /// The displayed result belongs to the selected custom molecule (same atoms in the same order
    /// and the same name).
    var resultMatchesCustom: Bool {
        guard molecule == .custom, let c = customMolecule, let out = result, out.input.isCustom,
              let atoms = out.input.atoms else { return false }
        return atoms.map(\.Z) == c.atoms.map(\.Z) && (out.input.title ?? "") == c.title
    }

    /// The displayed result belongs to the selected molecule (built-in or custom).
    var resultMatchesSelection: Bool {
        guard let out = result, out.isGeneral, isGeneral else { return false }
        return molecule == .custom ? resultMatchesCustom : out.input.molecule == molecule.rawValue
    }

    /// The converged optimisation (minimum or transition state) of the displayed result as a
    /// molecule that can be kept: its atoms, the level (method/basis set) and the energy.
    var optimizedMolecule: SavedMolecule? {
        guard let out = result, out.isGeneral, let opt = out.optimization, opt.converged,
              !opt.finalGeometry.atoms.isEmpty else { return nil }
        let title = out.input.isCustom ? (out.input.title ?? "") : out.moleculeNameSource
        var source = MoleculeChoice(rawValue: out.input.molecule)?.title ?? "custom molecule"
        if out.input.isCustom, resultMatchesCustom, let c = customMolecule { source = c.source }
        var m = SavedMolecule(title: title, atoms: opt.finalGeometry.atoms.map { CustomAtom($0) },
                              source: (opt.isTransitionState ? "transition-state search from " : "optimised from ") + source)
        m.charge = out.molecule?.charge ?? out.input.charge ?? 0
        m.level = MethodChoice.title(of: out.input.method) + "/" + out.input.basis
        m.transitionState = opt.isTransitionState
        m.energy = opt.finalEnergy
        return m
    }

    /// The optimised geometry of the displayed result is in the molecule menu already.
    var optimizedMoleculeSaved: Bool {
        guard let m = optimizedMolecule else { return false }
        return savedMolecules.contains { $0.sameStructure(as: m) }
    }

    /// Keeps the optimised geometry in the molecule menu (`select`: and selects it, the start of
    /// the next job).
    func saveOptimizedMolecule(select: Bool) {
        guard let m = optimizedMolecule else { return }
        if let existing = savedMolecules.first(where: { $0.sameStructure(as: m) }) {
            savedMessage = "In the molecule menu as “\(existing.menuTitle)”."
            if select { selectSaved(existing.id) }
            return
        }
        savedMolecules.append(m)
        persistMolecules()
        savedMessage = "Saved in the molecule menu as “\(m.menuTitle)”."
        if select { selectSaved(m.id) }
    }

    /// A point of a relaxed scan of a custom molecule as the start of a transition-state search:
    /// kept in the menu (with its level and point number) and selected, job TS optimisation +
    /// frequencies.
    func startCustomTransitionStateSearch(from point: MacQCScanPoint) {
        guard resultMatchesCustom, let c = customMolecule, let out = result else { return }
        var m = c
        m.id = UUID()
        m.atoms = point.geometry.atoms.map { CustomAtom($0) }
        m.level = "scan point \(point.index + 1), " + MethodChoice.title(of: out.input.method) + "/" + out.input.basis
        m.transitionState = false
        m.energy = point.energy
        m.source = "relaxed scan of " + c.source
        m.added = Date()
        m.warnings = []
        savedMolecules.append(m)
        persistMolecules()
        selectSaved(m.id)
        job = .tsFreq
    }

    /// Why the potential curve is not available for the selected molecule, or nil.
    var curveUnavailableReason: String? {
        molecule == .custom
            ? "The potential curve stretches the X–H bonds of the built-in molecules. For a custom molecule use the job Relaxed scan (any bond, angle or dihedral angle)."
            : nil
    }
}
