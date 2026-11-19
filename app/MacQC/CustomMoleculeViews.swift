//
//  CustomMoleculeViews.swift
//  MacQC
//
//  ΨLab 9.0: the sidebar views of custom molecules (the name drawn from its LaTeX source, the
//  editable name and charge, the coordinates in Å) and the offer to keep an optimised geometry
//  (molecule menu or .xyz file), shown after every converged optimisation or transition-state
//  search, with its level (method/basis set).
//

import SwiftUI

/// The selected custom molecule under the molecule menu: its name, formula, size and source, the
/// editable name (LaTeX) and charge, and the file actions.
struct CustomMoleculeHeader: View {
    @Bindable var model: AppModel
    let molecule: SavedMolecule
    @Binding var confirmRemove: Bool

    var body: some View {
        let m = molecule
        VStack(alignment: .leading, spacing: 4) {
            LaTeXText(m.nameSource, size: 17, weight: .semibold)
                .textSelection(.enabled)
            if let level = m.levelText {
                Text(levelLine(m, level: level))
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
            Text(sizeLine(m))
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        LabeledContent("Name (LaTeX)") {
            TextField("Name", text: $model.customTitle, prompt: Text(m.formula))
                .labelsHidden()
                .multilineTextAlignment(.trailing)
                .accessibilityLabel("Molecule name")
        }
        // 15.0: charge and multiplicity are set below the molecule menu (every molecule)
        ForEach(m.warnings, id: \.self) { w in
            Text(w)
                .font(.caption)
                .foregroundStyle(.orange)
        }
        HStack {
            Button("Open File…") { model.openStructureFile() }
                .help("Read another structure file (.xyz or Gaussian input); it is added to the molecule menu.")
            Button("Save as XYZ…") { model.exportXYZ(m) }
                .help("Save these coordinates as an .xyz file (Å; the comment line is the name and level).")
            Button("Save as GJF…") {
                model.exportGaussianInput(m, method: model.method.rawValue, basis: model.currentBasis,
                                          orbitals: model.orbitalsForBasis(model.currentBasis))
            }
            .help("Save these coordinates as a Gaussian input file (.gjf; route for the method and basis set of the sidebar, Cartesian coordinates in Å).")
            Button("Remove") { confirmRemove = true }
                .help("Remove this molecule from the molecule menu.")
        }
        .controlSize(.small)
        .disabled(model.isRunning || model.isScanning)
    }
}

extension CustomMoleculeHeader {
    fileprivate func levelLine(_ m: SavedMolecule, level: String) -> String {
        var text = "Geometry "
        text += m.transitionState ? "of the transition state at " : "optimised at "
        text += level
        if let e = m.energy {
            text += String(format: ", E = %.8f Eₕ", e)
        }
        return text
    }

    fileprivate func sizeLine(_ m: SavedMolecule) -> String {
        let atoms = m.atoms.count == 1 ? "1 atom" : "\(m.atoms.count) atoms"
        return "\(Fmt.formula(m.formula)) · \(atoms) · \(m.electrons) electrons · from \(m.source)"
    }
}

/// The Cartesian coordinates of a custom molecule in Å, as read from the file (or optimised).
struct CustomCoordinateTable: View {
    let molecule: SavedMolecule

    var body: some View {
        if molecule.atoms.count > 14 {
            ScrollView(.vertical) {
                grid
            }
            .frame(height: 280)
        } else {
            grid
        }
    }

    private var grid: some View {
        Grid(alignment: .trailing, horizontalSpacing: 10, verticalSpacing: 2) {
            GridRow {
                Text("Atom").gridColumnAlignment(.leading)
                Text("x (Å)")
                Text("y (Å)")
                Text("z (Å)")
            }
            .font(.caption.weight(.semibold))
            .foregroundStyle(.secondary)
            ForEach(Array(molecule.atoms.enumerated()), id: \.offset) { k, a in
                GridRow {
                    Text("\(a.symbol)\(k + 1)")
                    Text(String(format: "%.6f", a.x))
                    Text(String(format: "%.6f", a.y))
                    Text(String(format: "%.6f", a.z))
                }
            }
        }
        .font(.caption.monospacedDigit())
        .textSelection(.enabled)
    }
}

/// After a converged optimisation or transition-state search of the selected molecule: keep the
/// optimised coordinates in the molecule menu or save them as an .xyz file; the level of the
/// optimisation (method/basis set) and the energy go with them. For a custom molecule also "Use
/// optimised geometry" (keep and select, the start of the next job).
struct SaveOptimizedGeometryBox: View {
    @Bindable var model: AppModel

    var body: some View {
        if model.resultMatchesSelection, let m = model.optimizedMolecule {
            VStack(alignment: .leading, spacing: 6) {
                Text(offerText(m))
                    .font(.caption)
                    .foregroundStyle(.secondary)
                HStack {
                    Button(model.optimizedMoleculeSaved ? "In the Molecule Menu" : "Save in Molecule Menu") {
                        model.saveOptimizedMolecule(select: false)
                    }
                    .disabled(model.optimizedMoleculeSaved || model.isRunning)
                    .accessibilityLabel("Save optimised geometry in the molecule menu")
                    Button("Save as XYZ…") { model.exportXYZ(m) }
                        .accessibilityLabel("Save optimised geometry as XYZ")
                    Button("Save as GJF…") {
                        let level = model.gaussianLevel
                        model.exportGaussianInput(m, method: level.method, basis: level.basis, orbitals: level.orbitals)
                    }
                    .accessibilityLabel("Save optimised geometry as Gaussian input")
                    .help("A Gaussian input file (.gjf): the route for the level of the optimisation, Cartesian coordinates in Å.")
                }
                .controlSize(.small)
                if model.molecule == .custom {
                    Button(m.transitionState ? "Use Transition-State Geometry" : "Use Optimised Geometry") {
                        model.saveOptimizedMolecule(select: true)
                    }
                    .secondaryActionStyle()
                    .frame(maxWidth: .infinity)
                    .disabled(model.isRunning)
                    .help("Keeps the optimised coordinates in the molecule menu and selects them for the next job.")
                }
                if let message = model.savedMessage {
                    Text(message)
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
            }
        }
    }

    private func offerText(_ m: SavedMolecule) -> String {
        var text = m.transitionState ? "Transition state" : "Optimised geometry"
        text += " at " + (m.level ?? "")
        if let e = m.energy {
            text += String(format: " (E = %.8f Eₕ)", e)
        }
        text += ": save the coordinates to compute with them later."
        return text
    }
}
