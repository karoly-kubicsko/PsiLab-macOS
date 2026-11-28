//
//  MacQCApp.swift
//  MacQC
//
//  App entry point and menu commands.
//
//  23.0: every window, and every tab of a window (View ▸ Show Tab Bar, the + of the tab bar,
//  File ▸ New Tab), has its own model (AppModel with its NuclearModel): one tab can show an
//  electronic-structure calculation while another shows a nuclear-motion calculation, or two
//  nuclear-motion calculations side by side, each with its own inputs, results and running job.
//  The menu commands act on the model of the focused window or tab (@FocusedValue). 22.0 and
//  earlier had one model for the app, so a new tab showed the same calculation as the first.
//

import SwiftUI
import AppKit

@main
struct MacQCApp: App {
    var body: some Scene {
        WindowGroup(id: PsiLabWindow.id) {
            PsiLabWindow()
        }
        .defaultSize(width: 1180, height: 760)
        .commands {
            PsiLabCommands()
        }
    }
}

/// 23.0: the model of the focused window or tab, for the menu commands.
struct AppModelFocusedKey: FocusedValueKey {
    typealias Value = AppModel
}

extension FocusedValues {
    var appModel: AppModel? {
        get { self[AppModelFocusedKey.self] }
        set { self[AppModelFocusedKey.self] = newValue }
    }
}

/// 23.0: the content of a window or tab with its own model. The model is created once, when the
/// window appears (an initial value of @State would be built again whenever SwiftUI recreates this
/// view, and AppModel reads the startup examples).
struct PsiLabWindow: View {
    static let id = "psilab"
    @State private var model: AppModel?

    var body: some View {
        Group {
            if let model {
                ContentView(model: model)
                    .focusedSceneValue(\.appModel, model)
            } else {
                Color.clear
            }
        }
        .frame(minWidth: 1000, minHeight: 640)
        .onAppear {
            if model == nil { model = AppModel() }
        }
    }
}

/// The menu commands, for the model of the focused window or tab (disabled without one).
struct PsiLabCommands: Commands {
    @FocusedValue(\.appModel) private var focused: AppModel?
    @Environment(\.openWindow) private var openWindow

    var body: some Commands {
        // without a focused ΨLab window (all closed, or another window in front) every item but New
        // Window and New Tab is disabled; the placeholder model is never shown
        commands(focused ?? AppModel.placeholder, off: focused == nil)
    }

    @CommandsBuilder
    private func commands(_ model: AppModel, off: Bool) -> some Commands {
        // File > (23.0) New Window, New Tab, Open Structure File…, Open Output File… (13.0), Save
        // Optimised Geometry as XYZ… / as Gaussian Input… (13.0), Export Results as Text… (the
        // output file with its archive) / as JSON… / Potential Curve as Text…
        CommandGroup(replacing: .newItem) {
            // 23.0: a new window or a new tab of the front window, each with its own calculation
            Button("New Window") {
                openWindow(id: PsiLabWindow.id)
            }
            .keyboardShortcut("n", modifiers: [.command])

            Button("New Tab") {
                // the action of the + button of the tab bar; without a ΨLab window a new window
                if !NSApp.sendAction(#selector(NSResponder.newWindowForTab(_:)), to: nil, from: nil) {
                    openWindow(id: PsiLabWindow.id)
                }
            }
            .keyboardShortcut("t", modifiers: [.command])

            Divider()

            // 9.0: a custom molecule from a structure file (kept in the molecule menu)
            Button("Open Structure File…") {
                model.openStructureFile()
            }
            .keyboardShortcut("o", modifiers: [.command])
            .disabled(off || model.isRunning || model.isScanning)

            // 13.0: an output file (archive section) or a JSON export, every tab restored; 22.0: also
            // nuclear-motion output files (the nuclear-motion tabs)
            Button("Open Output File…") {
                model.openOutputFile()
            }
            .keyboardShortcut("o", modifiers: [.command, .shift])
            .disabled(off || model.isRunning || model.isScanning || model.nuclear.isRunning)

            Divider()

            Button("Save Optimised Geometry in Molecule Menu") {
                model.saveOptimizedMolecule(select: false)
            }
            .disabled(off || model.optimizedMolecule == nil || model.optimizedMoleculeSaved)

            Button("Save Optimised Geometry as XYZ…") {
                if let m = model.optimizedMolecule { model.exportXYZ(m) }
            }
            .disabled(off || model.optimizedMolecule == nil)

            Button("Save Optimised Geometry as Gaussian Input…") {
                if let m = model.optimizedMolecule {
                    let level = model.gaussianLevel
                    model.exportGaussianInput(m, method: level.method, basis: level.basis, orbitals: level.orbitals)
                }
            }
            .disabled(off || model.optimizedMolecule == nil)

            Divider()

            Button("Export Results as Text…") {
                if model.calculation == .nuclear {
                    model.nuclear.exportText()   // 20.0: the nuclear-motion output file
                } else {
                    model.exportText()
                }
            }
            .keyboardShortcut("e", modifiers: [.command, .shift])
            .disabled(off || model.calculation == .nuclear ? model.nuclear.result == nil : !model.canExportOrCopy)

            // 21.0: the eigenvalues, eigenvectors and/or Gaussian start vector (a sheet chooses)
            Button("Export Eigenvalues and Eigenvectors…") {
                model.nuclear.showExportSheet = true
            }
            .disabled(off || model.calculation != .nuclear || model.nuclear.result?.exportId == nil)

            Button("Export Results as JSON…") {
                if model.calculation == .nuclear {
                    model.nuclear.exportJSON()   // 23.0: the nuclear-motion result (its archive)
                } else {
                    model.exportJSON()
                }
            }
            .keyboardShortcut("e", modifiers: [.command, .shift, .option])
            .disabled(off || (model.calculation == .nuclear ? model.nuclear.result == nil : !model.canExportOrCopy))

            Button("Export Potential Curve as Text…") {
                model.exportCurve()
            }
            .disabled(off || model.curvePoints.isEmpty)
        }

        // View > the tabs and their pages (⌘1 … ⌘9)
        CommandGroup(before: .toolbar) {
            // 20.0: the kind of calculation (the switch at the top of the sidebar)
            Button("Electronic Structure Calculation") { model.calculation = .electronic }
                .keyboardShortcut("e", modifiers: [.command, .control])
                .disabled(off || model.isRunning || model.isScanning || model.nuclear.isRunning)
            Button("Nuclear-Motion Calculation") { model.calculation = .nuclear }
                .keyboardShortcut("n", modifiers: [.command, .control])
                .disabled(off || model.isRunning || model.isScanning || model.nuclear.isRunning)
            Divider()
            // 20.0: the electronic-structure tabs (disabled in the nuclear-motion interface)
            Group {
            Button("Summary") { model.show(.summary) }
                .keyboardShortcut("1", modifiers: [.command])
            Button("Hartree–Fock") { model.show(.methods, methods: .hartreeFock) }
                .keyboardShortcut("2", modifiers: [.command])
            Button("Coupled Cluster") { model.show(.methods, methods: .coupledCluster) }
                .keyboardShortcut("3", modifiers: [.command])
            Button("Full-CI") { model.show(.methods, methods: .fullCI) }
                .keyboardShortcut("4", modifiers: [.command])
            Button("Orbitals") { model.show(.methods, methods: .orbitals) }
                .keyboardShortcut("o", modifiers: [.command, .option])
            Button("Optimisation") { model.show(.optimization) }
                .keyboardShortcut("5", modifiers: [.command])
            Button("Frequencies") { model.show(.frequencies, vibration: .frequencies) }
                .keyboardShortcut("6", modifiers: [.command])
            Button("Thermochemistry") { model.show(.frequencies, vibration: .thermochemistry) }
                .keyboardShortcut("7", modifiers: [.command])
            Button("Relaxed Scan") { model.show(.potentialCurve, scan: .relaxed) }
                .keyboardShortcut("0", modifiers: [.command])
            Button("Potential Curve") { model.show(.potentialCurve, scan: .potentialCurve) }
                .keyboardShortcut("8", modifiers: [.command])
            Button("Output") { model.show(.output) }
                .keyboardShortcut("9", modifiers: [.command])
            }
            .disabled(off || model.calculation == .nuclear)
            Divider()
            // 20.0: the nuclear-motion tabs
            Group {
                Button("Energy Levels") { model.nuclear.selectedTab = .levels }
                Button("Averaged Geometry") { model.nuclear.selectedTab = .expectation }
                Button("Matrix Structure") { model.nuclear.selectedTab = .matrix }
                Button("Eigenvectors") { model.nuclear.selectedTab = .eigenvectors }
                Button("Convergence") { model.nuclear.selectedTab = .convergence }
                Button("Nuclear-Motion Output") { model.nuclear.selectedTab = .output }
            }
            .disabled(off || model.calculation != .nuclear)
            Divider()
        }

        // Edit > Copy Summary as Text
        CommandGroup(after: .pasteboard) {
            Button("Copy Summary as Text") {
                model.copySummary()
            }
            .keyboardShortcut("c", modifiers: [.command, .shift])
            .disabled(off || !model.canExportOrCopy)
        }
    }
}
