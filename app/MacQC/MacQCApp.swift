//
//  MacQCApp.swift
//  MacQC
//
//  App entry point and menu commands.
//

import SwiftUI

@main
struct MacQCApp: App {
    @State private var model = AppModel()

    var body: some Scene {
        WindowGroup {
            ContentView(model: model)
                .frame(minWidth: 1000, minHeight: 640)
        }
        .defaultSize(width: 1180, height: 760)
        .commands {
            // File > Open Structure File…, Open Output File… (13.0), Save Optimised Geometry as
            // XYZ… / as Gaussian Input… (13.0), Export Results as Text… (the output file with its
            // archive) / as JSON… / Potential Curve as Text… (replace the unused "New Window" item)
            CommandGroup(replacing: .newItem) {
                // 9.0: a custom molecule from a structure file (kept in the molecule menu)
                Button("Open Structure File…") {
                    model.openStructureFile()
                }
                .keyboardShortcut("o", modifiers: [.command])
                .disabled(model.isRunning || model.isScanning)

                // 13.0: an output file (archive section) or a JSON export, every tab restored
                Button("Open Output File…") {
                    model.openOutputFile()
                }
                .keyboardShortcut("o", modifiers: [.command, .shift])
                .disabled(model.isRunning || model.isScanning)

                Divider()

                Button("Save Optimised Geometry in Molecule Menu") {
                    model.saveOptimizedMolecule(select: false)
                }
                .disabled(model.optimizedMolecule == nil || model.optimizedMoleculeSaved)

                Button("Save Optimised Geometry as XYZ…") {
                    if let m = model.optimizedMolecule { model.exportXYZ(m) }
                }
                .disabled(model.optimizedMolecule == nil)

                Button("Save Optimised Geometry as Gaussian Input…") {
                    if let m = model.optimizedMolecule {
                        let level = model.gaussianLevel
                        model.exportGaussianInput(m, method: level.method, basis: level.basis, orbitals: level.orbitals)
                    }
                }
                .disabled(model.optimizedMolecule == nil)

                Divider()

                Button("Export Results as Text…") {
                    if model.calculation == .nuclear {
                        model.nuclear.exportText()   // 20.0: the nuclear-motion output file
                    } else {
                        model.exportText()
                    }
                }
                .keyboardShortcut("e", modifiers: [.command, .shift])
                .disabled(model.calculation == .nuclear ? model.nuclear.result == nil : !model.canExportOrCopy)

                Button("Export Results as JSON…") {
                    model.exportJSON()
                }
                .keyboardShortcut("e", modifiers: [.command, .shift, .option])
                .disabled(!model.canExportOrCopy)

                Button("Export Potential Curve as Text…") {
                    model.exportCurve()
                }
                .disabled(model.curvePoints.isEmpty)
            }

            // View > the tabs and their pages (⌘1 … ⌘9)
            CommandGroup(before: .toolbar) {
                // 20.0: the kind of calculation (the switch at the top of the sidebar)
                Button("Electronic Structure Calculation") { model.calculation = .electronic }
                    .keyboardShortcut("e", modifiers: [.command, .control])
                    .disabled(model.isRunning || model.isScanning || model.nuclear.isRunning)
                Button("Nuclear-Motion Calculation") { model.calculation = .nuclear }
                    .keyboardShortcut("n", modifiers: [.command, .control])
                    .disabled(model.isRunning || model.isScanning || model.nuclear.isRunning)
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
                .disabled(model.calculation == .nuclear)
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
                .disabled(model.calculation != .nuclear)
                Divider()
            }

            // Edit > Copy Summary as Text
            CommandGroup(after: .pasteboard) {
                Button("Copy Summary as Text") {
                    model.copySummary()
                }
                .keyboardShortcut("c", modifiers: [.command, .shift])
                .disabled(!model.canExportOrCopy)
            }
        }
    }
}
