//
//  SidebarView.swift
//  MacQC
//
//  Input sidebar: molecule, job (ΨLab 5.0; 7.0: relaxed scan and transition-state jobs),
//  geometry (the start of an optimisation, with the optimised values and a button that copies
//  them), method and basis set (HeH⁺: STO-nG or correlation-consistent), the scanned coordinate,
//  the transition-state options, the single point of jobs 5 and 9, thermochemistry and
//  optimisation settings, the size of the calculation (each step of the job, analytic or
//  numerical, with warnings and refusals), appearance, and the Run bar.
//

import SwiftUI

struct SidebarView: View {
    @Bindable var model: AppModel
    @State private var confirmRemove = false
    @AppStorage(AppearanceChoice.storageKey) private var appearanceRaw: String = AppearanceChoice.system.rawValue

    private static let numberFormat = FloatingPointFormatStyle<Double>()
        .precision(.fractionLength(0...6))

    var body: some View {
        Form {
            Section("Atom/Molecule") {
                // the atoms H–Ar, the molecules (16.0: in the order of MoleculeMenu.molecules), the
                // kept custom molecules (by name; 9.0) and "Custom molecule…", which opens a
                // structure file
                Picker("Atom/Molecule", selection: $model.moleculeSelection) {
                    Section("Atoms") {
                        ForEach(MoleculePresets.atoms) { p in
                            Text(p.title).tag(MoleculeSelection.preset(p.key))
                        }
                    }
                    Section("Molecules") {
                        ForEach(MoleculeMenu.molecules, id: \.self) { entry in
                            Text(MoleculeMenu.title(entry)).tag(entry)
                        }
                    }
                    if !model.savedMolecules.isEmpty {
                        Divider()
                        ForEach(model.savedMolecules) { m in
                            Text(m.menuTitle).tag(MoleculeSelection.saved(m.id))
                        }
                    }
                    Divider()
                    Text("Custom molecule…").tag(MoleculeSelection.openFile)
                }
                .disabled(model.isRunning || model.isScanning)
                .accessibilityLabel("Atom/Molecule")
                if model.molecule == .custom, let m = model.customMolecule {
                    CustomMoleculeHeader(model: model, molecule: m, confirmRemove: $confirmRemove)
                } else if model.molecule == .preset, let p = model.selectedPreset {
                    Text(p.note)
                        .font(.caption)
                        .foregroundStyle(.secondary)
                } else {
                    Text(model.molecule.subtitle)
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
                // 17.0: the switch is off by default and stays visible for HeH⁺ in STO-nG either way
                if model.molecule == .heh && BasisSets.isSzaboOstlund(model.hehBasis) {
                    Toggle("Szabo–Ostlund calculation", isOn: $model.hehSzaboOstlund)
                        .disabled(model.isRunning || model.isScanning)
                        .help("On: the Hartree–Fock and Full-CI program of Szabo and Ostlund (STO-nG Slater functions). Off: HeH⁺ in STO-nG as a general molecule (every method and job, charge and multiplicity).")
                }
                if model.isGeneral {
                    spinControls
                }
            }

            if model.isGeneral {
                jobSection
                generalSections
            } else {
                hehSections
            }

            Section("Appearance") {
                Picker("Appearance", selection: $appearanceRaw) {
                    ForEach(AppearanceChoice.allCases) { choice in
                        Text(choice.title).tag(choice.rawValue)
                    }
                }
                .pickerStyle(.segmented)
                .labelsHidden()
            }
        }
        .formStyle(.grouped)
        .safeAreaInset(edge: .bottom) {
            runBar
        }
        .task(id: model.estimateKey) {
            // Size and time of the calculation, refreshed shortly after an input changes.
            try? await Task.sleep(for: .milliseconds(200))
            guard !Task.isCancelled else { return }
            await model.refreshEstimate()
        }
        .alert("Structure file", isPresented: $model.showStructureError) {
            Button("OK", role: .cancel) {}
        } message: {
            Text(model.structureError ?? "The file could not be read.")
        }
        .confirmationDialog("Remove \(model.customMolecule?.menuTitle ?? "the molecule") from the molecule menu?",
                            isPresented: $confirmRemove) {
            Button("Remove", role: .destructive) { model.removeCustomMolecule() }
            Button("Cancel", role: .cancel) {}
        } message: {
            Text("Its coordinates are deleted from ΨLab's list (the structure file is not changed).")
        }
        .confirmationDialog("Long calculation", isPresented: $model.showLongRunConfirmation) {
            Button("Run") {
                Task { await model.run() }
            }
            Button("Cancel", role: .cancel) {}
        } message: {
            Text(model.longRunMessage + "\n\nYou can cancel it while it runs.")
        }
    }

    // MARK: HeH+ (Szabo–Ostlund)

    @ViewBuilder
    private var hehSections: some View {
        Section("Geometry") {
            LabeledContent("R (Å)") {
                numberField("R angstrom", value: $model.rAngstrom)
            }
        }

        Section("Basis") {
            basisPicker
            Text("STO-nG: the Slater functions of Szabo and Ostlund, Hartree–Fock and Full-CI. Any other basis set (6-31G … aug-cc-pCV5Z), or STO-nG with the Szabo–Ostlund calculation switched off, gives HeH⁺ every method and job type.")
                .font(.caption)
                .foregroundStyle(.secondary)
            LabeledContent("ζ (He)") {
                numberField("zeta He", value: $model.zeta1)
            }
            LabeledContent("ζ (H)") {
                numberField("zeta H", value: $model.zeta2)
            }
        }

        Section("Options") {
            Toggle("Exact erf (std::erf)", isOn: $model.exactErf)
                .help("Use the C++ standard-library erf instead of the Abramowitz–Stegun approximation of the original FORTRAN program.")
        }

        Section {
            Button("Reset to Szabo–Ostlund defaults") {
                model.resetToDefaults()
            }
            .secondaryActionStyle()
            .accessibilityLabel("Reset to Szabo–Ostlund defaults")
            .frame(maxWidth: .infinity)
            .disabled(model.isRunning)
        }
    }

    /// Basis-set menu of the selected molecule (HeH⁺: STO-1G, STO-2G, STO-3G first), 14.0: in
    /// sections by family (6-31G, 6-311G, def2, aug-cc-pVXZ, aug-cc-pCVXZ).
    private var basisPicker: some View {
        // 15.0: STO-1G, STO-2G and STO-3G are a family of the library (every molecule); HeH⁺ in
        // STO-nG is the Szabo–Ostlund calculation unless that option is off
        Picker("Basis set", selection: $model.currentBasis) {
            basisFamilySections(value: { $0 })
        }
        .disabled(model.isRunning || model.isScanning)
    }

    /// The basis sets of the library in one section per family, each tagged with `value(name)`.
    @ViewBuilder
    private func basisFamilySections<T: Hashable>(value: @escaping (String) -> T) -> some View {
        ForEach(BasisSets.families, id: \.self) { family in
            Section(family.title) {
                ForEach(family.names, id: \.self) { name in
                    Text(name).tag(value(name))
                }
            }
        }
    }

    /// The note under the basis-set menu of the general path (14.0: the Pople and def2 sets).
    private var basisCaption: String {
        let kind = BasisSets.isCartesianD(model.currentBasis)
            ? "\(model.currentBasis): Cartesian d functions (6D), as Gaussian uses the 6-31G family; the other sets are spherical (5D, 7F)."
            : "Gaussian functions from basissetexchange.org, spherical (5D, 7F) except the Cartesian d functions (6D) of the 6-31G family."
        switch model.molecule {
        case .heh:
            return kind + " He is not in the 6-311++G sets and has no core-valence sets (they fall back to aug-cc-pVXZ)."
                + (BasisSets.isSzaboOstlund(model.currentBasis) ? " STO-nG: switch the Szabo–Ostlund calculation on (Molecule) for the original HF and Full-CI program." : "")
        case .preset:
            return kind + " STO-1G covers H–Ca, STO-2G and STO-3G H–Xe; 6-31G H–Kr; aug-cc-pVXZ H–Ar."
        case .custom:
            return kind + " Elements: STO-1G H–Ca, STO-2G and STO-3G H–Xe; def2 H–Rn with effective core potentials from Rb on (TZVPPD, QZVPPD without Ce–Lu; mTZVP, mTZVPP also Th–Lr); 6-31G, 6-31G* H–Kr; 6-31+G… H–Ar; 6-311G, 6-311G* H–Ca, Ga–Kr and I (the diffuse ones up to Ca); aug-cc-pVXZ H–Ar and Sc–Kr (Sc–Zn up to QZ), aug-cc-pCVXZ Li–Ar."
        default:
            return kind + " The core-valence sets (CV) add tight functions for the 1s core."
        }
    }

    // MARK: General path

    private var jobSection: some View {
        Section("Job") {
            LabeledContent("Job") {
                Menu {
                    Picker("Job", selection: $model.job) {
                        ForEach(JobChoice.allCases) { choice in
                            Text("\(choice.number)   \(choice.title)").tag(choice)
                        }
                    }
                    .pickerStyle(.inline)
                    .labelsHidden()
                } label: {
                    Text(model.job.shortTitle)
                }
                .fixedSize()
                .accessibilityLabel("Job")
            }
            Text(model.job.title + ": " + model.job.explanation)
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .disabled(model.isRunning || model.isScanning)
    }

    /// The molecule's default geometry: experimental for H₂, LiH and H₂O.
    private var resetTitle: String {
        switch model.molecule {
        case .h2, .lih, .h2o: return "Experimental geometry"
        case .heh: return "Szabo–Ostlund R = 0.7743 Å"
        case .h2o2, .nh3, .custom, .preset: return "Default geometry"
        }
    }

    /// Symmetry of the geometry inputs.
    private var symmetryNote: String? {
        switch model.molecule {
        case .h2o:
            return "C₂ᵥ: both O–H bonds equal."
        case .h2o2:
            return "C₂: both O–H bonds and both O–O–H angles equal; dihedral 0° is cis (C₂ᵥ), 180° trans (C₂ₕ)."
        case .nh3:
            return "C₃ᵥ: three equal N–H bonds and H–N–H angles; 120° is planar (D₃ₕ, the transition state of the inversion)."
        default:
            return nil
        }
    }

    @ViewBuilder
    private var generalSections: some View {
        if model.molecule == .preset {
            presetGeometrySection
        } else if model.molecule == .custom {
            Section(model.job.usesOptimizer ? "Starting geometry (Å)" : "Geometry (Å)") {
                if let m = model.customMolecule {
                    CustomCoordinateTable(molecule: m)
                } else {
                    Text("Choose Custom molecule… in the molecule menu to read a structure file.")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
                SaveOptimizedGeometryBox(model: model)
            }
        } else {
            builtInGeometrySection
        }
        generalSectionsBelowGeometry
    }

    /// 16.0: the geometry parameters of a molecule of the menu (as for H₂O), none for an atom.
    @ViewBuilder
    private var presetGeometrySection: some View {
        let key = model.presetKey
        let params = PresetGeometry.parameters(key)
        Section(params.isEmpty || !model.job.usesOptimizer ? "Geometry" : "Starting geometry") {
            if params.isEmpty {
                Text("\(model.selectedPreset?.note ?? "A single atom"): no geometry parameters.")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            } else {
                ForEach(params.indices, id: \.self) { i in
                    LabeledContent("\(params[i].label) (\(params[i].unit))") {
                        numberField(params[i].label,
                                    value: Binding(get: { model.presetParameter(i) }, set: { model.setPresetParameter(i, $0) }))
                    }
                }
                if let note = PresetGeometry.note(key) {
                    Text(note)
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
                Button("Default geometry") {
                    model.resetPresetGeometry()
                }
                .secondaryActionStyle()
                .accessibilityLabel("Default geometry")
                .frame(maxWidth: .infinity)
                .disabled(model.isRunning || !model.presetGeometryChanged)
                if let o = model.optimizedPresetParameters {
                    VStack(alignment: .leading, spacing: 6) {
                        Text((o.transitionState ? "Transition state: " : "Optimised: ") + PresetGeometry.describe(key, o.values))
                            .font(.caption.monospacedDigit())
                            .textSelection(.enabled)
                        Button(o.transitionState ? "Use transition-state geometry" : "Use optimised geometry") {
                            model.useOptimizedPresetGeometry()
                        }
                        .secondaryActionStyle()
                        .accessibilityLabel("Use optimised geometry")
                        .frame(maxWidth: .infinity)
                        .disabled(!model.canUseOptimizedPresetGeometry || model.isRunning)
                    }
                }
            }
            SaveOptimizedGeometryBox(model: model)
        }
    }

    @ViewBuilder
    private var builtInGeometrySection: some View {
        Section(model.job.usesOptimizer ? "Starting geometry" : "Geometry") {
            LabeledContent("\(model.molecule.bondLabel) (Å)") {
                numberField("bond length angstrom", value: $model.bondLengthAngstrom)
            }
            if model.molecule == .h2o2 {
                LabeledContent("R(O–O) (Å)") {
                    numberField("O–O bond length angstrom", value: $model.ooBondLengthAngstrom)
                }
            }
            if let label = model.molecule.angleLabel {
                LabeledContent("\(label) (°)") {
                    numberField("angle", value: $model.angle)
                }
            }
            if model.molecule == .h2o2 {
                LabeledContent("Dihedral H–O–O–H (°)") {
                    numberField("dihedral", value: $model.dihedral)
                }
            }
            if let note = symmetryNote {
                Text(note)
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
            Button(resetTitle) {
                model.resetGeometry()
            }
            .secondaryActionStyle()
            .accessibilityLabel(resetTitle)
            .frame(maxWidth: .infinity)
            .disabled(model.isRunning)
            if let g = model.optimizedGeometry {
                let ts = model.result?.optimization?.isTransitionState == true
                VStack(alignment: .leading, spacing: 6) {
                    Text((ts ? "Transition state: " : "Optimised: ") + optimizedText(g))
                        .font(.caption.monospacedDigit())
                        .textSelection(.enabled)
                    Button(ts ? "Use transition-state geometry" : "Use optimised geometry") {
                        model.useOptimizedGeometry()
                    }
                    .secondaryActionStyle()
                    .accessibilityLabel("Use optimised geometry")
                    .frame(maxWidth: .infinity)
                    .disabled(!model.canUseOptimizedGeometry || model.isRunning)
                }
            }
            SaveOptimizedGeometryBox(model: model)
        }
    }

    @ViewBuilder
    private var generalSectionsBelowGeometry: some View {
        if model.job.hasScan {
            scanSection
        }
        if model.job.hasTransitionState {
            transitionStateSection
        }
        if (model.job.hasFrequencies || model.job.hasTransitionState) && model.hasAnalyticCorrelatedHessian {
            hessianSection
        }

        Section("Method and basis set") {
            Picker("Method", selection: $model.method) {
                ForEach(MethodChoice.families) { family in
                    Section(family.title) {
                        ForEach(family.methods) { choice in
                            Text(choice.title).tag(choice)
                        }
                    }
                }
            }
            .pickerStyle(.menu)
            Text(model.method.explanation + (model.method == .hf ? "" : " All electrons are correlated."))
                .font(.caption)
                .foregroundStyle(.secondary)
            // 15.0: the reference (RHF, UHF, ROHF)
            Picker("Reference", selection: $model.reference) {
                ForEach(ReferenceChoice.allCases) { r in
                    Text(r.title).tag(r)
                }
            }
            .pickerStyle(.menu)
            Text(referenceText)
                .font(.caption)
                .foregroundStyle(.secondary)
            basisPicker
            Text(basisCaption)
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .disabled(model.isRunning || model.isScanning)

        if model.job.hasSinglePoint {
            Section("Single-point energy") {
                Picker("SP method", selection: $model.spMethod) {
                    Text("Same (\(model.method.title))").tag(MethodChoice?.none)
                    ForEach(MethodChoice.families) { family in
                        Section(family.title) {
                            ForEach(family.methods) { choice in
                                Text(choice.title).tag(MethodChoice?.some(choice))
                            }
                        }
                    }
                }
                .pickerStyle(.menu)
                Picker("SP basis set", selection: $model.spBasis) {
                    Text("Same (\(model.currentBasis))").tag(String?.none)
                    basisFamilySections(value: { String?.some($0) })
                }
                .pickerStyle(.menu)
                Text((model.job.hasTransitionState ? "Energy at the transition state" : "Energy at the optimised geometry")
                     + "; E + ZPE, U, H and G add the thermal corrections of the frequency step.")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
            .disabled(model.isRunning || model.isScanning)
        }

        if model.job.hasFrequencies {
            Section("Thermochemistry") {
                LabeledContent("Temperature (K)") {
                    numberField("temperature", value: $model.temperature)
                }
                LabeledContent("Pressure (Pa)") {
                    numberField("pressure", value: $model.pressure)
                }
                HStack {
                    Button("1 atm") { model.pressure = 101325 }
                    Button("1 bar") { model.pressure = 100000 }
                    Spacer()
                }
                .controlSize(.small)
                Text("Ideal gas, rigid rotor, harmonic oscillators; CODATA 2018 constants, AME 2020 isotope masses.")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
            .disabled(model.isRunning || model.isScanning)
        }

        if model.job.usesOptimizer {
            Section(model.job.hasScan ? "Optimisation at each point" : "Optimisation") {
                Picker("Convergence", selection: $model.tightConvergence) {
                    Text("Normal").tag(false)
                    Text("Tight").tag(true)
                }
                .pickerStyle(.segmented)
                Stepper(value: $model.maxSteps, in: 5...500, step: 5) {
                    LabeledContent("Maximum steps", value: "\(model.maxSteps)")
                }
                Text(model.tightConvergence
                     ? "Max force 1.5e-5, RMS force 1.0e-5 E_h/bohr; max step 6.0e-5, RMS step 4.0e-5 bohr; energy change 1e-8 E_h."
                     : "Max force 4.5e-4, RMS force 3.0e-4 E_h/bohr; max step 1.8e-3, RMS step 1.2e-3 bohr; energy change 1e-6 E_h.")
                    .font(.caption)
                    .foregroundStyle(.secondary)
                if model.job.hasOptimization && !model.job.hasTransitionState {
                    optimizerControls
                }
                if model.hasAnalyticCorrelatedHessian && !model.job.hasTransitionState && model.optHessianMode == .automatic {
                    Toggle("Analytic initial Hessian", isOn: $model.analyticInitialHessian)
                    Text(model.analyticInitialHessian
                         ? "The first step uses the analytic \(model.method.title) Hessian (exact curvature; Lindh's model Hessian if it does not fit in memory); later steps update it (BFGS)."
                         : "The first step uses Lindh's model Hessian (no cost); later steps update it (BFGS).")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
            }
            .disabled(model.isRunning || model.isScanning)
        }

        Section("Size of the calculation") {
            EstimatePanel(model: model)
        }
    }

    // MARK: Relaxed scan and transition-state search (7.0)

    private var scanSection: some View {
        let s = model.scan
        let positions = model.estimate?.molecule.map { atomCoordinates($0.atoms) }
        let value = positions.flatMap { internalCoordinateValue($0, kind: s.kind, atoms: s.atoms) }
        let labels = model.atomLabels
        let names: [String] = s.usedAtoms.map { n in n >= 1 && n <= labels.count ? labels[n - 1] : "?\(n)" }
        let definition = "\(s.kind.rawValue)(\(names.joined(separator: "–")))"
        let unit = s.kind.unit
        let start: Double? = s.fromInput ? value : s.start
        var current = definition
        if let v = value {
            current += " = \(valueText(v, unit: unit)) at the input geometry"
        }
        var range = ""
        if let a = start {
            let end = a + Double(s.steps) * s.step
            range = "From \(valueText(a, unit: unit)) to \(valueText(end, unit: unit)). Gaussian input: \(s.line)"
        }
        return Section("Scanned coordinate") {
            Menu("Suggested scans") {
                ForEach(model.molecule.scanSuggestions) { suggestion in
                    Button(suggestion.title) {
                        model.scan = ScanSettings(suggestion)
                    }
                }
            }
            .accessibilityLabel("Suggested scans")
            Picker("Type", selection: $model.scan.kind) {
                ForEach(ScanKind.allCases) { kind in
                    Text(kind.title).tag(kind)
                }
            }
            .pickerStyle(.segmented)
            LabeledContent("Atoms") {
                HStack(spacing: 4) {
                    ForEach(0..<s.kind.atomCount, id: \.self) { k in
                        Picker("Atom \(k + 1)", selection: $model.scan.atoms[k]) {
                            ForEach(1...max(model.atomCount, 1), id: \.self) { n in
                                Text(n <= labels.count ? labels[n - 1] : "\(n)").tag(n)
                            }
                        }
                        .labelsHidden()
                        .fixedSize()
                    }
                }
            }
            Text(current)
                .font(.caption.monospacedDigit())
                .foregroundStyle(.secondary)
            Toggle("Start at the input geometry", isOn: $model.scan.fromInput)
            if !s.fromInput {
                LabeledContent("Start (\(unit))") {
                    numberField("scan start", value: $model.scan.start)
                }
            }
            Stepper(value: $model.scan.steps, in: 1...360) {
                LabeledContent("Steps", value: "\(s.steps) (\(s.steps + 1) points)")
            }
            LabeledContent("Step size (\(unit))") {
                numberField("scan step", value: $model.scan.step)
            }
            if !range.isEmpty {
                Text(range)
                    .font(.caption.monospacedDigit())
                    .foregroundStyle(.secondary)
                    .textSelection(.enabled)
            }
            if let problem = model.scanProblem {
                Label(problem, systemImage: "xmark.octagon.fill")
                    .font(.caption)
                    .foregroundStyle(.red)
            } else {
                Text(model.atomCount == 2
                     ? "A diatomic has nothing else to relax: the scan computes single points."
                     : "At each point the coordinate is frozen and all others are optimised, starting from the previous point.")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
        }
        .disabled(model.isRunning || model.isScanning)
    }

    private func valueText(_ v: Double, unit: String) -> String {
        unit == "Å" ? String(format: "%.4f Å", v) : String(format: "%.2f°", v)
    }

    /// ΨLab 10.0: analytic or numerical Hessians of Full CI, CCSD and CCSD(T).
    private var hessianSection: some View {
        Section("Hessian") {
            Picker("Second derivatives", selection: $model.analyticHessian) {
                Text("Analytic").tag(true)
                Text("Numerical").tag(false)
            }
            .pickerStyle(.segmented)
            .accessibilityLabel("Hessian")
            Text(model.analyticHessian ? analyticHessianText : numericalHessianText)
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .disabled(model.isRunning || model.isScanning)
    }

    private var analyticHessianText: String {
        let what = model.method == .fci ? "CI-vector" : "amplitude"
        let reference = model.effectiveReference
        var orbitals = "orbital response by CPHF"
        if reference == .uhf || reference == .rohf {
            orbitals = model.method == .fci
                ? "no orbital response (the Full-CI energy does not depend on the orbitals)"
                : "orbital response by coupled-perturbed \(reference.title)"
        }
        return "Analytic: \(what) responses along the internal directions and to an electric field (IR intensities), "
            + "\(orbitals); exact second derivatives. Numerical only if it does not fit in memory."
    }

    private var numericalHessianText: String {
        "Numerical: central differences of analytic gradients (step 0.005 bohr), IR intensities from the relaxed "
            + "dipole moments."
    }

    private var transitionStateSection: some View {
        Section("Transition-state search") {
            Picker("Hessian", selection: $model.tsHessianMode) {
                Text("CalcFC").tag(HessianMode.calcFC)
                Text("CalcAll").tag(HessianMode.calcAll)
                Text("RecalcFC").tag(HessianMode.recalcFC)
            }
            .pickerStyle(.segmented)
            if model.tsHessianMode == .recalcFC {
                Stepper(value: $model.recalcFC, in: 1...100) {
                    LabeledContent("Every N steps", value: "\(model.recalcFC)")
                }
            }
            Text(model.tsHessianMode == .calcAll
                 ? "CalcAll: the Hessian is computed at every step (exact, more expensive)."
                 : model.tsHessianMode == .recalcFC
                    ? "RecalcFC=\(model.recalcFC): the Hessian is computed at step 1 and every \(model.recalcFC) steps, updated (Bofill) in between."
                    : "CalcFC: the Hessian is computed at the first step, then updated (Bofill).")
                .font(.caption)
                .foregroundStyle(.secondary)
            Toggle("Eigenvalue test", isOn: $model.eigenTest)
            Text(model.eigenTest
                 ? "The search stops unless the Hessian has exactly one negative eigenvalue (Gaussian's default)."
                 : "No eigenvalue test (Gaussian's NoEigenTest): the search goes on with any number of negative eigenvalues.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .disabled(model.isRunning || model.isScanning)
    }

    private func optimizedText(_ g: GeometryInputs) -> String {
        var text = String(format: "%@ = %.6f Å", model.molecule.bondLabel, g.R * Units.bohrToAngstrom)
        if let roo = g.ROO {
            text += String(format: ", R(O–O) = %.6f Å", roo * Units.bohrToAngstrom)
        }
        if let a = g.angle, model.molecule.angleLabel != nil {
            text += String(format: ", angle %.4f°", a)
        }
        if let d = g.dihedral {
            text += String(format: ", dihedral %.3f°", d)
        }
        if !g.symmetric {
            text += " (not symmetric: means)"
        }
        return text
    }

    private func numberField(_ title: String, value: Binding<Double>) -> some View {
        TextField(title, value: value, format: Self.numberFormat)
            .labelsHidden()
            .multilineTextAlignment(.trailing)
            .monospacedDigit()
            .frame(width: 96)
    }

    // MARK: Charge, multiplicity, reference and optimiser (15.0)

    /// Charge and multiplicity of the selected molecule (steppers), with the electron count.
    @ViewBuilder
    private var spinControls: some View {
        Stepper(value: $model.charge, in: -10...20) {
            LabeledContent("Charge", value: model.charge > 0 ? "+\(model.charge)" : "\(model.charge)")
        }
        .disabled(model.isRunning || model.isScanning)
        .accessibilityLabel("Charge")
        Stepper(value: $model.multiplicity, in: 1...20) {
            LabeledContent("Multiplicity", value: "\(model.multiplicity) (\(AppModel.spinName(model.multiplicity)))")
        }
        .disabled(model.isRunning || model.isScanning)
        .accessibilityLabel("Multiplicity")
        if let problem = model.spinProblem {
            Label(problem, systemImage: "xmark.octagon.fill")
                .font(.caption)
                .foregroundStyle(.red)
        } else {
            Text("\(model.electronCount) electrons (\((model.electronCount + model.multiplicity - 1) / 2) α, \((model.electronCount - model.multiplicity + 1) / 2) β)")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        if model.charge != model.defaultCharge || model.multiplicity != model.defaultMultiplicity {
            Button("Charge \(model.defaultCharge > 0 ? "+" : "")\(model.defaultCharge), multiplicity \(model.defaultMultiplicity) (the molecule's)") {
                model.charge = model.defaultCharge
                model.multiplicity = model.defaultMultiplicity
            }
            .secondaryActionStyle()
            .controlSize(.small)
            .disabled(model.isRunning || model.isScanning)
        }
    }

    private var referenceText: String {
        if model.reference == .automatic {
            let open = model.multiplicity != 1
            return ReferenceChoice.automatic.explanation + " Here: " + (open ? (model.method == .fci ? "ROHF" : "UHF") : "RHF") + "."
        }
        return model.reference.explanation
    }

    /// Computed Hessians (Gaussian's CalcFC, CalcAll, RecalcFC=N) and the step (RFO or GEDIIS) of
    /// minimisations.
    @ViewBuilder
    private var optimizerControls: some View {
        Picker("Hessian", selection: $model.optHessianMode) {
            Text("Initial").tag(HessianMode.automatic)
            Text("CalcFC").tag(HessianMode.calcFC)
            Text("CalcAll").tag(HessianMode.calcAll)
            Text("RecalcFC").tag(HessianMode.recalcFC)
        }
        .pickerStyle(.segmented)
        if model.optHessianMode == .recalcFC {
            Stepper(value: $model.recalcFC, in: 1...100) {
                LabeledContent("Every N steps", value: "\(model.recalcFC)")
            }
        }
        Text(optHessianText)
            .font(.caption)
            .foregroundStyle(.secondary)
        Picker("Step", selection: $model.gediis) {
            Text("RFO").tag(false)
            Text("GEDIIS").tag(true)
        }
        .pickerStyle(.segmented)
        Text(model.gediis
             ? "Hybrid GEDIIS/RFO (Li and Frisch 2006, Gaussian's default): once the RMS force is below 0.01, the energy-represented DIIS interpolation of the last five points, then an RFO step from there; RFO when it is not downhill or too long."
             : "Rational-function optimisation (RFO) steps with a trust radius in redundant internal coordinates.")
            .font(.caption)
            .foregroundStyle(.secondary)
    }

    private var optHessianText: String {
        switch model.optHessianMode {
        case .automatic:
            return "The initial Hessian (analytic where possible, else Lindh's model), then BFGS updates."
        case .calcFC:
            return "CalcFC: the Hessian of the method is computed at the first step (analytic where available), then BFGS updates."
        case .calcAll:
            return "CalcAll: the Hessian is computed at every step (Newton steps; expensive)."
        case .recalcFC:
            return "RecalcFC=\(model.recalcFC): the Hessian is computed at step 1 and every \(model.recalcFC) steps, BFGS updates in between."
        }
    }

    // MARK: Run bar

    private var runBar: some View {
        GlassGroup(spacing: 10) {
            VStack(spacing: 8) {
                HStack(spacing: 8) {
                    Button {
                        model.requestRun()
                    } label: {
                        HStack(spacing: 8) {
                            if model.isRunning {
                                ProgressView()
                                    .controlSize(.small)
                            } else {
                                Image(systemName: "play.fill")
                            }
                            Text(model.isRunning ? "Running…" : "Run")
                                .fontWeight(.semibold)
                        }
                        .frame(maxWidth: .infinity)
                    }
                    .prominentActionStyle()
                    .accessibilityLabel("Run")
                    .controlSize(.large)
                    .keyboardShortcut("r", modifiers: [.command])
                    .disabled(model.isRunning || model.isScanning || model.refusal != nil || model.scanProblem != nil)

                    if model.isRunning {
                        Button("Cancel") {
                            model.cancelRun()
                        }
                        .secondaryActionStyle()
                        .controlSize(.large)
                        .accessibilityLabel("Cancel calculation")
                        .keyboardShortcut(".", modifiers: [.command])
                    }
                }

                if model.isRunning && model.isGeneral {
                    TimelineView(.periodic(from: .now, by: 1)) { context in
                        Text(model.progressDescription + "  ·  " + model.elapsedText(at: context.date))
                            .font(.caption.monospacedDigit())
                            .foregroundStyle(.secondary)
                            .frame(maxWidth: .infinity, alignment: .leading)
                            .lineLimit(2)
                    }
                }

                if let message = model.errorMessage {
                    Text(message)
                        .font(.caption)
                        .foregroundStyle(.red)
                        .frame(maxWidth: .infinity, alignment: .leading)
                        .lineLimit(4)
                } else if model.isGeneral && !model.isRunning, model.refusal != nil {
                    Label("Not possible with these settings (see Size of the calculation).", systemImage: "xmark.octagon.fill")
                        .font(.caption)
                        .foregroundStyle(.red)
                        .frame(maxWidth: .infinity, alignment: .leading)
                } else if model.isGeneral && !model.isRunning && model.isLongCalculation,
                          let seconds = model.estimatedSeconds {
                    Label("Long calculation: \(Fmt.duration(seconds)).", systemImage: "clock.badge.exclamationmark")
                        .font(.caption)
                        .foregroundStyle(.orange)
                        .frame(maxWidth: .infinity, alignment: .leading)
                }
            }
            .padding(12)
            .glassPanel(cornerRadius: 18)
        }
        .padding(.horizontal, 12)
        .padding(.bottom, 10)
    }
}

/// Basis functions, integrals, determinants, memory and time of the selected calculation,
/// with the backend's notes, warnings and refusal.
struct EstimatePanel: View {
    var model: AppModel

    var body: some View {
        if let out = model.estimate {
            // 16.0: an estimate of earlier inputs (e.g. the previous basis set) is dimmed while the
            // new one is computed
            let stale = model.estimatedKey != model.estimateKey
            VStack(alignment: .leading, spacing: 6) {
                if stale {
                    HStack(spacing: 8) {
                        ProgressView()
                            .controlSize(.small)
                        Text("Estimating…")
                            .foregroundStyle(.secondary)
                    }
                }
                VStack(alignment: .leading, spacing: 6) {
                    if let b = out.basis {
                        InfoRow("Basis functions", "\(b.functions)")
                        ForEach(Fmt.uniqueAtoms(b.atoms), id: \.symbol) { atom in
                            InfoRow("  \(atom.symbol)", (atom.basis == model.currentBasis ? "" : atom.basis + " ") + atom.contraction)
                        }
                        // 16.0: another set of the menu with exactly the same functions for these elements
                        let elements = Fmt.uniqueAtoms(b.atoms).map(\.symbol)
                        let same = BasisIdentity.sameFunctions(as: out.input.basis, elements: elements)
                        if !same.isEmpty {
                            Text("Same basis functions as \(same.joined(separator: ", ")) for \(elements.joined(separator: " and ")): "
                                 + "the Basis Set Exchange defines \(same.count == 1 ? "both sets" : "these sets") identically for "
                                 + "\(elements.count == 1 ? "this element" : "these elements").")
                                .font(.caption)
                                .foregroundStyle(.secondary)
                        }
                    }
                    if let e = out.estimate {
                        InfoRow("Two-electron integrals",
                                Fmt.count(e.uniqueIntegrals) + (e.integralsStored ? " (in memory)" : " (direct)"))
                        if let f = e.fci {
                            InfoRow("Determinants", Fmt.count(f.determinants))
                            InfoRow("Full-CI memory", Fmt.memory(f.memoryMB))  // 16.0: also when not possible
                        }
                        if let c = e.coupledCluster {
                            InfoRow("\(MethodChoice.title(of: c.method)) amplitudes", Fmt.count(c.amplitudes))
                            InfoRow("Orbitals (occ./virt.)", "\(c.occupied) / \(c.virtuals)")
                            InfoRow("\(MethodChoice.title(of: c.method)) memory", Fmt.memory(c.memoryMB))
                        }
                        // 16.0: the estimated times also of impossible jobs (marked as such)
                        let impossible = model.refusal != nil
                        if let j = out.job, j.type != "energy" {
                            ForEach(j.estimate.steps) { step in
                                VStack(alignment: .leading, spacing: 2) {
                                    InfoRow(step.name, Fmt.duration(step.seconds))
                                    Text("\(MethodChoice.title(of: step.method))/\(step.basis): \(step.how)")
                                        .font(.caption)
                                        .foregroundStyle(step.analytic ? Color.secondary : Color.orange)
                                }
                            }
                            InfoRow(impossible ? "Estimated time (job; not possible)" : "Estimated time (job)",
                                    Fmt.duration(j.estimate.totalSeconds) + " (\(e.threads) threads)")
                        } else {
                            InfoRow(impossible ? "Estimated time (not possible)" : "Estimated time",
                                    Fmt.duration(e.totalSeconds) + " (\(e.threads) threads)")
                        }
                    }
                    ForEach(out.notes ?? [], id: \.self) { note in
                        Text(note)
                            .font(.caption)
                            .foregroundStyle(.secondary)
                    }
                    ForEach(warnings(out), id: \.self) { warning in
                        Label(warning, systemImage: "clock.badge.exclamationmark")
                            .font(.caption)
                            .foregroundStyle(.orange)
                    }
                    if let refusal = model.refusal {
                        Label(refusal, systemImage: "xmark.octagon.fill")
                            .font(.caption)
                            .foregroundStyle(.red)
                            .textSelection(.enabled)
                    }
                }
                .opacity(stale ? 0.4 : 1)
            }
        } else {
            HStack(spacing: 8) {
                ProgressView()
                    .controlSize(.small)
                Text("Estimating…")
                    .foregroundStyle(.secondary)
            }
        }
    }

    /// Warnings of the energy and, for jobs, of the whole job (its "long calculation" warning
    /// replaces the one of a single energy).
    private func warnings(_ out: MacQCOutput) -> [String] {
        guard let j = out.job, j.type != "energy" else { return out.warnings ?? [] }
        let own = (out.warnings ?? []).filter { !$0.hasPrefix("Long calculation") }
        return own + j.estimate.warnings
    }
}
