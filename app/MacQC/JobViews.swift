//
//  JobViews.swift
//  MacQC
//
//  ΨLab 5.0 job results: the Optimisation tab (the molecule of every step, drawn and played step
//  by step (6.0; 12.0: in the 3D viewer), the steps table, energy and force charts, the geometry and internal coordinates
//  of each step, the optimised geometry with a button that copies it into the inputs; 7.0: the
//  transition-state search with the Hessian's eigenvalues and the followed mode of every step),
//  the Frequencies tab with its Frequencies and Thermochemistry pages (every normal mode animated
//  (6.0; 12.0: the selected mode in the 3D viewer), the table of modes, the IR spectrum; 7.0: the point group and the transition-state
//  check), and the job panel of the Summary tab.
//

import SwiftUI
import Charts

// MARK: - Optimisation

struct OptimizationView: View {
    @Bindable var model: AppModel
    @State private var selectedStep: Int?
    /// Play steps is on (initially off; on with the launch argument -psilabPlaySteps YES, which the
    /// UI smoke test uses: clicking the button by script takes too long on this page).
    @State private var playing = UserDefaults.standard.bool(forKey: "psilabPlaySteps")
    @State private var showStart = true

    var body: some View {
        if let out = model.result, let opt = out.optimization {
            ScrollView {
                VStack(alignment: .leading, spacing: 18) {
                    if let message = model.errorMessage {
                        ErrorBanner(message: message)
                    }
                    header(out, opt)
                    if !opt.history.isEmpty {
                        moleculePanel(opt)
                    }
                    resultPanel(out, opt)

                    SectionTitle("Steps")
                    stepsTable(opt)
                    charts(opt)

                    if let step = shownStep(opt) {
                        SectionTitle("Step \(step.step)")
                        stepDetail(step, opt: opt)
                    }
                }
                .padding(20)
            }
            .task(id: playing) {
                await playSteps()
            }
        } else {
            JobPlaceholder(model: model, title: "No optimisation", systemImage: "arrow.down.right.and.arrow.up.left",
                           text: "Choose a geometry optimisation or a transition-state optimisation (with or without frequencies) in the sidebar and press Run. The geometry inputs are the starting geometry.")
        }
    }

    private func shownStep(_ opt: MacQCOptimization) -> MacQCOptimizationStep? {
        if let s = selectedStep, let step = opt.history.first(where: { $0.step == s }) { return step }
        return opt.history.last
    }

    /// Position of the shown step in the history (the last step unless one was chosen).
    private func shownIndex(_ opt: MacQCOptimization) -> Int {
        if let s = selectedStep, let i = opt.history.firstIndex(where: { $0.step == s }) { return i }
        return max(opt.history.count - 1, 0)
    }

    private func select(_ index: Int, in opt: MacQCOptimization) {
        guard !opt.history.isEmpty else { return }
        playing = false
        selectedStep = opt.history[min(max(index, 0), opt.history.count - 1)].step
    }

    /// While Play is on: the next step every 0.7 s, the last one held for 1.6 s, then from the start.
    private func playSteps() async {
        guard playing else { return }
        while !Task.isCancelled {
            guard let opt = model.result?.optimization, opt.history.count > 1 else {
                playing = false
                return
            }
            let index = shownIndex(opt)
            let last = opt.history.count - 1
            try? await Task.sleep(for: .milliseconds(index == last ? 1600 : 700))
            if Task.isCancelled || !playing { return }
            selectedStep = opt.history[index == last ? 0 : index + 1].step
        }
    }

    /// The molecule at the shown step (to the scale of the whole optimisation), its energy, forces
    /// and geometry, and the buttons and slider that step through the optimisation or play it.
    private func moleculePanel(_ opt: MacQCOptimization) -> some View {
        let history = opt.history
        let index = shownIndex(opt)
        let step = history[index]
        // The start and every step, rotated onto each other (the backend's coordinates may turn).
        let aligned = alignedSequence(([opt.initialGeometry.atoms] + history.map { $0.geometry.atoms }).map { atomCoordinates($0) },
                                      atoms: step.geometry.atoms)
        let shown = atomsMoved(step.geometry.atoms, to: aligned[min(index + 1, aligned.count - 1)])
        let start = atomsMoved(opt.initialGeometry.atoms, to: aligned[0])
        return VStack(alignment: .leading, spacing: 10) {
            HStack(alignment: .top, spacing: 16) {
                // ΨLab 12.0: the 3D viewer (turn with the mouse, click atoms to measure)
                Molecule3DView(atoms: shown,
                               reference: showStart && index > 0 ? start : nil,
                               framing: aligned,
                               accessibilityText: "Molecule at this optimisation step")
                    .frame(minWidth: 260, maxWidth: .infinity, minHeight: 300, maxHeight: 300)
                VStack(alignment: .leading, spacing: 6) {
                    Text("Step \(step.step) of \(history.count)")
                        .font(.headline)
                    Grid(alignment: .leading, horizontalSpacing: 12, verticalSpacing: 4) {
                        GridRow {
                            Text("Energy").foregroundStyle(.secondary)
                            Text(Fmt.energy10(step.energy) + " E_h")
                        }
                        GridRow {
                            Text("ΔE").foregroundStyle(.secondary)
                            Text(step.deltaE.map { String(format: "%+.2e E_h", $0) } ?? "–")
                        }
                        GridRow {
                            Text("Max force").foregroundStyle(.secondary)
                            Text(step.criterion(0).map { String(format: "%.2e E_h/bohr", $0) } ?? "–")
                        }
                        if let ts = step.ts {
                            GridRow {
                                Text("Negative eigenvalues").foregroundStyle(.secondary)
                                Text("\(ts.negativeEigenvalues)")
                                    .foregroundStyle(ts.negativeEigenvalues == 1 ? Color.green : Color.orange)
                            }
                        }
                        ForEach(Array(Fmt.geometryRows(step.geometry).enumerated()), id: \.offset) { _, row in
                            GridRow {
                                Text(row.0).foregroundStyle(.secondary)
                                Text(row.1)
                            }
                        }
                    }
                    .monospacedDigit()
                    .textSelection(.enabled)
                    Label(step.converged ? "Converged" : "Not converged yet",
                          systemImage: step.converged ? "checkmark.circle.fill" : "circle")
                        .foregroundStyle(step.converged ? Color.green : Color.secondary)
                        .font(.callout)
                }
                .frame(width: 260, alignment: .leading)
            }
            HStack(spacing: 8) {
                Button {
                    select(0, in: opt)
                } label: {
                    Image(systemName: "backward.end.fill")
                }
                .help("First step")
                .accessibilityLabel("First step")
                Button {
                    select(index - 1, in: opt)
                } label: {
                    Image(systemName: "backward.frame.fill")
                }
                .help("Previous step")
                .accessibilityLabel("Previous step")
                .disabled(index == 0)
                Button {
                    playing.toggle()
                } label: {
                    Label(playing ? "Pause" : "Play steps", systemImage: playing ? "pause.fill" : "play.fill")
                }
                .help("Show the steps one after the other")
                .disabled(history.count < 2)
                Button {
                    select(index + 1, in: opt)
                } label: {
                    Image(systemName: "forward.frame.fill")
                }
                .help("Next step")
                .accessibilityLabel("Next step")
                .disabled(index == history.count - 1)
                Button {
                    select(history.count - 1, in: opt)
                } label: {
                    Image(systemName: "forward.end.fill")
                }
                .help("Last step")
                .accessibilityLabel("Last step")
                Slider(value: Binding(get: { Double(index) },
                                      set: { select(Int($0.rounded()), in: opt) }),
                       in: 0...Double(max(history.count - 1, 1)), step: 1)
                    .disabled(history.count < 2)
                    .accessibilityLabel("Optimisation step")
                Toggle("Starting geometry", isOn: $showStart)
                    .toggleStyle(.checkbox)
                    .help("Show the geometry of step 1 as a translucent molecule")
            }
            Text("The molecule of each optimisation step in 3D, the steps turned onto each other so that the molecule keeps its place: bond lengths in Å; translucent: the starting geometry. Drag to turn the molecule, pinch or ⌥-scroll to zoom; click two, three or four atoms to measure a distance, an angle or a dihedral angle. Choose a step with the buttons, the slider or a step number in the table, or press Play steps.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }

    private func header(_ out: MacQCOutput, _ opt: MacQCOptimization) -> some View {
        VStack(alignment: .leading, spacing: 4) {
            Text("\(out.moleculeTitle) · \(MethodChoice.title(of: out.input.method))/\(out.input.basis) · "
                 + (opt.isTransitionState ? "transition-state optimisation" : "geometry optimisation"))
                .font(.title2.weight(.semibold))
            Label(opt.converged
                  ? (opt.isTransitionState ? "Transition state (saddle point of order one) found after \(opt.steps) steps"
                                           : "Stationary point found after \(opt.steps) steps")
                  : opt.message,
                  systemImage: opt.converged ? "checkmark.seal.fill" : "exclamationmark.triangle.fill")
                .foregroundStyle(opt.converged ? Color.green : Color.orange)
        }
    }

    private func resultPanel(_ out: MacQCOutput, _ opt: MacQCOptimization) -> some View {
        VStack(alignment: .leading, spacing: 10) {
            Grid(alignment: .trailing, horizontalSpacing: 24, verticalSpacing: 5) {
                GridRow {
                    Text("Parameter").gridColumnAlignment(.leading)
                    Text("Start")
                    Text(opt.converged ? (opt.isTransitionState ? "Transition state" : "Optimised") : "Last step")
                    Text("Change")
                }
                .font(.subheadline.weight(.semibold))
                .foregroundStyle(.secondary)
                Divider().gridCellUnsizedAxes(.horizontal)
                ForEach(Array(opt.finalParameters.enumerated()), id: \.offset) { index, p in
                    let start = index < opt.initialParameters.count ? opt.initialParameters[index].value : p.value
                    GridRow {
                        Text(p.name)
                        Text(Fmt.parameterValue(start, angle: p.angle))
                        Text(Fmt.parameterValue(p.value, angle: p.angle))
                        Text(Fmt.parameterChange(p.value - start, angle: p.angle))
                    }
                }
                Divider().gridCellUnsizedAxes(.horizontal)
                GridRow {
                    Text("Energy (E_h)")
                    Text(Fmt.energy10(opt.initialEnergy))
                    Text(Fmt.energy10(opt.finalEnergy))
                    Text(String(format: "%+.10f", opt.finalEnergy - opt.initialEnergy))
                }
            }
            .monospacedDigit()
            .textSelection(.enabled)

            VStack(alignment: .leading, spacing: 4) {
                if out.input.molecule == "H2O" {
                    Text(opt.symmetricEnd
                         ? String(format: "C₂ᵥ symmetry kept (|R₁ − R₂| = %.1e bohr); all three internal coordinates were optimised without constraints.", opt.symmetryDeviation)
                         : String(format: "The optimised geometry is not C₂ᵥ (|R₁ − R₂| = %.1e bohr).", opt.symmetryDeviation))
                }
                if opt.isTransitionState {
                    Text(tsSummary(opt))
                }
                Text("Gradient: " + (opt.gradient ?? "–"))
                Text("Initial Hessian: \(opt.initialHessian); " + String(format: "%.1f s", opt.seconds))
                Text(String(format: "Criteria (%@): max force %.1e, RMS force %.1e E_h/bohr; max step %.1e, RMS step %.1e bohr; energy change %.0e E_h",
                            out.input.optConvergence, opt.criteria.maxForce, opt.criteria.rmsForce,
                            opt.criteria.maxDisplacement, opt.criteria.rmsDisplacement, opt.criteria.energy))
            }
            .font(.caption)
            .foregroundStyle(.secondary)

            if model.optimizedGeometry != nil {
                HStack {
                    Button {
                        model.useOptimizedGeometry()
                    } label: {
                        Label("Use optimised geometry as input", systemImage: "arrow.uturn.backward")
                    }
                    .secondaryActionStyle()
                    .disabled(!model.canUseOptimizedGeometry || model.isRunning)
                    if !model.canUseOptimizedGeometry {
                        Text("The geometry inputs already hold the optimised geometry.")
                            .font(.caption)
                            .foregroundStyle(.secondary)
                    }
                }
            }
            SaveOptimizedGeometryBox(model: model)
        }
        .dataPanel()
    }

    /// The Hessian, its evaluations and the eigenvalue test of a transition-state search.
    private func tsSummary(_ opt: MacQCOptimization) -> String {
        var parts = [opt.tsHessian == "calcall"
                     ? "Hessian computed at every step (CalcAll)"
                     : "Hessian computed at the first step (CalcFC), then Bofill updates"]
        if let n = opt.hessianEvaluations {
            parts.append("\(n) Hessian evaluation\(n == 1 ? "" : "s")")
        }
        if let n = opt.negativeEigenvalues {
            parts.append("\(n) negative eigenvalue\(n == 1 ? "" : "s") at the end")
        }
        if opt.eigenTestFailed == true {
            parts.append("the eigenvalue test failed (the Hessian must have exactly one negative eigenvalue)")
        }
        return parts.joined(separator: "; ") + "."
    }

    private func stepsTable(_ opt: MacQCOptimization) -> some View {
        let ts = opt.isTransitionState
        return VStack(alignment: .leading, spacing: 6) {
            Grid(alignment: .trailing, horizontalSpacing: 18, verticalSpacing: 4) {
                GridRow {
                    Text("Step").gridColumnAlignment(.leading)
                    Text("Energy (E_h)")
                    Text("ΔE (E_h)")
                    Text("Max force")
                    Text("RMS force")
                    Text("Max step")
                    Text("RMS step")
                    if ts {
                        Text("Neg. eig.")
                    }
                    Text("")
                }
                .font(.subheadline.weight(.semibold))
                .foregroundStyle(.secondary)
                Divider().gridCellUnsizedAxes(.horizontal)
                ForEach(opt.history) { s in
                    GridRow {
                        Button {
                            playing = false
                            selectedStep = s.step
                        } label: {
                            Text("\(s.step)")
                                .frame(minWidth: 22, alignment: .leading)
                        }
                        .buttonStyle(.link)
                        .accessibilityLabel("Optimisation step \(s.step)")
                        Text(Fmt.energy10(s.energy))
                        Text(s.deltaE.map { String(format: "%+.2e", $0) } ?? "–")
                        criterionText(s, 0)
                        criterionText(s, 1)
                        criterionText(s, 2)
                        criterionText(s, 3)
                        if ts {
                            Text(s.ts.map { "\($0.negativeEigenvalues)" } ?? "–")
                                .foregroundStyle(s.ts?.negativeEigenvalues == 1 ? Color.green : Color.orange)
                        }
                        Image(systemName: s.converged ? "checkmark.circle.fill" : "circle")
                            .foregroundStyle(s.converged ? Color.green : Color.secondary)
                    }
                    .fontWeight(shownStep(opt)?.step == s.step ? .semibold : .regular)
                }
            }
            .monospacedDigit()
            .textSelection(.enabled)
            Text(ts
                 ? "Green values meet their threshold; Neg. eig.: negative eigenvalues of the Hessian (one at a transition state). Click a step number to show its molecule above and its geometry, Hessian, followed mode and internal coordinates below."
                 : "Green values meet their threshold. Click a step number to show its molecule above and its geometry, internal coordinates and convergence table below.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }

    private func criterionText(_ s: MacQCOptimizationStep, _ index: Int) -> some View {
        let item = index < s.convergence.count ? s.convergence[index] : nil
        return Text(item?.value.map { String(format: "%.2e", $0) } ?? "–")
            .foregroundStyle(item?.met == true ? Color.green : Color.primary)
    }

    private func charts(_ opt: MacQCOptimization) -> some View {
        let shown = shownStep(opt)?.step ?? 0
        let final = opt.finalEnergy
        let energies = opt.history.map { (step: $0.step, value: $0.energy - final) }
        let forces = opt.history.compactMap { s in s.criterion(0).map { (step: s.step, value: max($0, 1.0e-12)) } }
        return HStack(alignment: .top, spacing: 16) {
            VStack(alignment: .leading, spacing: 6) {
                Text("Energy relative to the last step (mE_h)")
                    .font(.subheadline.weight(.semibold))
                Chart {
                    ForEach(energies, id: \.step) { p in
                        LineMark(x: .value("Step", p.step), y: .value("E - E_last (mE_h)", p.value * 1000))
                        PointMark(x: .value("Step", p.step), y: .value("E - E_last (mE_h)", p.value * 1000))
                    }
                    RuleMark(x: .value("Shown step", shown))
                        .foregroundStyle(Color.accentColor.opacity(0.6))
                        .lineStyle(StrokeStyle(lineWidth: 1, dash: [3, 3]))
                }
                .chartXAxisLabel("Step")
                .frame(height: 180)
            }
            .dataPanel()
            VStack(alignment: .leading, spacing: 6) {
                Text("Maximum force (E_h/bohr, logarithmic)")
                    .font(.subheadline.weight(.semibold))
                Chart {
                    ForEach(forces, id: \.step) { p in
                        LineMark(x: .value("Step", p.step), y: .value("Max force", p.value))
                            .foregroundStyle(Color.orange)
                        PointMark(x: .value("Step", p.step), y: .value("Max force", p.value))
                            .foregroundStyle(Color.orange)
                    }
                    RuleMark(x: .value("Shown step", shown))
                        .foregroundStyle(Color.accentColor.opacity(0.6))
                        .lineStyle(StrokeStyle(lineWidth: 1, dash: [3, 3]))
                    RuleMark(y: .value("Threshold", opt.criteria.maxForce))
                        .foregroundStyle(Color.green)
                        .lineStyle(StrokeStyle(lineWidth: 1, dash: [4, 3]))
                        .annotation(position: .top, alignment: .trailing) {
                            Text("threshold").font(.caption2).foregroundStyle(.green)
                        }
                }
                .chartYScale(domain: .automatic(includesZero: false), type: .log)
                .chartXAxisLabel("Step")
                .frame(height: 180)
            }
            .dataPanel()
        }
    }

    private func stepDetail(_ s: MacQCOptimizationStep, opt: MacQCOptimization) -> some View {
        VStack(alignment: .leading, spacing: 14) {
            HStack(alignment: .top, spacing: 16) {
                VStack(alignment: .leading, spacing: 6) {
                    Text("Geometry of this step (Å)")
                        .font(.subheadline.weight(.semibold))
                    AtomTable(atoms: s.geometry.atoms)
                    Text(Fmt.geometryLine(s.geometry))
                        .font(.caption.monospacedDigit())
                        .foregroundStyle(.secondary)
                }
                .dataPanel()
                VStack(alignment: .leading, spacing: 6) {
                    Text("Convergence")
                        .font(.subheadline.weight(.semibold))
                    Grid(alignment: .trailing, horizontalSpacing: 16, verticalSpacing: 4) {
                        GridRow {
                            Text("Item").gridColumnAlignment(.leading)
                            Text("Value")
                            Text("Threshold")
                            Text("Met")
                        }
                        .foregroundStyle(.secondary)
                        ForEach(s.convergence) { c in
                            GridRow {
                                Text(c.name)
                                Text(c.value.map { String(format: "%.2e", $0) } ?? "–")
                                Text(String(format: "%.2e", c.threshold))
                                Text(c.applicable ? (c.met ? "yes" : "no") : "–")
                                    .foregroundStyle(c.met ? Color.green : Color.primary)
                            }
                        }
                    }
                    .monospacedDigit()
                }
                .dataPanel()
            }
            if let ts = s.ts {
                tsPanel(ts, step: s)
            }
            VStack(alignment: .leading, spacing: 6) {
                Text("Internal coordinates (force = −dE/dq; step = change taken from this geometry)")
                    .font(.subheadline.weight(.semibold))
                Grid(alignment: .trailing, horizontalSpacing: 18, verticalSpacing: 4) {
                    GridRow {
                        Text("Name").gridColumnAlignment(.leading)
                        Text("Definition").gridColumnAlignment(.leading)
                        Text("Value")
                        Text("Force")
                        Text("Step")
                    }
                    .foregroundStyle(.secondary)
                    ForEach(s.coordinates) { q in
                        let angle = !q.definition.hasPrefix("R(")
                        GridRow {
                            Text(q.name)
                            Text(q.definition)
                            Text(angle ? String(format: "%.4f°", q.value * 180 / Double.pi)
                                       : String(format: "%.6f Å", q.value * Units.bohrToAngstrom))
                            Text(String(format: "%+.6f", -q.gradient))
                            Text(q.step.map { angle ? String(format: "%+.4f°", $0 * 180 / Double.pi)
                                                    : String(format: "%+.6f Å", $0 * Units.bohrToAngstrom) } ?? "–")
                        }
                    }
                }
                .monospacedDigit()
                .textSelection(.enabled)
                Text(stepCaption(s))
                    .font(.caption)
                    .foregroundStyle(.secondary)
                ForEach(s.warnings, id: \.self) { w in
                    Label(w, systemImage: "exclamationmark.triangle.fill")
                        .font(.caption)
                        .foregroundStyle(.orange)
                }
            }
            .dataPanel()
        }
    }

    /// The Hessian of a transition-state step: its eigenvalues, the followed mode and the shifts
    /// of the partitioned RFO step (Baker 1986).
    private func tsPanel(_ ts: MacQCTSStep, step s: MacQCOptimizationStep) -> some View {
        VStack(alignment: .leading, spacing: 6) {
            Text("Hessian and partitioned RFO step")
                .font(.subheadline.weight(.semibold))
            InfoRow("Hessian", ts.hessianComputed ? "computed at this step" : "Bofill update")
            InfoRow("Eigenvalues (E_h/bohr², E_h/rad²)",
                    ts.hessianEigenvalues.map { String(format: "%.5f", $0) }.joined(separator: "  "))
            InfoRow("Negative eigenvalues", "\(ts.negativeEigenvalues)")
            if let e = ts.followedEigenvalue {
                InfoRow("Followed mode: eigenvalue", String(format: "%.6f", e))
            }
            if let g = ts.followedGradient {
                InfoRow("Followed mode: gradient component", String(format: "%+.6f", g))
            }
            if let o = ts.followedOverlap {
                InfoRow("Overlap with the previous followed mode", String(format: "%.4f", o))
            }
            if let lp = ts.lambdaP, let ln = ts.lambdaN {
                InfoRow("Shifts λp (maximised), λn (minimised)", String(format: "%.6f, %.6f", lp, ln))
            }
            if let mode = ts.followedMode, !mode.isEmpty {
                InfoRow("Largest components of the followed mode", largestComponents(mode, names: s.coordinates.map { $0.definition }))
            }
            Text("The step goes uphill along the followed mode (the lowest at the first step, then the one most like the previous) and downhill along all others; it is scaled to the trust radius of 0.30 bohr/rad.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .monospacedDigit()
        .textSelection(.enabled)
        .dataPanel()
    }

    /// "D(H3-O1-O2-H4) 0.98, A(O2-O1-H3) −0.12, …" (at most four, by size).
    private func largestComponents(_ mode: [Double], names: [String]) -> String {
        let ranked = mode.enumerated().sorted { abs($0.element) > abs($1.element) }.prefix(4)
        return ranked.map { item in
            let name = item.offset < names.count ? names[item.offset] : "q\(item.offset + 1)"
            return String(format: "%@ %+.3f", name, item.element)
        }.joined(separator: ", ")
    }

    private func stepCaption(_ s: MacQCOptimizationStep) -> String {
        var parts = ["Hessian: \(s.hessianUpdate)"]
        if s.stepNorm > 0 {
            parts.append(String(format: "RFO step length %.4f (trust radius %.3f%@), predicted change %.2e E_h",
                                s.stepNorm, s.trustRadius, s.stepScaled ? ", scaled to it" : "", s.predictedChange))
        }
        if let r = s.ratio {
            parts.append(String(format: "previous step: actual/predicted energy change %.3f", r))
        }
        parts.append(String(format: "%.1f s", s.seconds))
        return parts.joined(separator: "; ")
    }
}

/// Cartesian coordinates of a geometry in Å.
struct AtomTable: View {
    var atoms: [MacQCAtom]

    var body: some View {
        Grid(alignment: .trailing, horizontalSpacing: 16, verticalSpacing: 3) {
            GridRow {
                Text("Atom").gridColumnAlignment(.leading)
                Text("x")
                Text("y")
                Text("z")
            }
            .foregroundStyle(.secondary)
            ForEach(Array(atoms.enumerated()), id: \.offset) { index, a in
                GridRow {
                    Text("\(a.symbol)\(index + 1)")
                    ForEach(0..<3, id: \.self) { k in
                        Text(String(format: "%.6f", (k < a.position.count ? a.position[k] : 0) * Units.bohrToAngstrom))
                    }
                }
            }
        }
        .monospacedDigit()
        .textSelection(.enabled)
    }
}

/// Placeholder of a job tab without a result: refusal, warning or how to get one.
struct JobPlaceholder: View {
    var model: AppModel
    var title: String
    var systemImage: String
    var text: String

    var body: some View {
        VStack(spacing: 14) {
            if let message = model.errorMessage {
                ErrorBanner(message: message)
                    .frame(maxWidth: 560)
            }
            if model.isGeneral, let refusal = model.refusal {
                StatusBanner(kind: .refusal, text: refusal)
                    .frame(maxWidth: 560)
            }
            ContentUnavailableView(title, systemImage: systemImage,
                                   description: Text(model.isGeneral ? text : "HeH⁺ in an STO-nG basis (Szabo–Ostlund) runs single-point energies. Choose a Gaussian basis set (6-31G … aug-cc-pCV5Z) for the other job types."))
        }
        .padding()
    }
}

// MARK: - Frequencies and thermochemistry

struct FrequenciesView: View {
    @Bindable var model: AppModel

    var body: some View {
        VStack(spacing: 0) {
            Picker("Page", selection: $model.vibrationPage) {
                ForEach(VibrationPage.allCases) { page in
                    Text(page.rawValue).tag(page)
                }
            }
            .pickerStyle(.segmented)
            .labelsHidden()
            .frame(maxWidth: 360)
            .padding(.top, 10)
            .padding(.bottom, 4)

            if let out = model.result, let f = out.frequencies {
                switch model.vibrationPage {
                case .frequencies:
                    FrequencyPage(model: model, out: out, f: f)
                case .thermochemistry:
                    if let t = out.thermochemistry {
                        ThermochemistryPage(model: model, out: out, t: t)
                    } else {
                        Spacer()
                    }
                }
            } else {
                JobPlaceholder(model: model, title: model.vibrationPage == .frequencies ? "No frequencies" : "No thermochemistry",
                               systemImage: "waveform.path",
                               text: "Choose a job with frequencies (Frequencies, Optimisation + frequencies, TS optimisation + frequencies, with or without a single-point energy) in the sidebar and press Run.")
                    .frame(maxHeight: .infinity)
            }
        }
    }
}

struct FrequencyPage: View {
    var model: AppModel
    var out: MacQCOutput
    var f: MacQCFrequencies
    @State private var selectedMode = 0
    @State private var animate = true
    @State private var amplitude = 1.0
    @State private var arrows = true

    var body: some View {
        ScrollView {
            VStack(alignment: .leading, spacing: 18) {
                if let message = model.errorMessage {
                    ErrorBanner(message: message)
                }
                VStack(alignment: .leading, spacing: 4) {
                    Text("\(out.moleculeTitle) · harmonic frequencies · \(MethodChoice.title(of: out.input.method))/\(out.input.basis)")
                        .font(.title2.weight(.semibold))
                    Text(geometryNote)
                        .foregroundStyle(.secondary)
                }
                if isTransitionStateJob {
                    transitionStateCheck
                }
                NotesPanel(notes: [], warnings: uniqueWarnings)

                SectionTitle("Normal modes")
                if !f.modes.isEmpty {
                    modeAnimations
                }
                modesTable
                if !f.modes.isEmpty {
                    irSpectrum
                }

                SectionTitle("How the Hessian was computed")
                VStack(alignment: .leading, spacing: 6) {
                    InfoRow("Hessian", f.analytic ? "analytic" : "numerical")
                    Text(f.description)
                        .font(.callout)
                        .textSelection(.enabled)
                    if let reason = f.numericalReason {
                        Text("Numerical because: " + reason)
                            .font(.caption)
                            .foregroundStyle(.secondary)
                    }
                    if f.gradientEvaluations > 0 {
                        InfoRow("Gradient evaluations", "\(f.gradientEvaluations)")
                    }
                    if f.energyEvaluations > 0 {
                        InfoRow("Energy evaluations", "\(f.energyEvaluations)")
                    }
                    if !f.coordinates.isEmpty {
                        InfoRow("Displaced coordinates", f.coordinates.joined(separator: ", "))
                    }
                    InfoRow("IR intensities", f.intensityDescription.isEmpty ? (f.intensities ? "yes" : "none") : f.intensityDescription)
                    InfoRow("Isotopes", zip(f.isotopes, f.masses).map { String(format: "%@ %.9f u", $0.0, $0.1) }.joined(separator: ", "))
                    InfoRow("Residual curvature (cm⁻¹)", f.lowFrequencies.map { String(format: "%.2f", $0) }.joined(separator: " "))
                }
                .dataPanel()
            }
            .padding(20)
        }
    }

    /// The warnings of the job and of the thermochemistry, each once (a left-out imaginary mode
    /// is reported by both).
    private var uniqueWarnings: [String] {
        var seen = Set<String>()
        return ((out.job?.warnings ?? []) + (out.thermochemistry?.warnings ?? [])).filter { seen.insert($0).inserted }
    }

    private var isTransitionStateJob: Bool {
        JobChoice(rawValue: out.input.job)?.hasTransitionState ?? false
    }

    /// Where the frequencies were computed, and the point group.
    private var geometryNote: String {
        var text: String
        if out.input.job == "freq" {
            text = "At the input geometry."
        } else if isTransitionStateJob {
            text = "At the transition state."
        } else {
            text = "At the optimised geometry."
        }
        if let g = f.pointGroup {
            text += " Point group \(Fmt.pointGroup(g))."
        }
        return text
    }

    /// Exactly one imaginary frequency confirms the transition state (Gaussian's check after
    /// Opt=TS Freq); its normal mode is the reaction coordinate.
    private var transitionStateCheck: some View {
        let imaginary = f.modes.enumerated().filter { $0.element.imaginary }
        let confirmed = imaginary.count == 1
        let text: String
        if confirmed, let m = imaginary.first {
            text = "One imaginary frequency, \(Fmt.wavenumber(m.element)) cm⁻¹" + (m.element.symmetry.isEmpty ? "" : " (\(m.element.symmetry))")
                + ": mode \(m.offset + 1) is the motion over the barrier."
        } else {
            text = "\(imaginary.count) imaginary frequencies: this is not a transition state (one is needed)."
        }
        return Label(text, systemImage: confirmed ? "checkmark.seal.fill" : "exclamationmark.triangle.fill")
            .foregroundStyle(confirmed ? Color.green : Color.orange)
            .font(.callout)
            .onAppear {
                if confirmed, let m = imaginary.first {
                    selectedMode = m.offset
                }
            }
    }

    /// ΨLab 12.0: the selected normal mode animated in the 3D viewer (turn it with the mouse), with
    /// Play/Pause, the amplitude and the displacement arrows, and a chooser of all modes; clicking
    /// a mode here, in the table or in the spectrum shows it.
    private var modeAnimations: some View {
        let index = min(max(selectedMode, 0), f.modes.count - 1)
        let m = f.modes[index]
        let atoms = f.geometry.atoms
        let equilibrium = atomCoordinates(atoms)
        let shortest = moleculeBonds(atoms).map { vectorDistance(equilibrium[$0.0], equilibrium[$0.1]) }.min() ?? 1.8
        let largest = m.displacement.map { vectorNorm($0) }.max() ?? 1
        // at amplitude 1 the largest atomic displacement is a quarter of the shortest bond
        let a = amplitude * 0.25 * shortest / max(largest, 1e-6)
        // 2 s per vibration at 1600 cm⁻¹, proportionally faster for higher wavenumbers
        let period = min(2.4, max(0.6, 2.0 * 1600.0 / max(abs(m.wavenumber), 1.0)))
        let vibration = Molecule3DVibration(displacement: m.displacement, amplitude: a, period: period,
                                            playing: animate, arrows: arrows)
        return VStack(alignment: .leading, spacing: 10) {
            HStack(spacing: 14) {
                Button {
                    animate.toggle()
                } label: {
                    Label(animate ? "Pause" : "Play", systemImage: animate ? "pause.fill" : "play.fill")
                }
                .help(animate ? "Stop the vibration" : "Animate the vibration")
                Toggle("Arrows", isOn: $arrows)
                    .toggleStyle(.checkbox)
                    .help("Arrows from the atoms towards the turning points")
                Text("Amplitude")
                    .foregroundStyle(.secondary)
                Slider(value: $amplitude, in: 0.4...2.0)
                    .frame(width: 140)
                    .accessibilityLabel("Amplitude")
                Spacer()
            }
            HStack(alignment: .top, spacing: 14) {
                VStack(alignment: .leading, spacing: 6) {
                    HStack(alignment: .firstTextBaseline, spacing: 8) {
                        Text("Mode \(index + 1)")
                            .font(.headline)
                        if !m.symmetry.isEmpty {
                            Text(m.symmetry)
                                .foregroundStyle(.secondary)
                        }
                        Spacer(minLength: 4)
                        Text(Fmt.wavenumber(m) + " cm⁻¹")
                            .font(.headline.monospacedDigit())
                            .foregroundStyle(m.imaginary ? Color.red : Color.primary)
                    }
                    Molecule3DView(atoms: atoms,
                                   reference: atoms,
                                   framing: [equilibrium,
                                             displacedCoordinates(equilibrium, m.displacement, a),
                                             displacedCoordinates(equilibrium, m.displacement, -a)],
                                   vibration: vibration,
                                   accessibilityText: "Animation of normal mode \(index + 1)")
                        .frame(minWidth: 300, maxWidth: .infinity, minHeight: 320, maxHeight: 320)
                    Text(m.irIntensity.map { String(format: "IR intensity %.2f km/mol", $0) } ?? "IR intensity –")
                        .font(.caption)
                        .monospacedDigit()
                        .foregroundStyle(.secondary)
                }
                modeChooser
                    .frame(width: 236)
            }
            Text("The atoms vibrate about the equilibrium geometry (translucent) along the selected normal mode, at a speed proportional to its wavenumber; the amplitude is exaggerated (at 1, the largest atomic displacement is a quarter of the shortest bond). Drag to turn the molecule, pinch or ⌥-scroll to zoom, click atoms to measure them at the equilibrium geometry. Choose a mode on the right, in the table or in the IR spectrum.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }

    /// All modes, two per row: number, symmetry, wavenumber and IR intensity.
    private var modeChooser: some View {
        ScrollView {
            LazyVGrid(columns: [GridItem(.flexible(), spacing: 6), GridItem(.flexible(), spacing: 6)], spacing: 6) {
                ForEach(f.modes.indices, id: \.self) { index in
                    modeTile(index)
                }
            }
            .padding(2)
        }
        .frame(maxHeight: 350)
    }

    private func modeTile(_ index: Int) -> some View {
        let m = f.modes[index]
        let selected = index == selectedMode
        return VStack(alignment: .leading, spacing: 2) {
            HStack(spacing: 4) {
                Text("Mode \(index + 1)")
                    .font(.caption.weight(.semibold))
                if !m.symmetry.isEmpty {
                    Text(m.symmetry)
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
            }
            Text(Fmt.wavenumber(m) + " cm⁻¹")
                .font(.callout.monospacedDigit())
                .foregroundStyle(m.imaginary ? Color.red : Color.primary)
            Text(m.irIntensity.map { String(format: "%.1f km/mol", $0) } ?? "–")
                .font(.caption2.monospacedDigit())
                .foregroundStyle(.secondary)
        }
        .padding(7)
        .frame(maxWidth: .infinity, alignment: .leading)
        .background(RoundedRectangle(cornerRadius: 8).fill(selected ? Color.accentColor.opacity(0.12) : Color.clear))
        .overlay(RoundedRectangle(cornerRadius: 8)
            .strokeBorder(selected ? Color.accentColor : Color.secondary.opacity(0.25), lineWidth: selected ? 2 : 1))
        .contentShape(Rectangle())
        .onTapGesture {
            selectedMode = index
        }
        .accessibilityElement(children: .combine)
        .accessibilityAddTraits(.isButton)
        .accessibilityLabel("Normal mode \(index + 1), \(Fmt.wavenumber(m)) per centimetre")
    }

    private var modesTable: some View {
        Grid(alignment: .trailing, horizontalSpacing: 22, verticalSpacing: 4) {
            GridRow {
                Text("Mode").gridColumnAlignment(.leading)
                Text("Symmetry").gridColumnAlignment(.leading)
                Text("Wavenumber (cm⁻¹)")
                Text("Reduced mass (u)")
                Text("Force constant (mdyn/Å)")
                Text("IR intensity (km/mol)")
            }
            .font(.subheadline.weight(.semibold))
            .foregroundStyle(.secondary)
            Divider().gridCellUnsizedAxes(.horizontal)
            ForEach(Array(f.modes.enumerated()), id: \.offset) { index, m in
                GridRow {
                    Button {
                        selectedMode = index
                    } label: {
                        Text("\(index + 1)").frame(minWidth: 22, alignment: .leading)
                    }
                    .buttonStyle(.link)
                    .accessibilityLabel("Normal mode \(index + 1)")
                    Text(m.symmetry.isEmpty ? "–" : m.symmetry)
                    Text(Fmt.wavenumber(m))
                        .foregroundStyle(m.imaginary ? Color.red : Color.primary)
                    Text(String(format: "%.6f", m.reducedMass))
                    Text(String(format: "%.6f", m.forceConstant))
                    Text(m.irIntensity.map { String(format: "%.4f", $0) } ?? "–")
                }
                .fontWeight(index == selectedMode ? .semibold : .regular)
            }
        }
        .monospacedDigit()
        .textSelection(.enabled)
        .dataPanel()
    }

    private var irSpectrum: some View {
        let sticks = f.modes.enumerated().compactMap { item -> (index: Int, w: Double, i: Double)? in
            guard let i = item.element.irIntensity, !item.element.imaginary else { return nil }
            return (item.offset, item.element.wavenumber, i)
        }
        return VStack(alignment: .leading, spacing: 8) {
            Text("IR stick spectrum")
                .font(.subheadline.weight(.semibold))
            if sticks.isEmpty {
                Text(f.intensities ? "No IR-active real modes." : "No IR intensities for this calculation.")
                    .foregroundStyle(.secondary)
                    .frame(maxWidth: .infinity, minHeight: 220)
            } else {
                Chart {
                    ForEach(sticks, id: \.index) { s in
                        RuleMark(x: .value("Wavenumber (cm⁻¹)", s.w), yStart: .value("Intensity", 0), yEnd: .value("Intensity", s.i))
                            .lineStyle(StrokeStyle(lineWidth: 3))
                            .foregroundStyle(s.index == selectedMode ? Color.orange : Color.blue)
                    }
                }
                .chartXScale(domain: 0...(max(4500, (sticks.map { $0.w }.max() ?? 4000) + 300)))
                .chartXAxisLabel("Wavenumber (cm⁻¹)")
                .chartYAxisLabel("km/mol")
                .frame(height: 220)
            }
        }
        .dataPanel()
    }
}

struct ThermochemistryPage: View {
    var model: AppModel
    var out: MacQCOutput
    var t: MacQCThermochemistry

    private func hartree(_ jPerMol: Double) -> Double { jPerMol / Units.joulePerMolPerHartree }

    var body: some View {
        ScrollView {
            VStack(alignment: .leading, spacing: 18) {
                if let message = model.errorMessage {
                    ErrorBanner(message: message)
                }
                VStack(alignment: .leading, spacing: 4) {
                    Text("\(out.moleculeTitle) · thermochemistry")
                        .font(.title2.weight(.semibold))
                    Text(String(format: "Ideal gas, rigid rotor, harmonic oscillators · T = %.2f K · p = %.2f Pa (%.6f atm, %.5f bar) · symmetry number %ld",
                                t.temperature, t.pressure, t.pressure / 101325, t.pressure / 1.0e5, t.symmetryNumber))
                        .foregroundStyle(.secondary)
                }
                NotesPanel(notes: [], warnings: t.warnings)

                SectionTitle("Energies")
                energyTable
                if let sp = out.singlePoint, let label = out.job?.finalEnergyLabel {
                    Text("The electronic energy is the single-point energy \(label) (\(Fmt.energy10(sp.energy)) E_h); the thermal corrections are from the frequency step (\(MethodChoice.title(of: out.input.method))/\(out.input.basis)).")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }

                SectionTitle("Entropy, heat capacity and thermal energy by contribution")
                contributionTable

                SectionTitle("Rotation and vibration")
                VStack(spacing: 8) {
                    InfoRow("Molecular mass", String(format: "%.6f u", t.mass))
                    InfoRow("Rotor", rotorDescription)
                    InfoRow("Moments of inertia (u bohr²)", t.moments.map { String(format: "%.6f", $0) }.joined(separator: "  "))
                    InfoRow("Rotational constants (GHz)", t.rotationalConstantsGHz.map { String(format: "%.6f", $0) }.joined(separator: "  "))
                    InfoRow("Rotational constants (cm⁻¹)", t.rotationalConstantsCm.map { String(format: "%.6f", $0) }.joined(separator: "  "))
                    InfoRow("Rotational temperatures (K)", t.rotationalTemperatures.map { String(format: "%.5f", $0) }.joined(separator: "  "))
                    InfoRow("Vibrational temperatures (K)", t.vibrationalTemperatures.map { String(format: "%.2f", $0) }.joined(separator: "  "))
                    if t.excludedModes > 0 {
                        InfoRow("Modes left out", "\(t.excludedModes) (imaginary or below 10 cm⁻¹)")
                    }
                }
                .dataPanel()
                Text("Constants: CODATA 2018; masses of the most abundant isotopes (AME 2020); electronic ground state a closed-shell singlet (q_el = 1).")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
            .padding(20)
        }
    }

    /// Linear, spherical, symmetric (NH₃) or asymmetric top, from the moments of inertia.
    private var rotorDescription: String {
        if t.linear { return "linear" }
        let m = t.moments.sorted()
        guard m.count == 3, let largest = m.last, largest > 0 else { return "nonlinear" }
        let equal = { (a: Double, b: Double) in abs(a - b) <= 1e-4 * largest }
        if equal(m[0], m[1]) && equal(m[1], m[2]) { return "nonlinear (spherical top)" }
        if equal(m[0], m[1]) { return "nonlinear (oblate symmetric top)" }
        if equal(m[1], m[2]) { return "nonlinear (prolate symmetric top)" }
        return "nonlinear (asymmetric top)"
    }

    private var energyTable: some View {
        let rows: [(label: String, value: Double)] = [
            ("Zero-point energy", t.corrections.zeroPoint),
            ("Thermal correction to U", t.corrections.thermalEnergy),
            ("Thermal correction to H = U + pV", t.corrections.enthalpy),
            ("Thermal correction to A = U − TS", t.corrections.helmholtz),
            ("Thermal correction to G = H − TS", t.corrections.gibbs),
        ]
        let totals: [(label: String, value: Double)] = [
            ("E (electronic)", t.electronicEnergy),
            ("E + ZPE", t.totals.zeroPoint),
            ("U", t.totals.thermalEnergy),
            ("H", t.totals.enthalpy),
            ("A", t.totals.helmholtz),
            ("G", t.totals.gibbs),
        ]
        return VStack(alignment: .leading, spacing: 6) {
        HStack(alignment: .top, spacing: 16) {
            Grid(alignment: .trailing, horizontalSpacing: 18, verticalSpacing: 4) {
                GridRow {
                    Text("Correction").gridColumnAlignment(.leading)
                    Text("E_h")
                    Text("kcal/mol")
                    Text("kJ/mol")
                }
                .foregroundStyle(.secondary)
                Divider().gridCellUnsizedAxes(.horizontal)
                ForEach(rows, id: \.label) { row in
                    GridRow {
                        Text(row.label)
                        Text(String(format: "%.8f", hartree(row.value)))
                        Text(String(format: "%.4f", row.value / 1000 / Units.joulePerCalorie))
                        Text(String(format: "%.4f", row.value / 1000))
                    }
                }
            }
            .monospacedDigit()
            .textSelection(.enabled)
            .dataPanel()
            Grid(alignment: .trailing, horizontalSpacing: 18, verticalSpacing: 4) {
                GridRow {
                    Text("Total").gridColumnAlignment(.leading)
                    Text("E_h")
                }
                .foregroundStyle(.secondary)
                Divider().gridCellUnsizedAxes(.horizontal)
                ForEach(totals, id: \.label) { row in
                    GridRow {
                        Text(row.label)
                        Text(Fmt.energy10(row.value))
                    }
                }
            }
            .monospacedDigit()
            .textSelection(.enabled)
            .dataPanel()
        }
        Text("H = U + pV; for the ideal gas pV = RT per mole (k_BT per molecule).")
            .font(.caption)
            .foregroundStyle(.secondary)
        }
    }

    private var contributionTable: some View {
        let terms: [(label: String, term: MacQCThermoTerm)] = [
            ("Electronic", t.terms.electronic),
            ("Translational", t.terms.translation),
            ("Rotational", t.terms.rotation),
            ("Vibrational", t.terms.vibration),
            ("Total", t.terms.total),
        ]
        return VStack(alignment: .leading, spacing: 6) {
            Grid(alignment: .trailing, horizontalSpacing: 18, verticalSpacing: 4) {
                GridRow {
                    Text("").gridColumnAlignment(.leading)
                    Text("E (kcal/mol)")
                    Text("C_V (cal/(mol K))")
                    Text("S (cal/(mol K))")
                    Text("S (J/(mol K))")
                    Text("q")
                }
                .foregroundStyle(.secondary)
                Divider().gridCellUnsizedAxes(.horizontal)
                ForEach(terms, id: \.label) { row in
                    GridRow {
                        Text(row.label)
                        Text(String(format: "%.4f", row.term.energy / 1000 / Units.joulePerCalorie))
                        Text(String(format: "%.4f", row.term.cv / Units.joulePerCalorie))
                        Text(String(format: "%.4f", row.term.entropy / Units.joulePerCalorie))
                        Text(String(format: "%.4f", row.term.entropy))
                        Text(String(format: "%.6e", row.term.q))
                    }
                    .fontWeight(row.label == "Total" ? .semibold : .regular)
                }
            }
            .monospacedDigit()
            .textSelection(.enabled)
            Text("The vibrational energy includes the zero-point energy; the vibrational partition function is counted from the bottom of the well.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }
}

// MARK: - Summary panel of a job

/// Job, final energy, optimised geometry, wavenumbers and thermochemistry (Summary tab).
struct JobSummaryPanel: View {
    var model: AppModel
    var out: MacQCOutput

    var body: some View {
        if let j = out.job, j.type != "energy" {
            VStack(alignment: .leading, spacing: 8) {
                InfoRow("Job", j.title)
                InfoRow("Termination", j.termination == "normal" ? "normal" : (j.message ?? j.termination))
                if let label = j.finalEnergyLabel, let e = j.finalEnergy {
                    InfoRow("Final energy (\(label))", Fmt.energy10(e) + " E_h")
                }
                if let opt = out.optimization {
                    InfoRow(opt.isTransitionState ? "Transition-state search" : "Optimisation",
                            opt.converged ? "converged in \(opt.steps) steps" : "not converged (\(opt.steps) steps)")
                    if opt.isTransitionState {
                        InfoRow("Hessian", Fmt.tsHessian(opt))
                    }
                    ForEach(opt.finalParameters) { p in
                        InfoRow("  " + p.name, Fmt.parameterValue(p.value, angle: p.angle))
                    }
                }
                if let scan = out.scan {
                    InfoRow("Scanned coordinate", scan.coordinate.definition + " (\(scan.coordinate.type))")
                    InfoRow("Points", "\(scan.points.count) of \(scan.steps + 1) · " + scan.message)
                    if let low = scan.points.first(where: { $0.index == scan.lowest }) {
                        InfoRow("Lowest point", Fmt.scanValue(low.value, angular: scan.isAngular) + " · " + Fmt.energy10(low.energy) + " E_h")
                        if let high = scan.points.first(where: { $0.index == scan.highest }), high.index != low.index {
                            InfoRow("Highest point", Fmt.scanValue(high.value, angular: scan.isAngular)
                                    + String(format: " · %.4f kcal/mol above the lowest", (high.energy - low.energy) * Units.kcalPerMolPerHartree))
                        }
                    }
                    ForEach(scan.extrema) { x in
                        InfoRow("  " + x.kind.capitalized + " (parabola)",
                                Fmt.scanValue(x.value, angular: scan.isAngular) + " · " + Fmt.energy10(x.energy) + " E_h")
                    }
                }
                if let f = out.frequencies {
                    InfoRow("Harmonic wavenumbers (cm⁻¹)", f.modes.map { Fmt.wavenumber($0) }.joined(separator: "  "))
                    InfoRow("Hessian", f.analytic ? "analytic" : "numerical (\(f.gradientEvaluations > 0 ? "\(f.gradientEvaluations) gradients" : "\(f.energyEvaluations) energies"))")
                    if let g = f.pointGroup {
                        InfoRow("Point group", Fmt.pointGroup(g))
                    }
                    if JobChoice(rawValue: out.input.job)?.hasTransitionState == true {
                        let n = f.modes.filter { $0.imaginary }.count
                        InfoRow("Transition state", n == 1 ? "confirmed: one imaginary frequency" : "not confirmed: \(n) imaginary frequencies")
                    }
                }
                if let t = out.thermochemistry {
                    InfoRow(String(format: "ZPE, H, G at %.2f K", t.temperature),
                            String(format: "%.8f, %@, %@ E_h", t.corrections.zeroPoint / Units.joulePerMolPerHartree,
                                   Fmt.energy10(t.totals.enthalpy), Fmt.energy10(t.totals.gibbs)))
                    InfoRow("Entropy", String(format: "%.4f J/(mol K)", t.terms.total.entropy))
                }
                if model.optimizedGeometry != nil {
                    Button {
                        model.useOptimizedGeometry()
                    } label: {
                        Label(out.optimization?.isTransitionState == true ? "Use transition-state geometry as input" : "Use optimised geometry as input",
                              systemImage: "arrow.uturn.backward")
                    }
                    .secondaryActionStyle()
                    .disabled(!model.canUseOptimizedGeometry || model.isRunning)
                }
                SaveOptimizedGeometryBox(model: model)
            }
            .dataPanel()
        }
    }
}
