//
//  NuclearDetailView.swift
//  MacQC
//
//  ΨLab 20.0: the result tabs of the nuclear-motion calculation: Levels (the energy levels, a
//  level diagram), Expectation values (the bond lengths and the angle of every level, the K
//  weights for J > 0), Convergence (the eigensolver's history, the DVR grids) and Output (the
//  output file, streamed while the job runs).
//

import SwiftUI
import Charts

struct NuclearDetailView: View {
    @Bindable var model: NuclearModel

    var body: some View {
        TabView(selection: $model.selectedTab) {
            NuclearLevelsView(model: model)
                .tabItem { Label("Levels", systemImage: "list.number") }
                .tag(NuclearTab.levels)
            NuclearExpectationView(model: model)
                .tabItem { Label("Geometry", systemImage: "ruler") }
                .tag(NuclearTab.expectation)
            NuclearMatrixView(model: model)
                .tabItem { Label("Matrix", systemImage: "square.grid.3x3") }
                .tag(NuclearTab.matrix)
            NuclearConvergenceView(model: model)
                .tabItem { Label("Convergence", systemImage: "chart.line.downtrend.xyaxis") }
                .tag(NuclearTab.convergence)
            NuclearOutputView(model: model)
                .tabItem { Label("Output", systemImage: "doc.plaintext") }
                .tag(NuclearTab.output)
        }
        .padding(.top, 4)
    }
}

/// Placeholder of a tab without a result.
private struct NuclearPlaceholder: View {
    var model: NuclearModel

    var body: some View {
        VStack(spacing: 14) {
            if let message = model.errorMessage {
                ErrorBanner(message: message)
                    .frame(maxWidth: 560)
            }
            if model.isRunning {
                ProgressView(model.progressLine.isEmpty ? "Running…" : model.progressLine)
                    .frame(maxWidth: 560)
            } else {
                ContentUnavailableView(
                    "No energy levels yet",
                    systemImage: "waveform.path.ecg",
                    description: Text("Choose the molecule, the basis (ns1, np) and jmax in the sidebar and press Run. jmax = 0 gives the vibrational levels, jmax > 0 the rovibrational levels of that J.")
                )
            }
        }
        .padding()
        .frame(maxWidth: .infinity, maxHeight: .infinity)
    }
}

/// Header of the result: molecule, J and the basis.
private struct NuclearHeader: View {
    var result: NuclearResult

    var body: some View {
        let input = result.input ?? NuclearInput()
        let iso = Isotopologue(rawValue: result.isotopologue ?? "") ?? input.isotopologue
        VStack(alignment: .leading, spacing: 4) {
            HStack(spacing: 0) {
                LaTeXText(iso.latex, size: 20, weight: .semibold)
                Text(input.J == 0 ? " · vibrational levels (J = 0)" : " · rovibrational levels, J = \(input.J), K = \(input.parityBlock == 0 ? 0 : 1) … \(input.J)")
                    .font(.title2.weight(.semibold))
            }
            Text("CVRQD surface · ns1 = \(input.radialPoints), np = \(input.angularPoints) · \(Fmt.count(Double(result.dimension ?? 0))) basis functions · \(NuclearEigensolver(rawValue: input.eigensolver.rawValue)?.title ?? "")")
                .foregroundStyle(.secondary)
        }
    }
}

// MARK: - Levels

struct NuclearLevelsView: View {
    var model: NuclearModel

    var body: some View {
        if let r = model.result, !r.levels.isEmpty {
            ScrollView {
                VStack(alignment: .leading, spacing: 18) {
                    if let message = model.errorMessage {
                        ErrorBanner(message: message)
                    }
                    NuclearHeader(result: r)
                    cards(r)
                    if r.converged == false {
                        Label("Not converged: the levels below are the best approximations after \(r.iterations ?? 0) iterations (see Convergence).",
                              systemImage: "exclamationmark.triangle.fill")
                            .foregroundStyle(.orange)
                    }
                    SectionTitle("Energy levels")
                    // the table at its own width (scrolls sideways in a narrow window), the diagram below
                    ScrollView(.horizontal) {
                        levelTable(r)
                            .fixedSize()
                            .dataPanel()
                    }
                    levelDiagram(r)
                        .frame(height: 320)
                        .dataPanel()
                    Text("Energies from the minimum of the potential energy surface; E − E₁ from the lowest level of this calculation (for jmax = 0: the vibrational band origins above the zero-point level). ‖Hx − Ex‖ is the residual of the eigenvector (E_h).")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
                .padding(20)
            }
        } else {
            NuclearPlaceholder(model: model)
        }
    }

    private func cards(_ r: NuclearResult) -> some View {
        let e1 = r.levels.first?.energy ?? 0
        let input = r.input ?? NuclearInput()
        return GlassGroup(spacing: 14) {
            HStack(spacing: 14) {
                EnergyCard(title: input.J == 0 ? "Zero-point energy" : "Lowest level", value: e1 * NuclearUnits.hartreeToCm, decimals: 4, unit: "cm⁻¹")
                if r.levels.count > 1 {
                    EnergyCard(title: "Second level − first", value: (r.levels[1].energy - e1) * NuclearUnits.hartreeToCm, decimals: 4, unit: "cm⁻¹")
                }
                EnergyCard(title: "Potential minimum (E_h)", value: r.potentialMinimum ?? 0, decimals: 10)
            }
        }
    }

    private func levelTable(_ r: NuclearResult) -> some View {
        let e1 = r.levels.first?.energy ?? 0
        let kValues = r.kValues ?? [0]
        let showK = kValues.count > 1
        return Grid(alignment: .trailing, horizontalSpacing: 18, verticalSpacing: 6) {
            GridRow {
                Text("n")
                Text("E (cm⁻¹)")
                Text("E − E₁ (cm⁻¹)")
                Text("E (E_h)")
                Text("‖Hx − Ex‖")
                if showK { Text("main K") }
            }
            .font(.caption.weight(.semibold))
            .foregroundStyle(.secondary)
            Divider()
            ForEach(Array(r.levels.enumerated()), id: \.offset) { k, l in
                GridRow {
                    Text("\(k + 1)")
                    Text(String(format: "%.4f", l.energy * NuclearUnits.hartreeToCm))
                    Text(String(format: "%.4f", (l.energy - e1) * NuclearUnits.hartreeToCm))
                    Text(String(format: "%.10f", l.energy))
                    Text(String(format: "%.1e", l.residual))
                        .foregroundStyle(l.residual < (r.input?.residualTolerance ?? 1e-8) ? Color.primary : Color.orange)
                    if showK, let i = l.kWeights.indices.max(by: { l.kWeights[$0] < l.kWeights[$1] }), i < kValues.count {
                        Text(String(format: "%d (%.0f %%)", kValues[i], 100 * l.kWeights[i]))
                    }
                }
                .monospacedDigit()
                .textSelection(.enabled)
            }
        }
    }

    /// A stick diagram of the levels above the first.
    private func levelDiagram(_ r: NuclearResult) -> some View {
        let e1 = r.levels.first?.energy ?? 0
        let points = r.levels.enumerated().map { (k: $0.offset + 1, e: ($0.element.energy - e1) * NuclearUnits.hartreeToCm) }
        return Chart(points, id: \.k) { p in
            RuleMark(xStart: .value("left", 0.15), xEnd: .value("right", 0.85), y: .value("E − E₁ (cm⁻¹)", p.e))
                .lineStyle(StrokeStyle(lineWidth: 2))
                .foregroundStyle(Color.accentColor)
                .annotation(position: .trailing, alignment: .leading, spacing: 4) {
                    Text("\(p.k)")
                        .font(.caption2.monospacedDigit())
                        .foregroundStyle(.secondary)
                }
        }
        .chartXScale(domain: 0...1)
        .chartYScale(domain: 0...max(1.0, (points.last?.e ?? 1) * 1.05))
        .chartXAxis(.hidden)
        .chartYAxisLabel("E − E₁ (cm⁻¹)")
        .accessibilityLabel("Energy level diagram")
    }
}

// MARK: - Expectation values

struct NuclearExpectationView: View {
    var model: NuclearModel

    var body: some View {
        if let r = model.result, !r.levels.isEmpty {
            ScrollView {
                VStack(alignment: .leading, spacing: 18) {
                    NuclearHeader(result: r)
                    SectionTitle("Vibrationally averaged geometry: expectation values of the bond lengths and the angle")
                    momentsTable(r)
                        .dataPanel()
                    Text("r₁ = r(O–H₁), r₂ = r(O–H₂) in Å; ⟨rⁿ⟩^(1/n) for n = 1, 2, 3, −1, −2, −3 (as the original program's expectval); the angle H–O–H in degrees. The DVR quadrature: |cᵢ|² at the grid points.")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                    if let kv = r.kValues, kv.count > 1 {
                        SectionTitle("Weights of the K blocks")
                        kTable(r, kv)
                            .dataPanel()
                    }
                }
                .padding(20)
            }
        } else {
            NuclearPlaceholder(model: model)
        }
    }

    private func momentsTable(_ r: NuclearResult) -> some View {
        let a = NuclearUnits.bohrToAngstrom
        func f(_ x: Double) -> String { String(format: "%.5f", x * a) }
        return Grid(alignment: .trailing, horizontalSpacing: 14, verticalSpacing: 6) {
            GridRow {
                Text("n")
                Text("bond")
                Text("⟨r⟩")
                Text("⟨r²⟩^½")
                Text("⟨r³⟩^⅓")
                Text("⟨r⁻¹⟩⁻¹")
                Text("⟨r⁻²⟩^−½")
                Text("⟨r⁻³⟩^−⅓")
                Text("⟨HOH⟩ (°)")
            }
            .font(.caption.weight(.semibold))
            .foregroundStyle(.secondary)
            Divider()
            ForEach(Array(r.levels.enumerated()), id: \.offset) { k, l in
                GridRow {
                    Text("\(k + 1)")
                    Text("r₁")
                    Text(f(l.r1.mean))
                    Text(f(l.r1.rms))
                    Text(f(l.r1.cubic))
                    Text(f(l.r1.inverse))
                    Text(f(l.r1.inverse2))
                    Text(f(l.r1.inverse3))
                    Text(String(format: "%.3f", l.theta))
                }
                .monospacedDigit()
                GridRow {
                    Text("")
                    Text("r₂")
                    Text(f(l.r2.mean))
                    Text(f(l.r2.rms))
                    Text(f(l.r2.cubic))
                    Text(f(l.r2.inverse))
                    Text(f(l.r2.inverse2))
                    Text(f(l.r2.inverse3))
                    Text("")
                }
                .monospacedDigit()
                .foregroundStyle(.secondary)
            }
        }
        .textSelection(.enabled)
    }

    private func kTable(_ r: NuclearResult, _ kv: [Int]) -> some View {
        Grid(alignment: .trailing, horizontalSpacing: 14, verticalSpacing: 6) {
            GridRow {
                Text("n")
                ForEach(kv, id: \.self) { k in Text("K = \(k)") }
            }
            .font(.caption.weight(.semibold))
            .foregroundStyle(.secondary)
            Divider()
            ForEach(Array(r.levels.enumerated()), id: \.offset) { n, l in
                GridRow {
                    Text("\(n + 1)")
                    ForEach(Array(l.kWeights.enumerated()), id: \.offset) { _, w in
                        Text(String(format: "%.5f", w))
                    }
                }
                .monospacedDigit()
            }
        }
        .textSelection(.enabled)
    }
}

// MARK: - Convergence

struct NuclearConvergenceView: View {
    var model: NuclearModel

    var body: some View {
        if let r = model.result, !r.levels.isEmpty {
            ScrollView {
                VStack(alignment: .leading, spacing: 18) {
                    NuclearHeader(result: r)
                    SectionTitle("Eigensolver")
                    stats(r)
                        .dataPanel()
                    if let h = r.history, h.count > 1 {
                        SectionTitle("History")
                        HStack(spacing: 18) {
                            lowestChart(h)
                                .frame(minHeight: 240)
                                .dataPanel()
                            changeChart(h, lanczos: r.input?.eigensolver == .lanczos)
                                .frame(minHeight: 240)
                                .dataPanel()
                        }
                    }
                    SectionTitle("Grids")
                    grids(r)
                        .dataPanel()
                }
                .padding(20)
            }
        } else {
            NuclearPlaceholder(model: model)
        }
    }

    private func stats(_ r: NuclearResult) -> some View {
        VStack(spacing: 8) {
            InfoRow("Method", r.solver ?? "")
            InfoRow("Converged", (r.converged ?? false) ? "yes" : "no")
            InfoRow("Iterations", "\(r.iterations ?? 0)")
            InfoRow("Products H x", Fmt.count(Double(r.matrixVectorProducts ?? 0)))
            if let li = r.linearIterations, li > 0 {
                InfoRow("Linear-solver iterations", Fmt.count(Double(li)))
            }
            InfoRow("Hamiltonian", String(format: "%@ × %@, %@ non-zero elements, %.2f s", Fmt.count(Double(r.dimension ?? 0)),
                                          Fmt.count(Double(r.dimension ?? 0)), Fmt.count(Double(r.nonzeros ?? 0)), r.hamiltonianSeconds ?? 0))
            InfoRow("Eigensolver time", String(format: "%.2f s", r.solverSeconds ?? 0))
            if let m = r.masses, m.count == 3, let mu1 = r.mu1, let mu2 = r.mu2 {
                InfoRow("Masses", String(format: "%.9f, %.9f, %.7f u", m[0], m[1], m[2]))
                InfoRow("Reduced masses", String(format: "μ₁ = %.4f, μ₂ = %.4f mₑ", mu1, mu2))
            }
            if let vmin = r.potentialMinimum, let rm = r.minimumR, let am = r.minimumAngle {
                InfoRow("Potential minimum", String(format: "%.12f E_h at r(OH) = %.5f Å, %.3f°", vmin, rm * NuclearUnits.bohrToAngstrom, am))
            }
            ForEach(r.notes ?? [], id: \.self) { note in
                Text(note)
                    .font(.caption)
                    .foregroundStyle(.orange)
                    .frame(maxWidth: .infinity, alignment: .leading)
            }
        }
    }

    private func lowestChart(_ h: [NuclearStep]) -> some View {
        Chart(h.filter { $0.lowest > 0 }, id: \.iteration) { s in
            LineMark(x: .value("Iteration", s.iteration), y: .value("Lowest level (cm⁻¹)", s.lowest * NuclearUnits.hartreeToCm))
            PointMark(x: .value("Iteration", s.iteration), y: .value("Lowest level (cm⁻¹)", s.lowest * NuclearUnits.hartreeToCm))
                .symbolSize(12)
        }
        .chartYScale(domain: .automatic(includesZero: false))
        .chartXAxisLabel("Iteration")
        .chartYAxisLabel("Lowest level (cm⁻¹)")
    }

    private func changeChart(_ h: [NuclearStep], lanczos: Bool) -> some View {
        let pts = h.compactMap { s -> (Int, Double)? in
            guard let c = s.change, c > 0 else { return nil }
            return (s.iteration, log10(c))
        }
        return Chart(pts, id: \.0) { p in
            LineMark(x: .value("Iteration", p.0), y: .value("log10", p.1))
            PointMark(x: .value("Iteration", p.0), y: .value("log10", p.1))
                .symbolSize(12)
        }
        .chartXAxisLabel("Iteration")
        .chartYAxisLabel(lanczos ? "log₁₀ change of the levels (E_h)" : "log₁₀ residual (E_h)")
    }

    private func grids(_ r: NuclearResult) -> some View {
        VStack(alignment: .leading, spacing: 10) {
            if let g = r.rGrid {
                Text("Radial points (R₁ and R₂), Å")
                    .font(.caption.weight(.semibold))
                    .foregroundStyle(.secondary)
                Text(g.map { String(format: "%.4f", $0 * NuclearUnits.bohrToAngstrom) }.joined(separator: "  "))
                    .font(.callout.monospacedDigit())
                    .textSelection(.enabled)
                    .fixedSize(horizontal: false, vertical: true)
            }
            if let c = r.cosGrid, let kv = r.kValues {
                ForEach(Array(zip(kv, c)), id: \.0) { k, q in
                    Text("Angular points, K = \(k): Jacobi angle θ (°)")
                        .font(.caption.weight(.semibold))
                        .foregroundStyle(.secondary)
                    Text(q.map { String(format: "%.2f", acos(max(-1, min(1, $0))) * 180 / .pi) }.joined(separator: "  "))
                        .font(.callout.monospacedDigit())
                        .textSelection(.enabled)
                        .fixedSize(horizontal: false, vertical: true)
                }
            }
        }
        .frame(maxWidth: .infinity, alignment: .leading)
    }
}

// MARK: - Output

struct NuclearOutputView: View {
    @Bindable var model: NuclearModel
    @State private var controller = OutputTextController()
    @State private var follow = true

    var body: some View {
        let text = model.outputText.isEmpty ? (model.result?.output ?? "") : model.outputText
        VStack(spacing: 0) {
            HStack(spacing: 10) {
                Text("Output file of the nuclear-motion calculation")
                    .foregroundStyle(.secondary)
                Spacer()
                if model.isRunning {
                    ProgressView()
                        .controlSize(.small)
                    Text("Running")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
                Toggle("Follow output", isOn: $follow)
                    .toggleStyle(.checkbox)
            }
            .padding(.horizontal, 14)
            .padding(.vertical, 8)
            Divider()
            if text.isEmpty {
                NuclearPlaceholder(model: model)
            } else {
                OutputTextView(text: text, generation: model.outputGeneration, follow: follow, controller: controller)
            }
        }
    }
}
