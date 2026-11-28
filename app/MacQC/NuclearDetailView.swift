//
//  NuclearDetailView.swift
//  MacQC
//
//  ΨLab 20.0/21.0: the result tabs of the nuclear-motion calculation: Levels (the energy levels, a
//  level diagram), Geometry (the expectation values of the bond lengths and the angle of every
//  level, the K weights for J > 0), Matrix, Eigenvectors (21.0, NuclearEigenvectorView),
//  Convergence (the eigensolver's history, 21.0 with the linear solves, the DVR grids) and Output
//  (the output file, streamed while the job runs).
//

import SwiftUI
import Charts

struct NuclearDetailView: View {
    @Bindable var model: NuclearModel

    var body: some View {
        TabView(selection: $model.selectedTab) {
            NuclearLevelsView(model: model)
                .modifier(NuclearOpenedFileInset(model: model))
                .tabItem { Label("Levels", systemImage: "list.number") }
                .tag(NuclearTab.levels)
            NuclearExpectationView(model: model)
                .modifier(NuclearOpenedFileInset(model: model))
                .tabItem { Label("Geometry", systemImage: "ruler") }
                .tag(NuclearTab.expectation)
            NuclearMatrixView(model: model)
                .modifier(NuclearOpenedFileInset(model: model))
                .tabItem { Label("Matrix", systemImage: "square.grid.3x3") }
                .tag(NuclearTab.matrix)
            NuclearEigenvectorView(model: model)
                .modifier(NuclearOpenedFileInset(model: model))
                .tabItem { Label("Eigenvectors", systemImage: "waveform") }
                .tag(NuclearTab.eigenvectors)
            NuclearConvergenceView(model: model)
                .modifier(NuclearOpenedFileInset(model: model))
                .tabItem { Label("Convergence", systemImage: "chart.line.downtrend.xyaxis") }
                .tag(NuclearTab.convergence)
            NuclearOutputView(model: model)
                .modifier(NuclearOpenedFileInset(model: model))
                .tabItem { Label("Output", systemImage: "doc.plaintext") }
                .tag(NuclearTab.output)
        }
        .padding(.top, 4)
        .sheet(isPresented: $model.showExportSheet) {
            NuclearExportSheet(model: model)
        }
    }
}

/// 22.0: the banner of an opened nuclear-motion output file above the content of a tab.
struct NuclearOpenedFileInset: ViewModifier {
    var model: NuclearModel

    func body(content: Content) -> some View {
        content.safeAreaInset(edge: .top, spacing: 0) {
            if let file = model.openedFile, model.result != nil, !file.example {
                OpenedFileBanner(file: file)
            }
        }
    }
}

/// Placeholder of a tab without a result.
struct NuclearPlaceholder: View {
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
            } else if model.loadingExample {
                // 23.0: the startup example (H₂¹⁶O, J = 10) is read from the app's resources
                ProgressView("Reading the example: H₂¹⁶O, J = 10, ns1 = np = 50, 100 levels…")
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
struct NuclearHeader: View {
    var result: NuclearResult

    var body: some View {
        let input = result.input ?? NuclearInput()
        let iso = Isotopologue(rawValue: result.isotopologue ?? "") ?? input.isotopologue
        let linear = input.eigensolver.usesLinearSolver ? " (" + input.linearSolver.title + ", " + input.preconditioner.title + ")" : ""
        VStack(alignment: .leading, spacing: 4) {
            HStack(spacing: 0) {
                LaTeXText(iso.latex, size: 20, weight: .semibold)
                Text(input.J == 0 ? " · vibrational levels (J = 0)" : " · rovibrational levels, J = \(input.J), K = \(input.parityBlock == 0 ? 0 : 1) … \(input.J)")
                    .font(.title2.weight(.semibold))
            }
            Text("CVRQD surface · \(input.coordinates.title) coordinates · ns1 = \(input.radialPoints), np = \(input.angularPoints) · \(Fmt.count(Double(result.dimension ?? 0))) basis functions · \(input.eigensolver.title)\(linear)")
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
                    HStack {
                        SectionTitle("Energy levels")
                        Spacer()
                        // 21.0: the eigenvalues, eigenvectors and/or Gaussian start vector as text
                        Button {
                            model.showExportSheet = true
                        } label: {
                            Label("Export…", systemImage: "square.and.arrow.up")
                        }
                        .disabled(r.exportId == nil)
                        .help("Export the eigenvalues, the eigenvectors and/or the Gaussian start vector as text")
                    }
                    // the table at its own width (scrolls sideways in a narrow window), the diagram below
                    ScrollView(.horizontal) {
                        levelTable(r)
                            .fixedSize()
                            .dataPanel()
                    }
                    levelDiagram(r)
                        .frame(height: 320)
                        .dataPanel()
                    Text("Energies from the minimum of the potential energy surface; E − E₁ from the lowest level of this calculation (for jmax = 0: the vibrational band origins above the zero-point level). ‖Hx − Ex‖ is the residual of the eigenvector (E_h); — when only the eigenvalues were computed by Lánczos.")
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
            // 21.0: the order of the output file and of the export: E (E_h) first
            GridRow {
                Text("n")
                Text("E (E_h)")
                Text("E (cm⁻¹)")
                Text("E − E₁ (cm⁻¹)")
                Text("‖Hx − Ex‖ (E_h)")
                if showK { Text("main K") }
            }
            .font(.caption.weight(.semibold))
            .foregroundStyle(.secondary)
            Divider()
            ForEach(Array(r.levels.enumerated()), id: \.offset) { k, l in
                GridRow {
                    Text("\(k + 1)")
                    Text(String(format: "%.10f", l.energy))
                    Text(String(format: "%.4f", l.energy * NuclearUnits.hartreeToCm))
                    Text(String(format: "%.4f", (l.energy - e1) * NuclearUnits.hartreeToCm))
                    if let res = l.residual {
                        Text(String(format: "%.1e", res))
                            .foregroundStyle(res < (r.input?.residualTolerance ?? 1e-8) ? Color.primary : Color.orange)
                    } else {
                        Text("—")
                            .foregroundStyle(.secondary)
                    }
                    if showK {
                        if let i = l.kWeights.indices.max(by: { l.kWeights[$0] < l.kWeights[$1] }), i < kValues.count {
                            Text(String(format: "%d (%.0f %%)", kValues[i], 100 * l.kWeights[i]))
                        } else {
                            Text("—")
                                .foregroundStyle(.secondary)
                        }
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

    /// 21.0: the table in a larger type (Task: the values were hard to read).
    private static let valueFont = Font.system(size: 15).monospacedDigit()
    private static let headerFont = Font.system(size: 13, weight: .semibold)
    /// 23.0: the headers of the expectation values: the quantity and its unit on one line, with
    /// real superscripts (LaTeX), in the size of the values (22.0: 13 pt on two lines)
    private static let momentHeaderSize: CGFloat = 15

    /// 22.0: a table across the whole width of the window (its columns share the free width);
    /// in a window too narrow for it, at its own width with a horizontal scroll bar.
    private func fullWidth<T: View>(_ table: T) -> some View {
        ViewThatFits(in: .horizontal) {
            table
                .frame(maxWidth: .infinity)
                .dataPanel()
            ScrollView(.horizontal) {
                table
                    .fixedSize()
                    .dataPanel()
            }
        }
    }

    var body: some View {
        if let r = model.result, !r.levels.isEmpty {
            ScrollView {
                VStack(alignment: .leading, spacing: 18) {
                    NuclearHeader(result: r)
                    if r.hasVectors == false {
                        ContentUnavailableView("Eigenvalues only", systemImage: "ruler",
                                               description: Text("The expectation values need the eigenvectors. Choose “Eigenvalues and Eigenvectors” under Levels in the sidebar and run again."))
                    } else if r.input?.expectationValues == false {
                        ContentUnavailableView("No expectation values", systemImage: "ruler",
                                               description: Text("Switch on “Expectation values” under Levels in the sidebar and run again."))
                    } else {
                        SectionTitle("Vibrationally averaged geometry: expectation values of the bond lengths and the angle")
                        fullWidth(momentsTable(r))
                        Text("r₁ = r(O–H₁), r₂ = r(O–H₂) in Å; ⟨rⁿ⟩^(1/n) for n = 1, 2, 3, −1, −2, −3 (as the original program's expectval); the angle H–O–H in degrees. The DVR quadrature: |cᵢ|² at the grid points.")
                            .font(.callout)
                            .foregroundStyle(.secondary)
                    }
                    if r.hasVectors != false, let kv = r.kValues, kv.count > 1 {
                        SectionTitle("Weights of the K blocks")
                        fullWidth(kTable(r, kv))
                    }
                }
                .padding(20)
            }
        } else {
            NuclearPlaceholder(model: model)
        }
    }

    /// 22.0: a cell of a value column: as wide as the column, which takes its share of the free
    /// width of a full-width table (fullWidth); the numbers right-aligned.
    private func cell(_ text: String) -> some View {
        Text(text)
            .frame(maxWidth: .infinity, alignment: .trailing)
    }

    /// 23.0: a header of the moments table, e.g. "⟨r^{-2}⟩^{-1/2} (Å)", on one line.
    private func header(_ latex: String) -> some View {
        LaTeXText(latex, size: Self.momentHeaderSize, weight: .semibold)
            .lineLimit(1)
            .fixedSize()
            .frame(maxWidth: .infinity, alignment: .trailing)
    }

    private func momentsTable(_ r: NuclearResult) -> some View {
        let a = NuclearUnits.bohrToAngstrom
        func f(_ x: Double) -> String { String(format: "%.5f", x * a) }
        // 23.0: the quantity and its unit on one line, 15 pt with superscripts (22.0: 13 pt, two
        // lines); 22.0: the seven value columns share the width of the window, n and the bond as
        // wide as their text; 12-pt column gaps (the 21.0 table fitted next to the sidebar with them)
        return Grid(alignment: .trailing, horizontalSpacing: 12, verticalSpacing: 9) {
            GridRow {
                Text("n")
                    .font(.system(size: Self.momentHeaderSize, weight: .semibold))
                Text("bond")
                    .font(.system(size: Self.momentHeaderSize, weight: .semibold))
                header("⟨r⟩ (Å)")
                header("⟨r^{2}⟩^{1/2} (Å)")
                header("⟨r^{3}⟩^{1/3} (Å)")
                header("⟨r^{−1}⟩^{−1} (Å)")
                header("⟨r^{−2}⟩^{−1/2} (Å)")
                header("⟨r^{−3}⟩^{−1/3} (Å)")
                header("⟨HOH⟩ (°)")
            }
            .foregroundStyle(.secondary)
            Divider()
            ForEach(Array(r.levels.enumerated()), id: \.offset) { k, l in
                GridRow {
                    Text("\(k + 1)")
                    Text("r₁")
                    cell(f(l.r1.mean))
                    cell(f(l.r1.rms))
                    cell(f(l.r1.cubic))
                    cell(f(l.r1.inverse))
                    cell(f(l.r1.inverse2))
                    cell(f(l.r1.inverse3))
                    cell(String(format: "%.3f", l.theta))
                }
                .font(Self.valueFont)
                GridRow {
                    Text("")
                    Text("r₂")
                    cell(f(l.r2.mean))
                    cell(f(l.r2.rms))
                    cell(f(l.r2.cubic))
                    cell(f(l.r2.inverse))
                    cell(f(l.r2.inverse2))
                    cell(f(l.r2.inverse3))
                    Text("")
                }
                .font(Self.valueFont)
                .foregroundStyle(.secondary)
            }
        }
        .textSelection(.enabled)
    }

    private func kTable(_ r: NuclearResult, _ kv: [Int]) -> some View {
        Grid(alignment: .trailing, horizontalSpacing: 16, verticalSpacing: 9) {
            GridRow {
                Text("n")
                ForEach(kv, id: \.self) { k in cell("K = \(k)") }
            }
            .font(Self.headerFont)
            .foregroundStyle(.secondary)
            Divider()
            ForEach(Array(r.levels.enumerated()), id: \.offset) { n, l in
                GridRow {
                    Text("\(n + 1)")
                    ForEach(Array(l.kWeights.enumerated()), id: \.offset) { _, w in
                        cell(String(format: "%.5f", w))
                    }
                }
                .font(Self.valueFont)
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
            ScrollViewReader { reader in
                ScrollView {
                    VStack(alignment: .leading, spacing: 18) {
                        NuclearHeader(result: r)
                        SectionTitle("Eigensolver")
                        stats(r)
                            .dataPanel()
                        if let h = r.history, h.count > 1 {
                            SectionTitle("History")
                            HStack(spacing: 18) {
                                lowestChart(h, perLevel: r.input?.eigensolver == .rayleighQuotient)
                                    .frame(minHeight: 240)
                                    .dataPanel()
                                changeChart(h, lanczos: r.input?.eigensolver == .lanczos)
                                    .frame(minHeight: 240)
                                    .dataPanel()
                            }
                        }
                        // 21.0: the linear solves, where the linear solver and the preconditioner show
                        let linear = (r.history ?? []).contains { $0.linearIterations > 0 }
                        if linear {
                            SectionTitle("Linear solves")
                                .id(Self.linearAnchor)
                            HStack(spacing: 18) {
                                linearIterationsChart(r.history ?? [])
                                    .frame(minHeight: 240)
                                    .dataPanel()
                                linearTraceChart(r)
                                    .frame(minHeight: 240)
                                    .dataPanel()
                            }
                            Text(linearCaption(r))
                                .font(.caption)
                                .foregroundStyle(.secondary)
                                .fixedSize(horizontal: false, vertical: true)
                        }
                        SectionTitle("Grids")
                        grids(r)
                            .dataPanel()
                    }
                    .padding(20)
                }
                .task {
                    // UI test (-psilabNuclearConvergenceScroll linear): open at the linear-solve charts
                    // (after the layout of the charts, on the main actor)
                    guard UserDefaults.standard.string(forKey: "psilabNuclearConvergenceScroll") == "linear" else { return }
                    try? await Task.sleep(for: .milliseconds(500))
                    reader.scrollTo(Self.linearAnchor, anchor: .top)
                }
            }
        } else {
            NuclearPlaceholder(model: model)
        }
    }

    private static let linearAnchor = "nuclear-convergence-linear-solves"

    private func linearCaption(_ r: NuclearResult) -> String {
        let input = r.input ?? NuclearInput()
        var text = "Left: the iterations of the linear solves of every step (\(input.linearSolver.longTitle), preconditioner \(input.preconditioner.title)). Right: the residual norm of a linear solve, iteration by iteration (MINRES: its estimate of the preconditioned residual), for the first and the last step."
        if input.eigensolver == .lanczos {
            text += " For Lánczos the outer iteration (the Krylov space of H⁻¹) does not depend on the linear solver or the preconditioner as long as every solve reaches the linear tolerance: the levels and their change (above) are the same for every choice, to rounding; the choices differ in the inner iterations, the products H x and the time."
        }
        if input.linearSolver == .bicg {
            text += " BiCG with r̃₀ = r₀ reproduces the iterates of CG for the symmetric H at twice the products H x."
        }
        return text
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
            InfoRow("Eigenvectors", (r.hasVectors ?? true) ? "computed" : "not kept (eigenvalues only)")
            InfoRow("Hamiltonian", String(format: "%@ × %@, %@ non-zero elements, %.2f s", Fmt.count(Double(r.dimension ?? 0)),
                                          Fmt.count(Double(r.dimension ?? 0)), Fmt.count(Double(r.nonzeros ?? 0)), r.hamiltonianSeconds ?? 0))
            InfoRow("Eigensolver time", String(format: "%.2f s", r.solverSeconds ?? 0))
            if let c = r.chebyshevInterval, c.count == 2 {
                InfoRow("Chebyshev interval", String(format: "eigenvalues of M⁻¹A in [%.3e, %.3e]", c[0], c[1]))
            }
            if let sv = r.startVector {
                InfoRow("Start vector", String(format: "Gaussian (Algorithm 5): μ = %.4f, %.4f bohr, %.2f°; σ = %.4f, %.4f bohr, %.2f°",
                                               sv.mu[0], sv.mu[1], sv.mu[2], sv.sigma[0], sv.sigma[1], sv.sigma[2]))
                InfoRow("Its quality", String(format: "Rayleigh quotient %.2f cm⁻¹%@", sv.rayleigh * NuclearUnits.hartreeToCm,
                                              sv.overlap.map { String(format: ", overlap |⟨q₁|x₁⟩| = %.4f with the lowest eigenvector", $0) } ?? ""))
            }
            if let m = r.masses, m.count == 3, let mu1 = r.mu1, let mu2 = r.mu2 {
                InfoRow("Masses", String(format: "%.9f, %.9f, %.7f u", m[0], m[1], m[2]))
                InfoRow("Coordinates", String(format: "%@: g₁ = %.6f, g₂ = %.6f (Sutcliffe–Tennyson); μ₁ = %.4f, μ₂ = %.4f mₑ",
                                              (r.input?.coordinates ?? .jacobi).title, r.g1 ?? 0, r.g2 ?? 0, mu1, mu2))
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

    private func lowestChart(_ h: [NuclearStep], perLevel: Bool) -> some View {
        let points = h.filter { $0.lowest > 0 }
        return Chart(points, id: \.iteration) { s in
            LineMark(x: .value("Iteration", s.iteration), y: .value("Level (cm⁻¹)", s.lowest * NuclearUnits.hartreeToCm))
                .foregroundStyle(by: .value("Level", perLevel ? "level \(s.level ?? 0)" : "lowest"))
            PointMark(x: .value("Iteration", s.iteration), y: .value("Level (cm⁻¹)", s.lowest * NuclearUnits.hartreeToCm))
                .symbolSize(12)
                .foregroundStyle(by: .value("Level", perLevel ? "level \(s.level ?? 0)" : "lowest"))
        }
        .chartYScale(domain: .automatic(includesZero: false))
        .chartLegend(perLevel ? .visible : .hidden)
        .chartXAxisLabel("Iteration")
        .chartYAxisLabel(perLevel ? "Level being refined (cm⁻¹)" : "Lowest level (cm⁻¹)")
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

    /// 21.0: the iterations of the linear solver(s) of every eigensolver step.
    private func linearIterationsChart(_ h: [NuclearStep]) -> some View {
        Chart(h, id: \.iteration) { s in
            BarMark(x: .value("Iteration", s.iteration), y: .value("Linear-solver iterations", s.linearIterations))
        }
        .chartXAxisLabel("Eigensolver iteration")
        .chartYAxisLabel("Linear-solver iterations")
        .accessibilityLabel("Linear-solver iterations per step")
    }

    /// 21.0: log₁₀ of the residual norm of the first and the last linear solve.
    private func linearTraceChart(_ r: NuclearResult) -> some View {
        struct P: Identifiable {
            let id: Int
            let series: String
            let k: Int
            let y: Double
        }
        var pts: [P] = []
        var n = 0
        for (name, t) in [("first solve (step \(r.firstSolve?.step ?? 0))", r.firstSolve), ("last solve (step \(r.lastSolve?.step ?? 0))", r.lastSolve)] {
            guard let t else { continue }
            for (k, v) in t.residuals.enumerated() where v > 0 && v.isFinite {
                pts.append(P(id: n, series: name, k: k, y: log10(v)))
                n += 1
            }
        }
        return Chart(pts) { p in
            LineMark(x: .value("Inner iteration", p.k), y: .value("log10 residual", p.y))
                .foregroundStyle(by: .value("Solve", p.series))
        }
        .chartXAxisLabel("Linear-solver iteration")
        .chartYAxisLabel("log₁₀ ‖b − A x‖")
        .accessibilityLabel("Residual of the linear solves")
    }

    private func grids(_ r: NuclearResult) -> some View {
        let coords = r.input?.coordinates ?? .jacobi
        return VStack(alignment: .leading, spacing: 10) {
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
                    Text("Angular points, K = \(k): \(coords.title) angle θ (°)")
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
