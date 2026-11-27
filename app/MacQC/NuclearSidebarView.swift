//
//  NuclearSidebarView.swift
//  MacQC
//
//  ΨLab 20.0/21.0: the sidebar of the nuclear-motion calculation. The basis of the discrete variable
//  representation is set by ns1 (Hermite points in R₁ and in R₂) and np (Legendre points in
//  cos θ) in Jacobi or (21.0) Radau coordinates, the rotation by jmax (J; jmax = 0: the vibrational
//  levels only), as in the original program; 21.0: the eigenvalues alone or with the eigenvectors
//  (the original's sv); the eigensolver (21.0 also inverse power and Rayleigh quotient iteration),
//  its linear solver (21.0 also Chebyshev and BiCG, for every iterative method), preconditioner,
//  the Gaussian start vector of Algorithm 5 and the tolerances in a section of their own.
//

import SwiftUI

struct NuclearSidebarView: View {
    @Bindable var model: NuclearModel

    private static let numberFormat = FloatingPointFormatStyle<Double>().precision(.fractionLength(0...6))

    var body: some View {
        Form {
            Section("Molecule") {
                Picker("Molecule", selection: $model.input.isotopologue) {
                    ForEach(Isotopologue.allCases) { iso in
                        Text(iso.title).tag(iso)
                    }
                }
                .accessibilityLabel("Isotopologue")
                MoleculeCaption(model.input.isotopologue.note)
                MoleculeCaption("CVRQD potential energy surface (Barletta et al., J. Chem. Phys. 125, 204307 (2006)); energies from its minimum")
            }

            Section("Basis (discrete variable representation)") {
                Picker("Coordinates", selection: $model.input.coordinates) {
                    ForEach(NuclearCoordinates.allCases) { c in
                        Text(c.title).tag(c)
                    }
                }
                .pickerStyle(.segmented)
                .accessibilityLabel("Coordinates")
                .help("The orthogonal coordinates of the Hamiltonian (nkoord of the original program: 1 = Jacobi, 2 = Radau)")
                Text(model.input.coordinates.explanation)
                    .font(.caption)
                    .foregroundStyle(.secondary)
                integerRow("Radial points ns1", help: model.input.coordinates == .jacobi
                           ? "Hermite DVR points in R₁ (O–H₁) and in R₂ (H₂ to the centre of mass of O–H₁)"
                           : "Hermite DVR points in R₁ and in R₂ (the Radau distances of H₁ and H₂ from the canonical point)",
                           value: $model.input.radialPoints, range: 2...120)
                integerRow("Angular points np", help: "(Associated) Legendre DVR points in cos θ, the \(model.input.coordinates.title) angle",
                           value: $model.input.angularPoints, range: 2...200)
                LabeledContent("Re (bohr)") {
                    numberField("Re bohr", value: $model.input.rCentre)
                }
                .help("Centre of the radial grids, in bohr as in the original program (Re = 2.000)")
                LabeledContent("Rmax (bohr)") {
                    numberField("Rmax bohr", value: $model.input.rHalfWidth)
                }
                .help("The outermost radial points lie at Re ± Rmax (bohr; the original: Rmax = 0.800)")
                Text(String(format: "Radial grid %.3f … %.3f bohr (%.3f … %.3f Å)",
                            model.input.rCentre - model.input.rHalfWidth, model.input.rCentre + model.input.rHalfWidth,
                            (model.input.rCentre - model.input.rHalfWidth) * NuclearUnits.bohrToAngstrom,
                            (model.input.rCentre + model.input.rHalfWidth) * NuclearUnits.bohrToAngstrom))
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
            .disabled(model.isRunning)

            Section("Rotation") {
                integerRow("jmax (J)", help: "Total angular momentum J; jmax = 0: vibrational levels only",
                           value: $model.input.J, range: 0...30)
                    // 21.0: J = 0 has only the block K = 0; the hidden picker kept K = 1 … J of an
                    // earlier J > 0 (the backend now ignores it for J = 0 as well)
                    .onChange(of: model.input.J) { _, newJ in
                        if newJ == 0 { model.input.parityBlock = 0 }
                    }
                if model.input.J > 0 {
                    Picker("K block", selection: $model.input.parityBlock) {
                        Text("K = 0 … \(model.input.J)").tag(0)
                        Text("K = 1 … \(model.input.J)").tag(1)
                    }
                    .pickerStyle(.segmented)
                    .accessibilityLabel("K block")
                    Text("Parity-adapted (Wang) functions of K, the body-fixed z axis along R₁: the two parity blocks of J (npar = \(model.input.parityBlock) of the original program).")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                } else {
                    Text("jmax = 0: the vibrational energy levels.")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
            }
            .disabled(model.isRunning)

            Section("Levels") {
                integerRow("Number of levels", help: "The lowest eigenvalues computed (nmax)", value: $model.input.levels, range: 1...200)
                // 21.0: the original's sv (0: eigenvalues only, 1: with the eigenvectors); one of the
                // two is always on, switching one off switches the other on
                Toggle("Eigenvalues", isOn: Binding(get: { !model.input.eigenvectors },
                                                    set: { model.input.eigenvectors = !$0 }))
                    .accessibilityLabel("Eigenvalues only")
                    .help("The eigenvalues only (sv = 0 of the original program): Lánczos forms no Ritz vectors and stops when the levels no longer change; no expectation values, no Eigenvectors tab")
                Toggle("Eigenvalues and Eigenvectors", isOn: Binding(get: { model.input.eigenvectors },
                                                                     set: { model.input.eigenvectors = $0 }))
                    .accessibilityLabel("Eigenvalues and eigenvectors")
                    .help("The eigenvectors as well (sv = 1): their residuals ‖Hx − Ex‖, the expectation values and the Eigenvectors tab")
                Toggle("Expectation values", isOn: $model.input.expectationValues)
                    .help("⟨r⟩, ⟨r²⟩^½, … of both O–H bonds and ⟨H–O–H⟩ for every level (needs the eigenvectors)")
                    .disabled(!model.input.eigenvectors)
            }
            .disabled(model.isRunning)

            Section("Eigensolver") {
                Picker("Method", selection: $model.input.eigensolver) {
                    ForEach(NuclearEigensolver.allCases) { s in
                        Text(s.title).tag(s)
                    }
                }
                .accessibilityLabel("Eigensolver")
                Text(model.input.eigensolver.explanation)
                    .font(.caption)
                    .foregroundStyle(.secondary)
                if model.input.eigensolver.usesLinearSolver {
                    Picker("Linear solver", selection: $model.input.linearSolver) {
                        ForEach(NuclearLinearSolver.allCases) { s in
                            Text(s.title).tag(s)
                        }
                    }
                    .pickerStyle(.segmented)
                    .accessibilityLabel("Linear solver")
                    .help(linearSolverHelp)
                    Picker("Preconditioner", selection: $model.input.preconditioner) {
                        ForEach(NuclearPreconditioner.allCases) { p in
                            Text(p.title).tag(p)
                        }
                    }
                    .pickerStyle(.segmented)
                    .accessibilityLabel("Preconditioner")
                    Text(preconditionerCaption)
                        .font(.caption)
                        .foregroundStyle(.secondary)
                    // 21.0: the start vector of Algorithm 5 (the original's psikiind)
                    Toggle("Gaussian start vector", isOn: $model.input.gaussianStart)
                        .accessibilityLabel("Gaussian start vector")
                        .help("Algorithm 5 of the thesis (psikiind): a product of Gaussians in r(O–H₁), r(O–H₂) and ∠HOH at the grid points, centred at the vibrationally averaged geometry (the minimum and a cubic anharmonic shift), with the widths of the harmonic ground state; it approximates the ground state (overlap about 0.997 for H₂O) and speeds up its convergence")
                    Text(model.input.gaussianStart
                         ? "Start vector: Gaussians in the bond lengths and the angle (Algorithm 5), close to the vibrational ground state."
                         : "Start vector: reproducible pseudo-random numbers.")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
                if model.input.eigensolver == .jacobiDavidson {
                    integerRow("Inner iterations", help: "Steps of the linear solver on the correction equation per iteration (0: the preconditioner alone)",
                               value: $model.input.jdInnerIterations, range: 0...200)
                }
                if model.input.eigensolver == .davidson {
                    integerRow("Inner iterations", help: "Steps of the linear solver on (H − θ) t = −r per iteration; 0: the classical correction t = −M⁻¹r with the preconditioner alone. Few steps only: an exact solve gives t = −u, nothing new (the reason for Jacobi–Davidson)",
                               value: $model.input.davidsonInnerIterations, range: 0...200)
                }
                if model.input.eigensolver != .dense {
                    LabeledContent("Eigenvalues (E_h)") {
                        toleranceField("Eigenvalue tolerance", value: $model.input.eigenvalueTolerance)
                    }
                    LabeledContent("Residual ‖Hx − Ex‖ (E_h)") {
                        toleranceField("Residual tolerance", value: $model.input.residualTolerance)
                    }
                    LabeledContent("Linear solves") {
                        toleranceField("Linear tolerance", value: $model.input.linearTolerance)
                    }
                    integerRow("Maximum iterations", help: "Lánczos steps, (Jacobi–)Davidson iterations or power/Rayleigh-quotient iterations (per level)",
                               value: $model.input.maxIterations, range: 10...20000)
                }
            }
            .disabled(model.isRunning)

            Section("Size of the calculation") {
                let e = model.estimate
                LabeledContent("Basis functions") {
                    Text("\(model.input.kValues.count) × \(model.input.radialPoints)² × \(model.input.angularPoints) = \(Fmt.count(Double(e.dimension ?? model.input.dimension)))")
                        .monospacedDigit()
                }
                if let nz = e.nonzeros {
                    LabeledContent("Non-zero elements", value: Fmt.count(Double(nz)))
                }
                if let mb = e.memoryMB {
                    LabeledContent("Memory", value: mb < 1024 ? String(format: "%.0f MB", mb) : String(format: "%.1f GB", mb / 1024))
                }
                if let refusal = e.refusal ?? e.error {
                    Label(refusal, systemImage: "xmark.octagon.fill")
                        .font(.caption)
                        .foregroundStyle(.red)
                }
            }
        }
        .formStyle(.grouped)
        .safeAreaInset(edge: .bottom) {
            runBar
        }
    }

    /// What the linear solver solves for the chosen eigensolver.
    private var linearSolverHelp: String {
        let system: String
        switch model.input.eigensolver {
        case .lanczos, .inversePower: system = "H x = b (H positive definite)"
        case .davidson: system = "(H − θ) t = −r, a few steps"
        case .jacobiDavidson: system = "the projected correction equation (I − QQᵀ)(H − θ)(I − QQᵀ) t = −r"
        case .rayleighQuotient: system = "(H − θ) y = x, nearly singular and indefinite"
        case .dense: system = "nothing"
        }
        return "Solves \(system). CG: preconditioned conjugate gradients; MINRES: minimal residual (indefinite systems); Chebyshev: Chebyshev iteration on the estimated eigenvalue interval of M⁻¹A (no inner products); BiCG: biconjugate gradients (for symmetric H the iterates of CG at twice the products)."
    }

    private var preconditionerCaption: String {
        switch model.input.eigensolver {
        case .lanczos, .inversePower:
            return "ILU: the symmetric Gauss–Seidel factors (D + L) D⁻¹ (D + U) of H, as in the original program; Jacobi: the diagonal D."
        default:
            return "ILU: the symmetric Gauss–Seidel factors of H − θ (D − θ on the diagonal; |D − θ| where the solver needs M positive definite), as the original's ILUprecondfordavidson; Jacobi: D − θ."
        }
    }

    private func integerRow(_ title: String, help: String, value: Binding<Int>, range: ClosedRange<Int>) -> some View {
        LabeledContent(title) {
            HStack(spacing: 4) {
                TextField(title, value: Binding(get: { value.wrappedValue },
                                                set: { value.wrappedValue = min(max($0, range.lowerBound), range.upperBound) }),
                          format: .number)
                    .labelsHidden()
                    .multilineTextAlignment(.trailing)
                    .monospacedDigit()
                    .frame(width: 64)
                Stepper(title, value: value, in: range)
                    .labelsHidden()
            }
        }
        .help(help)
        .accessibilityLabel(title)
    }

    private func numberField(_ title: String, value: Binding<Double>) -> some View {
        TextField(title, value: value, format: Self.numberFormat)
            .labelsHidden()
            .multilineTextAlignment(.trailing)
            .monospacedDigit()
            .frame(width: 96)
    }

    private func toleranceField(_ title: String, value: Binding<Double>) -> some View {
        TextField(title, value: value, format: .number.notation(.scientific).precision(.significantDigits(1...3)))
            .labelsHidden()
            .multilineTextAlignment(.trailing)
            .monospacedDigit()
            .frame(width: 96)
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
                    .disabled(model.isRunning || model.refusal != nil)

                    if model.isRunning {
                        Button("Cancel") {
                            model.cancel()
                        }
                        .secondaryActionStyle()
                        .controlSize(.large)
                        .accessibilityLabel("Cancel calculation")
                        .keyboardShortcut(".", modifiers: [.command])
                    }
                }
                if model.isRunning {
                    TimelineView(.periodic(from: .now, by: 1)) { context in
                        Text((model.progressLine.isEmpty ? "Starting…" : model.progressLine) + "  ·  " + model.elapsedText(at: context.date))
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
                }
            }
            .padding(12)
            .glassPanel(cornerRadius: 18)
        }
        .padding(.horizontal, 12)
        .padding(.bottom, 10)
    }
}
