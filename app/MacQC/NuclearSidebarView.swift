//
//  NuclearSidebarView.swift
//  MacQC
//
//  ΨLab 20.0: the sidebar of the nuclear-motion calculation. The basis of the discrete variable
//  representation is set by ns1 (Hermite points in R₁ and in R₂) and np (Legendre points in
//  cos θ), the rotation by jmax (J; jmax = 0: the vibrational levels only), as in the original
//  program; the eigensolver and its tolerances are in a section of their own.
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
                integerRow("Radial points ns1", help: "Hermite DVR points in R₁ (O–H₁) and in R₂ (H₂ to the centre of mass of O–H₁)",
                           value: $model.input.radialPoints, range: 2...120)
                integerRow("Angular points np", help: "(Associated) Legendre DVR points in cos θ, the Jacobi angle",
                           value: $model.input.angularPoints, range: 2...200)
                LabeledContent("Re (bohr)") {
                    numberField("Re bohr", value: $model.input.rCentre)
                }
                .help("Centre of the radial grids")
                LabeledContent("Rmax (bohr)") {
                    numberField("Rmax bohr", value: $model.input.rHalfWidth)
                }
                .help("The outermost radial points lie at Re ± Rmax")
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
                Toggle("Expectation values", isOn: $model.input.expectationValues)
                    .help("⟨r⟩, ⟨r²⟩^½, … of both O–H bonds and ⟨H–O–H⟩ for every level")
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
                if model.input.eigensolver == .lanczos {
                    Picker("Linear solver", selection: $model.input.linearSolver) {
                        ForEach(NuclearLinearSolver.allCases) { s in
                            Text(s.title).tag(s)
                        }
                    }
                    .pickerStyle(.segmented)
                    Picker("Preconditioner", selection: $model.input.preconditioner) {
                        ForEach(NuclearPreconditioner.allCases) { p in
                            Text(p.title).tag(p)
                        }
                    }
                    .pickerStyle(.segmented)
                    Text("ILU: the symmetric Gauss–Seidel factors (D + L) D⁻¹ (D + U) of H, as in the original program; Jacobi: the diagonal.")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
                if model.input.eigensolver == .jacobiDavidson {
                    integerRow("MINRES steps", help: "Steps of the correction equation per iteration",
                               value: $model.input.jdInnerIterations, range: 0...200)
                }
                if model.input.eigensolver != .dense {
                    LabeledContent("Eigenvalues (E_h)") {
                        toleranceField("Eigenvalue tolerance", value: $model.input.eigenvalueTolerance)
                    }
                    LabeledContent("Residual ‖Hx − Ex‖ (E_h)") {
                        toleranceField("Residual tolerance", value: $model.input.residualTolerance)
                    }
                    if model.input.eigensolver == .lanczos {
                        LabeledContent("Linear solves") {
                            toleranceField("Linear tolerance", value: $model.input.linearTolerance)
                        }
                    }
                    integerRow("Maximum iterations", help: "Lánczos steps or (Jacobi–)Davidson iterations",
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
