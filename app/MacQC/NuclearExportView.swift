//
//  NuclearExportView.swift
//  MacQC
//
//  ΨLab 21.0: the export of a nuclear-motion result: the eigenvalues (the level table of the output
//  file: n, E (E_h), E (cm⁻¹), E − E₁ (cm⁻¹), ‖Hx − Ex‖ (E_h)), the eigenvectors (the DVR coefficients
//  of every basis function, with its grid point, sign and weight) and/or the Gaussian start vector
//  of Algorithm 5, together or separately, as one text file (macqc_nuclear_export).
//

import SwiftUI

struct NuclearExportSheet: View {
    var model: NuclearModel
    @Environment(\.dismiss) private var dismiss
    @State private var eigenvalues = true
    @State private var eigenvectors = false
    @State private var startVector = false

    var body: some View {
        let vectorsAvailable = model.canExportVectors
        let startAvailable = model.canExportStartVector
        let r = model.result
        VStack(alignment: .leading, spacing: 14) {
            Text("Export Eigenvalues and Eigenvectors")
                .font(.title3.weight(.semibold))
            if let r {
                Text("\(r.isotopologue ?? "") · \(r.levels.count) levels · \(Fmt.count(Double(r.dimension ?? 0))) basis functions")
                    .foregroundStyle(.secondary)
            }
            VStack(alignment: .leading, spacing: 10) {
                Toggle(isOn: $eigenvalues) {
                    VStack(alignment: .leading, spacing: 2) {
                        Text("Eigenvalues")
                        Text("n, E (E_h), E (cm⁻¹), E − E₁ (cm⁻¹), ‖Hx − Ex‖ (E_h), as in the output file")
                            .font(.caption)
                            .foregroundStyle(.secondary)
                    }
                }
                Toggle(isOn: Binding(get: { eigenvectors && vectorsAvailable }, set: { eigenvectors = $0 })) {
                    VStack(alignment: .leading, spacing: 2) {
                        Text("Eigenvectors")
                        Text(vectorsAvailable
                             ? "the coefficients c_k of every DVR function (one row per basis function, with K, the grid point R₁, R₂, θ, its sign and weight; ψ = c/(s√w))"
                             : "not available: this calculation computed the eigenvalues only")
                            .font(.caption)
                            .foregroundStyle(.secondary)
                            .fixedSize(horizontal: false, vertical: true)
                    }
                }
                .disabled(!vectorsAvailable)
                Toggle(isOn: Binding(get: { startVector && startAvailable }, set: { startVector = $0 })) {
                    VStack(alignment: .leading, spacing: 2) {
                        Text("Gaussian start vector")
                        Text(startAvailable
                             ? "q₁ of Algorithm 5 in the same rows (a column before the eigenvectors when both are exported)"
                             : "not available: this calculation did not use the Gaussian start vector")
                            .font(.caption)
                            .foregroundStyle(.secondary)
                            .fixedSize(horizontal: false, vertical: true)
                    }
                }
                .disabled(!startAvailable)
            }
            .toggleStyle(.checkbox)
            if let r, eigenvectors && vectorsAvailable {
                let mb = Double(r.dimension ?? 0) * Double(r.levels.count + 1) * 23 / 1_048_576
                if mb > 50 {
                    Label(String(format: "About %.0f MB of text.", mb), systemImage: "exclamationmark.triangle")
                        .font(.caption)
                        .foregroundStyle(.orange)
                }
            }
            HStack {
                Spacer()
                Button("Cancel") { dismiss() }
                    .keyboardShortcut(.cancelAction)
                Button("Export…") {
                    let values = eigenvalues
                    let vectors = eigenvectors && vectorsAvailable
                    let start = startVector && startAvailable
                    dismiss()
                    model.exportLevels(eigenvalues: values, eigenvectors: vectors, startVector: start)
                }
                .keyboardShortcut(.defaultAction)
                .disabled(!(eigenvalues || (eigenvectors && vectorsAvailable) || (startVector && startAvailable)))
            }
        }
        .padding(20)
        .frame(width: 480)
    }
}
