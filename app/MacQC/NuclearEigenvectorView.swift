//
//  NuclearEigenvectorView.swift
//  MacQC
//
//  ΨLab 21.0: the Eigenvectors tab of the nuclear-motion calculation. The backend sends for every
//  level (with the eigenvectors, at most the first 60) the wavefunction ψ = c/(s√w) at the DVR grid
//  points (c the coefficient of a DVR function, s its sign, w its quadrature weight): cuts through
//  the point of the largest |ψ| in the planes (R₁, R₂), (R₁, θ) and (R₂, θ), and the probability
//  densities integrated over the third coordinate. Drawn as heat maps (a cell around every grid
//  point, to the midpoints between the points), with the one-dimensional densities ρ(R₁), ρ(R₂)
//  and ρ(θ) below. The nodes of ψ (its changes of sign) show the vibrational quantum numbers.
//

import SwiftUI
import Charts

private enum WavePlane: String, CaseIterable, Identifiable {
    case r1r2, r1Theta, r2Theta
    var id: String { rawValue }
    var title: String {
        switch self {
        case .r1r2: return "R₁ – R₂"
        case .r1Theta: return "R₁ – θ"
        case .r2Theta: return "R₂ – θ"
        }
    }
}

private enum WaveQuantity: String, CaseIterable, Identifiable {
    case psi, density
    var id: String { rawValue }
    var title: String { self == .psi ? "Wavefunction ψ (cut)" : "Probability density" }
}

/// A field on a (non-uniform) grid: columns at `xs`, rows at `ys`, value(ix, iy) = values[ix * ys.count + iy].
private struct WaveField {
    var xs: [Double]
    var ys: [Double]
    var values: [Double]
    var xLabel: String
    var yLabel: String
    var valueLabel: String
    var diverging: Bool
    var scale: Double          // largest |value|
    var marker: (x: Double, y: Double)?

    func value(_ ix: Int, _ iy: Int) -> Double { values[ix * ys.count + iy] }

    /// The edges of the cell of every point: the midpoints to its neighbours (in sorted order).
    static func edges(_ v: [Double], clamp: ClosedRange<Double>? = nil) -> [(Double, Double)] {
        let order = v.indices.sorted { v[$0] < v[$1] }
        var out = Array(repeating: (0.0, 0.0), count: v.count)
        for (p, idx) in order.enumerated() {
            let x = v[idx]
            let prev = p > 0 ? v[order[p - 1]] : nil
            let next = p + 1 < order.count ? v[order[p + 1]] : nil
            var lo = prev.map { (x + $0) / 2 } ?? (next.map { x - ($0 - x) / 2 } ?? x - 0.5)
            var hi = next.map { (x + $0) / 2 } ?? (prev.map { x + (x - $0) / 2 } ?? x + 0.5)
            if let c = clamp {
                lo = max(c.lowerBound, lo)
                hi = min(c.upperBound, hi)
            }
            out[idx] = (lo, hi)
        }
        return out
    }
}

struct NuclearEigenvectorView: View {
    var model: NuclearModel
    @State private var level = 0
    @State private var plane: WavePlane = .r1r2
    @State private var quantity: WaveQuantity = .psi
    @State private var hover = ""
    @Environment(\.colorScheme) private var colorScheme

    var body: some View {
        if let r = model.result, !r.levels.isEmpty {
            if r.hasVectors == false {
                VStack(alignment: .leading, spacing: 18) {
                    NuclearHeader(result: r)
                    ContentUnavailableView("Eigenvalues only", systemImage: "waveform",
                                           description: Text("This calculation computed the eigenvalues only. Choose “Eigenvalues and Eigenvectors” under Levels in the sidebar and run again to see the eigenvectors."))
                }
                .padding(20)
                .frame(maxWidth: .infinity, maxHeight: .infinity, alignment: .top)
            } else {
                content(r)
            }
        } else {
            NuclearPlaceholder(model: model)
        }
    }

    // MARK: Layout

    private func content(_ r: NuclearResult) -> some View {
        let available = r.levels.indices.filter { r.levels[$0].wavefunction != nil }
        let k = available.contains(level) ? level : (available.first ?? 0)
        return ScrollView {
            VStack(alignment: .leading, spacing: 16) {
                NuclearHeader(result: r)
                controls(r, available: available, shown: k)
                if let w = r.levels[k].wavefunction {
                    let field = makeField(r, w)
                    HStack(alignment: .top, spacing: 18) {
                        heatmap(field)
                            .frame(minWidth: 360, minHeight: 380)
                            .dataPanel()
                        info(r, k, w, field)
                            .frame(width: 250, alignment: .leading)
                            .dataPanel()
                    }
                    SectionTitle("One-dimensional probability densities")
                    marginals(r, w)
                    Text(caption(r, w))
                        .font(.caption)
                        .foregroundStyle(.secondary)
                        .fixedSize(horizontal: false, vertical: true)
                } else {
                    ContentUnavailableView("No picture of this level", systemImage: "waveform",
                                           description: Text("The eigenvectors are drawn for the first 60 levels."))
                }
            }
            .padding(20)
        }
    }

    private func controls(_ r: NuclearResult, available: [Int], shown: Int) -> some View {
        let e1 = r.levels.first?.energy ?? 0
        return ScrollView(.horizontal, showsIndicators: false) {
            HStack(spacing: 12) {
                Picker("Level", selection: Binding(get: { shown }, set: { level = $0 })) {
                    ForEach(available, id: \.self) { i in
                        Text(String(format: "%d  ·  %.2f cm⁻¹", i + 1, (r.levels[i].energy - e1) * NuclearUnits.hartreeToCm)).tag(i)
                    }
                }
                .fixedSize()
                .accessibilityLabel("Level")
                Stepper("Level", value: Binding(get: { shown }, set: { level = min(max($0, available.first ?? 0), available.last ?? 0) }),
                        in: (available.first ?? 0)...(available.last ?? 0))
                    .labelsHidden()
                Picker("Plane", selection: $plane) {
                    ForEach(WavePlane.allCases) { p in Text(p.title).tag(p) }
                }
                .pickerStyle(.segmented)
                .labelsHidden()
                .fixedSize()
                .accessibilityLabel("Plane")
                Picker("Quantity", selection: $quantity) {
                    ForEach(WaveQuantity.allCases) { q in Text(q.title).tag(q) }
                }
                .pickerStyle(.segmented)
                .labelsHidden()
                .fixedSize()
                .accessibilityLabel("Quantity")
            }
        }
    }

    // MARK: The field of the chosen plane and quantity

    private func thetaDegrees(_ r: NuclearResult, block: Int) -> [Double] {
        guard let c = r.cosGrid, block < c.count else { return [] }
        return c[block].map { acos(max(-1, min(1, $0))) * 180 / .pi }
    }

    private func makeField(_ r: NuclearResult, _ w: NuclearWavefunction) -> WaveField {
        let a = NuclearUnits.bohrToAngstrom
        let rs = (r.rGrid ?? []).map { $0 * a }
        let th = thetaDegrees(r, block: w.block)
        let psi = quantity == .psi
        let coords = r.input?.coordinates ?? .jacobi
        let rName = coords == .jacobi ? "" : " from K"
        let xs = rs
        var ys = rs
        var xLabel = "R₁\(rName) (Å)", yLabel = "R₂\(rName) (Å)"
        let raw: String
        var marker: (x: Double, y: Double)?
        switch plane {
        case .r1r2:
            raw = psi ? w.psiR1R2 : w.rhoR1R2
            if w.i < rs.count, w.j < rs.count { marker = (rs[w.i], rs[w.j]) }
        case .r1Theta:
            raw = psi ? w.psiR1Theta : w.rhoR1Theta
            ys = th
            yLabel = "θ (°)"
            if w.i < rs.count, w.l < th.count { marker = (rs[w.i], th[w.l]) }
        case .r2Theta:
            raw = psi ? w.psiR2Theta : w.rhoR2Theta
            xLabel = "R₂\(rName) (Å)"
            ys = th
            yLabel = "θ (°)"
            if w.j < rs.count, w.l < th.count { marker = (rs[w.j], th[w.l]) }
        }
        var values = NuclearWavefunction.floats(raw).map(Double.init)
        // densities per Å (and per degree): the backend's are per bohr and per radian
        if !psi {
            let f = plane == .r1r2 ? 1 / (a * a) : (1 / a) * (.pi / 180)
            values = values.map { $0 * f }
        }
        if values.count != xs.count * ys.count { values = Array(repeating: 0, count: xs.count * ys.count) }
        let scale = values.map(abs).max() ?? 1
        let valueLabel: String
        if psi {
            valueLabel = "ψ / ψ_max"
        } else {
            valueLabel = plane == .r1r2 ? "ρ (Å⁻²)" : "ρ (Å⁻¹ deg⁻¹)"
        }
        return WaveField(xs: xs, ys: ys, values: values, xLabel: xLabel, yLabel: yLabel, valueLabel: valueLabel,
                         diverging: psi, scale: scale > 0 ? scale : 1, marker: marker)
    }

    // MARK: Colours (diverging blue – gray – red for ψ; one blue ramp for ρ)

    private func rgb(_ hex: UInt32) -> (Double, Double, Double) {
        (Double((hex >> 16) & 0xff) / 255, Double((hex >> 8) & 0xff) / 255, Double(hex & 0xff) / 255)
    }

    private func mix(_ a: (Double, Double, Double), _ b: (Double, Double, Double), _ t: Double) -> Color {
        let t = max(0, min(1, t))
        return Color(red: a.0 + (b.0 - a.0) * t, green: a.1 + (b.1 - a.1) * t, blue: a.2 + (b.2 - a.2) * t)
    }

    private func color(_ v: Double, _ f: WaveField) -> Color {
        let dark = colorScheme == .dark
        let t = v / f.scale
        if f.diverging {
            let mid = rgb(dark ? 0x383838 : 0xefefef)
            let neg = rgb(dark ? 0x519ced : 0x1b64be)
            let pos = rgb(dark ? 0xeb6263 : 0xd23b3c)
            return t < 0 ? mix(mid, neg, -t) : mix(mid, pos, t)
        }
        let lo = rgb(dark ? 0x1f2630 : 0xf5f9ff)
        let hi = rgb(dark ? 0xb6d5f8 : 0x0e3d7e)
        return mix(lo, hi, t)
    }

    // MARK: Heat map

    private func heatmap(_ f: WaveField) -> some View {
        let xe = WaveField.edges(f.xs)
        let ye = WaveField.edges(f.ys, clamp: f.yLabel.hasPrefix("θ") ? 0...180 : nil)
        let xDomain = (xe.map(\.0).min() ?? 0)...(xe.map(\.1).max() ?? 1)
        let yDomain = (ye.map(\.0).min() ?? 0)...(ye.map(\.1).max() ?? 1)
        return VStack(alignment: .leading, spacing: 8) {
            Chart {
                // invisible corners: the axes; the cells are drawn below the plot
                PointMark(x: .value(f.xLabel, xDomain.lowerBound), y: .value(f.yLabel, yDomain.lowerBound)).opacity(0)
                PointMark(x: .value(f.xLabel, xDomain.upperBound), y: .value(f.yLabel, yDomain.upperBound)).opacity(0)
                if let m = f.marker {
                    PointMark(x: .value(f.xLabel, m.x), y: .value(f.yLabel, m.y))
                        .symbol(.cross)
                        .symbolSize(90)
                        .foregroundStyle(Color.primary)
                }
            }
            .chartXScale(domain: xDomain)
            .chartYScale(domain: yDomain)
            .chartXAxisLabel(f.xLabel)
            .chartYAxisLabel(f.yLabel)
            .chartBackground { proxy in
                GeometryReader { geo in
                    if let pf = proxy.plotFrame {
                        let frame = geo[pf]
                        Canvas { ctx, _ in
                            for ix in f.xs.indices {
                                guard let x0 = proxy.position(forX: xe[ix].0), let x1 = proxy.position(forX: xe[ix].1) else { continue }
                                for iy in f.ys.indices {
                                    guard let y0 = proxy.position(forY: ye[iy].0), let y1 = proxy.position(forY: ye[iy].1) else { continue }
                                    let rect = CGRect(x: min(x0, x1), y: min(y0, y1), width: abs(x1 - x0) + 0.5, height: abs(y1 - y0) + 0.5)
                                    ctx.fill(Path(rect), with: .color(self.color(f.value(ix, iy), f)))
                                }
                            }
                        }
                        .frame(width: frame.width, height: frame.height)
                        .offset(x: frame.minX, y: frame.minY)
                    }
                }
            }
            .chartOverlay { proxy in
                GeometryReader { geo in
                    Rectangle()
                        .fill(Color.clear)
                        .contentShape(Rectangle())
                        .onContinuousHover { phase in
                            switch phase {
                            case .active(let location):
                                guard let pf = proxy.plotFrame else { return }
                                let o = geo[pf].origin
                                guard let x = proxy.value(atX: location.x - o.x, as: Double.self),
                                      let y = proxy.value(atY: location.y - o.y, as: Double.self),
                                      let ix = xe.indices.first(where: { xe[$0].0 <= x && x < xe[$0].1 }),
                                      let iy = ye.indices.first(where: { ye[$0].0 <= y && y < ye[$0].1 }) else {
                                    hover = ""
                                    return
                                }
                                hover = String(format: "%@ = %.4f, %@ = %.4f:  %@ = %.4g", f.xLabel, f.xs[ix], f.yLabel, f.ys[iy],
                                               f.valueLabel, f.value(ix, iy))
                            case .ended:
                                hover = ""
                            }
                        }
                }
            }
            .accessibilityLabel(quantity == .psi ? "Wavefunction cut" : "Probability density")
            // the colour scale on its own line (in a narrow panel the readout squeezed its labels
            // to one character), the readout of the cell under the pointer below it
            colorBar(f)
            Text(hover.isEmpty ? "Point at a cell for its value; + marks the largest |ψ|." : hover)
                .font(.caption.monospacedDigit())
                .foregroundStyle(.secondary)
                .lineLimit(2)
                .frame(maxWidth: .infinity, alignment: .leading)
        }
    }

    private func colorBar(_ f: WaveField) -> some View {
        let stops = (0...20).map { k -> Gradient.Stop in
            let t = Double(k) / 20
            let v = f.diverging ? (2 * t - 1) * f.scale : t * f.scale
            return Gradient.Stop(color: color(v, f), location: t)
        }
        return HStack(spacing: 6) {
            Text(f.diverging ? "−1" : "0")
                .fixedSize()
            RoundedRectangle(cornerRadius: 3)
                .fill(LinearGradient(stops: stops, startPoint: .leading, endPoint: .trailing))
                .frame(minWidth: 80, idealWidth: 160, maxWidth: 200)
                .frame(height: 10)
            Text(f.diverging ? "+1" : String(format: "%.3g", f.scale))
                .fixedSize()
            Text(f.valueLabel)
                .foregroundStyle(.secondary)
                .fixedSize()
        }
        .font(.caption.monospacedDigit())
        .accessibilityLabel("Colour scale")
    }

    // MARK: Information and the one-dimensional densities

    private func info(_ r: NuclearResult, _ k: Int, _ w: NuclearWavefunction, _ f: WaveField) -> some View {
        let e1 = r.levels.first?.energy ?? 0
        let l = r.levels[k]
        let a = NuclearUnits.bohrToAngstrom
        let rs = r.rGrid ?? []
        let th = thetaDegrees(r, block: w.block)
        let kv = r.kValues ?? [0]
        return VStack(alignment: .leading, spacing: 8) {
            Text("Level \(k + 1)")
                .font(.title3.weight(.semibold))
            Text(String(format: "E = %.4f cm⁻¹", l.energy * NuclearUnits.hartreeToCm))
            Text(String(format: "E − E₁ = %.4f cm⁻¹", (l.energy - e1) * NuclearUnits.hartreeToCm))
            if let res = l.residual {
                Text(String(format: "‖Hx − Ex‖ = %.1e E_h", res))
                    .foregroundStyle(.secondary)
            }
            Divider()
            Text("Largest |ψ| at")
                .font(.caption.weight(.semibold))
                .foregroundStyle(.secondary)
            if w.i < rs.count, w.j < rs.count, w.l < th.count {
                Text(String(format: "R₁ = %.4f Å\nR₂ = %.4f Å\nθ = %.2f°", rs[w.i] * a, rs[w.j] * a, th[w.l]))
                    .monospacedDigit()
            }
            if kv.count > 1, w.block < kv.count {
                Text("K = \(kv[w.block])" + (w.block < l.kWeights.count ? String(format: " (weight %.3f)", l.kWeights[w.block]) : ""))
                Text("The cuts and the θ planes in this K block; R₁–R₂ density summed over K.")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
            if !l.kWeights.isEmpty, r.input?.expectationValues != false, l.r1.mean > 0 {
                Divider()
                Text(String(format: "⟨r₁⟩ = %.4f Å\n⟨r₂⟩ = %.4f Å\n⟨HOH⟩ = %.2f°", l.r1.mean * a, l.r2.mean * a, l.theta))
                    .monospacedDigit()
                    .foregroundStyle(.secondary)
            }
        }
        .textSelection(.enabled)
    }

    private struct Curve: Identifiable {
        let id: Int
        let x: Double
        let y: Double
    }

    private func marginals(_ r: NuclearResult, _ w: NuclearWavefunction) -> some View {
        let a = NuclearUnits.bohrToAngstrom
        let rs = r.rGrid ?? []
        let wr = r.rWeights ?? []
        let ns = rs.count
        let th = thetaDegrees(r, block: w.block)
        let np = th.count
        let d12 = NuclearWavefunction.floats(w.rhoR1R2).map(Double.init)
        let d1t = NuclearWavefunction.floats(w.rhoR1Theta).map(Double.init)
        var c1: [Curve] = [], c2: [Curve] = [], ct: [Curve] = []
        if wr.count == ns, d12.count == ns * ns {
            for i in 0..<ns {
                var s1 = 0.0, s2 = 0.0
                for j in 0..<ns {
                    s1 += d12[i * ns + j] * wr[j]
                    s2 += d12[j * ns + i] * wr[j]
                }
                c1.append(Curve(id: i, x: rs[i] * a, y: s1 / a))
                c2.append(Curve(id: i, x: rs[i] * a, y: s2 / a))
            }
        }
        if wr.count == ns, d1t.count == ns * np {
            for l in 0..<np {
                var s = 0.0
                for i in 0..<ns { s += d1t[i * np + l] * wr[i] }
                ct.append(Curve(id: l, x: th[l], y: s * .pi / 180))
            }
            ct.sort { $0.x < $1.x }
        }
        let rName = (r.input?.coordinates ?? .jacobi) == .jacobi ? "" : " from K"
        return HStack(spacing: 18) {
            curve(c1, x: "R₁\(rName) (Å)", y: "ρ(R₁) (Å⁻¹)")
            curve(c2, x: "R₂\(rName) (Å)", y: "ρ(R₂) (Å⁻¹)")
            curve(ct, x: "θ (°)", y: "ρ(θ) (deg⁻¹)")
        }
    }

    private func curve(_ c: [Curve], x: String, y: String) -> some View {
        Chart(c) { p in
            LineMark(x: .value(x, p.x), y: .value(y, p.y))
                .lineStyle(StrokeStyle(lineWidth: 2))
                .interpolationMethod(.monotone)
            PointMark(x: .value(x, p.x), y: .value(y, p.y))
                .symbolSize(16)
        }
        .chartXScale(domain: .automatic(includesZero: false))
        .chartXAxisLabel(x)
        .chartYAxisLabel(y)
        .frame(minHeight: 200)
        .dataPanel()
        .accessibilityLabel(y)
    }

    private func caption(_ r: NuclearResult, _ w: NuclearWavefunction) -> String {
        let coords = r.input?.coordinates ?? .jacobi
        var text = "\(coords.title) coordinates. "
        text += quantity == .psi
            ? "ψ at the grid points (the DVR coefficient over the sign and the square root of the quadrature weight of its DVR function), scaled to +1 at its largest value; the cut through that point (+). The overall sign of an eigenvector is arbitrary; the nodes (blue–red boundaries) count the vibrational quanta along each coordinate. "
            : "|ψ|² integrated over the third coordinate (in the R₁–R₂ plane also over the K blocks). "
        text += "Each cell reaches to the midpoints between the grid points (Hermite points in R, Legendre points in cos θ). The one-dimensional densities integrate the R₁–R₂ density over the other distance, and ρ(θ) over R₁ and R₂ (with the radial quadrature weights)."
        return text
    }
}
