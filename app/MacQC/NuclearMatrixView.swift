//
//  NuclearMatrixView.swift
//  MacQC
//
//  ΨLab 20.0: the structure of the sparse Hamiltonian of the nuclear-motion calculation. The
//  backend draws the matrix into at most 1024 × 1024 pixels (MatrixPicture): the whole matrix (a
//  pixel is a block of elements when the dimension is larger) and its leading rows and columns
//  (one pixel per element). Shown as the non-zero pattern (the share of non-zero elements of a
//  pixel) or as log₁₀|Hᵢⱼ|, as large as the tab allows and zoomable.
//

import SwiftUI
import AppKit

nonisolated struct NuclearMatrixPicture: Codable, Sendable, Equatable {
    var bins: Int
    var first: Int
    var span: Int
    var density: String     // base64, bins × bins bytes
    var magnitude: String   // base64, bins × bins bytes
    var logMin: Double
    var logMax: Double
}

private enum MatrixPart: String, CaseIterable, Identifiable {
    case whole, leading
    var id: String { rawValue }
}

private enum MatrixColouring: String, CaseIterable, Identifiable {
    case pattern, magnitude
    var id: String { rawValue }
    var title: String { self == .pattern ? "Non-zero" : "log₁₀ |Hᵢⱼ|" }
}

struct NuclearMatrixView: View {
    var model: NuclearModel
    @State private var part: MatrixPart = .whole
    @State private var colouring: MatrixColouring = .pattern
    @State private var zoom: CGFloat = 1
    @State private var pinchStart: CGFloat?
    @State private var hover: String = ""
    @State private var showInfo = false
    @Environment(\.colorScheme) private var colorScheme

    var body: some View {
        if let r = model.result, let whole = r.matrixWhole {
            let picture = (part == .leading ? r.matrixDetail : nil) ?? whole
            VStack(spacing: 0) {
                toolbar(r, picture)
                Divider()
                GeometryReader { geo in
                    let fit = max(200, min(geo.size.width, geo.size.height) - 48)
                    let side = fit * zoom
                    ScrollView([.horizontal, .vertical]) {
                        matrix(r, picture, side: side)
                            .padding(24)
                            .frame(minWidth: geo.size.width, minHeight: geo.size.height)
                    }
                    .gesture(MagnifyGesture()
                        .onChanged { value in
                            if pinchStart == nil { pinchStart = zoom }
                            zoom = min(16, max(1, (pinchStart ?? 1) * value.magnification))
                        }
                        .onEnded { _ in pinchStart = nil })
                }
                Divider()
                footer(r, picture)
            }
        } else if model.result != nil {
            ContentUnavailableView("No matrix picture", systemImage: "square.grid.3x3",
                                   description: Text("Run the calculation again to see the structure of its Hamiltonian."))
        } else {
            ContentUnavailableView("No Hamiltonian yet", systemImage: "square.grid.3x3",
                                   description: Text("The structure of the sparse Hamiltonian matrix appears here after a calculation: its non-zero elements, or the size of every element on a logarithmic scale."))
        }
    }

    // MARK: Controls

    private func toolbar(_ r: NuclearResult, _ p: NuclearMatrixPicture) -> some View {
        // scrolls sideways in a narrow window instead of widening it (which squeezed the sidebar)
        ScrollView(.horizontal, showsIndicators: false) {
            toolbarContent(r, p)
        }
    }

    private func toolbarContent(_ r: NuclearResult, _ p: NuclearMatrixPicture) -> some View {
        HStack(spacing: 12) {
            Picker("Part", selection: $part) {
                Text("Whole").tag(MatrixPart.whole)
                Text("First \(Fmt.count(Double(r.matrixDetail?.span ?? p.span)))").tag(MatrixPart.leading)
            }
            .pickerStyle(.segmented)
            .labelsHidden()
            .fixedSize()
            .accessibilityLabel("Part of the matrix")
            Picker("Colour", selection: $colouring) {
                ForEach(MatrixColouring.allCases) { c in
                    Text(c.title).tag(c)
                }
            }
            .pickerStyle(.segmented)
            .labelsHidden()
            .fixedSize()
            .accessibilityLabel("Colouring")
            Spacer(minLength: 8)
            Image(systemName: "minus.magnifyingglass")
                .foregroundStyle(.secondary)
            Slider(value: $zoom, in: 1...16)
                .frame(width: 100)
                .accessibilityLabel("Zoom")
            Image(systemName: "plus.magnifyingglass")
                .foregroundStyle(.secondary)
            Button("Fit") { zoom = 1 }
                .disabled(zoom == 1)
            Button {
                showInfo.toggle()
            } label: {
                Image(systemName: "info.circle")
            }
            .buttonStyle(.borderless)
            .help("How the matrix is ordered")
            .popover(isPresented: $showInfo, arrowEdge: .bottom) {
                info(r)
                    .padding(16)
                    .frame(width: 420)
            }
        }
        .padding(.horizontal, 14)
        .padding(.vertical, 8)
    }

    private func footer(_ r: NuclearResult, _ p: NuclearMatrixPicture) -> some View {
        let n = r.dimension ?? p.span
        let nz = r.nonzeros ?? 0
        let fill = n > 0 ? Double(nz) / (Double(n) * Double(n)) : 0
        let pixel = p.span > p.bins ? String(format: "1 px = %.1f × %.1f elements", Double(p.span) / Double(p.bins), Double(p.span) / Double(p.bins))
                                    : "1 px = 1 element"
        return HStack(spacing: 16) {
            Text("\(Fmt.count(Double(n)))² · \(Fmt.count(Double(nz))) non-zero (\(String(format: "%.3g", 100 * fill)) %, \(String(format: "%.1f", Double(nz) / Double(max(n, 1)))) per row) · \(pixel)")
                .foregroundStyle(.secondary)
                .lineLimit(1)
                .truncationMode(.middle)
                .minimumScaleFactor(0.85)
            Spacer(minLength: 8)
            if colouring == .magnitude {
                legend(p)
            }
            Text(hover)
                .monospacedDigit()
                .foregroundStyle(.secondary)
                .lineLimit(1)
                .frame(minWidth: 180, alignment: .trailing)
        }
        .font(.caption)
        .padding(.horizontal, 14)
        .padding(.vertical, 7)
    }

    private func legend(_ p: NuclearMatrixPicture) -> some View {
        HStack(spacing: 6) {
            Text(String(format: "10^%.1f", p.logMin))
            LinearGradient(colors: (0...6).map { Color(cgColor: Self.heat(Double($0) / 6)) }, startPoint: .leading, endPoint: .trailing)
                .frame(width: 120, height: 10)
                .clipShape(RoundedRectangle(cornerRadius: 2))
            Text(String(format: "10^%.1f E_h", p.logMax))
        }
    }

    private func info(_ r: NuclearResult) -> some View {
        let input = r.input ?? NuclearInput()
        let ns = input.radialPoints, np = input.angularPoints
        return VStack(alignment: .leading, spacing: 8) {
            Text("Order of the basis")
                .font(.headline)
            Text("Index = ((b · ns1 + i) · ns1 + j) · np + l, with b the K block, i the point in R₁, j the point in R₂ and l the angular point (the fastest).")
            Text("• Dense \(np) × \(np) blocks on the diagonal: j² (the angular kinetic energy) couples all angular points of one (i, j).")
            Text("• Diagonals \(np) apart within blocks of \(ns * np): the kinetic energy in R₂ couples j and j′.")
            Text("• Diagonals \(ns * np) apart: the kinetic energy in R₁ couples i and i′.")
            if (r.kValues?.count ?? 1) > 1 {
                Text("• Off-diagonal \(ns * ns * np) × \(ns * ns * np) blocks: the Coriolis coupling of K and K ± 1 (lines: the K blocks).")
            }
            Text("Non-zero elements per row: (ns1 − 1) + (ns1 − 1) + np = \(2 * ns + np - 2)\((r.kValues?.count ?? 1) > 1 ? ", plus np for each neighbouring K block" : "").")
        }
        .font(.callout)
        .fixedSize(horizontal: false, vertical: true)
    }

    // MARK: Picture

    private func matrix(_ r: NuclearResult, _ p: NuclearMatrixPicture, side: CGFloat) -> some View {
        let image = Self.image(p, colouring: colouring, dark: colorScheme == .dark)
        let nx = (r.input?.radialPoints ?? 0) * (r.input?.radialPoints ?? 0) * (r.input?.angularPoints ?? 0)
        let blocks = (r.kValues?.count ?? 1)
        return ZStack(alignment: .topLeading) {
            if let image {
                Image(decorative: image, scale: 1)
                    .resizable()
                    .interpolation(.none)
                    .frame(width: side, height: side)
            }
            // the K blocks
            if blocks > 1 && nx > 0 {
                Path { path in
                    for b in 1..<blocks {
                        let k = Double(b * nx - p.first) / Double(p.span)
                        guard k > 0 && k < 1 else { continue }
                        let x = side * CGFloat(k)
                        path.move(to: CGPoint(x: x, y: 0))
                        path.addLine(to: CGPoint(x: x, y: side))
                        path.move(to: CGPoint(x: 0, y: x))
                        path.addLine(to: CGPoint(x: side, y: x))
                    }
                }
                .stroke(Color.orange.opacity(0.8), style: StrokeStyle(lineWidth: 1, dash: [6, 4]))
                .frame(width: side, height: side)
                .allowsHitTesting(false)
            }
        }
        .frame(width: side, height: side)
        .overlay(Rectangle().strokeBorder(Color.secondary.opacity(0.5), lineWidth: 1))
        .onContinuousHover { phase in
            switch phase {
            case .active(let point):
                hover = Self.hoverText(p, at: point, side: side)
            case .ended:
                hover = ""
            }
        }
        .accessibilityLabel("Structure of the Hamiltonian matrix")
    }

    private static func hoverText(_ p: NuclearMatrixPicture, at point: CGPoint, side: CGFloat) -> String {
        guard side > 0, let density = Data(base64Encoded: p.density) else { return "" }
        let bi = min(p.bins - 1, max(0, Int(point.y / side * CGFloat(p.bins))))
        let bj = min(p.bins - 1, max(0, Int(point.x / side * CGFloat(p.bins))))
        func range(_ b: Int) -> (Int, Int) {
            // elements k with k * bins / span == b
            let lo = (b * p.span + p.bins - 1) / p.bins
            let hi = ((b + 1) * p.span + p.bins - 1) / p.bins - 1
            return (p.first + lo + 1, p.first + hi + 1)
        }
        let (r0, r1) = range(bi), (c0, c1) = range(bj)
        let d = Double(density[bi * p.bins + bj]) / 255
        let rows = r0 == r1 ? "row \(r0)" : "rows \(r0)–\(r1)"
        let cols = c0 == c1 ? "column \(c0)" : "columns \(c0)–\(c1)"
        let share = d == 0 ? "zero" : (r0 == r1 && c0 == c1 ? "non-zero" : String(format: "%.0f %% non-zero", max(1, 100 * d)))
        return "\(rows), \(cols): \(share)"
    }

    /// The colour map of log|H|: dark blue … teal … yellow (viridis-like).
    static func heat(_ t: Double) -> CGColor {
        let stops: [(Double, Double, Double)] = [(0.267, 0.005, 0.329), (0.231, 0.322, 0.545), (0.129, 0.569, 0.549),
                                                 (0.369, 0.788, 0.384), (0.993, 0.906, 0.144)]
        let x = max(0, min(1, t)) * Double(stops.count - 1)
        let k = min(stops.count - 2, Int(x))
        let f = x - Double(k)
        let a = stops[k], b = stops[k + 1]
        return CGColor(srgbRed: a.0 + (b.0 - a.0) * f, green: a.1 + (b.1 - a.1) * f, blue: a.2 + (b.2 - a.2) * f, alpha: 1)
    }

    /// The picture as an image: one pixel per bin.
    private static func image(_ p: NuclearMatrixPicture, colouring: MatrixColouring, dark: Bool) -> CGImage? {
        guard p.bins > 0, let density = Data(base64Encoded: p.density), let magnitude = Data(base64Encoded: p.magnitude),
              density.count == p.bins * p.bins, magnitude.count == p.bins * p.bins else { return nil }
        let n = p.bins
        var rgba = [UInt8](repeating: 0, count: n * n * 4)
        let background: (Double, Double, Double) = dark ? (0.11, 0.11, 0.12) : (1, 1, 1)
        // the accent of the pattern (blue), mixed with the background by the share of non-zero elements
        let ink: (Double, Double, Double) = dark ? (0.39, 0.66, 1.0) : (0.0, 0.32, 0.85)
        var palette = [(UInt8, UInt8, UInt8)](repeating: (0, 0, 0), count: 256)
        for v in 0..<256 {
            let c: (Double, Double, Double)
            if v == 0 {
                c = background
            } else if colouring == .pattern {
                let w = 0.35 + 0.65 * Double(v) / 255   // even a single element stays visible
                c = (background.0 + (ink.0 - background.0) * w, background.1 + (ink.1 - background.1) * w,
                     background.2 + (ink.2 - background.2) * w)
            } else {
                let col = heat(Double(v - 1) / 254).components ?? [0, 0, 0, 1]
                c = (Double(col[0]), Double(col[1]), Double(col[2]))
            }
            palette[v] = (UInt8(max(0, min(255, c.0 * 255))), UInt8(max(0, min(255, c.1 * 255))), UInt8(max(0, min(255, c.2 * 255))))
        }
        let source = colouring == .pattern ? density : magnitude
        source.withUnsafeBytes { (buf: UnsafeRawBufferPointer) in
            for q in 0..<(n * n) {
                let (r, g, b) = palette[Int(buf[q])]
                rgba[4 * q] = r
                rgba[4 * q + 1] = g
                rgba[4 * q + 2] = b
                rgba[4 * q + 3] = 255
            }
        }
        guard let provider = CGDataProvider(data: Data(rgba) as CFData),
              let space = CGColorSpace(name: CGColorSpace.sRGB) else { return nil }
        return CGImage(width: n, height: n, bitsPerComponent: 8, bitsPerPixel: 32, bytesPerRow: 4 * n, space: space,
                       bitmapInfo: CGBitmapInfo(rawValue: CGImageAlphaInfo.noneSkipLast.rawValue), provider: provider,
                       decode: nil, shouldInterpolate: false, intent: .defaultIntent)
    }
}
