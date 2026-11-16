//
//  OrbitalViews.swift
//  MacQC
//
//  ΨLab 8.0: the Orbitals page of the Methods tab. The molecular orbitals of the result (RHF
//  canonical orbitals: all occupied ones and up to ten virtual ones; Full CI natural orbitals;
//  the two Szabo–Ostlund orbitals of HeH⁺) as an energy-level diagram and a table, and the
//  chosen orbital either as isosurfaces ψ = ±c (blue: positive, red: negative; SceneKit, rotate
//  with the mouse) or as a contour plot in a plane through the molecule. The backend computes
//  the orbital on a grid and the surfaces (marching tetrahedra, vertices refined onto the exact
//  isovalue), so the drawings show the orbitals of the calculation exactly.
//

import SwiftUI
import SceneKit
import AppKit

struct OrbitalsView: View {
    @Bindable var model: AppModel

    var body: some View {
        if let out = model.result, let orbitals = out.orbitals, !orbitals.sets.isEmpty {
            OrbitalExplorer(model: model, out: out, orbitals: orbitals)
                .id(model.resultSerial)
        } else {
            VStack(spacing: 14) {
                if let message = model.errorMessage {
                    ErrorBanner(message: message)
                        .frame(maxWidth: 560)
                }
                ContentUnavailableView(
                    "No orbitals",
                    systemImage: "circle.hexagongrid",
                    description: Text("Run a calculation: the orbitals of its Hartree–Fock step (and the natural orbitals of Full CI) appear here, as isosurfaces or as a contour plot.")
                )
            }
            .padding()
        }
    }
}

/// How an orbital is shown.
enum OrbitalDisplay: String, CaseIterable, Identifiable {
    case surface = "Isosurface"
    case contour = "Contour plot"
    var id: String { rawValue }
}

/// A plane through the centre of the molecule, spanned by two of its principal axes.
enum OrbitalPlaneChoice: Int, CaseIterable, Identifiable {
    case axes12 = 0, axes13, axes23
    var id: Int { rawValue }
}

struct OrbitalExplorer: View {
    var model: AppModel
    var out: MacQCOutput
    var orbitals: MacQCOrbitals

    @State private var setIndex: Int
    @State private var selected: Int
    @State private var display: OrbitalDisplay
    @State private var logIsovalue: Double
    @State private var plane: OrbitalPlaneChoice
    @State private var showAtoms = true

    @State private var surfaces: (positive: OrbitalMesh, negative: OrbitalMesh)?
    @State private var surfaceSerial = 0
    @State private var contour: OrbitalContour?
    @State private var computing = false
    @State private var problem: String?
    @State private var seconds: Double?

    init(model: AppModel, out: MacQCOutput, orbitals: MacQCOrbitals) {
        self.model = model
        self.out = out
        self.orbitals = orbitals
        // Launch arguments (UI smoke test): -psilabOrbitalSet natural -psilabOrbital LUMO
        // -psilabOrbitalView contour -psilabIsovalue 0.03
        let defaults = UserDefaults.standard
        var chosenSet = 0
        if defaults.string(forKey: "psilabOrbitalSet")?.lowercased() == "natural",
           let k = orbitals.sets.firstIndex(where: { $0.isNatural }) {
            chosenSet = k
        }
        let list = orbitals.sets[chosenSet].orbitals
        var choice = list.firstIndex(where: { $0.label == "HOMO" }) ?? 0
        if let name = defaults.string(forKey: "psilabOrbital") {
            if let k = list.firstIndex(where: { $0.label.lowercased() == name.lowercased() }) {
                choice = k
            } else if let n = Int(name), let k = list.firstIndex(where: { $0.number == n }) {
                choice = k
            }
        }
        _setIndex = State(initialValue: chosenSet)
        _selected = State(initialValue: choice)
        _display = State(initialValue: defaults.string(forKey: "psilabOrbitalView")?.lowercased() == "contour" ? .contour : .surface)
        // planar molecules: the molecular plane; others: a plane that contains their unique axis
        _plane = State(initialValue: OrbitalGeometry(atoms: orbitals.atoms).isPlanar ? .axes12 : .axes13)
        let iso = defaults.double(forKey: "psilabIsovalue")
        _logIsovalue = State(initialValue: log10(iso > 0 ? min(max(iso, 0.002), 0.5) : 0.05))
    }

    private var orbitalSet: MacQCOrbitalSet { orbitals.sets[min(setIndex, orbitals.sets.count - 1)] }
    private var isovalue: Double { pow(10, logIsovalue) }
    private var orbital: MacQCOrbitalInfo? {
        selected < orbitalSet.orbitals.count ? orbitalSet.orbitals[selected] : nil
    }

    /// Changes when the drawing must be recomputed.
    private var drawingKey: String {
        "\(setIndex)|\(selected)|\(display.rawValue)|\(String(format: "%.4f", logIsovalue))|\(plane.rawValue)"
    }

    var body: some View {
        ScrollView {
            VStack(alignment: .leading, spacing: 16) {
                header
                HStack(alignment: .top, spacing: 16) {
                    VStack(alignment: .leading, spacing: 12) {
                        if !orbitalSet.isNatural {
                            OrbitalLevelDiagram(orbitalSet: orbitalSet, selected: $selected)
                                .frame(height: 250)
                                .dataPanel()
                        }
                        orbitalTable
                    }
                    .frame(width: 330)
                    VStack(alignment: .leading, spacing: 10) {
                        controls
                        drawing
                            .frame(minHeight: 420, maxHeight: 520)
                        caption
                    }
                    .frame(maxWidth: .infinity)
                }
            }
            .padding(20)
        }
        .task(id: drawingKey) {
            await compute()
        }
    }

    // MARK: Header and controls

    private var header: some View {
        VStack(alignment: .leading, spacing: 6) {
            Text("\(out.moleculeTitle) · \(orbitals.basis) · molecular orbitals")
                .font(.title2.weight(.semibold))
            if orbitals.sets.count > 1 {
                Picker("Orbitals", selection: Binding(get: { setIndex }, set: { k in
                    setIndex = k
                    let list = orbitals.sets[k].orbitals
                    selected = list.firstIndex(where: { $0.label == "HOMO" }) ?? 0
                })) {
                    ForEach(Array(orbitals.sets.enumerated()), id: \.offset) { k, s in
                        Text(s.title).tag(k)
                    }
                }
                .pickerStyle(.segmented)
                .frame(maxWidth: 460)
                .accessibilityLabel("Orbital set")
            } else {
                Text(orbitalSet.title)
                    .foregroundStyle(.secondary)
            }
        }
    }

    /// Two rows, so that they fit at the minimum window width.
    private var controls: some View {
        VStack(alignment: .leading, spacing: 8) {
            HStack(spacing: 12) {
                Picker("Display", selection: $display) {
                    ForEach(OrbitalDisplay.allCases) { d in
                        Text(d.rawValue).tag(d)
                    }
                }
                .pickerStyle(.segmented)
                .labelsHidden()
                .frame(width: 230)
                .accessibilityLabel("Orbital display")
                if display == .surface {
                    Toggle("Atoms", isOn: $showAtoms)
                        .toggleStyle(.checkbox)
                }
                if computing {
                    ProgressView()
                        .controlSize(.small)
                }
                Spacer(minLength: 0)
            }
            HStack(spacing: 10) {
                Text("Isovalue")
                    .foregroundStyle(.secondary)
                Slider(value: $logIsovalue, in: log10(0.002)...log10(0.5))
                    .frame(minWidth: 90, maxWidth: 140)
                    .accessibilityLabel("Isovalue")
                Text(String(format: "%.3f", isovalue))
                    .monospacedDigit()
                    .frame(width: 44, alignment: .leading)
                if display == .contour {
                    Picker("Plane", selection: $plane) {
                        ForEach(OrbitalPlaneChoice.allCases) { p in
                            Text(planeTitle(p)).tag(p)
                        }
                    }
                    .labelsHidden()
                    .frame(maxWidth: 250)
                    .accessibilityLabel("Contour plane")
                }
                Spacer(minLength: 0)
            }
        }
    }

    @ViewBuilder
    private var drawing: some View {
        if let problem {
            ErrorBanner(message: problem)
        } else if display == .surface {
            OrbitalSceneView(atoms: orbitals.atoms, surfaces: surfaces, serial: surfaceSerial, showAtoms: showAtoms)
                .background(Color(nsColor: .textBackgroundColor), in: RoundedRectangle(cornerRadius: 12))
                .clipShape(RoundedRectangle(cornerRadius: 12))
                .accessibilityLabel("Orbital isosurface")
        } else if let contour {
            OrbitalContourView(contour: contour, isovalue: isovalue)
                .background(Color(nsColor: .textBackgroundColor), in: RoundedRectangle(cornerRadius: 12))
                .accessibilityLabel("Orbital contour plot")
        } else {
            Color.clear
        }
    }

    private var caption: some View {
        let info = orbital
        var parts: [String] = []
        if let info {
            var name = "Orbital \(info.number)"
            if !info.label.isEmpty { name += " (\(info.label))" }
            parts.append(name)
            if let e = info.energy {
                parts.append(String(format: "ε = %.6f E_h = %.4f eV", e, e * 27.211386245981))
            }
            parts.append(orbitalSet.isNatural ? String(format: "occupation %.6f", info.occupation)
                                       : (info.occupation > 0 ? "occupied (2 electrons)" : "virtual"))
        }
        if display == .surface, let s = surfaces {
            parts.append("\(s.positive.triangleCount + s.negative.triangleCount) triangles")
        }
        if let t = seconds {
            parts.append(String(format: "%.2f s", t))
        }
        let detail = display == .surface
            ? "Blue: ψ = +\(String(format: "%.3f", isovalue)), red: ψ = −\(String(format: "%.3f", isovalue)) (a.u., bohr^−3/2). Drag to rotate, pinch or scroll to zoom. The overall sign of an orbital is arbitrary."
            : "Contours at ±\(String(format: "%.3f", isovalue)) × 2^k (blue solid: positive, red dashed: negative; the thick ones at the isovalue), the nodal line dotted; atoms projected onto the plane (paler: farther from it)."
        return VStack(alignment: .leading, spacing: 4) {
            Text(parts.joined(separator: " · "))
                .font(.callout.monospacedDigit())
                .textSelection(.enabled)
            Text(detail)
                .font(.caption)
                .foregroundStyle(.secondary)
        }
    }

    /// "Molecular plane", "Plane of N1, H2 (axes 1–3)" or "Principal axes 2–3": the nuclei within
    /// 0.05 bohr of the plane.
    private func planeTitle(_ p: OrbitalPlaneChoice) -> String {
        let g = OrbitalGeometry(atoms: orbitals.atoms)
        let axes: String
        switch p {
        case .axes12: axes = "axes 1–2"
        case .axes13: axes = "axes 1–3"
        case .axes23: axes = "axes 2–3"
        }
        let spec = g.planeSpec(p, points: 3)
        let inPlane = g.projectedAtoms(orbitals.atoms, onto: spec).enumerated()
            .filter { $0.element.distance < 0.05 }.map { "\($0.element.symbol)\($0.offset + 1)" }
        if inPlane.count == orbitals.atoms.count && orbitals.atoms.count > 2 { return "Molecular plane" }
        if inPlane.isEmpty { return "Principal " + axes }
        return "Plane of " + inPlane.joined(separator: ", ") + " (" + axes + ")"
    }

    // MARK: Table

    private var orbitalTable: some View {
        VStack(alignment: .leading, spacing: 4) {
            Grid(alignment: .trailing, horizontalSpacing: 12, verticalSpacing: 3) {
                GridRow {
                    Text("No.").gridColumnAlignment(.leading)
                    Text("").gridColumnAlignment(.leading)
                    Text(orbitalSet.isNatural ? "Occupation" : "ε (E_h)")
                    Text(orbitalSet.isNatural ? "" : "ε (eV)")
                    Text(orbitalSet.isNatural ? "" : "Occ.")
                }
                .font(.subheadline.weight(.semibold))
                .foregroundStyle(.secondary)
                Divider().gridCellUnsizedAxes(.horizontal)
                ForEach(Array(orbitalSet.orbitals.enumerated()), id: \.offset) { k, o in
                    GridRow {
                        Button {
                            selected = k
                        } label: {
                            Text("\(o.number)").frame(minWidth: 22, alignment: .leading)
                        }
                        .buttonStyle(.link)
                        .accessibilityLabel("Orbital \(o.number)")
                        Text(o.label)
                            .foregroundStyle(o.label == "HOMO" || o.label == "LUMO" ? Color.accentColor : Color.secondary)
                        if orbitalSet.isNatural {
                            Text(String(format: "%.6f", o.occupation))
                            Text("")
                            Text("")
                        } else {
                            Text(o.energy.map { String(format: "%.6f", $0) } ?? "–")
                            Text(o.energy.map { String(format: "%.3f", $0 * 27.211386245981) } ?? "–")
                            Text(String(format: "%.0f", o.occupation))
                        }
                    }
                    .fontWeight(k == selected ? .semibold : .regular)
                    .foregroundStyle(k == selected ? Color.primary : Color.primary.opacity(0.85))
                }
            }
            .monospacedDigit()
            .textSelection(.enabled)
            Text(orbitalSet.isNatural
                 ? "Natural orbitals: eigenvectors of the Full CI one-particle density matrix, by occupation."
                 : "All occupied orbitals and up to ten virtual ones. Click a number to draw the orbital.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }

    // MARK: Computation

    private func compute() async {
        guard selected < orbitalSet.coefficients.count else { return }
        computing = true
        problem = nil
        let start = Date()
        let orbitals = self.orbitals
        let setIndex = self.setIndex, selected = self.selected, iso = isovalue
        do {
            if display == .surface {
                let result = try await Task.detached(priority: .userInitiated) {
                    let json = try MacQCEngine.orbitalsJSON(orbitals)
                    return try MacQCEngine.orbitalSurfaces(json, set: setIndex, orbital: selected, isovalue: iso)
                }.value
                if Task.isCancelled { return }
                surfaces = result
                surfaceSerial += 1
            } else {
                let geometry = OrbitalGeometry(atoms: orbitals.atoms)
                let spec = geometry.planeSpec(plane, points: 141)
                let origin = spec.origin, u = spec.u, v = spec.v, n = spec.n
                let values = try await Task.detached(priority: .userInitiated) {
                    let json = try MacQCEngine.orbitalsJSON(orbitals)
                    return try MacQCEngine.orbitalPlane(json, set: setIndex, orbital: selected, origin: origin,
                                                        u: u, v: v, nu: n, nv: n)
                }.value
                if Task.isCancelled { return }
                contour = OrbitalContour(spec: spec, values: values, atoms: geometry.projectedAtoms(orbitals.atoms, onto: spec))
            }
            seconds = Date().timeIntervalSince(start)
        } catch {
            if !Task.isCancelled { problem = error.localizedDescription }
        }
        computing = false
    }
}

// MARK: - Energy levels

/// Orbital energies on a vertical axis (valence region: from 1.2 E_h below the HOMO to 0.8 E_h
/// above the LUMO; lower and higher orbitals are counted at the edges). Occupied levels carry
/// two electrons (↑↓); degenerate levels stand side by side. Click a level to draw it.
struct OrbitalLevelDiagram: View {
    var orbitalSet: MacQCOrbitalSet
    @Binding var selected: Int

    private struct Level {
        var indices: [Int]
        var energy: Double
    }

    private var window: (lo: Double, hi: Double) {
        let energies = orbitalSet.orbitals.compactMap { $0.energy }
        guard let homo = orbitalSet.orbitals.last(where: { $0.occupation > 0 })?.energy,
              let minE = energies.min(), let maxE = energies.max() else { return (-1, 1) }
        let lumo = orbitalSet.orbitals.first(where: { $0.occupation == 0 })?.energy ?? homo + 0.5
        return (max(minE, homo - 1.2) - 0.05, min(maxE, lumo + 0.8) + 0.05)
    }

    private var levels: [Level] {
        var out: [Level] = []
        for (k, o) in orbitalSet.orbitals.enumerated() {
            guard let e = o.energy else { continue }
            if let last = out.last, abs(last.energy - e) < 1e-4 {
                out[out.count - 1].indices.append(k)
            } else {
                out.append(Level(indices: [k], energy: e))
            }
        }
        return out
    }

    var body: some View {
        GeometryReader { geometry in
            let w = window
            let size = geometry.size
            let all = levels
            let shown = all.filter { $0.energy >= w.lo && $0.energy <= w.hi }
            let below = all.filter { $0.energy < w.lo }.reduce(0) { $0 + $1.indices.count }
            let above = all.filter { $0.energy > w.hi }.reduce(0) { $0 + $1.indices.count }
            ZStack(alignment: .topLeading) {
                Canvas { gc, sz in
                    let axisX = 46.0
                    var axis = Path()
                    axis.move(to: CGPoint(x: axisX, y: 8))
                    axis.addLine(to: CGPoint(x: axisX, y: Double(sz.height) - 8))
                    gc.stroke(axis, with: .color(.secondary.opacity(0.6)), lineWidth: 1)
                    for t in ticks(w) {
                        let y = yPosition(t, w, sz)
                        var tick = Path()
                        tick.move(to: CGPoint(x: axisX - 4, y: y))
                        tick.addLine(to: CGPoint(x: axisX, y: y))
                        gc.stroke(tick, with: .color(.secondary), lineWidth: 1)
                        gc.draw(Text(String(format: "%.1f", t)).font(.caption2.monospacedDigit()).foregroundStyle(.secondary),
                                at: CGPoint(x: axisX - 6, y: y), anchor: .trailing)
                    }
                    for level in shown {
                        let y = yPosition(level.energy, w, sz)
                        let n = level.indices.count
                        let width = min(46.0, (Double(sz.width) - axisX - 30) / Double(n) - 8)
                        for (slot, k) in level.indices.enumerated() {
                            let x0 = axisX + 18 + Double(slot) * (width + 8)
                            var line = Path()
                            line.move(to: CGPoint(x: x0, y: y))
                            line.addLine(to: CGPoint(x: x0 + width, y: y))
                            let isSelected = k == selected
                            let occupied = orbitalSet.orbitals[k].occupation > 0
                            let color: Color = isSelected ? Color.accentColor : (occupied ? Color.primary : Color.secondary)
                            gc.stroke(line, with: .color(color), lineWidth: isSelected ? 3.5 : 2)
                            if occupied {
                                gc.draw(Text("↑↓").font(.caption2), at: CGPoint(x: x0 + width / 2, y: y - 7), anchor: .center)
                            }
                        }
                        if let label = level.indices.map({ orbitalSet.orbitals[$0].label }).first(where: { $0 == "HOMO" || $0 == "LUMO" }) {
                            let x = axisX + 18 + Double(n) * (width + 8)
                            gc.draw(Text(label).font(.caption2.weight(.semibold)).foregroundStyle(Color.accentColor),
                                    at: CGPoint(x: x, y: y), anchor: .leading)
                        }
                    }
                    if below > 0 {
                        gc.draw(Text("\(below) lower").font(.caption2).foregroundStyle(.secondary),
                                at: CGPoint(x: axisX + 18, y: Double(sz.height) - 4), anchor: .bottomLeading)
                    }
                    if above > 0 {
                        gc.draw(Text("\(above) higher").font(.caption2).foregroundStyle(.secondary),
                                at: CGPoint(x: axisX + 18, y: 4), anchor: .topLeading)
                    }
                }
                .contentShape(Rectangle())
                .onTapGesture { location in
                    // the nearest level, then the slot within it
                    guard let level = shown.min(by: { abs(yPosition($0.energy, w, size) - Double(location.y)) <
                                                      abs(yPosition($1.energy, w, size) - Double(location.y)) }) else { return }
                    let n = level.indices.count
                    let width = min(46.0, (Double(size.width) - 46 - 30) / Double(n) - 8)
                    let slot = Int(((Double(location.x) - 46 - 18) / (width + 8)).rounded(.down))
                    selected = level.indices[min(max(slot, 0), n - 1)]
                }
                Text("E (E_h)")
                    .font(.caption2)
                    .foregroundStyle(.secondary)
                    .padding(.leading, 2)
            }
        }
        .accessibilityLabel("Orbital energy levels")
    }

    private func yPosition(_ e: Double, _ w: (lo: Double, hi: Double), _ size: CGSize) -> Double {
        let top = 16.0, bottom = Double(size.height) - 16
        return bottom - (e - w.lo) / max(w.hi - w.lo, 1e-6) * (bottom - top)
    }

    private func ticks(_ w: (lo: Double, hi: Double)) -> [Double] {
        let span = w.hi - w.lo
        let step = span > 2 ? 0.5 : (span > 1 ? 0.25 : 0.1)
        var t = (w.lo / step).rounded(.up) * step
        var out: [Double] = []
        while t <= w.hi + 1e-9 {
            out.append(t)
            t += step
        }
        return out
    }
}

// MARK: - Geometry of the drawings

/// Centre and principal axes of the nuclei (the planes of the contour plot, the camera).
struct OrbitalGeometry {
    var center: [Double] = [0, 0, 0]
    var axes: [[Double]] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]   // largest spread first
    var radius = 1.0                                         // largest distance of a nucleus from the centre
    var isPlanar = true

    init(atoms: [MacQCAtom]) {
        let positions = atomCoordinates(atoms)
        guard !positions.isEmpty else { return }
        let n = Double(positions.count)
        center = (0..<3).map { k in positions.map { $0[k] }.reduce(0, +) / n }
        let c = center
        var cov = [[Double]](repeating: [0, 0, 0], count: 3)
        for p in positions {
            for a in 0..<3 {
                for b in 0..<3 {
                    cov[a][b] += (p[a] - c[a]) * (p[b] - c[b])
                }
            }
        }
        let (values, vectors) = symmetricEigen(cov)
        let order = [0, 1, 2].sorted { values[$0] > values[$1] + 1e-12 || (abs(values[$0] - values[$1]) <= 1e-12 && $0 < $1) }
        var e = order.map { k -> [Double] in
            let v = vectors[k]
            var big = 0
            for i in 1..<3 where abs(v[i]) > abs(v[big]) + 1e-9 { big = i }
            return v[big] < 0 ? v.map { -$0 } : v
        }
        let normal = vectorCross(e[0], e[1])
        // A symmetric top (NH₃): the first axis is any direction in the plane of the two equal ones;
        // turn it towards the first nucleus off the unique axis, so that the planes of the contour
        // plot contain atoms (axes 1 and 3: the σᵥ plane through N1 and H2).
        if values[order[0]] - values[order[1]] < 1e-3 * max(values[order[0]], 1e-12) {
            for p in positions {
                let d = (0..<3).map { p[$0] - c[$0] }
                let along = vectorDot(d, normal)
                let inPlane = (0..<3).map { d[$0] - along * normal[$0] }
                let length = vectorNorm(inPlane)
                if length > 1e-6 {
                    e[0] = inPlane.map { $0 / length }
                    e[1] = vectorCross(normal, e[0])
                    break
                }
            }
        }
        e[2] = vectorCross(e[0], e[1])
        axes = e
        radius = positions.map { p in vectorDistance(p, c) }.max() ?? 1
        let off = positions.map { p in abs(vectorDot((0..<3).map { k in p[k] - c[k] }, e[2])) }.max() ?? 0
        isPlanar = off < 0.02
    }

    /// A square plane through the centre (n x n points), spanned by two principal axes, wide
    /// enough for the nuclei plus 4 bohr.
    func planeSpec(_ choice: OrbitalPlaneChoice, points n: Int) -> OrbitalPlaneSpec {
        let (a, b, normal): ([Double], [Double], [Double])
        switch choice {
        case .axes12: (a, b, normal) = (axes[0], axes[1], axes[2])
        case .axes13: (a, b, normal) = (axes[0], axes[2], axes[1])
        case .axes23: (a, b, normal) = (axes[1], axes[2], axes[0])
        }
        let half = radius + 4.0
        let step = 2 * half / Double(n - 1)
        let origin = (0..<3).map { center[$0] - half * a[$0] - half * b[$0] }
        return OrbitalPlaneSpec(origin: origin, u: a.map { $0 * step }, v: b.map { $0 * step }, n: n,
                                normal: normal, half: half, axisA: a, axisB: b, center: center)
    }
}

nonisolated struct OrbitalPlaneSpec: Sendable {
    var origin: [Double]
    var u: [Double]
    var v: [Double]
    var n: Int
    var normal: [Double]
    var half: Double          // bohr from the centre to an edge
    var axisA: [Double]
    var axisB: [Double]
    var center: [Double]
}

struct ProjectedAtom {
    var symbol: String
    var Z: Int
    var a: Double             // bohr along axisA from the centre
    var b: Double
    var distance: Double      // from the plane
}

struct OrbitalContour {
    var spec: OrbitalPlaneSpec
    var values: [Float]
    var atoms: [ProjectedAtom]
}

extension OrbitalGeometry {
    /// The nuclei projected onto a plane, with their distances from it.
    func projectedAtoms(_ atoms: [MacQCAtom], onto spec: OrbitalPlaneSpec) -> [ProjectedAtom] {
        atoms.map { atom in
            let d = (0..<3).map { k in (k < atom.position.count ? atom.position[k] : 0) - spec.center[k] }
            return ProjectedAtom(symbol: atom.symbol, Z: atom.Z, a: vectorDot(d, spec.axisA), b: vectorDot(d, spec.axisB),
                                 distance: abs(vectorDot(d, spec.normal)))
        }
    }
}

// MARK: - Contour plot

/// Contour lines of the orbital in a plane (marching squares, saddles resolved by the value at
/// the centre of the cell): ±c·2^k for k = −1 … 4, blue solid for positive values, red dashed
/// for negative ones, the lines at ±c thicker; the nodal line ψ = 0 dotted.
struct OrbitalContourView: View {
    var contour: OrbitalContour
    var isovalue: Double

    var body: some View {
        Canvas { gc, size in
            let spec = contour.spec
            let n = spec.n
            let side = min(Double(size.width), Double(size.height)) - 24
            let scale = side / (2 * spec.half)            // points per bohr
            let ox = (Double(size.width) - side) / 2, oy = (Double(size.height) - side) / 2
            func toScreen(_ i: Double, _ j: Double) -> CGPoint {
                CGPoint(x: ox + i / Double(n - 1) * side, y: oy + side - j / Double(n - 1) * side)
            }
            var frame = Path()
            frame.addRect(CGRect(x: ox, y: oy, width: side, height: side))
            gc.stroke(frame, with: .color(.secondary.opacity(0.4)), lineWidth: 1)
            var levels: [(Double, Bool)] = [(0, false)]
            for k in -1...4 {
                let c = isovalue * pow(2, Double(k))
                levels.append((c, k == 0))
                levels.append((-c, k == 0))
            }
            for (level, thick) in levels {
                let segments = marchingSquares(contour.values, n: n, level: Float(level))
                guard !segments.isEmpty else { continue }
                var path = Path()
                for s in segments {
                    path.move(to: toScreen(s.0, s.1))
                    path.addLine(to: toScreen(s.2, s.3))
                }
                let color: Color = level > 0 ? .blue : (level < 0 ? .red : .secondary)
                let dash: [CGFloat] = level > 0 ? [] : (level < 0 ? [5, 3] : [1.5, 3])
                gc.stroke(path, with: .color(color.opacity(level == 0 ? 0.8 : 1)),
                          style: StrokeStyle(lineWidth: thick ? 2.2 : 1, lineCap: .round, dash: dash))
            }
            // the nuclei, projected (paler when farther from the plane)
            for atom in contour.atoms {
                let i = (atom.a + spec.half) / (2 * spec.half) * Double(n - 1)
                let j = (atom.b + spec.half) / (2 * spec.half) * Double(n - 1)
                let p = toScreen(i, j)
                let r = AtomStyle.radius(atom.Z) * 0.55
                let alpha = max(0.25, 1 - atom.distance / 3)
                let rect = CGRect(x: Double(p.x) - r, y: Double(p.y) - r, width: 2 * r, height: 2 * r)
                gc.fill(Path(ellipseIn: rect), with: .color(AtomStyle.color(atom.Z).opacity(alpha)))
                gc.stroke(Path(ellipseIn: rect), with: .color(Color.primary.opacity(0.5 * alpha)), lineWidth: 1)
                gc.draw(Text(atom.symbol).font(.caption2.bold()).foregroundStyle(Color.primary.opacity(alpha)),
                        at: CGPoint(x: Double(p.x), y: Double(p.y) - r - 7), anchor: .center)
            }
            // scale bar: 1 Å
            let bar = 1.0 / Units.bohrToAngstrom * scale
            var path = Path()
            path.move(to: CGPoint(x: ox + 12, y: oy + side - 12))
            path.addLine(to: CGPoint(x: ox + 12 + bar, y: oy + side - 12))
            gc.stroke(path, with: .color(.primary), lineWidth: 2)
            gc.draw(Text("1 Å").font(.caption2), at: CGPoint(x: ox + 12 + bar / 2, y: oy + side - 22), anchor: .center)
        }
    }

    /// Segments (i0, j0, i1, j1) in grid coordinates where the values cross `level`.
    private func marchingSquares(_ f: [Float], n: Int, level: Float) -> [(Double, Double, Double, Double)] {
        var out: [(Double, Double, Double, Double)] = []
        func value(_ i: Int, _ j: Int) -> Float { f[j * n + i] }
        func cross(_ i0: Int, _ j0: Int, _ i1: Int, _ j1: Int) -> (Double, Double) {
            let a = value(i0, j0), b = value(i1, j1)
            let t = a == b ? 0.5 : Double((level - a) / (b - a))
            return (Double(i0) + t * Double(i1 - i0), Double(j0) + t * Double(j1 - j0))
        }
        for j in 0..<(n - 1) {
            for i in 0..<(n - 1) {
                // corners: 0 (i,j), 1 (i+1,j), 2 (i+1,j+1), 3 (i,j+1)
                let v0 = value(i, j), v1 = value(i + 1, j), v2 = value(i + 1, j + 1), v3 = value(i, j + 1)
                var code = 0
                if v0 >= level { code |= 1 }
                if v1 >= level { code |= 2 }
                if v2 >= level { code |= 4 }
                if v3 >= level { code |= 8 }
                if code == 0 || code == 15 { continue }
                let bottom = { cross(i, j, i + 1, j) }
                let right = { cross(i + 1, j, i + 1, j + 1) }
                let top = { cross(i, j + 1, i + 1, j + 1) }
                let left = { cross(i, j, i, j + 1) }
                func add(_ p: (Double, Double), _ q: (Double, Double)) { out.append((p.0, p.1, q.0, q.1)) }
                switch code {
                case 1, 14: add(left(), bottom())
                case 2, 13: add(bottom(), right())
                case 3, 12: add(left(), right())
                case 4, 11: add(right(), top())
                case 6, 9: add(bottom(), top())
                case 7, 8: add(left(), top())
                case 5, 10:
                    let centre = (v0 + v1 + v2 + v3) / 4
                    if (centre >= level) == (code == 5) {
                        add(left(), top())
                        add(bottom(), right())
                    } else {
                        add(left(), bottom())
                        add(right(), top())
                    }
                default: break
                }
            }
        }
        return out
    }
}

// MARK: - Isosurfaces (SceneKit)

/// The two lobes of an orbital and the molecule (balls and sticks, bohr), in a SceneKit view that
/// rotates with the mouse. The camera looks along the second principal axis from 30° above the
/// plane of the two longest ones, with the third axis (the normal of a planar molecule, the C₃
/// axis of NH₃) up: out-of-plane lobes appear above and below the molecule.
struct OrbitalSceneView: View {
    var atoms: [MacQCAtom]
    var surfaces: (positive: OrbitalMesh, negative: OrbitalMesh)?
    var serial: Int                 // changes with the surfaces
    var showAtoms: Bool

    @State private var scene = SCNScene()
    @State private var camera = SCNNode()
    @State private var prepared = false

    var body: some View {
        SceneView(scene: scene, pointOfView: camera, options: [.allowsCameraControl, .autoenablesDefaultLighting])
            .onAppear {
                prepare()
                update()
            }
            .onChange(of: serial) { update() }
            .onChange(of: showAtoms) { update() }
    }

    private func prepare() {
        guard !prepared else { return }
        prepared = true
        scene.background.contents = NSColor.textBackgroundColor
        let g = OrbitalGeometry(atoms: atoms)
        let cam = SCNCamera()
        cam.zNear = 0.1
        cam.zFar = 500
        cam.fieldOfView = 35
        camera.camera = cam
        let side = g.axes[1], normal = g.axes[2]
        let t = 30.0 * Double.pi / 180
        let back = (0..<3).map { cos(t) * side[$0] + sin(t) * normal[$0] }
        let up = (0..<3).map { -sin(t) * side[$0] + cos(t) * normal[$0] }
        let distance = (g.radius + 5.0) / tan(17.5 * Double.pi / 180)
        let eye = (0..<3).map { g.center[$0] + distance * back[$0] }
        camera.position = SCNVector3(x: CGFloat(eye[0]), y: CGFloat(eye[1]), z: CGFloat(eye[2]))
        camera.look(at: SCNVector3(x: CGFloat(g.center[0]), y: CGFloat(g.center[1]), z: CGFloat(g.center[2])),
                    up: SCNVector3(x: CGFloat(up[0]), y: CGFloat(up[1]), z: CGFloat(up[2])),
                    localFront: SCNVector3(x: 0, y: 0, z: -1))
        scene.rootNode.addChildNode(camera)
        let light = SCNNode()
        light.light = SCNLight()
        light.light?.type = .ambient
        light.light?.intensity = 350
        scene.rootNode.addChildNode(light)
    }

    private func update() {
        scene.rootNode.childNodes.filter { $0.name == "lobe" || $0.name == "molecule" }.forEach { $0.removeFromParentNode() }
        if showAtoms {
            scene.rootNode.addChildNode(moleculeNode())
        }
        if let s = surfaces {
            for (mesh, color) in [(s.positive, NSColor(calibratedRed: 0.18, green: 0.45, blue: 0.95, alpha: 1)),
                                  (s.negative, NSColor(calibratedRed: 0.93, green: 0.30, blue: 0.22, alpha: 1))]
            where mesh.triangleCount > 0 {
                let node = SCNNode(geometry: geometry(mesh, color: color))
                node.name = "lobe"
                scene.rootNode.addChildNode(node)
            }
        }
    }

    private func geometry(_ mesh: OrbitalMesh, color: NSColor) -> SCNGeometry {
        let stride = MemoryLayout<Float>.size * 3
        let positions = mesh.positions.withUnsafeBufferPointer { Data(buffer: $0) }
        let normals = mesh.normals.withUnsafeBufferPointer { Data(buffer: $0) }
        let indices = mesh.indices.withUnsafeBufferPointer { Data(buffer: $0) }
        let vertexSource = SCNGeometrySource(data: positions, semantic: .vertex, vectorCount: mesh.vertexCount,
                                             usesFloatComponents: true, componentsPerVector: 3,
                                             bytesPerComponent: MemoryLayout<Float>.size, dataOffset: 0, dataStride: stride)
        let normalSource = SCNGeometrySource(data: normals, semantic: .normal, vectorCount: mesh.vertexCount,
                                             usesFloatComponents: true, componentsPerVector: 3,
                                             bytesPerComponent: MemoryLayout<Float>.size, dataOffset: 0, dataStride: stride)
        let element = SCNGeometryElement(data: indices, primitiveType: .triangles, primitiveCount: mesh.triangleCount,
                                         bytesPerIndex: MemoryLayout<Int32>.size)
        let g = SCNGeometry(sources: [vertexSource, normalSource], elements: [element])
        let material = SCNMaterial()
        material.diffuse.contents = color
        material.specular.contents = NSColor(white: 1, alpha: 0.35)
        material.shininess = 0.4
        material.lightingModel = .blinn
        material.isDoubleSided = true
        g.materials = [material]
        return g
    }

    private func moleculeNode() -> SCNNode {
        let root = SCNNode()
        root.name = "molecule"
        let positions = atomCoordinates(atoms)
        for (i, atom) in atoms.enumerated() {
            let sphere = SCNSphere(radius: atom.Z <= 2 ? 0.32 : 0.48)
            sphere.segmentCount = 32
            let m = SCNMaterial()
            m.diffuse.contents = NSColor(AtomStyle.color(atom.Z))
            m.lightingModel = .blinn
            sphere.materials = [m]
            let node = SCNNode(geometry: sphere)
            node.position = SCNVector3(x: CGFloat(positions[i][0]), y: CGFloat(positions[i][1]), z: CGFloat(positions[i][2]))
            root.addChildNode(node)
        }
        for (i, j) in moleculeBonds(atoms) {
            let a = positions[i], b = positions[j]
            let d = (0..<3).map { b[$0] - a[$0] }
            let length = vectorNorm(d)
            guard length > 1e-6 else { continue }
            let cylinder = SCNCylinder(radius: 0.1, height: CGFloat(length))
            let m = SCNMaterial()
            m.diffuse.contents = NSColor(white: 0.6, alpha: 1)
            cylinder.materials = [m]
            let node = SCNNode(geometry: cylinder)
            node.position = SCNVector3(x: CGFloat((a[0] + b[0]) / 2), y: CGFloat((a[1] + b[1]) / 2), z: CGFloat((a[2] + b[2]) / 2))
            // rotate the cylinder's y axis onto the bond
            let y = [0.0, 1.0, 0.0]
            let unit = d.map { $0 / length }
            let axis = vectorCross(y, unit)
            let sine = vectorNorm(axis)
            let angle = atan2(sine, vectorDot(y, unit))
            if sine > 1e-9 {
                node.rotation = SCNVector4(x: CGFloat(axis[0] / sine), y: CGFloat(axis[1] / sine), z: CGFloat(axis[2] / sine),
                                           w: CGFloat(angle))
            } else if unit[1] < 0 {
                node.rotation = SCNVector4(x: 1, y: 0, z: 0, w: CGFloat(Double.pi))
            }
            root.addChildNode(node)
        }
        return root
    }
}
