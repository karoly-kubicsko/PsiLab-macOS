//
//  ScanViews.swift
//  MacQC
//
//  ΨLab 7.0: the Scan tab. Its Relaxed scan page shows the result of the job "Relaxed scan"
//  (Gaussian's Opt=ModRedundant with S): the energy along the scanned coordinate, the molecule at
//  every point (played point by point; 12.0: in the 3D viewer), the table of points with their constrained optimisations,
//  the minima and maxima, and a button that starts a transition-state search from a point. Its
//  Potential curve page is the curve of single-point energies along the bond length (versions
//  1–6, PotentialCurveView.swift).
//

import SwiftUI
import Charts

struct ScanTabView: View {
    @Bindable var model: AppModel

    var body: some View {
        VStack(spacing: 0) {
            Picker("Page", selection: $model.scanPage) {
                ForEach(ScanPage.allCases) { page in
                    Text(page.rawValue).tag(page)
                }
            }
            .pickerStyle(.segmented)
            .labelsHidden()
            .frame(maxWidth: 360)
            .padding(.top, 10)
            .padding(.bottom, 4)
            .accessibilityLabel("Scan page")

            switch model.scanPage {
            case .relaxed:
                RelaxedScanView(model: model)
            case .potentialCurve:
                PotentialCurveView(model: model)
            }
        }
    }
}

struct RelaxedScanView: View {
    @Bindable var model: AppModel
    @State private var selectedPoint: Int?
    /// Play points (initially off; on with the launch argument -psilabPlaySteps YES, as on the
    /// Optimisation tab).
    @State private var playing = UserDefaults.standard.bool(forKey: "psilabPlaySteps")

    var body: some View {
        if let out = model.result, let scan = out.scan, !scan.points.isEmpty {
            ScrollView {
                VStack(alignment: .leading, spacing: 18) {
                    if let message = model.errorMessage {
                        ErrorBanner(message: message)
                    }
                    header(out, scan)
                    moleculePanel(scan)
                    chart(scan)
                    SectionTitle("Points")
                    pointsTable(scan)
                    if !scan.extrema.isEmpty {
                        SectionTitle("Minima and maxima")
                        extremaPanel(scan)
                    }
                }
                .padding(20)
            }
            .task(id: playing) {
                await playPoints()
            }
        } else {
            JobPlaceholder(model: model, title: "No relaxed scan", systemImage: "chart.xyaxis.line",
                           text: "Choose the job Relaxed scan in the sidebar, pick the bond, angle or dihedral angle to scan, its start, the number of steps and the step size, and press Run. Each point is an optimisation with the scanned coordinate frozen.")
        }
    }

    // MARK: Selection and playing

    private func shownIndex(_ scan: MacQCScan) -> Int {
        if let s = selectedPoint, let i = scan.points.firstIndex(where: { $0.index == s }) { return i }
        if scan.highest >= 0, let i = scan.points.firstIndex(where: { $0.index == scan.highest }) { return i }
        return 0
    }

    private func select(_ position: Int, in scan: MacQCScan) {
        guard !scan.points.isEmpty else { return }
        playing = false
        selectedPoint = scan.points[min(max(position, 0), scan.points.count - 1)].index
    }

    /// While Play is on: the next point every 0.7 s, the last one held for 1.6 s, then from the start.
    private func playPoints() async {
        guard playing else { return }
        while !Task.isCancelled {
            guard let scan = model.result?.scan, scan.points.count > 1 else {
                playing = false
                return
            }
            let index = shownIndex(scan)
            let last = scan.points.count - 1
            try? await Task.sleep(for: .milliseconds(index == last ? 1600 : 700))
            if Task.isCancelled || !playing { return }
            selectedPoint = scan.points[index == last ? 0 : index + 1].index
        }
    }

    private func kcal(_ hartree: Double) -> Double { hartree * Units.kcalPerMolPerHartree }

    // MARK: Header

    private func header(_ out: MacQCOutput, _ scan: MacQCScan) -> some View {
        VStack(alignment: .leading, spacing: 4) {
            Text("\(out.moleculeTitle) · \(MethodChoice.title(of: out.input.method))/\(out.input.basis) · "
                 + (scan.relaxed ? "relaxed scan of " : "scan of ") + scan.coordinate.definition)
                .font(.title2.weight(.semibold))
            Text(String(format: "%@ from %@ in %ld steps of %@ (%ld points); ", scan.coordinate.type.capitalized,
                        Fmt.scanValue(scan.start, angular: scan.isAngular), scan.steps,
                        Fmt.scanValue(scan.stepSize, angular: scan.isAngular), scan.steps + 1)
                 + (scan.relaxed
                    ? "at each point the scanned coordinate is frozen and all other coordinates are optimised."
                    : "a diatomic has nothing else to relax, so the points are single-point energies."))
                .foregroundStyle(.secondary)
            Label(scan.message, systemImage: scan.converged ? "checkmark.seal.fill" : "exclamationmark.triangle.fill")
                .foregroundStyle(scan.converged ? Color.green : Color.orange)
        }
    }

    // MARK: Molecule

    private func moleculePanel(_ scan: MacQCScan) -> some View {
        let points = scan.points
        let index = shownIndex(scan)
        let point = points[index]
        let atoms = point.geometry.atoms
        let aligned = alignedSequence(points.map { atomCoordinates($0.geometry.atoms) }, atoms: atoms)
        let lowest = points.first(where: { $0.index == scan.lowest })?.energy ?? point.energy
        return VStack(alignment: .leading, spacing: 10) {
            // ΨLab 12.0: the 3D viewer; 15.0: the full width of the panel, the point's data below
            // the buttons
            Molecule3DView(atoms: atomsMoved(atoms, to: aligned[index]), framing: aligned,
                           accessibilityText: "Molecule at this scan point")
                .viewer3DFrame(Viewer3DSize.steps, minWidth: 480)
            HStack(spacing: 8) {
                Button {
                    select(0, in: scan)
                } label: {
                    Image(systemName: "backward.end.fill")
                }
                .help("First point")
                .accessibilityLabel("First point")
                Button {
                    select(index - 1, in: scan)
                } label: {
                    Image(systemName: "backward.frame.fill")
                }
                .help("Previous point")
                .accessibilityLabel("Previous point")
                .disabled(index == 0)
                Button {
                    playing.toggle()
                } label: {
                    Label(playing ? "Pause" : "Play points", systemImage: playing ? "pause.fill" : "play.fill")
                }
                .help("Show the points one after the other")
                .disabled(points.count < 2)
                Button {
                    select(index + 1, in: scan)
                } label: {
                    Image(systemName: "forward.frame.fill")
                }
                .help("Next point")
                .accessibilityLabel("Next point")
                .disabled(index == points.count - 1)
                Button {
                    select(points.count - 1, in: scan)
                } label: {
                    Image(systemName: "forward.end.fill")
                }
                .help("Last point")
                .accessibilityLabel("Last point")
                Slider(value: Binding(get: { Double(index) },
                                      set: { select(Int($0.rounded()), in: scan) }),
                       in: 0...Double(max(points.count - 1, 1)), step: 1)
                    .disabled(points.count < 2)
                    .accessibilityLabel("Scan point")
            }
            VStack(alignment: .leading, spacing: 6) {
                Text("Point \(index + 1) of \(points.count)")
                    .font(.headline)
                HStack(alignment: .top, spacing: 36) {
                    Grid(alignment: .leading, horizontalSpacing: 12, verticalSpacing: 4) {
                        GridRow {
                            Text(scan.coordinate.definition).foregroundStyle(.secondary)
                            Text(Fmt.scanValue(point.measured, angular: scan.isAngular))
                        }
                        GridRow {
                            Text("Energy").foregroundStyle(.secondary)
                            Text(Fmt.energy10(point.energy) + " E_h")
                        }
                        GridRow {
                            Text("Above the lowest").foregroundStyle(.secondary)
                            Text(String(format: "%.4f kcal/mol", kcal(point.energy - lowest)))
                        }
                        if scan.relaxed {
                            GridRow {
                                Text("Optimisation").foregroundStyle(.secondary)
                                Text("\(point.steps) step\(point.steps == 1 ? "" : "s")" + (point.converged ? "" : ", not converged"))
                                    .foregroundStyle(point.converged ? Color.primary : Color.orange)
                            }
                        }
                    }
                    Grid(alignment: .leading, horizontalSpacing: 12, verticalSpacing: 4) {
                        ForEach(Array(Fmt.geometryRows(point.geometry).enumerated()), id: \.offset) { _, row in
                            GridRow {
                                Text(row.0).foregroundStyle(.secondary)
                                Text(row.1)
                            }
                        }
                    }
                    Spacer(minLength: 0)
                }
                .monospacedDigit()
                .textSelection(.enabled)
                if scan.relaxed && atoms.count > 2 {
                    Button {
                        model.startTransitionStateSearch(from: point)
                    } label: {
                        Label("Start a TS search from this point", systemImage: "arrow.up.forward.circle")
                    }
                    .secondaryActionStyle()
                    .disabled(!model.scanMatchesMolecule || model.isRunning)
                    .help("Copies this geometry into the inputs and selects the job TS optimisation + frequencies (Gaussian's route: scan, then Opt=TS from the highest point).")
                }
                // 17.0: this point (or every point) as .xyz or Gaussian input
                if let out = model.result {
                    SaveStructureButtons(model: model, structure: pointStructure(point, of: scan),
                                         noun: "Point", method: out.input.method, basis: out.input.basis,
                                         allTag: "scan-points",
                                         all: { scan.points.compactMap { pointStructure($0, of: scan) } })
                }
            }
            Text("The optimised molecule at each scan point in 3D, the points turned onto each other (bond lengths in Å). Drag to turn the molecule, pinch or ⌥-scroll to zoom; click two, three or four atoms to measure a distance, an angle or a dihedral angle. The highest point is shown first. Choose a point with the buttons, the slider, the chart or the table, or press Play points.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }

    /// 17.0: a scan point as a structure to save: its optimised geometry, the scanned coordinate
    /// (ASCII, for the comment lines of .xyz and Gaussian files), the level and the energy.
    private func pointStructure(_ point: MacQCScanPoint, of scan: MacQCScan) -> SavedMolecule? {
        let position = (scan.points.firstIndex(where: { $0.index == point.index }) ?? 0) + 1
        let value = scan.isAngular ? String(format: "%.2f deg", point.measured) : String(format: "%.4f Angstrom", point.measured)
        let coordinate = scan.coordinate.definition.replacingOccurrences(of: "\u{2013}", with: "-")
        let label = (scan.relaxed ? "relaxed scan point " : "scan point ") + "\(position) of \(scan.points.count), "
            + coordinate + " = " + value
        return model.resultStructure(atoms: point.geometry.atoms, label: label,
                                     fileTag: String(format: "point%02ld", position), energy: point.energy)
    }

    // MARK: Chart

    private func chart(_ scan: MacQCScan) -> some View {
        let xLabel = scan.coordinate.definition + (scan.isAngular ? " (°)" : " (Å)")
        let lowest = scan.points.first(where: { $0.index == scan.lowest })?.energy ?? (scan.points.map { $0.energy }.min() ?? 0)
        let shown = scan.points[shownIndex(scan)]
        return VStack(alignment: .leading, spacing: 8) {
            Text("Energy relative to the lowest point (kcal/mol)")
                .font(.subheadline.weight(.semibold))
            Chart {
                ForEach(scan.points) { p in
                    LineMark(x: .value(xLabel, p.value), y: .value("ΔE (kcal/mol)", kcal(p.energy - lowest)))
                        .interpolationMethod(.catmullRom)
                    PointMark(x: .value(xLabel, p.value), y: .value("ΔE (kcal/mol)", kcal(p.energy - lowest)))
                        .foregroundStyle(p.converged ? Color.accentColor : Color.orange)
                        .symbolSize(p.index == shown.index ? 110 : 40)
                }
                ForEach(scan.extrema) { x in
                    PointMark(x: .value(xLabel, x.value), y: .value("ΔE (kcal/mol)", kcal(x.energy - lowest)))
                        .symbol(x.kind == "maximum" ? BasicChartSymbolShape.triangle : BasicChartSymbolShape.diamond)
                        .foregroundStyle(x.kind == "maximum" ? Color.red : Color.green)
                        .symbolSize(90)
                        .annotation(position: x.kind == "maximum" ? .top : .bottom) {
                            Text(x.kind == "maximum" ? "max" : "min")
                                .font(.caption2)
                                .foregroundStyle(x.kind == "maximum" ? Color.red : Color.green)
                        }
                }
                RuleMark(x: .value("Shown point", shown.value))
                    .foregroundStyle(Color.accentColor.opacity(0.6))
                    .lineStyle(StrokeStyle(lineWidth: 1, dash: [3, 3]))
            }
            .chartXAxisLabel(xLabel)
            .chartYAxisLabel("ΔE (kcal/mol)")
            .chartOverlay { proxy in
                GeometryReader { geometry in
                    Rectangle()
                        .fill(Color.clear)
                        .contentShape(Rectangle())
                        .onTapGesture { location in
                            guard let plotFrame = proxy.plotFrame else { return }
                            let origin = geometry[plotFrame].origin
                            guard let x = proxy.value(atX: location.x - origin.x, as: Double.self) else { return }
                            if let nearest = scan.points.enumerated().min(by: { abs($0.element.value - x) < abs($1.element.value - x) }) {
                                select(nearest.offset, in: scan)
                            }
                        }
                }
            }
            .frame(minHeight: 300)
            Text("Dots: the scan points (orange: optimisation not converged); triangle and diamond: maximum and minimum of the parabola through a point and its neighbours. Click the chart to show the nearest point.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }

    // MARK: Table

    private func pointsTable(_ scan: MacQCScan) -> some View {
        let lowest = scan.points.first(where: { $0.index == scan.lowest })?.energy ?? 0
        let shown = scan.points[shownIndex(scan)].index
        return VStack(alignment: .leading, spacing: 6) {
            Grid(alignment: .trailing, horizontalSpacing: 18, verticalSpacing: 4) {
                GridRow {
                    Text("Point").gridColumnAlignment(.leading)
                    Text("Target")
                    Text("Measured")
                    Text("Energy (E_h)")
                    Text("ΔE (kcal/mol)")
                    if scan.relaxed {
                        Text("Steps")
                        Text("Max force")
                    }
                    Text("")
                }
                .font(.subheadline.weight(.semibold))
                .foregroundStyle(.secondary)
                Divider().gridCellUnsizedAxes(.horizontal)
                ForEach(Array(scan.points.enumerated()), id: \.element.index) { position, p in
                    GridRow {
                        Button {
                            select(position, in: scan)
                        } label: {
                            Text("\(position + 1)")
                                .frame(minWidth: 22, alignment: .leading)
                        }
                        .buttonStyle(.link)
                        .accessibilityLabel("Scan point \(position + 1)")
                        Text(Fmt.scanValue(p.value, angular: scan.isAngular))
                        Text(Fmt.scanValue(p.measured, angular: scan.isAngular))
                        Text(Fmt.energy10(p.energy))
                        Text(String(format: "%.4f", kcal(p.energy - lowest)))
                        if scan.relaxed {
                            Text("\(p.steps)")
                            Text(p.maxForce.map { String(format: "%.2e", $0) } ?? "–")
                        }
                        Image(systemName: p.converged ? "checkmark.circle.fill" : "exclamationmark.triangle.fill")
                            .foregroundStyle(p.converged ? Color.green : Color.orange)
                    }
                    .fontWeight(p.index == shown ? .semibold : .regular)
                }
            }
            .monospacedDigit()
            .textSelection(.enabled)
            Text(String(format: "Target: the frozen value; measured: the value at the optimised geometry. %.1f s in all.", scan.seconds))
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }

    private func extremaPanel(_ scan: MacQCScan) -> some View {
        let lowest = scan.points.first(where: { $0.index == scan.lowest })?.energy ?? 0
        return VStack(alignment: .leading, spacing: 8) {
            ForEach(scan.extrema) { x in
                InfoRow("\(x.kind.capitalized) near point \(x.index + 1)",
                        Fmt.scanValue(x.value, angular: scan.isAngular) + " · " + Fmt.energy10(x.energy) + " E_h"
                        + String(format: " · %.4f kcal/mol", kcal(x.energy - lowest)))
            }
            Text("Vertex of the parabola through an interior point lower (minimum) or higher (maximum) than both neighbours. A maximum of a relaxed scan approximates a transition state: refine it with a TS optimisation (button above).")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .dataPanel()
    }
}
