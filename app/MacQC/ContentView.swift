//
//  ContentView.swift
//  MacQC
//
//  NavigationSplitView: input sidebar + tabbed detail area.
//

import SwiftUI
import AppKit

struct ContentView: View {
    @Bindable var model: AppModel
    @AppStorage(AppearanceChoice.storageKey) private var appearanceRaw: String = AppearanceChoice.system.rawValue

    var body: some View {
        NavigationSplitView {
            // 20.0: the kind of calculation at the top of the sidebar; each kind has its own
            // sidebar and tabs (the electronic-structure interface is unchanged)
            Group {
                if model.calculation == .electronic {
                    SidebarView(model: model)
                } else {
                    NuclearSidebarView(model: model.nuclear)
                }
            }
            .safeAreaInset(edge: .top, spacing: 0) {
                CalculationKindPicker(model: model)
            }
            // 16.0: wide enough for the segmented Hessian control (Initial | CalcFC | CalcAll | RecalcFC)
            .navigationSplitViewColumnWidth(min: 340, ideal: 340, max: 420)
        } detail: {
            if model.calculation == .electronic {
                DetailView(model: model)
            } else {
                NuclearDetailView(model: model.nuclear)
            }
        }
        .navigationTitle("ΨLab")
        .background(WindowFitter())
        .preferredColorScheme((AppearanceChoice(rawValue: appearanceRaw) ?? .system).colorScheme)
        .task {
            // 20.0: the nuclear-motion interface at launch (UI test): run it if asked, and no
            // electronic-structure calculation
            if model.calculation == .nuclear {
                if UserDefaults.standard.bool(forKey: "psilabNuclearAutoRun") { await model.nuclear.run() }
                return
            }
            // Show a result right away (the Szabo–Ostlund defaults, or the molecule given as a
            // launch argument). Large general calculations wait for Run.
            if model.result == nil && !model.openedAtLaunch {
                if model.isGeneral {
                    await model.refreshEstimate()
                    if model.refusal == nil && (!model.isLongCalculation || model.autoRunLongJobs) {
                        await model.run()
                    }
                } else {
                    await model.run()
                }
            }
        }
    }
}

/// 20.0: the two kinds of calculation at the top of the sidebar, drawn like the tab bar of the
/// detail column: a glass capsule, the chosen kind on a lighter capsule; locked while a
/// calculation runs.
struct CalculationKindPicker: View {
    @Bindable var model: AppModel
    @Namespace private var selection

    private var locked: Bool { model.isRunning || model.isScanning || model.nuclear.isRunning }

    var body: some View {
        HStack(spacing: 2) {
            ForEach(CalculationKind.allCases) { kind in
                let selected = model.calculation == kind
                Button {
                    withAnimation(.snappy(duration: 0.25)) {
                        model.calculation = kind
                    }
                } label: {
                    Text(kind.title)
                        .font(.system(size: 13, weight: selected ? .semibold : .medium))
                        .foregroundStyle(selected ? Color.primary : Color.secondary)
                        .lineLimit(1)
                        .frame(maxWidth: .infinity)
                        .padding(.vertical, 6)
                        .padding(.horizontal, 10)
                        .background {
                            if selected {
                                Capsule()
                                    .fill(Color.primary.opacity(0.13))
                                    .matchedGeometryEffect(id: "selected", in: selection)
                            }
                        }
                        .contentShape(Capsule())
                }
                .buttonStyle(.plain)
                .help(kind.help)
                .accessibilityLabel(kind.title)
                .accessibilityAddTraits(selected ? .isSelected : [])
            }
        }
        .padding(3)
        .glassCapsule()
        .disabled(locked)
        .opacity(locked ? 0.6 : 1)
        .accessibilityElement(children: .contain)
        .accessibilityLabel("Calculation")
        .padding(.horizontal, 14)
        .padding(.top, 6)
        .padding(.bottom, 4)
    }
}

/// Keeps the window on the screen: the default size (1180 × 760 pt) is larger than a 1024-pt
/// screen, and a restored frame can lie partly outside it. When the window appears, its frame
/// is clipped to the screen's visible area (the minimum width, 1000 pt, still fits).
struct WindowFitter: NSViewRepresentable {
    func makeNSView(context: Context) -> NSView { WindowFitterView() }
    func updateNSView(_ nsView: NSView, context: Context) {}
}

final class WindowFitterView: NSView {
    override func viewDidMoveToWindow() {
        super.viewDidMoveToWindow()
        fit()
        // SwiftUI may still apply the default size after the view is attached.
        Task { @MainActor [weak self] in self?.fit() }
    }

    private func fit() {
        guard let window, let screen = window.screen ?? NSScreen.main else { return }
        let visible = screen.visibleFrame
        var frame = window.frame
        frame.size.width = min(frame.width, visible.width)
        frame.size.height = min(frame.height, visible.height)
        frame.origin.x = min(max(frame.minX, visible.minX), visible.maxX - frame.width)
        frame.origin.y = min(max(frame.minY, visible.minY), visible.maxY - frame.height)
        if frame != window.frame {
            window.setFrame(frame, display: true)
        }
    }
}

/// Six tabs (they fit next to the 340-pt sidebar at the minimum window width of 1000 pt, about 600 pt): Summary,
/// Methods (Hartree–Fock, coupled cluster, Full CI), Optimisation, Frequencies (with the
/// thermochemistry page), Scan (relaxed scan and potential curve), Output.
struct DetailView: View {
    @Bindable var model: AppModel

    // 15.0: the TabView is the root of the detail column again (as in 12.0), so that macOS places
    // the tab bar in the window's toolbar. 13.0 had wrapped it in a VStack for the banner of an
    // opened output file, which moved the tabs below a title bar; the banner is now an inset of
    // every tab's content.
    var body: some View {
        TabView(selection: $model.selectedTab) {
            SummaryView(model: model)
                .modifier(OpenedFileInset(model: model))
                .tabItem { Label("Summary", systemImage: "list.bullet.rectangle") }
                .tag(DetailTab.summary)

            MethodsView(model: model)
                .modifier(OpenedFileInset(model: model))
                .tabItem { Label("Methods", systemImage: "square.grid.2x2") }
                .tag(DetailTab.methods)

            OptimizationView(model: model)
                .modifier(OpenedFileInset(model: model))
                .tabItem { Label("Optimisation", systemImage: "arrow.down.right.and.arrow.up.left") }
                .tag(DetailTab.optimization)

            FrequenciesView(model: model)
                .modifier(OpenedFileInset(model: model))
                .tabItem { Label("Frequencies", systemImage: "waveform.path") }
                .tag(DetailTab.frequencies)

            ScanTabView(model: model)
                .modifier(OpenedFileInset(model: model))
                .tabItem { Label("Scan", systemImage: "chart.xyaxis.line") }
                .tag(DetailTab.potentialCurve)

            OutputView(model: model)
                .modifier(OpenedFileInset(model: model))
                .tabItem { Label("Output", systemImage: "doc.plaintext") }
                .tag(DetailTab.output)
        }
        .padding(.top, 4)
    }
}

/// 15.0: the banner of an opened output file above the content of a tab; 18.0: none for the
/// start-up example (the H₂ Full-CI job of the app's resources).
struct OpenedFileInset: ViewModifier {
    @Bindable var model: AppModel

    func body(content: Content) -> some View {
        content.safeAreaInset(edge: .top, spacing: 0) {
            if let file = model.openedFile, model.result != nil, !file.example {
                OpenedFileBanner(file: file)
            }
        }
    }
}

// MARK: - Shared small views

/// 13.0: the results shown are those of an opened output file.
struct OpenedFileBanner: View {
    var file: OpenedOutputFile

    var body: some View {
        let detail = [file.program, file.written.isEmpty ? "" : "written " + file.written]
            .filter { !$0.isEmpty }.joined(separator: ", ")
        let text: String = "Results read from “\(file.name)”" + (detail.isEmpty ? "" : " (\(detail))")
            + ". Run calculates the inputs of the sidebar."
        Label {
            Text(text)
                .font(.callout)
                .textSelection(.enabled)
                .frame(maxWidth: .infinity, alignment: .leading)
        } icon: {
            Image(systemName: "doc.text.magnifyingglass")
                .foregroundStyle(Color.accentColor)
        }
        .padding(.horizontal, 14)
        .padding(.vertical, 6)
        .background(Color.accentColor.opacity(0.08))
        .accessibilityElement(children: .combine)
    }
}

/// Banner in the detail area before a result exists: why the current settings cannot be calculated (red), or that the
/// calculation will take long (orange).
struct StatusBanner: View {
    enum Kind {
        case refusal, warning
    }

    var kind: Kind
    var text: String

    var body: some View {
        let color: Color = kind == .refusal ? .red : .orange
        Label {
            VStack(alignment: .leading, spacing: 2) {
                Text(kind == .refusal ? "Not possible with these settings" : "Long calculation")
                    .font(.headline)
                Text(text)
                    .textSelection(.enabled)
            }
            .frame(maxWidth: .infinity, alignment: .leading)
        } icon: {
            Image(systemName: kind == .refusal ? "xmark.octagon.fill" : "clock.badge.exclamationmark")
                .foregroundStyle(color)
        }
        .padding(12)
        .background(color.opacity(0.12), in: RoundedRectangle(cornerRadius: 10))
        .overlay(RoundedRectangle(cornerRadius: 10).strokeBorder(color.opacity(0.35), lineWidth: 1))
    }
}

/// Red banner for errors returned by the backend (or by the wrapper).
struct ErrorBanner: View {
    var message: String

    var body: some View {
        Label {
            Text(message)
                .textSelection(.enabled)
                .frame(maxWidth: .infinity, alignment: .leading)
        } icon: {
            Image(systemName: "exclamationmark.triangle.fill")
                .foregroundStyle(.red)
        }
        .padding(12)
        .background(Color.red.opacity(0.12), in: RoundedRectangle(cornerRadius: 10))
        .overlay(RoundedRectangle(cornerRadius: 10).strokeBorder(Color.red.opacity(0.35), lineWidth: 1))
    }
}

/// Section heading used inside the detail tabs.
struct SectionTitle: View {
    var text: String

    init(_ text: String) {
        self.text = text
    }

    var body: some View {
        Text(text)
            .font(.title3.weight(.semibold))
    }
}

/// Placeholder when no result is available yet.
struct NoResultView: View {
    var model: AppModel

    var body: some View {
        VStack(spacing: 14) {
            if let message = model.errorMessage {
                ErrorBanner(message: message)
                    .frame(maxWidth: 560)
            }
            if model.isGeneral, let refusal = model.refusal {
                StatusBanner(kind: .refusal, text: refusal)
                    .frame(maxWidth: 560)
            } else if model.isGeneral && model.isLongCalculation {
                StatusBanner(kind: .warning, text: model.longRunMessage + " Run asks for confirmation; a running calculation can be cancelled.")
                    .frame(maxWidth: 560)
            }
            ContentUnavailableView(
                "No results yet",
                systemImage: "atom",
                description: Text(model.isGeneral
                                  ? "Choose the job, method and basis set in the sidebar and press Run."
                                  : "Press Run to start a Hartree–Fock and Full-CI calculation (Szabo–Ostlund).")
            )
        }
        .padding()
    }
}
