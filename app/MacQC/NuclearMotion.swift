//
//  NuclearMotion.swift
//  MacQC
//
//  ΨLab 20.0/21.0: nuclear-motion calculations (rovibrational energy levels of H₂¹⁶O, H₂¹⁷O, H₂¹⁸O
//  and D₂¹⁶O on the CVRQD potential energy surfaces), the second kind of calculation next to the
//  electronic structure. The backend (nuclear_motion.cpp, macqc_nuclear_run_job) builds the DVR
//  Hamiltonian of the original program cpp_DOPI_expval_array.cpp (21.0: in Jacobi or Radau
//  coordinates) and finds its lowest eigenvalues (and eigenvectors). This file holds the input and
//  result types and the model of the interface.
//

import Foundation
import Observation
import AppKit
import UniformTypeIdentifiers

/// The two kinds of calculation of the app (the switch at the top of the sidebar).
enum CalculationKind: String, CaseIterable, Identifiable {
    case electronic
    case nuclear

    var id: String { rawValue }

    var title: String {
        switch self {
        case .electronic: return "Electronic structure"
        case .nuclear: return "Nuclear motion"
        }
    }

    var help: String {
        switch self {
        case .electronic: return "Electronic structure calculation: Hartree–Fock, perturbation theory, CI, coupled cluster; energies, geometries, frequencies."
        case .nuclear: return "Nuclear-motion theory calculation: vibrational (jmax = 0) and rovibrational energy levels of the water isotopologues on the CVRQD potential energy surfaces."
        }
    }
}

nonisolated enum Isotopologue: String, CaseIterable, Identifiable, Codable, Sendable {
    case h216o = "H2(16)O"
    case h217o = "H2(17)O"
    case h218o = "H2(18)O"
    case d216o = "D2(16)O"

    var id: String { rawValue }

    /// LaTeX for the interface (LaTeXText).
    var latex: String {
        switch self {
        case .h216o: return "H$_2$$^{16}$O"
        case .h217o: return "H$_2$$^{17}$O"
        case .h218o: return "H$_2$$^{18}$O"
        case .d216o: return "D$_2$$^{16}$O"
        }
    }

    /// Plain text with Unicode sub- and superscripts (menus).
    var title: String {
        switch self {
        case .h216o: return "H₂¹⁶O"
        case .h217o: return "H₂¹⁷O"
        case .h218o: return "H₂¹⁸O"
        case .d216o: return "D₂¹⁶O"
        }
    }

    var note: String {
        switch self {
        case .h216o: return "10 electrons, X ^1A_1, closed shell, C_{2v}; nuclear masses m(H) = 1.00727647 u, m(^{16}O) = 15.990526 u"
        case .h217o: return "10 electrons, X ^1A_1, closed shell, C_{2v}; nuclear masses m(H) = 1.00727647 u, m(^{17}O) = 16.9947425 u"
        case .h218o: return "10 electrons, X ^1A_1, closed shell, C_{2v}; nuclear masses m(H) = 1.00727647 u, m(^{18}O) = 17.999161 u"
        case .d216o: return "10 electrons, X ^1A_1, closed shell, C_{2v}; nuclear masses m(D) = 2.013553214 u, m(^{16}O) = 15.990526 u"
        }
    }
}

nonisolated enum NuclearEigensolver: String, CaseIterable, Identifiable, Codable, Sendable {
    case lanczos = "lanczos"
    case jacobiDavidson = "jacobi-davidson"
    case davidson = "davidson"
    case inversePower = "inverse-power"          // 21.0 (original: invpower)
    case rayleighQuotient = "rayleigh-quotient"  // 21.0 (original: RQiteration)
    case dense = "dense"

    var id: String { rawValue }

    var title: String {
        switch self {
        case .lanczos: return "Lánczos"
        case .jacobiDavidson: return "Jacobi–Davidson"
        case .davidson: return "Davidson"
        case .inversePower: return "Inverse power iteration"
        case .rayleighQuotient: return "Rayleigh quotient iteration"
        case .dense: return "Dense"
        }
    }

    /// Iterative methods that solve linear systems (all but the dense diagonalisation).
    var usesLinearSolver: Bool { self != .dense }

    var explanation: String {
        switch self {
        case .lanczos:
            return "Shift-invert Lánczos with full reorthogonalisation: every step solves H x = q, so the smallest eigenvalues converge first (as in the original program)."
        case .jacobiDavidson:
            return "Jacobi–Davidson: the correction equation is solved approximately (the preconditioner of H − θ, then steps of the linear solver on the projected equation); the lowest levels one after the other."
        case .davidson:
            return "Davidson: the correction t = −M⁻¹r with the preconditioner M ≈ H − θ, or a few steps of the linear solver on (H − θ) t = −r."
        case .inversePower:
            return "Inverse power iteration (thesis §6.1, the original's invpower): x ← H⁻¹x, one linear solve per step; converges with the factor E₁/E₂. For several levels a block with Rayleigh–Ritz (inverse subspace iteration)."
        case .rayleighQuotient:
            return "Rayleigh quotient iteration (the original's RQiteration): (H − θ) y = x with θ = xᵀHx, cubic convergence; the levels one after the other with deflation. MINRES or BiCG solve the indefinite systems; with CG or Chebyshev the shift is θ − ‖r‖."
        case .dense:
            return "Dense diagonalisation of the whole matrix (Householder and QL): exact reference for small bases (at most 4000 functions)."
        }
    }
}

nonisolated enum NuclearLinearSolver: String, CaseIterable, Identifiable, Codable, Sendable {
    case pcg = "pcg"
    case minres = "minres"
    case chebyshev = "chebyshev"   // 21.0 (original: chebyshev)
    case bicg = "bicg"             // 21.0 (original: BiCGM)
    var id: String { rawValue }
    var title: String {
        switch self {
        case .pcg: return "CG"
        case .minres: return "MINRES"
        case .chebyshev: return "Chebyshev"
        case .bicg: return "BiCG"
        }
    }
    var longTitle: String {
        switch self {
        case .pcg: return "Preconditioned conjugate gradients"
        case .minres: return "MINRES"
        case .chebyshev: return "Chebyshev iteration"
        case .bicg: return "Biconjugate gradients"
        }
    }
}

/// 21.0: the coordinates of the DVR (the original's nkoord: 1 = Jacobi, 2 = Radau).
nonisolated enum NuclearCoordinates: String, CaseIterable, Identifiable, Codable, Sendable {
    case jacobi = "jacobi"
    case radau = "radau"
    var id: String { rawValue }
    var title: String { self == .jacobi ? "Jacobi" : "Radau" }
    /// R₁, R₂ and θ in words.
    var explanation: String {
        switch self {
        case .jacobi:
            return "R₁ = r(O–H₁), R₂ from the centre of mass of O–H₁ to H₂, θ the angle between them."
        case .radau:
            return "R₁, R₂ from the canonical point K (near O, on the line from O through the centre of mass) to H₁ and H₂, θ = ∠H₁KH₂; reduced masses m(H₁), m(H₂)."
        }
    }
}

nonisolated enum NuclearPreconditioner: String, CaseIterable, Identifiable, Codable, Sendable {
    case ilu = "ilu"
    case jacobi = "jacobi"
    case none = "none"
    var id: String { rawValue }
    var title: String {
        switch self {
        case .ilu: return "ILU"
        case .jacobi: return "Jacobi"
        case .none: return "None"
        }
    }
}

/// Input JSON of macqc_nuclear_run_job.
nonisolated struct NuclearInput: Codable, Sendable, Equatable {
    var isotopologue: Isotopologue = .h216o
    var J: Int = 0                     // jmax of the original program
    var parityBlock: Int = 0           // J > 0: 0 = K from 0 to J, 1 = K from 1 to J
    var radialPoints: Int = 20         // ns1
    var angularPoints: Int = 24        // np
    var rCentre: Double = 2.0          // bohr (Re)
    var rHalfWidth: Double = 0.8       // bohr (Rmax)
    var levels: Int = 10               // nmax
    var coordinates: NuclearCoordinates = .jacobi  // 21.0 (nkoord)
    var eigenvectors: Bool = true      // 21.0: false = eigenvalues only (sv = 0)
    var gaussianStart: Bool = false    // 21.0: the start vector of Algorithm 5 (psikiind)
    var eigensolver: NuclearEigensolver = .lanczos
    var linearSolver: NuclearLinearSolver = .pcg
    var preconditioner: NuclearPreconditioner = .ilu
    var eigenvalueTolerance: Double = 1e-10
    var residualTolerance: Double = 1e-8
    var linearTolerance: Double = 1e-12
    var maxIterations: Int = 1000
    var maxLinearIterations: Int = 10000
    var jdInnerIterations: Int = 15
    var davidsonInnerIterations: Int = 0   // 21.0
    var subspaceSize: Int = 0
    var expectationValues: Bool = true
    var threads: Int = 0

    /// K of the blocks (as the backend).
    var kValues: [Int] {
        if J == 0 { return [0] }
        return Array((parityBlock == 0 ? 0 : 1)...J)
    }

    /// Size of the Hamiltonian: blocks × ns1² × np.
    var dimension: Int { kValues.count * radialPoints * radialPoints * angularPoints }
}

nonisolated struct NuclearMoments: Codable, Sendable, Equatable {
    var mean: Double
    var rms: Double
    var cubic: Double
    var inverse: Double
    var inverse2: Double
    var inverse3: Double
}

nonisolated struct NuclearLevelResult: Codable, Sendable, Equatable {
    var energy: Double        // E_h from the minimum of the surface
    var residual: Double?     // nil: not known (Lánczos, eigenvalues only)
    var r1: NuclearMoments    // bohr
    var r2: NuclearMoments
    var theta: Double         // degrees
    var cosTheta: Double
    var kWeights: [Double]
    var wavefunction: NuclearWavefunction?  // 21.0
}

/// 21.0: an eigenvector on the DVR grid (backend: WavefunctionPlot). ψ at the grid points scaled to
/// a largest value of +1 at the point (block, i, j, l); cuts through that point in the planes
/// (R₁, R₂), (R₁, θ) and (R₂, θ), and the probability densities integrated over the third
/// coordinate (bohr⁻², bohr⁻¹ rad⁻¹). Base64 float32 arrays, row-major: R₁R₂ ns × ns (R₁ the row),
/// R₁θ and R₂θ ns × np (the distance the row).
nonisolated struct NuclearWavefunction: Codable, Sendable, Equatable {
    var block: Int
    var i: Int
    var j: Int
    var l: Int
    var psiR1R2: String
    var psiR1Theta: String
    var psiR2Theta: String
    var rhoR1R2: String
    var rhoR1Theta: String
    var rhoR2Theta: String

    static func floats(_ base64: String) -> [Float] {
        guard let data = Data(base64Encoded: base64) else { return [] }
        var out = [Float](repeating: 0, count: data.count / MemoryLayout<Float>.size)
        _ = out.withUnsafeMutableBytes { data.copyBytes(to: $0) }
        return out
    }
}

nonisolated struct NuclearStep: Codable, Sendable, Equatable {
    var iteration: Int
    var lowest: Double
    var change: Double?
    var converged: Int
    var linearIterations: Int
    var level: Int?            // 21.0: the level being refined (0: all at once)
    var linearResidual: Double?
    var seconds: Double?
}

/// 21.0: the residual norms of one linear solve (Convergence tab).
nonisolated struct NuclearLinearTrace: Codable, Sendable, Equatable {
    var step: Int
    var residuals: [Double]
}

/// 21.0: the Gaussian start vector of Algorithm 5.
nonisolated struct NuclearStartVector: Codable, Sendable, Equatable {
    var rayleigh: Double       // E_h
    var overlap: Double?       // |⟨q₁|x₁⟩|
    var mu: [Double]           // r₁, r₂ (bohr), angle (°)
    var sigma: [Double]
}

/// Result JSON of macqc_nuclear_run_job.
nonisolated struct NuclearResult: Codable, Sendable, Equatable {
    var error: String?
    var input: NuclearInput?
    var isotopologue: String?
    var masses: [Double]?
    var mu1: Double?
    var mu2: Double?
    var g1: Double?                         // 21.0
    var g2: Double?
    var hasVectors: Bool?                   // 21.0
    var startVector: NuclearStartVector?    // 21.0
    var chebyshevInterval: [Double]?        // 21.0
    var firstSolve: NuclearLinearTrace?     // 21.0
    var lastSolve: NuclearLinearTrace?
    var potentialMinimum: Double?
    var minimumR: Double?
    var minimumAngle: Double?
    var rGrid: [Double]?
    var rWeights: [Double]?                 // 21.0
    var kValues: [Int]?
    var cosGrid: [[Double]]?
    var dimension: Int?
    var nonzeros: Int?
    var hamiltonianSeconds: Double?
    var solverSeconds: Double?
    var iterations: Int?
    var matrixVectorProducts: Int?
    var linearIterations: Int?
    var converged: Bool?
    var solver: String?
    var history: [NuclearStep]?
    var levels: [NuclearLevelResult]
    var matrixWhole: NuclearMatrixPicture?     // 20.0: the structure of the Hamiltonian
    var matrixDetail: NuclearMatrixPicture?
    var notes: [String]?
    var output: String?
    var exportId: String?                   // 21.0: macqc_nuclear_export of this result
}

nonisolated struct NuclearEstimate: Codable, Sendable, Equatable {
    var error: String?
    var dimension: Int?
    var nonzeros: Int?
    var memoryMB: Double?
    var refusal: String?
}

/// 22.0: the base64 text of the matrix pictures in the archive of an output file: "rle:" and the
/// text with every run of n ≥ 5 'A' (zero bytes) written "~n." (backend: base64_zero_runs).
nonisolated enum NuclearBase64 {
    static func expand(_ s: String) -> String {
        guard s.hasPrefix("rle:") else { return s }
        var out = ""
        out.reserveCapacity(s.utf8.count * 4)
        var rest = Substring(s.dropFirst(4))
        while let tilde = rest.firstIndex(of: "~") {
            out.append(contentsOf: rest[..<tilde])
            let after = rest[rest.index(after: tilde)...]
            guard let dot = after.firstIndex(of: "."), let n = Int(after[..<dot]), n >= 0 else { return out }
            out += String(repeating: "A", count: n)
            rest = after[after.index(after: dot)...]
        }
        out.append(contentsOf: rest)
        return out
    }
}

extension NuclearMatrixPicture {
    /// 22.0: the pictures of an archive with their zero runs written out (NuclearBase64).
    mutating func expandZeroRuns() {
        density = NuclearBase64.expand(density)
        magnitude = NuclearBase64.expand(magnitude)
    }
}

nonisolated enum NuclearUnits {
    static let hartreeToCm = 219474.6313632     // CODATA 2018 (as the backend)
    static let bohrToAngstrom = 0.529177210903
}

nonisolated enum NuclearEngine {
    static func run(_ input: NuclearInput, job: MacQCJob?) async throws -> NuclearResult {
        let task = Task.detached(priority: .userInitiated) { () throws -> NuclearResult in
            let data = try JSONEncoder().encode(input)
            guard let json = String(data: data, encoding: .utf8) else { throw MacQCEngineError.encodingFailed }
            guard let raw = macqc_nuclear_run_job(json, job?.pointer) else { throw MacQCEngineError.backendReturnedNull }
            defer { macqc_free(raw) }
            let text = String(cString: raw)
            do {
                return try JSONDecoder().decode(NuclearResult.self, from: Data(text.utf8))
            } catch {
                throw MacQCEngineError.decodingFailed(String(describing: error))
            }
        }
        return try await task.value
    }

    /// 21.0: writes the eigenvalues, eigenvectors and/or Gaussian start vector of the last result
    /// (macqc_nuclear_export); nil on success, else the message.
    static func export(id: String, path: String, eigenvalues: Bool, eigenvectors: Bool, startVector: Bool) -> String? {
        let request: [String: Any] = ["exportId": id, "path": path, "eigenvalues": eigenvalues,
                                      "eigenvectors": eigenvectors, "startVector": startVector]
        guard let data = try? JSONSerialization.data(withJSONObject: request),
              let json = String(data: data, encoding: .utf8),
              let raw = macqc_nuclear_export(json) else { return "The export failed." }
        defer { macqc_free(raw) }
        guard let reply = try? JSONSerialization.jsonObject(with: Data(String(cString: raw).utf8)) as? [String: Any] else {
            return "The export failed."
        }
        return reply["error"] as? String
    }

    static func estimate(_ input: NuclearInput) -> NuclearEstimate {
        guard let data = try? JSONEncoder().encode(input), let json = String(data: data, encoding: .utf8),
              let raw = macqc_nuclear_estimate(json) else { return NuclearEstimate(refusal: "The estimate failed.") }
        defer { macqc_free(raw) }
        return (try? JSONDecoder().decode(NuclearEstimate.self, from: Data(String(cString: raw).utf8)))
            ?? NuclearEstimate(refusal: "The estimate could not be read.")
    }
}

/// Tabs of the nuclear-motion interface.
enum NuclearTab: Hashable {
    case levels, expectation, matrix, eigenvectors, convergence, output
}

/// The model of the nuclear-motion interface: its inputs, the running job and the result.
@MainActor
@Observable
final class NuclearModel {
    var input = NuclearInput()
    var result: NuclearResult?
    var errorMessage: String?
    var isRunning = false
    var progressLine = ""
    var outputText = ""
    var outputGeneration = 0
    var selectedTab: NuclearTab = .levels
    var runStarted: Date?
    var lastRunInput: NuclearInput?
    var showExportSheet = false        // 21.0: File ▸ Export Eigenvalues and Eigenvectors…
    /// 22.0: the result shown is that of an output file (File ▸ Open Output File…), not of a run.
    var openedFile: OpenedOutputFile?
    /// 23.0: the startup example is being read (NuclearPlaceholder says so).
    var loadingExample = false
    @ObservationIgnored private var runningJob: MacQCJob?

    /// `launchArguments`: take the UI test's launch arguments (23.0: the first window only).
    init(launchArguments: Bool = true) {
        guard launchArguments else { return }
        // launch arguments of the UI test: -psilabNuclearNs1 8 -psilabNuclearNp 10 -psilabNuclearJmax 0 …
        let d = UserDefaults.standard
        if d.object(forKey: "psilabNuclearNs1") != nil { input.radialPoints = max(2, d.integer(forKey: "psilabNuclearNs1")) }
        if d.object(forKey: "psilabNuclearNp") != nil { input.angularPoints = max(2, d.integer(forKey: "psilabNuclearNp")) }
        if d.object(forKey: "psilabNuclearJmax") != nil { input.J = max(0, d.integer(forKey: "psilabNuclearJmax")) }
        if d.object(forKey: "psilabNuclearLevels") != nil { input.levels = max(1, d.integer(forKey: "psilabNuclearLevels")) }
        if let s = d.string(forKey: "psilabNuclearSolver"), let e = NuclearEigensolver(rawValue: s) { input.eigensolver = e }
        // 21.0: -psilabNuclearCoordinates radau -psilabNuclearLinear minres -psilabNuclearGaussian YES …
        if let s = d.string(forKey: "psilabNuclearCoordinates"), let c = NuclearCoordinates(rawValue: s) { input.coordinates = c }
        if let s = d.string(forKey: "psilabNuclearLinear"), let l = NuclearLinearSolver(rawValue: s) { input.linearSolver = l }
        if d.object(forKey: "psilabNuclearGaussian") != nil { input.gaussianStart = d.bool(forKey: "psilabNuclearGaussian") }
        if d.object(forKey: "psilabNuclearEigenvectors") != nil { input.eigenvectors = d.bool(forKey: "psilabNuclearEigenvectors") }
        if let s = d.string(forKey: "psilabNuclearIsotopologue"),
           let i = Isotopologue.allCases.first(where: { $0.rawValue == s || $0.title == s }) { input.isotopologue = i }
    }

    /// The size of the calculation for the current inputs.
    var estimate: NuclearEstimate { NuclearEngine.estimate(input) }

    var refusal: String? {
        let e = estimate
        return e.refusal ?? e.error
    }

    func requestRun() {
        guard !isRunning else { return }
        Task { await run() }
    }

    func run() async {
        guard !isRunning, refusal == nil else { return }
        isRunning = true
        errorMessage = nil
        progressLine = ""
        outputText = ""
        outputGeneration += 1
        runStarted = Date()
        let input = self.input
        let handle = MacQCJob()
        runningJob = handle
        let poller = Task { [weak self] in
            var offset = 0
            while !Task.isCancelled {
                let line = handle.progressText()
                let piece = handle.output(from: offset)
                offset = max(offset, piece.total)
                guard let self else { return }
                self.progressLine = line
                if !piece.text.isEmpty { self.outputText += piece.text }
                try? await Task.sleep(for: .milliseconds(250))
            }
        }
        defer {
            poller.cancel()
            isRunning = false
            runningJob = nil
            runStarted = nil
            progressLine = ""
        }
        do {
            let r = try await NuclearEngine.run(input, job: handle)
            if let full = r.output, !full.isEmpty, full != outputText {
                outputText = full
                outputGeneration += 1
            }
            if r.error == "calculation cancelled" {
                errorMessage = "Calculation cancelled."
                return
            }
            if let e = r.error {
                errorMessage = e
                return
            }
            result = r
            lastRunInput = input
            openedFile = nil
        } catch {
            errorMessage = error.localizedDescription
        }
    }

    func cancel() {
        runningJob?.cancel()
    }

    /// 22.0: shows the result of an output file (its archive section) in every tab, as after the
    /// calculation; `text`: the file, shown in the Output tab. The inputs of the sidebar stay.
    func showOpenedResult(_ r: NuclearResult, text: String?, file: OpenedOutputFile) {
        var r = r
        r.matrixWhole?.expandZeroRuns()
        r.matrixDetail?.expandZeroRuns()
        if let text, !text.isEmpty { r.output = text }
        r.exportId = nil   // the eigenvectors themselves are not in the file
        result = r
        errorMessage = r.error
        outputText = r.output ?? ""
        outputGeneration += 1
        lastRunInput = r.input
        openedFile = file
        selectedTab = .levels
    }

    // MARK: 23.0: the startup example

    /// The output file of the startup example in the app's resources (examples/outputs): H₂¹⁶O on
    /// the CVRQD surface, J = 10 with the block K = 1 … 10, Radau coordinates, ns1 = np = 50
    /// (1 250 000 basis functions), the 100 lowest levels with their eigenvectors and expectation
    /// values by shift-invert Lánczos (conjugate gradients, ILU, the Gaussian start vector);
    /// 10 h 01 min on 12 threads with ΨLab 22.0.
    static let startupExampleName = "h216o-j10-radau-ns1-50-np-50"

    /// The example read once per launch (8.6 MB of text, decoded off the main thread) and shared by
    /// the models of all windows and tabs.
    private static var exampleTask: Task<(result: NuclearResult, text: String, header: ArchiveHeader?)?, Never>?

    /// The result, the text and the archive header of the startup example, nil if it cannot be read.
    nonisolated static func readExample(_ url: URL) -> (result: NuclearResult, text: String, header: ArchiveHeader?)? {
        guard let text = try? String(contentsOf: url, encoding: .utf8),
              let json = OutputArchive.json(in: text),
              let r = try? JSONDecoder().decode(NuclearResult.self, from: Data(json.utf8)) else { return nil }
        let header = try? JSONDecoder().decode(ArchiveHeader.self, from: Data(json.utf8))
        return (r, text, header)
    }

    /// 23.0: the nuclear-motion interface opens with the results of the startup example in every
    /// tab (Levels, Geometry, Matrix, Eigenvectors, Convergence, Output) and its inputs in the
    /// sidebar, so that Run repeats it (as the H₂ example of the electronic-structure interface,
    /// 17.0). Nothing is replaced if a result is shown or a calculation runs in the meantime.
    func showStartupExample() {
        guard result == nil,
              let url = Bundle.main.url(forResource: Self.startupExampleName, withExtension: "out") else { return }
        if Self.exampleTask == nil {
            Self.exampleTask = Task.detached(priority: .userInitiated) { NuclearModel.readExample(url) }
        }
        guard let task = Self.exampleTask else { return }
        loadingExample = true
        Task { [weak self] in
            let loaded = await task.value
            guard let self else { return }
            self.loadingExample = false
            guard let loaded, self.result == nil, !self.isRunning else { return }
            if let i = loaded.result.input { self.input = i }
            self.showOpenedResult(loaded.result, text: loaded.text,
                                  file: OpenedOutputFile(name: url.lastPathComponent,
                                                         program: loaded.header?.archive?.program ?? "",
                                                         written: loaded.header?.archive?.written ?? "",
                                                         example: true))
        }
    }

    /// File ▸ Export Results as Text…: the output file of the last calculation.
    func exportText() {
        guard let text = result?.output, !text.isEmpty else { return }
        let panel = NSSavePanel()
        panel.title = "Export Results as Text"
        panel.allowedContentTypes = [UTType.plainText]
        panel.canCreateDirectories = true
        let i = result?.input ?? input
        let name = i.isotopologue.rawValue.replacingOccurrences(of: "(", with: "").replacingOccurrences(of: ")", with: "")
        panel.nameFieldStringValue = "PsiLab-\(name)-J\(i.J)-\(i.coordinates.rawValue)-ns1-\(i.radialPoints)-np-\(i.angularPoints).out"
        guard panel.runModal() == .OK, let url = panel.url else { return }
        do {
            try text.write(to: url, atomically: true, encoding: .utf8)
        } catch {
            errorMessage = "The file could not be saved: \(error.localizedDescription)"
        }
    }

    /// 21.0: what the last result can export.
    var canExportVectors: Bool { result?.hasVectors == true }
    var canExportStartVector: Bool { result?.startVector != nil }

    /// 21.0: the eigenvalues (the level table of the output file), the eigenvectors and/or the
    /// Gaussian start vector of the last result, as one text file.
    func exportLevels(eigenvalues: Bool, eigenvectors: Bool, startVector: Bool) {
        guard let r = result, let id = r.exportId else {
            errorMessage = "Run the calculation again to export its levels."
            return
        }
        let panel = NSSavePanel()
        panel.title = "Export Eigenvalues and Eigenvectors"
        panel.allowedContentTypes = [UTType.plainText]
        panel.canCreateDirectories = true
        let i = r.input ?? input
        let name = i.isotopologue.rawValue.replacingOccurrences(of: "(", with: "").replacingOccurrences(of: ")", with: "")
        var parts: [String] = []
        if eigenvalues { parts.append("eigenvalues") }
        if eigenvectors { parts.append("eigenvectors") }
        if startVector { parts.append("startvector") }
        panel.nameFieldStringValue = "PsiLab-\(name)-J\(i.J)-\(i.coordinates.rawValue)-ns1-\(i.radialPoints)-np-\(i.angularPoints)-\(parts.joined(separator: "-")).txt"
        guard panel.runModal() == .OK, let url = panel.url else { return }
        if let message = NuclearEngine.export(id: id, path: url.path, eigenvalues: eigenvalues,
                                              eigenvectors: eigenvectors, startVector: startVector) {
            errorMessage = "The export failed: \(message)"
        } else {
            errorMessage = nil
        }
    }

    func elapsedText(at now: Date) -> String {
        guard let start = runStarted else { return "" }
        let s = Int(now.timeIntervalSince(start))
        return String(format: "%d:%02d", s / 60, s % 60)
    }
}
