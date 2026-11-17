//
//  Molecules.swift
//  MacQC
//
//  The molecules, job types, methods and basis sets offered in the sidebar (ΨLab 5.0).
//

import Foundation

/// HeH⁺ is the Szabo–Ostlund STO-nG calculation of versions 1 and 2, or (ΨLab 7.0) a general
/// calculation in a Gaussian basis set of the library; H₂, LiH, H₂O, H₂O₂ and NH₃ use the
/// Gaussian basis sets (correlation-consistent; 14.0 also Pople and def2).
nonisolated enum MoleculeChoice: String, CaseIterable, Identifiable, Sendable {
    case heh = "HeH+"
    case h2 = "H2"
    case lih = "LiH"
    case h2o = "H2O"
    case h2o2 = "H2O2"
    case nh3 = "NH3"
    /// ΨLab 9.0: a molecule from a structure file (the app keeps them; see CustomMolecules.swift).
    case custom = "custom"

    var id: String { rawValue }

    /// The built-in molecules, in the order of the molecule menu.
    static var builtIn: [MoleculeChoice] { allCases.filter { $0 != .custom } }

    /// Name in the user interface.
    var title: String {
        switch self {
        case .heh: return "HeH⁺"
        case .h2: return "H₂"
        case .lih: return "LiH"
        case .h2o: return "H₂O"
        case .h2o2: return "H₂O₂"
        case .nh3: return "NH₃"
        case .custom: return "Custom molecule"
        }
    }

    var electrons: Int {
        switch self {
        case .heh, .h2: return 2
        case .lih: return 4
        case .h2o, .nh3: return 10
        case .h2o2: return 18
        case .custom: return 0
        }
    }

    /// Element symbols in the order of the backend's atoms (He, H; O, H, H; O, O, H, H; N, H, H, H).
    var atomSymbols: [String] {
        switch self {
        case .heh: return ["He", "H"]
        case .h2: return ["H", "H"]
        case .lih: return ["Li", "H"]
        case .h2o: return ["O", "H", "H"]
        case .h2o2: return ["O", "O", "H", "H"]
        case .nh3: return ["N", "H", "H", "H"]
        case .custom: return []
        }
    }

    var atomCount: Int { atomSymbols.count }

    /// Atom labels with their numbers: "O1", "H2", …
    var atomLabels: [String] {
        atomSymbols.enumerated().map { "\($0.element)\($0.offset + 1)" }
    }

    /// Label of the bond length, e.g. "R(O–H)".
    var bondLabel: String {
        switch self {
        case .heh: return "R(He–H)"
        case .h2: return "R(H–H)"
        case .lih: return "R(Li–H)"
        case .h2o, .h2o2: return "R(O–H)"
        case .nh3: return "R(N–H)"
        case .custom: return "R"
        }
    }

    /// Label of the bond angle of the inputs (nil for diatomics).
    var angleLabel: String? {
        switch self {
        case .h2o: return "Angle H–O–H"
        case .h2o2: return "Angle O–O–H"
        case .nh3: return "Angle H–N–H"
        default: return nil
        }
    }

    /// Equilibrium bond length in bohr (H₂ 0.74144 Å, LiH 1.5949 Å, H₂O 0.9572 Å, H₂O₂ 0.9670 Å,
    /// NH₃ 1.012 Å, as in the backend); HeH⁺: the Szabo–Ostlund value 1.4632 bohr.
    var defaultBondLength: Double {
        switch self {
        case .heh: return 1.4632
        case .h2: return 0.74144 / Units.bohrToAngstrom
        case .lih: return 1.5949 / Units.bohrToAngstrom
        case .h2o: return 0.9572 / Units.bohrToAngstrom
        case .h2o2: return 0.9670 / Units.bohrToAngstrom
        case .nh3: return 1.012 / Units.bohrToAngstrom
        case .custom: return 1.0
        }
    }

    /// Default bond angle in degrees (H₂O 104.52°, H₂O₂ O–O–H 102.32°, NH₃ 106.7°; 0 for diatomics).
    var defaultAngle: Double {
        switch self {
        case .h2o: return 104.52
        case .h2o2: return 102.32
        case .nh3: return 106.7
        default: return 0
        }
    }

    static let defaultWaterAngle = 104.52
    /// H₂O₂: R(O–O) 1.4556 Å and the dihedral angle H–O–O–H 113.70° (the backend's defaults).
    static let defaultPeroxideOO = 1.4556 / Units.bohrToAngstrom
    static let defaultPeroxideDihedral = 113.70

    /// Allowed range of the input angle (the backend's limits).
    var angleRange: ClosedRange<Double> {
        switch self {
        case .nh3: return 60...120
        case .h2o2: return 30...175
        default: return 30...180
        }
    }

    /// Potential-curve range in bohr (H₂O, H₂O₂, NH₃: all X–H bonds stretched together at the
    /// chosen angles).
    var scanRange: ClosedRange<Double> {
        switch self {
        case .heh: return 0.8...5.0
        case .h2: return (0.40 / Units.bohrToAngstrom)...(3.00 / Units.bohrToAngstrom)
        case .lih: return (1.00 / Units.bohrToAngstrom)...(4.00 / Units.bohrToAngstrom)
        case .h2o, .h2o2, .nh3, .custom: return (0.70 / Units.bohrToAngstrom)...(2.00 / Units.bohrToAngstrom)
        }
    }

    /// Short note shown under the molecule picker.
    var subtitle: String {
        switch self {
        case .heh: return "2 electrons, charge +1: STO-nG (Szabo–Ostlund, HF and Full CI) or a Gaussian basis set (all methods and jobs)"
        case .h2: return "2 electrons, closed shell"
        case .lih: return "4 electrons, closed shell"
        case .h2o: return "10 electrons, closed shell (Full CI impossible)"
        case .h2o2: return "18 electrons, closed shell, C₂ (Full CI impossible)"
        case .nh3: return "10 electrons, closed shell, C₃ᵥ (Full CI impossible)"
        case .custom: return "A molecule from a structure file: .xyz (Å; the comment line is its name, LaTeX allowed) or Gaussian input (.gjf, .com)"
        }
    }

    /// Scan coordinates offered in the sidebar for this molecule (the first is the default).
    var scanSuggestions: [ScanSuggestion] {
        switch self {
        case .heh, .h2, .lih:
            return [ScanSuggestion(kind: .bond, atoms: [1, 2], title: "Bond " + bondLabel, steps: 10, step: 0.05, start: nil)]
        case .h2o:
            return [ScanSuggestion(kind: .angle, atoms: [2, 1, 3], title: "H–O–H bend from 100° to 180°", steps: 8, step: 10, start: 100),
                    ScanSuggestion(kind: .bond, atoms: [1, 2], title: "One O–H bond", steps: 10, step: 0.05, start: nil)]
        case .h2o2:
            return [ScanSuggestion(kind: .dihedral, atoms: [3, 1, 2, 4], title: "H–O–O–H torsion from 0° to 180°", steps: 18, step: 10, start: 0),
                    ScanSuggestion(kind: .bond, atoms: [1, 2], title: "O–O bond", steps: 10, step: 0.05, start: nil),
                    ScanSuggestion(kind: .angle, atoms: [2, 1, 3], title: "One O–O–H angle", steps: 8, step: 5, start: nil)]
        case .nh3:
            return [ScanSuggestion(kind: .dihedral, atoms: [2, 1, 3, 4], title: "Umbrella: D(H2–N1–H3–H4) from −110° to −180°", steps: 7, step: -10, start: -110),
                    ScanSuggestion(kind: .angle, atoms: [2, 1, 3], title: "One H–N–H angle", steps: 8, step: 2, start: nil),
                    ScanSuggestion(kind: .bond, atoms: [1, 2], title: "One N–H bond", steps: 10, step: 0.05, start: nil)]
        case .custom:
            return [ScanSuggestion(kind: .bond, atoms: [1, 2], title: "Bond between atoms 1 and 2", steps: 10, step: 0.05, start: nil),
                    ScanSuggestion(kind: .angle, atoms: [2, 1, 3], title: "Angle of atoms 2–1–3", steps: 8, step: 5, start: nil),
                    ScanSuggestion(kind: .dihedral, atoms: [4, 1, 2, 3], title: "Dihedral of atoms 4–1–2–3 from 0° to 180°", steps: 18, step: 10, start: 0)]
        }
    }

    /// The default scan of this molecule.
    var defaultScan: ScanSettings { scanSuggestions.first.map { ScanSettings($0) } ?? ScanSettings() }
}

/// Type of the scanned coordinate (Gaussian's B, A and D).
nonisolated enum ScanKind: String, CaseIterable, Identifiable, Sendable {
    case bond = "B"
    case angle = "A"
    case dihedral = "D"

    var id: String { rawValue }

    var title: String {
        switch self {
        case .bond: return "Bond"
        case .angle: return "Angle"
        case .dihedral: return "Dihedral"
        }
    }

    var atomCount: Int {
        switch self {
        case .bond: return 2
        case .angle: return 3
        case .dihedral: return 4
        }
    }

    var unit: String { self == .bond ? "Å" : "°" }
}

/// A scan offered in the sidebar's suggestions menu.
nonisolated struct ScanSuggestion: Identifiable, Sendable, Hashable {
    var kind: ScanKind
    var atoms: [Int]          // numbers from 1
    var title: String
    var steps: Int
    var step: Double          // Å or degrees
    var start: Double?        // nil: the value at the input geometry

    var id: String { title }
}

/// The scanned coordinate and the scan values of the sidebar (Gaussian's ModRedundant line
/// "B|A|D N1 N2 [N3 [N4]] S nsteps stepsize"; steps in Å or degrees).
nonisolated struct ScanSettings: Equatable, Sendable {
    var kind: ScanKind = .bond
    var atoms: [Int] = [1, 2, 3, 4]   // numbers from 1; the first kind.atomCount are used
    var steps: Int = 10
    var step: Double = 0.05           // Å or degrees
    var fromInput = true              // start at the value of the input geometry
    var start: Double = 0             // Å or degrees (when fromInput is false)

    init() {}

    init(_ s: ScanSuggestion) {
        kind = s.kind
        atoms = s.atoms + [1, 2, 3, 4].filter { !s.atoms.contains($0) }
        steps = s.steps
        step = s.step
        fromInput = s.start == nil
        start = s.start ?? 0
    }

    /// The atoms of the coordinate (numbers from 1).
    var usedAtoms: [Int] { Array(atoms.prefix(kind.atomCount)) }

    /// The ModRedundant line, e.g. "D 3 1 2 4 S 18 10".
    var line: String {
        let atomText = usedAtoms.map(String.init).joined(separator: " ")
        return "\(kind.rawValue) \(atomText) S \(steps) \(ScanSettings.number(step))"
    }

    /// Error in the choice of atoms (repeated or missing atoms), or nil.
    func atomProblem(atomCount: Int) -> String? {
        if atomCount < kind.atomCount {
            return "\(kind.title) needs \(kind.atomCount) atoms; this molecule has \(atomCount)."
        }
        let used = usedAtoms
        if used.contains(where: { $0 < 1 || $0 > atomCount }) { return "Choose atoms of the molecule." }
        if Set(used).count != used.count { return "An atom appears twice in the coordinate." }
        return nil
    }

    /// Shortest exact text of a number ("0.05", "10", "-2.5").
    static func number(_ x: Double) -> String {
        var t = String(format: "%.6f", x)
        while t.hasSuffix("0") { t.removeLast() }
        if t.hasSuffix(".") { t.removeLast() }
        return t == "-0" ? "0" : t
    }

    /// Parses a ModRedundant line ("[B|A|D] N1 N2 [N3 [N4]] S nsteps stepsize"; launch argument).
    init?(line: String) {
        let tokens = line.split(whereSeparator: { $0 == " " || $0 == "," }).map(String.init)
        guard let s = tokens.firstIndex(where: { $0.uppercased() == "S" }), tokens.count == s + 3,
              let n = Int(tokens[s + 1]), let size = Double(tokens[s + 2]) else { return nil }
        var head = Array(tokens[..<s])
        let kind: ScanKind
        if let first = head.first, let k = ScanKind(rawValue: first.uppercased()) {
            kind = k
            head.removeFirst()
        } else {
            switch head.count {
            case 2: kind = .bond
            case 3: kind = .angle
            case 4: kind = .dihedral
            default: return nil
            }
        }
        let numbers = head.compactMap { Int($0) }
        guard numbers.count == head.count, numbers.count == kind.atomCount else { return nil }
        self.init()
        self.kind = kind
        atoms = numbers + [1, 2, 3, 4].filter { !numbers.contains($0) }
        steps = n
        step = size
    }
}

/// Methods of the general path, in the order of the method menu. The raw value is the
/// backend's method name; all correlated methods start from RHF and correlate all electrons.
nonisolated enum MethodChoice: String, CaseIterable, Identifiable, Sendable {
    case hf = "RHF"
    case ccsd = "CCSD"
    case ccsdT = "CCSD(T)"
    case ccsdt = "CCSDT"
    case ccsdtq = "CCSDTQ"
    case fci = "FCI"

    var id: String { rawValue }

    /// Name in the user interface.
    var title: String {
        switch self {
        case .hf: return "HF"
        case .ccsd: return "CCSD"
        case .ccsdT: return "CCSD(T)"
        case .ccsdt: return "CCSDT"
        case .ccsdtq: return "CCSDTQ"
        case .fci: return "Full CI"
        }
    }

    /// ΨLab 10.0: an analytic Hessian beyond RHF (Full CI, CCSD, CCSD(T)).
    var hasAnalyticCorrelatedHessian: Bool {
        switch self {
        case .ccsd, .ccsdT, .fci: return true
        case .hf, .ccsdt, .ccsdtq: return false
        }
    }

    var isCoupledCluster: Bool {
        switch self {
        case .ccsd, .ccsdT, .ccsdt, .ccsdtq: return true
        case .hf, .fci: return false
        }
    }

    /// One-line description under the method menu.
    var explanation: String {
        switch self {
        case .hf: return "Restricted Hartree–Fock (no electron correlation)."
        case .ccsd: return "Coupled cluster with single and double excitations."
        case .ccsdT: return "CCSD plus the perturbative triples correction (T)."
        case .ccsdt: return "Coupled cluster with single, double and triple excitations."
        case .ccsdtq: return "Coupled cluster up to quadruple excitations."
        case .fci: return "Full configuration interaction: exact in the basis set."
        }
    }

    /// Accepts the raw value or the title (launch arguments), in any case.
    init?(name: String) {
        let n = name.lowercased()
        guard let m = MethodChoice.allCases.first(where: { $0.rawValue.lowercased() == n || $0.title.lowercased() == n })
        else { return nil }
        self = m
    }

    /// Title of a backend method name ("RHF" -> "HF").
    static func title(of backendName: String) -> String {
        MethodChoice(rawValue: backendName)?.title ?? backendName
    }
}

/// Job types (ΨLab 5.0 and 7.0), in the order of the job menu. The raw value is the backend's
/// job name; the number is the backend's job number.
nonisolated enum JobChoice: String, CaseIterable, Identifiable, Sendable {
    case energy = "energy"
    case opt = "opt"
    case freq = "freq"
    case optFreq = "opt+freq"
    case optFreqSP = "opt+freq+sp"
    case scan = "scan"
    case ts = "ts"
    case tsFreq = "ts+freq"
    case tsFreqSP = "ts+freq+sp"

    var id: String { rawValue }

    /// Full name (job menu, output file).
    var title: String {
        switch self {
        case .energy: return "Energy"
        case .opt: return "Geometry optimisation"
        case .freq: return "Frequencies"
        case .optFreq: return "Optimisation + frequencies"
        case .optFreqSP: return "Optimisation + frequencies + single-point energy"
        case .scan: return "Relaxed potential-energy scan"
        case .ts: return "Transition-state optimisation"
        case .tsFreq: return "TS optimisation + frequencies"
        case .tsFreqSP: return "TS optimisation + frequencies + single-point energy"
        }
    }

    /// Short name for the closed menu in the sidebar and for the result headers.
    var shortTitle: String {
        switch self {
        case .energy: return "Energy"
        case .opt: return "Optimisation"
        case .freq: return "Frequencies"
        case .optFreq: return "OPT + FREQ"
        case .optFreqSP: return "OPT + FREQ + SP"
        case .scan: return "Relaxed Scan"
        case .ts: return "TS Optimisation"
        case .tsFreq: return "TS + FREQ"
        case .tsFreqSP: return "TS + FREQ + SP"
        }
    }

    var number: Int { (JobChoice.allCases.firstIndex(of: self) ?? 0) + 1 }

    /// Optimisation to a minimum or to a transition state.
    var hasOptimization: Bool { self == .opt || self == .optFreq || self == .optFreqSP || hasTransitionState }
    var hasTransitionState: Bool { self == .ts || self == .tsFreq || self == .tsFreqSP }
    var hasFrequencies: Bool { self == .freq || self == .optFreq || self == .optFreqSP || self == .tsFreq || self == .tsFreqSP }
    var hasSinglePoint: Bool { self == .optFreqSP || self == .tsFreqSP }
    var hasScan: Bool { self == .scan }
    /// Uses the optimiser's convergence criteria and maximum number of steps.
    var usesOptimizer: Bool { hasOptimization || hasScan }

    /// One-line description under the job menu.
    var explanation: String {
        switch self {
        case .energy:
            return "Single-point energy at the input geometry."
        case .opt:
            return "Minimum-energy geometry from the input geometry (redundant internal coordinates, RFO steps, BFGS)."
        case .freq:
            return "Harmonic frequencies, IR intensities and thermochemistry at the input geometry (a warning if it is not a stationary point)."
        case .optFreq:
            return "Optimisation, then frequencies and thermochemistry at the optimised geometry."
        case .optFreqSP:
            return "Optimisation and frequencies, then a single-point energy with its own method and basis set; the final energies add the thermal corrections."
        case .scan:
            return "Gaussian's Opt=ModRedundant scan: the chosen bond, angle or dihedral is frozen at each value (start, start + step, …) and every other coordinate is optimised. Diatomics: single points."
        case .ts:
            return "Saddle point of order one (Gaussian's Opt=TS): Berny steps in redundant internal coordinates, a computed Hessian (CalcFC or CalcAll), partitioned RFO along the lowest mode, Bofill updates, trust radius 0.30."
        case .tsFreq:
            return "Transition-state search, then frequencies: the job checks that there is exactly one imaginary frequency and shows its normal mode."
        case .tsFreqSP:
            return "Transition-state search and frequencies, then a single-point energy at the transition state with its own method and basis set."
        }
    }

    /// Accepts the raw value, the title or the number (launch arguments), in any case.
    init?(name: String) {
        let n = name.lowercased()
        guard let j = JobChoice.allCases.first(where: {
            $0.rawValue == n || $0.title.lowercased() == n || "\($0.number)" == n
        }) else { return nil }
        self = j
    }
}

/// A family of basis sets: a section of the basis-set menus.
nonisolated struct BasisFamily: Sendable, Hashable {
    let title: String
    let names: [String]
}

nonisolated enum BasisSets {
    /// The Gaussian94 files in the app bundle's basis folder (from basissetexchange.org), by family
    /// in the order of the backend's library (basis.cpp). 14.0: the Pople 6-31G family (Cartesian
    /// d functions, 6D, as Gaussian uses them) and 6-311G family and the Karlsruhe def2 sets.
    static let families: [BasisFamily] = [
        BasisFamily(title: "Pople 6-31G (Cartesian d, 6D)",
                    names: ["6-31G", "6-31G*", "6-31+G", "6-31+G*", "6-31+G**", "6-31++G*", "6-31++G**"]),
        BasisFamily(title: "Pople 6-311G",
                    names: ["6-311G", "6-311G*", "6-311+G", "6-311+G*", "6-311+G**", "6-311++G*", "6-311++G**",
                            "6-311++G(2d,2p)", "6-311++G(3df,3pd)"]),
        BasisFamily(title: "Karlsruhe def2",
                    names: ["def2-SVP", "def2-mTZVP", "def2-mTZVPP", "def2-TZVP", "def2-TZVPP", "def2-TZVPPD",
                            "def2-QZVP", "def2-QZVPP", "def2-QZVPPD"]),
        BasisFamily(title: "Dunning aug-cc-pVXZ",
                    names: ["aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ", "aug-cc-pV5Z"]),
        BasisFamily(title: "Dunning aug-cc-pCVXZ (core-valence)",
                    names: ["aug-cc-pCVDZ", "aug-cc-pCVTZ", "aug-cc-pCVQZ", "aug-cc-pCV5Z"]),
    ]

    /// All 33 basis sets of the library in menu order.
    static let all: [String] = families.flatMap(\.names)

    /// HeH⁺ in the Slater-type basis of Szabo and Ostlund, expanded in 1, 2 or 3 Gaussians.
    static let szaboOstlund = ["STO-1G", "STO-2G", "STO-3G"]

    /// 14.0: the 6-31G family, whose d functions are six Cartesian functions (6D).
    static func isCartesianD(_ name: String) -> Bool {
        families[0].names.contains { $0.lowercased() == name.lowercased() }
    }

    /// The basis sets of a molecule's menu (HeH⁺: the STO-nG sets first).
    static func choices(for molecule: MoleculeChoice) -> [String] {
        molecule == .heh ? szaboOstlund + all : all
    }

    static func isSzaboOstlund(_ name: String) -> Bool {
        name.uppercased().hasPrefix("STO")
    }

    /// Folder with the .gbs files inside the app bundle ("" if missing; the backend then reports
    /// that it cannot read the file).
    static var directory: String {
        Bundle.main.resourceURL?.appendingPathComponent("basis", isDirectory: true).path ?? ""
    }
}

nonisolated enum MachineLimits {
    static var physicalMemoryMB: Double {
        Double(ProcessInfo.processInfo.physicalMemory) / 1_048_576
    }

    /// Two-electron integrals are kept in memory up to a quarter of the RAM (at most 8 GB).
    static var integralMemoryMB: Double {
        min(max(physicalMemoryMB / 4, 512), 8192)
    }

    /// Memory limit of Full CI and coupled cluster: 60 % of the RAM (at least 2 GB, at most
    /// 48 GB).
    static var fciMemoryMB: Double {
        min(max(physicalMemoryMB * 0.6, 2048), 49152)
    }
}
