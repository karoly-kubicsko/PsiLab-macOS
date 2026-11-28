//
//  AppModel.swift
//  MacQC
//
//  Observable application state: sidebar inputs (molecule, job type, method, basis set,
//  geometry, thermochemistry, optimisation, relaxed-scan and transition-state settings), size
//  estimate, the running job with its streamed output file, last result, potential-curve scan,
//  and the menu actions (export text / JSON, copy summary).
//

import SwiftUI
import AppKit
import UniformTypeIdentifiers

/// The six tabs of the detail area (ΨLab 5.0). Methods, Frequencies and Scan have pages.
enum DetailTab: Hashable {
    case summary, methods, optimization, frequencies, potentialCurve, output
}

/// Pages of the Methods tab (8.0: the molecular orbitals).
enum MethodsPage: String, CaseIterable, Identifiable {
    case hartreeFock = "Hartree–Fock"
    case coupledCluster = "Coupled cluster"
    case fullCI = "Full-CI"
    case orbitals = "Orbitals"
    var id: String { rawValue }
}

/// Pages of the Scan tab (ΨLab 7.0): the relaxed scan of a job, and the potential curve of
/// single-point energies along the bond length.
enum ScanPage: String, CaseIterable, Identifiable {
    case relaxed = "Relaxed scan"
    case potentialCurve = "Potential curve"
    var id: String { rawValue }
}

/// The inputs of a geometry (bohr and degrees): the start of the next job when copied from an
/// optimised or scanned geometry.
struct GeometryInputs: Equatable {
    var R: Double
    var angle: Double?
    var ROO: Double?
    var dihedral: Double?
    /// The geometry has the symmetry of the inputs (otherwise R and the angle are means).
    var symmetric: Bool
}

/// Pages of the Frequencies tab.
enum VibrationPage: String, CaseIterable, Identifiable {
    case frequencies = "Frequencies"
    case thermochemistry = "Thermochemistry"
    var id: String { rawValue }
}

/// One point of the potential-curve scan: the (R)HF energy and that of the correlated method
/// (Full CI or coupled cluster; missing in an HF-only scan).
struct CurvePoint: Identifiable, Sendable {
    var R: Double
    var hf: Double
    var correlated: Double?
    var id: Double { R }
}

/// Flattened (method, R, energy) sample for Swift Charts.
struct CurveSample: Identifiable {
    let id = UUID()
    var method: String
    var R: Double
    var energy: Double
}

/// Minimum of one curve (refined by a parabola through the three lowest grid points).
struct CurveMinimum: Identifiable {
    var method: String
    var R: Double
    var energy: Double
    var refined: Bool
    var id: String { method }
}

@MainActor
@Observable
final class AppModel {

    // MARK: Kind of calculation (20.0)

    /// Electronic structure or nuclear motion: the switch at the top of the sidebar; the nuclear-
    /// motion calculation has its own sidebar, tabs and model (NuclearModel).
    var calculation: CalculationKind = .electronic
    let nuclear: NuclearModel

    // MARK: Sidebar inputs

    static let defaultR = 1.4632
    static let defaultZeta1 = 2.0925
    static let defaultZeta2 = 1.24
    static let defaultGaussians = 3

    var molecule: MoleculeChoice = .heh

    // HeH+: Szabo–Ostlund (STO-1G, STO-2G, STO-3G) or (7.0) a Gaussian basis set of the library
    var R: Double = AppModel.defaultR             // bohr (the backend's unit), also of the general HeH+
    var zeta1: Double = AppModel.defaultZeta1   // He
    var zeta2: Double = AppModel.defaultZeta2   // H
    var hehBasis: String = "STO-\(AppModel.defaultGaussians)G"
    var exactErf: Bool = false

    /// Number of Gaussians of the Szabo–Ostlund basis (reads and writes `hehBasis`).
    var nGaussians: Int {
        get { Int(hehBasis.dropFirst(4).prefix(1)) ?? AppModel.defaultGaussians }
        set { hehBasis = "STO-\(newValue)G" }
    }

    // General path: H2, LiH, H2O, H2O2, NH3 and HeH+ in a Gaussian basis set (14.0: Pople, def2, Dunning)
    var basis: String = "aug-cc-pVDZ"
    var method: MethodChoice = .hf
    var bondLengths: [MoleculeChoice: Double] = [:]   // bohr; missing = experimental value
    var angles: [MoleculeChoice: Double] = [:]        // degrees; missing = default angle
    var ooBondLength: Double = MoleculeChoice.defaultPeroxideOO       // bohr, H2O2
    var dihedral: Double = MoleculeChoice.defaultPeroxideDihedral     // degrees, H2O2

    // ΨLab 9.0: custom molecules (structure files), kept in Application Support/ΨLab/Molecules.json
    var savedMolecules: [SavedMolecule] = MoleculeStore.load()
    var customID: UUID?                    // the selected one (molecule == .custom)
    var structureError: String?            // why the last structure file could not be read
    var showStructureError = false
    var savedMessage: String?              // after keeping an optimised geometry

    // Relaxed scan and transition-state search (7.0)
    var scans: [MoleculeChoice: ScanSettings] = [:]   // missing = the molecule's default scan
    var tsHessianMode: HessianMode = .calcFC          // 15.0: CalcFC, CalcAll or RecalcFC=N
    var eigenTest = true                              // exactly one negative eigenvalue
    var analyticHessian = true                        // 10.0: Full CI, CCSD, CCSD(T) Hessians analytic
    var analyticInitialHessian = true                 // 10.0: minimisations start from the analytic Hessian
    // 15.0: computed Hessians of minimisations (Opt=CalcFC, CalcAll, RecalcFC=N), N, and the step
    var optHessianMode: HessianMode = .automatic
    var recalcFC: Int = 3
    var gediis = false                                // hybrid GEDIIS/RFO instead of RFO

    /// 7.0–14.0: the Hessian of the transition-state search at every step.
    var tsCalcAll: Bool {
        get { tsHessianMode == .calcAll }
        set { tsHessianMode = newValue ? .calcAll : .calcFC }
    }

    // 15.0: the selected atom or molecule of the menu (molecule == .preset), the charge and the
    // multiplicity chosen in the sidebar for each molecule (missing: the molecule's own), the
    // reference, and HeH⁺ in STO-nG as the Szabo–Ostlund calculation (or a general one)
    var presetKey: String = "O2"
    /// 16.0: the geometry parameters of the molecules of the menu (Å, degrees; missing: the
    /// backend's geometry), see PresetGeometry
    var presetValues: [String: [Double]] = [:]
    var charges: [String: Int] = [:]
    var multiplicities: [String: Int] = [:]
    var reference: ReferenceChoice = .automatic

    /// 15.0: the reference the calculation uses (the choice, or Gaussian's default: RHF for
    /// closed-shell singlets, UHF for open shells, ROHF for open-shell Full CI).
    var effectiveReference: ReferenceChoice {
        if reference != .automatic { return reference }
        if multiplicity == 1 { return .rhf }
        return method == .fci ? .rohf : .uhf
    }

    /// 15.0: an analytic Hessian beyond the SCF for the method on this reference. UHF and ROHF:
    /// every correlated method (CCSDT and CCSDTQ included) and Full CI; RHF: all but CCSDT, CCSDTQ.
    var hasAnalyticCorrelatedHessian: Bool {
        if method == .hf { return false }
        return effectiveReference == .rhf ? method.hasAnalyticCorrelatedHessian : true
    }
    /// 17.0: off by default (HeH⁺ in STO-nG is then a general molecule); on with the switch in the
    /// sidebar or the launch argument -psilabSzaboOstlund YES.
    var hehSzaboOstlund = false

    // Job (ΨLab 5.0)
    var job: JobChoice = .energy
    var spMethod: MethodChoice? = nil      // single point of job 5; nil = the job's method
    var spBasis: String? = nil             // nil = the job's basis set
    var temperature: Double = 298.15       // K
    var pressure: Double = 101325          // Pa
    var tightConvergence = false
    var maxSteps = 50

    /// The bond length in angstrom; reads and writes `R` (bohr).
    var rAngstrom: Double {
        get { R * Units.bohrToAngstrom }
        set { R = newValue / Units.bohrToAngstrom }
    }

    /// Bond length (bohr) of the selected molecule (HeH⁺: `R`, in both basis-set families).
    var bondLength: Double {
        get { molecule == .heh ? R : (bondLengths[molecule] ?? molecule.defaultBondLength) }
        set {
            if molecule == .heh { R = newValue } else { bondLengths[molecule] = newValue }
        }
    }

    var bondLengthAngstrom: Double {
        get { bondLength * Units.bohrToAngstrom }
        set { bondLength = newValue / Units.bohrToAngstrom }
    }

    /// The bond angle of the selected molecule (H–O–H, O–O–H or H–N–H) in degrees.
    var angle: Double {
        get { angles[molecule] ?? molecule.defaultAngle }
        set { angles[molecule] = newValue }
    }

    var ooBondLengthAngstrom: Double {
        get { ooBondLength * Units.bohrToAngstrom }
        set { ooBondLength = newValue / Units.bohrToAngstrom }
    }

    /// The scan of the selected molecule.
    var scan: ScanSettings {
        get { scans[molecule] ?? molecule.defaultScan }
        set { scans[molecule] = newValue }
    }

    /// The basis set of the selected molecule's menu (HeH⁺: `hehBasis`).
    var currentBasis: String {
        get { molecule == .heh ? hehBasis : basis }
        set {
            if molecule == .heh { hehBasis = newValue } else { basis = newValue }
        }
    }

    /// The general path: every molecule except HeH⁺ in STO-nG (Szabo–Ostlund; 15.0: unless that
    /// option is off).
    var isGeneral: Bool { molecule != .heh || !BasisSets.isSzaboOstlund(hehBasis) || !hehSzaboOstlund }

    // MARK: Charge and multiplicity (15.0)

    /// The key of the selected molecule for its charge and multiplicity ("O2", "H2O", a custom
    /// molecule's identifier).
    var spinKey: String {
        switch molecule {
        case .preset: return "preset:" + presetKey
        case .custom: return "custom:" + (customID?.uuidString ?? "")
        default: return molecule.rawValue
        }
    }

    /// The molecule's own charge (HeH⁺, H₂⁺, H₃⁺, H₅⁺: +1; custom: from its file).
    var defaultCharge: Int {
        switch molecule {
        case .heh: return 1
        case .preset: return selectedPreset?.charge ?? 0
        case .custom: return customMolecule?.charge ?? 0
        default: return 0
        }
    }

    /// The molecule's own multiplicity (the ground states of the atoms and of O₂, NO, NO₂ …).
    var defaultMultiplicity: Int {
        switch molecule {
        case .preset: return selectedPreset?.multiplicity ?? 1
        case .custom: return customMolecule?.multiplicity ?? 1
        default: return 1
        }
    }

    /// Nuclear charges of the selected molecule.
    var nuclearCharge: Int {
        switch molecule {
        case .preset: return (selectedPreset?.symbols ?? []).reduce(0) { $0 + Elements.atomicNumber($1) }
        case .custom: return (customMolecule?.atoms ?? []).reduce(0) { $0 + $1.Z }
        default: return molecule.atomSymbols.reduce(0) { $0 + Elements.atomicNumber($1) }
        }
    }

    var charge: Int {
        get { molecule == .custom ? defaultCharge : (charges[spinKey] ?? defaultCharge) }
        set {
            if molecule == .custom { customCharge = newValue } else { charges[spinKey] = newValue == defaultCharge ? nil : newValue }
        }
    }

    var multiplicity: Int {
        get { molecule == .custom ? defaultMultiplicity : (multiplicities[spinKey] ?? defaultMultiplicity) }
        set {
            if molecule == .custom {
                updateCustom { $0.multiplicity = newValue }
            } else {
                multiplicities[spinKey] = newValue == defaultMultiplicity ? nil : newValue
            }
        }
    }

    var electronCount: Int { nuclearCharge - charge }

    /// Why charge and multiplicity do not fit together (nil if they do).
    var spinProblem: String? {
        let n = electronCount
        if n < 1 { return "No electrons: the charge is too high." }
        let unpaired = multiplicity - 1
        if unpaired > n { return "Multiplicity \(multiplicity) needs at least \(unpaired) electrons; there are \(n)." }
        if (n - unpaired) % 2 != 0 {
            return "\(n) electrons and multiplicity \(multiplicity) do not fit: an \(n % 2 == 0 ? "even" : "odd") number of electrons needs an \(n % 2 == 0 ? "odd" : "even") multiplicity."
        }
        return nil
    }

    /// Name of a multiplicity ("singlet", "doublet", …).
    static func spinName(_ m: Int) -> String {
        let names = ["singlet", "doublet", "triplet", "quartet", "quintet", "sextet", "septet", "octet", "nonet", "decet"]
        return m >= 1 && m <= names.count ? names[m - 1] : "multiplicity \(m)"
    }

    // MARK: Size estimate (general molecules)

    /// Output of an `estimateOnly` call for the current inputs (nil for HeH+).
    var estimate: MacQCOutput?
    /// 16.0: the `estimateKey` the displayed estimate was made for (it is out of date while a new
    /// one is computed)
    var estimatedKey = ""

    /// Changes whenever an input that affects the estimate changes.
    var estimateKey: String {
        guard isGeneral else { return "HeH+" }
        let s = scan
        let custom = molecule == .custom
            ? (customMolecule.map { "\($0.id)|\($0.charge)|\($0.multiplicity)|\($0.atoms.hashValue)" } ?? "none") : ""
        return "\(molecule.rawValue)|\(custom)|\(currentBasis)|\(method.rawValue)|\(bondLength)|\(angle)|\(ooBondLength)|\(dihedral)|"
            + "\(job.rawValue)|\(spMethod?.rawValue ?? "")|\(spBasis ?? "")|\(temperature)|\(pressure)|"
            + "\(tightConvergence)|\(maxSteps)|\(s.line)|\(s.fromInput)|\(s.start)|\(tsHessianMode.rawValue)|\(eigenTest)|\(analyticHessian)|\(analyticInitialHessian)|"
            + "\(presetKey)|\(presetValues[presetKey] ?? [])|\(charge)|\(multiplicity)|\(reference.rawValue)|\(hehSzaboOstlund)|\(optHessianMode.rawValue)|\(recalcFC)|\(gediis)"
    }

    /// The backend refuses the calculation (reason), or the inputs are invalid.
    var refusal: String? {
        guard isGeneral, let e = estimate else { return nil }
        if let r = e.job?.estimate.refusal { return r }
        if let r = e.refusal { return r }
        if e.estimate == nil, let message = e.error { return message }
        return nil
    }

    /// Refusal of a single energy with these settings (potential curves run energies only).
    var energyRefusal: String? {
        guard isGeneral, let e = estimate else { return nil }
        if let r = e.refusal { return r }
        if e.estimate == nil, let message = e.error { return message }
        return nil
    }

    /// The whole job is estimated to take longer than 30 s.
    var isLongCalculation: Bool {
        estimate?.job?.estimate.longCalculation ?? estimate?.estimate?.longCalculation ?? false
    }

    /// Estimated time of the whole job (seconds).
    var estimatedSeconds: Double? {
        estimate?.job?.estimate.totalSeconds ?? estimate?.estimate?.totalSeconds
    }

    /// The single-point method and basis set of jobs 5 and 9 (the job's own when not chosen).
    var effectiveSPMethod: MethodChoice { spMethod ?? method }
    var effectiveSPBasis: String { spBasis ?? currentBasis }

    /// Problem of the scan settings that the app reports before asking the backend (atoms), or nil.
    var scanProblem: String? {
        guard job == .scan else { return nil }
        return scan.atomProblem(atomCount: atomCount)
    }

    // MARK: Results and the running job

    var selectedTab: DetailTab = .summary
    var methodsPage: MethodsPage = .hartreeFock
    var scanPage: ScanPage = .potentialCurve
    var vibrationPage: VibrationPage = .frequencies
    var isRunning = false
    var result: MacQCOutput? {
        didSet { resultSerial += 1 }
    }
    /// Changes with every new result (the orbital page starts afresh).
    var resultSerial = 0
    /// 13.0: the output file (or JSON export) whose results are shown, nil for a calculation of
    /// this session.
    var openedFile: OpenedOutputFile?
    /// An output file was requested at launch (-psilabOpenOutput): no calculation at launch.
    var openedAtLaunch = false
    var errorMessage: String?
    var progress = MacQCProgress()
    /// The backend's description of the current step ("Optimisation step 4 · E = … · max force …").
    var progressLine = ""
    var runStarted: Date?
    var showLongRunConfirmation = false

    /// The output file of the latest run (streamed while it runs; also after a cancelled or
    /// failed run). `outputGeneration` changes when the text is replaced rather than extended.
    var outputText = ""
    var outputGeneration = 0

    /// Starts the job at launch even when it is estimated to be long (launch argument
    /// -psilabAutoRun YES, used by the UI smoke test).
    var autoRunLongJobs = false

    @ObservationIgnored private var runningJob: MacQCJob?

    // MARK: Potential curve

    /// Step choices in angstrom (converted to bohr for the scan grid).
    static let scanStepChoices: [Double] = [0.025, 0.05, 0.075, 0.1, 0.15]

    var scanStep: Double = 0.05   // angstrom
    var isScanning = false
    var scanProgress: Double = 0
    var scanMessage: String?
    var curvePoints: [CurvePoint] = []
    var curveInputDescription: String?   // parameters the scan was computed with
    var curveHFLabel = "HF"               // "HF" for HeH+, "RHF" for the other molecules
    var curveCorrelatedLabel = "FCI"      // "FCI", "CCSD", "CCSD(T)", "CCSDT" or "CCSDTQ"
    var showLongScanConfirmation = false

    @ObservationIgnored private var scanTask: Task<Void, Never>?
    @ObservationIgnored private var scanJob: MacQCJob?

    // MARK: Start-up

    /// 23.0: one model per window and per tab of a window (MacQCApp). The first model takes the
    /// launch arguments of ΨLab (the UI test's -psilab…); every later window or tab opens with the
    /// startup examples (H₂ in the electronic-structure interface, H₂¹⁶O with J = 10 in the
    /// nuclear-motion interface), as the app does without launch arguments.
    private static var launchArgumentsTaken = false

    /// 23.0: the models of all open windows and tabs (weak), so that the kept molecules (one file,
    /// Molecules.json) are the same in every molecule menu: a model that saves the list hands it to
    /// the others, so that a later save in another window does not drop it.
    private struct WeakModel {
        weak var model: AppModel?
    }
    private static var liveModels: [WeakModel] = []

    /// 23.0: the model of the menu commands while no ΨLab window is focused (they are disabled
    /// then; it is never shown and reads no example).
    static let placeholder = AppModel(placeholder: true)

    convenience init() {
        self.init(placeholder: false)
    }

    init(placeholder: Bool) {
        let first = !placeholder && !AppModel.launchArgumentsTaken
        if first { AppModel.launchArgumentsTaken = true }
        nuclear = NuclearModel(launchArguments: first)
        guard !placeholder else { return }
        AppModel.liveModels.removeAll { $0.model == nil }
        AppModel.liveModels.append(WeakModel(model: self))
        if first { applyLaunchArguments() }
        // 17.0: without launch arguments of ΨLab the app opens with the results of H₂,
        // Full CI/aug-cc-pCV5Z//Full CI/aug-cc-pCVQZ (opt+freq+sp), read from the bundled output file.
        if !first || !AppModel.hasPsiLabLaunchArguments {
            showStartupExample()
            // 23.0: the nuclear-motion interface with the results of its example (H₂¹⁶O, J = 10);
            // -psilabNuclearExample NO (UI test): the default inputs instead, whose Run is short
            let skip = first && UserDefaults.standard.object(forKey: "psilabNuclearExample") != nil
                && !UserDefaults.standard.bool(forKey: "psilabNuclearExample")
            if !skip { nuclear.showStartupExample() }
        }
    }

    /// 23.0: the kept molecules saved by `source`: the same list in the menus of the other windows.
    static func moleculesSaved(_ list: [SavedMolecule], by source: AppModel) {
        liveModels.removeAll { $0.model == nil }
        for w in liveModels {
            guard let m = w.model, m !== source else { continue }
            m.savedMolecules = list
        }
    }

    /// The launch arguments of ΨLab (the first window only, 23.0).
    private func applyLaunchArguments() {
        // Launch arguments (also used by the UI test): -psilabMolecule H2O -psilabBasis aug-cc-pVTZ
        // -psilabMethod FCI (or HF, CCSD, CCSD(T), CCSDT, CCSDTQ)
        let defaults = UserDefaults.standard
        if let m = defaults.string(forKey: "psilabMolecule"),
           let choice = MoleculeChoice.allCases.first(where: { $0.rawValue.lowercased() == m.lowercased() }) {
            molecule = choice
        } else if let m = defaults.string(forKey: "psilabMolecule"),
                  let p = MoleculePresets.all.first(where: { $0.key.lowercased() == m.lowercased() }) {
            // 15.0: -psilabMolecule O2 (an atom or molecule of the menu)
            presetKey = p.key
            molecule = .preset
        }
        // 15.0: -psilabCharge 1 -psilabMultiplicity 3 -psilabReference ROHF -psilabOptHessian recalcfc
        // -psilabRecalcFC 2 -psilabGEDIIS YES
        if defaults.object(forKey: "psilabCharge") != nil { charge = defaults.integer(forKey: "psilabCharge") }
        if defaults.object(forKey: "psilabMultiplicity") != nil {
            let m = defaults.integer(forKey: "psilabMultiplicity")
            if m >= 1 { multiplicity = m }
        }
        if let r = defaults.string(forKey: "psilabReference"),
           let choice = ReferenceChoice.allCases.first(where: { $0.rawValue.lowercased() == r.lowercased() }) {
            reference = choice
        }
        if let h = defaults.string(forKey: "psilabOptHessian"), let mode = HessianMode(rawValue: h.lowercased()) {
            optHessianMode = mode
        }
        if defaults.object(forKey: "psilabRecalcFC") != nil {
            recalcFC = max(1, defaults.integer(forKey: "psilabRecalcFC"))
        }
        if defaults.object(forKey: "psilabGEDIIS") != nil { gediis = defaults.bool(forKey: "psilabGEDIIS") }
        // ΨLab 9.0: -psilabStructure "<text of an .xyz or .gjf file>" -psilabStructureName co2.xyz
        // (the text itself: the sandboxed app reads files only through the file panel)
        if let text = defaults.string(forKey: "psilabStructure") {
            addStructure(text: text, fileName: defaults.string(forKey: "psilabStructureName") ?? "structure.xyz", select: true)
        } else if molecule == .custom {
            if let last = savedMolecules.last { selectSaved(last.id) } else { molecule = .heh }
        }
        // HeH+: -psilabBasis STO-2G (Szabo–Ostlund) or aug-cc-pVTZ (general path); 14.0: any set, e.g. "6-31G*"
        if let b = defaults.string(forKey: "psilabBasis"),
           let name = BasisSets.choices(for: molecule).first(where: { $0.lowercased() == b.lowercased() }) {
            currentBasis = name
        }
        if let m = defaults.string(forKey: "psilabMethod"),
           let choice = MethodChoice(name: m) {
            method = choice
        }
        // ΨLab 5.0: -psilabJob opt+freq -psilabSPMethod "CCSD(T)" -psilabSPBasis aug-cc-pVTZ
        // -psilabR 1.0 (angstrom) -psilabAngle 110 -psilabTab output -psilabAutoRun YES
        if let j = defaults.string(forKey: "psilabJob"), let choice = JobChoice(name: j) {
            job = choice
        }
        if let m = defaults.string(forKey: "psilabSPMethod"), let choice = MethodChoice(name: m) {
            spMethod = choice
        }
        if let b = defaults.string(forKey: "psilabSPBasis"),
           let name = BasisSets.all.first(where: { $0.lowercased() == b.lowercased() }) {
            spBasis = name
        }
        if defaults.object(forKey: "psilabR") != nil {
            let r = defaults.double(forKey: "psilabR")
            if r > 0 { bondLength = r / Units.bohrToAngstrom }
        }
        if defaults.object(forKey: "psilabAngle") != nil {
            let a = defaults.double(forKey: "psilabAngle")
            if a > 0 { angle = a }
        }
        // ΨLab 7.0: -psilabROO 1.40 (Å) -psilabDihedral 160 -psilabScan "D 3 1 2 4 S 18 10"
        // -psilabScanStart 0 -psilabTSHessian calcall -psilabEigenTest NO -psilabHessian numerical
        // 10.0: -psilabInitialHessian model
        if defaults.object(forKey: "psilabROO") != nil {
            let r = defaults.double(forKey: "psilabROO")
            if r > 0 { ooBondLength = r / Units.bohrToAngstrom }
        }
        if defaults.object(forKey: "psilabDihedral") != nil {
            dihedral = defaults.double(forKey: "psilabDihedral")
        }
        if let line = defaults.string(forKey: "psilabScan"), var s = ScanSettings(line: line) {
            if defaults.object(forKey: "psilabScanStart") != nil {
                s.fromInput = false
                s.start = defaults.double(forKey: "psilabScanStart")
            }
            scan = s
        }
        if let h = defaults.string(forKey: "psilabTSHessian") {
            tsCalcAll = h.lowercased() == "calcall"
        }
        if defaults.object(forKey: "psilabEigenTest") != nil {
            eigenTest = defaults.bool(forKey: "psilabEigenTest")
        }
        if let h = defaults.string(forKey: "psilabHessian") {
            analyticHessian = h.lowercased() != "numerical"
        }
        if let h = defaults.string(forKey: "psilabInitialHessian") {
            analyticInitialHessian = h.lowercased() != "model"
        }
        if let t = defaults.string(forKey: "psilabTab") {
            switch t.lowercased() {
            case "methods": selectedTab = .methods
            case "orbitals":
                selectedTab = .methods
                methodsPage = .orbitals
            case "optimisation", "optimization": selectedTab = .optimization
            case "frequencies": selectedTab = .frequencies
            case "thermochemistry":
                selectedTab = .frequencies
                vibrationPage = .thermochemistry
            case "curve", "potential curve":
                selectedTab = .potentialCurve
                scanPage = .potentialCurve
            case "scan", "relaxed scan":
                selectedTab = .potentialCurve
                scanPage = .relaxed
            case "output": selectedTab = .output
            default: selectedTab = .summary
            }
        }
        autoRunLongJobs = defaults.bool(forKey: "psilabAutoRun")
        // 20.0: -psilabCalculation nuclear opens the nuclear-motion interface
        if defaults.string(forKey: "psilabCalculation") == "nuclear" { calculation = .nuclear }
        // 13.0: -psilabOpenOutput path: open an output file (or a JSON export) at launch; no
        // calculation then, also when the file cannot be read (its error stays visible). The app is
        // sandboxed: the path must lie in its container.
        if let path = defaults.string(forKey: "psilabOpenOutput"), !path.isEmpty {
            openedAtLaunch = true
            loadOutputFile(from: URL(fileURLWithPath: (path as NSString).expandingTildeInPath),
                           showSummary: defaults.string(forKey: "psilabTab") == nil)
        }
        // 17.0: -psilabSzaboOstlund YES: HeH⁺ in STO-nG as the Szabo–Ostlund calculation (off by default)
        if defaults.object(forKey: "psilabSzaboOstlund") != nil {
            hehSzaboOstlund = defaults.bool(forKey: "psilabSzaboOstlund")
        }
    }

    /// 17.0: the app was launched with arguments of ΨLab (-psilabMolecule, -psilabOpenOutput, …;
    /// the UI test uses them); otherwise it opens with the startup example.
    static var hasPsiLabLaunchArguments: Bool {
        // 23.0: -psilabNuclearExample only chooses the inputs of the nuclear-motion interface
        ProcessInfo.processInfo.arguments.dropFirst().contains { $0.hasPrefix("-psilab") && $0 != "-psilabNuclearExample" }
    }

    /// 17.0: the startup example: H₂ optimised and its frequencies computed with Full CI/aug-cc-pCVQZ,
    /// then a Full CI/aug-cc-pCV5Z single point at the optimised geometry (tight convergence, hybrid
    /// GEDIIS/RFO steps). The sidebar holds the inputs of that job, so that Run repeats it; its
    /// results are read from the output file in the app's resources (no calculation at launch).
    func showStartupExample() {
        molecule = .h2
        method = .fci
        basis = "aug-cc-pCVQZ"
        job = .optFreqSP
        spMethod = nil
        spBasis = "aug-cc-pCV5Z"
        tightConvergence = true
        gediis = true
        bondLengths[.h2] = nil   // 0.74144 Å, the starting geometry of the example
        guard let url = Bundle.main.url(forResource: AppModel.startupExampleName, withExtension: "out") else { return }
        openedAtLaunch = true
        loadOutputFile(from: url, example: true)
    }

    /// The output file of the startup example in the app's resources (examples/outputs).
    static let startupExampleName = "h2-fullci-aug-cc-pcv5z-opt-aug-cc-pcvqz-optfreqsp"

    // MARK: Inputs

    func currentInput(printLevel: Int = 2) -> MacQCInput {
        var input = MacQCInput()
        input.printLevel = printLevel
        if isGeneral && molecule == .custom {
            // ΨLab 9.0: the atoms of the selected custom molecule (Å)
            input.molecule = MoleculeChoice.custom.rawValue
            input.R = 0
            if let c = customMolecule {
                input.atoms = c.atoms
                input.charge = c.charge
                input.multiplicity = c.multiplicity
                input.title = c.title
            } else {
                input.atoms = []
            }
        }
        if isGeneral {
            if molecule == .preset {
                // 16.0: an atom or molecule of the menu: the backend's geometry, or the coordinates of
                // the changed geometry parameters
                input.molecule = presetKey
                input.R = 0
                input.atoms = presetAtoms
            } else if molecule != .custom {
                input.molecule = molecule.rawValue
                input.R = bondLength
            }
            // 15.0: charge and multiplicity of every molecule (sent when they differ from the
            // molecule's own; a custom molecule always sends its own), the reference
            if molecule != .custom {
                input.charge = charge == defaultCharge ? nil : charge
                input.multiplicity = multiplicity == defaultMultiplicity ? nil : multiplicity
            }
            input.reference = reference == .automatic ? nil : reference.rawValue
            if molecule == .heh && BasisSets.isSzaboOstlund(hehBasis) { input.szaboOstlund = false }
            input.angle = molecule.angleLabel == nil ? MoleculeChoice.defaultWaterAngle : angle
            if molecule == .h2o2 {
                input.ROO = ooBondLength
                input.dihedral = dihedral
            }
            input.basis = currentBasis
            input.method = method.rawValue
            input.basisDirectory = BasisSets.directory
            input.threads = 0
            input.memoryMB = MachineLimits.integralMemoryMB
            input.fciMemoryMB = MachineLimits.fciMemoryMB
            input.job = job.rawValue
            if job.hasSinglePoint {
                input.spMethod = spMethod?.rawValue ?? ""
                input.spBasis = spBasis ?? ""
            }
            input.temperature = temperature
            input.pressure = pressure
            input.optConvergence = tightConvergence ? "tight" : "normal"
            input.maxSteps = maxSteps
            if job == .scan {
                let s = scan
                input.scan = s.line
                input.scanStart = s.fromInput ? nil : s.start
            }
            input.tsHessian = tsHessianMode == .automatic ? "calcfc" : tsHessianMode.rawValue
            input.optHessian = job.hasTransitionState ? "auto" : optHessianMode.rawValue
            let recalc = job.hasTransitionState ? tsHessianMode == .recalcFC : optHessianMode == .recalcFC
            input.recalcFC = recalc ? recalcFC : 0
            input.optAlgorithm = gediis && !job.hasTransitionState ? "gediis" : "rfo"
            input.eigenTest = eigenTest
            input.hessian = analyticHessian ? "analytic" : "numerical"
            // "auto": the analytic Hessian (HF, Full CI, CCSD, CCSD(T)) where possible, else Lindh's model
            input.initialHessian = hasAnalyticCorrelatedHessian && !analyticInitialHessian ? "model" : "auto"
        } else {
            input.molecule = MoleculeChoice.heh.rawValue
            input.basis = hehBasis   // STO-1G, STO-2G or STO-3G: the Szabo–Ostlund calculation
            input.R = R
            input.zeta1 = zeta1
            input.zeta2 = zeta2
            input.ZA = 2
            input.ZB = 1
            input.nGaussians = nGaussians
            input.exactErf = exactErf
        }
        return input
    }

    func resetToDefaults() {
        R = AppModel.defaultR
        zeta1 = AppModel.defaultZeta1
        zeta2 = AppModel.defaultZeta2
        nGaussians = AppModel.defaultGaussians
        exactErf = false
    }

    /// Equilibrium geometry of the selected molecule (the backend's defaults).
    func resetGeometry() {
        if molecule == .heh {
            R = AppModel.defaultR
        } else {
            bondLengths[molecule] = nil
        }
        angles[molecule] = nil
        if molecule == .h2o2 {
            ooBondLength = MoleculeChoice.defaultPeroxideOO
            dihedral = MoleculeChoice.defaultPeroxideDihedral
        }
    }

    /// The inputs of a geometry of the selected molecule: R (the mean X–H distance), the mean
    /// angle, and for H₂O₂ R(O–O) and the dihedral angle.
    func geometryInputs(_ g: MacQCGeometry, symmetric: Bool) -> GeometryInputs? {
        guard let r = g.R else { return nil }
        return GeometryInputs(R: r, angle: g.angle, ROO: g.ROO, dihedral: g.dihedral, symmetric: symmetric)
    }

    /// The optimised geometry (minimum or transition state) of the displayed result, if it is for
    /// the selected molecule and converged.
    var optimizedGeometry: GeometryInputs? {
        guard isGeneral, molecule != .custom, molecule != .preset, let out = result, out.input.molecule == molecule.rawValue, out.input.isGeneral,
              let opt = out.optimization, opt.converged else { return nil }
        let symmetric = opt.symmetricEnd || (out.input.molecule != "H2O" && out.input.molecule != "H2O2" && out.input.molecule != "NH3")
        return geometryInputs(opt.finalGeometry, symmetric: symmetric)
    }

    /// Copying the geometry changes the inputs.
    func differsFromInputs(_ g: GeometryInputs) -> Bool {
        if abs(g.R - bondLength) > 1e-9 { return true }
        if molecule.angleLabel != nil, let a = g.angle, abs(a - angle) > 1e-9 { return true }
        if molecule == .h2o2 {
            if let r = g.ROO, abs(r - ooBondLength) > 1e-9 { return true }
            if let d = g.dihedral, abs(d - dihedral) > 1e-9 { return true }
        }
        return false
    }

    /// The optimised geometry differs from the inputs (so copying it changes something).
    var canUseOptimizedGeometry: Bool {
        guard let g = optimizedGeometry else { return false }
        return differsFromInputs(g)
    }

    /// Copies a geometry into the geometry inputs (the start of the next job).
    func useGeometry(_ g: GeometryInputs) {
        bondLength = g.R
        if molecule.angleLabel != nil, let a = g.angle {
            angle = min(max(a, molecule.angleRange.lowerBound), molecule.angleRange.upperBound)
        }
        if molecule == .h2o2 {
            if let r = g.ROO { ooBondLength = r }
            if let d = g.dihedral { dihedral = d }
        }
    }

    /// Copies the optimised geometry into the geometry inputs (the start of the next job).
    func useOptimizedGeometry() {
        guard let g = optimizedGeometry else { return }
        useGeometry(g)
    }

    /// A point of the displayed relaxed scan as the start of a transition-state search
    /// (Gaussian's usual route: scan, then Opt=TS from the highest point): copies its geometry
    /// into the inputs and selects the job TS optimisation + frequencies.
    func startTransitionStateSearch(from point: MacQCScanPoint) {
        if molecule == .custom {
            startCustomTransitionStateSearch(from: point)
            return
        }
        guard let out = result, out.input.molecule == molecule.rawValue,
              let g = geometryInputs(point.geometry, symmetric: true) else { return }
        useGeometry(g)
        job = .tsFreq
    }

    /// The displayed scan belongs to the selected molecule (its points can be copied).
    var scanMatchesMolecule: Bool {
        guard let out = result, out.scan != nil else { return false }
        if molecule == .custom { return resultMatchesCustom }
        return out.input.molecule == molecule.rawValue && out.input.isGeneral && isGeneral
    }

    // MARK: Estimate

    /// Asks the backend for the size, time and memory of the current input (no calculation).
    func refreshEstimate() async {
        guard isGeneral else {
            estimate = nil
            return
        }
        let key = estimateKey
        var input = currentInput(printLevel: 0)
        input.estimateOnly = true
        do {
            let output = try await MacQCEngine.run(input)
            if key == estimateKey {   // inputs unchanged meanwhile
                estimate = output
                estimatedKey = key
            }
        } catch {
            if key == estimateKey {
                estimate = nil
                estimatedKey = key
            }
        }
    }

    /// Text of the estimate for a confirmation dialog.
    var longRunMessage: String {
        if let j = estimate?.job, !j.estimate.warnings.isEmpty {
            return j.estimate.warnings.joined(separator: "\n")
        }
        return estimate?.warnings?.joined(separator: "\n") ?? "This calculation may take a long time."
    }

    // MARK: Single calculation

    /// Run button: asks for confirmation first when the estimate warns of a long calculation.
    func requestRun() {
        guard !isRunning else { return }
        if isGeneral && (refusal != nil || scanProblem != nil) { return }
        if isGeneral && isLongCalculation {
            showLongRunConfirmation = true
        } else {
            Task { await run() }
        }
    }

    func run() async {
        guard !isRunning else { return }
        isRunning = true
        errorMessage = nil
        savedMessage = nil
        progress = MacQCProgress()
        progressLine = ""
        runStarted = Date()
        let input = currentInput(printLevel: 2)
        if input.isGeneral {
            outputText = ""
            outputGeneration += 1
        }
        let handle = MacQCJob()
        runningJob = handle
        // Poll the backend's progress and the new part of the output file while the job runs.
        let poller = Task { [weak self] in
            var offset = 0
            while !Task.isCancelled {
                let p = handle.progress()
                let line = handle.progressText()
                let piece = handle.output(from: offset)
                offset = max(offset, piece.total)
                guard let self else { return }
                self.progress = p
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
            let output = try await MacQCEngine.run(input, job: handle)
            // The complete output file (also of a cancelled or failed job).
            if input.isGeneral {
                let full = output.output ?? handle.output(from: 0).text
                if !full.isEmpty && full != outputText {
                    outputText = full
                    outputGeneration += 1
                }
            }
            if output.error == "calculation cancelled" {
                errorMessage = "Calculation cancelled."   // keep the previous result visible
                return
            }
            result = output
            openedFile = nil   // 13.0: a result of this session (a cancelled or failed run keeps the file's)
            errorMessage = output.error
            // UI smoke test (-psilabSaveOptimized YES): keep the optimised geometry
            if UserDefaults.standard.bool(forKey: "psilabSaveOptimized") && optimizedMolecule != nil {
                saveOptimizedMolecule(select: false)
            }
            if output.scan != nil {
                scanPage = .relaxed
            }
        } catch {
            // Keep the previous result visible; show why the new run failed.
            errorMessage = error.localizedDescription
        }
    }

    /// The progress line under the Run button.
    var progressDescription: String {
        progressLine.isEmpty ? progress.text : progressLine
    }

    func cancelRun() {
        runningJob?.cancel()
    }

    /// Elapsed time of the running calculation, e.g. "1:05".
    func elapsedText(at now: Date) -> String {
        guard let start = runStarted else { return "" }
        let s = Int(now.timeIntervalSince(start))
        return String(format: "%d:%02d", s / 60, s % 60)
    }

    // MARK: Potential-curve scan

    var scanRange: ClosedRange<Double> { molecule.scanRange }

    /// Grid of R values over the molecule's scan range (inclusive) with the selected step.
    func scanGrid() -> [Double] {
        let step = scanStep / Units.bohrToAngstrom   // grid is in bohr
        let start = scanRange.lowerBound, end = scanRange.upperBound
        guard step > 0 else { return [] }
        var values: [Double] = []
        let count = Int(((end - start) / step + 1.0e-9).rounded(.down))
        for i in 0...count {
            values.append(start + Double(i) * step)
        }
        if let last = values.last, last < end - 1.0e-9 {
            values.append(end)
        }
        return values
    }

    /// Estimated time of a scan (general molecules), from the single-point estimate plus
    /// about 0.2 s per point for setting up each calculation (basis files, start-up).
    var scanSeconds: Double {
        guard isGeneral, let e = estimate?.estimate else { return 0 }
        return Double(scanGrid().count) * (e.totalSeconds + 0.2)
    }

    /// Run scan button: confirmation first when the scan would take more than 30 s.
    func requestScan() {
        guard !isScanning else { return }
        if let reason = curveUnavailableReason {
            scanMessage = reason
            return
        }
        if isGeneral && energyRefusal != nil {
            scanMessage = energyRefusal
            return
        }
        if isGeneral && scanSeconds > 30 {
            showLongScanConfirmation = true
        } else {
            startScan()
        }
    }

    func startScan() {
        guard !isScanning else { return }
        scanTask = Task { await self.runScan() }
    }

    func cancelScan() {
        scanTask?.cancel()
        scanJob?.cancel()
    }

    private func runScan() async {
        isScanning = true
        scanProgress = 0
        scanMessage = nil
        defer {
            isScanning = false
            scanTask = nil
            scanJob = nil
        }

        let grid = scanGrid()
        var points: [CurvePoint] = []
        var failures: [String] = []
        var template = currentInput(printLevel: 0)
        template.job = "energy"   // a curve is a series of single-point energies, whatever the job
        template.scan = ""
        template.scanStart = nil
        let general = template.isGeneral

        for (index, r) in grid.enumerated() {
            if Task.isCancelled {
                scanMessage = "Scan cancelled."
                break
            }
            var input = template
            input.R = r
            let job = MacQCJob()
            scanJob = job
            do {
                let output = try await MacQCEngine.run(input, job: job)
                if output.error == "calculation cancelled" {
                    scanMessage = "Scan cancelled."
                    break
                }
                if general, let rhf = output.rhf, output.error == nil {
                    points.append(CurvePoint(R: r, hf: rhf.totalEnergy, correlated: output.correlatedEnergy))
                } else if !general, let hf = output.hf, let fci = output.fci, output.error == nil {
                    points.append(CurvePoint(R: r, hf: hf.totalEnergy, correlated: fci.totalEnergy))
                } else {
                    failures.append(String(format: "R = %.3f: %@", r, output.error ?? "no result"))
                }
            } catch {
                failures.append(String(format: "R = %.3f: %@", r, error.localizedDescription))
            }
            scanProgress = Double(index + 1) / Double(grid.count)
        }

        curvePoints = points
        curveHFLabel = general ? "RHF" : "HF"
        curveCorrelatedLabel = general ? (template.method == "FCI" ? "FCI" : template.method) : "FCI"
        if general {
            let choice = MoleculeChoice(rawValue: template.molecule)
            var fixed = ""
            if let label = choice?.angleLabel {
                fixed = String(format: " · %@ %.2f°", label.lowercased(), template.angle)
            }
            if template.molecule == "H2O2", let roo = template.ROO, let d = template.dihedral {
                fixed += String(format: ", R(O–O) %.4f Å, dihedral %.2f°", roo * Units.bohrToAngstrom, d)
            }
            curveInputDescription = String(
                format: "%@ · %@ · %@%@", choice?.title ?? template.molecule, template.basis,
                template.method == "RHF" ? "RHF" : "RHF and " + MethodChoice.title(of: template.method), fixed)
        } else {
            curveInputDescription = String(
                format: "STO-%ldG, ζ(He) = %.4f, ζ(H) = %.4f, %@ erf",
                template.nGaussians, template.zeta1, template.zeta2,
                template.exactErf ? "exact" : "A&S")
        }
        if !failures.isEmpty {
            scanMessage = "\(failures.count) point(s) failed: " + failures.prefix(3).joined(separator: "; ")
        }
    }

    var curveSamples: [CurveSample] {
        var samples: [CurveSample] = []
        for p in curvePoints {
            samples.append(CurveSample(method: curveHFLabel, R: p.R, energy: p.hf))
        }
        for p in curvePoints {
            if let e = p.correlated {
                samples.append(CurveSample(method: curveCorrelatedLabel, R: p.R, energy: e))
            }
        }
        return samples
    }

    var curveMinima: [CurveMinimum] {
        var minima: [CurveMinimum] = []
        if let m = AppModel.minimum(of: curvePoints.map { ($0.R, $0.hf) }, method: curveHFLabel) {
            minima.append(m)
        }
        let correlatedSamples = curvePoints.compactMap { p in p.correlated.map { (p.R, $0) } }
        if let m = AppModel.minimum(of: correlatedSamples, method: curveCorrelatedLabel) {
            minima.append(m)
        }
        return minima
    }

    /// Lowest grid point; if it is interior and the grid is equally spaced, refined by the
    /// vertex of the parabola through it and its two neighbours.
    static func minimum(of samples: [(Double, Double)], method: String) -> CurveMinimum? {
        guard !samples.isEmpty else { return nil }
        var best = 0
        for i in samples.indices where samples[i].1 < samples[best].1 {
            best = i
        }
        var result = CurveMinimum(method: method, R: samples[best].0, energy: samples[best].1,
                                  refined: false)
        if best > 0 && best < samples.count - 1 {
            let (r0, e0) = samples[best - 1]
            let (r1, e1) = samples[best]
            let (r2, e2) = samples[best + 1]
            let h1 = r1 - r0
            let h2 = r2 - r1
            let curvature = e0 - 2.0 * e1 + e2
            if abs(h1 - h2) < 1.0e-9 && curvature > 0 {
                let shift = -0.5 * h1 * (e2 - e0) / curvature
                let energy = e1 - (e2 - e0) * (e2 - e0) / (8.0 * curvature)
                if abs(shift) <= h1 && energy.isFinite {
                    result.R = r1 + shift
                    result.energy = energy
                    result.refined = true
                }
            }
        }
        return result
    }

    // MARK: Summary text / export

    func summaryText() -> String? {
        guard let out = result else { return nil }
        if out.isGeneral { return generalSummaryText(out) }
        guard let hf = out.hf else { return nil }
        var lines: [String] = []
        lines.append("ΨLab – HeH+ (STO-\(out.input.nGaussians)G)")
        lines.append(String(format: "R = %@, zeta(He) = %.4f, zeta(H) = %.4f, %@ erf",
                            Fmt.lengthBoth(out.input.R), out.input.zeta1, out.input.zeta2,
                            out.input.exactErf ? "exact" : "Abramowitz–Stegun"))
        lines.append("")
        lines.append("Hartree–Fock")
        lines.append("  converged:            \(hf.converged ? "yes" : "no") (\(hf.iterations) iterations)")
        lines.append("  electronic energy:    " + Fmt.energyWithUnit(hf.electronicEnergy))
        lines.append("  nuclear repulsion:    " + Fmt.energyWithUnit(hf.nuclearRepulsion))
        lines.append("  total energy:         " + Fmt.energyWithUnit(hf.totalEnergy))
        if hf.orbitalEnergies.count >= 2 {
            lines.append("  orbital energies:     " + hf.orbitalEnergies.map { Fmt.energy($0) }.joined(separator: ", ") + " E_h")
        }
        if let fci = out.fci {
            lines.append("")
            lines.append("Full-CI")
            lines.append("  total energy:         " + Fmt.energyWithUnit(fci.totalEnergy))
            lines.append("  correlation energy:   " + Fmt.energyWithUnit(fci.correlationEnergy))
            for (i, s) in fci.states.enumerated() {
                lines.append(String(format: "  state %ld: %@ E_h  <S^2> = %.4f  %@",
                                    i + 1, Fmt.energy(s.energy), s.S2, s.spin))
            }
            for w in fci.singletWeights {
                lines.append(String(format: "  weight %@: %.6f", w.configuration, w.weight))
            }
        }
        if let message = out.error {
            lines.append("")
            lines.append("Error: \(message)")
        }
        return lines.joined(separator: "\n")
    }

    private func generalSummaryText(_ out: MacQCOutput) -> String? {
        guard let mol = out.molecule, let rhf = out.rhf else { return nil }
        var lines: [String] = []
        lines.append("ΨLab – \(out.moleculeTitle), \(out.input.method)/\(out.basis?.name ?? out.input.basis)")
        lines.append(Fmt.moleculeGeometry(mol, plain: true))
        if let b = out.basis {
            let atoms = Fmt.uniqueAtoms(b.atoms).map { "\($0.symbol) \($0.contraction)" }.joined(separator: ", ")
            lines.append("\(b.functions) basis functions (\(atoms))")
        }
        for note in out.notes ?? [] {
            lines.append("Note: \(note)")
        }
        lines.append("")
        lines.append("Hartree–Fock (\(out.referenceLabel))")
        lines.append("  converged:            \(rhf.converged ? "yes" : "no") (\(rhf.iterations) iterations)")
        lines.append("  total energy:         " + Fmt.energy10(rhf.totalEnergy) + " E_h")
        lines.append("  electronic energy:    " + Fmt.energy10(rhf.electronicEnergy) + " E_h")
        lines.append("  nuclear repulsion:    " + Fmt.energy10(rhf.nuclearRepulsion) + " E_h")
        if let homo = rhf.homo, let lumo = rhf.lumo {
            lines.append(String(format: "  HOMO, LUMO:           %.8f, %.8f E_h", homo, lumo))
        }
        if let f = out.fullCI {
            lines.append("")
            lines.append("Full-CI")
            lines.append(String(format: "  determinants:         %@ (%ld orbitals, %ld electrons)",
                                Fmt.count(f.determinants), f.orbitals, f.alphaElectrons + f.betaElectrons))
            lines.append("  converged:            \(f.converged ? "yes" : "no") (\(f.iterations) Davidson iterations)")
            lines.append("  total energy:         " + Fmt.energy10(f.totalEnergy) + " E_h")
            lines.append("  correlation energy:   " + Fmt.energy10(f.correlationEnergy) + " E_h")
            lines.append(String(format: "  <S^2>:                %.6f", f.S2))
            lines.append(String(format: "  RHF determinant:      |C_0|^2 = %.6f", f.hfWeight))
        }
        if let c = out.coupledCluster {
            lines.append("")
            lines.append("Coupled cluster \(c.method)")
            lines.append("  converged:            \(c.converged ? "yes" : "no") (\(c.iterations) iterations)")
            if let t = c.triplesCorrection {
                lines.append("  CCSD energy:          " + Fmt.energy10(c.ccEnergy) + " E_h")
                lines.append("  (T) correction:       " + Fmt.energy10(t) + " E_h")
            }
            lines.append("  total energy:         " + Fmt.energy10(c.totalEnergy) + " E_h")
            lines.append("  correlation energy:   " + Fmt.energy10(c.totalCorrelation) + " E_h")
            lines.append(String(format: "  T1 diagnostic:        %.6f", c.t1Diagnostic))
        }
        if let j = out.job, j.type != "energy" {
            lines.append("")
            lines.append("Job: \(j.title)")
            if let label = j.finalEnergyLabel, let e = j.finalEnergy {
                lines.append("  final energy (\(label)): " + Fmt.energy10(e) + " E_h")
            }
            if let opt = out.optimization {
                let what = opt.isTransitionState ? "TS search:" : "optimisation:"
                let padding = String(repeating: " ", count: max(1, 22 - what.count))
                let state = opt.converged ? "converged" : "NOT converged"
                lines.append("  \(what)\(padding)\(state) in \(opt.steps) steps")
                if opt.isTransitionState {
                    lines.append("  Hessian:              " + Fmt.tsHessian(opt))
                }
                for p in opt.finalParameters {
                    lines.append("  " + Fmt.parameter(p))
                }
            }
            if let scan = out.scan {
                lines.append("  relaxed scan of \(scan.coordinate.definition): \(scan.points.count) points, " + scan.message)
                for p in scan.points {
                    lines.append(String(format: "    %3ld  %12.6f %@  %18.10f E_h  %10.4f kcal/mol%@", p.index + 1, p.value,
                                        scan.isAngular ? "deg" : "A  ", p.energy,
                                        p.relativeEnergy * Units.joulePerMolPerHartree / 1000 / Units.joulePerCalorie,
                                        p.converged ? "" : "  not converged"))
                }
                for x in scan.extrema {
                    lines.append(String(format: "    %@ near point %ld: %.6f %@, E = %.10f E_h (parabola)", x.kind,
                                        x.index + 1, x.value, scan.isAngular ? "deg" : "A", x.energy))
                }
            }
            if let f = out.frequencies {
                let list = f.modes.map { Fmt.wavenumber($0) }.joined(separator: ", ")
                lines.append("  harmonic wavenumbers: " + list + " cm-1 (" + (f.analytic ? "analytic Hessian" : "numerical Hessian") + ")")
                if let g = f.pointGroup {
                    lines.append("  point group:          " + g)
                }
                let imaginary = f.modes.filter { $0.imaginary }.count
                if JobChoice(rawValue: out.input.job)?.hasTransitionState == true {
                    lines.append("  transition state:     " + (imaginary == 1 ? "confirmed (one imaginary frequency)" : "NOT confirmed (\(imaginary) imaginary frequencies)"))
                }
            }
            if let t = out.thermochemistry {
                lines.append(String(format: "  thermochemistry at %.2f K, %.0f Pa:", t.temperature, t.pressure))
                lines.append("    ZPE " + String(format: "%.10f E_h", t.corrections.zeroPoint / Units.joulePerMolPerHartree))
                lines.append("    E + ZPE " + Fmt.energy10(t.totals.zeroPoint) + " E_h, H " + Fmt.energy10(t.totals.enthalpy)
                             + " E_h, G " + Fmt.energy10(t.totals.gibbs) + " E_h")
                lines.append(String(format: "    S %.4f J/(mol K)", t.terms.total.entropy))
            }
        }
        for warning in out.warnings ?? [] {
            lines.append("Warning: \(warning)")
        }
        if let message = out.error {
            lines.append("")
            lines.append("Error: \(message)")
        }
        return lines.joined(separator: "\n")
    }

    func copySummary() {
        guard let text = summaryText() else { return }
        let pasteboard = NSPasteboard.general
        pasteboard.clearContents()
        pasteboard.setString(text, forType: .string)
    }

    /// Plain-text export. H₂, LiH, H₂O: the output file of the job, exactly as in the Output
    /// tab and in psilab's output file. HeH⁺: summary, the FORTRAN-style report and, if a scan
    /// was run, the potential-curve table with its minima.
    func fullTextReport() -> String? {
        if let out = result, out.isGeneral, let file = out.output, !file.isEmpty {
            return file
        }
        guard let out = result, let summary = summaryText() else { return nil }
        let rule = String(repeating: "=", count: 78)
        var parts: [String] = [summary]
        if !out.report.isEmpty {
            let title = out.isGeneral
                ? "DETAILED REPORT"
                : "DETAILED REPORT (Hartree–Fock as in Szabo & Ostlund, Appendix B; then Full CI)"
            parts.append(rule + "\n" + title + "\n" + rule + "\n" + out.report)
        }
        if let curve = curveTableText() {
            parts.append(rule + "\n" + curve)
        }
        return parts.joined(separator: "\n\n") + "\n"
    }

    /// The potential-curve scan as a fixed-width table, or nil if no scan has been run.
    func curveTableText() -> String? {
        guard !curvePoints.isEmpty else { return nil }
        var lines: [String] = []
        lines.append("POTENTIAL CURVE")
        if let description = curveInputDescription {
            lines.append(description)
        }
        lines.append("")
        let hfHead = curveHFLabel == "HF" ? " E_HF (E_h)" : "E_RHF (E_h)"
        let hasCorrelated = curvePoints.contains { $0.correlated != nil }
        let corrHead = String(repeating: " ", count: max(0, 11 - curveCorrelatedLabel.count)) + "E_\(curveCorrelatedLabel) (E_h)"
        lines.append("    R (bohr)       R (Å)     \(hfHead)" + (hasCorrelated ? "  " + corrHead : ""))
        for p in curvePoints {
            var line = String(format: "%12.4f  %12.4f  %14.8f", p.R, p.R * Units.bohrToAngstrom, p.hf)
            if let e = p.correlated {
                line += String(format: "  %14.8f", e)
            }
            lines.append(line)
        }
        lines.append("")
        for m in curveMinima {
            lines.append(String(format: "%@ minimum (%@): R = %@, E = %.8f E_h",
                                m.method, m.refined ? "parabolic fit" : "grid point",
                                Fmt.lengthBoth(m.R), m.energy))
        }
        return lines.joined(separator: "\n")
    }

    /// File name stem for the exports, e.g. "PsiLab-HeH+-STO3G-R1.4632" or
    /// "PsiLab-H2O-CCSD(T)-aug-cc-pVDZ-opt-R1.8088".
    private func exportName(_ out: MacQCOutput) -> String {
        if out.isGeneral {
            let method = out.input.method
            let job = out.input.job == "energy" ? "" : "-" + out.input.job
            if out.input.isCustom {
                let name = (out.molecule?.formula ?? "custom").replacingOccurrences(of: " ", with: "")
                return String(format: "PsiLab-%@-%@-%@%@", name, method, out.input.basis, job)
            }
            return String(format: "PsiLab-%@-%@-%@%@-R%.4f", out.input.molecule, method, out.input.basis, job, out.input.R)
        }
        return String(format: "PsiLab-HeH+-STO%ldG-R%.4f", out.input.nGaussians, out.input.R)
    }

    func exportText() {
        // 13.0: with the archive of the results (the file can be opened again)
        guard let out = result, let text = exportableText() else { return }

        let panel = NSSavePanel()
        panel.title = "Export Results as Text"
        panel.allowedContentTypes = [UTType.plainText]
        panel.canCreateDirectories = true
        panel.nameFieldStringValue = exportName(out) + (out.isGeneral ? ".out" : ".txt")
        guard panel.runModal() == .OK, let url = panel.url else { return }
        do {
            try text.write(to: url, atomically: true, encoding: .utf8)
        } catch {
            errorMessage = "Could not write the file: \(error.localizedDescription)"
        }
    }

    func exportJSON() {
        guard let out = result else { return }

        let encoder = JSONEncoder()
        encoder.outputFormatting = [.prettyPrinted, .sortedKeys, .withoutEscapingSlashes]
        let data: Data
        do {
            // 17.0: the file in the backend's form ("Full CI"; the interface writes Full-CI)
            let encoded = try encoder.encode(out)
            data = Data(InterfaceText.backendForm(String(decoding: encoded, as: UTF8.self)).utf8)
        } catch {
            errorMessage = "Could not encode the results: \(error.localizedDescription)"
            return
        }

        let panel = NSSavePanel()
        panel.title = "Export Results as JSON"
        panel.allowedContentTypes = [UTType.json]
        panel.canCreateDirectories = true
        panel.nameFieldStringValue = exportName(out) + ".json"
        guard panel.runModal() == .OK, let url = panel.url else { return }
        do {
            try data.write(to: url, options: .atomic)
        } catch {
            errorMessage = "Could not write the file: \(error.localizedDescription)"
        }
    }

    var canExportOrCopy: Bool { result?.hf != nil || result?.rhf != nil || (result?.output?.isEmpty == false) }

    /// File > Export Potential Curve as Text… (the scan table with its minima).
    func exportCurve() {
        guard let text = curveTableText() else { return }
        let panel = NSSavePanel()
        panel.title = "Export Potential Curve as Text"
        panel.allowedContentTypes = [UTType.plainText]
        panel.canCreateDirectories = true
        panel.nameFieldStringValue = "PsiLab-curve-\(molecule.rawValue).txt"
        guard panel.runModal() == .OK, let url = panel.url else { return }
        do {
            try (text + "\n").write(to: url, atomically: true, encoding: .utf8)
        } catch {
            errorMessage = "Could not write the file: \(error.localizedDescription)"
        }
    }

    /// Selects a tab (and its page) from the View menu.
    func show(_ tab: DetailTab, methods: MethodsPage? = nil, vibration: VibrationPage? = nil, scan: ScanPage? = nil) {
        if let m = methods { methodsPage = m }
        if let v = vibration { vibrationPage = v }
        if let s = scan { scanPage = s }
        selectedTab = tab
    }
}
