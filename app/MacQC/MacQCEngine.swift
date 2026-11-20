//
//  MacQCEngine.swift
//  MacQC
//
//  Swift wrapper around the C API of the C++ backend (macqc_run / macqc_free)
//  and Codable models that mirror the JSON format documented in README.md.
//
//  All model types are marked `nonisolated` so they work unchanged whether or not
//  the project uses "Default Actor Isolation = MainActor" (the Xcode 26 template
//  default) and can be decoded on a background thread.
//

import Foundation

// MARK: - Input

/// Input JSON of `macqc_run`. HeH⁺ with `basis` "STO-1G", "STO-2G" or "STO-3G" selects the
/// Szabo–Ostlund calculation (keys R … printLevel); H₂, LiH, H₂O, H₂O₂, NH₃ and HeH⁺ with a
/// Gaussian basis set of the library the general path (R, basis, angle, method, job, …). The
/// backend echoes only the keys of the selected path (the Szabo–Ostlund echo has no "molecule"),
/// so every key is optional when decoding.
nonisolated struct MacQCInput: Codable, Sendable, Equatable {
    var molecule: String = "HeH+"
    var R: Double = 1.4632          // bohr (H2O, H2O2: O–H; NH3: N–H)
    var zeta1: Double = 2.0925
    var zeta2: Double = 1.24
    var ZA: Double = 2
    var ZB: Double = 1
    var nGaussians: Int = 3
    var exactErf: Bool = false
    var maxIter: Int = 25
    var convergence: Double = 1.0e-4
    var printLevel: Int = 2
    // General molecules
    var basis: String = "aug-cc-pVDZ"
    var angle: Double = 104.52      // degrees: H–O–H (H2O), O–O–H (H2O2), H–N–H (NH3)
    var method: String = "HF"       // 15.0: "HF" (14.0 "RHF"), "MP2" … "MP5", "CISD" … "QCISDTQ", "CCSD" … "CCSDTQ", "FCI"
    var basisDirectory: String = ""
    var threads: Int = 0
    var memoryMB: Double = 1024
    var fciMemoryMB: Double = 4096
    var estimateOnly: Bool = false
    // ΨLab 5.0 job types (general molecules)
    var job: String = "energy"      // "energy", "opt", "freq", "opt+freq", "opt+freq+sp"
    var spMethod: String = ""       // opt+freq+sp ("" = method)
    var spBasis: String = ""        // opt+freq+sp ("" = basis)
    var temperature: Double = 298.15
    var pressure: Double = 101325
    var optConvergence: String = "normal"
    var maxSteps: Int = 50
    // ΨLab 7.0: H2O2, the relaxed scan and the transition-state search
    var ROO: Double? = nil          // bohr, O–O of H2O2 (nil: the default)
    var dihedral: Double? = nil     // degrees, H–O–O–H of H2O2 (nil: the default)
    var scan: String = ""           // job "scan": ModRedundant line, e.g. "D 3 1 2 4 S 18 10.0"
    var scanStart: Double? = nil    // Å or degrees (nil: the value at the input geometry)
    var tsHessian: String = "calcfc" // "calcfc" (first step) or "calcall" (every step)
    var eigenTest: Bool = true      // the TS search needs exactly one negative eigenvalue
    // ΨLab 10.0: Hessians of Full-CI, CCSD and CCSD(T): "analytic" or "numerical"
    var hessian: String = "analytic"
    var initialHessian: String = "auto"   // minimisations: "auto", "analytic" or "model"
    // ΨLab 9.0: molecule "custom" (a structure file): atoms in Å, charge, multiplicity, name (LaTeX)
    var atoms: [CustomAtom]? = nil
    var charge: Int? = nil
    var multiplicity: Int? = nil
    var title: String? = nil
    // ΨLab 15.0: the reference (nil: automatic), HeH⁺ in STO-nG as a general molecule
    // (szaboOstlund false), computed Hessians of minimisations, RecalcFC=N, the step (rfo, gediis)
    var reference: String? = nil
    var szaboOstlund: Bool? = nil
    var optHessian: String = "auto"
    var recalcFC: Int = 0
    var optAlgorithm: String = "rfo"

    init() {}

    private enum CodingKeys: String, CodingKey {
        case molecule, R, zeta1, zeta2, ZA, ZB, nGaussians, exactErf, maxIter, convergence, printLevel
        case basis, angle, method, basisDirectory, threads, memoryMB, fciMemoryMB, estimateOnly
        case job, spMethod, spBasis, temperature, pressure, optConvergence, maxSteps
        case ROO, dihedral, scan, scanStart, tsHessian, eigenTest, hessian, initialHessian
        case atoms, charge, multiplicity, title
        case reference, szaboOstlund, optHessian, recalcFC, optAlgorithm
    }

    init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        let d = MacQCInput()
        molecule = try c.decodeIfPresent(String.self, forKey: .molecule) ?? d.molecule
        R = try c.decodeIfPresent(Double.self, forKey: .R) ?? d.R
        zeta1 = try c.decodeIfPresent(Double.self, forKey: .zeta1) ?? d.zeta1
        zeta2 = try c.decodeIfPresent(Double.self, forKey: .zeta2) ?? d.zeta2
        ZA = try c.decodeIfPresent(Double.self, forKey: .ZA) ?? d.ZA
        ZB = try c.decodeIfPresent(Double.self, forKey: .ZB) ?? d.ZB
        nGaussians = try c.decodeIfPresent(Int.self, forKey: .nGaussians) ?? d.nGaussians
        exactErf = try c.decodeIfPresent(Bool.self, forKey: .exactErf) ?? d.exactErf
        maxIter = try c.decodeIfPresent(Int.self, forKey: .maxIter) ?? d.maxIter
        convergence = try c.decodeIfPresent(Double.self, forKey: .convergence) ?? d.convergence
        printLevel = try c.decodeIfPresent(Int.self, forKey: .printLevel) ?? d.printLevel
        basis = try c.decodeIfPresent(String.self, forKey: .basis) ?? d.basis
        angle = try c.decodeIfPresent(Double.self, forKey: .angle) ?? d.angle
        method = try c.decodeIfPresent(String.self, forKey: .method) ?? d.method
        basisDirectory = try c.decodeIfPresent(String.self, forKey: .basisDirectory) ?? d.basisDirectory
        threads = try c.decodeIfPresent(Int.self, forKey: .threads) ?? d.threads
        memoryMB = try c.decodeIfPresent(Double.self, forKey: .memoryMB) ?? d.memoryMB
        fciMemoryMB = try c.decodeIfPresent(Double.self, forKey: .fciMemoryMB) ?? d.fciMemoryMB
        estimateOnly = try c.decodeIfPresent(Bool.self, forKey: .estimateOnly) ?? d.estimateOnly
        job = try c.decodeIfPresent(String.self, forKey: .job) ?? d.job
        spMethod = try c.decodeIfPresent(String.self, forKey: .spMethod) ?? d.spMethod
        spBasis = try c.decodeIfPresent(String.self, forKey: .spBasis) ?? d.spBasis
        temperature = try c.decodeIfPresent(Double.self, forKey: .temperature) ?? d.temperature
        pressure = try c.decodeIfPresent(Double.self, forKey: .pressure) ?? d.pressure
        optConvergence = try c.decodeIfPresent(String.self, forKey: .optConvergence) ?? d.optConvergence
        maxSteps = try c.decodeIfPresent(Int.self, forKey: .maxSteps) ?? d.maxSteps
        ROO = try c.decodeIfPresent(Double.self, forKey: .ROO)
        dihedral = try c.decodeIfPresent(Double.self, forKey: .dihedral)
        scan = try c.decodeIfPresent(String.self, forKey: .scan) ?? d.scan
        scanStart = try c.decodeIfPresent(Double.self, forKey: .scanStart)
        tsHessian = try c.decodeIfPresent(String.self, forKey: .tsHessian) ?? d.tsHessian
        eigenTest = try c.decodeIfPresent(Bool.self, forKey: .eigenTest) ?? d.eigenTest
        hessian = try c.decodeIfPresent(String.self, forKey: .hessian) ?? d.hessian
        initialHessian = try c.decodeIfPresent(String.self, forKey: .initialHessian) ?? d.initialHessian
        atoms = try c.decodeIfPresent([CustomAtom].self, forKey: .atoms)
        charge = try c.decodeIfPresent(Int.self, forKey: .charge)
        multiplicity = try c.decodeIfPresent(Int.self, forKey: .multiplicity)
        title = try c.decodeIfPresent(String.self, forKey: .title)
        reference = try c.decodeIfPresent(String.self, forKey: .reference)
        szaboOstlund = try c.decodeIfPresent(Bool.self, forKey: .szaboOstlund)
        optHessian = try c.decodeIfPresent(String.self, forKey: .optHessian) ?? d.optHessian
        recalcFC = try c.decodeIfPresent(Int.self, forKey: .recalcFC) ?? d.recalcFC
        optAlgorithm = try c.decodeIfPresent(String.self, forKey: .optAlgorithm) ?? d.optAlgorithm
        // The echo of the Szabo–Ostlund calculation has no "molecule" and no "basis".
        if !c.contains(.molecule) {
            basis = "STO-\(nGaussians)G"
        }
    }

    /// True for the general Gaussian-basis path: H₂, LiH, H₂O, H₂O₂, NH₃, and HeH⁺ in a
    /// Gaussian basis set of the library (false: HeH⁺ in STO-nG, Szabo–Ostlund).
    var isGeneral: Bool { molecule != "HeH+" || !basis.uppercased().hasPrefix("STO") || szaboOstlund == false }

    /// ΨLab 9.0: a custom molecule (atoms from a structure file).
    var isCustom: Bool { molecule == "custom" }
}

// MARK: - Output

/// The 16 unique integrals, named exactly as in the FORTRAN program.
nonisolated struct MacQCIntegrals: Codable, Sendable, Equatable {
    var S12: Double
    var T11: Double
    var T12: Double
    var T22: Double
    var V11A: Double
    var V12A: Double
    var V22A: Double
    var V11B: Double
    var V12B: Double
    var V22B: Double
    var V1111: Double
    var V2111: Double
    var V2121: Double
    var V2211: Double
    var V2221: Double
    var V2222: Double

    /// Name/value pairs in the order used by the printed report.
    var namedValues: [(name: String, value: Double)] {
        [
            ("S12", S12), ("T11", T11), ("T12", T12), ("T22", T22),
            ("V11A", V11A), ("V12A", V12A), ("V22A", V22A),
            ("V11B", V11B), ("V12B", V12B), ("V22B", V22B),
            ("V1111", V1111), ("V2111", V2111), ("V2121", V2121),
            ("V2211", V2211), ("V2221", V2221), ("V2222", V2222),
        ]
    }
}

nonisolated struct MacQCSCFIteration: Codable, Sendable, Equatable, Identifiable {
    var iteration: Int
    var electronicEnergy: Double
    var deltaP: Double

    var id: Int { iteration }
}

nonisolated struct MacQCHF: Codable, Sendable, Equatable {
    var converged: Bool
    var iterations: Int
    var electronicEnergy: Double
    var nuclearRepulsion: Double
    var totalEnergy: Double
    var orbitalEnergies: [Double]
    var C: [[Double]]
    var P: [[Double]]
    var S: [[Double]]
    var X: [[Double]]
    var H: [[Double]]
    var F: [[Double]]
    var G: [[Double]]
    var PS: [[Double]]
    var integrals: MacQCIntegrals
    var scfHistory: [MacQCSCFIteration]
}

nonisolated struct MacQCFCIState: Codable, Sendable, Equatable {
    var energy: Double
    var S2: Double
    var spin: String
    var vector: [Double]
}

nonisolated struct MacQCSingletWeight: Codable, Sendable, Equatable {
    var configuration: String
    var weight: Double
}

nonisolated struct MacQCFCI: Codable, Sendable, Equatable {
    var totalEnergy: Double
    var electronicEnergy: Double
    var correlationEnergy: Double
    var states: [MacQCFCIState]
    var groundStateVector: [Double]
    var determinants: [String]
    var singletWeights: [MacQCSingletWeight]
}

// MARK: - Output of the general molecules (H2, LiH, H2O)

nonisolated struct MacQCAtom: Codable, Sendable, Equatable {
    var symbol: String
    var Z: Int
    var position: [Double]   // bohr
    var ecp: MacQCEcp?       // 14.0: in the orbitals object, the atom's effective core potential
}

/// 14.0: an effective core potential U = U_L + Σ_l |lm⟩[U_l − U_L]⟨lm| with Gaussian terms
/// c·exp(−a r²), each term [a, c].
nonisolated struct MacQCEcp: Codable, Sendable, Equatable {
    var core: Int                    // core electrons replaced
    var L: Int                       // angular momentum of the local part
    var local: [[Double]]            // U_L
    var projected: [[[Double]]]      // U_l − U_L, l = 0 … L−1
}

nonisolated struct MacQCMoleculeInfo: Codable, Sendable, Equatable {
    var name: String
    var electrons: Int
    var charge: Int
    var R: Double            // bohr
    var angle: Double?       // degrees (H2O, H2O2, NH3)
    var ROO: Double?         // bohr (H2O2)
    var dihedral: Double?    // degrees (H2O2)
    var nuclearRepulsion: Double
    var atoms: [MacQCAtom]
    // ΨLab 9.0
    var formula: String?     // Hill order, e.g. "CO2"
    var custom: Bool?
    var title: String?       // custom molecules: the name from the structure file (LaTeX)
    var coreElectrons: Int?  // 14.0: replaced by effective core potentials (0 or nil: none)
}

nonisolated struct MacQCAtomBasis: Codable, Sendable, Equatable {
    var symbol: String
    var basis: String
    var contraction: String
    var primitives: String
    var functions: Int
    var ecp: String?          // 14.0: "ECP28" (28 core electrons replaced), nil: all electrons
    var ecpElectrons: Int?
}

nonisolated struct MacQCBasisInfo: Codable, Sendable, Equatable {
    var name: String
    var functions: Int
    var shells: Int
    var maxL: Int
    var atoms: [MacQCAtomBasis]
    /// 14.0 (nil in older results: spherical): the d shells are Cartesian (6D, the 6-31G family),
    /// and the kind of functions in words, e.g. "spherical harmonics (5D, 7F)".
    var cartesianD: Bool?
    var functionTypes: String?
}

nonisolated struct MacQCFCIEstimate: Codable, Sendable, Equatable {
    var orbitals: Int
    var alphaElectrons: Int
    var betaElectrons: Int
    var alphaStrings: Double
    var betaStrings: Double
    var determinants: Double
    var memoryMB: Double
    var seconds: Double
    var feasible: Bool
    var reason: String?
}

nonisolated struct MacQCEstimate: Codable, Sendable, Equatable {
    var basisFunctions: Int
    var uniqueIntegrals: Double
    var integralMemoryMB: Double
    var integralsStored: Bool
    var threads: Int
    var rhfSeconds: Double
    var totalSeconds: Double
    var longCalculation: Bool
    var fci: MacQCFCIEstimate?
    var coupledCluster: MacQCCCEstimate?
}

nonisolated struct MacQCCCEstimate: Codable, Sendable, Equatable {
    var method: String
    var engine: String
    var occupied: Int
    var virtuals: Int
    var amplitudes: Double
    var memoryMB: Double
    var ladderIntegralsStored: Bool
    var iterations: Int
    var iterationSeconds: Double
    var triplesSeconds: Double
    var seconds: Double
    var feasible: Bool
    var reason: String?
}

nonisolated struct MacQCRHFIteration: Codable, Sendable, Equatable, Identifiable {
    var iteration: Int
    var totalEnergy: Double
    var deltaE: Double
    var gradient: Double
    var diis: Bool

    var id: Int { iteration }
}

nonisolated struct MacQCRHF: Codable, Sendable, Equatable {
    var converged: Bool
    var iterations: Int
    var totalEnergy: Double
    var electronicEnergy: Double
    var nuclearRepulsion: Double
    var kineticEnergy: Double
    var nuclearAttractionEnergy: Double
    var electronRepulsionEnergy: Double
    var virialRatio: Double
    var orbitalEnergies: [Double]
    var occupiedOrbitals: Int
    var homo: Double?
    var lumo: Double?
    var smallestOverlapEigenvalue: Double
    var removedFunctions: Int
    var integralsStored: Bool
    var uniqueIntegrals: Double
    var threads: Int
    var integralSeconds: Double
    var scfSeconds: Double
    var scfHistory: [MacQCRHFIteration]
    // ΨLab 15.0: the reference (RHF, UHF, ROHF), spin and ⟨S²⟩ (open shells)
    var reference: String?
    var multiplicity: Int?
    var alphaElectrons: Int?
    var betaElectrons: Int?
    var S2: Double?
    var S2Exact: Double?
}

/// ΨLab 15.0: a method of the spin-orbital engine (MP2–MP5, CISD … CISDTQ, QCISD … QCISDTQ, and
/// coupled cluster on open-shell references): energies, the energy terms (E(2), E(3), (T), …).
nonisolated struct MacQCCorrelatedTerm: Codable, Sendable, Equatable {
    var label: String
    var value: Double
}

nonisolated struct MacQCCorrelated: Codable, Sendable, Equatable {
    var method: String
    var label: String?            // with the reference prefix, e.g. "UMP2"
    var reference: String?
    var converged: Bool
    var iterations: Int
    var referenceEnergy: Double
    var correlationEnergy: Double
    var totalEnergy: Double
    var t1Diagnostic: Double?
    var terms: [MacQCCorrelatedTerm]?
    var seconds: Double?

    var title: String { label ?? method }
}

nonisolated struct MacQCDeterminant: Codable, Sendable, Equatable {
    var coefficient: Double
    var alpha: [Int]
    var beta: [Int]
    var excitation: Int
}

nonisolated struct MacQCDavidsonIteration: Codable, Sendable, Equatable, Identifiable {
    var iteration: Int
    var totalEnergy: Double
    var residual: Double
    var subspace: Int

    var id: Int { iteration }
}

nonisolated struct MacQCFullCI: Codable, Sendable, Equatable {
    var converged: Bool
    var iterations: Int
    var totalEnergy: Double
    var electronicEnergy: Double
    var correlationEnergy: Double
    var S2: Double
    var hfWeight: Double
    var orbitals: Int
    var alphaElectrons: Int
    var betaElectrons: Int
    var alphaStrings: Double
    var betaStrings: Double
    var determinants: Double
    var leadingDeterminants: [MacQCDeterminant]
    var naturalOccupations: [Double]
    var history: [MacQCDavidsonIteration]
    var transformSeconds: Double
    var ciSeconds: Double
    // ΨLab 15.0: spin adaptation (configuration state functions of the requested spin)
    var csfs: Double?
    var spin: Double?
}

nonisolated struct MacQCCCAmplitude: Codable, Sendable, Equatable {
    var value: Double
    var rank: Int
    var occupied: [Int]
    var virtuals: [Int]
    var spins: String
}

nonisolated struct MacQCCCIteration: Codable, Sendable, Equatable, Identifiable {
    var iteration: Int
    var totalEnergy: Double
    var correlationEnergy: Double
    var deltaE: Double
    var residual: Double

    var id: Int { iteration }
}

/// Coupled-cluster result (CCSD, CCSD(T), CCSDT, CCSDTQ).
nonisolated struct MacQCCoupledCluster: Codable, Sendable, Equatable {
    var method: String
    var engine: String
    var converged: Bool
    var iterations: Int
    var occupied: Int
    var virtuals: Int
    var referenceEnergy: Double
    var correlationEnergy: Double      // without (T)
    var ccEnergy: Double               // reference + correlation
    var triplesCorrection: Double?     // CCSD(T) only
    var totalEnergy: Double
    var t1Diagnostic: Double
    var amplitudeCounts: [Double]      // index = rank
    var amplitudeNorms: [Double]
    var leadingAmplitudes: [MacQCCCAmplitude]
    var history: [MacQCCCIteration]
    var ladder: String?
    var transformSeconds: Double
    var ccSeconds: Double
    var triplesSeconds: Double?

    /// Correlation energy including (T).
    var totalCorrelation: Double { totalEnergy - referenceEnergy }
}

// MARK: - Jobs of ΨLab 5.0 (optimisation, frequencies, thermochemistry)

/// A geometry: atoms (bohr) and the parameters of the inputs (R in bohr; H₂O: both O–H
/// distances, their mean and the H–O–H angle in degrees; H₂O₂: also R(O–O), both O–O–H angles
/// and the dihedral angle; NH₃: the three N–H distances and H–N–H angles, their means and the
/// pyramidality 360° − ΣH–N–H).
nonisolated struct MacQCGeometry: Codable, Sendable, Equatable {
    var atoms: [MacQCAtom]
    var R: Double?
    var R1: Double?
    var R2: Double?
    var R3: Double?
    var ROO: Double?
    var angle: Double?
    var angle1: Double?
    var angle2: Double?
    var angle3: Double?
    var dihedral: Double?
    var pyramidal: Double?
}

nonisolated struct MacQCJobStep: Codable, Sendable, Equatable, Identifiable {
    var name: String
    var method: String
    var basis: String
    var how: String
    var analytic: Bool
    var evaluations: Int
    var evaluationKind: String      // "steps", "gradients", "energies", "Hessians"
    var seconds: Double
    var memoryMB: Double

    var id: String { name }
}

nonisolated struct MacQCJobEstimate: Codable, Sendable, Equatable {
    var steps: [MacQCJobStep]
    var totalSeconds: Double
    var memoryMB: Double
    var longCalculation: Bool
    var refusal: String?
    var warnings: [String]
    var notes: [String]
}

nonisolated struct MacQCJobInfo: Codable, Sendable, Equatable {
    var type: String
    var title: String
    var termination: String          // "normal", "error", "refused", "cancelled", "estimate"
    var completed: Bool
    var message: String?
    var estimate: MacQCJobEstimate
    var warnings: [String]
    var seconds: Double
    var initialGeometry: MacQCGeometry?
    var finalGeometry: MacQCGeometry?
    var finalEnergy: Double?
    var finalEnergyLabel: String?
    var singlePointMethod: String?
    var singlePointBasis: String?
}

nonisolated struct MacQCParameter: Codable, Sendable, Equatable, Identifiable {
    var name: String
    var value: Double                // bohr or degrees
    var angle: Bool

    var id: String { name }
}

nonisolated struct MacQCInternalCoordinate: Codable, Sendable, Equatable, Identifiable {
    var name: String
    var definition: String
    var value: Double                // bohr or radians
    var gradient: Double             // E_h/bohr or E_h/rad
    var step: Double?

    var id: String { name }
}

nonisolated struct MacQCConvergenceItem: Codable, Sendable, Equatable, Identifiable {
    var name: String
    var value: Double?
    var threshold: Double
    var met: Bool
    var applicable: Bool

    var id: String { name }
}

/// A step of the transition-state search (ΨLab 7.0): the Hessian's eigenvalues (internal
/// coordinates), the followed mode and the shifts of the partitioned RFO step.
nonisolated struct MacQCTSStep: Codable, Sendable, Equatable {
    var hessianEigenvalues: [Double]
    var negativeEigenvalues: Int
    var followedEigenvalue: Double?
    var followedGradient: Double?
    var followedOverlap: Double?
    var followedMode: [Double]?
    var lambdaP: Double?
    var lambdaN: Double?
    var hessianComputed: Bool
}

nonisolated struct MacQCOptimizationStep: Codable, Sendable, Equatable, Identifiable {
    var step: Int
    var energy: Double
    var deltaE: Double?
    var geometry: MacQCGeometry
    var coordinates: [MacQCInternalCoordinate]
    var convergence: [MacQCConvergenceItem]
    var converged: Bool
    var hessianUpdate: String
    var trustRadius: Double
    var stepNorm: Double
    var stepScaled: Bool
    var predictedChange: Double
    var ratio: Double?
    var warnings: [String]
    var seconds: Double
    var ts: MacQCTSStep?             // transition-state search only

    var id: Int { step }

    /// Value of a convergence item by position: 0 max force, 1 RMS force, 2 max step, 3 RMS step.
    func criterion(_ index: Int) -> Double? {
        index < convergence.count ? convergence[index].value : nil
    }
}

nonisolated struct MacQCCriteria: Codable, Sendable, Equatable {
    var maxForce: Double
    var rmsForce: Double
    var maxDisplacement: Double
    var rmsDisplacement: Double
    var energy: Double
}

nonisolated struct MacQCOptimization: Codable, Sendable, Equatable {
    var converged: Bool
    var steps: Int
    var message: String
    var initialHessian: String
    var gradient: String?
    var criteria: MacQCCriteria
    var initialEnergy: Double
    var finalEnergy: Double
    var initialGeometry: MacQCGeometry
    var finalGeometry: MacQCGeometry
    var initialParameters: [MacQCParameter]
    var finalParameters: [MacQCParameter]
    var symmetricStart: Bool
    var symmetricEnd: Bool
    var symmetryDeviation: Double
    var history: [MacQCOptimizationStep]
    var seconds: Double
    // ΨLab 7.0: transition-state search
    var transitionState: Bool?
    var tsHessian: String?           // "calcfc" or "calcall"
    var hessianEvaluations: Int?
    var negativeEigenvalues: Int?    // of the last Hessian
    var eigenTestFailed: Bool?

    var isTransitionState: Bool { transitionState ?? false }
}

nonisolated struct MacQCNormalMode: Codable, Sendable, Equatable {
    var wavenumber: Double           // cm⁻¹, negative = imaginary
    var imaginary: Bool
    var symmetry: String
    var reducedMass: Double          // u
    var forceConstant: Double        // mdyn/Å
    var irIntensity: Double?         // km/mol
    var displacement: [[Double]]     // normalised Cartesian displacement per atom
}

nonisolated struct MacQCFrequencies: Codable, Sendable, Equatable {
    var analytic: Bool
    var description: String
    var numericalReason: String?
    var gradientEvaluations: Int
    var energyEvaluations: Int
    var coordinates: [String]
    var intensities: Bool
    var intensityDescription: String
    var dipole: [Double]?
    var linear: Bool
    var pointGroup: String?          // "C2v", "C2h", "D3h", … (four atoms; nil when not determined)
    var vibrations: Int
    var masses: [Double]
    var isotopes: [String]
    var geometry: MacQCGeometry
    var energy: Double
    var lowFrequencies: [Double]
    var modes: [MacQCNormalMode]
    var seconds: Double
}

nonisolated struct MacQCThermoTerm: Codable, Sendable, Equatable {
    var energy: Double               // J/mol
    var entropy: Double              // J/(mol K)
    var cv: Double                   // J/(mol K)
    var q: Double
}

nonisolated struct MacQCThermoTerms: Codable, Sendable, Equatable {
    var translation: MacQCThermoTerm
    var rotation: MacQCThermoTerm
    var vibration: MacQCThermoTerm
    var electronic: MacQCThermoTerm
    var total: MacQCThermoTerm
}

/// Thermal corrections (J/mol) or totals (E_h).
nonisolated struct MacQCThermoValues: Codable, Sendable, Equatable {
    var zeroPoint: Double
    var thermalEnergy: Double
    var enthalpy: Double
    var gibbs: Double
    var helmholtz: Double
}

nonisolated struct MacQCThermochemistry: Codable, Sendable, Equatable {
    var temperature: Double
    var pressure: Double
    var symmetryNumber: Int
    var linear: Bool
    var mass: Double
    var moments: [Double]                    // u bohr²
    var rotationalConstantsGHz: [Double]
    var rotationalConstantsCm: [Double]
    var rotationalTemperatures: [Double]
    var vibrationalTemperatures: [Double]
    var excludedModes: Int
    var warnings: [String]
    var terms: MacQCThermoTerms
    var zeroPointEnergy: Double              // J/mol
    var corrections: MacQCThermoValues       // J/mol
    var electronicEnergy: Double             // E_h
    var totals: MacQCThermoValues            // E_h
}

nonisolated struct MacQCSinglePoint: Codable, Sendable, Equatable {
    var method: String
    var energy: Double
    var converged: Bool
    var rhfEnergy: Double
    var correlationEnergy: Double
    var geometry: MacQCGeometry
    var warnings: [String]
    var seconds: Double
}

// MARK: - Relaxed scan (ΨLab 7.0)

nonisolated struct MacQCScanCoordinate: Codable, Sendable, Equatable {
    var type: String                 // "bond", "angle", "dihedral"
    var letter: String               // "B", "A", "D"
    var atoms: [Int]                 // atom numbers from 1
    var definition: String           // "D(H3–O1–O2–H4)"
    var unit: String                 // "angstrom" or "degrees"
}

nonisolated struct MacQCScanOptimizationStep: Codable, Sendable, Equatable {
    var step: Int
    var energy: Double
    var maxForce: Double?
}

nonisolated struct MacQCScanPoint: Codable, Sendable, Equatable, Identifiable {
    var index: Int
    var value: Double                // target, Å or degrees
    var measured: Double             // at the optimised geometry
    var energy: Double
    var relativeEnergy: Double       // E_h above the lowest point
    var converged: Bool
    var steps: Int
    var maxForce: Double?
    var geometry: MacQCGeometry
    var parameters: [MacQCParameter]
    var optimizationSteps: [MacQCScanOptimizationStep]
    var warnings: [String]
    var seconds: Double

    var id: Int { index }
}

nonisolated struct MacQCScanExtremum: Codable, Sendable, Equatable, Identifiable {
    var kind: String                 // "minimum" or "maximum"
    var index: Int                   // scan point next to it
    var value: Double                // parabola through the point and its neighbours (Å or degrees)
    var energy: Double

    var id: String { "\(kind)\(index)" }
}

/// A relaxed potential-energy scan (Gaussian Opt=ModRedundant with S): one constrained
/// optimisation per value of the scanned coordinate.
nonisolated struct MacQCScan: Codable, Sendable, Equatable {
    var coordinate: MacQCScanCoordinate
    var start: Double
    var stepSize: Double
    var steps: Int
    var inputValue: Double
    var relaxed: Bool                // false for a diatomic (single points)
    var converged: Bool
    var message: String
    var lowest: Int
    var highest: Int
    var points: [MacQCScanPoint]
    var extrema: [MacQCScanExtremum]
    var seconds: Double

    var isAngular: Bool { coordinate.unit == "degrees" }
    var unitSymbol: String { isAngular ? "°" : "Å" }
}

// MARK: - Molecular orbitals (ΨLab 8.0)

/// A contracted shell of real solid-harmonic Gaussians or (14.0) of Cartesian Gaussians
/// (coefficients include all normalisation).
nonisolated struct MacQCShell: Codable, Sendable, Equatable {
    var l: Int
    var atom: Int
    var center: [Double]             // bohr
    var exponents: [Double]
    var coefficients: [Double]
    var first: Int                   // first basis function (m = -l … l, or xx, xy, xz, yy, yz, zz)
    /// 14.0: false for a Cartesian shell (the d shells of the 6-31G family, 6D); nil (older
    /// results) or true: spherical. Passed back to the backend with the orbitals.
    var pure: Bool?

    var isCartesian: Bool { pure == false }
    var size: Int { isCartesian ? (l + 1) * (l + 2) / 2 : 2 * l + 1 }
}

nonisolated struct MacQCOrbitalInfo: Codable, Sendable, Equatable, Identifiable {
    var number: Int                  // 1-based, in energy order (canonical) or occupation order (natural)
    var energy: Double?              // E_h (canonical orbitals)
    var occupation: Double           // 2 or 0 (RHF), natural occupation (Full-CI)
    var label: String                // "HOMO", "HOMO-1", "LUMO", "LUMO+2" or ""

    var id: Int { number }
}

nonisolated struct MacQCOrbitalSet: Codable, Sendable, Equatable {
    var kind: String                 // "canonical" or "natural"
    var title: String
    var orbitals: [MacQCOrbitalInfo]
    var coefficients: [[Double]]     // per orbital, the AO coefficients

    var isNatural: Bool { kind == "natural" }
}

/// The "orbitals" object of the output: atoms, basis and orbital sets (RHF canonical orbitals,
/// all occupied and up to ten virtual ones; Full-CI natural orbitals; for HeH⁺ in STO-nG the
/// two Szabo–Ostlund orbitals). Passed back to the backend to draw an orbital.
nonisolated struct MacQCOrbitals: Codable, Sendable, Equatable {
    var basis: String
    var functions: Int
    var atoms: [MacQCAtom]
    var shells: [MacQCShell]
    var sets: [MacQCOrbitalSet]
}

/// An isosurface: positions (bohr) and unit normals, three floats per vertex; three indices per
/// triangle, counter-clockwise seen from outside.
nonisolated struct OrbitalMesh: Sendable, Equatable {
    var positions: [Float] = []
    var normals: [Float] = []
    var indices: [Int32] = []

    var vertexCount: Int { positions.count / 3 }
    var triangleCount: Int { indices.count / 3 }
}

/// Output JSON of `macqc_run`. HeH+: `hf` / `fci` (null when the backend reports an error
/// before producing them). H2, LiH, H2O: `molecule`, `basis`, `estimate`, `notes`,
/// `warnings`, `refusal`, `rhf`, `fullCI`, `coupledCluster`.
nonisolated struct MacQCOutput: Codable, Sendable, Equatable {
    var input: MacQCInput
    var error: String?
    var hf: MacQCHF?
    var fci: MacQCFCI?
    var report: String
    var molecule: MacQCMoleculeInfo?
    var basis: MacQCBasisInfo?
    var estimate: MacQCEstimate?
    var notes: [String]?
    var warnings: [String]?
    var refusal: String?
    var rhf: MacQCRHF?
    var fullCI: MacQCFullCI?
    var coupledCluster: MacQCCoupledCluster?
    var correlated: MacQCCorrelated?   // 15.0: the spin-orbital engine's methods
    var reference: String?             // 15.0: RHF, UHF or ROHF
    // ΨLab 5.0 (H2, LiH, H2O): the job and its output file; rhf / fullCI / coupledCluster /
    // molecule / basis describe the final energy calculation of the job.
    var job: MacQCJobInfo?
    var optimization: MacQCOptimization?
    var frequencies: MacQCFrequencies?
    var thermochemistry: MacQCThermochemistry?
    var singlePoint: MacQCSinglePoint?
    var scan: MacQCScan?             // ΨLab 7.0: job "scan"
    var orbitals: MacQCOrbitals?     // ΨLab 8.0: molecular orbitals of the final calculation
    var output: String?

    /// Energy of the correlated method of this result (Full-CI or coupled cluster).
    var correlatedEnergy: Double? { coupledCluster?.totalEnergy ?? correlated?.totalEnergy ?? fullCI?.totalEnergy }

    /// 15.0: the reference of the SCF ("RHF", "UHF", "ROHF").
    var referenceLabel: String { rhf?.reference ?? reference ?? "RHF" }

    /// Title of the method of this result: "HF", "CCSD(T)", "Full-CI", …
    var methodTitle: String { MethodChoice.title(of: input.method) }

    /// True for the general path (all molecules except HeH⁺ in STO-nG).
    var isGeneral: Bool { input.isGeneral }

    /// ΨLab 9.0: the name of the molecule as LaTeX source: a custom molecule's name (or its
    /// formula), "H$_2$O" for the built-in molecules.
    var moleculeNameSource: String {
        if input.isCustom {
            let t = (input.title ?? molecule?.title ?? "").trimmingCharacters(in: .whitespaces)
            if !t.isEmpty { return t }
            let f = molecule?.formula ?? "?"
            return f.replacingOccurrences(of: #"(\d+)"#, with: #"\$_{$1}\$"#, options: .regularExpression)
        }
        // 15.0: the atoms and molecules of the menu (Unicode names, e.g. "H₅⁺", "NO₂")
        if let p = MoleculePresets.find(input.molecule) { return p.title }
        let name = molecule?.name ?? input.molecule
        if name == "HeH+" { return "HeH$^+$" }
        return name.replacingOccurrences(of: #"(\d+)"#, with: #"\$_{$1}\$"#, options: .regularExpression)
    }

    /// The name in plain text with Unicode subscripts ("CO₂", "H₂O", "HeH⁺").
    var moleculeTitle: String { LaTeX.plain(moleculeNameSource) }

    /// 17.0: the level of the result as Level1//Level2 when a single point followed the
    /// optimisation ("Full-CI/aug-cc-pCV5Z//Full-CI/aug-cc-pCVQZ": the single-point energy at
    /// Level1 on the geometry optimised at Level2), else method/basis set.
    var compositeLevel: String {
        let level = MethodChoice.title(of: input.method) + "/" + input.basis
        guard let sp = singlePoint else { return level }
        let spBasis = job?.singlePointBasis ?? (input.spBasis.isEmpty ? input.basis : input.spBasis)
        return MethodChoice.title(of: sp.method) + "/" + spBasis + "//" + level
    }
}

// MARK: - Jobs (cancellation and progress)

/// Progress of a running calculation, as reported by the backend.
nonisolated struct MacQCProgress: Sendable, Equatable {
    var stage: Int = 0        // 0 idle, 1 integrals, 2 SCF, 3 transformation, 4 Full-CI, 5 done,
                              // 6 coupled cluster, 7 (T) correction (iteration = percent);
                              // 5.0: 8 gradient, 9 Hessian, 10 CPHF, 11 Lambda, 12 optimisation,
                              // 13 numerical Hessian; 7.0: 14 scan point (iteration = point,
                              // residual = number of points)
    var iteration: Int = 0
    var energy: Double = 0
    var residual: Double = 0

    var text: String {
        switch stage {
        case 1:
            return "Integrals…"
        case 2:
            return String(format: "SCF iteration %ld · E = %.8f E_h · gradient %.1e", iteration, energy, residual)
        case 3:
            return "Transforming the integrals to the RHF orbitals…"
        case 4:
            if iteration == 0 { return "Full-CI: setting up the string Hamiltonian…" }
            return String(format: "Full-CI · Davidson iteration %ld · E = %.8f E_h · residual %.1e",
                          iteration, energy, residual)
        case 5:
            return "Done"
        case 6:
            if iteration == 0 { return "Coupled cluster: setting up…" }
            return String(format: "Coupled cluster iteration %ld · E = %.8f E_h · max |Δt| %.1e",
                          iteration, energy, residual)
        case 7:
            return String(format: "(T) correction · %ld %% done", iteration)
        case 8:
            return "Gradient…"
        case 9:
            // ΨLab 10.0: the parts of an analytic Hessian (residual: perturbation number)
            let k = Int(residual)
            switch iteration {
            case 1: return k > 0 ? "Hessian · derivative integrals (coordinate \(k))" : "Hessian · derivative integrals…"
            case 2: return "Hessian · orbital Lagrangian and Z-vector…"
            case 3: return "Hessian · second-derivative integrals…"
            case 4: return "Hessian · orbital response terms (perturbation \(k))"
            case 5: return "Hessian · Z-vector integral terms (perturbation \(k))"
            case 6: return "Hessian · response contractions (perturbation \(k))"
            case 7: return "Hessian · amplitude responses (perturbation \(k))"
            default: return "Hessian…"
            }
        case 10:
            return String(format: "Orbital response · iteration %ld · residual %.1e", iteration, residual)
        case 11:
            return String(format: "Lambda equations · iteration %ld · max change %.1e", iteration, residual)
        case 12:
            return String(format: "Optimisation step %ld", iteration)
        case 13:
            return String(format: "Numerical Hessian · displacement %ld", iteration)
        case 14:
            return String(format: "Scan point %ld of %ld", iteration, Int(residual))
        default:
            return "Starting…"
        }
    }
}

/// A backend job handle: lets another thread cancel a running calculation and read its progress.
nonisolated final class MacQCJob: @unchecked Sendable {
    private let handle: OpaquePointer?

    init() {
        handle = macqc_job_create()
    }

    deinit {
        macqc_job_free(handle)
    }

    func cancel() {
        macqc_job_cancel(handle)
    }

    func progress() -> MacQCProgress {
        var stage: Int32 = 0
        var iteration: Int32 = 0
        var energy = 0.0
        var residual = 0.0
        macqc_job_progress(handle, &stage, &iteration, &energy, &residual)
        return MacQCProgress(stage: Int(stage), iteration: Int(iteration), energy: energy, residual: residual)
    }

    /// The output file written so far from byte `offset` on, and the total length in bytes
    /// (the next call passes that total as its offset).
    func output(from offset: Int) -> (text: String, total: Int) {
        var total: Int64 = 0
        guard let raw = macqc_job_output(handle, Int64(offset), &total) else { return ("", offset) }
        defer { macqc_free(raw) }
        return (String(cString: raw), Int(total))
    }

    /// One line about the current step of the job ("Optimisation step 4 · E = … · max force …").
    func progressText() -> String {
        guard let raw = macqc_job_progress_text(handle) else { return "" }
        defer { macqc_free(raw) }
        return String(cString: raw)
    }

    fileprivate var pointer: OpaquePointer? { handle }
}

// MARK: - Engine

nonisolated enum MacQCEngineError: LocalizedError {
    case encodingFailed
    case backendReturnedNull
    case decodingFailed(String)
    case orbital(String)
    case structure(String)

    var errorDescription: String? {
        switch self {
        case .encodingFailed:
            return "The input could not be encoded as JSON."
        case .backendReturnedNull:
            return "The backend returned no result (out of memory?)."
        case .decodingFailed(let detail):
            return "The backend output could not be decoded: \(detail)"
        case .orbital(let message):
            return "The orbital could not be drawn: \(message)"
        case .structure(let message):
            return "The structure file could not be read: \(message)"
        }
    }
}

nonisolated enum MacQCEngine {

    /// ΨLab 9.0: reads a structure file (.xyz, .gjf, .com) given as text; `fileName` decides the
    /// format by its extension (otherwise it is recognised from the text). Throws the backend's
    /// explanation for a file that is not a valid structure.
    static func parseStructure(_ text: String, fileName: String) throws -> ParsedStructure {
        guard let raw = macqc_parse_structure(text, fileName) else { throw MacQCEngineError.backendReturnedNull }
        defer { macqc_free(raw) }
        let json = String(cString: raw)
        let parsed: ParsedStructure
        do {
            parsed = try JSONDecoder().decode(ParsedStructure.self, from: Data(json.utf8))
        } catch {
            throw MacQCEngineError.decodingFailed(String(describing: error))
        }
        if let message = parsed.error { throw MacQCEngineError.structure(message) }
        return parsed
    }

    /// Synchronous call into the C++ backend. Blocks the calling thread, so do not call it
    /// from the main thread for long scans; use `run(_:)` instead.
    static func runSync(_ input: MacQCInput, job: MacQCJob? = nil) throws -> MacQCOutput {
        let encoder = JSONEncoder()
        let inputData = try encoder.encode(input)
        guard let inputJSON = String(data: inputData, encoding: .utf8) else {
            throw MacQCEngineError.encodingFailed
        }

        // macqc_run_job never throws and returns a malloc'ed, NUL-terminated UTF-8 string.
        guard let raw = macqc_run_job(inputJSON, job?.pointer) else {
            throw MacQCEngineError.backendReturnedNull
        }
        defer { macqc_free(raw) }
        let outputText = String(cString: raw)

        do {
            return try JSONDecoder().decode(MacQCOutput.self, from: Data(outputText.utf8))
        } catch {
            throw MacQCEngineError.decodingFailed(String(describing: error))
        }
    }

    /// Runs the calculation on a background thread.
    static func run(_ input: MacQCInput, job: MacQCJob? = nil) async throws -> MacQCOutput {
        let task = Task.detached(priority: .userInitiated) {
            try MacQCEngine.runSync(input, job: job)
        }
        return try await task.value
    }

    // MARK: Orbitals (ΨLab 8.0)

    /// The orbitals object as the JSON text the backend reads.
    static func orbitalsJSON(_ orbitals: MacQCOrbitals) throws -> String {
        let data = try JSONEncoder().encode(orbitals)
        guard let text = String(data: data, encoding: .utf8) else { throw MacQCEngineError.encodingFailed }
        return text
    }

    private static func copyMesh(_ m: macqc_mesh) -> OrbitalMesh {
        var mesh = OrbitalMesh()
        let nv = Int(m.vertex_count) * 3, ni = Int(m.triangle_count) * 3
        if let p = m.positions, nv > 0 { mesh.positions = Array(UnsafeBufferPointer<Float>(start: p, count: nv)) }
        if let n = m.normals, nv > 0 { mesh.normals = Array(UnsafeBufferPointer<Float>(start: n, count: nv)) }
        if let i = m.indices, ni > 0 { mesh.indices = Array(UnsafeBufferPointer<Int32>(start: i, count: ni)) }
        return mesh
    }

    /// The isosurfaces psi = +isovalue and psi = -isovalue of an orbital (blocks: call it off
    /// the main thread).
    static func orbitalSurfaces(_ json: String, set setIndex: Int, orbital: Int, isovalue: Double,
                                maxPoints: Int = 400_000) throws -> (positive: OrbitalMesh, negative: OrbitalMesh) {
        var positive = macqc_mesh()
        var negative = macqc_mesh()
        if let error = macqc_orbital_surface(json, Int32(setIndex), Int32(orbital), isovalue, Int32(maxPoints),
                                             &positive, &negative) {
            defer { macqc_free(error) }
            throw MacQCEngineError.orbital(String(cString: error))
        }
        defer {
            macqc_mesh_free(&positive)
            macqc_mesh_free(&negative)
        }
        return (copyMesh(positive), copyMesh(negative))
    }

    /// The orbital at origin + i u + j v (bohr), i < nu fastest, j < nv.
    static func orbitalPlane(_ json: String, set setIndex: Int, orbital: Int, origin: [Double], u: [Double], v: [Double],
                             nu: Int, nv: Int) throws -> [Float] {
        var values = [Float](repeating: 0, count: nu * nv)
        if let error = macqc_orbital_plane(json, Int32(setIndex), Int32(orbital), origin, u, v, Int32(nu), Int32(nv), &values) {
            defer { macqc_free(error) }
            throw MacQCEngineError.orbital(String(cString: error))
        }
        return values
    }
}
