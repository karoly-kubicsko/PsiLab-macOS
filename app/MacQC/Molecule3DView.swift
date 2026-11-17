//
//  Molecule3DView.swift
//  MacQC
//
//  ΨLab 12.0: a 3D molecule viewer in the manner of PyMOL, rendered with SceneKit. It shows the
//  molecule of the Summary, every optimisation step, every scan point and every normal mode:
//
//  • representations: ball and stick, sticks (licorice), spheres (van der Waals radii) and lines;
//    bonds are split in the middle and coloured by their two atoms (CPK colours);
//  • the mouse turns the molecule about its centre (arcball), a pinch or ⌥-scroll zooms (a plain
//    scroll scrolls the page);
//    Reset view returns to the starting orientation, which is the view of the former 2D drawings
//    (the plane of a planar molecule, the heaviest atom on top);
//  • depth cueing (atoms farther away fade into the background), black, white or window-coloured
//    background, atom labels (element, element and number), bond lengths in Å;
//  • clicking atoms selects up to four of them and measures, like PyMOL's picking: the distance of
//    two (Å), the angle of three, the dihedral angle of four (degrees), with dashed lines;
//  • a translucent reference geometry (the start of an optimisation, the equilibrium of a mode);
//  • normal modes: the atoms oscillate along the mode, r(t) = r₀ + a sin(2πt/T) d, as keyframe
//    animations that SceneKit plays on its render thread (no redraw of the SwiftUI page per
//    frame), with arrows to the turning points;
//  • the view can be saved as a PNG image.
//
//  The scene is in bohr, like the coordinates of the backend.
//

import SwiftUI
import SceneKit
import AppKit
import QuartzCore
import UniformTypeIdentifiers
import simd

// MARK: - Display settings

/// How atoms and bonds are drawn (PyMOL's "show as").
nonisolated enum MoleculeRepresentation: String, CaseIterable, Identifiable, Sendable {
    case ballAndStick, sticks, spheres, lines

    var id: String { rawValue }

    var title: String {
        switch self {
        case .ballAndStick: return "Ball and stick"
        case .sticks: return "Sticks"
        case .spheres: return "Spheres (van der Waals)"
        case .lines: return "Lines"
        }
    }

    /// Radius of an atom (bohr).
    func atomRadius(_ Z: Int) -> Double {
        let angstrom: Double
        switch self {
        case .ballAndStick: angstrom = 0.24 * Elements.vanDerWaalsRadius(Z)
        case .sticks: angstrom = 0.20
        case .spheres: angstrom = Elements.vanDerWaalsRadius(Z)
        case .lines: angstrom = 0.07
        }
        return angstrom / Units.bohrToAngstrom
    }

    /// Radius of a bond (bohr; 0: bonds are not drawn).
    var bondRadius: Double {
        let angstrom: Double
        switch self {
        case .ballAndStick: angstrom = 0.11
        case .sticks: angstrom = 0.20
        case .spheres: angstrom = 0
        case .lines: angstrom = 0.045
        }
        return angstrom / Units.bohrToAngstrom
    }
}

nonisolated enum MoleculeLabelChoice: String, CaseIterable, Identifiable, Sendable {
    case hidden, elements, numbered

    var id: String { rawValue }

    var title: String {
        switch self {
        case .hidden: return "No atom labels"
        case .elements: return "Element symbols"
        case .numbered: return "Symbols and atom numbers"
        }
    }
}

nonisolated enum MoleculeBackgroundChoice: String, CaseIterable, Identifiable, Sendable {
    case automatic, black, white

    var id: String { rawValue }

    var title: String {
        switch self {
        case .automatic: return "Background as the window"
        case .black: return "Black background"
        case .white: return "White background"
        }
    }
}

struct Molecule3DStyle: Equatable {
    var representation: MoleculeRepresentation = .ballAndStick
    var labels: MoleculeLabelChoice = .hidden
    var bondLengths = true
    var depthCue = true
    var darkBackground = false
    /// 0 … 1: the grey of the background.
    var backgroundWhite = 1.0
}

/// A normal mode for the viewer: the atoms move to r₀ ± amplitude · displacement.
struct Molecule3DVibration: Equatable {
    var displacement: [[Double]]
    /// Bohr per unit of the displacement at the turning points.
    var amplitude: Double
    /// Seconds per vibration.
    var period: Double
    var playing: Bool
    var arrows: Bool
}

struct MoleculeBond: Equatable, Hashable {
    var a: Int
    var b: Int
}

/// Everything the scene shows; the controller changes the scene only where this changed.
struct Molecule3DConfiguration: Equatable {
    var Z: [Int] = []
    var symbols: [String] = []
    /// Bohr; for a normal mode, the equilibrium geometry.
    var positions: [[Double]] = []
    var bonds: [MoleculeBond] = []
    /// A translucent second geometry of the same atoms (bohr).
    var reference: [[Double]]?
    /// The geometries the starting view must show (all steps of an optimisation, …).
    var framing: [[[Double]]] = []
    var vibration: Molecule3DVibration?
    var style = Molecule3DStyle()
    /// Picked atoms (indices from 0), in the order of picking.
    var selection: [Int] = []
}

// MARK: - SwiftUI view

/// The 3D viewer with its display menu, Reset view and Save image buttons, and the measurement of
/// the picked atoms. The display settings are kept between launches and shared by all viewers.
struct Molecule3DView: View {
    var atoms: [MacQCAtom]
    /// A translucent second geometry of the same atoms (start of an optimisation, equilibrium).
    var reference: [MacQCAtom]? = nil
    /// Geometries the starting view must contain (default: the atoms).
    var framing: [[[Double]]] = []
    var vibration: Molecule3DVibration? = nil
    var accessibilityText = "Molecule in 3D"

    @State private var controller = Molecule3DController()
    @State private var selection: [Int] = []
    @AppStorage("viewer3D.representation") private var representation: MoleculeRepresentation = .ballAndStick
    @AppStorage("viewer3D.labels") private var labels: MoleculeLabelChoice = .hidden
    @AppStorage("viewer3D.bondLengths") private var bondLengths = true
    @AppStorage("viewer3D.depthCue") private var depthCue = true
    @AppStorage("viewer3D.background") private var background: MoleculeBackgroundChoice = .automatic
    @Environment(\.colorScheme) private var colorScheme

    var body: some View {
        MoleculeSceneRepresentable(controller: controller, configuration: configuration) { picked in
            pick(picked)
        }
        .accessibilityLabel(accessibilityText)
        .overlay(alignment: .topTrailing) {
            toolbar
                .padding(6)
        }
        .overlay(alignment: .bottomLeading) {
            if let text = measurement {
                Text(text)
                    .font(.callout.monospacedDigit())
                    .padding(.horizontal, 8)
                    .padding(.vertical, 4)
                    .background(.regularMaterial, in: RoundedRectangle(cornerRadius: 6))
                    .padding(8)
                    .textSelection(.enabled)
            }
        }
        .clipShape(RoundedRectangle(cornerRadius: 8))
        .overlay(RoundedRectangle(cornerRadius: 8).strokeBorder(Color.secondary.opacity(0.3), lineWidth: 1))
        .onChange(of: atoms.map(\.Z)) {
            selection = []
        }
    }

    private var isDark: Bool {
        switch background {
        case .automatic: return colorScheme == .dark
        case .black: return true
        case .white: return false
        }
    }

    private var configuration: Molecule3DConfiguration {
        let positions = atomCoordinates(atoms)
        var c = Molecule3DConfiguration()
        c.Z = atoms.map(\.Z)
        c.symbols = atoms.map(\.symbol)
        c.positions = positions
        c.bonds = moleculeBonds(atoms).map { MoleculeBond(a: $0.0, b: $0.1) }
        if let reference, reference.count == atoms.count {
            c.reference = atomCoordinates(reference)
        }
        c.framing = framing.isEmpty ? [positions] : framing
        c.vibration = vibration
        c.style = Molecule3DStyle(representation: representation, labels: labels, bondLengths: bondLengths,
                                  depthCue: depthCue, darkBackground: isDark,
                                  backgroundWhite: isDark ? (background == .black ? 0 : 0.11) : 1)
        c.selection = selection.filter { $0 < atoms.count }
        return c
    }

    /// PyMOL-like picking: a click adds an atom (a fifth starts again), a click on a picked atom
    /// removes it, a click beside the molecule clears the selection.
    private func pick(_ index: Int?) {
        guard let i = index else {
            selection = []
            return
        }
        if let k = selection.firstIndex(of: i) {
            selection.remove(at: k)
        } else if selection.count >= 4 {
            selection = [i]
        } else {
            selection.append(i)
        }
    }

    /// The position, distance, angle or dihedral angle of the picked atoms (for a normal mode: at
    /// the equilibrium geometry).
    private var measurement: String? {
        let picked = selection.filter { $0 < atoms.count }
        guard !picked.isEmpty else { return nil }
        let positions = atomCoordinates(atoms)
        let names = picked.map { "\(atoms[$0].symbol)\($0 + 1)" }.joined(separator: "–")
        let numbers = picked.map { $0 + 1 }
        switch picked.count {
        case 1:
            let p = positions[picked[0]].map { $0 * Units.bohrToAngstrom }
            return String(format: "%@  (%.4f, %.4f, %.4f) Å", names, p[0], p[1], p[2])
        case 2:
            return internalCoordinateValue(positions, kind: .bond, atoms: numbers)
                .map { String(format: "Distance %@  %.4f Å", names, $0) }
        case 3:
            return internalCoordinateValue(positions, kind: .angle, atoms: numbers)
                .map { String(format: "Angle %@  %.2f°", names, $0) }
        default:
            return internalCoordinateValue(positions, kind: .dihedral, atoms: numbers)
                .map { String(format: "Dihedral %@  %.2f°", names, $0) } ?? "Dihedral \(names): undefined (three atoms in a line)"
        }
    }

    private var toolbar: some View {
        HStack(spacing: 4) {
            Menu {
                Picker("Representation", selection: $representation) {
                    ForEach(MoleculeRepresentation.allCases) { r in
                        Text(r.title).tag(r)
                    }
                }
                .pickerStyle(.inline)
                Picker("Labels", selection: $labels) {
                    ForEach(MoleculeLabelChoice.allCases) { l in
                        Text(l.title).tag(l)
                    }
                }
                .pickerStyle(.inline)
                Toggle("Bond lengths (Å)", isOn: $bondLengths)
                Toggle("Depth cueing", isOn: $depthCue)
                Picker("Background", selection: $background) {
                    ForEach(MoleculeBackgroundChoice.allCases) { b in
                        Text(b.title).tag(b)
                    }
                }
                .pickerStyle(.inline)
            } label: {
                Image(systemName: "paintpalette")
            }
            .menuStyle(.borderlessButton)
            .menuIndicator(.hidden)
            .fixedSize()
            .help("Representation, labels, bond lengths, depth cueing and background")
            .accessibilityLabel("Display options")
            Button {
                controller.resetView()
            } label: {
                Image(systemName: "scope")
            }
            .buttonStyle(.borderless)
            .help("Reset the view: the starting orientation and size")
            .accessibilityLabel("Reset view")
            Button {
                controller.saveImage()
            } label: {
                Image(systemName: "camera")
            }
            .buttonStyle(.borderless)
            .help("Save the view as a PNG image")
            .accessibilityLabel("Save image")
        }
        .padding(.horizontal, 7)
        .padding(.vertical, 4)
        .background(.regularMaterial, in: RoundedRectangle(cornerRadius: 7))
    }
}

/// The SceneKit view inside SwiftUI.
struct MoleculeSceneRepresentable: NSViewRepresentable {
    let controller: Molecule3DController
    let configuration: Molecule3DConfiguration
    let onPick: (Int?) -> Void

    func makeCoordinator() -> Molecule3DController {
        controller
    }

    func makeNSView(context: Context) -> MoleculeSCNView {
        let view = MoleculeSCNView(frame: NSRect(x: 0, y: 0, width: 400, height: 260), options: nil)
        controller.attach(view)
        controller.onPick = onPick
        controller.apply(configuration)
        return view
    }

    func updateNSView(_ view: MoleculeSCNView, context: Context) {
        controller.onPick = onPick
        controller.apply(configuration)
    }

    static func dismantleNSView(_ view: MoleculeSCNView, coordinator: Molecule3DController) {
        coordinator.detach()
    }
}

/// SCNView that reports clicks (not drags) for picking atoms. The viewers sit on scrolling
/// pages: a plain scroll scrolls the page, ⌥-scroll (or a pinch) zooms the molecule.
final class MoleculeSCNView: SCNView {
    var onClick: ((CGPoint) -> Void)?

    override func scrollWheel(with event: NSEvent) {
        if event.modifierFlags.contains(.option) || event.modifierFlags.contains(.command) {
            super.scrollWheel(with: event)
        } else {
            nextResponder?.scrollWheel(with: event)
        }
    }

    func installClickRecognizer() {
        let recognizer = NSClickGestureRecognizer(target: self, action: #selector(clicked(_:)))
        recognizer.delaysPrimaryMouseButtonEvents = false
        addGestureRecognizer(recognizer)
    }

    @objc private func clicked(_ recognizer: NSClickGestureRecognizer) {
        guard recognizer.state == .ended else { return }
        onClick?(recognizer.location(in: self))
    }
}

// MARK: - Scene

/// Builds and updates the SceneKit scene of one viewer.
final class Molecule3DController {
    var onPick: ((Int?) -> Void)?

    private weak var view: MoleculeSCNView?
    private var scene: SCNScene?
    private let cameraNode = SCNNode()
    private let moleculeRoot = SCNNode()
    private let ghostRoot = SCNNode()
    private let arrowRoot = SCNNode()
    private let annotationRoot = SCNNode()
    private var atomNodes: [SCNNode] = []
    private var bondNodes: [(MoleculeBond, SCNNode)] = []
    private var last: Molecule3DConfiguration?
    private var cameraKey: [Double] = []
    /// Centre and radius (bohr) of the region the starting view shows.
    private var center = SIMD3<Float>(0, 0, 0)
    private var radius: Float = 4
    private var depthTimer: Timer?
    private var depthCueDistance: Float = -1
    private var materials: [String: SCNMaterial] = [:]
    private var labelImages: [String: CGImage] = [:]

    // MARK: Life cycle

    func attach(_ view: MoleculeSCNView) {
        // (again after SwiftUI made a new view: start from an empty scene)
        depthTimer?.invalidate()
        last = nil
        cameraKey = []
        cameraNode.childNodes.forEach { $0.removeFromParentNode() }
        for node in [cameraNode, moleculeRoot, ghostRoot, arrowRoot, annotationRoot] {
            node.removeFromParentNode()
        }
        let scene = SCNScene()
        self.scene = scene
        self.view = view
        cameraNode.camera = {
            let camera = SCNCamera()
            camera.fieldOfView = 30
            camera.zNear = 0.05
            camera.zFar = 2000
            return camera
        }()
        // a key light from the upper left and a weak fill light from the lower right, both turning
        // with the camera (the lighting stays the same while the molecule is turned), and ambient
        let key = SCNNode()
        key.light = SCNLight()
        key.light?.type = .directional
        key.light?.intensity = 1150
        key.eulerAngles = SCNVector3(x: -0.45, y: -0.5, z: 0)
        cameraNode.addChildNode(key)
        let fill = SCNNode()
        fill.light = SCNLight()
        fill.light?.type = .directional
        fill.light?.intensity = 320
        fill.eulerAngles = SCNVector3(x: 0.5, y: 0.6, z: 0)
        cameraNode.addChildNode(fill)
        let ambient = SCNNode()
        ambient.light = SCNLight()
        ambient.light?.type = .ambient
        ambient.light?.intensity = 420
        scene.rootNode.addChildNode(ambient)
        scene.rootNode.addChildNode(cameraNode)
        ghostRoot.opacity = 0.28
        for node in [moleculeRoot, ghostRoot, arrowRoot, annotationRoot] {
            scene.rootNode.addChildNode(node)
        }
        view.scene = scene
        view.pointOfView = cameraNode
        view.allowsCameraControl = true
        view.defaultCameraController.interactionMode = .orbitArcball
        view.defaultCameraController.inertiaEnabled = true
        view.autoenablesDefaultLighting = false
        view.antialiasingMode = .multisampling4X
        view.preferredFramesPerSecond = 60
        view.rendersContinuously = false
        view.installClickRecognizer()
        view.onClick = { [weak self] point in
            self?.click(at: point)
        }
        // Depth cueing follows the distance of the camera (zooming changes it): checked 10 times
        // a second, also while the mouse is dragged.
        let timer = Timer(timeInterval: 0.1, repeats: true) { [weak self] _ in
            MainActor.assumeIsolated {
                self?.updateDepthCue(force: false)
            }
        }
        RunLoop.main.add(timer, forMode: .common)
        depthTimer = timer
    }

    func detach() {
        depthTimer?.invalidate()
        depthTimer = nil
        view?.isPlaying = false
    }

    // MARK: Updates

    func apply(_ c: Molecule3DConfiguration) {
        guard view != nil, scene != nil, c != last else { return }
        let previous = last
        last = c
        let structureChanged = previous.map {
            $0.Z != c.Z || $0.symbols != c.symbols || $0.bonds != c.bonds || $0.style.representation != c.style.representation
                || $0.style.labels != c.style.labels || $0.style.darkBackground != c.style.darkBackground
        } ?? true
        let positionsChanged = previous?.positions != c.positions
        // a new optimisation step or scan point: the atoms glide there
        let glide = !structureChanged && positionsChanged && c.vibration == nil && previous?.vibration == nil

        SCNTransaction.begin()
        SCNTransaction.animationDuration = glide ? 0.3 : 0
        if structureChanged {
            buildMolecule(c)
        }
        placeMolecule(c)
        SCNTransaction.commit()

        SCNTransaction.begin()
        SCNTransaction.animationDuration = 0
        if structureChanged || previous?.reference != c.reference {
            buildGhost(c)
        }
        if structureChanged || positionsChanged || previous?.vibration != c.vibration {
            buildVibration(c)
        }
        if structureChanged || positionsChanged || previous?.selection != c.selection
            || previous?.style.bondLengths != c.style.bondLengths || (previous?.vibration == nil) != (c.vibration == nil) {
            buildAnnotations(c)
        }
        applyBackground(c)
        let key = Molecule3DController.cameraKey(c)
        if key != cameraKey {
            cameraKey = key
            resetView()
        }
        SCNTransaction.commit()
        updateDepthCue(force: true)
    }

    /// The starting view changes only with another molecule or another first geometry (the start of
    /// an optimisation, the equilibrium of the modes), not with the step or mode shown.
    private static func cameraKey(_ c: Molecule3DConfiguration) -> [Double] {
        let first = c.framing.first ?? c.positions
        return c.Z.map { Double($0) } + first.flatMap { $0.map { ($0 * 1000).rounded() } }
    }

    // MARK: Atoms and bonds

    private func buildMolecule(_ c: Molecule3DConfiguration) {
        moleculeRoot.childNodes.forEach { $0.removeFromParentNode() }
        atomNodes = []
        bondNodes = []
        let rep = c.style.representation
        for (i, Z) in c.Z.enumerated() {
            let node = SCNNode(geometry: sphere(Z, radius: rep.atomRadius(Z)))
            node.name = "atom:\(i)"
            if c.style.labels != .hidden {
                let text = c.style.labels == .numbered ? "\(c.symbols[i])\(i + 1)" : c.symbols[i]
                node.addChildNode(label(text, height: 0.55, lift: 0, toward: Float(rep.atomRadius(Z)) + 0.05,
                                        dark: c.style.darkBackground, plate: false))
            }
            moleculeRoot.addChildNode(node)
            atomNodes.append(node)
        }
        if rep.bondRadius > 0 {
            for bond in c.bonds where bond.a < c.Z.count && bond.b < c.Z.count {
                let node = SCNNode()
                let first = SCNNode(geometry: halfCylinder(c.Z[bond.a], radius: rep.bondRadius))
                first.simdPosition = SIMD3<Float>(0, 0.25, 0)
                let second = SCNNode(geometry: halfCylinder(c.Z[bond.b], radius: rep.bondRadius))
                second.simdPosition = SIMD3<Float>(0, 0.75, 0)
                node.addChildNode(first)
                node.addChildNode(second)
                moleculeRoot.addChildNode(node)
                bondNodes.append((bond, node))
            }
        }
    }

    private func placeMolecule(_ c: Molecule3DConfiguration) {
        for (i, node) in atomNodes.enumerated() where i < c.positions.count {
            node.simdPosition = vector(c.positions[i])
        }
        for (bond, node) in bondNodes where bond.a < c.positions.count && bond.b < c.positions.count {
            let (q, length) = Molecule3DController.bondFrame(vector(c.positions[bond.a]), vector(c.positions[bond.b]))
            node.simdPosition = vector(c.positions[bond.a])
            node.simdOrientation = q
            node.simdScale = SIMD3<Float>(1, max(length, 1e-4), 1)
        }
    }

    /// The translucent second geometry (same atoms and bonds).
    private func buildGhost(_ c: Molecule3DConfiguration) {
        ghostRoot.childNodes.forEach { $0.removeFromParentNode() }
        guard let reference = c.reference, reference.count == c.Z.count else { return }
        let rep = c.style.representation
        for (i, Z) in c.Z.enumerated() {
            let node = SCNNode(geometry: sphere(Z, radius: rep.atomRadius(Z)))
            node.simdPosition = vector(reference[i])
            ghostRoot.addChildNode(node)
        }
        guard rep.bondRadius > 0 else { return }
        for bond in c.bonds where bond.a < c.Z.count && bond.b < c.Z.count {
            let a = vector(reference[bond.a]), b = vector(reference[bond.b])
            let (q, length) = Molecule3DController.bondFrame(a, b)
            let node = SCNNode()
            node.simdPosition = a
            node.simdOrientation = q
            node.simdScale = SIMD3<Float>(1, max(length, 1e-4), 1)
            let first = SCNNode(geometry: halfCylinder(c.Z[bond.a], radius: rep.bondRadius))
            first.simdPosition = SIMD3<Float>(0, 0.25, 0)
            let second = SCNNode(geometry: halfCylinder(c.Z[bond.b], radius: rep.bondRadius))
            second.simdPosition = SIMD3<Float>(0, 0.75, 0)
            node.addChildNode(first)
            node.addChildNode(second)
            ghostRoot.addChildNode(node)
        }
    }

    // MARK: Normal modes

    /// Keyframe animations of the atoms and bonds along the mode (60 samples of one period, played
    /// in a loop by SceneKit), and the arrows from the atoms towards the turning points.
    private func buildVibration(_ c: Molecule3DConfiguration) {
        atomNodes.forEach { $0.removeAnimation(forKey: "vibration") }
        bondNodes.forEach { $0.1.removeAnimation(forKey: "vibration") }
        arrowRoot.childNodes.forEach { $0.removeFromParentNode() }
        guard let v = c.vibration, v.displacement.count == c.positions.count else {
            view?.isPlaying = false
            view?.rendersContinuously = false
            return
        }
        let samples = 60
        let times = (0...samples).map { NSNumber(value: Double($0) / Double(samples)) }
        let frames: [[[Double]]] = (0...samples).map { k in
            displacedCoordinates(c.positions, v.displacement, v.amplitude * sin(2 * Double.pi * Double(k) / Double(samples)))
        }
        if v.playing {
            for (i, node) in atomNodes.enumerated() where i < c.positions.count {
                let animation = CAKeyframeAnimation(keyPath: "position")
                animation.values = frames.map { NSValue(scnVector3: scnVector(vector($0[i]))) }
                animation.keyTimes = times
                animation.duration = v.period
                animation.repeatCount = .infinity
                node.addAnimation(animation, forKey: "vibration")
            }
            for (bond, node) in bondNodes {
                let animation = CAKeyframeAnimation(keyPath: "transform")
                animation.values = frames.map { f in
                    NSValue(scnMatrix4: SCNMatrix4FromMat4(Molecule3DController.bondMatrix(vector(f[bond.a]), vector(f[bond.b]))))
                }
                animation.keyTimes = times
                animation.duration = v.period
                animation.repeatCount = .infinity
                node.addAnimation(animation, forKey: "vibration")
            }
        }
        view?.isPlaying = v.playing
        view?.rendersContinuously = v.playing
        guard v.arrows else { return }
        let rep = c.style.representation
        let shaft = Float(max(0.05, 0.35 * max(rep.bondRadius, 0.12)))
        for (i, d) in v.displacement.enumerated() where i < c.positions.count {
            let move = vector(d) * Float(v.amplitude)
            let length = simd_length(move)
            guard length > 0.04 else { continue }
            let direction = move / length
            let start = vector(c.positions[i]) + direction * Float(rep.atomRadius(c.Z[i]) * 0.9)
            arrowRoot.addChildNode(arrow(from: start, direction: direction, length: max(length, 0.3), radius: shaft))
        }
    }

    // MARK: Labels, selection and measurements

    private func buildAnnotations(_ c: Molecule3DConfiguration) {
        annotationRoot.childNodes.forEach { $0.removeFromParentNode() }
        for node in atomNodes {
            node.childNodes.filter { $0.name == "halo" }.forEach { $0.removeFromParentNode() }
        }
        let rep = c.style.representation
        // bond lengths (not for a vibrating mode), a little above the middle of each bond
        if c.style.bondLengths && c.vibration == nil {
            for bond in c.bonds where bond.a < c.positions.count && bond.b < c.positions.count {
                let a = vector(c.positions[bond.a]), b = vector(c.positions[bond.b])
                let r = Double(simd_length(b - a)) * Units.bohrToAngstrom
                let node = label(String(format: "%.4f Å", r), height: 0.42, lift: 0.42,
                                 toward: Float(max(rep.bondRadius, 0.1)) + 0.05, dark: c.style.darkBackground, plate: true)
                node.simdPosition = (a + b) / 2
                annotationRoot.addChildNode(node)
            }
        }
        // picked atoms: a pink shell (PyMOL's pk1 … pk4) and yellow dashes between them
        let picked = c.selection.filter { $0 < atomNodes.count }
        for i in picked {
            let r = rep.atomRadius(c.Z[i]) * 1.25 + 0.18
            let halo = SCNNode(geometry: haloSphere(radius: r))
            halo.name = "halo"
            atomNodes[i].addChildNode(halo)
        }
        if picked.count >= 2 {
            for k in 0..<(picked.count - 1) {
                annotationRoot.addChildNode(dashes(from: vector(c.positions[picked[k]]), to: vector(c.positions[picked[k + 1]])))
            }
        }
    }

    private func click(at point: CGPoint) {
        guard let view else { return }
        let hits = view.hitTest(point, options: [SCNHitTestOption.ignoreHiddenNodes: true])
        for hit in hits {
            var node: SCNNode? = hit.node
            while let n = node {
                if let name = n.name, name.hasPrefix("atom:"), let i = Int(name.dropFirst(5)) {
                    onPick?(i)
                    return
                }
                node = n.parent
            }
        }
        onPick?(nil)
    }

    // MARK: Camera, background, depth cueing

    /// The starting view: the orientation of the former 2D drawings (MoleculeFrame: the plane of a
    /// planar molecule, the heaviest atom on top), at a distance that shows all framing geometries.
    func resetView() {
        guard let c = last, let view else { return }
        let geometries = c.framing.isEmpty ? [c.positions] : c.framing
        guard !geometries.flatMap({ $0 }).isEmpty else { return }
        var heavy: Int?
        for (i, Z) in c.Z.enumerated() where heavy == nil || Z > c.Z[heavy!] { heavy = i }
        let frame = MoleculeFrame(geometries: geometries, heavyAtom: heavy)
        let right = vector(frame.axisH)
        let up = vector(frame.flipVertical ? frame.axisV.map { -$0 } : frame.axisV)
        let toward = vector(frame.flipVertical ? frame.axisD.map { -$0 } : frame.axisD)
        let middle = (0..<3).map { frame.origin[$0] + frame.centerH * frame.axisH[$0] + frame.centerV * frame.axisV[$0] }
        center = vector(middle)
        let largestAtom = Float(c.Z.map { c.style.representation.atomRadius($0) }.max() ?? 0.6)
        var halfWidth: Float = 0.5, halfHeight: Float = 0.5, front: Float = 0, reach: Float = 0.5
        for p in geometries.flatMap({ $0 }) {
            let d = vector(p) - center
            halfWidth = max(halfWidth, abs(simd_dot(d, right)))
            halfHeight = max(halfHeight, abs(simd_dot(d, up)))
            front = max(front, simd_dot(d, toward))
            reach = max(reach, simd_length(d))
        }
        let margin = largestAtom + 0.6
        halfWidth += margin
        halfHeight += margin
        radius = reach + largestAtom
        let bounds = view.bounds.size
        let aspect = bounds.width > 10 && bounds.height > 10 ? Float(bounds.width / bounds.height) : 1.8
        let t = Float(tan(15.0 * Double.pi / 180))
        let distance = 1.08 * max(halfHeight / t, halfWidth / (t * min(max(aspect, 0.7), 3))) + front + largestAtom
        cameraNode.simdPosition = center + distance * toward
        cameraNode.simdLook(at: center, up: up, localFront: SIMD3<Float>(0, 0, -1))
        view.defaultCameraController.target = scnVector(center)
        depthCueDistance = -1
        updateDepthCue(force: true)
    }

    private func applyBackground(_ c: Molecule3DConfiguration) {
        let color = NSColor(white: CGFloat(c.style.backgroundWhite), alpha: 1)
        view?.backgroundColor = color
        scene?.background.contents = color
        scene?.fogColor = color
    }

    /// Linear fog from just in front of the molecule's centre to well behind it: the farthest atoms
    /// are about half faded into the background (PyMOL's depth_cue).
    private func updateDepthCue(force: Bool) {
        guard let scene, let c = last else { return }
        guard c.style.depthCue else {
            if scene.fogEndDistance != 0 {
                scene.fogStartDistance = 0
                scene.fogEndDistance = 0
            }
            return
        }
        let distance = simd_length(cameraNode.simdWorldPosition - center)
        guard force || abs(distance - depthCueDistance) > 0.01 * max(radius, 1) else { return }
        depthCueDistance = distance
        scene.fogDensityExponent = 1
        scene.fogStartDistance = CGFloat(max(distance - 0.25 * radius, 0.01))
        scene.fogEndDistance = CGFloat(distance + 2.2 * radius)
    }

    // MARK: Image

    /// Saves what the viewer shows as a PNG image (at the resolution of the screen).
    func saveImage() {
        guard let view else { return }
        let image = view.snapshot()
        guard let cg = image.cgImage(forProposedRect: nil, context: nil, hints: nil),
              let data = NSBitmapImageRep(cgImage: cg).representation(using: .png, properties: [:]) else { return }
        let panel = NSSavePanel()
        panel.allowedContentTypes = [.png]
        panel.nameFieldStringValue = "molecule.png"
        guard panel.runModal() == .OK, let url = panel.url else { return }
        try? data.write(to: url)
    }

    // MARK: Geometry

    private func material(_ Z: Int) -> SCNMaterial {
        let key = "Z\(Z)"
        if let m = materials[key] { return m }
        let (r, g, b) = Elements.cpk(Z)
        let m = SCNMaterial()
        m.lightingModel = .physicallyBased
        m.diffuse.contents = NSColor(srgbRed: CGFloat(r), green: CGFloat(g), blue: CGFloat(b), alpha: 1)
        m.metalness.contents = 0.0
        m.roughness.contents = 0.38
        materials[key] = m
        return m
    }

    private func sphere(_ Z: Int, radius: Double) -> SCNSphere {
        let s = SCNSphere(radius: CGFloat(radius))
        s.segmentCount = 48
        s.materials = [material(Z)]
        return s
    }

    /// Half a bond: radius r, height 0.5 along y (the bond node scales y to the bond length).
    private func halfCylinder(_ Z: Int, radius: Double) -> SCNCylinder {
        let c = SCNCylinder(radius: CGFloat(radius), height: 0.5)
        c.radialSegmentCount = 28
        c.materials = [material(Z)]
        return c
    }

    private func haloSphere(radius: Double) -> SCNSphere {
        let s = SCNSphere(radius: CGFloat(radius))
        s.segmentCount = 32
        let m = SCNMaterial()
        m.lightingModel = .constant
        m.diffuse.contents = NSColor(srgbRed: 1.0, green: 0.42, blue: 0.75, alpha: 1)
        m.transparency = 0.38
        m.writesToDepthBuffer = false
        s.materials = [m]
        return s
    }

    /// Yellow dashes between two picked atoms (PyMOL's distance object).
    private func dashes(from a: SIMD3<Float>, to b: SIMD3<Float>) -> SCNNode {
        let root = SCNNode()
        let d = b - a
        let length = simd_length(d)
        guard length > 1e-3 else { return root }
        let u = d / length
        let q = Molecule3DController.rotationFromY(to: u)
        let m = SCNMaterial()
        m.lightingModel = .constant
        m.diffuse.contents = NSColor(srgbRed: 1.0, green: 0.86, blue: 0.1, alpha: 1)
        let dash: Float = 0.22, gap: Float = 0.14
        var s: Float = 0
        while s < length {
            let e = min(s + dash, length)
            let cylinder = SCNCylinder(radius: 0.045, height: CGFloat(e - s))
            cylinder.radialSegmentCount = 10
            cylinder.materials = [m]
            let node = SCNNode(geometry: cylinder)
            node.simdOrientation = q
            node.simdPosition = a + u * ((s + e) / 2)
            root.addChildNode(node)
            s += dash + gap
        }
        return root
    }

    private func arrow(from start: SIMD3<Float>, direction u: SIMD3<Float>, length: Float, radius: Float) -> SCNNode {
        let m = SCNMaterial()
        m.lightingModel = .physicallyBased
        m.diffuse.contents = NSColor.systemOrange
        m.roughness.contents = 0.45
        let head = min(0.45 * length, 7 * radius)
        let node = SCNNode()
        node.simdPosition = start
        node.simdOrientation = Molecule3DController.rotationFromY(to: u)
        let shaft = SCNCylinder(radius: CGFloat(radius), height: CGFloat(max(length - head, 0.01)))
        shaft.radialSegmentCount = 16
        shaft.materials = [m]
        let shaftNode = SCNNode(geometry: shaft)
        shaftNode.simdPosition = SIMD3<Float>(0, (length - head) / 2, 0)
        let cone = SCNCone(topRadius: 0, bottomRadius: CGFloat(2.6 * radius), height: CGFloat(head))
        cone.radialSegmentCount = 20
        cone.materials = [m]
        let coneNode = SCNNode(geometry: cone)
        coneNode.simdPosition = SIMD3<Float>(0, length - head / 2, 0)
        node.addChildNode(shaftNode)
        node.addChildNode(coneNode)
        return node
    }

    /// A text that always faces the viewer: `height` bohr high, `lift` bohr above its point on the
    /// screen and `toward` bohr nearer to the viewer (in front of its own atom or bond). With a
    /// plate: on a translucent rounded rectangle and drawn over everything (bond lengths).
    private func label(_ text: String, height: Float, lift: Float, toward: Float, dark: Bool, plate: Bool) -> SCNNode {
        let holder = SCNNode()
        let billboard = SCNBillboardConstraint()
        billboard.freeAxes = .all
        holder.constraints = [billboard]
        guard let image = labelImage(text, dark: dark, plate: plate) else { return holder }
        let width = height * Float(image.width) / Float(max(image.height, 1))
        let plane = SCNPlane(width: CGFloat(width), height: CGFloat(height))
        let m = SCNMaterial()
        m.lightingModel = .constant
        m.diffuse.contents = image
        m.isDoubleSided = true
        m.writesToDepthBuffer = false
        if plate {
            m.readsFromDepthBuffer = false
        }
        plane.materials = [m]
        let node = SCNNode(geometry: plane)
        node.simdPosition = SIMD3<Float>(0, lift, toward)
        node.renderingOrder = plate ? 2000 : 1000
        holder.addChildNode(node)
        return holder
    }

    /// The text of a label as an image (26-point text, sharp when zoomed in); cached.
    private func labelImage(_ text: String, dark: Bool, plate: Bool) -> CGImage? {
        let key = "\(text)|\(dark)|\(plate)"
        if let image = labelImages[key] { return image }
        let font = NSFont.systemFont(ofSize: 26, weight: .semibold)
        let foreground = dark ? NSColor.white : NSColor.black
        let string = NSAttributedString(string: text, attributes: [.font: font, .foregroundColor: foreground])
        let size = string.size()
        let padX: CGFloat = plate ? 10 : 4, padY: CGFloat = plate ? 4 : 2
        let width = Int((size.width + 2 * padX).rounded(.up)), height = Int((size.height + 2 * padY).rounded(.up))
        guard let space = CGColorSpace(name: CGColorSpace.sRGB),
              let context = CGContext(data: nil, width: width, height: height, bitsPerComponent: 8, bytesPerRow: 0,
                                      space: space, bitmapInfo: CGImageAlphaInfo.premultipliedLast.rawValue) else { return nil }
        let graphics = NSGraphicsContext(cgContext: context, flipped: false)
        NSGraphicsContext.saveGraphicsState()
        NSGraphicsContext.current = graphics
        let rect = NSRect(x: 0, y: 0, width: width, height: height)
        if plate {
            (dark ? NSColor(white: 0.12, alpha: 0.78) : NSColor(white: 1, alpha: 0.82)).setFill()
            NSBezierPath(roundedRect: rect, xRadius: 9, yRadius: 9).fill()
        } else {
            // a soft outline in the background colour keeps the symbol readable on any atom
            let outline = NSAttributedString(string: text, attributes: [
                .font: font, .strokeColor: dark ? NSColor.black : NSColor.white, .strokeWidth: 9.0,
            ])
            outline.draw(at: NSPoint(x: padX, y: padY))
        }
        string.draw(at: NSPoint(x: padX, y: padY))
        NSGraphicsContext.restoreGraphicsState()
        let image = context.makeImage()
        if labelImages.count > 400 { labelImages.removeAll() }
        labelImages[key] = image
        return image
    }

    // MARK: Vector helpers

    private func vector(_ p: [Double]) -> SIMD3<Float> {
        SIMD3<Float>(Float(p.count > 0 ? p[0] : 0), Float(p.count > 1 ? p[1] : 0), Float(p.count > 2 ? p[2] : 0))
    }

    private func scnVector(_ v: SIMD3<Float>) -> SCNVector3 {
        SCNVector3(x: CGFloat(v.x), y: CGFloat(v.y), z: CGFloat(v.z))
    }

    /// The rotation that turns the y axis onto the unit vector u.
    nonisolated static func rotationFromY(to u: SIMD3<Float>) -> simd_quatf {
        let y = SIMD3<Float>(0, 1, 0)
        if simd_dot(y, u) < -0.999_999 {
            return simd_quatf(angle: Float.pi, axis: SIMD3<Float>(1, 0, 0))
        }
        return simd_quatf(from: y, to: u)
    }

    /// Orientation and length of a bond node from a to b.
    nonisolated static func bondFrame(_ a: SIMD3<Float>, _ b: SIMD3<Float>) -> (simd_quatf, Float) {
        let d = b - a
        let length = simd_length(d)
        guard length > 1e-6 else { return (simd_quatf(angle: 0, axis: SIMD3<Float>(0, 1, 0)), 0) }
        return (rotationFromY(to: d / length), length)
    }

    /// The transform of a bond node: at a, y axis along the bond, scaled to its length.
    nonisolated static func bondMatrix(_ a: SIMD3<Float>, _ b: SIMD3<Float>) -> simd_float4x4 {
        let (q, length) = bondFrame(a, b)
        var m = simd_float4x4(q)
        m.columns.1 = m.columns.1 * max(length, 1e-4)
        m.columns.3 = SIMD4<Float>(a.x, a.y, a.z, 1)
        return m
    }
}
