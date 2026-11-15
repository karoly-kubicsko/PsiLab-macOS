//
//  MoleculeDrawing.swift
//  MacQC
//
//  ΨLab 6.0: drawings of the molecule. The Optimisation tab draws the geometry of every
//  optimisation step (bond lengths and the H–O–H angle written on it, the starting geometry
//  dashed) and plays the steps one after the other; the Frequencies tab animates every normal
//  mode. ΨLab 7.0: the Scan tab draws the molecule at every point of a relaxed scan. Planar
//  molecules (H₂, LiH, H₂O, HeH⁺, planar NH₃ and H₂O₂) are drawn in their plane; non-planar ones
//  (H₂O₂, NH₃) in a tilted view along their principal axes, the nearer atoms drawn last and a
//  little larger. Successive geometries are rotated onto each other first, so that a sequence
//  shows how the molecule moves, not how its coordinates turn.
//

import SwiftUI
import AppKit

// MARK: - Geometry helpers

/// Cartesian coordinates (bohr) of the atoms, three per atom.
func atomCoordinates(_ atoms: [MacQCAtom]) -> [[Double]] {
    atoms.map { a in (0..<3).map { $0 < a.position.count ? a.position[$0] : 0 } }
}

/// The atoms with other coordinates (bohr).
func atomsMoved(_ atoms: [MacQCAtom], to positions: [[Double]]) -> [MacQCAtom] {
    atoms.enumerated().map { i, a in
        MacQCAtom(symbol: a.symbol, Z: a.Z, position: i < positions.count ? positions[i] : a.position)
    }
}

/// Index of the heaviest atom (the first of equal ones).
func heaviestAtom(_ atoms: [MacQCAtom]) -> Int? {
    guard !atoms.isEmpty else { return nil }
    var hub = 0
    for (i, a) in atoms.enumerated() where a.Z > atoms[hub].Z { hub = i }
    return hub
}

/// Bonds of the molecules of ΨLab, fixed by the atom numbers (so they stay the same during an
/// optimisation or a scan, whatever the distances): with one atom heavier than hydrogen, every
/// other atom bonded to it (HeH⁺, LiH, H₂O, NH₃); with two (H₂O₂), the bond between them and the
/// hydrogen atoms bonded to them in turn (H3 to O1, H4 to O2); H₂: H–H.
func moleculeBonds(_ atoms: [MacQCAtom]) -> [(Int, Int)] {
    guard atoms.count > 1 else { return [] }
    let heavy = atoms.indices.filter { atoms[$0].Z > 1 }
    if heavy.count >= 2 {
        var bonds: [(Int, Int)] = []
        for a in 0..<heavy.count {
            for b in (a + 1)..<heavy.count {
                bonds.append((heavy[a], heavy[b]))
            }
        }
        let light = atoms.indices.filter { atoms[$0].Z <= 1 }
        for (k, h) in light.enumerated() {
            let partner = heavy[k % heavy.count]
            bonds.append((min(partner, h), max(partner, h)))
        }
        return bonds
    }
    guard let hub = heaviestAtom(atoms) else { return [] }
    return atoms.indices.filter { $0 != hub }.map { (min(hub, $0), max(hub, $0)) }
}

func vectorDistance(_ a: [Double], _ b: [Double]) -> Double {
    var s = 0.0
    for k in 0..<min(a.count, b.count) { s += (a[k] - b[k]) * (a[k] - b[k]) }
    return s.squareRoot()
}

func vectorNorm(_ a: [Double]) -> Double {
    a.reduce(0) { $0 + $1 * $1 }.squareRoot()
}

func vectorDot(_ a: [Double], _ b: [Double]) -> Double {
    var s = 0.0
    for k in 0..<min(a.count, b.count) { s += a[k] * b[k] }
    return s
}

func vectorCross(_ a: [Double], _ b: [Double]) -> [Double] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

/// Positions moved along a displacement: r + factor · d (missing displacements count as zero).
func displacedCoordinates(_ positions: [[Double]], _ displacement: [[Double]], _ factor: Double) -> [[Double]] {
    positions.enumerated().map { i, p in
        let d = i < displacement.count && displacement[i].count == 3 ? displacement[i] : [0, 0, 0]
        return (0..<3).map { p[$0] + factor * d[$0] }
    }
}

/// Eigenvalues and eigenvectors (`vectors[k]` belongs to `values[k]`) of a small symmetric
/// matrix, by cyclic Jacobi rotations.
func symmetricEigen(_ a: [[Double]]) -> (values: [Double], vectors: [[Double]]) {
    let n = a.count
    var m = a
    var v = (0..<n).map { i in (0..<n).map { j in i == j ? 1.0 : 0.0 } }
    for _ in 0..<100 {
        var off = 0.0
        for i in 0..<n {
            for j in (i + 1)..<n { off += m[i][j] * m[i][j] }
        }
        if off < 1e-22 { break }
        for p in 0..<n {
            for q in (p + 1)..<n {
                if abs(m[p][q]) < 1e-300 { continue }
                let theta = (m[q][q] - m[p][p]) / (2 * m[p][q])
                let t = (theta >= 0 ? 1.0 : -1.0) / (abs(theta) + (theta * theta + 1).squareRoot())
                let c = 1 / (t * t + 1).squareRoot()
                let s = t * c
                for k in 0..<n {
                    let mkp = m[k][p], mkq = m[k][q]
                    m[k][p] = c * mkp - s * mkq
                    m[k][q] = s * mkp + c * mkq
                }
                for k in 0..<n {
                    let mpk = m[p][k], mqk = m[q][k]
                    m[p][k] = c * mpk - s * mqk
                    m[q][k] = s * mpk + c * mqk
                }
                for k in 0..<n {
                    let vkp = v[k][p], vkq = v[k][q]
                    v[k][p] = c * vkp - s * vkq
                    v[k][q] = s * vkp + c * vkq
                }
            }
        }
    }
    return ((0..<n).map { m[$0][$0] }, (0..<n).map { k in (0..<n).map { v[$0][k] } })
}

/// `x` rotated and translated onto `y` (least squares with weights, Horn's quaternion method,
/// J. Opt. Soc. Am. A 4, 629 (1987)); a proper rotation, so no mirror images.
func alignedGeometry(_ x: [[Double]], onto y: [[Double]], weights: [Double]) -> [[Double]] {
    let n = x.count
    guard n > 1, y.count == n else { return x }
    let w = (0..<n).map { $0 < weights.count ? weights[$0] : 1.0 }
    let total = w.reduce(0, +)
    guard total > 0 else { return x }
    let cx = (0..<3).map { k in (0..<n).map { w[$0] * x[$0][k] }.reduce(0, +) / total }
    let cy = (0..<3).map { k in (0..<n).map { w[$0] * y[$0][k] }.reduce(0, +) / total }
    var S = [[Double]](repeating: [0, 0, 0], count: 3)
    for i in 0..<n {
        for a in 0..<3 {
            for b in 0..<3 {
                S[a][b] += w[i] * (x[i][a] - cx[a]) * (y[i][b] - cy[b])
            }
        }
    }
    let (xx, xy, xz) = (S[0][0], S[0][1], S[0][2])
    let (yx, yy, yz) = (S[1][0], S[1][1], S[1][2])
    let (zx, zy, zz) = (S[2][0], S[2][1], S[2][2])
    let N: [[Double]] = [[xx + yy + zz, yz - zy, zx - xz, xy - yx],
                         [yz - zy, xx - yy - zz, xy + yx, zx + xz],
                         [zx - xz, xy + yx, -xx + yy - zz, yz + zy],
                         [xy - yx, zx + xz, yz + zy, -xx - yy + zz]]
    let (values, vectors) = symmetricEigen(N)
    var best = 0
    for k in 1..<4 where values[k] > values[best] { best = k }
    let q = vectors[best]
    let (q0, q1, q2, q3) = (q[0], q[1], q[2], q[3])
    let R: [[Double]] = [[q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3, 2 * (q1 * q2 - q0 * q3), 2 * (q1 * q3 + q0 * q2)],
                         [2 * (q1 * q2 + q0 * q3), q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3, 2 * (q2 * q3 - q0 * q1)],
                         [2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1), q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3]]
    return x.map { p in
        let d = (0..<3).map { p[$0] - cx[$0] }
        return (0..<3).map { r in vectorDot(R[r], d) + cy[r] }
    }
}

/// A sequence of geometries of one molecule (optimisation steps, scan points), each rotated onto
/// the one before it with the nuclear charges as weights: the heavy atoms keep their places and
/// the drawing shows how the molecule changes.
func alignedSequence(_ geometries: [[[Double]]], atoms: [MacQCAtom]) -> [[[Double]]] {
    guard var previous = geometries.first else { return [] }
    let weights = atoms.map { Double($0.Z) }
    var result = [previous]
    for g in geometries.dropFirst() {
        let a = alignedGeometry(g, onto: previous, weights: weights)
        result.append(a)
        previous = a
    }
    return result
}

// MARK: - Projection and drawing

/// The view on a set of geometries of one molecule. Planar sets are drawn in their plane: the
/// principal axis along which the atoms spread most horizontal, the second vertical. Non-planar
/// sets in a tilted view: a symmetric top (NH₃) with its unique axis nearly vertical, seen 25°
/// from above; others (H₂O₂) with the longest axis horizontal, seen 30° off the plane of the two
/// longest axes. The heaviest atom is put on top (H₂O is drawn as the usual inverted V); turning
/// the view upside down also turns its depth, so a drawing is never a mirror image.
struct MoleculeFrame {
    var origin: [Double] = [0, 0, 0]
    var axisH: [Double] = [0, 0, 1]
    var axisV: [Double] = [1, 0, 0]
    var axisD: [Double] = [0, 1, 0]
    var flipVertical = false
    var centerH = 0.0
    var centerV = 0.0
    var spreadH = 1.0
    var spreadV = 0.0
    /// Largest distance of an atom from the plane of the drawing (bohr; 0 when planar).
    var depthRange = 0.0

    init(geometries: [[[Double]]], heavyAtom: Int?) {
        let all = geometries.flatMap { $0 }
        guard !all.isEmpty else { return }
        let count = Double(all.count)
        let c = (0..<3).map { k in all.map { $0[k] }.reduce(0, +) / count }
        var cov = [[Double]](repeating: [0, 0, 0], count: 3)
        for p in all {
            for a in 0..<3 {
                for b in 0..<3 {
                    cov[a][b] += (p[a] - c[a]) * (p[b] - c[b]) / count
                }
            }
        }
        let (values, vectors) = symmetricEigen(cov)
        let order = [0, 1, 2].sorted { values[$0] > values[$1] + 1e-12 || (abs(values[$0] - values[$1]) <= 1e-12 && $0 < $1) }
        // Each axis with its largest component positive (LiH: Li on the left, as before).
        var e = order.map { k -> [Double] in
            let v = vectors[k]
            var big = 0
            for i in 1..<3 where abs(v[i]) > abs(v[big]) + 1e-9 { big = i }
            return v[big] < 0 ? v.map { -$0 } : v
        }
        e[2] = vectorCross(e[0], e[1])
        origin = c
        let thickness = all.map { p in vectorDot((0..<3).map { k in p[k] - c[k] }, e[2]) }
        let planar = (thickness.max() ?? 0) - (thickness.min() ?? 0) < 0.05
        axisH = e[0]
        if planar {
            axisV = e[1]
        } else if values[order[1]] > 0.8 * values[order[0]] {
            let b = 25.0 * Double.pi / 180
            axisV = (0..<3).map { cos(b) * e[2][$0] + sin(b) * e[1][$0] }
        } else {
            let a = 30.0 * Double.pi / 180
            axisV = (0..<3).map { cos(a) * e[1][$0] + sin(a) * e[2][$0] }
        }
        axisD = vectorCross(axisH, axisV)
        let ah = axisH, av = axisV
        let h = all.map { p in vectorDot((0..<3).map { k in p[k] - c[k] }, ah) }
        let v = all.map { p in vectorDot((0..<3).map { k in p[k] - c[k] }, av) }
        let lh = h.min() ?? 0, hh = h.max() ?? 0, lv = v.min() ?? 0, hv = v.max() ?? 0
        centerH = (hh + lh) / 2
        centerV = (hv + lv) / 2
        spreadH = hh - lh
        spreadV = hv - lv
        let axis = axisD
        depthRange = all.map { p in abs(vectorDot((0..<3).map { k in p[k] - c[k] }, axis)) }.max() ?? 0
        if planar { depthRange = 0 }
        if let heavy = heavyAtom, let first = geometries.first, heavy < first.count, first.count > 1 {
            let up = first.map { p in vectorDot((0..<3).map { k in p[k] - c[k] }, av) }
            let others = first.indices.filter { $0 != heavy }.map { up[$0] }
            let mean = others.reduce(0, +) / Double(others.count)
            flipVertical = up[heavy] < mean - 1e-6
        }
    }

    private func horizontal(_ p: [Double]) -> Double {
        vectorDot((0..<3).map { p[$0] - origin[$0] }, axisH)
    }

    private func vertical(_ p: [Double]) -> Double {
        vectorDot((0..<3).map { p[$0] - origin[$0] }, axisV)
    }

    /// Distance towards the viewer (bohr).
    func depth(_ p: [Double]) -> Double {
        let d = vectorDot((0..<3).map { p[$0] - origin[$0] }, axisD)
        return flipVertical ? -d : d
    }

    /// Size factor of an atom: up to 12 % larger in front, smaller behind (1 for planar views).
    func depthScale(_ p: [Double]) -> Double {
        guard depthRange > 1e-6 else { return 1 }
        return 1 + 0.12 * depth(p) / depthRange
    }

    /// Pixels per bohr so that the region fits into `size` with `padding` pixels around it.
    func scale(for size: CGSize, padding: Double) -> Double {
        let w = max(Double(size.width) - 2 * padding, 10)
        let h = max(Double(size.height) - 2 * padding, 10)
        var s = w / max(spreadH, 0.5)
        if spreadV > 1e-6 { s = min(s, h / spreadV) }
        return min(s, 150)
    }

    func point(_ p: [Double], size: CGSize, scale: Double) -> CGPoint {
        let x = Double(size.width) / 2 + (horizontal(p) - centerH) * scale
        let dy = (vertical(p) - centerV) * scale
        return CGPoint(x: x, y: Double(size.height) / 2 + (flipVertical ? dy : -dy))
    }
}

/// Colours and sizes of the atoms (CPK colours: H white, He cyan, Li violet, N blue, O red).
enum AtomStyle {
    static func color(_ Z: Int) -> Color {
        switch Z {
        case 1: return Color(white: 0.93)
        case 2: return Color(red: 0.75, green: 0.95, blue: 1.0)
        case 3: return Color(red: 0.75, green: 0.45, blue: 1.0)
        case 7: return Color(red: 0.19, green: 0.31, blue: 0.97)
        case 8: return Color(red: 0.92, green: 0.16, blue: 0.12)
        default: return Color.gray
        }
    }

    static func radius(_ Z: Int) -> Double {
        switch Z {
        case 1, 2: return 11
        case 3: return 17
        default: return 15
        }
    }

    static func labelColor(_ Z: Int) -> Color {
        Z <= 2 ? Color.black : Color.white
    }
}

/// Draws atoms, bonds, outlines, arrows and labels into a Canvas for one MoleculeFrame.
struct MoleculeSketch {
    let frame: MoleculeFrame
    let size: CGSize
    let scale: Double

    init(frame: MoleculeFrame, size: CGSize, padding: Double) {
        self.frame = frame
        self.size = size
        self.scale = frame.scale(for: size, padding: padding)
    }

    func point(_ p: [Double]) -> CGPoint {
        frame.point(p, size: size, scale: scale)
    }

    func bonds(_ gc: GraphicsContext, _ positions: [[Double]], _ bonds: [(Int, Int)], color: Color,
               width: Double, dashed: Bool = false) {
        for (i, j) in bonds where i < positions.count && j < positions.count {
            var path = Path()
            path.move(to: point(positions[i]))
            path.addLine(to: point(positions[j]))
            gc.stroke(path, with: .color(color),
                      style: StrokeStyle(lineWidth: width, lineCap: .round, dash: dashed ? [5, 4] : []))
        }
    }

    /// The atoms, the farthest first (non-planar views: nearer atoms a little larger).
    func atoms(_ gc: GraphicsContext, _ atoms: [MacQCAtom], _ positions: [[Double]]) {
        let order = atoms.indices.filter { $0 < positions.count }
            .sorted { frame.depth(positions[$0]) < frame.depth(positions[$1]) }
        for i in order {
            let atom = atoms[i]
            let p = point(positions[i])
            let r = AtomStyle.radius(atom.Z) * frame.depthScale(positions[i])
            let rect = CGRect(x: Double(p.x) - r, y: Double(p.y) - r, width: 2 * r, height: 2 * r)
            let highlight = CGPoint(x: Double(p.x) - 0.35 * r, y: Double(p.y) - 0.35 * r)
            gc.fill(Path(ellipseIn: rect),
                    with: .radialGradient(Gradient(colors: [Color.white, AtomStyle.color(atom.Z)]),
                                          center: highlight, startRadius: 0, endRadius: 1.3 * r))
            gc.stroke(Path(ellipseIn: rect), with: .color(Color.primary.opacity(0.55)), lineWidth: 1)
            gc.draw(Text(atom.symbol).font(.caption.bold()).foregroundStyle(AtomStyle.labelColor(atom.Z)), at: p)
        }
    }

    /// Dashed circles: a reference geometry (equilibrium, or the start of an optimisation).
    func outlines(_ gc: GraphicsContext, _ atoms: [MacQCAtom], _ positions: [[Double]]) {
        for (i, atom) in atoms.enumerated() where i < positions.count {
            let p = point(positions[i])
            let r = AtomStyle.radius(atom.Z) * frame.depthScale(positions[i])
            let rect = CGRect(x: Double(p.x) - r, y: Double(p.y) - r, width: 2 * r, height: 2 * r)
            gc.stroke(Path(ellipseIn: rect), with: .color(Color.secondary.opacity(0.8)),
                      style: StrokeStyle(lineWidth: 1.2, dash: [4, 3]))
        }
    }

    func arrow(_ gc: GraphicsContext, from a: [Double], to b: [Double], color: Color) {
        let base = point(a)
        let tip = point(b)
        let dx = Double(tip.x - base.x)
        let dy = Double(tip.y - base.y)
        let length = (dx * dx + dy * dy).squareRoot()
        guard length > 6 else { return }
        let ux = dx / length
        let uy = dy / length
        let tx = Double(tip.x)
        let ty = Double(tip.y)
        var path = Path()
        path.move(to: base)
        path.addLine(to: tip)
        path.move(to: CGPoint(x: tx - 9 * ux + 5 * uy, y: ty - 9 * uy - 5 * ux))
        path.addLine(to: tip)
        path.addLine(to: CGPoint(x: tx - 9 * ux - 5 * uy, y: ty - 9 * uy + 5 * ux))
        gc.stroke(path, with: .color(color), style: StrokeStyle(lineWidth: 2.2, lineCap: .round, lineJoin: .round))
    }

    /// A small text on a background, centred at `p`.
    func label(_ gc: GraphicsContext, _ string: String, at p: CGPoint) {
        let text = gc.resolve(Text(string).font(.caption.monospacedDigit()).foregroundStyle(Color.primary))
        let s = text.measure(in: CGSize(width: 300, height: 40))
        let rect = CGRect(x: Double(p.x) - Double(s.width) / 2 - 4, y: Double(p.y) - Double(s.height) / 2 - 2,
                          width: Double(s.width) + 8, height: Double(s.height) + 4)
        gc.fill(Path(roundedRect: rect, cornerRadius: 4),
                with: .color(Color(nsColor: .controlBackgroundColor).opacity(0.85)))
        gc.draw(text, at: p, anchor: .center)
    }
}

// MARK: - Optimisation step

/// The molecule at one optimisation step or scan point, drawn to the scale of the whole sequence
/// (playing the steps shows how it moves), with its bond lengths (Å) and, for H₂O, the H–O–H
/// angle; optionally the starting geometry as dashed outlines.
struct OptimizationGeometryCanvas: View {
    var atoms: [MacQCAtom]
    var start: [MacQCAtom]?
    var frame: MoleculeFrame
    var label = "Molecule at this optimisation step"

    var body: some View {
        Canvas { gc, size in
            draw(gc, size: size)
        }
        .accessibilityLabel(label)
    }

    private func draw(_ gc: GraphicsContext, size: CGSize) {
        guard !atoms.isEmpty else { return }
        let sketch = MoleculeSketch(frame: frame, size: size, padding: 40)
        let positions = atomCoordinates(atoms)
        let bonds = moleculeBonds(atoms)
        if let start, start.count == atoms.count {
            let startPositions = atomCoordinates(start)
            sketch.bonds(gc, startPositions, moleculeBonds(start), color: Color.secondary.opacity(0.7), width: 1.5, dashed: true)
            sketch.outlines(gc, start, startPositions)
        }
        sketch.bonds(gc, positions, bonds, color: Color.secondary, width: 5)
        if bonds.count == 2, let hub = heaviestAtom(atoms) {
            angleMark(gc, sketch, positions, hub: hub, ends: bonds.map { $0.0 == hub ? $0.1 : $0.0 })
        }
        sketch.atoms(gc, atoms, positions)
        // bond lengths beside each bond, on the side away from the centre of the molecule
        let points = positions.map { sketch.point($0) }
        let cx = points.map { Double($0.x) }.reduce(0, +) / Double(points.count)
        let cy = points.map { Double($0.y) }.reduce(0, +) / Double(points.count)
        for (i, j) in bonds {
            let a = points[i], b = points[j]
            let dx = Double(b.x - a.x), dy = Double(b.y - a.y)
            let length = max((dx * dx + dy * dy).squareRoot(), 1)
            var nx = -dy / length, ny = dx / length
            let mx = (Double(a.x) + Double(b.x)) / 2, my = (Double(a.y) + Double(b.y)) / 2
            let side = nx * (cx - mx) + ny * (cy - my)
            if side > 1e-6 || (abs(side) <= 1e-6 && ny > 0) {
                nx = -nx
                ny = -ny
            }
            let r = vectorDistance(positions[i], positions[j]) * Units.bohrToAngstrom
            sketch.label(gc, String(format: "%.4f Å", r), at: CGPoint(x: mx + 18 * nx, y: my + 18 * ny))
        }
    }

    /// Arc and value of the angle at the central atom.
    private func angleMark(_ gc: GraphicsContext, _ sketch: MoleculeSketch, _ positions: [[Double]], hub: Int, ends: [Int]) {
        guard ends.count == 2 else { return }
        let u = (0..<3).map { positions[ends[0]][$0] - positions[hub][$0] }
        let v = (0..<3).map { positions[ends[1]][$0] - positions[hub][$0] }
        let cosine = (0..<3).map { u[$0] * v[$0] }.reduce(0, +) / max(vectorNorm(u) * vectorNorm(v), 1e-12)
        let degrees = acos(min(1, max(-1, cosine))) * 180 / Double.pi
        let c = sketch.point(positions[hub])
        let p1 = sketch.point(positions[ends[0]]), p2 = sketch.point(positions[ends[1]])
        let a1 = atan2(Double(p1.y - c.y), Double(p1.x - c.x))
        let a2 = atan2(Double(p2.y - c.y), Double(p2.x - c.x))
        var diff = a2 - a1
        while diff > Double.pi { diff -= 2 * Double.pi }
        while diff <= -Double.pi { diff += 2 * Double.pi }
        let radius = 30.0
        var arc = Path()
        for k in 0...24 {
            let t = a1 + diff * Double(k) / 24
            let q = CGPoint(x: Double(c.x) + radius * cos(t), y: Double(c.y) + radius * sin(t))
            if k == 0 { arc.move(to: q) } else { arc.addLine(to: q) }
        }
        gc.stroke(arc, with: .color(Color.accentColor), lineWidth: 1.5)
        let mid = a1 + diff / 2
        sketch.label(gc, String(format: "%.2f°", degrees),
                     at: CGPoint(x: Double(c.x) + 52 * cos(mid), y: Double(c.y) + 52 * sin(mid)))
    }
}

// MARK: - Normal mode

/// One normal mode as an animation: the atoms oscillate about the equilibrium geometry (dashed)
/// along the mode's Cartesian displacement, at a speed proportional to the wavenumber;
/// optionally with arrows to the turning points.
struct NormalModeCanvas: View {
    var atoms: [MacQCAtom]
    var displacement: [[Double]]
    var wavenumber: Double
    var animate: Bool
    /// 1: the largest atomic displacement is a quarter of the shortest bond.
    var amplitude: Double
    var arrows: Bool

    /// Seconds per vibration: 2 s at 1600 cm⁻¹, proportionally faster for higher wavenumbers.
    private var period: Double {
        min(2.4, max(0.6, 2.0 * 1600.0 / max(abs(wavenumber), 1.0)))
    }

    var body: some View {
        TimelineView(.animation(minimumInterval: 1.0 / 30.0, paused: !animate)) { context in
            let phase = animate ? sin(context.date.timeIntervalSinceReferenceDate * 2 * Double.pi / period) : 0
            Canvas { gc, size in
                draw(gc, size: size, phase: phase)
            }
        }
    }

    /// Bohr per unit of the normalised displacement.
    private func bohrPerUnit(_ positions: [[Double]]) -> Double {
        let shortest = moleculeBonds(atoms).map { vectorDistance(positions[$0.0], positions[$0.1]) }.min() ?? 1.8
        let largest = displacement.map { vectorNorm($0) }.max() ?? 1
        return amplitude * 0.25 * shortest / max(largest, 1e-6)
    }

    private func draw(_ gc: GraphicsContext, size: CGSize, phase: Double) {
        guard !atoms.isEmpty else { return }
        let equilibrium = atomCoordinates(atoms)
        let a = bohrPerUnit(equilibrium)
        let frame = MoleculeFrame(geometries: [displacedCoordinates(equilibrium, displacement, a),
                                               displacedCoordinates(equilibrium, displacement, -a)],
                                  heavyAtom: heaviestAtom(atoms))
        let sketch = MoleculeSketch(frame: frame, size: size, padding: 28)
        let bonds = moleculeBonds(atoms)
        sketch.bonds(gc, equilibrium, bonds, color: Color.secondary.opacity(0.5), width: 1.5, dashed: true)
        sketch.outlines(gc, atoms, equilibrium)
        if arrows {
            let turning = displacedCoordinates(equilibrium, displacement, a)
            for i in equilibrium.indices {
                sketch.arrow(gc, from: equilibrium[i], to: turning[i], color: Color.orange)
            }
        }
        let now = displacedCoordinates(equilibrium, displacement, a * phase)
        sketch.bonds(gc, now, bonds, color: Color.secondary, width: 5)
        sketch.atoms(gc, atoms, now)
    }
}

// MARK: - Internal coordinates

/// Value of a bond (Å), angle or dihedral angle (degrees; IUPAC sign, positive clockwise seen
/// along the middle bond, as in the backend) for atom numbers from 1; nil if an atom is missing
/// or the dihedral angle is undefined.
func internalCoordinateValue(_ positions: [[Double]], kind: ScanKind, atoms numbers: [Int]) -> Double? {
    let n = numbers.prefix(kind.atomCount).map { $0 - 1 }
    guard n.count == kind.atomCount, n.allSatisfy({ $0 >= 0 && $0 < positions.count }), Set(n).count == n.count else { return nil }
    let p = n.map { positions[$0] }
    func minus(_ a: [Double], _ b: [Double]) -> [Double] { (0..<3).map { a[$0] - b[$0] } }
    switch kind {
    case .bond:
        return vectorDistance(p[0], p[1]) * Units.bohrToAngstrom
    case .angle:
        let u = minus(p[0], p[1]), v = minus(p[2], p[1])
        let c = vectorDot(u, v) / max(vectorNorm(u) * vectorNorm(v), 1e-12)
        return acos(min(1, max(-1, c))) * 180 / Double.pi
    case .dihedral:
        let b1 = minus(p[1], p[0]), b2 = minus(p[2], p[1]), b3 = minus(p[3], p[2])
        let n1 = vectorCross(b1, b2), n2 = vectorCross(b2, b3)
        guard vectorNorm(n1) > 1e-8, vectorNorm(n2) > 1e-8 else { return nil }
        return atan2(vectorNorm(b2) * vectorDot(b1, n2), vectorDot(n1, n2)) * 180 / Double.pi
    }
}
