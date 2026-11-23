//
//  MoleculeDrawing.swift
//  MacQC
//
//  ΨLab 6.0: geometry of the molecule drawings: bonds, coordinates, alignment of the geometries of
//  an optimisation or a scan onto each other (so that a sequence shows how the molecule moves,
//  not how its coordinates turn), internal coordinates. ΨLab 12.0: the molecules are drawn by the
//  3D viewer (Molecule3DView.swift); MoleculeFrame gives its starting orientation, that of the
//  former 2D drawings: planar molecules (H₂, LiH, H₂O, HeH⁺, planar NH₃ and H₂O₂) in their
//  plane, non-planar ones (H₂O₂, NH₃) in a tilted view along their principal axes.
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
    // ΨLab 9.0: custom molecules (any atoms in any order) are bonded by distance
    let builtInPatterns: [[Int]] = [[2, 1], [1, 1], [3, 1], [8, 1, 1], [8, 8, 1, 1], [7, 1, 1, 1]]
    if !builtInPatterns.contains(atoms.map(\.Z)) { return distanceBonds(atoms) }
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

/// Bonds of a custom molecule: atoms closer than 1.25 times the sum of their covalent radii (the
/// backend uses 1.3 for its internal coordinates; 19.0: H–H 1.45 times); every atom keeps at least a bond to its
/// nearest neighbour, so that no atom floats alone in the drawing of a stretched geometry.
func distanceBonds(_ atoms: [MacQCAtom]) -> [(Int, Int)] {
    let p = atomCoordinates(atoms)
    var bonds: [(Int, Int)] = []
    var bonded = Array(repeating: false, count: atoms.count)
    for i in 0..<atoms.count {
        for j in (i + 1)..<atoms.count {
            // 19.0: H–H pairs up to 1.45 times (0.90 Å): the three-centre bonds of H₃⁺ (0.87 Å)
            let factor = atoms[i].Z == 1 && atoms[j].Z == 1 ? 1.45 : 1.25
            let limit = factor * (Elements.covalentRadiusBohr(atoms[i].Z) + Elements.covalentRadiusBohr(atoms[j].Z))
            if vectorDistance(p[i], p[j]) < limit {
                bonds.append((i, j))
                bonded[i] = true
                bonded[j] = true
            }
        }
    }
    for i in 0..<atoms.count where !bonded[i] && atoms.count > 1 {
        var nearest = i == 0 ? 1 : 0
        for j in 0..<atoms.count where j != i && vectorDistance(p[i], p[j]) < vectorDistance(p[i], p[nearest]) {
            nearest = j
        }
        bonds.append((min(i, nearest), max(i, nearest)))
    }
    return bonds
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

// MARK: - Orientation

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
        if planar, let heavy = heavyAtom, let first = geometries.first, heavy < first.count,
           values[order[0]] - values[order[1]] <= 1e-3 * max(values[order[0]], 1e-12) {
            // 19.0: a planar symmetric top (H₃⁺): the two axes in the plane are degenerate, so any
            // pair would do; the heaviest (first) atom straight above the centre, as in the sidebar
            // (H1 at the apex, H2–H3 the horizontal base)
            let d = (0..<3).map { first[heavy][$0] - c[$0] }
            let inPlane = (0..<3).map { d[$0] - vectorDot(d, e[2]) * e[2][$0] }
            let n = vectorNorm(inPlane)
            if n > 1e-3 {
                axisV = inPlane.map { $0 / n }
                axisH = vectorCross(axisV, e[2])
            } else {
                axisV = e[1]
            }
        } else if planar {
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
}

/// Colours and sizes of the atoms (CPK colours: H white, He cyan, Li violet, C grey, N blue, O red,
/// …; 9.0: every element, see ElementData.swift).
enum AtomStyle {
    static func color(_ Z: Int) -> Color {
        let (r, g, b) = Elements.cpk(Z)
        return Color(red: r, green: g, blue: b)
    }

    static func radius(_ Z: Int) -> Double {
        switch Z {
        case 1, 2: return 11
        case 3: return 17
        case 4...10: return 15
        default: return min(19, max(13, 8 + 7 * (Elements.covalentRadiusBohr(Z) * Units.bohrToAngstrom)))
        }
    }

    /// Black on light atoms, white on dark ones.
    static func labelColor(_ Z: Int) -> Color {
        let (r, g, b) = Elements.cpk(Z)
        return 0.299 * r + 0.587 * g + 0.114 * b > 0.6 ? Color.black : Color.white
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
