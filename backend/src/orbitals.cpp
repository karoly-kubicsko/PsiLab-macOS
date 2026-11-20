// Molecular orbitals for display (ΨLab 8.0): orbital sets, values on grids, isosurfaces.
#include "macqc/orbitals.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>

#include "macqc/ci.h"
#include "macqc/gto_integrals.h"
#include "macqc/parallel.h"
#include "macqc/rhf.h"

namespace macqc {

namespace {

constexpr double kPiOrbitals = 3.14159265358979323846;

std::string frontier_label(int k, int homo) {
    // k, homo: 0-based orbital numbers
    if (k == homo) return "HOMO";
    if (k == homo + 1) return "LUMO";
    if (k < homo) return "HOMO-" + std::to_string(homo - k);
    return "LUMO+" + std::to_string(k - homo - 1);
}

// x^a y^b z^c for a + b + c = l in the order of solid_harmonic_matrix: (l,0,0), (l-1,1,0),
// (l-1,0,1), (l-2,2,0), ... (0,0,l).
void cartesian_monomials(int l, double x, double y, double z, double* out) {
    double px[kMaxAngularMomentum + 1], py[kMaxAngularMomentum + 1], pz[kMaxAngularMomentum + 1];
    px[0] = py[0] = pz[0] = 1.0;
    for (int k = 1; k <= l; ++k) {
        px[k] = px[k - 1] * x;
        py[k] = py[k - 1] * y;
        pz[k] = pz[k - 1] * z;
    }
    int idx = 0;
    for (int i = 0; i <= l; ++i)          // i = l - (x exponent)
        for (int c = 0; c <= i; ++c) {    // z exponent
            const int a = l - i, b = i - c;
            out[idx++] = px[a] * py[b] * pz[c];
        }
}

}  // namespace

OrbitalData make_orbital_data(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf,
                              const FCIResultGeneral* fci, int virtuals) {
    OrbitalData d;
    d.atoms = mol.atoms;
    d.basis_name = basis.name;
    d.n_basis = basis.n_functions;
    d.shells = basis.shells;
    const std::size_t nbf = rhf.C.rows(), nmo = rhf.C.cols();
    if (nbf != static_cast<std::size_t>(basis.n_functions) || nmo == 0) return d;
    const int nocc = rhf.n_occupied;
    const int last = std::min<int>(static_cast<int>(nmo), nocc + std::max(0, virtuals));
    OrbitalSet can;
    can.kind = "canonical";
    can.title = "RHF canonical orbitals";
    for (int k = 0; k < last; ++k) {
        OrbitalInfo o;
        o.number = k + 1;
        o.energy = k < static_cast<int>(rhf.orbital_energies.size()) ? rhf.orbital_energies[k] : 0.0;
        o.has_energy = k < static_cast<int>(rhf.orbital_energies.size());
        o.occupation = k < nocc ? 2.0 : 0.0;
        o.label = frontier_label(k, nocc - 1);
        can.orbitals.push_back(o);
        std::vector<double> c(nbf);
        for (std::size_t mu = 0; mu < nbf; ++mu) c[mu] = rhf.C(mu, static_cast<std::size_t>(k));
        can.coefficients.push_back(std::move(c));
    }
    d.sets.push_back(std::move(can));
    if (fci && fci->natural_orbitals.rows() == nmo && fci->natural_orbitals.cols() > 0) {
        const Matrix& U = fci->natural_orbitals;  // MO x NO
        OrbitalSet nat;
        nat.kind = "natural";
        nat.title = "Full CI natural orbitals";
        const int nno = std::min<int>(last, static_cast<int>(U.cols()));
        for (int k = 0; k < nno; ++k) {
            OrbitalInfo o;
            o.number = k + 1;
            o.occupation = k < static_cast<int>(fci->natural_occupations.size()) ? fci->natural_occupations[k] : 0.0;
            o.label = "";
            nat.orbitals.push_back(o);
            std::vector<double> c(nbf, 0.0);
            for (std::size_t mu = 0; mu < nbf; ++mu) {
                double s = 0.0;
                for (std::size_t p = 0; p < nmo; ++p) s += rhf.C(mu, p) * U(p, static_cast<std::size_t>(k));
                c[mu] = s;
            }
            // sign convention: the largest coefficient positive
            std::size_t big = 0;
            for (std::size_t mu = 1; mu < nbf; ++mu)
                if (std::fabs(c[mu]) > std::fabs(c[big])) big = mu;
            if (c[big] < 0.0)
                for (double& x : c) x = -x;
            nat.coefficients.push_back(std::move(c));
        }
        d.sets.push_back(std::move(nat));
    }
    return d;
}

OrbitalData make_sto_orbital_data(int n_gaussians, double R, double zeta1, double zeta2, double ZA, double ZB,
                                  const std::vector<double>& C, const std::vector<double>& energies) {
    if (n_gaussians < 1 || n_gaussians > 3) throw std::invalid_argument("STO-nG: n must be 1, 2 or 3");
    if (C.size() != 4) throw std::invalid_argument("STO-nG orbitals: C must be 2 x 2");
    // Szabo & Ostlund, Appendix B: STO-nG expansions of a 1s Slater function with zeta = 1.
    static const double coef[3][3] = {{1.000000, 0.678914, 0.444635},
                                      {0.000000, 0.430129, 0.535328},
                                      {0.000000, 0.000000, 0.154329}};
    static const double expon[3][3] = {{0.270950, 0.151623, 0.109818},
                                       {0.000000, 0.851819, 0.405771},
                                       {0.000000, 0.000000, 2.227660}};
    OrbitalData d;
    d.basis_name = "STO-" + std::to_string(n_gaussians) + "G";
    d.n_basis = 2;
    auto atom = [](double Z, double z) {
        Atom a;
        a.Z = static_cast<int>(std::lround(Z));
        a.symbol = a.Z == 2 ? "He" : (a.Z == 1 ? "H" : "X");
        a.position = {0.0, 0.0, z};
        return a;
    };
    d.atoms = {atom(ZA, 0.0), atom(ZB, R)};
    const int col = n_gaussians - 1;
    for (int center = 0; center < 2; ++center) {
        Shell s;
        s.l = 0;
        s.atom = center;
        s.center = {0.0, 0.0, center == 0 ? 0.0 : R};
        s.first_function = center;
        const double zeta = center == 0 ? zeta1 : zeta2;
        for (int i = 0; i < n_gaussians; ++i) {
            const double a = expon[i][col] * zeta * zeta;
            s.exponents.push_back(a);
            // as in the program: the coefficient times the normalisation of exp(-a r^2)
            s.coefficients.push_back(coef[i][col] * std::pow(2.0 * a / kPiOrbitals, 0.75));
        }
        d.shells.push_back(s);
    }
    OrbitalSet set;
    set.kind = "canonical";
    set.title = "Hartree–Fock orbitals (Szabo–Ostlund)";
    for (int k = 0; k < 2; ++k) {
        OrbitalInfo o;
        o.number = k + 1;
        o.has_energy = k < static_cast<int>(energies.size());
        o.energy = o.has_energy ? energies[k] : 0.0;
        o.occupation = k == 0 ? 2.0 : 0.0;
        o.label = k == 0 ? "HOMO" : "LUMO";
        set.orbitals.push_back(o);
        set.coefficients.push_back({C[0 * 2 + k], C[1 * 2 + k]});  // row-major C(mu, k)
    }
    d.sets.push_back(std::move(set));
    return d;
}

void basis_function_values(const std::vector<Shell>& shells, int n_basis, const std::array<double, 3>& r,
                           double* values) {
    std::fill(values, values + n_basis, 0.0);
    double mono[(kMaxAngularMomentum + 1) * (kMaxAngularMomentum + 2) / 2];
    for (const Shell& s : shells) {
        const double x = r[0] - s.center[0], y = r[1] - s.center[1], z = r[2] - s.center[2];
        const double r2 = x * x + y * y + z * z;
        double radial = 0.0;
        for (std::size_t p = 0; p < s.exponents.size(); ++p) {
            const double e = s.exponents[p] * r2;
            if (e < 60.0) radial += s.coefficients[p] * std::exp(-e);
        }
        if (radial == 0.0) continue;
        const int l = s.l;
        double* out = values + s.first_function;
        if (l == 0) {
            out[0] = radial * solid_harmonic_matrix(0)[0];
            continue;
        }
        cartesian_monomials(l, x, y, z, mono);
        // the functions of the shell: real solid harmonics, or (14.0) normalised Cartesian components
        const ShellTransform& T = shell_transform(s);
        for (int m = 0; m < T.size(); ++m) {
            double v = 0.0;
            for (const auto& [c, f] : T.rows[static_cast<std::size_t>(m)]) v += f * mono[c];
            out[m] = v * radial;
        }
    }
}

double orbital_value(const OrbitalData& d, const std::vector<double>& c, const std::array<double, 3>& r) {
    std::vector<double> phi(static_cast<std::size_t>(d.n_basis));
    basis_function_values(d.shells, d.n_basis, r, phi.data());
    double v = 0.0;
    for (std::size_t mu = 0; mu < phi.size() && mu < c.size(); ++mu) v += c[mu] * phi[mu];
    return v;
}

std::vector<double> orbital_grid(const OrbitalData& d, const std::vector<double>& c, const GridSpec& g, int threads) {
    if (c.size() != static_cast<std::size_t>(d.n_basis))
        throw std::invalid_argument("orbital_grid: the coefficients do not match the basis");
    std::vector<double> out(g.size(), 0.0);
    // Shells whose coefficients are all zero in this orbital contribute nothing.
    std::vector<Shell> used;
    std::vector<int> first;
    int n_used = 0;
    for (const Shell& s : d.shells) {
        bool any = false;
        for (int m = 0; m < s.size(); ++m) any = any || c[static_cast<std::size_t>(s.first_function + m)] != 0.0;
        if (!any) continue;
        Shell t = s;
        first.push_back(s.first_function);
        t.first_function = n_used;
        n_used += s.size();
        used.push_back(t);
    }
    std::vector<double> cu(static_cast<std::size_t>(n_used));
    for (std::size_t k = 0; k < used.size(); ++k)
        for (int m = 0; m < used[k].size(); ++m)
            cu[static_cast<std::size_t>(used[k].first_function + m)] = c[static_cast<std::size_t>(first[k] + m)];
    const int nth = std::max(1, std::min(resolve_threads(threads), g.n[2] * g.n[1]));
    std::atomic<int> next{0};
    const int rows = g.n[2] * g.n[1];
    run_parallel(nth, [&](int) {
        std::vector<double> phi(static_cast<std::size_t>(std::max(n_used, 1)));
        for (int row; (row = next++) < rows;) {
            const int k = row / g.n[1], j = row % g.n[1];
            for (int i = 0; i < g.n[0]; ++i) {
                basis_function_values(used, n_used, g.point(i, j, k), phi.data());
                double v = 0.0;
                for (int mu = 0; mu < n_used; ++mu) v += cu[static_cast<std::size_t>(mu)] * phi[static_cast<std::size_t>(mu)];
                out[g.index(i, j, k)] = v;
            }
        }
    });
    return out;
}

Mesh isosurface(const std::vector<double>& f, const GridSpec& g, double level) {
    Mesh mesh;
    const int nx = g.n[0], ny = g.n[1], nz = g.n[2];
    if (nx < 2 || ny < 2 || nz < 2 || f.size() != g.size()) return mesh;
    const double hx = g.axes[0][0], hy = g.axes[1][1], hz = g.axes[2][2];
    // -grad f at a grid node (central differences, one-sided at the faces)
    auto minus_gradient = [&](int i, int j, int k) {
        auto value = [&](int a, int b, int c) { return f[g.index(a, b, c)]; };
        std::array<double, 3> gr{};
        const int ia = std::max(i - 1, 0), ib = std::min(i + 1, nx - 1);
        const int ja = std::max(j - 1, 0), jb = std::min(j + 1, ny - 1);
        const int ka = std::max(k - 1, 0), kb = std::min(k + 1, nz - 1);
        gr[0] = -(value(ib, j, k) - value(ia, j, k)) / ((ib - ia) * hx);
        gr[1] = -(value(i, jb, k) - value(i, ja, k)) / ((jb - ja) * hy);
        gr[2] = -(value(i, j, kb) - value(i, j, ka)) / ((kb - ka) * hz);
        return gr;
    };
    std::unordered_map<std::uint64_t, int> edge_vertex;
    const std::uint64_t total = static_cast<std::uint64_t>(g.size());
    auto vertex_on_edge = [&](std::size_t a, std::size_t b, int ia, int ja, int ka, int ib, int jb, int kb) {
        const std::uint64_t key = std::min<std::uint64_t>(a, b) * total + std::max<std::uint64_t>(a, b);
        auto it = edge_vertex.find(key);
        if (it != edge_vertex.end()) return it->second;
        const double fa = f[a], fb = f[b];
        double t = (fa == fb) ? 0.5 : (level - fa) / (fb - fa);
        t = std::min(1.0, std::max(0.0, t));
        const std::array<double, 3> pa = g.point(ia, ja, ka), pb = g.point(ib, jb, kb);
        const std::array<double, 3> na = minus_gradient(ia, ja, ka), nb = minus_gradient(ib, jb, kb);
        std::array<double, 3> n{};
        double len = 0.0;
        for (int d = 0; d < 3; ++d) {
            mesh.positions.push_back(static_cast<float>(pa[d] + t * (pb[d] - pa[d])));
            n[d] = na[d] + t * (nb[d] - na[d]);
            len += n[d] * n[d];
        }
        len = std::sqrt(len);
        for (int d = 0; d < 3; ++d) mesh.normals.push_back(static_cast<float>(len > 0.0 ? n[d] / len : 0.0));
        const int id = static_cast<int>(mesh.positions.size() / 3 - 1);
        edge_vertex.emplace(key, id);
        return id;
    };
    // Six tetrahedra around the diagonal 0-7 of each cube (corner = dx + 2 dy + 4 dz).
    static const int tets[6][4] = {{0, 1, 3, 7}, {0, 3, 2, 7}, {0, 2, 6, 7}, {0, 6, 4, 7}, {0, 4, 5, 7}, {0, 5, 1, 7}};
    auto add_triangle = [&](int a, int b, int c) {
        if (a == b || b == c || a == c) return;
        const float* pa = &mesh.positions[3 * static_cast<std::size_t>(a)];
        const float* pb = &mesh.positions[3 * static_cast<std::size_t>(b)];
        const float* pc = &mesh.positions[3 * static_cast<std::size_t>(c)];
        const double u[3] = {pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]};
        const double v[3] = {pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]};
        const double cr[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
        if (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2] < 1e-24) return;  // zero area
        double outward = 0.0;
        for (int d = 0; d < 3; ++d)
            outward += cr[d] * (mesh.normals[3 * static_cast<std::size_t>(a) + d] + mesh.normals[3 * static_cast<std::size_t>(b) + d] +
                                mesh.normals[3 * static_cast<std::size_t>(c) + d]);
        mesh.indices.push_back(a);
        if (outward >= 0.0) {
            mesh.indices.push_back(b);
            mesh.indices.push_back(c);
        } else {
            mesh.indices.push_back(c);
            mesh.indices.push_back(b);
        }
    };
    for (int k = 0; k + 1 < nz; ++k)
        for (int j = 0; j + 1 < ny; ++j)
            for (int i = 0; i + 1 < nx; ++i) {
                int ci[8], cj[8], ck[8];
                std::size_t id[8];
                double v[8];
                bool any_in = false, any_out = false;
                for (int corner = 0; corner < 8; ++corner) {
                    ci[corner] = i + (corner & 1);
                    cj[corner] = j + ((corner >> 1) & 1);
                    ck[corner] = k + ((corner >> 2) & 1);
                    id[corner] = g.index(ci[corner], cj[corner], ck[corner]);
                    v[corner] = f[id[corner]];
                    (v[corner] >= level ? any_in : any_out) = true;
                }
                if (!any_in || !any_out) continue;
                for (const auto& t : tets) {
                    int in[4], out[4], nin = 0, nout = 0;
                    for (int q = 0; q < 4; ++q) (v[t[q]] >= level ? in[nin++] : out[nout++]) = t[q];
                    if (nin == 0 || nout == 0) continue;
                    auto vtx = [&](int a, int b) {
                        return vertex_on_edge(id[a], id[b], ci[a], cj[a], ck[a], ci[b], cj[b], ck[b]);
                    };
                    if (nin == 1) {
                        add_triangle(vtx(in[0], out[0]), vtx(in[0], out[1]), vtx(in[0], out[2]));
                    } else if (nin == 3) {
                        add_triangle(vtx(out[0], in[0]), vtx(out[0], in[1]), vtx(out[0], in[2]));
                    } else {
                        // a quadrilateral: in0-out0, in0-out1, in1-out1, in1-out0
                        const int a = vtx(in[0], out[0]), b = vtx(in[0], out[1]), c = vtx(in[1], out[1]),
                                  e = vtx(in[1], out[0]);
                        add_triangle(a, b, c);
                        add_triangle(a, c, e);
                    }
                }
            }
    return mesh;
}

namespace {

// Moves every vertex onto sign * psi = isovalue of the exact orbital (up to three Newton steps
// along its normal, each at most half a grid spacing) and recomputes its normal from the
// gradient of the orbital (central differences): the grid's linear interpolation is poor where
// the orbital changes fast, e.g. near a nucleus.
void refine_surface(Mesh& m, const OrbitalData& d, const std::vector<double>& c, double sign, double isovalue, double h,
                    int threads) {
    const std::size_t nv = m.vertex_count();
    if (nv == 0) return;
    const int nth = std::max(1, std::min<int>(resolve_threads(threads), static_cast<int>(nv / 64) + 1));
    std::atomic<std::size_t> next{0};
    const std::size_t chunk = 256;
    run_parallel(nth, [&](int) {
        std::vector<double> phi(static_cast<std::size_t>(std::max(d.n_basis, 1)));
        auto f = [&](const std::array<double, 3>& p) {
            basis_function_values(d.shells, d.n_basis, p, phi.data());
            double v = 0.0;
            for (std::size_t mu = 0; mu < phi.size(); ++mu) v += c[mu] * phi[mu];
            return sign * v;
        };
        for (std::size_t start; (start = next.fetch_add(chunk)) < nv;) {
            for (std::size_t vi = start; vi < std::min(nv, start + chunk); ++vi) {
                std::array<double, 3> p{}, n{};
                for (int k = 0; k < 3; ++k) {
                    p[k] = m.positions[3 * vi + k];
                    n[k] = m.normals[3 * vi + k];
                }
                const double delta = 0.02 * h;
                for (int iter = 0; iter < 3; ++iter) {
                    std::array<double, 3> a = p, b = p;
                    for (int k = 0; k < 3; ++k) {
                        a[k] += delta * n[k];
                        b[k] -= delta * n[k];
                    }
                    const double f0 = f(p) - isovalue;
                    const double slope = (f(a) - f(b)) / (2.0 * delta);
                    if (!(std::fabs(slope) > 1e-12)) break;
                    const double t = std::max(-0.5 * h, std::min(0.5 * h, -f0 / slope));
                    for (int k = 0; k < 3; ++k) p[k] += t * n[k];
                    if (std::fabs(t) < 1e-7) break;
                }
                std::array<double, 3> g{};
                double len = 0.0;
                for (int k = 0; k < 3; ++k) {
                    std::array<double, 3> a = p, b = p;
                    a[k] += delta;
                    b[k] -= delta;
                    g[k] = -(f(a) - f(b)) / (2.0 * delta);  // outward: sign * psi decreases
                    len += g[k] * g[k];
                }
                len = std::sqrt(len);
                for (int k = 0; k < 3; ++k) {
                    m.positions[3 * vi + k] = static_cast<float>(p[k]);
                    if (len > 0.0) m.normals[3 * vi + k] = static_cast<float>(g[k] / len);
                }
            }
        }
    });
}

}  // namespace

OrbitalSurfaces orbital_surfaces(const OrbitalData& d, const std::vector<double>& c, double isovalue, int max_points,
                                 int threads) {
    if (!(isovalue > 0.0) || !std::isfinite(isovalue)) throw std::invalid_argument("the isovalue must be positive");
    if (d.atoms.empty()) throw std::invalid_argument("orbital_surfaces: no atoms");
    OrbitalSurfaces s;
    std::array<double, 3> lo{}, hi{};
    for (int k = 0; k < 3; ++k) {
        lo[k] = std::numeric_limits<double>::infinity();
        hi[k] = -std::numeric_limits<double>::infinity();
    }
    for (const Atom& a : d.atoms)
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], a.position[k]);
            hi[k] = std::max(hi[k], a.position[k]);
        }
    // coarse search for the region |psi| >= isovalue
    const double coarse = 0.5, reach = 12.0;
    GridSpec cg;
    for (int k = 0; k < 3; ++k) {
        cg.origin[k] = lo[k] - reach;
        cg.n[k] = static_cast<int>(std::ceil((hi[k] - lo[k] + 2.0 * reach) / coarse)) + 1;
        cg.axes[k] = {0.0, 0.0, 0.0};
        cg.axes[k][k] = coarse;
    }
    const std::vector<double> cv = orbital_grid(d, c, cg, threads);
    std::array<double, 3> blo{}, bhi{};
    for (int k = 0; k < 3; ++k) {
        blo[k] = std::numeric_limits<double>::infinity();
        bhi[k] = -std::numeric_limits<double>::infinity();
    }
    bool found = false;
    for (int kk = 0; kk < cg.n[2]; ++kk)
        for (int jj = 0; jj < cg.n[1]; ++jj)
            for (int ii = 0; ii < cg.n[0]; ++ii)
                if (std::fabs(cv[cg.index(ii, jj, kk)]) >= 0.5 * isovalue) {
                    found = true;
                    const std::array<double, 3> p = cg.point(ii, jj, kk);
                    for (int k = 0; k < 3; ++k) {
                        blo[k] = std::min(blo[k], p[k]);
                        bhi[k] = std::max(bhi[k], p[k]);
                    }
                }
    if (!found) {
        // nothing reaches half the isovalue on the coarse grid: a box around the atoms
        for (int k = 0; k < 3; ++k) {
            blo[k] = lo[k] - 2.0;
            bhi[k] = hi[k] + 2.0;
        }
    }
    const double margin = 1.0;
    double volume = 1.0;
    for (int k = 0; k < 3; ++k) {
        blo[k] = std::max(blo[k] - margin, cg.origin[k]);
        bhi[k] = std::min(bhi[k] + margin, cg.origin[k] + (cg.n[k] - 1) * coarse);
        volume *= bhi[k] - blo[k];
    }
    const double target = std::max(1000, max_points);
    double h = std::max(0.04, std::cbrt(volume / target));
    GridSpec g;
    for (;;) {
        for (int k = 0; k < 3; ++k) {
            g.origin[k] = blo[k];
            g.n[k] = std::max(2, static_cast<int>(std::ceil((bhi[k] - blo[k]) / h)) + 1);
            g.axes[k] = {0.0, 0.0, 0.0};
            g.axes[k][k] = h;
        }
        if (static_cast<double>(g.size()) <= target || h > 1.0) break;
        h *= 1.02;
    }
    std::vector<double> v = orbital_grid(d, c, g, threads);
    for (double x : v) s.max_abs = std::max(s.max_abs, std::fabs(x));
    s.positive = isosurface(v, g, isovalue);
    for (double& x : v) x = -x;
    s.negative = isosurface(v, g, isovalue);
    s.grid = g;
    refine_surface(s.positive, d, c, 1.0, isovalue, h, threads);
    refine_surface(s.negative, d, c, -1.0, isovalue, h, threads);
    return s;
}

}  // namespace macqc
