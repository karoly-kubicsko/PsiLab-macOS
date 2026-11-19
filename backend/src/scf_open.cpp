// ΨLab 15.0: open-shell Hartree-Fock. Unrestricted (UHF, Pople-Nesbet) and restricted open-shell
// (ROHF, Roothaan effective Fock matrix) determinants for any spin multiplicity, with the same
// integrals, orthogonalisation, guesses and DIIS as the closed-shell RHF of rhf.cpp.
//
//   UHF:  F_a = H + J(P_a + P_b) - K(P_a),  F_b = H + J(P_a + P_b) - K(P_b),
//         E = 1/2 tr[P_a (H + F_a)] + 1/2 tr[P_b (H + F_b)] + V_nn.
//   ROHF: the same energy with P_a = C_(docc+socc) C', P_b = C_docc C' of one set of orbitals,
//         made stationary by the Roothaan effective Fock matrix (MO basis, blocks c = doubly,
//         o = singly occupied, v = virtual)
//             F_eff = (F_a + F_b)/2 except F_eff(c,o) = F_b(c,o), F_eff(o,v) = F_a(o,v),
//         which commutes with the occupations n = 2, 1, 0 exactly at the ROHF solution: the DIIS
//         error vector is X' (F_eff P S - S P F_eff) X with P = P_a + P_b. The orbital energies
//         are the eigenvalues of F_eff (PySCF's convention, rohf.get_roothaan_fock).
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "macqc/gto_integrals.h"
#include "macqc/job.h"
#include "macqc/rhf.h"

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

Matrix mm(const Matrix& a, const Matrix& b) {
    const std::size_t n = a.rows(), m = b.cols(), kk = a.cols();
    Matrix c(n, m);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < kk; ++k) {
            const double aik = a(i, k);
            if (aik == 0.0) continue;
            for (std::size_t j = 0; j < m; ++j) c(i, j) += aik * b(k, j);
        }
    return c;
}

Matrix mtm(const Matrix& a, const Matrix& b) {  // a' b
    const std::size_t n = a.cols(), m = b.cols(), kk = a.rows();
    Matrix c(n, m);
    for (std::size_t k = 0; k < kk; ++k)
        for (std::size_t i = 0; i < n; ++i) {
            const double aki = a(k, i);
            if (aki == 0.0) continue;
            for (std::size_t j = 0; j < m; ++j) c(i, j) += aki * b(k, j);
        }
    return c;
}

double dot(const Matrix& a, const Matrix& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.rows(); ++i)
        for (std::size_t j = 0; j < a.cols(); ++j) s += a(i, j) * b(i, j);
    return s;
}

Matrix add(const Matrix& a, const Matrix& b, double fb = 1.0) {
    Matrix c(a.rows(), a.cols());
    for (std::size_t i = 0; i < a.rows(); ++i)
        for (std::size_t j = 0; j < a.cols(); ++j) c(i, j) = a(i, j) + fb * b(i, j);
    return c;
}

Matrix scaled(const Matrix& a, double f) {
    Matrix c = a;
    for (std::size_t i = 0; i < a.rows(); ++i)
        for (std::size_t j = 0; j < a.cols(); ++j) c(i, j) *= f;
    return c;
}

bool solve_linear(std::vector<double> a, std::vector<double> b, std::size_t n, std::vector<double>& x) {
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t piv = col;
        for (std::size_t r = col + 1; r < n; ++r)
            if (std::fabs(a[r * n + col]) > std::fabs(a[piv * n + col])) piv = r;
        if (std::fabs(a[piv * n + col]) < 1.0e-14) return false;
        if (piv != col) {
            for (std::size_t c = 0; c < n; ++c) std::swap(a[piv * n + c], a[col * n + c]);
            std::swap(b[piv], b[col]);
        }
        for (std::size_t r = col + 1; r < n; ++r) {
            const double f = a[r * n + col] / a[col * n + col];
            if (f == 0.0) continue;
            for (std::size_t c = col; c < n; ++c) a[r * n + c] -= f * a[col * n + c];
            b[r] -= f * b[col];
        }
    }
    x.assign(n, 0.0);
    for (std::size_t r = n; r-- > 0;) {
        double s = b[r];
        for (std::size_t c = r + 1; c < n; ++c) s -= a[r * n + c] * x[c];
        x[r] = s / a[r * n + r];
    }
    return true;
}

// P = C_occ C_occ' (one spin).
Matrix spin_density(const Matrix& C, int nocc) {
    const std::size_t n = C.rows();
    Matrix P(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j <= i; ++j) {
            double s = 0.0;
            for (int k = 0; k < nocc; ++k) s += C(i, static_cast<std::size_t>(k)) * C(j, static_cast<std::size_t>(k));
            P(i, j) = P(j, i) = s;
        }
    return P;
}

// The orthonormalised first `nocc` columns of guess orbitals (Loewdin in the metric S), or an
// empty matrix if they do not fit.
Matrix orthonormal_occupied(const Matrix* guess, int nocc, const Matrix& S) {
    const std::size_t n = S.rows();
    if (!guess || guess->rows() != n || guess->cols() < static_cast<std::size_t>(nocc)) return Matrix();
    if (nocc == 0) return Matrix(n, 0);
    Matrix Co(n, static_cast<std::size_t>(nocc));
    for (std::size_t i = 0; i < n; ++i)
        for (int k = 0; k < nocc; ++k) Co(i, static_cast<std::size_t>(k)) = (*guess)(i, static_cast<std::size_t>(k));
    const Matrix M = mtm(Co, mm(S, Co));
    const EigenResult me = symmetric_eigen(M);
    if (me.values.empty() || me.values.front() < 1.0e-6) return Matrix();
    const std::size_t no = static_cast<std::size_t>(nocc);
    Matrix Mh(no, no);
    for (std::size_t a = 0; a < no; ++a)
        for (std::size_t b = 0; b < no; ++b) {
            double s = 0.0;
            for (std::size_t k = 0; k < no; ++k)
                s += me.vectors(a, k) * me.vectors(b, k) / std::sqrt(me.values[k]);
            Mh(a, b) = s;
        }
    return mm(Co, Mh);
}

// Common set-up: integrals, orthogonalisation, Fock builder.
struct OpenSetup {
    std::size_t n = 0, nmo = 0;
    Matrix S, T, V, H, X;
};

OpenSetup setup_integrals(const Molecule& mol, const BasisSet& basis, const RHFOptions& opt, RHFResult& r, int nmin) {
    OpenSetup s;
    s.n = static_cast<std::size_t>(basis.n_functions);
    if (s.n < static_cast<std::size_t>(nmin)) throw std::runtime_error("the basis set has fewer functions than occupied orbitals");
    s.S = overlap_matrix(basis);
    s.T = kinetic_matrix(basis);
    s.V = nuclear_attraction_matrix(basis, mol);
    s.H = add(s.T, s.V);
    if (opt.field[0] != 0.0 || opt.field[1] != 0.0 || opt.field[2] != 0.0) {
        const std::vector<Matrix> D = dipole_matrices(basis, {0.0, 0.0, 0.0});
        for (int k = 0; k < 3; ++k) {
            const double f = opt.field[static_cast<std::size_t>(k)];
            if (f == 0.0) continue;
            s.H = add(s.H, D[static_cast<std::size_t>(k)], f);
            for (const Atom& a : mol.atoms) r.nuclear_repulsion -= f * a.nuclear_charge() * a.position[static_cast<std::size_t>(k)];
        }
    }
    r.S = s.S;
    r.H = s.H;
    const EigenResult se = symmetric_eigen(s.S);
    r.smallest_overlap_eigenvalue = se.values.empty() ? 0.0 : se.values.front();
    std::vector<std::size_t> keep;
    for (std::size_t k = 0; k < s.n; ++k)
        if (se.values[k] >= opt.lindep_tol) keep.push_back(k);
    s.nmo = keep.size();
    r.removed_functions = static_cast<int>(s.n - s.nmo);
    r.n_mo = static_cast<int>(s.nmo);
    if (s.nmo < static_cast<std::size_t>(nmin)) throw std::runtime_error("too many linearly dependent basis functions");
    s.X = Matrix(s.n, s.nmo);
    for (std::size_t k = 0; k < s.nmo; ++k) {
        const double f = 1.0 / std::sqrt(se.values[keep[k]]);
        for (std::size_t i = 0; i < s.n; ++i) s.X(i, k) = se.vectors(i, keep[k]) * f;
    }
    return s;
}

// DIIS: extrapolation coefficients of the stored error vectors (sets of matrices), or empty.
std::vector<double> diis_coefficients(std::vector<std::vector<Matrix>>& errs) {
    while (errs.size() >= 2) {
        const std::size_t m = errs.size(), dim = m + 1;
        std::vector<double> B(dim * dim, 0.0), rhs(dim, 0.0), c;
        double scale = 0.0;
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j <= i; ++j) {
                double v = 0.0;
                for (std::size_t k = 0; k < errs[i].size(); ++k) v += dot(errs[i][k], errs[j][k]);
                B[i * dim + j] = B[j * dim + i] = v;
                if (i == j) scale = std::max(scale, v);
            }
        if (scale > 0.0)
            for (std::size_t i = 0; i < m; ++i)
                for (std::size_t j = 0; j < m; ++j) B[i * dim + j] /= scale;
        for (std::size_t i = 0; i < m; ++i) B[i * dim + m] = B[m * dim + i] = -1.0;
        rhs[m] = -1.0;
        if (solve_linear(B, rhs, dim, c)) {
            c.resize(m);
            return c;
        }
        errs.erase(errs.begin());
        return {};  // caller drops its oldest Fock matrices as well
    }
    return {};
}

double max_abs(const Matrix& m) {
    double g = 0.0;
    for (std::size_t i = 0; i < m.rows(); ++i)
        for (std::size_t j = 0; j < m.cols(); ++j) g = std::max(g, std::fabs(m(i, j)));
    return g;
}

// <S^2> of a determinant with alpha orbitals Ca and beta orbitals Cb (first na, nb occupied).
double determinant_s2(const Matrix& Ca, const Matrix& Cb, int na, int nb, const Matrix& S) {
    const double sz = 0.5 * (na - nb);
    double ov = 0.0;
    if (na > 0 && nb > 0) {
        const Matrix SCb = mm(S, Cb);
        for (int i = 0; i < na; ++i)
            for (int j = 0; j < nb; ++j) {
                double s = 0.0;
                for (std::size_t m = 0; m < S.rows(); ++m) s += Ca(m, static_cast<std::size_t>(i)) * SCb(m, static_cast<std::size_t>(j));
                ov += s * s;
            }
    }
    return sz * (sz + 1.0) + nb - ov;
}

// The energy components of the result from the converged densities and Fock matrices.
void finish(RHFResult& r, const OpenSetup& s, const Matrix& Pa, const Matrix& Pb, const Matrix& Fa, const Matrix& Fb) {
    r.Pa = Pa;
    r.Pb = Pb;
    r.Fa = Fa;
    r.Fb = Fb;
    r.P = add(Pa, Pb);
    r.kinetic_energy = dot(r.P, s.T);
    r.nuclear_attraction_energy = dot(r.P, s.V);
    r.two_electron_energy = 0.5 * (dot(Pa, add(Fa, s.H, -1.0)) + dot(Pb, add(Fb, s.H, -1.0)));
    const double potential = r.total_energy - r.kinetic_energy;
    r.virial_ratio = r.kinetic_energy != 0.0 ? -potential / r.kinetic_energy : 0.0;
}

}  // namespace

std::pair<int, int> alpha_beta_electrons(int electrons, int multiplicity) {
    if (multiplicity < 1) throw std::runtime_error("the spin multiplicity must be at least 1");
    const int twice_s = multiplicity - 1;
    if (electrons <= 0) throw std::runtime_error("the molecule has no electrons");
    if (twice_s > electrons)
        throw std::runtime_error(std::to_string(electrons) + " electron" + (electrons == 1 ? "" : "s") +
                                 " cannot have spin multiplicity " + std::to_string(multiplicity) +
                                 " (at most " + std::to_string(electrons + 1) + ")");
    if ((electrons - twice_s) % 2 != 0)
        throw std::runtime_error(std::to_string(electrons) + " electrons cannot have spin multiplicity " +
                                 std::to_string(multiplicity) + ": an " +
                                 (electrons % 2 ? "odd number of electrons needs an even multiplicity (doublet, quartet, ...)"
                                                : "even number of electrons needs an odd multiplicity (singlet, triplet, ...)"));
    return {(electrons + twice_s) / 2, (electrons - twice_s) / 2};
}

std::string canonical_reference_name(const std::string& name) {
    std::string t;
    for (char c : name)
        if (c != ' ' && c != '-' && c != '_') t += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (t == "RHF" || t == "R") return "RHF";
    if (t == "UHF" || t == "U") return "UHF";
    if (t == "ROHF" || t == "RO") return "ROHF";
    return "";
}

RHFResult run_uhf(const Molecule& mol, const BasisSet& basis, const RHFOptions& opt, const RHFProgress& progress) {
    check_ecps(mol, basis);
    const auto ab = alpha_beta_electrons(mol.electrons(), opt.multiplicity);
    const int na = ab.first, nb = ab.second;
    report_progress(JobStage::Integrals);
    check_cancelled();
    RHFResult r;
    r.reference = "UHF";
    r.multiplicity = opt.multiplicity;
    r.n_alpha = na;
    r.n_beta = nb;
    r.n_occupied = na;
    r.n_basis = basis.n_functions;
    r.nuclear_repulsion = mol.nuclear_repulsion();
    const auto t0 = Clock::now();
    const OpenSetup s = setup_integrals(mol, basis, opt, r, na);
    const std::size_t n = s.n, nmo = s.nmo;
    FockBuilder fock(basis, opt.threads, opt.memory_limit_mb);
    r.eri_stored = fock.stored();
    r.unique_integrals = fock.unique_integrals();
    r.threads = fock.threads();
    r.integral_seconds = seconds_since(t0);
    const auto t1 = Clock::now();

    auto diagonalize = [&](const Matrix& F, Matrix& C, std::vector<double>& eps) {
        const EigenResult e = symmetric_eigen(mtm(s.X, mm(F, s.X)));
        C = mm(s.X, e.vectors);
        eps = e.values;
    };

    Matrix Ca, Cb;
    std::vector<double> ea, eb;
    diagonalize(s.H, Ca, ea);
    Cb = Ca;
    eb = ea;
    Matrix Pa = spin_density(Ca, na), Pb = spin_density(Cb, nb);
    r.guess = "core";
    const Matrix Goa = orthonormal_occupied(opt.guess_orbitals, na, s.S);
    const Matrix Gob = orthonormal_occupied(opt.guess_orbitals_b ? opt.guess_orbitals_b : opt.guess_orbitals, nb, s.S);
    if (Goa.rows() && Gob.rows()) {
        Pa = spin_density(Goa, na);
        Pb = spin_density(Gob, nb);
        r.guess_from_orbitals = true;
        r.guess = "orbitals";
    } else if (opt.guess == "sad") {
        const Matrix P = sad_guess_density(mol, basis, opt.threads, opt.memory_limit_mb);
        Pa = scaled(P, 0.5);
        Pb = Pa;
        r.guess = "sad";
    } else if (opt.guess != "core") {
        throw std::invalid_argument("unknown SCF guess '" + opt.guess + "' (sad or core)");
    }

    // Level shift of the virtual orbitals while the orbital gradient is large (open shells with
    // degenerate partly filled shells, e.g. the B, C, O, F atoms, would otherwise oscillate).
    const double shift = opt.level_shift > 0.0 ? opt.level_shift : 0.3;
    std::vector<Matrix> hist_fa, hist_fb;
    std::vector<std::vector<Matrix>> hist_e;
    double e_prev = 0.0;
    Matrix Fa, Fb;
    bool shifting = na != nb;
    for (int it = 1; it <= opt.max_iter; ++it) {
        check_cancelled();
        std::vector<Matrix> J, K;
        fock.build_jk({Pa, Pb}, J, K);
        const Matrix Jt = add(J[0], J[1]);
        Fa = add(add(s.H, Jt), K[0], -1.0);
        Fb = add(add(s.H, Jt), K[1], -1.0);
        const double e_el = 0.5 * (dot(Pa, add(s.H, Fa)) + dot(Pb, add(s.H, Fb)));
        const double e_tot = e_el + r.nuclear_repulsion;
        const Matrix FPSa = mm(Fa, mm(Pa, s.S)), FPSb = mm(Fb, mm(Pb, s.S));
        const Matrix erra = mtm(s.X, mm(add(FPSa, FPSa.transposed(), -1.0), s.X));
        const Matrix errb = mtm(s.X, mm(add(FPSb, FPSb.transposed(), -1.0), s.X));
        const double gmax = std::max(max_abs(erra), max_abs(errb));

        RHFIteration rec;
        rec.iteration = it;
        rec.energy = e_tot;
        rec.delta_e = it == 1 ? 0.0 : e_tot - e_prev;
        rec.gradient = gmax;
        r.electronic_energy = e_el;
        r.total_energy = e_tot;
        r.iterations = it;
        report_progress(JobStage::SCF, it, e_tot, gmax);
        const bool done = it > 1 && std::fabs(e_tot - e_prev) < opt.energy_tol && gmax < opt.gradient_tol;
        e_prev = e_tot;
        if (done) {
            r.converged = true;
            r.history.push_back(rec);
            if (progress) progress(rec);
            break;
        }
        if (gmax < 1.0e-2) shifting = false;

        Matrix Fxa = Fa, Fxb = Fb;
        hist_fa.push_back(Fa);
        hist_fb.push_back(Fb);
        hist_e.push_back({erra, errb});
        if (static_cast<int>(hist_fa.size()) > opt.diis_size) {
            hist_fa.erase(hist_fa.begin());
            hist_fb.erase(hist_fb.begin());
            hist_e.erase(hist_e.begin());
        }
        if (!shifting) {
            while (hist_e.size() >= 2) {
                const std::size_t before = hist_e.size();
                const std::vector<double> c = diis_coefficients(hist_e);
                if (!c.empty()) {
                    Fxa = Matrix(n, n);
                    Fxb = Matrix(n, n);
                    for (std::size_t i = 0; i < c.size(); ++i) {
                        Fxa = add(Fxa, hist_fa[i], c[i]);
                        Fxb = add(Fxb, hist_fb[i], c[i]);
                    }
                    rec.diis = true;
                    break;
                }
                if (hist_e.size() < before) {
                    hist_fa.erase(hist_fa.begin());
                    hist_fb.erase(hist_fb.begin());
                }
            }
        } else {
            // F + mu (S - S P S): the virtual orbitals of each spin raised by mu
            const Matrix SPSa = mm(s.S, mm(Pa, s.S)), SPSb = mm(s.S, mm(Pb, s.S));
            Fxa = add(Fa, add(s.S, SPSa, -1.0), shift);
            Fxb = add(Fb, add(s.S, SPSb, -1.0), shift);
        }
        r.history.push_back(rec);
        if (progress) progress(rec);
        diagonalize(Fxa, Ca, ea);
        diagonalize(Fxb, Cb, eb);
        Pa = spin_density(Ca, na);
        Pb = spin_density(Cb, nb);
    }
    if (r.converged) {
        diagonalize(Fa, Ca, ea);
        diagonalize(Fb, Cb, eb);
    }
    r.C = Ca;
    r.Cb = Cb;
    r.orbital_energies = ea;
    r.orbital_energies_b = eb;
    finish(r, s, Pa, Pb, Fa, Fb);
    r.s2 = determinant_s2(Ca, Cb, na, nb, s.S);
    const double S_spin = 0.5 * (opt.multiplicity - 1);
    r.s2_exact = S_spin * (S_spin + 1.0);
    r.scf_seconds = seconds_since(t1);
    (void)nmo;
    return r;
}

RHFResult run_rohf(const Molecule& mol, const BasisSet& basis, const RHFOptions& opt, const RHFProgress& progress) {
    check_ecps(mol, basis);
    const auto ab = alpha_beta_electrons(mol.electrons(), opt.multiplicity);
    const int na = ab.first, nb = ab.second;
    report_progress(JobStage::Integrals);
    check_cancelled();
    RHFResult r;
    r.reference = "ROHF";
    r.multiplicity = opt.multiplicity;
    r.n_alpha = na;
    r.n_beta = nb;
    r.n_occupied = na;
    r.n_basis = basis.n_functions;
    r.nuclear_repulsion = mol.nuclear_repulsion();
    const auto t0 = Clock::now();
    const OpenSetup s = setup_integrals(mol, basis, opt, r, na);
    const std::size_t n = s.n, nmo = s.nmo;
    FockBuilder fock(basis, opt.threads, opt.memory_limit_mb);
    r.eri_stored = fock.stored();
    r.unique_integrals = fock.unique_integrals();
    r.threads = fock.threads();
    r.integral_seconds = seconds_since(t0);
    const auto t1 = Clock::now();

    // occupation numbers of the MOs: 2 (closed), 1 (open), 0
    auto occupation = [&](std::size_t p) { return p < static_cast<std::size_t>(nb) ? 2 : (p < static_cast<std::size_t>(na) ? 1 : 0); };
    auto diagonalize = [&](const Matrix& F, Matrix& C, std::vector<double>& eps) {
        const EigenResult e = symmetric_eigen(mtm(s.X, mm(F, s.X)));
        C = mm(s.X, e.vectors);
        eps = e.values;
    };
    // The Roothaan effective Fock matrix of orbitals C (MO basis) and in the AO basis, S C F C' S.
    auto effective = [&](const Matrix& C, const Matrix& Fa, const Matrix& Fb, double shift, Matrix& Fmo) {
        const Matrix Fam = mtm(C, mm(Fa, C)), Fbm = mtm(C, mm(Fb, C));
        Fmo = Matrix(nmo, nmo);
        for (std::size_t p = 0; p < nmo; ++p)
            for (std::size_t q = 0; q < nmo; ++q) {
                const int np = occupation(p), nq = occupation(q);
                double v = 0.5 * (Fam(p, q) + Fbm(p, q));
                if ((np == 2 && nq == 1) || (np == 1 && nq == 2)) v = Fbm(p, q);
                else if ((np == 1 && nq == 0) || (np == 0 && nq == 1)) v = Fam(p, q);
                Fmo(p, q) = v;
            }
        if (shift > 0.0)
            for (std::size_t p = 0; p < nmo; ++p) {
                const int np = occupation(p);
                if (np == 1) Fmo(p, p) += 0.5 * shift;
                else if (np == 0) Fmo(p, p) += shift;
            }
        const Matrix SC = mm(s.S, C);
        return mm(SC, mm(Fmo, SC.transposed()));
    };

    Matrix C;
    std::vector<double> eps;
    diagonalize(s.H, C, eps);
    Matrix Pa = spin_density(C, na), Pb = spin_density(C, nb);
    r.guess = "core";
    // starting orbitals: the given ones, orthonormalised as a whole set
    const Matrix Go = orthonormal_occupied(opt.guess_orbitals, na, s.S);
    if (Go.rows()) {
        Pa = spin_density(Go, na);
        Pb = spin_density(Go, nb);
        r.guess_from_orbitals = true;
        r.guess = "orbitals";
    } else if (opt.guess == "sad") {
        const Matrix P = sad_guess_density(mol, basis, opt.threads, opt.memory_limit_mb);
        Pa = scaled(P, 0.5);
        Pb = Pa;
        r.guess = "sad";
    } else if (opt.guess != "core") {
        throw std::invalid_argument("unknown SCF guess '" + opt.guess + "' (sad or core)");
    }
    // ROHF needs orbitals for the effective Fock matrix: from the guess density, the natural
    // orbitals of P_a + P_b (the occupied ones first).
    {
        const Matrix P = add(Pa, Pb);
        const Matrix M = mtm(s.X, mm(s.S, mm(P, mm(s.S, s.X))));
        EigenResult e = symmetric_eigen(M);
        Matrix V(nmo, nmo);
        for (std::size_t k = 0; k < nmo; ++k)
            for (std::size_t i = 0; i < nmo; ++i) V(i, k) = e.vectors(i, nmo - 1 - k);  // descending occupation
        C = mm(s.X, V);
        Pa = spin_density(C, na);
        Pb = spin_density(C, nb);
    }

    const double shift = opt.level_shift > 0.0 ? opt.level_shift : 0.3;
    std::vector<Matrix> hist_f;
    std::vector<std::vector<Matrix>> hist_e;
    double e_prev = 0.0;
    Matrix Fa, Fb, Feff_mo;
    bool shifting = na != nb;
    for (int it = 1; it <= opt.max_iter; ++it) {
        check_cancelled();
        std::vector<Matrix> J, K;
        fock.build_jk({Pa, Pb}, J, K);
        const Matrix Jt = add(J[0], J[1]);
        Fa = add(add(s.H, Jt), K[0], -1.0);
        Fb = add(add(s.H, Jt), K[1], -1.0);
        const double e_el = 0.5 * (dot(Pa, add(s.H, Fa)) + dot(Pb, add(s.H, Fb)));
        const double e_tot = e_el + r.nuclear_repulsion;
        const Matrix Feff = effective(C, Fa, Fb, 0.0, Feff_mo);
        const Matrix P = add(Pa, Pb);
        const Matrix FPS = mm(Feff, mm(P, s.S));
        const Matrix err = mtm(s.X, mm(add(FPS, FPS.transposed(), -1.0), s.X));
        const double gmax = max_abs(err);

        RHFIteration rec;
        rec.iteration = it;
        rec.energy = e_tot;
        rec.delta_e = it == 1 ? 0.0 : e_tot - e_prev;
        rec.gradient = gmax;
        r.electronic_energy = e_el;
        r.total_energy = e_tot;
        r.iterations = it;
        report_progress(JobStage::SCF, it, e_tot, gmax);
        const bool done = it > 1 && std::fabs(e_tot - e_prev) < opt.energy_tol && gmax < opt.gradient_tol;
        e_prev = e_tot;
        if (done) {
            r.converged = true;
            r.history.push_back(rec);
            if (progress) progress(rec);
            break;
        }
        if (gmax < 1.0e-2) shifting = false;
        Matrix Fx = Feff;
        hist_f.push_back(Feff);
        hist_e.push_back({err});
        if (static_cast<int>(hist_f.size()) > opt.diis_size) {
            hist_f.erase(hist_f.begin());
            hist_e.erase(hist_e.begin());
        }
        if (!shifting) {
            while (hist_e.size() >= 2) {
                const std::size_t before = hist_e.size();
                const std::vector<double> c = diis_coefficients(hist_e);
                if (!c.empty()) {
                    Fx = Matrix(n, n);
                    for (std::size_t i = 0; i < c.size(); ++i) Fx = add(Fx, hist_f[i], c[i]);
                    rec.diis = true;
                    break;
                }
                if (hist_e.size() < before) hist_f.erase(hist_f.begin());
            }
        } else {
            Matrix tmp;
            Fx = effective(C, Fa, Fb, shift, tmp);
        }
        r.history.push_back(rec);
        if (progress) progress(rec);
        diagonalize(Fx, C, eps);
        Pa = spin_density(C, na);
        Pb = spin_density(C, nb);
    }
    if (r.converged) {
        // canonical ROHF orbitals: eigenvectors of the effective Fock matrix of the solution
        const Matrix Feff = effective(C, Fa, Fb, 0.0, Feff_mo);
        diagonalize(Feff, C, eps);
    }
    r.C = C;
    r.orbital_energies = eps;
    finish(r, s, Pa, Pb, Fa, Fb);
    r.s2 = determinant_s2(C, C, na, nb, s.S);
    const double S_spin = 0.5 * (opt.multiplicity - 1);
    r.s2_exact = S_spin * (S_spin + 1.0);
    r.scf_seconds = seconds_since(t1);
    return r;
}

RHFResult run_scf(const Molecule& mol, const BasisSet& basis, const RHFOptions& opt, const RHFProgress& progress) {
    const std::string ref = canonical_reference_name(opt.reference.empty() ? "RHF" : opt.reference);
    if (ref == "UHF") return run_uhf(mol, basis, opt, progress);
    if (ref == "ROHF") return run_rohf(mol, basis, opt, progress);
    if (ref != "RHF") throw std::invalid_argument("unknown reference '" + opt.reference + "' (RHF, UHF or ROHF)");
    if (opt.multiplicity != 1)
        throw std::runtime_error("restricted Hartree-Fock (RHF) needs a closed-shell singlet; use UHF or ROHF for "
                                 "multiplicity " + std::to_string(opt.multiplicity));
    return run_rhf(mol, basis, opt, progress);
}

}  // namespace macqc
