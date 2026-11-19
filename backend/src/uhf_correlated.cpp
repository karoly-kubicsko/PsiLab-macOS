// ΨLab 15.0: analytic gradients of correlated methods on a UHF reference (open shells).
//
// The method supplies its Lagrangian densities in the spatial orbitals of each spin
// (so_methods.h): E = sum_s D_s h_s + sum G_aa (pq|rs)_aa + G_bb (pq|rs)_bb + G_ab (pq|rs)_ab + V_nn.
// With the orbital Lagrangians
//   X^a_tp = 2 sum_q h^a_tq D^a_qp + 4 sum_qrs (tq|rs)_aa G_aa,pqrs + 2 sum_qrs (tq|rs)_ab G_ab,pqrs,
//   X^b_tr = 2 sum_s h^b_ts D^b_sr + 4 sum (tq|rs)_bb G_bb + 2 sum_pqs (pq|ts)_ab G_ab,pqrs,
// (G_aa, G_bb 8-fold, G_ab 4-fold symmetrised), the occupied-virtual orbital response of both spins
// is handled by one Z-vector of the coupled-perturbed UHF equations,
//   (e_a - e_i) Z^s_ai + [J(DZ_a + DZ_b) - K(DZ_s)]_ai = X^s_ai - X^s_ia,  DZ_s = sym(C_v Z^s C_o'),
// and the gradient is
//   dE/dx = sum D1 h^x + sum Gamma (mn|ls)^x + sum W S^x + V_nn^x
// with the relaxed one-particle density D1 = sum_s C_s D_s C_s' - DZ_a - DZ_b, the AO two-particle
// density of G_aa, G_bb, G_ab and the Z-vector terms (Coulomb -DZ_t P_t, exchange +DZ_s P_s), and the
// energy-weighted density w_s = -X^s/2 + (X^s_ai - X^s_ia)/2 + Z^s_ai e_i (virtual-occupied block)
// + C_o' [J(DZ_t) - K(DZ_s)] C_o (occupied block), W = sum_s C_s w_s C_s' symmetrised. For a closed
// shell (alpha = beta orbitals) this is the RHF expression of correlated_gradient.h.
#include "macqc/uhf_correlated.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

#include "macqc/gto_derivatives.h"
#include "macqc/gto_integrals.h"
#include "macqc/job.h"
#include "macqc/scf_derivatives.h"

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

Matrix mul(const Matrix& A, bool ta, const Matrix& B, bool tb, int threads) {
    const std::size_t m = ta ? A.cols() : A.rows();
    const std::size_t k = ta ? A.rows() : A.cols();
    const std::size_t n = tb ? B.rows() : B.cols();
    Matrix C(m, n);
    if (m && n && k) gemm(ta, tb, m, n, k, 1.0, A.data(), A.cols(), B.data(), B.cols(), 0.0, C.data(), n, threads);
    return C;
}

Matrix columns(const Matrix& C, std::size_t c0, std::size_t nc) {
    Matrix out(C.rows(), nc);
    for (std::size_t i = 0; i < C.rows(); ++i)
        for (std::size_t j = 0; j < nc; ++j) out(i, j) = C(i, c0 + j);
    return out;
}

Matrix symmetrized(const Matrix& A) {
    Matrix S(A.rows(), A.cols());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) S(i, j) = 0.5 * (A(i, j) + A(j, i));
    return S;
}

void add(Matrix& y, double a, const Matrix& x) {
    for (std::size_t i = 0; i < y.rows(); ++i)
        for (std::size_t j = 0; j < y.cols(); ++j) y(i, j) += a * x(i, j);
}

// (pq|rs) of the spin pair (sp, sq) from the MO integrals of the reference
struct UInts {
    const MOIntegrals* aa;
    const MOIntegrals* bb;
    const MOIntegrals* ab;
    double operator()(int sp, int sq, std::size_t p, std::size_t q, std::size_t r, std::size_t s) const {
        if (sp == 0 && sq == 0) return (*aa)(p, q, r, s);
        if (sp == 1 && sq == 1) return (*bb)(p, q, r, s);
        if (sp == 0) return ab->mixed ? ab->mixed_value(p, q, r, s) : (*ab)(p, q, r, s);
        return ab->mixed ? ab->mixed_value(r, s, p, q) : (*ab)(r, s, p, q);
    }
};

// 8-fold symmetrisation of a full n^4 array into a packed one
Packed4 pack8(const std::vector<double>& G, std::size_t n) {
    Packed4 P(n);
    for (std::size_t p = 0; p < n; ++p)
        for (std::size_t q = 0; q < n; ++q)
            for (std::size_t r = 0; r < n; ++r)
                for (std::size_t s = 0; s < n; ++s) {
                    const double v = G[((p * n + q) * n + r) * n + s];
                    if (v != 0.0) P.add_raw(p, q, r, s, v);
                }
    return P;
}

// 4-fold symmetrisation (p <-> q, r <-> s) of the alpha-beta density, full layout
std::vector<double> sym4(const std::vector<double>& G, std::size_t na, std::size_t nb) {
    std::vector<double> out(G.size(), 0.0);
    for (std::size_t p = 0; p < na; ++p)
        for (std::size_t q = 0; q < na; ++q)
            for (std::size_t r = 0; r < nb; ++r)
                for (std::size_t s = 0; s < nb; ++s) {
                    auto at = [&](std::size_t a, std::size_t b, std::size_t c, std::size_t d) {
                        return G[((a * na + b) * nb + c) * nb + d];
                    };
                    out[((p * na + q) * nb + r) * nb + s] = 0.25 * (at(p, q, r, s) + at(q, p, r, s) + at(p, q, s, r) + at(q, p, s, r));
                }
    return out;
}

// AO two-particle density of the alpha-beta part (bra Ca, ket Cb), symmetrised over the exchange of
// the pairs, in the packed layout of FockBuilder::packed().
void back_transform_mixed(const std::vector<double>& Gab, const Matrix& Ca, const Matrix& Cb, std::vector<double>& out,
                          int nt) {
    const std::size_t N = Ca.rows(), na = Ca.cols(), nb = Cb.cols();
    const std::size_t NP = N * (N + 1) / 2, npb = nb * (nb + 1) / 2;
    std::vector<double> half(NP * npb, 0.0);  // [mu nu][rs]
    Matrix M(na, na);
    for (std::size_t r = 0, rs = 0; r < nb; ++r)
        for (std::size_t s = 0; s <= r; ++s, ++rs) {
            check_cancelled();
            const double w = r == s ? 1.0 : 2.0;  // r s and s r (G symmetric in r, s)
            for (std::size_t p = 0; p < na; ++p)
                for (std::size_t q = 0; q < na; ++q) M(p, q) = w * Gab[((p * na + q) * nb + r) * nb + s];
            const Matrix A = mul(Ca, false, mul(M, false, Ca, true, nt), false, nt);
            for (std::size_t m = 0, mn = 0; m < N; ++m)
                for (std::size_t k = 0; k <= m; ++k, ++mn) half[mn * npb + rs] = A(m, k);
        }
    // full AO pair matrix Gamma(mn, ls), then the symmetrised packed value
    std::vector<double> full(NP * NP, 0.0);
    Matrix M2(nb, nb);
    for (std::size_t mn = 0; mn < NP; ++mn) {
        check_cancelled();
        const double* h = half.data() + mn * npb;
        for (std::size_t r = 0, rs = 0; r < nb; ++r)
            for (std::size_t s = 0; s <= r; ++s, ++rs) {
                const double x = (r == s ? 1.0 : 0.5) * h[rs];
                M2(r, s) = M2(s, r) = x;
            }
        const Matrix B = mul(Cb, false, mul(M2, false, Cb, true, nt), false, nt);
        for (std::size_t l = 0, ls = 0; l < N; ++l)
            for (std::size_t s = 0; s <= l; ++s, ++ls) full[mn * NP + ls] = B(l, s);
    }
    for (std::size_t mn = 0; mn < NP; ++mn)
        for (std::size_t ls = 0; ls <= mn; ++ls) out[mn * (mn + 1) / 2 + ls] += 0.5 * (full[mn * NP + ls] + full[ls * NP + mn]);
}

}  // namespace

UHFCorrelatedGradient uhf_correlated_gradient(const Molecule& mol, const BasisSet& basis, const RHFResult& scf,
                                             const MOIntegrals& aa, const MOIntegrals& bb, const MOIntegrals& ab,
                                             const SODensities& d, int threads, double memory_mb,
                                             const CPHFCallback& progress) {
    if (scf.reference != "UHF") throw std::invalid_argument("uhf_correlated_gradient: UHF reference needed");
    if (scf.removed_functions > 0) throw std::runtime_error("analytic correlated gradient with removed near-linear dependencies");
    const auto t0 = Clock::now();
    const int nt = resolve_threads(threads);
    const int nat = static_cast<int>(mol.atoms.size());
    const Matrix* Cs[2] = {&scf.C, &scf.Cb};
    const std::vector<double>* es[2] = {&scf.orbital_energies, &scf.orbital_energies_b};
    const std::size_t N = scf.C.rows();
    const std::size_t n[2] = {scf.C.cols(), scf.Cb.cols()};
    const std::size_t o[2] = {static_cast<std::size_t>(scf.n_alpha), static_cast<std::size_t>(scf.n_beta)};
    const UInts g{&aa, &bb, &ab};
    UHFCorrelatedGradient res;

    // symmetrised densities
    Matrix Ds[2] = {symmetrized(d.D[0]), symmetrized(d.D[1])};
    const Packed4 Gaa = pack8(d.Gaa, n[0]), Gbb = pack8(d.Gbb, n[1]);
    const std::vector<double> Gab = sym4(d.Gab, n[0], n[1]);
    // the energy from the densities
    {
        double e = scf.nuclear_repulsion;
        for (int s = 0; s < 2; ++s) {
            const Matrix& h = s ? bb.h : aa.h;
            for (std::size_t p = 0; p < n[s]; ++p)
                for (std::size_t q = 0; q < n[s]; ++q) e += Ds[s](p, q) * h(p, q);
        }
        for (std::size_t p = 0; p < n[0]; ++p)
            for (std::size_t q = 0; q < n[0]; ++q)
                for (std::size_t r = 0; r < n[0]; ++r)
                    for (std::size_t s = 0; s < n[0]; ++s) e += Gaa.get(p, q, r, s) * aa(p, q, r, s);
        for (std::size_t p = 0; p < n[1]; ++p)
            for (std::size_t q = 0; q < n[1]; ++q)
                for (std::size_t r = 0; r < n[1]; ++r)
                    for (std::size_t s = 0; s < n[1]; ++s) e += Gbb.get(p, q, r, s) * bb(p, q, r, s);
        for (std::size_t p = 0; p < n[0]; ++p)
            for (std::size_t q = 0; q < n[0]; ++q)
                for (std::size_t r = 0; r < n[1]; ++r)
                    for (std::size_t s = 0; s < n[1]; ++s)
                        e += Gab[((p * n[0] + q) * n[1] + r) * n[1] + s] * g(0, 1, p, q, r, s);
        res.density_energy = e;
    }

    // ---- orbital Lagrangians
    report_progress(JobStage::Gradient, 1);
    Matrix X[2];
    for (int s = 0; s < 2; ++s) {
        const Matrix& h = s ? bb.h : aa.h;
        X[s] = mul(h, false, Ds[s], false, nt);
        for (std::size_t t = 0; t < n[s]; ++t)
            for (std::size_t p = 0; p < n[s]; ++p) X[s](t, p) *= 2.0;
    }
    for (int s = 0; s < 2; ++s) {
        const std::size_t ns = n[s];
        const Packed4& Gs = s ? Gbb : Gaa;
        for (std::size_t t = 0; t < ns; ++t) {
            check_cancelled();
            for (std::size_t p = 0; p < ns; ++p) {
                double v = 0.0;
                for (std::size_t q = 0; q < ns; ++q)
                    for (std::size_t r = 0; r < ns; ++r)
                        for (std::size_t u = 0; u < ns; ++u) v += g(s, s, t, q, r, u) * Gs.get(p, q, r, u);
                X[s](t, p) += 4.0 * v;
            }
        }
    }
    for (std::size_t t = 0; t < n[0]; ++t) {
        check_cancelled();
        for (std::size_t p = 0; p < n[0]; ++p) {
            double v = 0.0;
            for (std::size_t q = 0; q < n[0]; ++q)
                for (std::size_t r = 0; r < n[1]; ++r)
                    for (std::size_t u = 0; u < n[1]; ++u) v += g(0, 1, t, q, r, u) * Gab[((p * n[0] + q) * n[1] + r) * n[1] + u];
            X[0](t, p) += 2.0 * v;
        }
    }
    for (std::size_t t = 0; t < n[1]; ++t) {
        check_cancelled();
        for (std::size_t r = 0; r < n[1]; ++r) {
            double v = 0.0;
            for (std::size_t p = 0; p < n[0]; ++p)
                for (std::size_t q = 0; q < n[0]; ++q)
                    for (std::size_t u = 0; u < n[1]; ++u) v += g(0, 1, p, q, t, u) * Gab[((p * n[0] + q) * n[1] + r) * n[1] + u];
            X[1](t, r) += 2.0 * v;
        }
    }
    for (int s = 0; s < 2; ++s)
        for (std::size_t p = 0; p < n[s]; ++p)
            for (std::size_t q = 0; q < p; ++q)
                if ((p < o[s]) == (q < o[s]))
                    res.orbital_lagrangian_oo_vv = std::max(res.orbital_lagrangian_oo_vv, std::fabs(X[s](p, q) - X[s](q, p)));

    // ---- Z-vector (coupled-perturbed UHF)
    const auto t1 = Clock::now();
    FockBuilder fock(basis, nt, memory_mb);
    const std::size_t v[2] = {n[0] - o[0], n[1] - o[1]};
    Matrix Co[2], Cv[2];
    for (int s = 0; s < 2; ++s) {
        Co[s] = columns(*Cs[s], 0, o[s]);
        Cv[s] = columns(*Cs[s], o[s], v[s]);
    }
    std::vector<double> rhs, diag;
    for (int s = 0; s < 2; ++s)
        for (std::size_t a = 0; a < v[s]; ++a)
            for (std::size_t i = 0; i < o[s]; ++i) {
                rhs.push_back(X[s](o[s] + a, i) - X[s](i, o[s] + a));
                diag.push_back((*es[s])[o[s] + a] - (*es[s])[i]);
            }
    auto unpack = [&](const std::vector<double>& z, Matrix* Z) {
        std::size_t pos = 0;
        for (int s = 0; s < 2; ++s) {
            Z[s] = Matrix(v[s], o[s]);
            for (std::size_t a = 0; a < v[s]; ++a)
                for (std::size_t i = 0; i < o[s]; ++i) Z[s](a, i) = z[pos++];
        }
    };
    auto densities_of = [&](const Matrix* Z, Matrix* DZ) {
        for (int s = 0; s < 2; ++s) DZ[s] = symmetrized(mul(mul(Cv[s], false, Z[s], false, nt), false, Co[s], true, nt));
    };
    auto apply = [&](const std::vector<std::vector<double>>& batch) {
        std::vector<Matrix> dens;
        std::vector<std::array<Matrix, 2>> Zs(batch.size());
        for (std::size_t k = 0; k < batch.size(); ++k) {
            Matrix Z[2], DZ[2];
            unpack(batch[k], Z);
            densities_of(Z, DZ);
            Zs[k] = {Z[0], Z[1]};
            // P^U = C_v U C_o' + transpose = 2 DZ
            for (int s = 0; s < 2; ++s) {
                Matrix P = DZ[s];
                for (std::size_t i = 0; i < N; ++i)
                    for (std::size_t j = 0; j < N; ++j) P(i, j) *= 2.0;
                dens.push_back(P);
            }
        }
        std::vector<Matrix> J, K;
        fock.build_jk(dens, J, K);
        std::vector<std::vector<double>> out(batch.size());
        for (std::size_t k = 0; k < batch.size(); ++k) {
            Matrix Jt = J[2 * k];
            add(Jt, 1.0, J[2 * k + 1]);
            for (int s = 0; s < 2; ++s) {
                Matrix G = Jt;
                add(G, -1.0, K[2 * k + static_cast<std::size_t>(s)]);
                const Matrix Gvo = mul(mul(Cv[s], true, G, false, nt), false, Co[s], false, nt);
                for (std::size_t a = 0; a < v[s]; ++a)
                    for (std::size_t i = 0; i < o[s]; ++i)
                        out[k].push_back(((*es[s])[o[s] + a] - (*es[s])[i]) * Zs[k][static_cast<std::size_t>(s)](a, i) + Gvo(a, i));
            }
        }
        return out;
    };
    const LinearBatchSolution zs = solve_symmetric_batch(apply, diag, {rhs}, 1.0e-10, 100, progress);
    res.zvector_converged = zs.converged;
    res.zvector_iterations = zs.iterations;
    res.zvector_residual = zs.residual;
    Matrix Z[2], DZ[2];
    unpack(zs.x[0], Z);
    densities_of(Z, DZ);
    Matrix GDZ[2];
    {
        std::vector<Matrix> J, K;
        fock.build_jk({DZ[0], DZ[1]}, J, K);
        Matrix Jt = J[0];
        add(Jt, 1.0, J[1]);
        for (int s = 0; s < 2; ++s) {
            GDZ[s] = Jt;
            add(GDZ[s], -1.0, K[static_cast<std::size_t>(s)]);
        }
    }
    res.zvector_seconds = seconds_since(t1);

    // ---- energy-weighted density, relaxed density
    Matrix Wao(N, N), D1(N, N);
    for (int s = 0; s < 2; ++s) {
        const std::size_t ns = n[s], os = o[s];
        Matrix w(ns, ns);
        for (std::size_t t = 0; t < ns; ++t)
            for (std::size_t p = 0; p < ns; ++p) w(t, p) = -0.5 * X[s](t, p);
        for (std::size_t a = 0; a < v[s]; ++a)
            for (std::size_t i = 0; i < os; ++i)
                w(os + a, i) += 0.5 * (X[s](os + a, i) - X[s](i, os + a)) + Z[s](a, i) * (*es[s])[i];
        const Matrix Goo = mul(mul(Co[s], true, GDZ[s], false, nt), false, Co[s], false, nt);
        for (std::size_t k = 0; k < os; ++k)
            for (std::size_t l = 0; l < os; ++l) w(k, l) += Goo(k, l);
        add(Wao, 1.0, symmetrized(mul(mul(*Cs[s], false, w, false, nt), false, *Cs[s], true, nt)));
        add(D1, 1.0, mul(mul(*Cs[s], false, Ds[s], false, nt), false, *Cs[s], true, nt));
        add(D1, -1.0, DZ[s]);
    }
    D1 = symmetrized(D1);
    res.relaxed_density_ao = D1;
    res.relaxed_dipole = dipole_moment(mol, basis, D1);

    // ---- AO two-particle density
    const auto t2 = Clock::now();
    report_progress(JobStage::Gradient, 2);
    std::vector<double> gao = back_transform_packed(Gaa, scf.C, nt);
    {
        const std::vector<double> gb = back_transform_packed(Gbb, scf.Cb, nt);
        for (std::size_t k = 0; k < gao.size(); ++k) gao[k] += gb[k];
    }
    back_transform_mixed(Gab, scf.C, scf.Cb, gao, nt);
    res.backtransform_seconds = seconds_since(t2);
    TwoParticleDensity gamma;
    gamma.packed = &gao;
    {
        SeparableDensity c;
        c.A = DZ[0];
        add(c.A, 1.0, DZ[1]);
        c.B = scf.P;
        c.coulomb = -1.0;
        c.exchange = 0.0;
        gamma.separable.push_back(c);
        for (int s = 0; s < 2; ++s) {
            SeparableDensity x;
            x.A = DZ[s];
            x.B = s ? scf.Pb : scf.Pa;
            x.coulomb = 0.0;
            x.exchange = 1.0;
            gamma.separable.push_back(x);
        }
    }

    // ---- derivative integrals
    const auto t3 = Clock::now();
    report_progress(JobStage::Gradient, 3);
    res.nuclear = nuclear_repulsion_gradient(mol);
    res.one_electron = core_hamiltonian_gradient(basis, mol, D1, nt);
    res.two_electron = eri_gradient(basis, nat, gamma, nt);
    res.overlap = overlap_gradient(basis, nat, Wao, nt);
    res.integral_seconds = seconds_since(t3);
    res.gradient.assign(static_cast<std::size_t>(3 * nat), 0.0);
    for (std::size_t x = 0; x < res.gradient.size(); ++x)
        res.gradient[x] = res.nuclear[x] + res.one_electron[x] + res.two_electron[x] + res.overlap[x];
    res.seconds = seconds_since(t0);
    return res;
}

}  // namespace macqc
