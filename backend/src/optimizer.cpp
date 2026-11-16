#include "macqc/optimizer.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <stdexcept>

#include "macqc/gto_derivatives.h"
#include "macqc/hessian.h"
#include "macqc/job.h"
#include "macqc/rhf_derivatives.h"

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }
constexpr double kPi = 3.14159265358979323846;

std::string format(const char* fmt, ...) {
    char buf[1024];
    va_list args;
    va_start(args, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, args);
    va_end(args);
    return buf;
}

int period_row(int Z) { return Z <= 2 ? 0 : (Z <= 10 ? 1 : 2); }

double distance(const Molecule& m, int a, int b) {
    double s = 0.0;
    for (std::size_t k = 0; k < 3; ++k) {
        const double d = m.atoms[static_cast<std::size_t>(a)].position[k] - m.atoms[static_cast<std::size_t>(b)].position[k];
        s += d * d;
    }
    return std::sqrt(s);
}

Matrix mat_mul(const Matrix& A, const Matrix& B) {
    Matrix C(A.rows(), B.cols());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t k = 0; k < A.cols(); ++k) {
            const double a = A(i, k);
            if (a == 0.0) continue;
            for (std::size_t j = 0; j < B.cols(); ++j) C(i, j) += a * B(k, j);
        }
    return C;
}

double norm(const std::vector<double>& v) {
    double s = 0.0;
    for (double x : v) s += x * x;
    return std::sqrt(s);
}

// Symmetric difference |R(O-H1) - R(O-H2)| of a triatomic with the heavy atom first.
double c2v_deviation(const Molecule& m) {
    if (m.atoms.size() != 3 || m.custom) return 0.0;
    return std::fabs(distance(m, 0, 1) - distance(m, 0, 2));
}

}  // namespace

ConvergenceCriteria convergence_criteria(bool tight) {
    ConvergenceCriteria c;
    if (tight) {
        c.max_force = 1.5e-5;
        c.rms_force = 1.0e-5;
        c.max_displacement = 6.0e-5;
        c.rms_displacement = 4.0e-5;
        c.energy = 1.0e-8;
    }
    return c;
}

std::vector<double> lindh_model_hessian(const Molecule& m, const InternalCoordinates& ic) {
    // Lindh, Bernhardsson, Karlström & Malmqvist, Chem. Phys. Lett. 241, 423 (1995).
    static const double alpha[3][3] = {{1.0000, 0.3949, 0.3949}, {0.3949, 0.2800, 0.2800}, {0.3949, 0.2800, 0.2800}};
    static const double rref[3][3] = {{1.35, 2.10, 2.53}, {2.10, 2.87, 3.40}, {2.53, 3.40, 3.40}};
    auto rho = [&](int a, int b) {
        const int ra = period_row(m.atoms[static_cast<std::size_t>(a)].Z), rb = period_row(m.atoms[static_cast<std::size_t>(b)].Z);
        const double r = distance(m, a, b);
        return std::exp(alpha[ra][rb] * (rref[ra][rb] * rref[ra][rb] - r * r));
    };
    std::vector<double> h;
    for (const auto& c : ic.coordinates()) {
        double k = 0.0;
        if (c.type == InternalType::Stretch) k = 0.45 * rho(c.atoms[0], c.atoms[1]);
        else if (c.type == InternalType::Dihedral)
            k = 0.005 * rho(c.atoms[0], c.atoms[1]) * rho(c.atoms[1], c.atoms[2]) * rho(c.atoms[2], c.atoms[3]);
        else k = 0.15 * rho(c.atoms[0], c.atoms[1]) * rho(c.atoms[1], c.atoms[2]);
        h.push_back(std::max(k, 0.005));
    }
    return h;
}

std::vector<double> schlegel_model_hessian(const Molecule& m, const InternalCoordinates& ic) {
    // Schlegel, Theor. Chim. Acta 66, 333 (1984): F_r = 1.734 / (r - B)^3; bends 0.160 (H) / 0.250.
    static const double Btab[3][3] = {{-0.244, 0.352, 0.660}, {0.352, 1.085, 1.522}, {0.660, 1.522, 2.068}};
    std::vector<double> h;
    for (const auto& c : ic.coordinates()) {
        if (c.type == InternalType::Stretch) {
            const int ra = period_row(m.atoms[static_cast<std::size_t>(c.atoms[0])].Z),
                      rb = period_row(m.atoms[static_cast<std::size_t>(c.atoms[1])].Z);
            const double d = distance(m, c.atoms[0], c.atoms[1]) - Btab[ra][rb];
            h.push_back(std::max(0.005, 1.734 / (d * d * d)));
        } else if (c.type == InternalType::Dihedral) {
            // F = A - B (R - R_cov) of the central bond, A = 0.0023, B = 0.07
            const double r = distance(m, c.atoms[1], c.atoms[2]);
            const double rcov = covalent_radius(m.atoms[static_cast<std::size_t>(c.atoms[1])].Z) +
                                covalent_radius(m.atoms[static_cast<std::size_t>(c.atoms[2])].Z);
            h.push_back(std::max(0.0023, 0.0023 - 0.07 * (r - rcov)));
        } else {
            const bool hydrogen = m.atoms[static_cast<std::size_t>(c.atoms[0])].Z == 1 ||
                                  m.atoms[static_cast<std::size_t>(c.atoms[2])].Z == 1;
            h.push_back(hydrogen ? 0.160 : 0.250);
        }
    }
    return h;
}

std::vector<GeometricParameter> geometric_parameters(const Molecule& mol) {
    std::vector<GeometricParameter> out;
    if (mol.atoms.size() < 2) return out;  // 9.0: a single atom (custom molecule)
    const InternalCoordinates ic(mol);
    const std::vector<double> q = ic.values(mol);
    for (std::size_t i = 0; i < ic.size(); ++i) {
        const InternalCoordinate& c = ic.coordinates()[i];
        if (c.type == InternalType::Stretch) {
            out.push_back({c.definition(mol), q[i], false});
        } else if (c.type == InternalType::Bend || c.type == InternalType::Dihedral) {
            out.push_back({c.definition(mol), q[i] * 180.0 / kPi, true});
        } else if (c.component == 0) {
            InternalCoordinate b = c;
            b.type = InternalType::Bend;
            const std::string def = "A(" + mol.atoms[static_cast<std::size_t>(c.atoms[0])].symbol + std::to_string(c.atoms[0] + 1) +
                                    "-" + mol.atoms[static_cast<std::size_t>(c.atoms[1])].symbol + std::to_string(c.atoms[1] + 1) +
                                    "-" + mol.atoms[static_cast<std::size_t>(c.atoms[2])].symbol + std::to_string(c.atoms[2] + 1) + ")";
            // angle from the bond vectors
            double u[3], v[3], uu = 0, vv = 0, uv = 0;
            for (std::size_t k = 0; k < 3; ++k) {
                u[k] = mol.atoms[static_cast<std::size_t>(c.atoms[0])].position[k] - mol.atoms[static_cast<std::size_t>(c.atoms[1])].position[k];
                v[k] = mol.atoms[static_cast<std::size_t>(c.atoms[2])].position[k] - mol.atoms[static_cast<std::size_t>(c.atoms[1])].position[k];
                uu += u[k] * u[k];
                vv += v[k] * v[k];
                uv += u[k] * v[k];
            }
            out.push_back({def, std::acos(std::max(-1.0, std::min(1.0, uv / std::sqrt(uu * vv)))) * 180.0 / kPi, true});
        }
    }
    return out;
}

namespace {

struct PrimitiveConstraint {
    InternalCoordinate c;
    double target = 0.0;
};

// A bend frozen at 180 degrees becomes its two linear-bend components, both 0.
std::vector<PrimitiveConstraint> expand_constraints(const std::vector<CoordinateConstraint>& cs) {
    std::vector<PrimitiveConstraint> out;
    for (const auto& c : cs) {
        if (c.coordinate.type == InternalType::Bend && std::fabs(c.target - kPi) < 1.0e-9) {
            for (int k = 0; k < 2; ++k) {
                PrimitiveConstraint p;
                p.c = c.coordinate;
                p.c.type = InternalType::LinearBend;
                p.c.component = k;
                p.target = 0.0;
                out.push_back(p);
            }
        } else {
            out.push_back({c.coordinate, c.target});
        }
    }
    return out;
}

std::vector<InternalCoordinate> required_coordinates(const std::vector<PrimitiveConstraint>& pcs) {
    std::vector<InternalCoordinate> r;
    for (const auto& p : pcs) r.push_back(p.c);
    return r;
}

// Projector onto the space that keeps the constraints, and the change that drives them to
// their targets: P~ = P - P C (C' P C)^- C' P, dq_c = P C (C' P C)^- (target - q_c).
struct ConstraintProjection {
    Matrix Pt;
    std::vector<double> dq_c;
    std::vector<int> index;
    double deviation = 0.0;
};

ConstraintProjection project_constraints(const InternalCoordinates& ic, const std::vector<double>& q, const Matrix& P,
                                         const std::vector<PrimitiveConstraint>& pcs, double scale = 1.0) {
    ConstraintProjection r;
    const std::size_t nq = q.size(), nc = pcs.size();
    r.Pt = P;
    r.dq_c.assign(nq, 0.0);
    if (nc == 0) return r;
    std::vector<double> delta(nc);
    for (std::size_t k = 0; k < nc; ++k) {
        const int i = ic.find(pcs[k].c);
        if (i < 0) throw std::runtime_error("internal error: a constrained coordinate is missing from the coordinate set");
        r.index.push_back(i);
        double d = pcs[k].target - q[static_cast<std::size_t>(i)];
        if (pcs[k].c.type == InternalType::Dihedral) d = wrap_angle(d);
        delta[k] = d;
        r.deviation = std::max(r.deviation, std::fabs(d));
    }
    Matrix M(nc, nc);
    for (std::size_t a = 0; a < nc; ++a)
        for (std::size_t b = 0; b < nc; ++b)
            M(a, b) = P(static_cast<std::size_t>(r.index[a]), static_cast<std::size_t>(r.index[b]));
    const Matrix Mi = generalized_inverse(M, 1.0e-10);
    // PC (nq x nc) = columns of P
    Matrix PCM(nq, nc);  // P C M^-
    for (std::size_t i = 0; i < nq; ++i)
        for (std::size_t b = 0; b < nc; ++b) {
            double s = 0.0;
            for (std::size_t a = 0; a < nc; ++a) s += P(i, static_cast<std::size_t>(r.index[a])) * Mi(a, b);
            PCM(i, b) = s;
        }
    for (std::size_t i = 0; i < nq; ++i)
        for (std::size_t j = 0; j < nq; ++j) {
            double s = 0.0;
            for (std::size_t b = 0; b < nc; ++b) s += PCM(i, b) * P(static_cast<std::size_t>(r.index[b]), j);
            r.Pt(i, j) -= s;
        }
    for (std::size_t i = 0; i < nq; ++i)
        for (std::size_t b = 0; b < nc; ++b) r.dq_c[i] += PCM(i, b) * scale * delta[b];
    return r;
}

struct Projectors {
    Matrix B, G, Gi, P;
    int rank = 0;
};

// H_q = G^- B (H_x - K) B' G^-, K = sum_i g_q,i d2q_i/dx2 (exact also away from a stationary point).
Matrix internal_hessian(const InternalCoordinates& ic, const Molecule& x, const Matrix& Hx, const Matrix& B,
                        const Matrix& Gi, const std::vector<double>& gq0) {
    const std::size_t nq = ic.size(), n3 = 3 * x.atoms.size();
    Matrix Hc = Hx;
    const Matrix K = ic.b_derivative_contracted(x, gq0);
    for (std::size_t a = 0; a < n3; ++a)
        for (std::size_t b = 0; b < n3; ++b) Hc(a, b) -= K(a, b);
    const Matrix GB = mat_mul(Gi, B);
    const Matrix T = mat_mul(GB, Hc);
    Matrix H(nq, nq);
    for (std::size_t i = 0; i < nq; ++i)
        for (std::size_t j = 0; j < nq; ++j) {
            double s = 0.0;
            for (std::size_t k = 0; k < n3; ++k) s += T(i, k) * GB(j, k);
            H(i, j) = s;
        }
    for (std::size_t i = 0; i < nq; ++i)
        for (std::size_t j = 0; j < i; ++j) H(i, j) = H(j, i) = 0.5 * (H(i, j) + H(j, i));
    return H;
}

Projectors projectors(const InternalCoordinates& ic, const Molecule& x) {
    Projectors p;
    const std::size_t nq = ic.size(), n3 = 3 * x.atoms.size();
    p.B = ic.wilson_b(x);
    p.G = Matrix(nq, nq);
    for (std::size_t i = 0; i < nq; ++i)
        for (std::size_t j = 0; j < nq; ++j) {
            double s = 0.0;
            for (std::size_t k = 0; k < n3; ++k) s += p.B(i, k) * p.B(j, k);
            p.G(i, j) = s;
        }
    p.Gi = generalized_inverse(p.G, 1.0e-8, &p.rank);
    p.P = mat_mul(p.G, p.Gi);
    return p;
}

}  // namespace

Molecule apply_constraints(const Molecule& mol, const std::vector<CoordinateConstraint>& constraints, double* residual) {
    // A bend set to 180 degrees goes to 179 degrees as a bend first (no change of the other
    // coordinates), then to 180 through its linear-bend components, which are redundant with the
    // other coordinates and would otherwise move them.
    std::vector<CoordinateConstraint> first = constraints;
    bool two_phases = false;
    for (auto& c : first)
        if (c.coordinate.type == InternalType::Bend && std::fabs(c.target - kPi) < 1.0e-9 &&
            internal_value(mol, c.coordinate) < 179.0 * kPi / 180.0) {
            c.target = 179.0 * kPi / 180.0;
            two_phases = true;
        }
    const Molecule start = two_phases ? apply_constraints(mol, first, nullptr) : mol;
    const std::vector<PrimitiveConstraint> pcs = expand_constraints(constraints);
    const std::vector<InternalCoordinate> required = required_coordinates(pcs);
    Molecule x = start;
    double dev = 0.0;
    for (int pass = 0; pass < 400; ++pass) {
        const InternalCoordinates ic(x, required);
        const std::vector<double> q = ic.values(x);
        const Projectors pr = projectors(ic, x);
        // largest remaining change per kind -> substep scale
        double ds = 0.0, da = 0.0;
        dev = 0.0;
        for (const auto& pc : pcs) {
            const int i = ic.find(pc.c);
            if (i < 0) throw std::runtime_error("internal error: a constrained coordinate is missing");
            double d = pc.target - q[static_cast<std::size_t>(i)];
            if (pc.c.type == InternalType::Dihedral) d = wrap_angle(d);
            dev = std::max(dev, std::fabs(d));
            if (pc.c.type == InternalType::Stretch) ds = std::max(ds, std::fabs(d));
            else da = std::max(da, std::fabs(d));
        }
        if (dev < 1.0e-10) break;
        double scale = 1.0;
        if (ds > 0.2) scale = std::min(scale, 0.2 / ds);
        if (da > 20.0 * kPi / 180.0) scale = std::min(scale, 20.0 * kPi / 180.0 / da);
        const ConstraintProjection cp = project_constraints(ic, q, pr.P, pcs, scale);
        const BackTransform bt = back_transform(ic, x, cp.dq_c);
        x = bt.molecule;
    }
    if (residual) *residual = dev;
    return x;
}

OptimizationResult optimize_geometry(const Molecule& start, const BasisSet& basis, const std::string& method,
                                     const ElectronicOptions& eo, const OptimizationOptions& oo,
                                     const ElectronicCallbacks& cb, const OptimizationCallback& on_step) {
    const auto t0 = Clock::now();
    OptimizationResult res;
    res.criteria = convergence_criteria(oo.tight);
    res.initial = start;
    res.initial_parameters = geometric_parameters(start);
    res.symmetric_start = start.atoms.size() == 3 && !start.custom && c2v_deviation(start) < 1.0e-8;
    const ConvergenceCriteria& cc = res.criteria;
    const std::vector<PrimitiveConstraint> pcs = expand_constraints(oo.constraints);
    const std::vector<InternalCoordinate> required = required_coordinates(pcs);
    const bool constrained = !pcs.empty();
    const bool ts = oo.transition_state;
    if (ts && constrained) throw std::invalid_argument("a transition-state search with frozen coordinates is not available");
    if (ts && oo.ts_hessian != "calcfc" && oo.ts_hessian != "calcall")
        throw std::invalid_argument("ts_hessian must be calcfc or calcall");
    res.transition_state = ts;
    res.ts_hessian = ts ? oo.ts_hessian : "";
    std::vector<double> followed_prev;  // the mode followed at the previous step (transition states)
    const double constraint_threshold = 0.01 * cc.max_displacement;

    Molecule x = start;
    InternalCoordinates ic(x, required);
    Matrix H;
    bool need_initial = true;
    std::vector<double> q_prev, g_prev;
    double e_prev = 0.0, predicted_prev = 0.0, step_prev = 0.0;
    double trust = oo.trust_radius;
    Matrix orbitals;
    const Matrix* guess = oo.guess_orbitals;
    int removed0 = -1;
    std::vector<double> negative_direction;  // of the initial analytic Hessian, if it had a negative eigenvalue
    double negative_curvature = 0.0;
    bool first_build = true;

    for (int step = 1; step <= oo.max_steps; ++step) {
        check_cancelled();
        const auto ts_start = Clock::now();
        OptimizationStep rec;
        rec.step = step;
        rec.geometry = x;
        if (oo.on_step_start) oo.on_step_start(step, x);
        const BasisSet bx = basis_at_geometry(basis, x);
        GradientResult g = run_gradient(x, bx, method, eo, guess, cb);
        rec.warnings = g.point.warnings;
        if (!g.point.converged) {
            res.message = format("The optimisation stopped at step %d: the %s calculation did not converge.", step,
                                 method.c_str());
            res.final_point = g;
            res.final_geometry = x;
            res.final_energy = g.point.energy;
            res.steps = step;
            break;
        }
        if (removed0 < 0) {
            removed0 = g.point.rhf.removed_functions;
        } else if (g.point.rhf.removed_functions != removed0) {
            res.removed_functions_changed = true;
            rec.warnings.push_back(format(
                "The number of removed near-linear dependencies changed from %d to %d: the potential energy surface is "
                "not continuous here.",
                removed0, g.point.rhf.removed_functions));
            removed0 = g.point.rhf.removed_functions;
        }
        orbitals = g.point.rhf.C;
        guess = &orbitals;
        rec.energy = g.point.energy;
        rec.delta_e = step == 1 ? 0.0 : rec.energy - e_prev;
        rec.gradient_description = g.description;
        rec.cartesian_gradient = g.gradient;

        // ---- internal coordinates, gradient
        if (ic.needs_rebuild(x)) {
            ic = InternalCoordinates(x, required);
            need_initial = true;
            q_prev.clear();
            rec.warnings.push_back("The internal coordinates were rebuilt (a bend crossed the linear threshold).");
        }
        const std::vector<double> q = ic.values(x);
        const std::size_t nq = q.size(), n3 = 3 * x.atoms.size();
        const Projectors pr = projectors(ic, x);
        const Matrix& B = pr.B;
        const Matrix& Gi = pr.Gi;
        const Matrix& P = pr.P;
        const std::vector<double> gq0 = internal_gradient(B, Gi, g.gradient);
        std::vector<double> gq(nq, 0.0);  // projected gradient of the full surface (Hessian updates)
        for (std::size_t i = 0; i < nq; ++i)
            for (std::size_t j = 0; j < nq; ++j) gq[i] += P(i, j) * gq0[j];
        const ConstraintProjection cp = project_constraints(ic, q, P, pcs);
        const Matrix& Pt = cp.Pt;  // = P without constraints
        std::vector<double> gt(nq, 0.0);  // gradient within the constrained space (step, convergence)
        for (std::size_t i = 0; i < nq; ++i)
            for (std::size_t j = 0; j < nq; ++j) gt[i] += Pt(i, j) * gq[j];
        for (const auto& c : ic.coordinates()) {
            rec.names.push_back(c.name());
            rec.definitions.push_back(c.definition(x));
        }
        rec.q = q;
        rec.gq = constrained ? gt : gq;
        rec.constrained = cp.index;
        rec.constraint_deviation = cp.deviation;

        // ---- Hessian: computed (transition states), initial (carried over, analytic or model) or updated
        const bool ts_compute = ts && (need_initial || q_prev.size() != nq || oo.ts_hessian == "calcall");
        if (ts_compute) {
            HessianOptions ho;
            ho.reference = &g;
            ho.dipoles = false;
            ho.force_numerical = oo.numerical_hessian;
            const std::string ctx = format("Step %d · Hessian (%s)", step, oo.ts_hessian == "calcall" ? "CalcAll" : "CalcFC");
            report_progress_text(ctx);
            const HessianResult hr = run_hessian(x, bx, method, eo, ho, cb, &orbitals, [&](int i, int n, const std::string& what) {
                report_progress_text(ctx + format(": displacement %d of %d (%s)", i, n, what.c_str()));
            });
            H = internal_hessian(ic, x, hr.hessian, B, Gi, gq0);
            ++res.hessian_evaluations;
            rec.hessian_computed = true;
            rec.hessian_description = hr.description;
            rec.hessian_update = std::string("computed Hessian (") + (oo.ts_hessian == "calcall" ? "CalcAll" : "CalcFC") +
                                 "): " + hr.description + ", transformed to internal coordinates";
            if (res.initial_hessian.empty()) res.initial_hessian = std::string("computed (") +
                                                                    (oo.ts_hessian == "calcall" ? "CalcAll, every step" : "CalcFC") +
                                                                    "): " + hr.description;
            if (!q_prev.empty() && q_prev.size() == nq && std::fabs(predicted_prev) > 1.0e-14) {
                rec.ratio = rec.delta_e / predicted_prev;
                rec.has_ratio = true;
            }
            need_initial = false;
        } else if (need_initial || q_prev.size() != nq) {
            bool carried = false;
            if (first_build && oo.initial_hessian_matrix && oo.hessian_coordinates.size() == nq &&
                oo.initial_hessian_matrix->rows() == nq) {
                carried = true;
                for (std::size_t i = 0; i < nq; ++i)
                    if (!oo.hessian_coordinates[i].same_as(ic.coordinates()[i])) carried = false;
            }
            std::string kind = oo.initial_hessian;
            if (kind == "auto")
                kind = (method == "RHF" && bx.n_functions <= oo.analytic_hessian_max_functions) ? "analytic" : "lindh";
            // 10.0: "analytic" also for Full CI, CCSD and CCSD(T) (the analytic correlated Hessian);
            // "auto" keeps the model Hessian for them (an analytic correlated Hessian costs several gradients)
            if (kind == "analytic" && !analytic_hessian_method(method)) kind = "lindh";
            H = Matrix(nq, nq);
            if (carried) {
                H = *oo.initial_hessian_matrix;
                rec.hessian_update = "initial Hessian: " + oo.initial_hessian_label;
                res.initial_hessian = oo.initial_hessian_label;
            } else if (kind == "analytic") {
                Matrix Hc;
                std::string label = "analytic RHF Hessian";
                if (method == "RHF") {
                    Hc = rhf_hessian(x, bx, g.point.rhf, eo.threads, eo.memory_mb, 1e-8).hessian;
                } else {
                    HessianOptions ho;
                    ho.reference = &g;
                    ho.dipoles = false;
                    report_progress_text(format("Step %d · initial Hessian", step));
                    const HessianResult hr = run_hessian(x, bx, method, eo, ho, cb, &orbitals);
                    Hc = hr.hessian;
                    label = (hr.analytic ? "analytic " : "numerical ") + method + " Hessian";
                    ++res.hessian_evaluations;
                }
                const Matrix K = ic.b_derivative_contracted(x, gq0);
                for (std::size_t a = 0; a < n3; ++a)
                    for (std::size_t b = 0; b < n3; ++b) Hc(a, b) -= K(a, b);
                // H_q = G^- B (H_x - K) B' G^-
                const Matrix GB = mat_mul(Gi, B);
                const Matrix T = mat_mul(GB, Hc);
                for (std::size_t i = 0; i < nq; ++i)
                    for (std::size_t j = 0; j < nq; ++j) {
                        double s = 0.0;
                        for (std::size_t k = 0; k < n3; ++k) s += T(i, k) * GB(j, k);
                        H(i, j) = s;
                    }
                rec.hessian_update = "initial Hessian: " + label + " in internal coordinates";
                res.initial_hessian = label;
                // A minimisation needs a positive definite Hessian: negative eigenvalues (in the
                // non-redundant space) are replaced by their magnitude; the most negative direction
                // is kept, to leave a saddle point along it if the gradient there vanishes.
                const Matrix PHP = mat_mul(mat_mul(P, H), P);
                const EigenResult eh = symmetric_eigen(PHP);
                Matrix Hn(nq, nq);
                int negative = 0;
                for (std::size_t k = 0; k < nq; ++k) {
                    double lam = eh.values[k];
                    double proj = 0.0;
                    for (std::size_t i = 0; i < nq; ++i)
                        for (std::size_t j = 0; j < nq; ++j) proj += eh.vectors(i, k) * P(i, j) * eh.vectors(j, k);
                    if (proj < 0.5) continue;  // redundant direction
                    if (lam < 0.0) {
                        ++negative;
                        if (negative_direction.empty()) {
                            negative_direction.assign(nq, 0.0);
                            for (std::size_t i = 0; i < nq; ++i) negative_direction[i] = eh.vectors(i, k);
                            negative_curvature = lam;
                        }
                        lam = std::max(-lam, 0.005);
                    }
                    lam = std::max(lam, 0.005);
                    for (std::size_t i = 0; i < nq; ++i)
                        for (std::size_t j = 0; j < nq; ++j) Hn(i, j) += lam * eh.vectors(i, k) * eh.vectors(j, k);
                }
                H = Hn;
                if (negative)
                    rec.hessian_update += format("; %d negative eigenvalue%s (lowest %.4f) replaced by %s magnitude for the "
                                                 "minimisation",
                                                 negative, negative == 1 ? "" : "s", negative_curvature,
                                                 negative == 1 ? "its" : "their");
            } else {
                const std::vector<double> d = kind == "schlegel" ? schlegel_model_hessian(x, ic) : lindh_model_hessian(x, ic);
                for (std::size_t i = 0; i < nq; ++i) H(i, i) = d[i];
                rec.hessian_update = std::string("initial Hessian: ") +
                                     (kind == "schlegel" ? "Schlegel model Hessian" : "Lindh model Hessian");
                res.initial_hessian = kind == "schlegel" ? "Schlegel model Hessian" : "Lindh model Hessian";
            }
            need_initial = false;
        } else {
            const std::vector<double> s = ic.difference(q, q_prev);
            std::vector<double> y(nq);
            for (std::size_t i = 0; i < nq; ++i) y[i] = gq[i] - g_prev[i];
            double sy = 0.0;
            for (std::size_t i = 0; i < nq; ++i) sy += s[i] * y[i];
            std::vector<double> Hs(nq, 0.0);
            for (std::size_t i = 0; i < nq; ++i)
                for (std::size_t j = 0; j < nq; ++j) Hs[i] += H(i, j) * s[j];
            double sHs = 0.0;
            for (std::size_t i = 0; i < nq; ++i) sHs += s[i] * Hs[i];
            if (ts) {
                // Bofill: phi * Murtagh-Sargent (SR1) + (1 - phi) * Powell symmetric Broyden
                std::vector<double> xi(nq);
                for (std::size_t i = 0; i < nq; ++i) xi[i] = y[i] - Hs[i];
                double xs = 0.0, xx = 0.0, ss = 0.0;
                for (std::size_t i = 0; i < nq; ++i) {
                    xs += xi[i] * s[i];
                    xx += xi[i] * xi[i];
                    ss += s[i] * s[i];
                }
                if (ss > 1.0e-20 && xx > 1.0e-24) {
                    const double phi = xs * xs / (xx * ss);
                    const bool sr1 = std::fabs(xs) > 1.0e-12 * std::sqrt(xx * ss);
                    for (std::size_t i = 0; i < nq; ++i)
                        for (std::size_t j = 0; j < nq; ++j) {
                            const double ms = sr1 ? xi[i] * xi[j] / xs : 0.0;
                            const double psb = (xi[i] * s[j] + s[i] * xi[j]) / ss - xs * s[i] * s[j] / (ss * ss);
                            H(i, j) += (sr1 ? phi : 0.0) * ms + (1.0 - (sr1 ? phi : 0.0)) * psb;
                        }
                    rec.hessian_update = format("Bofill update (phi = %.4f: Murtagh-Sargent weight)", sr1 ? phi : 0.0);
                } else {
                    rec.hessian_update = "Bofill update skipped (no change of the gradient)";
                }
            } else if (sy > 1.0e-12 * norm(s) * norm(y) && sy > 0.0 && sHs > 1.0e-14) {
                for (std::size_t i = 0; i < nq; ++i)
                    for (std::size_t j = 0; j < nq; ++j) H(i, j) += y[i] * y[j] / sy - Hs[i] * Hs[j] / sHs;
                rec.hessian_update = format("BFGS update (s.y = %.3e)", sy);
            } else {
                rec.hessian_update = format("BFGS update skipped (curvature condition s.y = %.3e not positive)", sy);
            }
            // trust radius from the quality of the previous step (minimisations only; fixed for
            // transition states, as in Gaussian)
            if (std::fabs(predicted_prev) > 1.0e-14) {
                rec.ratio = rec.delta_e / predicted_prev;
                rec.has_ratio = true;
                if (!ts) {
                    if (rec.ratio < 0.25) trust = std::max(oo.min_trust, 0.25 * step_prev);
                    else if (rec.ratio > 0.75 && step_prev > 0.8 * trust) trust = std::min(oo.max_trust, 2.0 * trust);
                }
            }
            if (rec.delta_e > 1.0e-6 && !ts)
                rec.warnings.push_back(format("The energy rose by %.2e E_h; the trust radius is reduced.", rec.delta_e));
        }
        first_build = false;

        std::vector<double> dq(nq, 0.0);
        bool eigen_test_failed = false;
        if (ts) {
            // ---- partitioned RFO (Baker 1986) in the eigenvectors of the projected Hessian
            rec.ts_step = true;
            const Matrix PHP = mat_mul(mat_mul(P, H), P);
            const EigenResult eh = symmetric_eigen(PHP);
            std::vector<std::size_t> modes;  // non-redundant eigenvectors, ascending eigenvalues
            for (std::size_t k = 0; k < nq; ++k) {
                double proj = 0.0;
                for (std::size_t i = 0; i < nq; ++i)
                    for (std::size_t j = 0; j < nq; ++j) proj += eh.vectors(i, k) * P(i, j) * eh.vectors(j, k);
                if (proj >= 0.5) modes.push_back(k);
            }
            int negative = 0;
            for (std::size_t k : modes) {
                rec.hessian_eigenvalues.push_back(eh.values[k]);
                if (eh.values[k] < 0.0) ++negative;
            }
            rec.negative_eigenvalues = negative;
            // the followed mode: the lowest, then the one most like the previous followed mode
            std::size_t kf = modes.empty() ? 0 : modes[0];
            double best = -1.0;
            if (!followed_prev.empty() && followed_prev.size() == nq) {
                for (std::size_t k : modes) {
                    double ov = 0.0;
                    for (std::size_t i = 0; i < nq; ++i) ov += eh.vectors(i, k) * followed_prev[i];
                    if (std::fabs(ov) > best) {
                        best = std::fabs(ov);
                        kf = k;
                    }
                }
            }
            std::vector<double> vk(nq);
            double sign = 1.0;
            if (!followed_prev.empty() && followed_prev.size() == nq) {
                double ov = 0.0;
                for (std::size_t i = 0; i < nq; ++i) ov += eh.vectors(i, kf) * followed_prev[i];
                if (ov < 0.0) sign = -1.0;
            }
            for (std::size_t i = 0; i < nq; ++i) vk[i] = sign * eh.vectors(i, kf);
            rec.followed_mode = vk;
            rec.followed_overlap = best < 0.0 ? 1.0 : best;
            auto component = [&](std::size_t k) {
                double f = 0.0;
                for (std::size_t i = 0; i < nq; ++i) f += eh.vectors(i, k) * gq[i];
                return f;
            };
            const double bk = eh.values[kf], Fk = component(kf);
            const double lp = 0.5 * bk + 0.5 * std::sqrt(bk * bk + 4.0 * Fk * Fk);
            std::vector<std::size_t> others;
            for (std::size_t k : modes)
                if (k != kf) others.push_back(k);
            const std::size_t m = others.size();
            Matrix Aug(m + 1, m + 1);
            std::vector<double> Fo(m);
            for (std::size_t a = 0; a < m; ++a) {
                Fo[a] = component(others[a]);
                Aug(a, a) = eh.values[others[a]];
                Aug(a, m) = Aug(m, a) = Fo[a];
            }
            const double ln = symmetric_eigen(Aug).values[0];
            rec.followed_eigenvalue = bk;
            rec.followed_gradient = Fk;
            rec.lambda_p = lp;
            rec.lambda_n = ln;
            rec.rfo_lambda = ln;
            if (std::fabs(bk - lp) > 1.0e-14)
                for (std::size_t i = 0; i < nq; ++i) dq[i] += -Fk / (bk - lp) * eh.vectors(i, kf);
            for (std::size_t a = 0; a < m; ++a) {
                const double den = eh.values[others[a]] - ln;
                if (std::fabs(den) > 1.0e-14)
                    for (std::size_t i = 0; i < nq; ++i) dq[i] += -Fo[a] / den * eh.vectors(i, others[a]);
            }
            followed_prev = vk;
            if (oo.eigen_test && negative != 1) {
                eigen_test_failed = true;
                rec.warnings.push_back(format(
                    "Eigenvalue test: the Hessian has %d negative eigenvalue%s; a transition state (first-order saddle point) "
                    "needs exactly one.",
                    negative, negative == 1 ? "" : "s"));
            }
        } else {
        // ---- projected Hessian and RFO step (within the constrained space)
        Matrix Ht(nq, nq);
        {
            const Matrix PH = mat_mul(mat_mul(Pt, H), Pt);
            for (std::size_t i = 0; i < nq; ++i)
                for (std::size_t j = 0; j < nq; ++j) Ht(i, j) = PH(i, j) + 1000.0 * ((i == j ? 1.0 : 0.0) - Pt(i, j));
        }
        Matrix A(nq + 1, nq + 1);
        for (std::size_t i = 0; i < nq; ++i) {
            for (std::size_t j = 0; j < nq; ++j) A(i, j) = Ht(i, j);
            A(i, nq) = A(nq, i) = gt[i];
        }
        const EigenResult ea = symmetric_eigen(A);
        const double last = ea.vectors(nq, 0);
        rec.rfo_lambda = ea.values[0];
        if (std::fabs(last) > 1.0e-8) {
            for (std::size_t i = 0; i < nq; ++i) dq[i] = ea.vectors(i, 0) / last;
        } else {
            for (std::size_t i = 0; i < nq; ++i) dq[i] = -gt[i];
            rec.warnings.push_back("The RFO step was degenerate; a steepest-descent step is used.");
        }
        {
            std::vector<double> pdq(nq, 0.0);
            for (std::size_t i = 0; i < nq; ++i)
                for (std::size_t j = 0; j < nq; ++j) pdq[i] += Pt(i, j) * dq[j];
            dq = pdq;
        }
        if (step == 1 && !negative_direction.empty()) {
            std::vector<double> nd(nq, 0.0);
            for (std::size_t i = 0; i < nq; ++i)
                for (std::size_t j = 0; j < nq; ++j) nd[i] += Pt(i, j) * negative_direction[j];
            double gv = 0.0;
            for (std::size_t i = 0; i < nq; ++i) gv += gt[i] * nd[i];
            if (std::fabs(gv) < 1.0e-5 && norm(nd) > 0.5) {
                // Saddle point along the negative-curvature direction with no gradient there
                // (e.g. a linear H2O started at 180 degrees): leave it along that direction.
                for (std::size_t i = 0; i < nq; ++i) dq[i] += 0.1 * nd[i];
                rec.warnings.push_back(format(
                    "The starting geometry is a saddle point along a direction with negative curvature (%.4f E_h/unit^2) "
                    "and no gradient; the first step follows that direction.",
                    negative_curvature));
            }
        }
        }  // minimisation step
        if (ts) trust = oo.trust_radius;  // fixed for transition states
        double sn = norm(dq);
        if (sn > trust) {
            for (double& d : dq) d *= trust / sn;
            rec.step_scaled = true;
            sn = trust;
        }
        // drive the frozen coordinates to their targets
        for (std::size_t i = 0; i < nq; ++i) dq[i] += cp.dq_c[i];
        if (constrained) sn = norm(dq);
        rec.trust_radius = trust;
        rec.step_norm = sn;
        double pred = 0.0;
        {
            const Matrix PHP = mat_mul(mat_mul(P, H), P);
            for (std::size_t i = 0; i < nq; ++i) {
                pred += gq[i] * dq[i];
                for (std::size_t j = 0; j < nq; ++j) pred += 0.5 * dq[i] * PHP(i, j) * dq[j];
            }
        }
        rec.predicted_change = pred;
        rec.dq = dq;

        // ---- convergence
        double mf = 0.0, rf = 0.0, md = 0.0, rd = 0.0;
        for (std::size_t i = 0; i < nq; ++i) {
            mf = std::max(mf, std::fabs(gt[i]));
            rf += gt[i] * gt[i];
            md = std::max(md, std::fabs(dq[i]));
            rd += dq[i] * dq[i];
        }
        const double nd = static_cast<double>(std::max<std::size_t>(1, nq));
        rf = std::sqrt(rf / nd);
        rd = std::sqrt(rd / nd);
        rec.convergence = {{"Maximum force", mf, cc.max_force, mf < cc.max_force, true},
                           {"RMS force", rf, cc.rms_force, rf < cc.rms_force, true},
                           {"Maximum displacement", md, cc.max_displacement, md < cc.max_displacement, true},
                           {"RMS displacement", rd, cc.rms_displacement, rd < cc.rms_displacement, true},
                           {"Energy change", std::fabs(rec.delta_e), cc.energy,
                            step > 1 && std::fabs(rec.delta_e) < cc.energy, step > 1}};
        if (constrained)
            rec.convergence.push_back({"Constraint deviation", cp.deviation, constraint_threshold,
                                       cp.deviation < constraint_threshold, true});
        const bool four = rec.convergence[0].met && rec.convergence[1].met && rec.convergence[2].met &&
                          rec.convergence[3].met;
        const bool energy_ok = step == 1 || rec.convergence[4].met;
        const bool constraint_ok = !constrained || cp.deviation < constraint_threshold;
        const bool tiny_forces = mf < 0.01 * cc.max_force && rf < 0.01 * cc.rms_force;
        rec.converged = ((four && energy_ok) || tiny_forces) && constraint_ok && !eigen_test_failed;
        report_progress(JobStage::Optimization, step, rec.energy, mf);
        if (ts) res.final_negative_eigenvalues = rec.negative_eigenvalues;
        if (eigen_test_failed) {
            rec.dq.assign(nq, 0.0);
            rec.seconds = seconds_since(ts_start);
            res.history.push_back(rec);
            if (on_step) on_step(res.history.back());
            res.eigen_test_failed = true;
            res.steps = step;
            res.final_point = std::move(g);
            res.final_geometry = x;
            res.final_energy = rec.energy;
            res.final_coordinates = ic.coordinates();
            res.final_hessian = H;
            res.message = format(
                "The transition-state search stopped at step %d: the eigenvalue test failed (%d negative eigenvalues of "
                "the Hessian, one needed). Start nearer to the transition state (e.g. from the highest point of a relaxed "
                "scan), compute the Hessian at every step (CalcAll), or switch the test off (Gaussian: NoEigenTest).",
                step, rec.negative_eigenvalues);
            break;
        }

        if (rec.converged) {
            rec.dq.assign(nq, 0.0);
            rec.seconds = seconds_since(ts_start);
            res.history.push_back(rec);
            if (on_step) on_step(res.history.back());
            res.converged = true;
            res.steps = step;
            res.final_point = std::move(g);
            res.final_geometry = x;
            res.final_energy = rec.energy;
            res.final_constraint_deviation = cp.deviation;
            res.message = constrained ? format("Constrained stationary point found after %d step%s.", step, step == 1 ? "" : "s")
                          : ts ? format("Transition state found after %d step%s (%d negative eigenvalue of the Hessian, %.6f).",
                                        step, step == 1 ? "" : "s", rec.negative_eigenvalues, rec.followed_eigenvalue)
                               : format("Stationary point found after %d step%s.", step, step == 1 ? "" : "s");
            res.final_coordinates = ic.coordinates();
            res.final_hessian = H;
            break;
        }
        // ---- take the step
        const BackTransform bt = back_transform(ic, x, dq);
        rec.backtransform_iterations = bt.iterations;
        rec.backtransform_residual = bt.residual;
        if (!bt.converged)
            rec.warnings.push_back(format("The back-transformation to Cartesian coordinates stopped after %d iterations "
                                          "(residual %.1e).",
                                          bt.iterations, bt.residual));
        rec.seconds = seconds_since(ts_start);
        res.history.push_back(rec);
        if (on_step) on_step(res.history.back());
        q_prev = q;
        g_prev = gq;
        e_prev = rec.energy;
        predicted_prev = pred;
        step_prev = sn;
        res.final_point = std::move(g);
        res.final_geometry = x;
        res.final_energy = rec.energy;
        res.final_constraint_deviation = cp.deviation;
        res.final_coordinates = ic.coordinates();
        res.final_hessian = H;
        res.steps = step;
        if (step == oo.max_steps) {
            res.message = format("The optimisation did not converge within %d steps.", oo.max_steps);
            break;
        }
        x = bt.molecule;
    }
    res.final_parameters = geometric_parameters(res.final_geometry);
    res.initial_energy = res.history.empty() ? 0.0 : res.history.front().energy;
    res.symmetry_deviation = c2v_deviation(res.final_geometry);
    res.symmetric_end = res.final_geometry.atoms.size() == 3 && !res.final_geometry.custom && res.symmetry_deviation < 1.0e-6;
    res.seconds = seconds_since(t0);
    return res;
}

}  // namespace macqc
