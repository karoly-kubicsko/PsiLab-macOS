#include "macqc/molecular.h"

#include "macqc/gradient.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <stdexcept>

#include "macqc/cc.h"
#include "macqc/gto_integrals.h"
#include "macqc/mo_integrals.h"

namespace macqc {

namespace {

double distance_bohr(const Atom& a, const Atom& b) {
    double s = 0.0;
    for (std::size_t k = 0; k < 3; ++k) s += (a.position[k] - b.position[k]) * (a.position[k] - b.position[k]);
    return std::sqrt(s);
}

std::string format(const char* fmt, ...) {
    char buf[1024];
    va_list args;
    va_start(args, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, args);
    va_end(args);
    return buf;
}

inline int nc(int l) { return (l + 1) * (l + 2) / 2; }

// Work units of one shell quartet class (la >= lb, lc >= ld as in the engine).
struct ClassCost {
    double per_primitive = 0.0;  // vertical recurrence + accumulation, per primitive quartet
    double per_quartet = 0.0;    // horizontal recurrences + spherical transformation
};

ClassCost class_cost(int la, int lb, int lc, int ld) {
    const int Lab = la + lb, Lcd = lc + ld, L = Lab + Lcd;
    double vrr = 0.0;
    for (int E = 1; E <= Lab; ++E) vrr += nc(E) * (L - E + 1);
    for (int F = 1; F <= Lcd; ++F) {
        const int emin = std::max(0, la - (Lcd - F));
        double se = 0.0;
        for (int E = emin; E <= Lab; ++E) se += nc(E);
        vrr += (Lcd - F + 1) * nc(F) * se;
    }
    double nk = 0.0;
    for (int F = lc; F <= Lcd; ++F) nk += nc(F);
    double acc = 0.0;
    for (int E = la; E <= Lab; ++E) acc += nc(E) * nk;
    double post = 0.0;
    for (int j = 1; j <= lb; ++j)
        for (int e = la; e <= la + lb - j; ++e) post += nc(e) * nc(j) * nk;
    post += 2.0 * nc(la) * nc(lb) * nk;
    const double nbra = (2 * la + 1) * (2 * lb + 1);
    for (int j = 1; j <= ld; ++j)
        for (int f = lc; f <= lc + ld - j; ++f) post += nc(f) * nc(j) * nbra;
    post += 2.0 * nc(lc) * nc(ld) * nbra;
    ClassCost c;
    c.per_primitive = L == 0 ? 8.0 : vrr + acc + 40.0;
    c.per_quartet = post + 150.0;
    return c;
}

// Seconds per work unit and per stored integral in a Fock build (one thread), calibrated
// with H2O and LiH in aug-cc-pV5Z on a 2.1 GHz x86-64 core (Apple silicon cores are faster).
constexpr double kSecondsPerUnit = 1.14e-9;
constexpr double kSecondsPerStoredIntegral = 3.5e-9;

std::string angstrom_bohr(double bohr) {
    return format("%.6f \xC3\x85 (%.6f bohr)", bohr * kBohrInAngstrom, bohr);
}

std::string human_seconds(double s) {
    if (s < 1.0) return "less than a second";
    if (s < 90.0) return format("about %.0f s", s);
    if (s < 5400.0) return format("about %.0f min", s / 60.0);
    return format("about %.1f h", s / 3600.0);
}

std::string thousands(std::size_t v) {
    std::string s = std::to_string(v), out;
    int k = 0;
    for (auto it = s.rbegin(); it != s.rend(); ++it) {
        if (k && k % 3 == 0) out += ',';
        out += *it;
        ++k;
    }
    return std::string(out.rbegin(), out.rend());
}

std::string thousands_d(double v) {
    if (v >= 1.0e15) return format("%.2e", v);
    return thousands(static_cast<std::size_t>(std::llround(v)));
}

std::string orbital_list(const std::vector<int>& o) {
    std::string s;
    for (std::size_t k = 0; k < o.size(); ++k) s += (k ? " " : "") + std::to_string(o[k]);
    return s;
}

std::string fci_section_impl(const FCIResultGeneral& f, double rhf_energy, int level) {
    const int mult = f.alpha - f.beta + 1;
    static const char* spin_names[] = {"singlet", "doublet", "triplet", "quartet", "quintet", "sextet", "septet"};
    const std::string spin_name = mult >= 1 && mult <= 7 ? spin_names[mult - 1] : format("multiplicity %d", mult);
    std::string s = format("\nFull CI         all electrons correlated, M_S = %s determinants, lowest %s state\n",
                           f.alpha == f.beta ? "0" : format("%g", 0.5 * (f.alpha - f.beta)).c_str(), spin_name.c_str());
    s += format("    %d electrons (%d alpha, %d beta) in %d orbitals: %s alpha x %s beta strings = %s determinants\n",
                f.alpha + f.beta, f.alpha, f.beta, f.orbitals, thousands(f.alpha_strings).c_str(),
                thousands(f.beta_strings).c_str(), thousands(f.determinants).c_str());
    if (f.spin_adapted)
        s += format("    spin-adapted: %s configuration state functions with S = %g (Weyl-Paldus), the Davidson\n"
                    "    vectors projected onto S = %g (Loewdin projector from S^2), sigma = H c in the determinants\n",
                    thousands(static_cast<std::size_t>(std::llround(f.csfs))).c_str(), f.spin, f.spin);
    s += format("    integral transformation to the %s orbitals: %.2f s\n", f.alpha == f.beta ? "RHF" : "ROHF",
                f.transform_seconds);
    s += f.alpha == f.beta ? "    Davidson method: diagonal preconditioner, singlet projection (C = C^T),\n"
                           : "    Davidson method: diagonal preconditioner,\n";
    s += "    converged when |H c - E c| < 1.0e-06 and |dE| < 1.0e-10 E_h\n";
    if (level >= 2) {
        s += "    iter        total energy (E_h)        residual   subspace\n";
        for (const auto& it : f.history)
            s += format("    %4d  %24.12f  %12.3e  %6d\n", it.iteration, it.energy, it.residual, it.subspace);
    }
    s += f.converged ? format("    Davidson converged in %d iterations (%.2f s)\n", f.iterations, f.ci_seconds)
                     : format("    Davidson NOT converged after %d iterations\n", f.iterations);
    s += "\nFull CI energies (E_h)\n";
    s += format("    RHF energy                %20.10f\n", rhf_energy);
    s += format("    Full CI energy            %20.10f\n", f.total_energy);
    s += format("    correlation energy        %20.10f\n", f.correlation_energy);
    s += format("    <S^2>                     %20.10f\n", f.s2);
    s += format("    weight of the RHF determinant |C_0|^2  %.8f\n", f.hf_weight);
    s += "\nLeading determinants (occupied MOs, alpha | beta; excitation level relative to RHF)\n";
    for (const auto& d : f.leading)
        s += format("    %+12.8f   %-14s | %-14s  %d\n", d.coefficient, orbital_list(d.alpha).c_str(),
                    orbital_list(d.beta).c_str(), d.excitation);
    const std::size_t nshow = level >= 2 ? std::min<std::size_t>(f.natural_occupations.size(), 20)
                                         : std::min<std::size_t>(f.natural_occupations.size(), 8);
    s += format("\nNatural orbital occupation numbers (largest %zu; they sum to %d)\n", nshow, f.alpha + f.beta);
    for (std::size_t k = 0; k < nshow; ++k)
        s += format("  %12.8f", f.natural_occupations[k]) + ((k % 5 == 4 || k + 1 == nshow) ? "\n" : "");
    return s;
}

std::string cc_section_impl(const CCResult& c, int level) {
    std::string s = "\nCoupled cluster " + c.method + ", all electrons correlated\n";
    s += format("    %d occupied and %d virtual orbitals (%d electrons); %s implementation\n", c.occupied, c.virtuals,
                2 * c.occupied, c.engine.c_str());
    std::string counts;
    const char* names[] = {"", "T1", "T2", "T3", "T4"};
    for (std::size_t k = 1; k < c.amplitude_counts.size() && k < 5; ++k)
        counts += format("%s%s %s", counts.empty() ? "" : ", ", names[k], thousands_d(c.amplitude_counts[k]).c_str());
    s += "    amplitudes: " + counts +
         (c.engine.find("spin-orbital") != std::string::npos ? " (distinct spin-orbital sets)\n" : " (spatial)\n");
    if (!c.ladder.empty()) s += "    (vv|vv) ladder: " + c.ladder + "\n";
    s += format("    integral transformation: %.2f s; Jacobi updates with DIIS, converged when\n", c.transform_seconds);
    s += format("    max |t_new - t| < %.1e and |dE| < %.1e E_h\n", c.amplitude_tol, c.energy_tol);
    if (level >= 2) {
        s += "    iter        total energy (E_h)     correlation (E_h)        dE (E_h)   max |dt|\n";
        for (const auto& it : c.history)
            s += format("    %4d  %24.12f  %20.12f  %14.3e  %9.2e\n", it.iteration, it.energy, it.correlation, it.delta_e,
                        it.residual);
    }
    s += c.converged ? format("    converged in %d iterations (%.2f s)\n", c.iterations, c.cc_seconds)
                     : format("    NOT converged after %d iterations\n", c.iterations);
    if (c.has_triples) s += format("    (T) correction: %.2f s\n", c.triples_seconds);
    s += "\n" + c.method + " energies (E_h)\n";
    s += format("    RHF energy                %20.10f\n", c.reference_energy);
    const std::string base = c.has_triples ? std::string("CCSD") : c.method;
    s += format("    %-26s%20.10f\n", (base + " correlation energy").c_str(), c.correlation_energy);
    s += format("    %-26s%20.10f\n", (base + " energy").c_str(), c.cc_energy);
    if (c.has_triples) {
        s += format("    (T) correction            %20.10f\n", c.triples_correction);
        s += format("    CCSD(T) correlation energy%20.10f\n", c.correlation_energy + c.triples_correction);
        s += format("    CCSD(T) energy            %20.10f\n", c.total_energy);
    }
    s += format("    T1 diagnostic             %20.6f\n", c.t1_diagnostic);
    std::string norms;
    for (std::size_t k = 1; k < c.rank_norms.size() && k < 5; ++k)
        norms += format("%s|%s| = %.6f", norms.empty() ? "" : ", ", names[k], c.rank_norms[k]);
    if (!norms.empty()) s += "    amplitude norms: " + norms + "\n";
    s += "\nLargest amplitudes (occupied -> virtual MOs, spin a = alpha, b = beta)\n";
    for (const auto& a : c.leading)
        s += format("    T%d  %+12.8f   %-12s -> %-12s %s\n", a.rank, a.value, orbital_list(a.occupied).c_str(),
                    orbital_list(a.virtuals).c_str(), a.spins.c_str());
    return s;
}

}  // namespace

std::string fci_report_section(const FCIResultGeneral& f, double rhf_energy, int level) {
    return fci_section_impl(f, rhf_energy, level);
}

std::string cc_report_section(const CCResult& c, int level) { return cc_section_impl(c, level); }

std::string so_report_section(const SOResult& r, int level) {
    std::string s;
    const std::string label = method_label(r.method, r.reference);
    const bool mp = r.method.rfind("MP", 0) == 0;
    const bool ci = r.method.rfind("CI", 0) == 0;
    const bool qci = r.method.rfind("QCI", 0) == 0;
    s += "\n" + label + (mp ? " (Moller-Plesset perturbation theory, H0 = Fock operator, all electrons correlated)"
                          : ci ? " (configuration interaction from the reference, intermediate normalisation, all "
                                 "electrons correlated)"
                          : qci ? " (quadratic configuration interaction, all electrons correlated)"
                                : " (coupled cluster, all electrons correlated)") +
         "\n";
    s += format("  %s on the %s reference; spin-orbital engine (dense spin blocks, verified equations)\n",
                r.converged ? "converged" : "NOT converged", r.reference.c_str());
    s += format("  reference energy (%s)        %20.10f E_h\n", r.reference.c_str(), r.reference_energy);
    double running = r.reference_energy;
    for (const auto& t : r.terms) {
        if (mp) {
            running += t.value;
            s += format("  %-28s %20.10f E_h   MP%c total %20.10f E_h\n", (t.label + " correction").c_str(), t.value,
                        t.label[2], running);
        } else {
            s += format("  %-28s %20.10f E_h\n", t.label.c_str(), t.value);
        }
    }
    s += format("  correlation energy           %20.10f E_h\n", r.correlation_energy);
    s += format("  total energy %-15s %20.10f E_h\n", label.c_str(), r.total_energy);
    if (level >= 1) {
        s += "  amplitudes (unique, in spin orbitals):";
        for (const auto& a : r.amplitudes) s += format(" %s %s (|t| %.3e);", a.name.c_str(), thousands_d(a.count).c_str(), a.norm);
        s += "\n";
        if (!mp) s += format("  T1 diagnostic %.6f\n", r.t1_diagnostic);
        s += format("  %d iterations, %.2f s\n", r.iterations, r.seconds);
    }
    return s;
}

namespace {
std::string compact_upper(const std::string& method) {
    std::string m;
    for (char c : method)
        if (c != '-' && c != ' ' && c != '_') m += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return m;
}

std::string correlated_name(const std::string& m) {
    static const char* names[] = {"CCSD", "CCSD(T)", "CCSDT", "CCSDTQ", "FCI", "MP2", "MP3", "MP4", "MP5", "CISD",
                                  "CISD(T)", "CISDT", "CISDTQ", "QCISD", "QCISD(T)", "QCISDT", "QCISDTQ"};
    for (const char* n : names)
        if (m == n) return n;
    if (m == "CCSD[T]") return "CCSD(T)";
    if (m == "QCISD[T]") return "QCISD(T)";
    if (m == "FULLCI") return "FCI";
    if (m == "MP4(SDTQ)" || m == "MP4SDTQ") return "MP4";
    return "";
}
}  // namespace

MethodSpec parse_method_spec(const std::string& text) {
    const std::string m = compact_upper(text);
    if (m == "HF" || m == "SCF") return {"RHF", ""};
    if (m == "RHF" || m == "UHF" || m == "ROHF") return {"RHF", m};
    std::string c = correlated_name(m);
    if (!c.empty()) return {c, ""};
    if (m.rfind("RO", 0) == 0 && !(c = correlated_name(m.substr(2))).empty()) return {c, "ROHF"};
    if (m.rfind("U", 0) == 0 && !(c = correlated_name(m.substr(1))).empty()) return {c, "UHF"};
    if (m.rfind("R", 0) == 0 && !(c = correlated_name(m.substr(1))).empty()) return {c, "RHF"};
    return {"", ""};
}

std::string effective_reference(const std::string& method, const std::string& requested, int multiplicity) {
    const std::string r = requested.empty() ? "" : canonical_reference_name(requested);
    if (!requested.empty() && r.empty()) throw std::runtime_error("reference must be RHF, UHF or ROHF");
    if (r == "RHF" && multiplicity != 1)
        throw std::runtime_error(format("restricted Hartree-Fock (RHF) needs a closed-shell singlet; multiplicity %d "
                                        "needs UHF or ROHF",
                                        multiplicity));
    if (!r.empty()) return r;
    if (multiplicity == 1) return "RHF";
    return method == "FCI" ? "ROHF" : "UHF";
}

std::string method_label(const std::string& method, const std::string& reference) {
    if (method == "RHF") return reference.empty() ? "RHF" : reference;
    if (reference == "UHF") return "U" + method;
    if (reference == "ROHF") return "RO" + method;
    return method;
}

std::string canonical_method_name(const std::string& method) { return parse_method_spec(method).method; }

bool coupled_cluster_method(const std::string& canonical, CCMethod* method) {
    CCMethod mm;
    if (canonical == "CCSD") mm = CCMethod::CCSD;
    else if (canonical == "CCSD(T)") mm = CCMethod::CCSD_T;
    else if (canonical == "CCSDT") mm = CCMethod::CCSDT;
    else if (canonical == "CCSDTQ") mm = CCMethod::CCSDTQ;
    else return false;
    if (method) *method = mm;
    return true;
}

std::string validate_molecular(const MolecularInput& in) {
    const std::string m = canonical_molecule_name(in.molecule);
    if (m.empty() && !find_molecule_preset(in.molecule))
        return "molecule must be HeH+, H2, LiH, H2O, H2O2, NH3, an atom H-Ar, H2+, H3+, H5+, N2, O2, F2, Cl2, NO, NO2, CO, CO2, N2O, H2CO, CH4 "
               "or custom";
    if (in.multiplicity < 1 || in.multiplicity > 20) return "the multiplicity must be between 1 and 20";
    if (m == "custom") {
        if (in.atoms.empty()) return "a custom molecule needs atoms (read from an .xyz or .gjf file)";
        if (in.atoms.size() > 200) return "a custom molecule may have at most 200 atoms";
        if (in.charge < -20 || in.charge > 20) return "the charge must be between -20 and +20";
        int nel = -in.charge;
        for (const Atom& a : in.atoms) {
            if (a.Z < 1 || a.Z > kMaxAtomicNumber) return "unknown element in the custom molecule (Z = 1 to 118)";
            nel += a.Z;
        }
        if (nel < 1) return format("the molecule has %d electrons with charge %+d; it needs at least one", nel, in.charge);
        if (in.multiplicity < 1 || in.multiplicity > 20) return "the multiplicity must be between 1 and 20";
    }
    {
        // 15.0: the electrons must fit the multiplicity
        int nuclear = 0, charge = in.charge_given ? in.charge : 0;
        if (m == "custom") {
            for (const Atom& a : in.atoms) nuclear += a.Z;
            charge = in.charge;
        } else if (const MoleculePreset* pre = m.empty() ? find_molecule_preset(in.molecule) : nullptr) {
            for (const auto& a : pre->atoms) nuclear += atomic_number(a.first);
            if (!in.charge_given) charge = pre->charge;
        } else {
            nuclear = m == "H2" ? 2 : m == "LiH" ? 4 : (m == "H2O" || m == "NH3") ? 10 : m == "H2O2" ? 18 : 3;
            if (m == "HeH+" && !in.charge_given) charge = 1;
        }
        const int mult = in.multiplicity_given || m == "custom" ? in.multiplicity
                         : (m.empty() && find_molecule_preset(in.molecule)) ? find_molecule_preset(in.molecule)->multiplicity
                                                                             : 1;
        try {
            alpha_beta_electrons(nuclear - charge, mult);
        } catch (const std::exception& e) {
            return e.what();
        }
    }
    if (canonical_basis_name(in.basis).empty()) {
        std::string list;
        for (const auto& b : supported_basis_sets()) list += (list.empty() ? "" : ", ") + b;
        return "unknown basis set '" + in.basis + "' (available: " + list + ")";
    }
    if (canonical_method_name(in.method).empty())
        return "method must be HF (RHF, UHF, ROHF), MP2, MP3, MP4, MP5, CISD, CISD(T), CISDT, CISDTQ, QCISD, QCISD(T), "
               "QCISDT, QCISDTQ, CCSD, CCSD(T), CCSDT, CCSDTQ or FCI (with U or RO for the reference, e.g. UMP2)";
    if (!in.reference.empty() && canonical_reference_name(in.reference).empty())
        return "reference must be RHF, UHF or ROHF";
    if (!std::isfinite(in.fci_memory_mb) || in.fci_memory_mb < 64.0)
        return "fciMemoryMB (memory limit of Full CI and coupled cluster) must be at least 64";
    if (!std::isfinite(in.R) || in.R < 0.0) return "R must be a positive, finite number (bohr)";
    if (in.R > 0.0 && in.R < 0.3) return "R must be at least 0.3 bohr (0.16 \xC3\x85)";
    if (in.R > 50.0) return "R must be at most 50 bohr (26.5 \xC3\x85)";
    if (m == "H2O" && !(std::isfinite(in.angle) && in.angle >= 30.0 && in.angle <= 180.0))
        return "the H-O-H angle must be between 30 and 180 degrees";
    if (m == "H2O2") {
        if (!(std::isfinite(in.angle) && in.angle >= 30.0 && in.angle <= 175.0))
            return "the O-O-H angle must be between 30 and 175 degrees";
        if (!std::isfinite(in.r_oo) || in.r_oo < 0.0 || (in.r_oo > 0.0 && in.r_oo < 0.3) || in.r_oo > 50.0)
            return "R(O-O) must be between 0.3 and 50 bohr";
        if (!(std::isfinite(in.dihedral) && in.dihedral >= -360.0 && in.dihedral <= 360.0))
            return "the H-O-O-H dihedral angle must be between -360 and 360 degrees";
    }
    if (m == "NH3" && !(std::isfinite(in.angle) && in.angle >= 60.0 && in.angle <= 120.0))
        return "the H-N-H angle must be between 60 and 120 degrees (120: planar)";
    if (in.threads < 0 || in.threads > 64) return "threads must be between 0 (automatic) and 64";
    if (!std::isfinite(in.memory_mb) || in.memory_mb < 0.0) return "memoryMB must be >= 0";
    if (in.max_iter < 1 || in.max_iter > 1000) return "scfMaxIter must be between 1 and 1000";
    if (!(in.energy_tol > 0.0) || !(in.gradient_tol > 0.0)) return "SCF convergence thresholds must be positive";
    if (in.print_level < 0 || in.print_level > 2) return "printLevel must be 0, 1 or 2";
    return "";
}

int expected_scf_iterations(const std::string& molecule) {
    // Observed with the core-Hamiltonian guess and DIIS near equilibrium (7, 10, 13), plus one.
    if (molecule == "H2" || molecule == "HeH+") return 8;
    if (molecule == "LiH") return 11;
    if (molecule == "H2O2") return 16;
    if (molecule == "custom") return 18;
    return 14;
}

double estimate_rhf_seconds(const BasisSet& basis, bool eri_stored, int threads, int iterations) {
    const ShellPairList pairs = make_shell_pairs(basis);
    // Group the pairs by (la, lb, number of primitive pairs) to keep the double loop cheap.
    struct Group {
        int la, lb;
        double count = 0.0, prims = 0.0;
    };
    std::vector<Group> groups;
    for (const ShellPair& p : pairs.pairs) {
        const int la = basis.shells[p.a].l, lb = basis.shells[p.b].l;
        const double k = static_cast<double>(p.prims.size());
        if (k == 0.0) continue;
        auto it = std::find_if(groups.begin(), groups.end(),
                               [&](const Group& g) { return g.la == la && g.lb == lb; });
        if (it == groups.end()) {
            groups.push_back({la, lb, 0.0, 0.0});
            it = groups.end() - 1;
        }
        it->count += 1.0;
        it->prims += k;
    }
    // Sum over unordered pairs of pairs: ~ half of all ordered combinations.
    double units = 0.0;
    for (const Group& g : groups)
        for (const Group& h : groups) {
            const ClassCost c = class_cost(g.la, g.lb, h.la, h.lb);
            units += 0.5 * (g.prims * h.prims * c.per_primitive + g.count * h.count * c.per_quartet);
        }
    const double t = std::max(1, threads);
    const double one_pass = units * kSecondsPerUnit / t;
    if (eri_stored) {
        const double n = static_cast<double>(unique_eri_count(static_cast<std::size_t>(basis.n_functions)));
        return one_pass + iterations * n * kSecondsPerStoredIntegral / t;
    }
    return iterations * one_pass;
}

MolecularResult prepare_molecular(const MolecularInput& input) {
    const std::string problem = validate_molecular(input);
    if (!problem.empty()) throw std::runtime_error(problem);
    MolecularResult r;
    r.input = input;
    // 15.0: presets of the molecule menu run as custom molecules
    if (canonical_molecule_name(input.molecule).empty()) {
        const MoleculePreset* pre = find_molecule_preset(input.molecule);
        r.input.preset = pre->key;
        if (r.input.atoms.empty())
            for (const auto& a : pre->atoms) {
                Atom at;
                at.Z = atomic_number(a.first);
                for (int k = 0; k < 3; ++k) at.position[static_cast<std::size_t>(k)] = a.second[static_cast<std::size_t>(k)] / kBohrInAngstrom;
                r.input.atoms.push_back(at);
            }
        if (!input.charge_given) r.input.charge = pre->charge;
        if (!input.multiplicity_given) r.input.multiplicity = pre->multiplicity;
        if (r.input.title.empty()) r.input.title = pre->title;
        r.input.molecule = "custom";
    } else {
        r.input.molecule = canonical_molecule_name(input.molecule);
    }
    r.input.basis = canonical_basis_name(input.basis);
    {
        const MethodSpec ms = parse_method_spec(input.method);
        r.input.method = ms.method;
        if (r.input.reference.empty()) r.input.reference = ms.reference;
        else r.input.reference = canonical_reference_name(r.input.reference);
    }
    if (r.input.molecule == "custom") {
        try {
            r.molecule = make_custom_molecule(r.input.atoms, r.input.charge, r.input.title);
        } catch (const std::invalid_argument& e) {
            throw std::runtime_error(e.what());
        }
        r.input.atoms = r.molecule.atoms;
        r.input.R = r.molecule.atoms.size() >= 2 ? distance_bohr(r.molecule.atoms[0], r.molecule.atoms[1]) : 0.0;
    } else {
        if (!(r.input.R > 0.0)) r.input.R = default_bond_length(r.input.molecule);
        if (r.input.molecule == "H2O2" && !(r.input.r_oo > 0.0)) r.input.r_oo = default_oo_bond_length();
        r.molecule = make_molecule(r.input.molecule, r.input.R, r.input.angle, r.input.r_oo, r.input.dihedral);
        if (input.charge_given) r.molecule.charge = input.charge;  // 15.0: e.g. H2O+
        r.input.charge = r.molecule.charge;
    }
    // 15.0: electrons and spin
    try {
        alpha_beta_electrons(r.molecule.electrons(), r.input.multiplicity);
        r.reference = effective_reference(r.input.method, r.input.reference, r.input.multiplicity);
    } catch (const std::runtime_error& e) {
        throw std::runtime_error(e.what());
    }
    r.basis = build_basis(r.molecule, r.input.basis, r.input.basis_directory);
    r.molecule = attach_ecps(r.molecule, r.basis);  // 14.0: effective core potentials
    for (const auto& a : r.basis.atoms)
        if (!a.note.empty() && std::find(r.notes.begin(), r.notes.end(), a.note) == r.notes.end())
            r.notes.push_back(a.note);

    CostEstimate& e = r.estimate;
    e.basis_functions = r.basis.n_functions;
    e.shells = static_cast<int>(r.basis.shells.size());
    e.max_l = r.basis.max_l;
    e.unique_eris = unique_eri_count(static_cast<std::size_t>(r.basis.n_functions));
    e.eri_memory_mb = static_cast<double>(e.unique_eris) * 8.0 / (1024.0 * 1024.0);
    e.eri_stored = r.input.memory_mb > 0.0 && e.eri_memory_mb <= r.input.memory_mb;
    e.threads = resolve_threads(r.input.threads);
    e.rhf_seconds =
        estimate_rhf_seconds(r.basis, e.eri_stored, e.threads, expected_scf_iterations(r.molecule.name));
    e.total_seconds = e.rhf_seconds;
    const auto ab = alpha_beta_electrons(r.molecule.electrons(), r.input.multiplicity);
    if (uses_so_engine(r.input.method, r.reference)) {
        const int nbf = r.basis.n_functions;
        const int on[2] = {ab.first, ab.second}, vn[2] = {nbf - ab.first, nbf - ab.second};
        const SOSize sz = estimate_so_method(r.input.method, on, vn);
        CostEstimate::SOCost c;
        c.method = method_label(r.input.method, r.reference);
        c.amplitudes = sz.amplitudes;
        c.memory_bytes = sz.memory_bytes;
        // dense contractions at about 1.5 GFlop/s per thread (portable GEMM, measured), 20 iterations
        const double iterations = r.input.method.rfind("MP", 0) == 0 ? 2.0 * (r.input.method[2] - '1') : 20.0;
        c.seconds = iterations * sz.flops_per_iteration / (1.5e9 * e.threads) +
                    2.0 * estimate_rhf_seconds(r.basis, false, e.threads, 1);
        if (vn[0] <= 0 && vn[1] <= 0) {
            c.feasible = false;
            c.reason = "there are no virtual orbitals to correlate in this basis set";
        } else if (c.memory_bytes / 1048576.0 > r.input.fci_memory_mb) {
            c.feasible = false;
            c.reason = format("%s would need about %.1f GB of memory (limit %.1f GB): %s amplitudes in dense spin blocks.",
                              c.method.c_str(), c.memory_bytes / 1073741824.0, r.input.fci_memory_mb / 1024.0,
                              thousands_d(c.amplitudes).c_str());
        } else if (c.seconds > kMaxCCSeconds) {
            c.feasible = false;
            c.reason = format("%s would take about %s; ΨLab refuses calculations above three days.", c.method.c_str(),
                              human_seconds(c.seconds).c_str());
        }
        e.so = c;
        e.total_seconds += c.seconds;
        if (!c.feasible) {
            r.refusal = c.reason + " Use a smaller basis set or a lower-order method.";
            return r;
        }
        if (e.total_seconds > kLongCalculationSeconds)
            r.warnings.push_back(format("Long calculation: %s with %s amplitudes after %s; about %.1f GB of memory. "
                                        "Estimated time %s on %d threads.",
                                        c.method.c_str(), thousands_d(c.amplitudes).c_str(), r.reference.c_str(),
                                        c.memory_bytes / 1073741824.0, human_seconds(e.total_seconds).c_str(), e.threads));
        return r;
    }
    if (r.input.method == "FCI") {
        const int nel = r.molecule.electrons();
        const double pass = estimate_rhf_seconds(r.basis, false, e.threads, 1);
        e.fci = estimate_fci(r.basis.n_functions, r.basis.n_functions, ab.first, ab.second, e.threads,
                             r.input.fci_memory_mb, pass);
        (void)nel;
        e.total_seconds += e.fci->seconds;
        if (!e.fci->feasible) {
            r.refusal = e.fci->reason;
            if (r.molecule.custom)
                r.refusal += " Full CI is feasible only for very few electrons in small basis sets; use coupled cluster.";
            else if (r.molecule.name == "H2O" || r.molecule.name == "H2O2" || r.molecule.name == "NH3")
                r.refusal += " Full CI for " + r.molecule.name + " is impossible in every basis set offered; use RHF.";
            else if (r.molecule.name == "LiH")
                r.refusal += " For LiH, Full CI is possible up to aug-cc-pCVTZ (11 million determinants).";
            return r;
        }
        if (e.total_seconds > kLongCalculationSeconds)
            r.warnings.push_back(format(
                "Long calculation: Full CI with %s determinants in %d orbitals after RHF; about %.1f GB of "
                "memory. Estimated time %s on %d threads.",
                thousands(static_cast<std::size_t>(e.fci->determinants)).c_str(), e.fci->orbitals,
                e.fci->memory_bytes / (1024.0 * 1024.0 * 1024.0), human_seconds(e.total_seconds).c_str(), e.threads));
        return r;
    }
    CCMethod ccm;
    if (coupled_cluster_method(r.input.method, &ccm)) {
        const int nel = r.molecule.electrons();
        const double pass = estimate_rhf_seconds(r.basis, false, e.threads, 1);
        e.cc = estimate_cc(ccm, nel / 2, r.basis.n_functions, r.basis.n_functions, static_cast<double>(e.unique_eris),
                           pass, e.threads, r.input.fci_memory_mb, r.input.memory_mb);
        e.total_seconds += e.cc->seconds;
        const int rank = cc_max_rank(ccm);
        if (nel == 2)
            r.notes.push_back(ccm == CCMethod::CCSD ? "For two electrons CCSD is exact (identical to Full CI)."
                              : ccm == CCMethod::CCSD_T
                                  ? "For two electrons CCSD is exact (identical to Full CI); the (T) correction is zero."
                                  : "For two electrons CCSD is exact (identical to Full CI); triple and quadruple "
                                    "excitations vanish, so CCSDT and CCSDTQ give the CCSD energy.");
        else if (nel == 4 && rank == 4)
            r.notes.push_back("For four electrons CCSDTQ is exact (identical to Full CI).");
        if (!e.cc->feasible) {
            r.refusal = e.cc->reason;
            if (rank >= 3) r.refusal += " Use a smaller basis set, or CCSD(T).";
            else r.refusal += " Use a smaller basis set.";
            return r;
        }
        if (e.total_seconds > kLongCalculationSeconds)
            r.warnings.push_back(format(
                "Long calculation: %s with %s amplitudes (%d occupied, %d virtual orbitals) after RHF; about %.1f GB "
                "of memory. Estimated time %s on %d threads.",
                e.cc->method.c_str(), thousands_d(e.cc->amplitudes).c_str(), e.cc->occupied, e.cc->virtuals,
                e.cc->memory_bytes / (1024.0 * 1024.0 * 1024.0), human_seconds(e.total_seconds).c_str(), e.threads));
        return r;
    }
    if (e.rhf_seconds > kLongCalculationSeconds)
        r.warnings.push_back(format(
            "Long calculation: %d basis functions (up to %c functions); %s two-electron integrals %s. "
            "Estimated time %s on %d threads.",
            e.basis_functions, angular_letter(e.max_l), thousands(e.unique_eris).c_str(),
            e.eri_stored ? format("kept in memory (%.0f MB)", e.eri_memory_mb).c_str()
                         : format("recomputed in every SCF iteration (storing them would need %.1f GB)",
                                  e.eri_memory_mb / 1024.0)
                               .c_str(),
            human_seconds(e.rhf_seconds).c_str(), e.threads));
    return r;
}

void run_molecular(MolecularResult& r, const RHFProgress& progress, const FCIProgress& fci_progress,
                   const CCProgress& cc_progress) {
    if (!r.refusal.empty()) throw std::runtime_error(r.refusal);
    if (r.reference != "RHF" || uses_so_engine(r.input.method, r.reference)) {
        // 15.0: open shells and the methods of the spin-orbital engine
        ElectronicOptions eo;
        eo.threads = r.input.threads;
        eo.memory_mb = r.input.memory_mb;
        eo.fci_memory_mb = r.input.fci_memory_mb;
        eo.scf_max_iter = r.input.max_iter;
        eo.scf_energy_tol = r.input.energy_tol;
        eo.scf_gradient_tol = r.input.gradient_tol;
        eo.cc_amplitude_tol = 1.0e-8;
        eo.reference = r.reference;
        eo.multiplicity = r.input.multiplicity;
        ElectronicCallbacks cb;
        cb.scf = progress;
        cb.fci = fci_progress;
        cb.cc = cc_progress;
        SinglePoint sp = run_single_point(r.molecule, r.basis, r.input.method, eo, nullptr, cb);
        r.rhf = sp.rhf;
        r.so = sp.so;
        r.fci = sp.fci;
        r.cc = sp.cc;
        if (r.rhf->removed_functions > 0)
            r.warnings.push_back(format("%d near-linearly dependent combination(s) of basis functions removed "
                                        "(overlap eigenvalues below %.0e).",
                                        r.rhf->removed_functions, 1.0e-7));
        for (const auto& w : sp.warnings) r.warnings.push_back(w);
        r.report = molecular_report(r);
        return;
    }
    RHFOptions opt;
    opt.max_iter = r.input.max_iter;
    opt.energy_tol = r.input.energy_tol;
    opt.gradient_tol = r.input.gradient_tol;
    opt.threads = r.input.threads;
    opt.memory_limit_mb = r.input.memory_mb;
    r.rhf = run_rhf(r.molecule, r.basis, opt, progress);
    if (r.rhf->removed_functions > 0)
        r.warnings.push_back(format("%d near-linearly dependent combination(s) of basis functions removed "
                                    "(overlap eigenvalues below %.0e).",
                                    r.rhf->removed_functions, opt.lindep_tol));
    if (!r.rhf->converged)
        r.warnings.push_back(format("The SCF did not converge in %d iterations.", r.rhf->iterations));
    if (r.input.method == "FCI" && r.rhf->converged) {
        const auto t0 = std::chrono::steady_clock::now();
        const MOIntegrals mo = transform_to_mo(r.basis, r.rhf->H, r.rhf->C, r.input.threads);
        const double transform = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        FCIOptions fo;
        fo.threads = r.input.threads;
        fo.memory_limit_mb = r.input.fci_memory_mb;
        const int nocc = r.rhf->n_occupied;
        r.fci = run_fci_determinants(mo, nocc, nocc, r.rhf->nuclear_repulsion, r.rhf->total_energy, fo, fci_progress);
        r.fci->transform_seconds = transform;
        if (!r.fci->converged)
            r.warnings.push_back(format("The Full CI (Davidson) did not converge in %d iterations.", r.fci->iterations));
    }
    CCMethod ccm;
    if (coupled_cluster_method(r.input.method, &ccm) && r.rhf->converged) {
        CCOptions co;
        co.threads = r.input.threads;
        co.memory_limit_mb = r.input.fci_memory_mb;
        co.eri_memory_mb = r.input.memory_mb;
        const int nocc = r.rhf->n_occupied;
        const auto t0 = std::chrono::steady_clock::now();
        if (cc_max_rank(ccm) == 2) {
            const CCIntegrals ints = cc_integrals_from_ao(r.basis, r.rhf->H, r.rhf->C, nocc, r.input.threads);
            auto ladder = make_ao_ladder(r.basis, r.rhf->C, nocc, r.input.threads, r.input.memory_mb);
            const double transform = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            r.cc = run_ccsd_closed_shell(ints, *ladder, r.rhf->total_energy, ccm == CCMethod::CCSD_T, co, cc_progress);
            r.cc->transform_seconds = transform;
        } else {
            const MOIntegrals mo = transform_to_mo(r.basis, r.rhf->H, r.rhf->C, r.input.threads);
            const double transform = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            r.cc = run_cc_spin_orbital(mo, nocc, cc_max_rank(ccm), r.rhf->nuclear_repulsion, co, cc_progress);
            r.cc->transform_seconds = transform;
            // Report against the SCF energy (the same determinant; equal to ~1e-10 E_h).
            r.cc->reference_energy = r.rhf->total_energy;
            r.cc->cc_energy = r.cc->reference_energy + r.cc->correlation_energy;
            r.cc->total_energy = r.cc->cc_energy;
        }
        if (!r.cc->converged)
            r.warnings.push_back(format("The %s iterations did not converge in %d iterations.", r.cc->method.c_str(),
                                        r.cc->iterations));
        if (ccm == CCMethod::CCSD_T && !r.cc->converged)
            r.warnings.push_back("The (T) correction was not computed because CCSD did not converge.");
    }
    r.report = molecular_report(r);
}

std::string molecular_report(const MolecularResult& r) {
    const int level = r.input.print_level;
    if (level <= 0) return "";
    const std::string rule(78, '=');
    const Molecule& mol = r.molecule;
    std::string s;
    const std::string scf_name = r.reference == "UHF" ? "unrestricted Hartree-Fock"
                                 : r.reference == "ROHF" ? "restricted open-shell Hartree-Fock"
                                                          : "restricted Hartree-Fock";
    s += std::string("\xCE\xA8Lab \xE2\x80\x93 ") + scf_name +
         (r.input.method == "FCI" ? " and Full CI"
                                  : (r.input.method == "RHF" ? "" : " and " + method_label(r.input.method, r.reference))) +
         "\n" + rule + "\n";
    s += format("Molecule        %s, %d electrons%s, %s (charge %d, multiplicity %d)\n", mol.display_name().c_str(),
                mol.electrons(),
                mol.core_electrons() ? format(" (and %d core electrons in ECPs)", mol.core_electrons()).c_str() : "",
                r.input.multiplicity == 1 ? "closed shell" : "open shell", mol.charge, r.input.multiplicity);
    if (mol.custom)
        s += "Title           " + (mol.title.empty() ? std::string("(none)") : mol.title) + "\n" +
             format("Geometry        %zu atoms from Cartesian coordinates\n", mol.atoms.size());
    else if (mol.name == "H2O2")
        s += "Geometry        R(O-O) = " + angstrom_bohr(r.input.r_oo) + ", R(O-H) = " + angstrom_bohr(r.input.R) +
             format(", angle O-O-H = %.4f deg, dihedral H-O-O-H = %.4f deg\n", r.input.angle, r.input.dihedral);
    else if (mol.name == "NH3")
        s += "Geometry        R(N-H) = " + angstrom_bohr(r.input.R) + format(", angle H-N-H = %.4f deg (C3v)\n", r.input.angle);
    else if (mol.name == "H2O")
        s += "Geometry        R(O-H) = " + angstrom_bohr(r.input.R) + format(", angle H-O-H = %.4f deg\n", r.input.angle);
    else if (mol.name == "LiH")
        s += "Geometry        R(Li-H) = " + angstrom_bohr(r.input.R) + "\n";
    else if (mol.name == "HeH+")
        s += "Geometry        R(He-H) = " + angstrom_bohr(r.input.R) + "\n";
    else
        s += "Geometry        R(H-H) = " + angstrom_bohr(r.input.R) + "\n";
    s += "Cartesian coordinates (bohr)\n";
    for (const Atom& a : mol.atoms)
        s += format("    %-2s  %16.10f  %16.10f  %16.10f\n", a.symbol.c_str(), a.position[0], a.position[1],
                    a.position[2]);
    s += format("Nuclear repulsion energy  %20.10f E_h\n\n", mol.nuclear_repulsion());

    const BasisSet& b = r.basis;
    s += "Basis set       " + b.name + " (" + basis_function_types(b.name) + "; basissetexchange.org, Gaussian94 format)\n";
    int width = 13;  // 14.0: names up to "6-311++G(3df,3pd)"
    for (const auto& a : b.atoms) width = std::max(width, static_cast<int>(a.basis.size()));
    std::vector<std::string> seen;
    for (const auto& a : b.atoms) {
        if (std::find(seen.begin(), seen.end(), a.symbol) != seen.end()) continue;
        seen.push_back(a.symbol);
        s += format("    %-2s  %-*s %-16s -> %-16s %4d functions%s\n", a.symbol.c_str(), width, a.basis.c_str(),
                    a.primitives.c_str(), a.contraction.c_str(), a.functions,
                    a.ecp.empty() ? "" : (", " + a.ecp).c_str());
    }
    for (const auto& n : r.notes) s += "    Note: " + n + "\n";
    if (b.ecp_electrons > 0)
        s += format("    Effective core potentials (ECPn: n core electrons replaced): %d core electrons; the\n"
                    "    nuclear charges are Z less the core electrons\n",
                    b.ecp_electrons);
    s += format("    %d basis functions, %zu shells, highest angular momentum %c (l = %d)\n", b.n_functions,
                b.shells.size(), angular_letter(b.max_l), b.max_l);
    const CostEstimate& e = r.estimate;
    s += format("    two-electron integrals: %s unique, ", thousands(e.unique_eris).c_str());
    s += e.eri_stored ? format("kept in memory (%.1f MB)\n", e.eri_memory_mb)
                      : format("integral-direct (storing them would need %.1f GB)\n", e.eri_memory_mb / 1024.0);
    if (e.fci) {
        const FCISize& f = *e.fci;
        s += format("    Full CI: %d electrons (%d alpha, %d beta) in %d orbitals, %s x %s = %s determinants,\n"
                    "             about %.0f MB of memory, estimated time %s\n",
                    f.alpha + f.beta, f.alpha, f.beta, f.orbitals, thousands_d(f.alpha_strings).c_str(),
                    thousands_d(f.beta_strings).c_str(), thousands_d(f.determinants).c_str(),
                    f.memory_bytes / (1024.0 * 1024.0), human_seconds(e.total_seconds).c_str());
    }
    if (e.cc) {
        const CCSize& c = *e.cc;
        s += format("    %s: %d occupied, %d virtual orbitals, %s amplitudes, about %.0f MB of memory,\n"
                    "             estimated time %s\n",
                    c.method.c_str(), c.occupied, c.virtuals, thousands_d(c.amplitudes).c_str(),
                    c.memory_bytes / (1024.0 * 1024.0), human_seconds(e.total_seconds).c_str());
    }
    if (e.so)
        s += format("    %s: %s amplitudes (dense spin blocks), about %.0f MB of memory, estimated time %s\n",
                    e.so->method.c_str(), thousands_d(e.so->amplitudes).c_str(), e.so->memory_bytes / 1048576.0,
                    human_seconds(e.total_seconds).c_str());
    for (const auto& w : r.warnings) s += "Warning: " + w + "\n";
    if (!r.refusal.empty()) s += "NOT POSSIBLE: " + r.refusal + "\n";
    if (!r.rhf) return s;

    const RHFResult& h = *r.rhf;
    s += format("    overlap matrix: smallest eigenvalue %.3e; %s\n", h.smallest_overlap_eigenvalue,
                h.removed_functions ? format("%d functions removed", h.removed_functions).c_str()
                                    : "no functions removed");
    const char* guess = h.guess == "sad"        ? "superposition of atomic densities (SAD)"
                        : h.guess == "orbitals" ? "orbitals of the previous geometry"
                                                : "core-Hamiltonian guess";
    s += format("\nSCF             %s, DIIS (up to 8 Fock matrices)\n", guess);
    s += format("    converged when |dE| < %.1e E_h and max |X'(FPS - SPF)X| < %.1e\n", r.input.energy_tol,
                r.input.gradient_tol);
    if (level >= 2) {
        s += "    iter        total energy (E_h)        dE (E_h)      gradient\n";
        for (const auto& it : h.history)
            s += format("    %4d  %24.12f  %14.3e  %12.3e%s\n", it.iteration, it.energy, it.delta_e, it.gradient,
                        it.diis ? "  DIIS" : "");
    }
    s += h.converged ? format("    SCF converged in %d iterations\n", h.iterations)
                     : format("    SCF NOT converged after %d iterations\n", h.iterations);
    s += format("    time: integrals %.2f s, SCF %.2f s, %d threads\n\n", h.integral_seconds, h.scf_seconds, h.threads);

    if (h.open_shell())
        s += format("Spin            %d alpha, %d beta electrons; <S^2> = %.6f (S(S+1) = %.6f)\n\n", h.n_alpha, h.n_beta,
                    h.s2, h.s2_exact);
    s += "Energies (E_h)\n";
    s += format("    total energy (%s)%*s%20.10f\n", r.reference.c_str(), static_cast<int>(11 - r.reference.size()), "",
                h.total_energy);
    s += format("    electronic energy         %20.10f\n", h.electronic_energy);
    s += format("    nuclear repulsion         %20.10f\n", h.nuclear_repulsion);
    s += format("    kinetic energy            %20.10f\n", h.kinetic_energy);
    s += format("    nuclear attraction        %20.10f\n", h.nuclear_attraction_energy);
    s += format("    electron repulsion        %20.10f\n", h.two_electron_energy);
    s += format("    virial ratio -V/T         %20.10f\n\n", h.virial_ratio);

    auto orbitals = [&](const std::vector<double>& eps, int nocc, const std::string& title) {
        std::string t;
        const int nshow = level >= 2 ? static_cast<int>(eps.size()) : std::min(static_cast<int>(eps.size()), nocc + 5);
        t += format("%s (E_h), %d occupied (*)%s\n", title.c_str(), nocc,
                    nshow < static_cast<int>(eps.size()) ? ", lowest virtual orbitals" : "");
        for (int k = 0; k < nshow; ++k) {
            t += format("  %4d%s %15.8f", k + 1, k < nocc ? "*" : " ", eps[k]);
            t += (k % 4 == 3 || k + 1 == nshow) ? "\n" : "";
        }
        if (nocc >= 1 && nocc < static_cast<int>(eps.size()))
            t += format("    HOMO %.8f E_h, LUMO %.8f E_h, gap %.8f E_h\n", eps[nocc - 1], eps[nocc], eps[nocc] - eps[nocc - 1]);
        return t;
    };
    if (h.reference == "UHF") {
        s += orbitals(h.orbital_energies, h.n_alpha, "Alpha orbital energies");
        s += orbitals(h.orbital_energies_b, h.n_beta, "Beta orbital energies");
    } else if (h.reference == "ROHF") {
        s += orbitals(h.orbital_energies, h.n_alpha,
                      format("ROHF orbital energies (Roothaan effective Fock matrix; %d doubly, %d singly occupied)",
                             h.n_beta, h.n_alpha - h.n_beta));
    } else {
        s += orbitals(h.orbital_energies, h.n_occupied, "Orbital energies");
    }
    if (r.fci) s += fci_section_impl(*r.fci, r.rhf->total_energy, level);
    if (r.cc) s += cc_section_impl(*r.cc, level);
    if (r.so) s += so_report_section(*r.so, level);
    return s;
}

}  // namespace macqc
