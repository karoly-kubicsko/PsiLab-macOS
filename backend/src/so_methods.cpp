// ΨLab 15.0: the correlated methods of the spin-orbital engine (see so_methods.h).
#include "macqc/so_methods.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

#include "macqc/gto_integrals.h"
#include "macqc/job.h"

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;
using so::BlockTensor;
using so::Spaces;

struct TermSpec {
    char kind;
    int target;
    const char* block;
    double coef;
    int nclusters;
    int cl[4][3];
    const char* hletters;
    const char* steps[4];
    const char* zletters;
    const char* cosets;
};

const TermSpec kTerms[] = {
#include "macqc/so_diagrams.inc"
};

const char* kPart = "ABCD";
const char* kHole = "IJKL";

std::string target_letters(int m) { return std::string(kPart).substr(0, static_cast<std::size_t>(m)) + std::string(kHole).substr(0, static_cast<std::size_t>(m)); }
std::string amp_modes(int m) { return std::string(static_cast<std::size_t>(m), 'v') + std::string(static_cast<std::size_t>(m), 'o'); }
double factorial(int m) {
    double f = 1.0;
    for (int k = 2; k <= m; ++k) f *= k;
    return f;
}

struct MethodInfo {
    std::string family;  // "CC", "QCI", "CI", "MP"
    int rank = 2;        // CC/QCI/CI: highest rank; MP: order n
    bool triples = false;  // (T)
    double triples_singles = 1.0;
};

bool method_info(const std::string& m, MethodInfo& info) {
    static const std::map<std::string, MethodInfo> table = {
        {"MP2", {"MP", 2, false, 0}},       {"MP3", {"MP", 3, false, 0}},       {"MP4", {"MP", 4, false, 0}},
        {"MP5", {"MP", 5, false, 0}},       {"CISD", {"CI", 2, false, 0}},      {"CISD(T)", {"CI", 2, true, 2.0}},
        {"CISDT", {"CI", 3, false, 0}},     {"CISDTQ", {"CI", 4, false, 0}},    {"QCISD", {"QCI", 2, false, 0}},
        {"QCISD(T)", {"QCI", 2, true, 2.0}}, {"QCISDT", {"QCI", 3, false, 0}},  {"QCISDTQ", {"QCI", 4, false, 0}},
        {"CCSD", {"CC", 2, false, 0}},      {"CCSD(T)", {"CC", 2, true, 1.0}},  {"CCSDT", {"CC", 3, false, 0}},
        {"CCSDTQ", {"CC", 4, false, 0}},
    };
    const auto it = table.find(m);
    if (it == table.end()) return false;
    info = it->second;
    return true;
}

// allowed spin patterns of an amplitude of rank m (beta particles = beta holes)
std::vector<unsigned> amplitude_keys(int m) {
    std::vector<unsigned> keys;
    const unsigned all = 1u << (2 * m);
    for (unsigned k = 0; k < all; ++k) {
        int bp = 0, bh = 0;
        for (int i = 0; i < m; ++i) {
            bp += (k >> i) & 1u;
            bh += (k >> (m + i)) & 1u;
        }
        if (bp == bh) keys.push_back(k);
    }
    return keys;
}

}  // namespace

// ------------------------------------------------------------------------------------------------

double SpatialHamiltonian::g(int sp, int sq, std::size_t p, std::size_t q, std::size_t r, std::size_t s) const {
    if (sp == 0 && sq == 0) return (*aa)(p, q, r, s);
    if (sp == 1 && sq == 1) return (*bb)(p, q, r, s);
    if (sp == 0) return ab->mixed ? ab->mixed_value(p, q, r, s) : (*ab)(p, q, r, s);
    return ab->mixed ? ab->mixed_value(r, s, p, q) : (*ab)(r, s, p, q);
}

bool so_engine_method(const std::string& canonical) {
    MethodInfo i;
    return method_info(canonical, i);
}

int so_method_max_rank(const std::string& canonical) {
    MethodInfo i;
    if (!method_info(canonical, i)) return 0;
    if (i.family == "MP") return i.rank == 5 ? 4 : (i.rank == 4 ? 3 : 2);
    return i.triples ? 3 : i.rank;
}

// ------------------------------------------------------------------------------------------------

struct Var {
    std::string name;   // "T2", "T2(1)", "C3", "T3[T]"
    int rank = 2;
    int order = 0;      // MP order (0 otherwise)
    int stage = 0;
    bool ci = false;    // CI: denominators shifted by E_c
    int input = -1, residual = -1, lambda = -1;
};

class SOMethodState {
public:
    std::string method, reference;
    MethodInfo info;
    Spaces sp;
    so::Graph g;
    std::unique_ptr<so::Engine> eng;
    std::vector<Var> vars;
    std::vector<BlockTensor> t, lambda;
    std::vector<double> eo[2], ev[2];  // orbital-energy diagonals (denominators)
    int nstages = 1;
    int hin[4] = {-1, -1, -1, -1};  // oo, ov, vo, vv
    std::map<std::string, int> vin;  // "oo|oo" ...
    int delta = -1;
    int fock[4] = {-1, -1, -1, -1};
    int eref = -1, energy = -1, ec = -1;
    std::vector<int> stage_energy;
    std::vector<std::pair<std::string, int>> energy_terms;
    int n[2] = {0, 0}, nocc[2] = {0, 0};
    BlockTensor hval[4];
    std::map<std::string, BlockTensor> vval;
    bool lambda_done = false;
    SODensities dens;
};

namespace {

const char* kHBlocks[4] = {"oo", "ov", "vo", "vv"};
const char* kVParts[3] = {"oo", "vo", "vv"};

int fock_index(const std::string& ul) {
    for (int k = 0; k < 4; ++k)
        if (ul == kHBlocks[k]) return k;
    throw std::logic_error("fock block " + ul);
}

// The Hamiltonian block node of a term ("F:o|v" -> f_ov node, "V:vo|vv" -> input).
int hblock(SOMethodState& s, const std::string& block) {
    const std::string body = block.substr(2);
    const auto bar = body.find('|');
    const std::string u = body.substr(0, bar), l = body.substr(bar + 1);
    if (block[0] == 'F') return s.fock[fock_index(u + l)];
    return s.vin.at(u + "|" + l);
}

int build_term(SOMethodState& s, const TermSpec& t, const std::vector<int>& amps) {
    int z = hblock(s, t.block);
    for (int k = 0; k < t.nclusters; ++k) z = s.g.einsum(t.steps[k], z, amps[static_cast<std::size_t>(k)]);
    if (t.target == 0) return s.g.sum({{t.coef, z}});
    const std::string target = target_letters(t.target);
    std::vector<std::pair<std::string, std::string>> perms;
    std::vector<double> coefs;
    std::istringstream in(t.cosets);
    std::string tok;
    while (in >> tok) {
        coefs.push_back(t.coef * (tok.back() == '-' ? -1.0 : 1.0));
        perms.emplace_back(tok.substr(0, tok.size() - 1), target);
    }
    return s.g.permsum(z, perms, coefs);
}

void build_hamiltonian_nodes(SOMethodState& s) {
    so::Graph& g = s.g;
    for (int k = 0; k < 4; ++k) s.hin[k] = g.input(std::string("h:") + kHBlocks[k], kHBlocks[k]);
    for (const char* u : kVParts)
        for (const char* l : kVParts) {
            const std::string key = std::string(u) + "|" + l;
            s.vin[key] = g.input("V:" + key, std::string(u) + l);
        }
    s.delta = g.input("delta", "oo", false);
    const int foo = g.einsum("kilj,ij->kl", s.vin["oo|oo"], s.delta);
    s.fock[0] = g.sum({{1.0, s.hin[0]}, {1.0, foo}});
    s.fock[1] = g.sum({{1.0, s.hin[1]}, {1.0, g.einsum("ijak,jk->ia", s.vin["oo|vo"], s.delta)}});
    s.fock[2] = g.sum({{1.0, s.hin[2]}, {1.0, g.einsum("ajik,jk->ai", s.vin["vo|oo"], s.delta)}});
    s.fock[3] = g.sum({{1.0, s.hin[3]}, {1.0, g.einsum("ajbk,jk->ab", s.vin["vo|vo"], s.delta)}});
    const int e1 = g.einsum("ij,ij->", s.hin[0], s.delta);
    const int e2 = g.einsum("kl,kl->", foo, s.delta);
    s.eref = g.sum({{1.0, e1}, {0.5, e2}});
}

int add_var(SOMethodState& s, const std::string& name, int rank, int order, int stage, bool ci) {
    Var v;
    v.name = name;
    v.rank = rank;
    v.order = order;
    v.stage = stage;
    v.ci = ci;
    v.input = s.g.input("t:" + name, amp_modes(rank));
    v.lambda = s.g.input("l:" + name, amp_modes(rank));
    s.vars.push_back(v);
    return static_cast<int>(s.vars.size()) - 1;
}

int full_dot(so::Graph& g, int a, int b, int rank) {
    const std::string l = target_letters(rank);
    return g.einsum(l + "," + l + "->", a, b);
}

// Residuals and energies of the CC / QCI / CI families.
void build_cluster_method(SOMethodState& s) {
    const MethodInfo& mi = s.info;
    const int nmax = mi.rank;
    const bool ci = mi.family == "CI";
    std::vector<int> amp(5, -1);
    std::vector<int> vidx(5, -1);
    for (int m = 1; m <= nmax; ++m) {
        vidx[static_cast<std::size_t>(m)] = add_var(s, (ci ? "C" : "T") + std::to_string(m), m, 0, 0, ci);
        amp[static_cast<std::size_t>(m)] = s.vars[static_cast<std::size_t>(vidx[static_cast<std::size_t>(m)])].input;
    }
    auto included = [&](const TermSpec& t, int m) {
        for (int k = 0; k < t.nclusters; ++k)
            if (t.cl[k][0] > nmax) return false;
        if (mi.family == "CC") return t.kind == 'C';
        if (mi.family == "CI") return (t.kind == 'C' && t.nclusters <= 1) || t.kind == 'D';
        // QCI
        if (t.kind != 'C') return false;
        if (t.nclusters <= 1) return true;
        if (t.nclusters != 2 || m == 0) return false;
        int r1 = t.cl[0][0], r2 = t.cl[1][0];
        if (r1 > r2) std::swap(r1, r2);
        if (r1 == 2 && r2 == 2 && m == 2) return true;
        if (std::min(r1, r2) == 2 && std::max(r1, r2) == m) return true;
        if (nmax == 4 && m == 4 && r1 == 2 && r2 == 2) return true;
        if (m == 1 && r1 == 1 && r2 == 2) return true;
        return false;
    };
    std::vector<std::vector<std::pair<double, int>>> res(static_cast<std::size_t>(nmax) + 1);
    std::vector<std::pair<double, int>> ecorr;
    for (const TermSpec& t : kTerms) {
        if (t.target > nmax) continue;
        if (!included(t, t.target)) continue;
        if (t.target == 0 && mi.family != "CC" && t.nclusters > 1) continue;
        std::vector<int> a;
        for (int k = 0; k < t.nclusters; ++k) a.push_back(amp[static_cast<std::size_t>(t.cl[k][0])]);
        const int node = build_term(s, t, a);
        if (t.target == 0) ecorr.emplace_back(1.0, node);
        else res[static_cast<std::size_t>(t.target)].emplace_back(1.0, node);
    }
    const int ec = s.g.sum(ecorr);
    if (ci) {
        s.ec = ec;
        for (int m = 1; m <= nmax; ++m)
            res[static_cast<std::size_t>(m)].emplace_back(-1.0, s.g.scale(ec, amp[static_cast<std::size_t>(m)]));
    }
    for (int m = 1; m <= nmax; ++m)
        s.vars[static_cast<std::size_t>(vidx[static_cast<std::size_t>(m)])].residual = s.g.sum(res[static_cast<std::size_t>(m)]);
    const int ebase = s.g.sum({{1.0, s.eref}, {1.0, ec}});
    s.energy_terms.emplace_back(ci ? "E_c (CI correlation)" : "correlation", ec);
    s.stage_energy.push_back(ebase);
    s.energy = ebase;
    s.nstages = 1;
    if (mi.triples) {
        // (T): t3 with <T|[F_oo + F_vv, T3] + (V T2)_C|0> = 0, E_T = <t3|X> + f <t3|Y>
        const int v3 = add_var(s, ci ? "C3[T]" : "T3[T]", 3, 0, 1, false);
        const int t3 = s.vars[static_cast<std::size_t>(v3)].input;
        std::vector<std::pair<double, int>> lin, X, Y;
        for (const TermSpec& t : kTerms) {
            if (t.target != 3) continue;
            const std::string b = t.block;
            if (t.kind == 'C' && t.nclusters == 1 && t.cl[0][0] == 3 && (b == "F:o|o" || b == "F:v|v"))
                lin.emplace_back(1.0, build_term(s, t, {t3}));
            else if (t.kind == 'C' && t.nclusters == 1 && t.cl[0][0] == 2 && b[0] == 'V')
                X.emplace_back(1.0, build_term(s, t, {amp[2]}));
            else if (t.kind == 'D' && t.nclusters == 1 && t.cl[0][0] <= 2)
                Y.emplace_back(1.0, build_term(s, t, {amp[static_cast<std::size_t>(t.cl[0][0])]}));
        }
        const int xn = s.g.sum(X), yn = s.g.sum(Y);
        lin.emplace_back(1.0, xn);
        s.vars[static_cast<std::size_t>(v3)].residual = s.g.sum(lin);
        const double w = 1.0 / 36.0;
        const int e4 = s.g.sum({{w, full_dot(s.g, t3, xn, 3)}});
        const int e5 = s.g.sum({{w * mi.triples_singles, full_dot(s.g, t3, yn, 3)}});
        s.energy_terms.emplace_back("E[4]_T", e4);
        s.energy_terms.emplace_back("E[5]_ST", e5);
        s.energy = s.g.sum({{1.0, ebase}, {1.0, e4}, {1.0, e5}});
        s.stage_energy.push_back(s.energy);
        s.nstages = 2;
    }
}

// Moller-Plesset: amplitudes order by order.
void build_mp(SOMethodState& s) {
    const int n = s.info.rank;
    // ranks needed at each order (see so_methods.h)
    std::vector<std::set<int>> need(static_cast<std::size_t>(n));
    need[static_cast<std::size_t>(n - 1)] = {1, 2};
    for (int j = n - 2; j >= 1; --j) {
        const int top = std::min({j + 1, *need[static_cast<std::size_t>(j + 1)].rbegin() + 2, 4});
        for (int r = 1; r <= top; ++r) need[static_cast<std::size_t>(j)].insert(r);
    }
    std::map<std::pair<int, int>, int> amp;  // (order, rank) -> input
    std::map<std::pair<int, int>, int> var;
    for (int j = 1; j <= n - 1; ++j)
        for (int r : need[static_cast<std::size_t>(j)]) {
            const int v = add_var(s, "T" + std::to_string(r) + "(" + std::to_string(j) + ")", r, j, j - 1, false);
            var[{j, r}] = v;
            amp[{j, r}] = s.vars[static_cast<std::size_t>(v)].input;
        }
    std::map<std::pair<int, int>, std::vector<std::pair<double, int>>> res;
    std::map<int, std::vector<std::pair<double, int>>> eorder;
    for (const TermSpec& t : kTerms) {
        if (t.kind != 'C') continue;
        const std::string b = t.block;
        const int ov = (b == "F:o|o" || b == "F:v|v") ? 0 : 1;
        const int k = t.nclusters;
        // all assignments of orders (>= 1) to the clusters
        std::vector<int> ord(static_cast<std::size_t>(k), 1);
        const int maxorder = n;  // energies up to order n
        std::function<void(int)> rec = [&](int c) {
            if (c == k) {
                int total = ov;
                for (int o : ord) total += o;
                std::vector<int> a;
                for (int q = 0; q < k; ++q) {
                    const auto it = amp.find({ord[static_cast<std::size_t>(q)], t.cl[q][0]});
                    if (it == amp.end()) return;
                    a.push_back(it->second);
                }
                if (t.target == 0) {
                    if (total >= 2 && total <= maxorder) eorder[total].emplace_back(1.0, build_term(s, t, a));
                } else if (total >= 1 && total <= n - 1 && var.count({total, t.target})) {
                    res[{total, t.target}].emplace_back(1.0, build_term(s, t, a));
                }
                return;
            }
            for (int o = 1; o <= n - 1; ++o) {
                ord[static_cast<std::size_t>(c)] = o;
                rec(c + 1);
            }
        };
        rec(0);
    }
    for (auto& kv : var) s.vars[static_cast<std::size_t>(kv.second)].residual = s.g.sum(res[kv.first]);
    std::vector<std::pair<double, int>> etot = {{1.0, s.eref}};
    for (int order = 2; order <= n; ++order) {
        const int e = s.g.sum(eorder[order]);
        s.energy_terms.emplace_back("E(" + std::to_string(order) + ")", e);
        etot.emplace_back(1.0, e);
        // the stage of order j - 1 completes the energy up to order j
        s.stage_energy.push_back(s.g.sum(etot));
    }
    s.energy = s.stage_energy.back();
    s.nstages = n - 1;
}

// ---- inputs ---------------------------------------------------------------------------------

std::size_t spatial_index(char space, int spin, std::size_t i, const int* nocc) {
    return space == 'o' ? i : static_cast<std::size_t>(nocc[spin]) + i;
}

using GFunc = std::function<double(int, int, std::size_t, std::size_t, std::size_t, std::size_t)>;

BlockTensor make_h_block(const std::string& ul, const Matrix* h, const Spaces& sp, const int* nocc) {
    BlockTensor b(ul);
    for (int s = 0; s < 2; ++s) {
        const unsigned key = s ? 3u : 0u;
        Tensor& t = b.block(key, sp);
        const std::size_t d0 = sp.dim(ul[0], s), d1 = sp.dim(ul[1], s);
        for (std::size_t p = 0; p < d0; ++p)
            for (std::size_t q = 0; q < d1; ++q)
                t[p * d1 + q] = h[s](spatial_index(ul[0], s, p, nocc), spatial_index(ul[1], s, q, nocc));
    }
    return b;
}

BlockTensor make_v_block(const std::string& modes, const GFunc& g, const Spaces& sp, const int* nocc) {
    BlockTensor b(modes);
    for (unsigned key = 0; key < 16; ++key) {
        const int s0 = key & 1u, s1 = (key >> 1) & 1u, s2 = (key >> 2) & 1u, s3 = (key >> 3) & 1u;
        if (s0 + s1 != s2 + s3) continue;
        const bool coul = s0 == s2 && s1 == s3, exch = s0 == s3 && s1 == s2;
        if (!coul && !exch) continue;
        Tensor& t = b.block(key, sp);
        const std::size_t d0 = sp.dim(modes[0], s0), d1 = sp.dim(modes[1], s1), d2 = sp.dim(modes[2], s2),
                          d3 = sp.dim(modes[3], s3);
        std::size_t pos = 0;
        for (std::size_t P = 0; P < d0; ++P) {
            const std::size_t p = spatial_index(modes[0], s0, P, nocc);
            for (std::size_t Q = 0; Q < d1; ++Q) {
                const std::size_t q = spatial_index(modes[1], s1, Q, nocc);
                for (std::size_t R = 0; R < d2; ++R) {
                    const std::size_t r = spatial_index(modes[2], s2, R, nocc);
                    for (std::size_t S = 0; S < d3; ++S, ++pos) {
                        const std::size_t sidx = spatial_index(modes[3], s3, S, nocc);
                        double v = 0.0;
                        if (coul) v += g(s0, s1, p, r, q, sidx);
                        if (exch) v -= g(s0, s1, p, sidx, q, r);
                        t[pos] = v;
                    }
                }
            }
        }
    }
    return b;
}

BlockTensor make_delta(const Spaces& sp) {
    BlockTensor d("oo");
    for (int s = 0; s < 2; ++s) {
        Tensor& t = d.block(s ? 3u : 0u, sp);
        const std::size_t o = sp.dim('o', s);
        for (std::size_t i = 0; i < o; ++i) t[i * o + i] = 1.0;
    }
    return d;
}

void set_hamiltonian(SOMethodState& s, const Matrix* h, const GFunc& g) {
    for (int k = 0; k < 4; ++k) {
        s.hval[k] = make_h_block(kHBlocks[k], h, s.sp, s.nocc);
        s.eng->set_input(s.hin[k], s.hval[k]);
    }
    for (auto& kv : s.vin) {
        const std::string modes = kv.first.substr(0, 2) + kv.first.substr(3);
        s.vval[kv.first] = make_v_block(modes, g, s.sp, s.nocc);
        s.eng->set_input(kv.second, s.vval[kv.first]);
    }
    s.eng->set_input(s.delta, make_delta(s.sp));
}

// ---- amplitude arithmetic -------------------------------------------------------------------

BlockTensor zero_amplitude(int m, const Spaces& sp) {
    BlockTensor b(amp_modes(m));
    for (unsigned k : amplitude_keys(m)) b.block(k, sp);
    return b;
}

// out = -R / (D - shift), D = sum eps_v - sum eps_o, on the allowed blocks of rank m; returns max |out|
double jacobi_step(const BlockTensor& R, int m, const Spaces& sp, const std::vector<double>* eo,
                   const std::vector<double>* ev, double shift, BlockTensor& out) {
    out = zero_amplitude(m, sp);
    double mx = 0.0;
    for (auto& kv : out.blocks()) {
        const Tensor* r = R.find(kv.first);
        if (!r) continue;
        const unsigned key = kv.first;
        const std::vector<std::size_t>& dims = kv.second.dims();
        std::vector<std::size_t> idx(dims.size(), 0);
        for (std::size_t pos = 0; pos < kv.second.size(); ++pos) {
            double d = -shift;
            for (int k = 0; k < m; ++k) d += ev[(key >> k) & 1u][idx[static_cast<std::size_t>(k)]];
            for (int k = 0; k < m; ++k) d -= eo[(key >> (m + k)) & 1u][idx[static_cast<std::size_t>(m + k)]];
            if (std::fabs(d) < 1.0e-8) d = d < 0 ? -1.0e-8 : 1.0e-8;
            const double v = -(*r)[pos] / d;
            kv.second[pos] = v;
            mx = std::max(mx, std::fabs(v));
            for (std::size_t q = dims.size(); q-- > 0;) {
                if (++idx[q] < dims[q]) break;
                idx[q] = 0;
            }
        }
    }
    return mx;
}

double amplitude_count(int m, const Spaces& sp) {
    double n = 0.0;
    for (unsigned k : amplitude_keys(m)) {
        double b = 1.0;
        for (std::size_t d : so::block_dims(amp_modes(m), k, sp)) b *= static_cast<double>(d);
        n += b;
    }
    return n;
}

void flatten(const std::vector<const BlockTensor*>& xs, std::vector<double>& out) {
    out.clear();
    for (const BlockTensor* x : xs)
        for (const auto& kv : x->blocks()) out.insert(out.end(), kv.second.vec().begin(), kv.second.vec().end());
}

void unflatten(const std::vector<double>& in, const std::vector<BlockTensor*>& xs) {
    std::size_t pos = 0;
    for (BlockTensor* x : xs)
        for (auto& kv : x->blocks()) {
            std::copy(in.begin() + static_cast<std::ptrdiff_t>(pos),
                      in.begin() + static_cast<std::ptrdiff_t>(pos + kv.second.size()), kv.second.vec().begin());
            pos += kv.second.size();
        }
}

// DIIS (Pulay) on flattened vectors; returns the extrapolated vector or the last one.
class Diis {
public:
    explicit Diis(int size) : size_(size) {}
    std::vector<double> extrapolate(const std::vector<double>& x, const std::vector<double>& err) {
        xs_.push_back(x);
        es_.push_back(err);
        if (static_cast<int>(xs_.size()) > size_) {
            xs_.erase(xs_.begin());
            es_.erase(es_.begin());
        }
        const std::size_t m = xs_.size();
        if (m < 2) return x;
        const std::size_t N = m + 1;
        std::vector<double> B(N * N, 0.0), rhs(N, 0.0);
        double scale = 0.0;
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j <= i; ++j) {
                double v = 0.0;
                for (std::size_t k = 0; k < es_[i].size(); ++k) v += es_[i][k] * es_[j][k];
                B[i * N + j] = B[j * N + i] = v;
                scale = std::max(scale, v);
            }
        if (scale <= 0.0) return x;
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) B[i * N + j] /= scale;
        for (std::size_t i = 0; i < m; ++i) B[i * N + m] = B[m * N + i] = -1.0;
        rhs[m] = -1.0;
        // Gaussian elimination
        std::vector<double> a = B, b = rhs;
        for (std::size_t c = 0; c < N; ++c) {
            std::size_t piv = c;
            for (std::size_t r = c + 1; r < N; ++r)
                if (std::fabs(a[r * N + c]) > std::fabs(a[piv * N + c])) piv = r;
            if (std::fabs(a[piv * N + c]) < 1e-14) {
                xs_.erase(xs_.begin());
                es_.erase(es_.begin());
                return x;
            }
            if (piv != c) {
                for (std::size_t q = 0; q < N; ++q) std::swap(a[piv * N + q], a[c * N + q]);
                std::swap(b[piv], b[c]);
            }
            for (std::size_t r = c + 1; r < N; ++r) {
                const double f = a[r * N + c] / a[c * N + c];
                for (std::size_t q = c; q < N; ++q) a[r * N + q] -= f * a[c * N + q];
                b[r] -= f * b[c];
            }
        }
        std::vector<double> cvec(N, 0.0);
        for (std::size_t r = N; r-- > 0;) {
            double sum = b[r];
            for (std::size_t q = r + 1; q < N; ++q) sum -= a[r * N + q] * cvec[q];
            cvec[r] = sum / a[r * N + r];
        }
        std::vector<double> out(x.size(), 0.0);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t k = 0; k < out.size(); ++k) out[k] += cvec[i] * xs_[i][k];
        return out;
    }

private:
    int size_;
    std::vector<std::vector<double>> xs_, es_;
};

std::vector<double> diagonal_of(const BlockTensor& f, int spin, const Spaces& sp) {
    const Tensor* t = f.find(spin ? 3u : 0u);
    const std::size_t d = sp.dim(f.modes()[0], spin);
    std::vector<double> e(d, 0.0);
    if (t)
        for (std::size_t i = 0; i < d; ++i) e[i] = (*t)[i * d + i];
    return e;
}

double sorted_norm(const BlockTensor& t, int m) {
    const double f = factorial(m);
    return std::sqrt(t.dot(t)) / f;
}

}  // namespace

// ------------------------------------------------------------------------------------------------

SOSize estimate_so_method(const std::string& method, const int nocc[2], const int nvir[2]) {
    SOSize s;
    MethodInfo mi;
    if (!method_info(method, mi)) return s;
    Spaces sp;
    for (int k = 0; k < 2; ++k) {
        sp.o[k] = nocc[k];
        sp.v[k] = nvir[k];
    }
    auto amp_size = [&](int m) {
        double n = 0.0;
        for (unsigned k : amplitude_keys(m)) {
            double b = 1.0;
            for (std::size_t d : so::block_dims(amp_modes(m), k, sp)) b *= static_cast<double>(d);
            n += b;
        }
        return n;
    };
    std::vector<int> ranks;
    if (mi.family == "MP") {
        const int n = mi.rank;
        ranks = {1, 2};
        if (n >= 3) ranks.insert(ranks.end(), {1, 2});
        if (n >= 4) ranks.insert(ranks.end(), {1, 2, 3});
        if (n >= 5) ranks.insert(ranks.end(), {1, 2, 3, 4});
    } else {
        for (int m = 1; m <= mi.rank; ++m) ranks.push_back(m);
        if (mi.triples) ranks.push_back(3);
    }
    double amps = 0.0, largest = 0.0;
    for (int r : ranks) {
        const double a = amp_size(r);
        amps += a;
        largest = std::max(largest, a);
    }
    s.amplitudes = amps;
    const double o = nocc[0] + nocc[1], v = nvir[0] + nvir[1];
    // Hamiltonian blocks (dense spin blocks, about 3/8 of v^4 for vvvv) and the stored amplitudes,
    // residuals, multipliers, DIIS vectors and the largest intermediates (a few times the largest
    // amplitude tensor, or o^2 v^2 x v^2 / o^2 for the ladder)
    const double vvvv = 0.375 * v * v * v * v, other = (o + v) * (o + v) * (o + v) * (o + v) / 16.0 * 6.0;
    s.memory_bytes = 8.0 * (vvvv + other + 14.0 * amps + 6.0 * largest);
    const int top = ranks.empty() ? 2 : *std::max_element(ranks.begin(), ranks.end());
    s.flops_per_iteration = std::pow(o, top) * std::pow(v, top + 2) * 0.1;
    return s;
}

SOResult run_so_method(const std::string& method, const SpatialHamiltonian& H, const std::string& reference,
                       const SOOptions& opt, const SOProgress& progress) {
    const auto t0 = Clock::now();
    auto st = std::make_shared<SOMethodState>();
    SOMethodState& s = *st;
    if (!method_info(method, s.info)) throw std::invalid_argument("the spin-orbital engine has no method " + method);
    s.method = method;
    s.reference = reference;
    for (int k = 0; k < 2; ++k) {
        s.n[k] = H.n[k];
        s.nocc[k] = H.nocc[k];
        s.sp.o[k] = H.nocc[k];
        s.sp.v[k] = H.n[k] - H.nocc[k];
    }
    if (s.sp.o[0] + s.sp.o[1] == 0) throw std::runtime_error(method + " needs electrons");
    set_tensor_threads(opt.threads);
    build_hamiltonian_nodes(s);
    if (s.info.family == "MP") build_mp(s);
    else build_cluster_method(s);
    s.eng = std::make_unique<so::Engine>(s.g, s.sp);
    report_progress(JobStage::Transformation);
    set_hamiltonian(s, H.h, [&H](int a, int b, std::size_t p, std::size_t q, std::size_t r, std::size_t x) {
        return H.g(a, b, p, q, r, x);
    });
    // Fock diagonal for the denominators
    s.eng->forward({s.fock[0], s.fock[3], s.eref});
    for (int k = 0; k < 2; ++k) {
        s.eo[k] = diagonal_of(s.eng->value(s.fock[0]), k, s.sp);
        s.ev[k] = diagonal_of(s.eng->value(s.fock[3]), k, s.sp);
    }
    SOResult res;
    res.method = method;
    res.reference = reference;
    res.reference_energy = s.eng->value(s.eref).scalar() + H.nuclear_repulsion;
    for (Var& v : s.vars) {
        s.t.push_back(zero_amplitude(v.rank, s.sp));
        s.lambda.push_back(zero_amplitude(v.rank, s.sp));
        s.eng->set_input(v.input, s.t.back());
        s.eng->set_input(v.lambda, s.lambda.back());
    }
    report_progress(JobStage::CoupledCluster);
    double e_last = 0.0;
    bool all_converged = true;
    for (int stage = 0; stage < s.nstages; ++stage) {
        std::vector<std::size_t> vi;
        for (std::size_t k = 0; k < s.vars.size(); ++k)
            if (s.vars[k].stage == stage) vi.push_back(k);
        if (vi.empty()) continue;
        std::vector<int> outs, volatile_in;
        for (std::size_t k : vi) {
            outs.push_back(s.vars[k].residual);
            volatile_in.push_back(s.vars[k].input);
        }
        const int enode = s.stage_energy[static_cast<std::size_t>(std::min<int>(stage, static_cast<int>(s.stage_energy.size()) - 1))];
        outs.push_back(enode);
        if (s.ec >= 0) outs.push_back(s.ec);
        double vec = 0.0;
        for (std::size_t k : vi) vec += amplitude_count(s.vars[k].rank, s.sp);
        const int dsize = std::min(opt.diis_size, static_cast<int>(0.25 * opt.memory_mb * 1048576.0 / (16.0 * vec)));
        Diis diis(dsize);
        double e_prev = 0.0;
        bool converged = false;
        std::string label = s.info.family == "MP" ? "MP order " + std::to_string(stage + 1) + " amplitudes"
                                                  : (stage == 0 ? "amplitudes" : "(T) amplitudes");
        for (int it = 1; it <= opt.max_iter; ++it) {
            check_cancelled();
            for (std::size_t k : vi) s.eng->set_input(s.vars[k].input, s.t[k]);
            s.eng->forward(outs, false, volatile_in);
            const double e = s.eng->value(enode).scalar() + H.nuclear_repulsion;
            const double shift = s.ec >= 0 ? s.eng->value(s.ec).scalar() : 0.0;
            double dmax = 0.0;
            std::vector<BlockTensor> steps(vi.size());
            for (std::size_t q = 0; q < vi.size(); ++q) {
                const std::size_t k = vi[q];
                dmax = std::max(dmax, jacobi_step(s.eng->value(s.vars[k].residual), s.vars[k].rank, s.sp, s.eo, s.ev,
                                                  s.vars[k].ci ? shift : 0.0, steps[q]));
            }
            SOIteration rec;
            rec.iteration = it;
            rec.energy = e;
            rec.delta_e = it == 1 ? 0.0 : e - e_prev;
            rec.residual = dmax;
            rec.stage = label;
            res.history.push_back(rec);
            res.iterations += 1;
            report_progress(JobStage::CoupledCluster, it, e, dmax);
            if (progress) progress(rec);
            e_last = e;
            if (dmax < opt.amplitude_tol && (it > 1 || dmax == 0.0) && std::fabs(e - e_prev) < opt.energy_tol) {
                converged = true;
                break;
            }
            e_prev = e;
            std::vector<double> x, err;
            std::vector<const BlockTensor*> cur;
            std::vector<BlockTensor*> curw;
            for (std::size_t q = 0; q < vi.size(); ++q) {
                s.t[vi[q]].axpy(1.0, steps[q]);
                cur.push_back(&s.t[vi[q]]);
                curw.push_back(&s.t[vi[q]]);
            }
            std::vector<const BlockTensor*> sc;
            for (auto& b : steps) sc.push_back(&b);
            if (dsize >= 2) {
                flatten(cur, x);
                flatten(sc, err);
                unflatten(diis.extrapolate(x, err), curw);
            }
        }
        if (!converged) all_converged = false;
        for (std::size_t k : vi) s.eng->set_input(s.vars[k].input, s.t[k]);
    }
    // final energy and its components
    std::vector<int> outs = {s.energy, s.eref};
    for (const auto& et : s.energy_terms) outs.push_back(et.second);
    s.eng->forward(outs, false);
    res.total_energy = s.eng->value(s.energy).scalar() + H.nuclear_repulsion;
    res.correlation_energy = res.total_energy - res.reference_energy;
    for (const auto& et : s.energy_terms) res.terms.push_back({et.first, s.eng->value(et.second).scalar()});
    res.converged = all_converged;
    (void)e_last;
    for (std::size_t k = 0; k < s.vars.size(); ++k) {
        SOAmplitudeInfo a;
        a.name = s.vars[k].name;
        a.rank = s.vars[k].rank;
        a.norm = sorted_norm(s.t[k], a.rank);
        double cnt = 0.0;
        for (const auto& kv : s.t[k].blocks()) cnt += static_cast<double>(kv.second.size());
        a.count = cnt / (factorial(a.rank) * factorial(a.rank));
        res.amplitudes.push_back(a);
        if (a.rank == 1 && s.vars[k].order <= 1) {
            const double nel = s.sp.o[0] + s.sp.o[1];
            res.t1_diagnostic = std::sqrt(s.t[k].dot(s.t[k]) / nel);
        }
    }
    res.seconds = std::chrono::duration<double>(Clock::now() - t0).count();
    res.state = st;
    return res;
}

// ---- Lambda equations and densities -------------------------------------------------------------

namespace {

// gradient with respect to the unique amplitudes: m!^2 times the antisymmetric part
BlockTensor sorted_gradient(const BlockTensor& g, int m) {
    if (g.empty()) return g;
    BlockTensor x = so::antisymmetrize(g);
    const double f = factorial(m);
    x.scale(f * f);
    return x;
}

// Adjoints of the method's graph for the Lagrangian with the current multipliers: seeds E = 1 and
// R_v = lambda_v / m!^2 for the variables in `which` (all if empty).
std::vector<std::pair<int, BlockTensor>> lagrangian_seeds(const SOMethodState& s, const std::vector<std::size_t>& which,
                                                         bool with_energy) {
    std::vector<std::pair<int, BlockTensor>> seeds;
    if (with_energy) seeds.emplace_back(s.energy, BlockTensor::make_scalar(1.0));
    for (std::size_t k : which) {
        BlockTensor l = s.lambda[k];
        const double f = factorial(s.vars[k].rank);
        l.scale(1.0 / (f * f));
        seeds.emplace_back(s.vars[k].residual, std::move(l));
    }
    return seeds;
}

// Spatial densities from the adjoints of the Hamiltonian inputs.
void collect_densities(const SOMethodState& s, const std::function<const BlockTensor&(int)>& adj, SODensities& d) {
    for (int k = 0; k < 2; ++k) {
        d.n[k] = s.n[k];
        d.D[k] = Matrix(static_cast<std::size_t>(s.n[k]), static_cast<std::size_t>(s.n[k]));
    }
    const std::size_t na = static_cast<std::size_t>(s.n[0]), nb = static_cast<std::size_t>(s.n[1]);
    d.Gaa.assign(na * na * na * na, 0.0);
    d.Gbb.assign(nb * nb * nb * nb, 0.0);
    d.Gab.assign(na * na * nb * nb, 0.0);
    for (int k = 0; k < 4; ++k) {
        const BlockTensor& a = adj(s.hin[k]);
        const std::string ul = kHBlocks[k];
        for (const auto& kv : a.blocks()) {
            const int spin = kv.first == 0u ? 0 : (kv.first == 3u ? 1 : -1);
            if (spin < 0) continue;
            const std::size_t d1 = s.sp.dim(ul[1], spin);
            for (std::size_t p = 0; p < s.sp.dim(ul[0], spin); ++p)
                for (std::size_t q = 0; q < d1; ++q)
                    d.D[spin](spatial_index(ul[0], spin, p, s.nocc), spatial_index(ul[1], spin, q, s.nocc)) += kv.second[p * d1 + q];
        }
    }
    auto add = [&](int s1, int s2, std::size_t p, std::size_t q, std::size_t r, std::size_t x, double v) {
        if (s1 == 0 && s2 == 0) d.Gaa[((p * na + q) * na + r) * na + x] += v;
        else if (s1 == 1 && s2 == 1) d.Gbb[((p * nb + q) * nb + r) * nb + x] += v;
        else if (s1 == 0) d.Gab[((p * na + q) * nb + r) * nb + x] += v;
        else d.Gab[((r * na + x) * nb + p) * nb + q] += v;
    };
    for (const auto& vk : s.vin) {
        const BlockTensor& a = adj(vk.second);
        const std::string modes = vk.first.substr(0, 2) + vk.first.substr(3);
        for (const auto& kv : a.blocks()) {
            const unsigned key = kv.first;
            const int s0 = key & 1u, s1 = (key >> 1) & 1u, s2 = (key >> 2) & 1u, s3 = (key >> 3) & 1u;
            const bool coul = s0 == s2 && s1 == s3, exch = s0 == s3 && s1 == s2;
            const std::vector<std::size_t>& dims = kv.second.dims();
            std::size_t pos = 0;
            for (std::size_t P = 0; P < dims[0]; ++P) {
                const std::size_t p = spatial_index(modes[0], s0, P, s.nocc);
                for (std::size_t Q = 0; Q < dims[1]; ++Q) {
                    const std::size_t q = spatial_index(modes[1], s1, Q, s.nocc);
                    for (std::size_t R = 0; R < dims[2]; ++R) {
                        const std::size_t r = spatial_index(modes[2], s2, R, s.nocc);
                        for (std::size_t S = 0; S < dims[3]; ++S, ++pos) {
                            const double v = kv.second[pos];
                            if (v == 0.0) continue;
                            const std::size_t x = spatial_index(modes[3], s3, S, s.nocc);
                            if (coul) add(s0, s1, p, r, q, x, v);
                            if (exch) add(s0, s1, p, x, q, r, -v);
                        }
                    }
                }
            }
        }
    }
}

}  // namespace

SODensities so_densities(SOResult& result, const SOOptions& opt, const SOProgress& progress) {
    SOMethodState& s = *result.state;
    std::vector<std::size_t> all;
    for (std::size_t k = 0; k < s.vars.size(); ++k) all.push_back(k);
    for (std::size_t k = 0; k < s.vars.size(); ++k) s.eng->set_input(s.vars[k].input, s.t[k]);
    // values of everything (residuals and energy at the solution)
    std::vector<int> outs = {s.energy};
    for (const Var& v : s.vars) outs.push_back(v.residual);
    if (s.ec >= 0) outs.push_back(s.ec);
    s.eng->forward(outs);
    const double shift = s.ec >= 0 ? s.eng->value(s.ec).scalar() : 0.0;
    SODensities d;
    d.lambda_converged = true;
    report_progress(JobStage::Lambda);
    for (int stage = s.nstages - 1; stage >= 0; --stage) {
        std::vector<std::size_t> vi, later;
        for (std::size_t k = 0; k < s.vars.size(); ++k) {
            if (s.vars[k].stage == stage) vi.push_back(k);
            else if (s.vars[k].stage > stage) later.push_back(k);
        }
        if (vi.empty()) continue;
        // constant part: E and the multipliers of later stages
        s.eng->reverse(lagrangian_seeds(s, later, true));
        std::vector<BlockTensor> base(vi.size());
        for (std::size_t q = 0; q < vi.size(); ++q) base[q] = s.eng->adjoint(s.vars[vi[q]].input);
        double vec = 0.0;
        for (std::size_t k : vi) vec += amplitude_count(s.vars[k].rank, s.sp);
        const int dsize = std::min(opt.diis_size, static_cast<int>(0.25 * opt.memory_mb * 1048576.0 / (16.0 * vec)));
        Diis diis(dsize);
        bool converged = false;
        double worst = 0.0;
        for (int it = 1; it <= opt.max_iter; ++it) {
            check_cancelled();
            s.eng->reverse(lagrangian_seeds(s, vi, false));
            worst = 0.0;
            std::vector<BlockTensor> steps(vi.size());
            for (std::size_t q = 0; q < vi.size(); ++q) {
                const std::size_t k = vi[q];
                BlockTensor g = s.eng->adjoint(s.vars[k].input);
                g.axpy(1.0, base[q]);
                const BlockTensor G = sorted_gradient(g, s.vars[k].rank);
                worst = std::max(worst, jacobi_step(G, s.vars[k].rank, s.sp, s.eo, s.ev, s.vars[k].ci ? shift : 0.0, steps[q]));
            }
            SOIteration rec;
            rec.iteration = it;
            rec.residual = worst;
            rec.stage = "Lambda";
            report_progress(JobStage::Lambda, it, 0.0, worst);
            if (progress) progress(rec);
            d.lambda_iterations += 1;
            if (worst < opt.amplitude_tol && it > 1) {
                converged = true;
                break;
            }
            std::vector<double> x, err;
            std::vector<const BlockTensor*> cur, sc;
            std::vector<BlockTensor*> curw;
            for (std::size_t q = 0; q < vi.size(); ++q) {
                s.lambda[vi[q]].axpy(1.0, steps[q]);
                cur.push_back(&s.lambda[vi[q]]);
                curw.push_back(&s.lambda[vi[q]]);
                sc.push_back(&steps[q]);
            }
            if (dsize >= 2) {
                flatten(cur, x);
                flatten(sc, err);
                unflatten(diis.extrapolate(x, err), curw);
            }
        }
        d.lambda_residual = std::max(d.lambda_residual, worst);
        if (!converged) d.lambda_converged = false;
    }
    for (std::size_t k = 0; k < s.vars.size(); ++k) s.eng->set_input(s.vars[k].lambda, s.lambda[k]);
    s.lambda_done = true;
    // densities: adjoints of the Hamiltonian for the Lagrangian E + sum lambda R
    s.eng->reverse(lagrangian_seeds(s, all, true));
    collect_densities(s, [&](int id) -> const BlockTensor& { return s.eng->adjoint(id); }, d);
    d.energy = s.eng->value(s.energy).scalar();
    s.dens = d;
    return d;
}

CorrelatedDensity so_closed_shell_density(const SODensities& d) {
    const std::size_t n = static_cast<std::size_t>(d.n[0]);
    CorrelatedDensity c;
    c.D = Matrix(n, n);
    for (std::size_t p = 0; p < n; ++p)
        for (std::size_t q = 0; q < n; ++q) c.D(p, q) = d.D[0](p, q) + d.D[1](p, q);
    c.G = Packed4(n);
    for (std::size_t p = 0; p < n; ++p)
        for (std::size_t q = 0; q < n; ++q)
            for (std::size_t r = 0; r < n; ++r)
                for (std::size_t x = 0; x < n; ++x) {
                    const std::size_t i = ((p * n + q) * n + r) * n + x;
                    const double v = d.Gaa[i] + d.Gbb[i] + d.Gab[i];
                    if (v != 0.0) c.G.add_raw(p, q, r, x, v);
                }
    return c;
}

SOSecondOrder so_second_order(SOResult& result, const std::vector<HamiltonianDerivative>& pert, const SOOptions& opt,
                              const SOProgress& progress) {
    SOMethodState& s = *result.state;
    if (!s.lambda_done) throw std::logic_error("so_second_order: the Lambda equations have not been solved");
    SOSecondOrder out;
    const std::size_t np = pert.size(), nv = s.vars.size();
    for (std::size_t k = 0; k < nv; ++k) {
        s.eng->set_input(s.vars[k].input, s.t[k]);
        s.eng->set_input(s.vars[k].lambda, s.lambda[k]);
    }
    std::vector<int> outs = {s.energy};
    for (const Var& v : s.vars) outs.push_back(v.residual);
    if (s.ec >= 0) outs.push_back(s.ec);
    s.eng->forward(outs, true);
    const double shift = s.ec >= 0 ? s.eng->value(s.ec).scalar() : 0.0;
    std::vector<std::vector<BlockTensor>> tx(np, std::vector<BlockTensor>(nv)), ltt(np, std::vector<BlockTensor>(nv));
    out.dD.resize(np);
    report_progress(JobStage::Hessian, 3);
    for (std::size_t x = 0; x < np; ++x) {
        // the perturbed Hamiltonian as input tangents
        std::vector<std::pair<int, BlockTensor>> base;
        for (int k = 0; k < 4; ++k) base.emplace_back(s.hin[k], make_h_block(kHBlocks[k], pert[x].h, s.sp, s.nocc));
        if (pert[x].g)
            for (const auto& kv : s.vin)
                base.emplace_back(kv.second, make_v_block(kv.first.substr(0, 2) + kv.first.substr(3), pert[x].g, s.sp, s.nocc));
        for (std::size_t k = 0; k < nv; ++k) tx[x][k] = zero_amplitude(s.vars[k].rank, s.sp);
        for (int stage = 0; stage < s.nstages; ++stage) {
            std::vector<std::size_t> vi;
            for (std::size_t k = 0; k < nv; ++k)
                if (s.vars[k].stage == stage) vi.push_back(k);
            if (vi.empty()) continue;
            std::vector<int> routs;
            for (std::size_t k : vi) routs.push_back(s.vars[k].residual);
            double vec = 0.0;
            for (std::size_t k : vi) vec += amplitude_count(s.vars[k].rank, s.sp);
            const int dsize = std::min(opt.diis_size, static_cast<int>(0.25 * opt.memory_mb * 1048576.0 / (16.0 * vec)));
            Diis diis(dsize);
            bool converged = false;
            double worst = 0.0;
            for (int it = 1; it <= opt.max_iter; ++it) {
                check_cancelled();
                std::vector<std::pair<int, BlockTensor>> tin = base;
                for (std::size_t k = 0; k < nv; ++k)
                    if (s.vars[k].stage <= stage) tin.emplace_back(s.vars[k].input, tx[x][k]);
                s.eng->tangent(tin, routs);
                worst = 0.0;
                std::vector<BlockTensor> steps(vi.size());
                for (std::size_t q = 0; q < vi.size(); ++q) {
                    const std::size_t k = vi[q];
                    worst = std::max(worst, jacobi_step(s.eng->tangent_of(s.vars[k].residual), s.vars[k].rank, s.sp, s.eo,
                                                        s.ev, s.vars[k].ci ? shift : 0.0, steps[q]));
                }
                ++out.iterations;
                SOIteration rec;
                rec.iteration = it;
                rec.residual = worst;
                rec.stage = "amplitude response " + std::to_string(x + 1) + " of " + std::to_string(np);
                report_progress(JobStage::Hessian, static_cast<int>(x + 1), 0.0, worst);
                if (progress) progress(rec);
                if (worst < opt.amplitude_tol && it > 1) {
                    converged = true;
                    break;
                }
                std::vector<double> xv, err;
                std::vector<const BlockTensor*> cur, sc;
                std::vector<BlockTensor*> curw;
                for (std::size_t q = 0; q < vi.size(); ++q) {
                    tx[x][vi[q]].axpy(1.0, steps[q]);
                    cur.push_back(&tx[x][vi[q]]);
                    curw.push_back(&tx[x][vi[q]]);
                    sc.push_back(&steps[q]);
                }
                if (dsize >= 2) {
                    flatten(cur, xv);
                    flatten(sc, err);
                    unflatten(diis.extrapolate(xv, err), curw);
                }
            }
            out.residual = std::max(out.residual, worst);
            if (!converged) out.converged = false;
        }
        // forward over reverse with the amplitude response only: dD[t^x] and L_tt t^x
        std::vector<std::pair<int, BlockTensor>> tin;
        for (std::size_t k = 0; k < nv; ++k) tin.emplace_back(s.vars[k].input, tx[x][k]);
        s.eng->tangent(tin, outs);
        std::vector<std::size_t> all;
        for (std::size_t k = 0; k < nv; ++k) all.push_back(k);
        s.eng->reverse(lagrangian_seeds(s, all, true), true);
        collect_densities(s, [&](int id) -> const BlockTensor& { return s.eng->adjoint_tangent(id); }, out.dD[x]);
        for (std::size_t k = 0; k < nv; ++k) ltt[x][k] = s.eng->adjoint_tangent(s.vars[k].input);
    }
    out.pair = Matrix(np, np);
    for (std::size_t x = 0; x < np; ++x)
        for (std::size_t y = 0; y < np; ++y) {
            double v = 0.0;
            for (std::size_t k = 0; k < nv; ++k) v += tx[y][k].dot(ltt[x][k]);
            out.pair(x, y) = v;
        }
    for (std::size_t x = 0; x < np; ++x)
        for (std::size_t y = 0; y < x; ++y) out.pair(x, y) = out.pair(y, x) = 0.5 * (out.pair(x, y) + out.pair(y, x));
    return out;
}

}  // namespace macqc
