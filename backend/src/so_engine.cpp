// ΨLab 15.0: block-sparse spin-orbital tensors and a computation graph with automatic
// differentiation (see so_engine.h).
#include "macqc/so_engine.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <stdexcept>

#include "macqc/job.h"

namespace macqc {
namespace so {

namespace {

struct ParsedSpec {
    std::string a, b, c;
};

ParsedSpec parse2(const std::string& spec) {
    const auto comma = spec.find(',');
    const auto arrow = spec.find("->");
    if (comma == std::string::npos || arrow == std::string::npos || comma > arrow)
        throw std::invalid_argument("so::einsum: bad spec '" + spec + "'");
    return {spec.substr(0, comma), spec.substr(comma + 1, arrow - comma - 1), spec.substr(arrow + 2)};
}

std::pair<std::string, std::string> parse1(const std::string& spec) {
    const auto arrow = spec.find("->");
    if (arrow == std::string::npos) throw std::invalid_argument("so::transpose: bad spec '" + spec + "'");
    return {spec.substr(0, arrow), spec.substr(arrow + 2)};
}

std::string result_modes(const ParsedSpec& p, const std::string& ma, const std::string& mb) {
    std::string out;
    for (char ch : p.c) {
        const auto ia = p.a.find(ch);
        if (ia != std::string::npos) {
            out += ma[ia];
            continue;
        }
        const auto ib = p.b.find(ch);
        if (ib == std::string::npos) throw std::invalid_argument("so::einsum: result letter not in the operands");
        out += mb[ib];
    }
    return out;
}

}  // namespace

std::vector<std::size_t> block_dims(const std::string& modes, unsigned key, const Spaces& sp) {
    std::vector<std::size_t> d(modes.size());
    for (std::size_t k = 0; k < modes.size(); ++k) d[k] = sp.dim(modes[k], static_cast<int>((key >> k) & 1u));
    return d;
}

Tensor& BlockTensor::block(unsigned key, const Spaces& sp) {
    auto it = blocks_.find(key);
    if (it != blocks_.end()) return it->second;
    return blocks_.emplace(key, Tensor(block_dims(modes_, key, sp), 0.0)).first->second;
}

const Tensor* BlockTensor::find(unsigned key) const {
    const auto it = blocks_.find(key);
    return it == blocks_.end() ? nullptr : &it->second;
}

void BlockTensor::axpy(double a, const BlockTensor& x) {
    if (x.empty() || a == 0.0) return;
    if (modes_.empty() && blocks_.empty()) modes_ = x.modes_;
    if (modes_ != x.modes_) throw std::logic_error("BlockTensor::axpy: modes differ");
    for (const auto& kv : x.blocks_) {
        auto it = blocks_.find(kv.first);
        if (it == blocks_.end()) {
            Tensor t = kv.second;
            if (a != 1.0) t *= a;
            blocks_.emplace(kv.first, std::move(t));
        } else {
            it->second.axpy(a, kv.second);
        }
    }
}

void BlockTensor::scale(double a) {
    for (auto& kv : blocks_) kv.second *= a;
}

double BlockTensor::dot(const BlockTensor& x) const {
    double s = 0.0;
    for (const auto& kv : blocks_) {
        const Tensor* o = x.find(kv.first);
        if (o) s += kv.second.dot(*o);
    }
    return s;
}

double BlockTensor::max_abs() const {
    double m = 0.0;
    for (const auto& kv : blocks_)
        for (double v : kv.second.vec()) m = std::max(m, std::fabs(v));
    return m;
}

std::size_t BlockTensor::size() const {
    std::size_t n = 0;
    for (const auto& kv : blocks_) n += kv.second.size();
    return n;
}

double BlockTensor::scalar() const {
    const Tensor* t = find(0u);
    return t && t->size() ? (*t)[0] : 0.0;
}

BlockTensor BlockTensor::make_scalar(double v) {
    BlockTensor s("");
    Tensor t(std::vector<std::size_t>{}, v);
    s.blocks_.emplace(0u, std::move(t));
    return s;
}

void einsum_acc(const Spaces& sp, const std::string& spec, const BlockTensor& a, const BlockTensor& b, BlockTensor& c,
                double alpha) {
    if (a.empty() || b.empty() || alpha == 0.0) return;
    const ParsedSpec p = parse2(spec);
    if (p.a.size() != a.rank() || p.b.size() != b.rank()) throw std::invalid_argument("so::einsum: rank mismatch in '" + spec + "'");
    const std::string modes = result_modes(p, a.modes(), b.modes());
    if (c.modes().empty() && c.empty()) c = BlockTensor(modes);
    if (c.modes() != modes) throw std::invalid_argument("so::einsum: accumulator modes in '" + spec + "'");
    // letter positions
    std::vector<int> b_in_a(p.b.size(), -1), c_from_a(p.c.size(), -1), c_from_b(p.c.size(), -1);
    for (std::size_t j = 0; j < p.b.size(); ++j) {
        const auto ia = p.a.find(p.b[j]);
        if (ia != std::string::npos) b_in_a[j] = static_cast<int>(ia);
    }
    for (std::size_t k = 0; k < p.c.size(); ++k) {
        const auto ia = p.a.find(p.c[k]);
        if (ia != std::string::npos) c_from_a[k] = static_cast<int>(ia);
        else c_from_b[k] = static_cast<int>(p.b.find(p.c[k]));
    }
    for (const auto& ka : a.blocks()) {
        for (const auto& kb : b.blocks()) {
            bool ok = true;
            for (std::size_t j = 0; j < p.b.size() && ok; ++j)
                if (b_in_a[j] >= 0 && (((ka.first >> b_in_a[j]) & 1u) != ((kb.first >> j) & 1u))) ok = false;
            if (!ok) continue;
            unsigned key = 0;
            for (std::size_t k = 0; k < p.c.size(); ++k) {
                const unsigned bit = c_from_a[k] >= 0 ? (ka.first >> c_from_a[k]) & 1u : (kb.first >> c_from_b[k]) & 1u;
                key |= bit << k;
            }
            if (ka.second.size() == 0 || kb.second.size() == 0) continue;
            Tensor& out = c.block(key, sp);
            if (out.size() == 0) continue;
            macqc::einsum_acc(spec, ka.second, kb.second, out, alpha);
        }
    }
}

void transpose_acc(const BlockTensor& a, const std::string& from, const std::string& to, BlockTensor& c, double alpha) {
    if (a.empty() || alpha == 0.0) return;
    if (from.size() != a.rank() || to.size() != a.rank()) throw std::invalid_argument("so::transpose: rank mismatch");
    std::vector<int> src(to.size());
    std::string modes(to.size(), ' ');
    for (std::size_t k = 0; k < to.size(); ++k) {
        const auto i = from.find(to[k]);
        if (i == std::string::npos) throw std::invalid_argument("so::transpose: letters differ");
        src[k] = static_cast<int>(i);
        modes[k] = a.modes()[i];
    }
    if (c.modes().empty() && c.empty()) c = BlockTensor(modes);
    if (c.modes() != modes) throw std::invalid_argument("so::transpose: accumulator modes");
    const std::string spec = from + "->" + to;
    for (const auto& kv : a.blocks()) {
        unsigned key = 0;
        for (std::size_t k = 0; k < to.size(); ++k) key |= ((kv.first >> src[k]) & 1u) << k;
        Tensor t = macqc::einsum(spec, kv.second, alpha);
        auto it = c.blocks().find(key);
        if (it == c.blocks().end()) c.blocks().emplace(key, std::move(t));
        else it->second.axpy(1.0, t);
    }
}

// ---- graph ----------------------------------------------------------------------------------

int Graph::add(Node n, const std::string& key) {
    if (!key.empty()) {
        const auto it = cse_.find(key);
        if (it != cse_.end()) return it->second;
    }
    nodes_.push_back(std::move(n));
    const int id = static_cast<int>(nodes_.size()) - 1;
    if (!key.empty()) cse_.emplace(key, id);
    return id;
}

int Graph::input(const std::string& name, const std::string& modes, bool variable) {
    const auto it = inputs_.find(name);
    if (it != inputs_.end()) return it->second;
    Node n;
    n.op = Op::Input;
    n.name = name;
    n.modes = modes;
    n.variable = variable;
    const int id = add(std::move(n), "");
    inputs_.emplace(name, id);
    return id;
}

int Graph::find_input(const std::string& name) const {
    const auto it = inputs_.find(name);
    return it == inputs_.end() ? -1 : it->second;
}

static std::string coef_key(double c) {
    char buf[40];
    std::snprintf(buf, sizeof buf, "%.17g", c);
    return buf;
}

int Graph::einsum(const std::string& spec, int a, int b, double coef) {
    const ParsedSpec p = parse2(spec);
    Node n;
    n.op = Op::Einsum;
    n.spec = spec;
    n.args = {a, b};
    n.coefs = {coef};
    n.modes = result_modes(p, node(a).modes, node(b).modes);
    if (p.a.size() != node(a).modes.size() || p.b.size() != node(b).modes.size())
        throw std::invalid_argument("Graph::einsum: rank mismatch in '" + spec + "'");
    return add(std::move(n), "E|" + spec + "|" + std::to_string(a) + "|" + std::to_string(b) + "|" + coef_key(coef));
}

int Graph::transpose(int a, const std::string& from, const std::string& to, double coef) {
    if (from == to && coef == 1.0) return a;
    Node n;
    n.op = Op::Transpose;
    n.spec = from + "->" + to;
    n.args = {a};
    n.coefs = {coef};
    const std::string& ma = node(a).modes;
    n.modes.assign(to.size(), ' ');
    for (std::size_t k = 0; k < to.size(); ++k) {
        const auto i = from.find(to[k]);
        if (i == std::string::npos || from.size() != ma.size()) throw std::invalid_argument("Graph::transpose: letters");
        n.modes[k] = ma[i];
    }
    return add(std::move(n), "T|" + n.spec + "|" + std::to_string(a) + "|" + coef_key(coef));
}

int Graph::sum(const std::vector<std::pair<double, int>>& terms) {
    std::vector<std::pair<double, int>> t;
    for (const auto& x : terms)
        if (x.second >= 0 && x.first != 0.0) t.push_back(x);
    if (t.empty()) return -1;
    if (t.size() == 1 && t[0].first == 1.0) return t[0].second;
    Node n;
    n.op = Op::Sum;
    n.modes = node(t[0].second).modes;
    std::string key = "S";
    for (const auto& x : t) {
        if (node(x.second).modes != n.modes) throw std::invalid_argument("Graph::sum: modes differ");
        n.args.push_back(x.second);
        n.coefs.push_back(x.first);
        key += "|" + std::to_string(x.second) + ":" + coef_key(x.first);
    }
    return add(std::move(n), key);
}

int Graph::scale(int scalar, int a, double coef) {
    if (!node(scalar).modes.empty()) throw std::invalid_argument("Graph::scale: not a scalar");
    Node n;
    n.op = Op::Scale;
    n.args = {scalar, a};
    n.coefs = {coef};
    n.modes = node(a).modes;
    return add(std::move(n), "X|" + std::to_string(scalar) + "|" + std::to_string(a) + "|" + coef_key(coef));
}

int Graph::permsum(int a, const std::vector<std::pair<std::string, std::string>>& perms, const std::vector<double>& coefs) {
    if (perms.empty()) return -1;
    if (perms.size() == 1) return transpose(a, perms[0].first, perms[0].second, coefs[0]);
    Node n;
    n.op = Op::Permsum;
    n.args = {a};
    n.coefs = coefs;
    std::string key = "P|" + std::to_string(a);
    const std::string& ma = node(a).modes;
    for (std::size_t k = 0; k < perms.size(); ++k) {
        n.spec += (k ? ";" : "") + perms[k].first + "->" + perms[k].second;
        key += "|" + perms[k].first + ">" + perms[k].second + ":" + coef_key(coefs[k]);
    }
    const std::string& to = perms[0].second;
    const std::string& from = perms[0].first;
    n.modes.assign(to.size(), ' ');
    for (std::size_t k = 0; k < to.size(); ++k) n.modes[k] = ma[from.find(to[k])];
    return add(std::move(n), key);
}

namespace {
std::vector<std::pair<std::string, std::string>> split_perms(const std::string& spec) {
    std::vector<std::pair<std::string, std::string>> out;
    std::size_t start = 0;
    while (start <= spec.size()) {
        const auto semi = spec.find(';', start);
        const std::string one = spec.substr(start, semi == std::string::npos ? std::string::npos : semi - start);
        out.push_back(parse1(one));
        if (semi == std::string::npos) break;
        start = semi + 1;
    }
    return out;
}
}  // namespace

// ---- engine -----------------------------------------------------------------------------------

Engine::Engine(const Graph& g, const Spaces& sp)
    : g_(g), sp_(sp), val_(g.size()), tan_(g.size()), adj_(g.size()), adt_(g.size()), consumers_(g.size()),
      operand_(g.size(), 0), anc_inputs_(g.size()), stamp_in_(g.size(), 0), stamp_node_(g.size(), -1),
      valid_(g.size(), 0) {
    for (std::size_t i = 0; i < g.size(); ++i) {
        const Node& n = g.node(static_cast<int>(i));
        for (int a : n.args) consumers_[static_cast<std::size_t>(a)].push_back(static_cast<int>(i));
        if (n.op == Op::Einsum || n.op == Op::Scale)
            for (int a : n.args) operand_[static_cast<std::size_t>(a)] = 1;
        if (n.op == Op::Input) {
            anc_inputs_[i] = {static_cast<int>(i)};
        } else {
            std::vector<int> u;
            for (int a : n.args) u.insert(u.end(), anc_inputs_[static_cast<std::size_t>(a)].begin(), anc_inputs_[static_cast<std::size_t>(a)].end());
            std::sort(u.begin(), u.end());
            u.erase(std::unique(u.begin(), u.end()), u.end());
            anc_inputs_[i] = std::move(u);
        }
    }
}

void Engine::set_input(int id, const BlockTensor& v) {
    const Node& n = g_.node(id);
    if (n.op != Op::Input) throw std::logic_error("Engine::set_input: not an input");
    if (!v.empty() && v.modes() != n.modes) throw std::logic_error("Engine::set_input: modes of " + n.name);
    val_[static_cast<std::size_t>(id)] = v.empty() ? BlockTensor(n.modes) : v;
    stamp_in_[static_cast<std::size_t>(id)] = ++clock_;
    valid_[static_cast<std::size_t>(id)] = 1;
}

std::vector<int> Engine::needed(const std::vector<int>& outputs) const {
    std::vector<char> need(g_.size(), 0);
    std::vector<int> stack;
    for (int o : outputs)
        if (o >= 0) stack.push_back(o);
    while (!stack.empty()) {
        const int id = stack.back();
        stack.pop_back();
        if (need[static_cast<std::size_t>(id)]) continue;
        need[static_cast<std::size_t>(id)] = 1;
        for (int a : g_.node(id).args) stack.push_back(a);
    }
    std::vector<int> order;
    for (std::size_t i = 0; i < need.size(); ++i)
        if (need[i]) order.push_back(static_cast<int>(i));
    return order;
}

bool Engine::up_to_date(int id) const {
    const std::size_t i = static_cast<std::size_t>(id);
    if (!valid_[i]) return false;
    if (g_.node(id).op == Op::Input) return true;
    for (int a : anc_inputs_[i])
        if (stamp_in_[static_cast<std::size_t>(a)] > stamp_node_[i]) return false;
    return true;
}

// one use of the value of `a` is over: release it when nobody needs it any more
void Engine::use(int a) {
    const std::size_t i = static_cast<std::size_t>(a);
    if (--remaining_[i] > 0) return;
    if (g_.node(a).op == Op::Input || is_output_[i] || protect_[i]) return;
    if (keep_ && operand_[i]) return;
    val_[i].clear();
    valid_[i] = 0;
}

void Engine::eval(int id) {
    const std::size_t i = static_cast<std::size_t>(id);
    if (in_set_[i]) return;
    in_set_[i] = 1;
    const Node& n = g_.node(id);
    if (n.op == Op::Input) {
        if (val_[i].modes() != n.modes) val_[i] = BlockTensor(n.modes);
        return;
    }
    if (keep_)
        for (int a : n.args) eval(a);  // a later reverse pass needs the operands: (re)compute them
    if (up_to_date(id)) {
        // cached: the uses of the arguments through this node are over
        for (int a : n.args) {
            if (in_set_[static_cast<std::size_t>(a)]) use(a);
            else --remaining_[static_cast<std::size_t>(a)];
        }
        return;
    }
    check_cancelled();
    BlockTensor out(n.modes);
    switch (n.op) {
        case Op::Einsum:
            eval(n.args[0]);
            eval(n.args[1]);
            einsum_acc(sp_, n.spec, val_[static_cast<std::size_t>(n.args[0])], val_[static_cast<std::size_t>(n.args[1])], out, n.coefs[0]);
            use(n.args[0]);
            use(n.args[1]);
            break;
        case Op::Transpose: {
            eval(n.args[0]);
            const auto ft = parse1(n.spec);
            transpose_acc(val_[static_cast<std::size_t>(n.args[0])], ft.first, ft.second, out, n.coefs[0]);
            use(n.args[0]);
            break;
        }
        case Op::Sum:
            for (std::size_t k = 0; k < n.args.size(); ++k) {
                eval(n.args[k]);
                out.axpy(n.coefs[k], val_[static_cast<std::size_t>(n.args[k])]);
                use(n.args[k]);
            }
            break;
        case Op::Scale:
            eval(n.args[0]);
            eval(n.args[1]);
            out.axpy(n.coefs[0] * val_[static_cast<std::size_t>(n.args[0])].scalar(), val_[static_cast<std::size_t>(n.args[1])]);
            use(n.args[0]);
            use(n.args[1]);
            break;
        case Op::Permsum: {
            eval(n.args[0]);
            const auto perms = split_perms(n.spec);
            for (std::size_t k = 0; k < perms.size(); ++k)
                transpose_acc(val_[static_cast<std::size_t>(n.args[0])], perms[k].first, perms[k].second, out, n.coefs[k]);
            use(n.args[0]);
            break;
        }
        default: break;
    }
    val_[i] = std::move(out);
    valid_[i] = 1;
    stamp_node_[i] = ++clock_;
}

void Engine::forward(const std::vector<int>& outputs, bool keep, const std::vector<int>& volatile_inputs) {
    keep_ = keep;
    const std::vector<int> order = needed(outputs);
    const std::size_t N = g_.size();
    remaining_.assign(N, 0);
    in_set_.assign(N, 0);
    is_output_.assign(N, 0);
    protect_.assign(N, 0);
    std::vector<char> need(N, 0);
    for (int id : order) need[static_cast<std::size_t>(id)] = 1;
    for (int id : order)
        for (int a : g_.node(id).args) ++remaining_[static_cast<std::size_t>(a)];
    for (int o : outputs)
        if (o >= 0) is_output_[static_cast<std::size_t>(o)] = 1;
    if (!volatile_inputs.empty()) {
        // keep what does not depend on the inputs that change between calls
        std::vector<char> vol(N, 0);
        for (int v : volatile_inputs) vol[static_cast<std::size_t>(v)] = 1;
        std::vector<char> constant(N, 0);
        for (int id : order) {
            bool dep = false;
            for (int a : anc_inputs_[static_cast<std::size_t>(id)])
                if (vol[static_cast<std::size_t>(a)]) dep = true;
            constant[static_cast<std::size_t>(id)] = !dep;
        }
        // keep the constant nodes that changing nodes consume (the frontier), not their subtrees
        for (int id : order) {
            if (constant[static_cast<std::size_t>(id)]) continue;
            for (int a : g_.node(id).args)
                if (constant[static_cast<std::size_t>(a)]) protect_[static_cast<std::size_t>(a)] = 1;
        }
    }
    for (int o : outputs)
        if (o >= 0) eval(o);
}

void Engine::tangent(const std::vector<std::pair<int, BlockTensor>>& input_tangents, const std::vector<int>& outputs) {
    for (auto& t : tan_) t.clear();
    std::vector<char> is_set(g_.size(), 0);
    for (const auto& it : input_tangents) {
        tan_[static_cast<std::size_t>(it.first)] = it.second;
        is_set[static_cast<std::size_t>(it.first)] = 1;
    }
    const std::vector<int> order = needed(outputs);
    std::vector<int> rem(g_.size(), 0);
    std::vector<char> outp(g_.size(), 0);
    for (int id : order)
        for (int a : g_.node(id).args) ++rem[static_cast<std::size_t>(a)];
    for (int o : outputs)
        if (o >= 0) outp[static_cast<std::size_t>(o)] = 1;
    auto done_with = [&](int a) {
        const std::size_t i = static_cast<std::size_t>(a);
        if (--rem[i] > 0 || outp[i] || operand_[i] || g_.node(a).op == Op::Input) return;
        tan_[i].clear();
    };
    for (int id : order) {
        check_cancelled();
        const Node& n = g_.node(id);
        BlockTensor& out = tan_[static_cast<std::size_t>(id)];
        if (n.op == Op::Input) {
            if (!is_set[static_cast<std::size_t>(id)]) out = BlockTensor(n.modes);
            continue;
        }
        out = BlockTensor(n.modes);
        switch (n.op) {
            case Op::Einsum: {
                const int a = n.args[0], b = n.args[1];
                einsum_acc(sp_, n.spec, tan_[static_cast<std::size_t>(a)], value(b), out, n.coefs[0]);
                einsum_acc(sp_, n.spec, value(a), tan_[static_cast<std::size_t>(b)], out, n.coefs[0]);
                break;
            }
            case Op::Transpose: {
                const auto ft = parse1(n.spec);
                transpose_acc(tan_[static_cast<std::size_t>(n.args[0])], ft.first, ft.second, out, n.coefs[0]);
                break;
            }
            case Op::Sum:
                for (std::size_t k = 0; k < n.args.size(); ++k) out.axpy(n.coefs[k], tan_[static_cast<std::size_t>(n.args[k])]);
                break;
            case Op::Scale: {
                const double s = value(n.args[0]).scalar();
                const double ds = tan_[static_cast<std::size_t>(n.args[0])].scalar();
                out.axpy(n.coefs[0] * ds, value(n.args[1]));
                out.axpy(n.coefs[0] * s, tan_[static_cast<std::size_t>(n.args[1])]);
                break;
            }
            case Op::Permsum: {
                const auto perms = split_perms(n.spec);
                for (std::size_t k = 0; k < perms.size(); ++k)
                    transpose_acc(tan_[static_cast<std::size_t>(n.args[0])], perms[k].first, perms[k].second, out, n.coefs[k]);
                break;
            }
            default: break;
        }
        for (int a : n.args) done_with(a);
    }
}

void Engine::settle(int a, bool wt) {
    const std::size_t i = static_cast<std::size_t>(a);
    if (remaining_[i] <= 0) return;
    if (--remaining_[i] == 0) propagate(a, wt);
}

void Engine::propagate(int id, bool wt) {
    const std::size_t r = static_cast<std::size_t>(id);
    if (processed_[r]) return;
    processed_[r] = 1;
    const Node& n = g_.node(id);
    if (n.op == Op::Input) return;
    check_cancelled();
    const bool zero = adj_[r].empty() && adt_[r].empty();
    const double coef = n.coefs.empty() ? 1.0 : n.coefs[0];
    if (zero) {
        for (int a : n.args) settle(a, wt);
        return;
    }
    const BlockTensor cbar = std::move(adj_[r]);
    const BlockTensor cbt = std::move(adt_[r]);
    adj_[r].clear();
    adt_[r].clear();
    switch (n.op) {
        case Op::Einsum: {
            const ParsedSpec p = parse2(n.spec);
            const std::string spa = p.c + "," + p.b + "->" + p.a;
            const std::string spb = p.a + "," + p.c + "->" + p.b;
            const std::size_t a = static_cast<std::size_t>(n.args[0]), b = static_cast<std::size_t>(n.args[1]);
            einsum_acc(sp_, spa, cbar, val_[b], adj_[a], coef);
            einsum_acc(sp_, spb, val_[a], cbar, adj_[b], coef);
            if (wt) {
                einsum_acc(sp_, spa, cbt, val_[b], adt_[a], coef);
                einsum_acc(sp_, spa, cbar, tan_[b], adt_[a], coef);
                einsum_acc(sp_, spb, tan_[a], cbar, adt_[b], coef);
                einsum_acc(sp_, spb, val_[a], cbt, adt_[b], coef);
            }
            settle(n.args[0], wt);
            settle(n.args[1], wt);
            break;
        }
        case Op::Transpose: {
            const auto ft = parse1(n.spec);
            const std::size_t a = static_cast<std::size_t>(n.args[0]);
            transpose_acc(cbar, ft.second, ft.first, adj_[a], coef);
            if (wt) transpose_acc(cbt, ft.second, ft.first, adt_[a], coef);
            settle(n.args[0], wt);
            break;
        }
        case Op::Sum:
            for (std::size_t k = 0; k < n.args.size(); ++k) {
                const std::size_t a = static_cast<std::size_t>(n.args[k]);
                adj_[a].axpy(n.coefs[k], cbar);
                if (wt) adt_[a].axpy(n.coefs[k], cbt);
                settle(n.args[k], wt);  // depth first: the term is propagated before the next one
            }
            break;
        case Op::Scale: {
            const std::size_t s = static_cast<std::size_t>(n.args[0]), a = static_cast<std::size_t>(n.args[1]);
            const double sv = val_[s].scalar();
            adj_[a].axpy(coef * sv, cbar);
            const double sb = coef * val_[a].dot(cbar);
            if (sb != 0.0) adj_[s].axpy(1.0, BlockTensor::make_scalar(sb));
            if (wt) {
                const double ds = tan_[s].scalar();
                adt_[a].axpy(coef * ds, cbar);
                adt_[a].axpy(coef * sv, cbt);
                const double st = coef * (tan_[a].dot(cbar) + val_[a].dot(cbt));
                if (st != 0.0) adt_[s].axpy(1.0, BlockTensor::make_scalar(st));
            }
            settle(n.args[0], wt);
            settle(n.args[1], wt);
            break;
        }
        case Op::Permsum: {
            const auto perms = split_perms(n.spec);
            const std::size_t a = static_cast<std::size_t>(n.args[0]);
            for (std::size_t k = 0; k < perms.size(); ++k) {
                transpose_acc(cbar, perms[k].second, perms[k].first, adj_[a], n.coefs[k]);
                if (wt) transpose_acc(cbt, perms[k].second, perms[k].first, adt_[a], n.coefs[k]);
            }
            settle(n.args[0], wt);
            break;
        }
        default: break;
    }
}

void Engine::reverse(const std::vector<std::pair<int, BlockTensor>>& seeds, bool with_tangents,
                     const std::vector<std::pair<int, BlockTensor>>& seed_tangents) {
    for (auto& a : adj_) a.clear();
    for (auto& a : adt_) a.clear();
    std::vector<int> roots;
    for (const auto& s : seeds) {
        adj_[static_cast<std::size_t>(s.first)].axpy(1.0, s.second);
        roots.push_back(s.first);
    }
    if (with_tangents)
        for (const auto& s : seed_tangents) {
            adt_[static_cast<std::size_t>(s.first)].axpy(1.0, s.second);
            roots.push_back(s.first);
        }
    // pending contributions: consumers within the ancestors of the seeds
    const std::vector<int> order = needed(roots);
    const std::size_t N = g_.size();
    std::vector<char> in(N, 0);
    for (int id : order) in[static_cast<std::size_t>(id)] = 1;
    remaining_.assign(N, 0);
    processed_.assign(N, 0);
    for (int id : order)
        for (int a : g_.node(id).args) ++remaining_[static_cast<std::size_t>(a)];
    // nodes without consumers among the ancestors (the seeds, normally) start the propagation
    for (std::size_t k = order.size(); k-- > 0;) {
        const int id = order[k];
        if (remaining_[static_cast<std::size_t>(id)] == 0 && !processed_[static_cast<std::size_t>(id)]) propagate(id, with_tangents);
    }
}

void Engine::release() {
    for (std::size_t i = 0; i < g_.size(); ++i) {
        if (g_.node(static_cast<int>(i)).op == Op::Input) continue;
        val_[i].clear();
        tan_[i].clear();
        valid_[i] = 0;
    }
}

BlockTensor antisymmetrize(const BlockTensor& t) {
    const std::size_t r = t.rank();
    if (r == 0 || r % 2) return t;
    const std::size_t m = r / 2;
    const std::string base = std::string("abcdefgh").substr(0, r);
    auto pass = [&](const BlockTensor& in, std::size_t offset) {
        std::vector<int> perm(m);
        std::iota(perm.begin(), perm.end(), 0);
        BlockTensor out;
        do {
            int sign = 1;
            for (std::size_t i = 0; i < m; ++i)
                for (std::size_t j = i + 1; j < m; ++j)
                    if (perm[i] > perm[j]) sign = -sign;
            std::string to = base;
            for (std::size_t i = 0; i < m; ++i) to[offset + i] = base[offset + static_cast<std::size_t>(perm[i])];
            transpose_acc(in, base, to, out, static_cast<double>(sign));
        } while (std::next_permutation(perm.begin(), perm.end()));
        return out;
    };
    BlockTensor x = pass(pass(t, 0), m);
    double f = 1.0;
    for (std::size_t k = 2; k <= m; ++k) f *= static_cast<double>(k);
    x.scale(1.0 / (f * f));
    return x;
}

}  // namespace so
}  // namespace macqc
