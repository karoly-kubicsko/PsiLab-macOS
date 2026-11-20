// ΨLab 15.0: dense spin-orbital tensors with automatic differentiation, the engine of the
// correlated methods MP2-MP5, CISD ... CISDTQ, QCISD ... QCISDTQ and the open-shell coupled cluster.
//
// Tensors. Every index is occupied ('o') or virtual ('v') and runs over the spin orbitals of that
// space, alpha first; a BlockTensor stores one dense Tensor per spin pattern (bit k = spin of index
// k), only the patterns that are not zero (e.g. t2 has the blocks aaaa, abab, abba, baab, baba,
// bbbb). einsum contracts the blocks whose spins agree on the shared indices with the dense einsum
// of tensor.h (matrix products).
//
// Graphs. A method is written once as a graph of operations on BlockTensors: inputs (the
// Hamiltonian blocks, amplitudes, multipliers), einsum, transposition, linear combination and the
// product with a scalar node. Identical operations are created once (common subexpressions). The
// Engine evaluates a graph for given inputs and differentiates it:
//   - reverse mode (vector-Jacobian products): adjoints of all nodes for adjoints of some
//     (Lambda equations, densities as derivatives with respect to the Hamiltonian);
//   - forward mode (Jacobian-vector products): tangents of all nodes for tangents of the inputs
//     (amplitude responses);
//   - forward over reverse: the tangents of the adjoints (Hessian-vector products, density
//     responses of the analytic Hessians).
// All operations are multilinear, so every derivative is exact.
#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "macqc/tensor.h"

namespace macqc {
namespace so {

// Numbers of occupied and virtual orbitals of each spin (0 alpha, 1 beta).
struct Spaces {
    int o[2] = {0, 0};
    int v[2] = {0, 0};
    std::size_t dim(char space, int spin) const {
        return static_cast<std::size_t>(space == 'o' ? o[spin] : v[spin]);
    }
};

class BlockTensor {
public:
    BlockTensor() = default;
    explicit BlockTensor(std::string modes) : modes_(std::move(modes)) {}

    const std::string& modes() const { return modes_; }
    std::size_t rank() const { return modes_.size(); }
    std::map<unsigned, Tensor>& blocks() { return blocks_; }
    const std::map<unsigned, Tensor>& blocks() const { return blocks_; }
    bool empty() const { return blocks_.empty(); }
    // The block of a spin pattern, created (zero) if absent.
    Tensor& block(unsigned key, const Spaces& sp);
    const Tensor* find(unsigned key) const;
    void clear() { blocks_.clear(); }

    // this += a x (x must have the same modes; missing blocks are created)
    void axpy(double a, const BlockTensor& x);
    void scale(double a);
    double dot(const BlockTensor& x) const;  // sum over all elements
    double max_abs() const;
    std::size_t size() const;                // stored elements
    // Scalar value of a rank-0 tensor (0 if empty).
    double scalar() const;
    static BlockTensor make_scalar(double v);

private:
    std::string modes_;
    std::map<unsigned, Tensor> blocks_;
};

// Dimensions of the dense block of a spin pattern.
std::vector<std::size_t> block_dims(const std::string& modes, unsigned key, const Spaces& sp);

// c += alpha * einsum(spec, a, b) over the spin blocks ("ijab,abcd->ijcd"; every letter of an
// operand must occur in the other operand or in the result; no repeated letters).
void einsum_acc(const Spaces& sp, const std::string& spec, const BlockTensor& a, const BlockTensor& b, BlockTensor& c,
                double alpha = 1.0);
// c += alpha * a relabelled: the index with letter from[k] of a becomes the index with the same
// letter in `to` (a transposition).
void transpose_acc(const BlockTensor& a, const std::string& from, const std::string& to, BlockTensor& c,
                   double alpha = 1.0);

// ---- graphs -----------------------------------------------------------------------------

enum class Op { Input, Einsum, Transpose, Sum, Scale, Permsum };

struct Node {
    Op op = Op::Input;
    std::string modes;
    std::string name;            // inputs
    bool variable = true;        // inputs: differentiated (false: constant, e.g. the identity)
    std::string spec;            // Einsum: "ab,bc->ac"; Transpose: "from->to"
    std::vector<int> args;
    std::vector<double> coefs;   // Sum: one per argument; Einsum, Transpose, Scale: coefs[0]
};

class Graph {
public:
    int input(const std::string& name, const std::string& modes, bool variable = true);
    int einsum(const std::string& spec, int a, int b, double coef = 1.0);
    int transpose(int a, const std::string& from, const std::string& to, double coef = 1.0);
    int sum(const std::vector<std::pair<double, int>>& terms);  // terms may be empty: returns -1
    int scale(int scalar, int a, double coef = 1.0);             // coef * scalar * a
    // sum_k coef_k * a relabelled from_k -> to_k in one node (antisymmetrisers; no copies are kept)
    int permsum(int a, const std::vector<std::pair<std::string, std::string>>& perms, const std::vector<double>& coefs);
    // Linear combination accumulator: adds coef * node to the list of `key`.
    const Node& node(int id) const { return nodes_[static_cast<std::size_t>(id)]; }
    std::size_t size() const { return nodes_.size(); }
    int find_input(const std::string& name) const;

private:
    int add(Node n, const std::string& key);
    std::vector<Node> nodes_;
    std::map<std::string, int> cse_;
    std::map<std::string, int> inputs_;
};

// Evaluation and derivatives of a graph.
//
// Memory and time: forward() evaluates depth first and adds every term of a linear combination to
// it as soon as the term is known, releasing the term; values that a later reverse or tangent pass
// needs (operands of products) are kept only with `keep`. A node whose inputs have not been set
// since it was computed is not computed again (cached), and with `volatile_inputs` the nodes that do
// not depend on those inputs are kept for the next call (e.g. the source terms of the (T)
// equations while the triples are iterated). reverse() propagates depth first as well, so the
// adjoint of a sum is not copied to all its terms at once.
class Engine {
public:
    Engine(const Graph& g, const Spaces& sp);

    void set_input(int id, const BlockTensor& v);
    const BlockTensor& value(int id) const { return val_[static_cast<std::size_t>(id)]; }
    void forward(const std::vector<int>& outputs, bool keep = true, const std::vector<int>& volatile_inputs = {});
    // Tangents of the nodes needed for `outputs` for the given input tangents (others zero); the
    // values must have been computed by forward(..., keep = true).
    void tangent(const std::vector<std::pair<int, BlockTensor>>& input_tangents, const std::vector<int>& outputs);
    const BlockTensor& tangent_of(int id) const { return tan_[static_cast<std::size_t>(id)]; }
    // Adjoints of the inputs for the seeds (node, adjoint); with `with_tangents`, also the tangents of
    // the adjoints (forward over reverse): the seeds' tangents and the node tangents of the last
    // tangent() call. The seeds' ancestors must have been computed by forward(..., keep = true).
    void reverse(const std::vector<std::pair<int, BlockTensor>>& seeds, bool with_tangents = false,
                 const std::vector<std::pair<int, BlockTensor>>& seed_tangents = {});
    const BlockTensor& adjoint(int id) const { return adj_[static_cast<std::size_t>(id)]; }
    const BlockTensor& adjoint_tangent(int id) const { return adt_[static_cast<std::size_t>(id)]; }
    // Drops values, tangents and adjoints of intermediate nodes (memory).
    void release();
    const Spaces& spaces() const { return sp_; }

private:
    std::vector<int> needed(const std::vector<int>& outputs) const;
    bool up_to_date(int id) const;
    void eval(int id);
    void use(int a);
    void propagate(int id, bool wt);
    void settle(int a, bool wt);
    const Graph& g_;
    Spaces sp_;
    std::vector<BlockTensor> val_, tan_, adj_, adt_;
    std::vector<std::vector<int>> consumers_;
    std::vector<char> operand_;             // operand of an einsum or scale node
    std::vector<std::vector<int>> anc_inputs_;
    std::vector<long> stamp_in_, stamp_node_;
    std::vector<char> valid_;
    long clock_ = 0;
    // per call
    std::vector<int> remaining_;
    std::vector<char> in_set_, is_output_, protect_, processed_;
    bool keep_ = true;
};

// Antisymmetrised projection over the particle indices (first m) and the hole indices (last m)
// of a rank-2m tensor: (1/m!^2) sum over signed permutations of both groups.
BlockTensor antisymmetrize(const BlockTensor& t);

}  // namespace so
}  // namespace macqc
