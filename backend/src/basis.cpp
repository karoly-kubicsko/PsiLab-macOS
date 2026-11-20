#include "macqc/basis.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <set>
#include <sstream>
#include <stdexcept>

namespace macqc {

namespace {

std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

std::string upper(std::string s) {
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

std::string lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::vector<std::string> split(const std::string& s) {
    std::istringstream in(s);
    std::vector<std::string> t;
    std::string w;
    while (in >> w) t.push_back(w);
    return t;
}

[[noreturn]] void fail(int line, const std::string& msg) {
    throw std::runtime_error("line " + std::to_string(line) + ": " + msg);
}

// Fortran-style real ("1.5D+01", "0.0297400", "2.0E-3").
bool parse_real(std::string t, double& value) {
    for (char& c : t)
        if (c == 'D' || c == 'd') c = 'E';
    if (t.empty()) return false;
    char* end = nullptr;
    value = std::strtod(t.c_str(), &end);
    return end && *end == '\0' && std::isfinite(value);
}

bool parse_int(const std::string& t, int& value) {
    if (t.empty()) return false;
    char* end = nullptr;
    const long v = std::strtol(t.c_str(), &end, 10);
    if (!end || *end != '\0' || v < -1000000 || v > 1000000) return false;
    value = static_cast<int>(v);
    return true;
}

int shell_l(const std::string& type) {
    static const std::string letters = "SPDFGHIK";
    if (type.size() == 1) {
        const auto k = letters.find(type[0]);
        if (k != std::string::npos) return static_cast<int>(k);
    }
    return -1;
}

double double_factorial(int n) {  // n!! for n >= -1
    double r = 1.0;
    for (int k = n; k > 1; k -= 2) r *= k;
    return r;
}

// The basis-set library (14.0: 33 sets in five families; 15.0: 36 in six, with STO-nG), in display order. `cartesian_d`:
// the d shells are six Cartesian functions (6D), as Gaussian uses the 6-31G family by default
// and as the Basis Set Exchange marks these shells; every other shell of every set is
// spherical (5D, 7F), Gaussian's default for all other sets.
const BasisSetInfo kBases[] = {
    // 15.0: the minimal sets (STO-2G and STO-3G from the BSE, STO-1G fitted by backend/tools/sto_ng.py)
    {"STO-1G", "sto-1g.gbs", "Pople STO-nG", false},
    {"STO-2G", "sto-2g.gbs", "Pople STO-nG", false},
    {"STO-3G", "sto-3g.gbs", "Pople STO-nG", false},
    {"6-31G", "6-31g.gbs", "Pople 6-31G", true},
    {"6-31G*", "6-31gs.gbs", "Pople 6-31G", true},
    {"6-31+G", "6-31pg.gbs", "Pople 6-31G", true},
    {"6-31+G*", "6-31pgs.gbs", "Pople 6-31G", true},
    {"6-31+G**", "6-31pgss.gbs", "Pople 6-31G", true},
    {"6-31++G*", "6-31ppgs.gbs", "Pople 6-31G", true},
    {"6-31++G**", "6-31ppgss.gbs", "Pople 6-31G", true},
    {"6-311G", "6-311g.gbs", "Pople 6-311G", false},
    {"6-311G*", "6-311gs.gbs", "Pople 6-311G", false},
    {"6-311+G", "6-311pg.gbs", "Pople 6-311G", false},
    {"6-311+G*", "6-311pgs.gbs", "Pople 6-311G", false},
    {"6-311+G**", "6-311pgss.gbs", "Pople 6-311G", false},
    {"6-311++G*", "6-311ppgs.gbs", "Pople 6-311G", false},
    {"6-311++G**", "6-311ppgss.gbs", "Pople 6-311G", false},
    {"6-311++G(2d,2p)", "6-311ppg_2d_2p_.gbs", "Pople 6-311G", false},
    {"6-311++G(3df,3pd)", "6-311ppg_3df_3pd_.gbs", "Pople 6-311G", false},
    {"def2-SVP", "def2-svp.gbs", "Karlsruhe def2", false},
    {"def2-mTZVP", "def2-mtzvp.gbs", "Karlsruhe def2", false},
    {"def2-mTZVPP", "def2-mtzvpp.gbs", "Karlsruhe def2", false},
    {"def2-TZVP", "def2-tzvp.gbs", "Karlsruhe def2", false},
    {"def2-TZVPP", "def2-tzvpp.gbs", "Karlsruhe def2", false},
    {"def2-TZVPPD", "def2-tzvppd.gbs", "Karlsruhe def2", false},
    {"def2-QZVP", "def2-qzvp.gbs", "Karlsruhe def2", false},
    {"def2-QZVPP", "def2-qzvpp.gbs", "Karlsruhe def2", false},
    {"def2-QZVPPD", "def2-qzvppd.gbs", "Karlsruhe def2", false},
    {"aug-cc-pVDZ", "aug-cc-pvdz.gbs", "Dunning aug-cc-pVXZ", false},
    {"aug-cc-pVTZ", "aug-cc-pvtz.gbs", "Dunning aug-cc-pVXZ", false},
    {"aug-cc-pVQZ", "aug-cc-pvqz.gbs", "Dunning aug-cc-pVXZ", false},
    {"aug-cc-pV5Z", "aug-cc-pv5z.gbs", "Dunning aug-cc-pVXZ", false},
    {"aug-cc-pCVDZ", "aug-cc-pcvdz.gbs", "Dunning aug-cc-pCVXZ", false},
    {"aug-cc-pCVTZ", "aug-cc-pcvtz.gbs", "Dunning aug-cc-pCVXZ", false},
    {"aug-cc-pCVQZ", "aug-cc-pcvqz.gbs", "Dunning aug-cc-pCVXZ", false},
    {"aug-cc-pCV5Z", "aug-cc-pcv5z.gbs", "Dunning aug-cc-pCVXZ", false},
};

// Spelling-independent key of a basis-set name: lower case without blanks; Gaussian's
// polarisation notation (d) = *, (d,p) = ** at the end of a Pople name; "def2svp" = "def2-svp".
std::string name_key(const std::string& name) {
    std::string k;
    for (char c : lower(trim(name)))
        if (c != ' ' && c != '\t') k += c;
    auto ends_with = [&](const std::string& t) { return k.size() >= t.size() && k.compare(k.size() - t.size(), t.size(), t) == 0; };
    if (k.rfind("6-31", 0) == 0) {
        if (ends_with("g(d,p)")) k = k.substr(0, k.size() - 5) + "**";
        else if (ends_with("g(d)")) k = k.substr(0, k.size() - 3) + "*";
    }
    if (k.rfind("def2", 0) == 0 && k.size() > 4 && k[4] != '-') k.insert(4, "-");
    return k;
}

// aug-cc-pCVXZ -> aug-cc-pVXZ (the valence set used for hydrogen), else "".
std::string valence_partner(const std::string& canonical) {
    if (canonical.rfind("aug-cc-pCV", 0) == 0) return "aug-cc-pV" + canonical.substr(10);
    return "";
}

// "H-Ar", "H, Li-Ca", "H-Ca, Ga-Kr": the elements of a basis-set file as ranges of atomic numbers.
std::string element_ranges(const BasisFile& f) {
    std::string out;
    auto it = f.elements.begin();
    while (it != f.elements.end()) {
        const int first = it->first;
        int last = first;
        auto next = std::next(it);
        while (next != f.elements.end() && next->first == last + 1) {
            last = next->first;
            ++next;
        }
        out += (out.empty() ? "" : ", ") + element_symbol(first) + (last == first ? "" : "–" + element_symbol(last));
        it = next;
    }
    return out;
}

// "ECP28": the potential for reports (its number of core electrons).
std::string format_ecp(const AtomEcp& e) { return "ECP" + std::to_string(e.core_electrons); }

std::string read_text_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot read the basis-set file " + path);
    return std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
}

}  // namespace

char angular_letter(int l) {
    static const char letters[] = "spdfghik";
    return (l >= 0 && l < 8) ? letters[l] : '?';
}

BasisFile parse_gaussian94(const std::string& text) {
    BasisFile file;
    // the lines (without the byte-order mark), numbered from 1
    std::vector<std::string> lines;
    {
        std::istringstream in(text);
        std::string raw;
        while (std::getline(in, raw)) {
            if (lines.empty() && raw.size() >= 3 && raw.compare(0, 3, "\xEF\xBB\xBF") == 0) raw.erase(0, 3);
            lines.push_back(raw);
        }
    }
    std::size_t i = 0;  // next line
    int line = 0;       // number of the line last read
    bool in_header = true;
    // the next line that is not blank and not a comment, trimmed ("" at the end of the text)
    auto next_line = [&]() -> std::string {
        while (i < lines.size()) {
            const std::string t = trim(lines[i++]);
            line = static_cast<int>(i);
            if (t.empty()) continue;
            if (t[0] == '!') {
                if (in_header) file.header += t + "\n";
                continue;
            }
            in_header = false;
            return t;
        }
        line = static_cast<int>(lines.size());
        return "";
    };
    // 14.0: "<El>-ECP L n_core", then L + 1 parts (the local part U_L first, then U_l - U_L for
    // l = 0 .. L-1): a label line, the number of terms, "n exponent coefficient" per term.
    auto read_ecp = [&](const std::vector<std::string>& head, int Z) {
        int L = 0, ncore = 0;
        if (head.size() != 3 || !parse_int(head[1], L) || !parse_int(head[2], ncore))
            fail(line, "expected an ECP line such as 'I-ECP 3 28'");
        if (L < 0 || L > 6) fail(line, "ECP angular momentum " + head[1] + " is beyond the supported maximum (l = 6)");
        if (ncore < 0 || ncore >= Z || ncore % 2 != 0) fail(line, "invalid number of ECP core electrons " + head[2]);
        auto ecp = std::make_shared<AtomEcp>();
        ecp->core_electrons = ncore;
        ecp->local_l = L;
        ecp->projected.resize(static_cast<std::size_t>(L));
        for (int part = 0; part <= L; ++part) {
            const std::string label = next_line();
            if (label.empty()) fail(line, "the file ends inside an ECP");
            int nterm = 0;
            const std::string count = next_line();
            if (!parse_int(count, nterm) || nterm < 0 || nterm > 100)
                fail(line, "expected the number of ECP terms, got '" + count + "'");
            std::vector<EcpTerm>& terms = part == 0 ? ecp->local : ecp->projected[static_cast<std::size_t>(part - 1)];
            for (int k = 0; k < nterm; ++k) {
                const auto t = split(next_line());
                EcpTerm term;
                if (t.size() != 3 || !parse_int(t[0], term.power) || !parse_real(t[1], term.exponent) ||
                    !parse_real(t[2], term.coefficient))
                    fail(line, "expected an ECP term 'n exponent coefficient'");
                if (term.power != 2)
                    fail(line, "ECP terms r^(n-2) exp(-a r^2) with n = " + t[0] +
                                   " are not supported (ΨLab supports the Gaussian terms, n = 2)");
                if (!(term.exponent > 0.0)) fail(line, "ECP exponents must be positive");
                if (term.coefficient != 0.0) terms.push_back(term);
            }
        }
        ElementBasis& eb = file.elements[Z];
        eb.Z = Z;
        if (eb.ecp) fail(line, "element " + element_symbol(Z) + " has two ECPs");
        eb.ecp = ecp;
    };
    for (;;) {
        std::string s = next_line();
        if (s.empty()) break;
        if (s.rfind("****", 0) == 0) continue;
        // Element line: "<symbol> 0" (a leading '-' is allowed).
        auto t = split(s);
        int zero = -1;
        std::string sym = t[0];
        if (!sym.empty() && sym[0] == '-') sym.erase(0, 1);
        if (t.size() != 2 || !parse_int(t[1], zero) || zero != 0)
            fail(line, "expected an element line such as 'O 0', got '" + s + "'");
        const int Z = atomic_number(sym);
        if (Z == 0) fail(line, "unknown or unsupported element '" + t[0] + "'");
        const int element_line = line;
        s = next_line();
        if (s.empty()) throw std::runtime_error("basis file ends without '****' after the last element");
        t = split(s);
        const std::string first = upper(t[0]);
        if (first.size() > 4 && first.compare(first.size() - 4, 4, "-ECP") == 0) {
            read_ecp(t, Z);
            continue;
        }
        if (file.elements.count(Z) && !file.elements[Z].shells.empty())
            fail(element_line, "element " + sym + " is defined twice");
        ElementBasis* current = &file.elements[Z];
        current->Z = Z;
        // Shell lines "<type> <nprim> <scale>" with their primitives, up to "****".
        for (;;) {
            if (s.rfind("****", 0) == 0) break;
            int nprim = 0;
            double scale = 1.0;
            if (t.size() != 3 || !parse_int(t[1], nprim) || !parse_real(t[2], scale))
                fail(line, "expected a shell line such as 'S 3 1.00', got '" + s + "'");
            const std::string type = upper(t[0]);
            const bool sp = (type == "SP" || type == "L");
            const int l = sp ? 0 : shell_l(type);
            if (l < 0) fail(line, "unknown shell type '" + t[0] + "'");
            if (l > kMaxAngularMomentum - 1)
                fail(line, std::string("angular momentum ") + angular_letter(l) +
                               " is beyond the supported maximum (h, l = 5)");
            if (nprim < 1 || nprim > 100) fail(line, "invalid number of primitives " + t[1]);
            if (!(scale > 0.0)) fail(line, "invalid scale factor " + t[2]);
            ShellData a, b;
            a.l = l;
            b.l = 1;
            for (int k = 0; k < nprim; ++k) {
                if (i >= lines.size()) fail(line, "file ends inside a shell");
                const std::string raw = trim(lines[i++]);
                line = static_cast<int>(i);
                const auto p = split(raw);
                double e = 0.0, c1 = 0.0, c2 = 0.0;
                if (p.size() != (sp ? 3u : 2u) || !parse_real(p[0], e) || !parse_real(p[1], c1) ||
                    (sp && !parse_real(p[2], c2)))
                    fail(line, std::string("expected ") + (sp ? "exponent and two coefficients" : "exponent and coefficient") +
                                   ", got '" + raw + "'");
                if (!(e > 0.0)) fail(line, "exponents must be positive");
                e *= scale * scale;
                a.exponents.push_back(e);
                a.coefficients.push_back(c1);
                if (sp) {
                    b.exponents.push_back(e);
                    b.coefficients.push_back(c2);
                }
            }
            current->shells.push_back(a);
            if (sp) current->shells.push_back(b);
            s = next_line();
            if (s.empty()) throw std::runtime_error("basis file ends without '****' after the last element");
            t = split(s);
        }
    }
    if (file.elements.empty()) throw std::runtime_error("no element found in the basis-set file");
    for (const auto& [Z, eb] : file.elements)
        if (eb.shells.empty())
            throw std::runtime_error("element " + element_symbol(Z) + " has an ECP but no basis functions");
    return file;
}

const std::vector<BasisSetInfo>& basis_set_library() {
    static const std::vector<BasisSetInfo> library(std::begin(kBases), std::end(kBases));
    return library;
}

const std::vector<std::string>& supported_basis_sets() {
    static const std::vector<std::string> names = [] {
        std::vector<std::string> v;
        for (const auto& b : kBases) v.push_back(b.name);
        return v;
    }();
    return names;
}

const BasisSetInfo* basis_set_info(const std::string& name) {
    const std::string k = name_key(name);
    for (const auto& b : kBases)
        if (k == name_key(b.name)) return &b;
    return nullptr;
}

std::string canonical_basis_name(const std::string& name) {
    const BasisSetInfo* b = basis_set_info(name);
    return b ? b->name : "";
}

std::string basis_file_name(const std::string& canonical_name) {
    if (const BasisSetInfo* b = basis_set_info(canonical_name)) return b->file;
    return lower(canonical_name) + ".gbs";
}

std::string basis_function_types(const std::string& name) {
    const BasisSetInfo* b = basis_set_info(name);
    if (b && b->cartesian_d) return "Cartesian d functions (6D), spherical f functions (7F)";
    return "spherical harmonics (5D, 7F)";
}

std::string default_basis_directory() {
    if (const char* env = std::getenv("PSILAB_BASIS_DIR"))
        if (*env) return env;
#ifdef MACQC_BASIS_DIR
    return MACQC_BASIS_DIR;
#else
    return "basis";
#endif
}

std::vector<double> normalized_coefficients(int l, const std::vector<double>& exps,
                                            const std::vector<double>& coefs) {
    const std::size_t n = exps.size();
    std::vector<double> c(n);
    const double pi = 3.14159265358979323846;
    for (std::size_t i = 0; i < n; ++i)
        c[i] = coefs[i] * std::pow(2.0 * exps[i] / pi, 0.75) * std::pow(4.0 * exps[i], 0.5 * l) /
               std::sqrt(double_factorial(2 * l - 1));
    // <phi|phi> with x^l-normalised primitives: sum d_i d_j (2 sqrt(a_i a_j)/(a_i + a_j))^(l + 3/2)
    double norm = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            norm += coefs[i] * coefs[j] *
                    std::pow(2.0 * std::sqrt(exps[i] * exps[j]) / (exps[i] + exps[j]), l + 1.5);
    if (!(norm > 0.0)) throw std::runtime_error("a contracted shell has zero norm");
    const double scale = 1.0 / std::sqrt(norm);
    for (double& x : c) x *= scale;
    return c;
}

BasisSet build_basis(const Molecule& mol, const std::string& basis_name,
                     const std::map<std::string, BasisFile>& files) {
    const std::string name = canonical_basis_name(basis_name);
    if (name.empty()) throw std::runtime_error("unknown basis set '" + basis_name + "'");
    auto main_file = files.find(name);
    if (main_file == files.end()) throw std::runtime_error("basis set " + name + " was not loaded");

    BasisSet basis;
    basis.name = name;
    basis.cartesian_d = basis_set_info(name)->cartesian_d;
    basis.source = main_file->second.header;
    std::set<std::string> extra_sources;
    for (std::size_t a = 0; a < mol.atoms.size(); ++a) {
        const Atom& atom = mol.atoms[a];
        AtomBasisInfo info;
        info.symbol = atom.symbol;
        info.basis = name;
        const ElementBasis* eb = nullptr;
        auto it = main_file->second.elements.find(atom.Z);
        if (it != main_file->second.elements.end()) {
            eb = &it->second;
        } else {
            const std::string partner = valence_partner(name);
            if (atom.Z <= 2 && !partner.empty()) {
                auto pf = files.find(partner);
                if (pf == files.end()) throw std::runtime_error("basis set " + partner + " was not loaded");
                auto pit = pf->second.elements.find(atom.Z);
                if (pit == pf->second.elements.end())
                    throw std::runtime_error(partner + " has no functions for " + atom.symbol);
                eb = &pit->second;
                info.basis = partner;
                info.note = atom.Z == 1 ? name + " is not defined for H (hydrogen has no core electrons); the " + partner +
                                              " functions are used for H, as is customary."
                                        : name + " is not defined for He (the 1s shell of helium is its valence shell); "
                                                 "the " + partner + " functions are used for He, as for H.";
                extra_sources.insert(partner);
            } else {
                std::string msg = name + " is not available for " + atom.symbol +
                                  ": the Basis Set Exchange does not define it for this element.";
                const std::string family = basis_set_info(name)->family;
                const bool dunning = family.rfind("Dunning", 0) == 0;
                if (!dunning) {
                    // 14.0: the Pople and def2 sets, every element for which the BSE defines them
                    msg = name + " is not available for " + atom.symbol + ": ΨLab has it for " +
                          element_ranges(main_file->second) + " (the elements for which the Basis Set Exchange "
                          "defines it).";
                    if (family == "Karlsruhe def2" && atom.Z >= 58 && atom.Z <= 71)
                        msg += " (The diffuse def2-TZVPPD and def2-QZVPPD are not defined for the lanthanides "
                               "Ce–Lu; use def2-TZVPP or def2-QZVPP.)";
                    else if (family == "Karlsruhe def2")
                        msg += " (The def2 sets end at Rn; def2-mTZVP and def2-mTZVPP also cover Th–Lr.)";
                    else if (family == "Pople 6-31G")
                        msg += " (The diffuse 6-31+G and 6-31++G sets are defined for H–Ar, 6-31G and 6-31G* for "
                               "H–Kr; the def2 sets cover H–Rn.)";
                    else
                        msg += " (6-31G and 6-31G* cover H–Kr, the def2 sets H–Rn.)";
                } else if (name.find("pCV") != std::string::npos && atom.Z > 18)
                    msg += " (The core-valence sets are defined for Li–Ar; use the aug-cc-pVXZ set.)";
                else if (atom.Z == 19 || atom.Z == 20)
                    msg += " (The Basis Set Exchange has no aug-cc-pVXZ sets for K and Ca.)";
                else if (atom.Z >= 21 && atom.Z <= 30 && name == "aug-cc-pV5Z")
                    msg += " (aug-cc-pV5Z of Sc–Zn contains i functions, l = 6; ΨLab's integrals go up to h functions. "
                           "Use aug-cc-pVQZ or a smaller set.)";
                else if (atom.Z > 36)
                    msg += " (Elements beyond Kr need effective core potentials: ΨLab has the aug-cc-pVXZ sets for "
                           "H–Ar and Sc–Kr and, 14.0, the def2 sets with their effective core potentials up to Rn: "
                           "use def2-SVP ... def2-QZVPPD.)";
                if (name == "aug-cc-pCV5Z" && atom.Z == 3)
                    msg += " (The BSE has aug-cc-pCV5Z only for B–Ar.) For " + mol.display_name() +
                           " choose aug-cc-pCVQZ or aug-cc-pV5Z.";
                throw std::runtime_error(msg);
            }
        }
        // 14.0: the effective core potential of the atom (from the file of its functions)
        basis.ecps.push_back(eb->ecp);
        if (eb->ecp) {
            const AtomEcp& e = *eb->ecp;
            info.ecp_electrons = e.core_electrons;
            info.ecp = format_ecp(e);
            basis.ecp_electrons += e.core_electrons;
        }
        int count[kMaxAngularMomentum] = {0};
        std::set<double> prims[kMaxAngularMomentum];
        for (const ShellData& sd : eb->shells) {
            Shell sh;
            sh.l = sd.l;
            sh.atom = static_cast<int>(a);
            sh.center = atom.position;
            sh.exponents = sd.exponents;
            sh.coefficients = normalized_coefficients(sd.l, sd.exponents, sd.coefficients);
            sh.pure = !(basis.cartesian_d && sd.l == 2);  // 14.0: 6D for the 6-31G family
            sh.first_function = basis.n_functions;
            basis.n_functions += sh.size();
            basis.n_primitives += static_cast<int>(sd.exponents.size()) * sh.size();
            basis.max_l = std::max(basis.max_l, sd.l);
            info.functions += sh.size();
            ++count[sd.l];
            for (double e : sd.exponents) prims[sd.l].insert(e);
            basis.shells.push_back(std::move(sh));
        }
        std::string con = "[", pri = "(";
        for (int l = 0; l < kMaxAngularMomentum; ++l) {
            if (count[l]) con += std::to_string(count[l]) + angular_letter(l);
            if (!prims[l].empty()) pri += std::to_string(prims[l].size()) + angular_letter(l);
        }
        info.contraction = con + "]";
        info.primitives = pri + ")";
        basis.atoms.push_back(info);
    }
    for (const auto& s : extra_sources) basis.source += files.at(s).header;
    return basis;
}

// A basis set without the per-atom list (built by hand, e.g. in tests) has no ECPs.
Molecule attach_ecps(const Molecule& mol, const BasisSet& basis) {
    if (!basis.ecps.empty() && basis.ecps.size() != mol.atoms.size())
        throw std::logic_error("attach_ecps: the basis set was built for another molecule");
    Molecule m = mol;
    for (std::size_t a = 0; a < m.atoms.size(); ++a) m.atoms[a].ecp = basis.ecps.empty() ? nullptr : basis.ecps[a];
    return m;
}

void check_ecps(const Molecule& mol, const BasisSet& basis) {
    bool ok = basis.ecps.empty() || basis.ecps.size() == mol.atoms.size();
    for (std::size_t a = 0; ok && a < mol.atoms.size(); ++a)
        ok = mol.atoms[a].ecp == (basis.ecps.empty() ? nullptr : basis.ecps[a]);
    if (!ok)
        throw std::logic_error("internal error: the effective core potentials of the molecule are not those of "
                               "the basis set (attach_ecps)");
}

BasisSet build_basis(const Molecule& mol, const std::string& basis_name, const std::string& directory) {
    const std::string name = canonical_basis_name(basis_name);
    if (name.empty()) {
        std::string list;
        for (const auto& b : supported_basis_sets()) list += (list.empty() ? "" : ", ") + b;
        throw std::runtime_error("unknown basis set '" + basis_name + "' (available: " + list + ")");
    }
    std::string dir = directory.empty() ? default_basis_directory() : directory;
    if (!dir.empty() && dir.back() != '/') dir += '/';
    std::map<std::string, BasisFile> files;
    auto load = [&](const std::string& n) {
        const std::string path = dir + basis_file_name(n);
        std::string text;
        try {
            text = read_text_file(path);
        } catch (const std::exception&) {
            throw std::runtime_error("cannot read the basis-set file " + path +
                                     " (set the basis directory, e.g. PSILAB_BASIS_DIR)");
        }
        try {
            files[n] = parse_gaussian94(text);
        } catch (const std::exception& e) {
            throw std::runtime_error(path + ", " + e.what());
        }
    };
    load(name);
    const std::string partner = valence_partner(name);
    bool need_partner = false;
    for (const Atom& at : mol.atoms)
        if (at.Z <= 2 && !files[name].elements.count(at.Z)) need_partner = true;
    if (need_partner && !partner.empty()) load(partner);
    return build_basis(mol, name, files);
}

}  // namespace macqc
