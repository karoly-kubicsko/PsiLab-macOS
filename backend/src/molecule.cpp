#include "macqc/molecule.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

namespace macqc {

namespace {

const char* const kSymbols[] = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
    "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
    "Fl", "Mc", "Lv", "Ts", "Og"};

std::string compact_lower(const std::string& s) {
    std::string out;
    for (char c : s)
        if (!std::isspace(static_cast<unsigned char>(c)) && c != '_' && c != '-')
            out += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return out;
}

}  // namespace

std::string element_symbol(int Z) {
    if (Z < 1 || Z > kMaxAtomicNumber) throw std::invalid_argument("element Z = " + std::to_string(Z) + " not supported");
    return kSymbols[Z - 1];
}

int atomic_number(const std::string& symbol) {
    const std::string s = compact_lower(symbol);
    if (!s.empty() && std::all_of(s.begin(), s.end(), [](char c) { return std::isdigit(static_cast<unsigned char>(c)); })) {
        const long z = s.size() <= 3 ? std::stol(s) : 0;
        return z >= 1 && z <= kMaxAtomicNumber ? static_cast<int>(z) : 0;
    }
    for (int z = 1; z <= kMaxAtomicNumber; ++z)
        if (s == compact_lower(kSymbols[z - 1])) return z;
    return 0;
}

std::string canonical_molecule_name(const std::string& name) {
    const std::string s = compact_lower(name);
    if (s == "h2") return "H2";
    if (s == "lih") return "LiH";
    if (s == "h2o" || s == "water") return "H2O";
    if (s == "h2o2" || s == "hooh" || s == "hydrogenperoxide" || s == "peroxide") return "H2O2";
    if (s == "nh3" || s == "ammonia" || s == "h3n") return "NH3";
    if (s == "heh+" || s == "heh") return "HeH+";
    if (s == "custom" || s == "custommolecule" || s == "xyz" || s == "cartesian") return "custom";
    return "";
}

namespace {

std::vector<MoleculePreset> build_presets() {
    std::vector<MoleculePreset> v;
    // atoms H-Ar: ground-state terms and multiplicities
    const char* syms[18] = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"};
    const int mult[18] = {2, 1, 2, 1, 2, 3, 4, 3, 2, 1, 2, 1, 2, 3, 4, 3, 2, 1};
    const char* term[18] = {"2S", "1S", "2S", "1S", "2P", "3P", "4S", "3P", "2P", "1S",
                            "2S", "1S", "2P", "3P", "4S", "3P", "2P", "1S"};
    for (int k = 0; k < 18; ++k) {
        MoleculePreset p;
        p.key = syms[k];
        p.title = std::string(syms[k]) + " atom";
        p.multiplicity = mult[k];
        p.state = term[k];
        p.atoms = {{syms[k], {0.0, 0.0, 0.0}}};
        p.source = "ground state of the atom";
        p.atom = true;
        v.push_back(p);
    }
    auto diatomic = [&](const char* key, const char* title, const char* a, const char* b, double r, int charge, int m,
                        const char* state, const char* source) {
        MoleculePreset p;
        p.key = key;
        p.title = title;
        p.charge = charge;
        p.multiplicity = m;
        p.state = state;
        p.atoms = {{a, {0.0, 0.0, -0.5 * r}}, {b, {0.0, 0.0, 0.5 * r}}};
        p.source = source;
        v.push_back(p);
    };
    diatomic("H2+", "H$_2^+$", "H", "H", 1.9972 * kBohrInAngstrom, 1, 2, "X 2Sigma_g+",
             "r_e = 1.9972 bohr (exact Born-Oppenheimer equilibrium)");
    {
        // H3+: equilateral triangle, r_e = 1.6500 bohr (high-accuracy potential energy surfaces)
        MoleculePreset p;
        p.key = "H3+";
        p.title = "H$_3^+$";
        p.charge = 1;
        p.multiplicity = 1;
        p.state = "X 1A1'";
        const double a = 1.6500 * kBohrInAngstrom, rc = a / std::sqrt(3.0);
        for (int k = 0; k < 3; ++k) {
            const double t = (90.0 + 120.0 * k) * 3.14159265358979323846 / 180.0;
            p.atoms.push_back({"H", {rc * std::cos(t), rc * std::sin(t), 0.0}});
        }
        p.source = "D3h, r_e = 1.6500 bohr";
        v.push_back(p);
    }
    {
        // H5+: the C2v minimum, a proton shared between an H2 and an H3+ apex (MP2/aug-cc-pVTZ, ΨLab 15.0)
        MoleculePreset p;
        p.key = "H5+";
        p.title = "H$_5^+$";
        p.charge = 1;
        p.multiplicity = 1;
        p.state = "X 1A1 (C2v)";
#include "macqc/h5plus_geometry.inc"
        p.source = "C2v minimum optimised at MP2/aug-cc-pVTZ";
        v.push_back(p);
    }
    diatomic("N2", "N$_2$", "N", "N", 1.09768, 0, 1, "X 1Sigma_g+", "r_e (Huber and Herzberg)");
    diatomic("O2", "O$_2$", "O", "O", 1.20752, 0, 3, "X 3Sigma_g-", "r_e (Huber and Herzberg)");
    diatomic("F2", "F$_2$", "F", "F", 1.41193, 0, 1, "X 1Sigma_g+", "r_e (Huber and Herzberg)");
    diatomic("Cl2", "Cl$_2$", "Cl", "Cl", 1.9879, 0, 1, "X 1Sigma_g+", "r_e (Huber and Herzberg)");
    {
        MoleculePreset p;
        p.key = "NO";
        p.title = "NO";
        p.multiplicity = 2;
        p.state = "X 2Pi";
        p.atoms = {{"N", {0.0, 0.0, 0.0}}, {"O", {0.0, 0.0, 1.15077}}};
        p.source = "r_e (Huber and Herzberg)";
        v.push_back(p);
    }
    {
        // NO2: r(N-O) = 1.1934 A, angle O-N-O = 134.1 degrees (microwave r_e)
        MoleculePreset p;
        p.key = "NO2";
        p.title = "NO$_2$";
        p.multiplicity = 2;
        p.state = "X 2A1";
        const double r = 1.1934, half = 0.5 * 134.1 * 3.14159265358979323846 / 180.0;
        p.atoms = {{"N", {0.0, 0.0, 0.0}},
                   {"O", {r * std::sin(half), 0.0, -r * std::cos(half)}},
                   {"O", {-r * std::sin(half), 0.0, -r * std::cos(half)}}};
        p.source = "r_e = 1.1934 A, 134.1 degrees (microwave spectroscopy)";
        v.push_back(p);
    }
    // ΨLab 15.0: closed-shell molecules of the structure files in examples/ (coordinates in
    // angstrom, as in the files; experimental geometries)
    auto molecule = [&](const char* key, const char* title, const char* state,
                        std::vector<std::pair<std::string, std::array<double, 3>>> atoms, const char* source) {
        MoleculePreset p;
        p.key = key;
        p.title = title;
        p.multiplicity = 1;
        p.state = state;
        p.atoms = std::move(atoms);
        p.source = source;
        v.push_back(p);
    };
    molecule("CO", "CO", "X 1Sigma+", {{"C", {0.0, 0.0, 0.0}}, {"O", {0.0, 0.0, 1.1283}}},
             "r = 1.1283 A (examples/carbon-monoxide.xyz)");
    molecule("CO2", "CO$_2$", "X 1Sigma_g+",
             {{"C", {0.0, 0.0, 0.0}}, {"O", {0.0, 0.0, 1.16}}, {"O", {0.0, 0.0, -1.16}}},
             "linear, r(C-O) = 1.160 A (examples/carbon-dioxide.xyz)");
    molecule("N2O", "N$_2$O", "X 1Sigma+",
             {{"N", {0.0, 0.0, -1.1273}}, {"N", {0.0, 0.0, 0.0}}, {"O", {0.0, 0.0, 1.1851}}},
             "linear N-N-O, r(N-N) = 1.1273 A, r(N-O) = 1.1851 A (examples/dinitrogen-monoxide.xyz)");
    molecule("H2CO", "H$_2$CO", "X 1A1",
             {{"C", {0.0, 0.0, 0.0}},
              {"O", {0.0, 0.0, 1.2033}},
              {"H", {0.0, 0.934192, -0.581709}},
              {"H", {0.0, -0.934192, -0.581709}}},
             "C2v, r(C-O) = 1.2033 A, r(C-H) = 1.1005 A, H-C-H 116.2 degrees (examples/formaldehyde.xyz)");
    molecule("CH4", "CH$_4$", "X 1A1",
             {{"C", {0.0, 0.0, 0.0}},
              {"H", {0.62758, 0.62758, 0.62758}},
              {"H", {-0.62758, -0.62758, 0.62758}},
              {"H", {-0.62758, 0.62758, -0.62758}},
              {"H", {0.62758, -0.62758, -0.62758}}},
             "Td, r(C-H) = 1.0870 A (examples/methane.xyz)");
    return v;
}

}  // namespace

const std::vector<MoleculePreset>& molecule_presets() {
    static const std::vector<MoleculePreset> v = build_presets();
    return v;
}

const MoleculePreset* find_molecule_preset(const std::string& name) {
    std::string s = compact_lower(name);
    const std::string atom_suffix = "atom";
    if (s.size() > atom_suffix.size() && s.compare(s.size() - atom_suffix.size(), atom_suffix.size(), atom_suffix) == 0)
        s = s.substr(0, s.size() - atom_suffix.size());
    const std::string plus = "plus";
    if (s.size() > plus.size() && s.compare(s.size() - plus.size(), plus.size(), plus) == 0)
        s = s.substr(0, s.size() - plus.size()) + "+";
    for (const MoleculePreset& p : molecule_presets())
        if (compact_lower(p.key) == s) return &p;
    static const std::pair<const char*, const char*> names[] = {
        {"nitrogen", "N2"}, {"oxygen", "O2"}, {"fluorine", "F2"}, {"chlorine", "Cl2"},
        {"nitricoxide", "NO"}, {"nitrogenmonoxide", "NO"}, {"nitrogendioxide", "NO2"},
        {"carbonmonoxide", "CO"}, {"carbondioxide", "CO2"}, {"dinitrogenmonoxide", "N2O"}, {"nitrousoxide", "N2O"},
        {"formaldehyde", "H2CO"}, {"methanal", "H2CO"}, {"ch2o", "H2CO"}, {"methane", "CH4"}};
    for (const auto& n : names)
        if (s == n.first)
            for (const MoleculePreset& p : molecule_presets())
                if (p.key == n.second) return &p;
    return nullptr;
}

double default_bond_length(const std::string& molecule) {
    const std::string m = canonical_molecule_name(molecule);
    if (m == "H2") return 0.74144 / kBohrInAngstrom;
    if (m == "LiH") return 1.5949 / kBohrInAngstrom;
    if (m == "H2O") return 0.9572 / kBohrInAngstrom;
    if (m == "H2O2") return 0.9670 / kBohrInAngstrom;
    if (m == "NH3") return 1.012 / kBohrInAngstrom;
    if (m == "HeH+") return 1.4632;
    throw std::invalid_argument("unknown molecule '" + molecule + "'");
}

double default_oo_bond_length() { return kDefaultPeroxideOOAngstrom / kBohrInAngstrom; }

double default_angle(const std::string& molecule) {
    const std::string m = canonical_molecule_name(molecule);
    if (m == "H2O2") return kDefaultPeroxideAngle;
    if (m == "NH3") return kDefaultAmmoniaAngle;
    return kDefaultWaterAngle;
}

int Molecule::electrons() const {
    int n = -charge;
    for (const Atom& a : atoms) n += a.nuclear_charge();  // 14.0: without ECP core electrons
    return n;
}

int Molecule::core_electrons() const {
    int n = 0;
    for (const Atom& a : atoms) n += a.core_electrons();
    return n;
}

double Molecule::nuclear_repulsion() const {
    double e = 0.0;
    for (std::size_t i = 0; i < atoms.size(); ++i)
        for (std::size_t j = 0; j < i; ++j) {
            const double dx = atoms[i].position[0] - atoms[j].position[0];
            const double dy = atoms[i].position[1] - atoms[j].position[1];
            const double dz = atoms[i].position[2] - atoms[j].position[2];
            e += atoms[i].nuclear_charge() * atoms[j].nuclear_charge() / std::sqrt(dx * dx + dy * dy + dz * dz);
        }
    return e;
}

std::string Molecule::formula() const {
    // Hill order: C, then H, then the others alphabetically; without carbon all alphabetically
    std::vector<std::pair<std::string, int>> counts;
    for (const Atom& a : atoms) {
        const std::string s = a.symbol.empty() ? element_symbol(a.Z) : a.symbol;
        auto it = std::find_if(counts.begin(), counts.end(), [&](const auto& c) { return c.first == s; });
        if (it == counts.end()) counts.emplace_back(s, 1);
        else ++it->second;
    }
    const bool carbon = std::any_of(counts.begin(), counts.end(), [](const auto& c) { return c.first == "C"; });
    auto rank = [&](const std::string& s) { return carbon ? (s == "C" ? 0 : (s == "H" ? 1 : 2)) : 2; };
    std::sort(counts.begin(), counts.end(), [&](const auto& a, const auto& b) {
        if (rank(a.first) != rank(b.first)) return rank(a.first) < rank(b.first);
        return a.first < b.first;
    });
    std::string f;
    for (const auto& c : counts) f += c.first + (c.second > 1 ? std::to_string(c.second) : "");
    if (charge == 1) f += "+";
    else if (charge == -1) f += "-";
    else if (charge != 0) f += " " + std::to_string(std::abs(charge)) + (charge > 0 ? "+" : "-");
    return f;
}

std::string Molecule::display_name() const { return custom ? formula() : name; }

Molecule make_custom_molecule(const std::vector<Atom>& in, int charge, const std::string& title) {
    if (in.empty()) throw std::invalid_argument("the molecule has no atoms");
    Molecule mol;
    mol.name = "custom";
    mol.custom = true;
    mol.title = title;
    mol.charge = charge;
    for (std::size_t i = 0; i < in.size(); ++i) {
        Atom a = in[i];
        if (a.Z < 1 || a.Z > kMaxAtomicNumber)
            throw std::invalid_argument("atom " + std::to_string(i + 1) + ": unknown element (Z = 1 to 118)");
        a.symbol = element_symbol(a.Z);
        for (double x : a.position)
            if (!std::isfinite(x)) throw std::invalid_argument("atom " + std::to_string(i + 1) + ": the coordinates must be finite numbers");
        for (std::size_t j = 0; j < mol.atoms.size(); ++j) {
            double d2 = 0.0;
            for (std::size_t k = 0; k < 3; ++k) d2 += (a.position[k] - mol.atoms[j].position[k]) * (a.position[k] - mol.atoms[j].position[k]);
            if (std::sqrt(d2) < 0.3)
                throw std::invalid_argument("atoms " + std::to_string(j + 1) + " and " + std::to_string(i + 1) +
                                            " are closer than 0.3 bohr (0.16 \xC3\x85)");
        }
        mol.atoms.push_back(a);
    }
    return mol;
}

Molecule make_molecule(const std::string& name, double R, double angle_degrees, double r_oo, double dihedral) {
    const std::string m = canonical_molecule_name(name);
    if (!(std::isfinite(R) && R > 0.0))
        throw std::invalid_argument("the bond length must be a positive, finite number");
    if (R < 0.3) throw std::invalid_argument("the bond length must be at least 0.3 bohr (0.16 \xC3\x85)");
    if (R > 50.0) throw std::invalid_argument("the bond length must be at most 50 bohr (26.5 \xC3\x85)");
    Molecule mol;
    mol.name = m;
    auto atom = [](int Z, double x, double y, double z) {
        Atom a;
        a.Z = Z;
        a.symbol = element_symbol(Z);
        a.position = {x, y, z};
        return a;
    };
    if (m == "H2") {
        mol.atoms = {atom(1, 0.0, 0.0, -0.5 * R), atom(1, 0.0, 0.0, 0.5 * R)};
    } else if (m == "LiH") {
        mol.atoms = {atom(3, 0.0, 0.0, 0.0), atom(1, 0.0, 0.0, R)};
    } else if (m == "H2O") {
        if (!(std::isfinite(angle_degrees) && angle_degrees >= 30.0 && angle_degrees <= 180.0))
            throw std::invalid_argument("the H-O-H angle must be between 30 and 180 degrees");
        const double half = 0.5 * angle_degrees * 3.14159265358979323846 / 180.0;
        const double x = R * std::sin(half), z = R * std::cos(half);
        mol.atoms = {atom(8, 0.0, 0.0, 0.0), atom(1, x, 0.0, z), atom(1, -x, 0.0, z)};
    } else if (m == "H2O2") {
        if (!(r_oo > 0.0)) r_oo = default_oo_bond_length();
        if (!std::isfinite(r_oo) || r_oo < 0.3 || r_oo > 50.0)
            throw std::invalid_argument("the O-O bond length must be between 0.3 and 50 bohr");
        if (!(std::isfinite(angle_degrees) && angle_degrees >= 30.0 && angle_degrees <= 175.0))
            throw std::invalid_argument("the O-O-H angle must be between 30 and 175 degrees");
        if (!(std::isfinite(dihedral) && dihedral >= -360.0 && dihedral <= 360.0))
            throw std::invalid_argument("the H-O-O-H dihedral angle must be between -360 and 360 degrees");
        const double deg = 3.14159265358979323846 / 180.0;
        const double c = std::cos(angle_degrees * deg), s = std::sin(angle_degrees * deg);
        const double h = 0.5 * dihedral * deg;
        // H3 and H4 related by the C2 rotation about z: (x, y, z) -> (-x, -y, z); the dihedral
        // angle H3-O1-O2-H4 has the IUPAC sign (positive: clockwise seen along O1 -> O2)
        mol.atoms = {atom(8, -0.5 * r_oo, 0.0, 0.0), atom(8, 0.5 * r_oo, 0.0, 0.0),
                     atom(1, -0.5 * r_oo + R * c, R * s * std::sin(h), R * s * std::cos(h)),
                     atom(1, 0.5 * r_oo - R * c, -R * s * std::sin(h), R * s * std::cos(h))};
    } else if (m == "HeH+") {
        // 7.0: HeH+ through the general path (correlation-consistent basis sets): He at the origin,
        // H at z = R, charge +1
        mol.atoms = {atom(2, 0.0, 0.0, 0.0), atom(1, 0.0, 0.0, R)};
        mol.charge = 1;
    } else if (m == "NH3") {
        if (!(std::isfinite(angle_degrees) && angle_degrees >= 60.0 && angle_degrees <= 120.0))
            throw std::invalid_argument("the H-N-H angle must be between 60 and 120 degrees (120: planar)");
        // H_k = R (sin t cos p_k, sin t sin p_k, cos t), p_k = 0, 120, 240 degrees, with
        // cos(angle) = (3 cos^2 t - 1) / 2
        const double ca = std::cos(angle_degrees * 3.14159265358979323846 / 180.0);
        double c2 = std::max(0.0, (2.0 * ca + 1.0) / 3.0);
        if (c2 < 1.0e-14) c2 = 0.0;  // 120 degrees: exactly planar (cos 120 deg rounds to -0.5000000000000002)
        const double ct = std::sqrt(c2), st = std::sqrt(std::max(0.0, 1.0 - c2));
        mol.atoms = {atom(7, 0.0, 0.0, 0.0)};
        for (int k = 0; k < 3; ++k) {
            const double p = k * 2.0 * 3.14159265358979323846 / 3.0;
            mol.atoms.push_back(atom(1, R * st * std::cos(p), R * st * std::sin(p), R * ct));
        }
    } else {
        throw std::invalid_argument("unknown molecule '" + name + "' (HeH+, H2, LiH, H2O, H2O2 or NH3)");
    }
    return mol;
}

}  // namespace macqc
