#include "macqc/structure_file.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

namespace macqc {

namespace {

std::string lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::string trim(const std::string& s) {
    std::size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}

// Lines without line terminators (LF, CRLF or CR), without a UTF-8 byte-order mark.
std::vector<std::string> split_lines(const std::string& text) {
    std::string t = text;
    if (t.size() >= 3 && static_cast<unsigned char>(t[0]) == 0xEF && static_cast<unsigned char>(t[1]) == 0xBB &&
        static_cast<unsigned char>(t[2]) == 0xBF)
        t.erase(0, 3);
    std::vector<std::string> lines;
    std::string cur;
    for (std::size_t i = 0; i < t.size(); ++i) {
        const char c = t[i];
        if (c == '\r' || c == '\n') {
            lines.push_back(cur);
            cur.clear();
            if (c == '\r' && i + 1 < t.size() && t[i + 1] == '\n') ++i;
        } else {
            cur += c;
        }
    }
    if (!cur.empty()) lines.push_back(cur);
    return lines;
}

// Fields separated by blanks, tabs or commas.
std::vector<std::string> fields(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : line) {
        if (std::isspace(static_cast<unsigned char>(c)) || c == ',') {
            if (!cur.empty()) out.push_back(cur);
            cur.clear();
        } else {
            cur += c;
        }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

bool blank(const std::string& line) { return trim(line).empty(); }

// A floating-point number (also with a Fortran D exponent); false if the field is not one.
bool to_number(const std::string& f, double* x) {
    std::string s = f;
    for (char& c : s)
        if (c == 'd' || c == 'D') c = 'e';
    if (s.empty()) return false;
    char* end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    if (end == s.c_str() || *end != '\0' || !std::isfinite(v)) return false;
    *x = v;
    return true;
}

bool to_integer(const std::string& f, long* n) {
    if (f.empty()) return false;
    char* end = nullptr;
    const long v = std::strtol(f.c_str(), &end, 10);
    if (end == f.c_str() || *end != '\0') return false;
    *n = v;
    return true;
}

std::string at_line(std::size_t i) { return "line " + std::to_string(i + 1) + ": "; }

// Element of an atom label: "C", "c", "C1", "C12", "Cl", "C(Fragment=1)", "C-C_3", "6".
int element_of_label(const std::string& label, std::size_t line) {
    std::string s = label;
    const auto paren = s.find('(');
    if (paren != std::string::npos) s = s.substr(0, paren);
    const std::string ls = lower(s);
    if (ls == "bq" || (ls.size() > 3 && ls.compare(ls.size() - 3, 3, "-bq") == 0))
        throw std::runtime_error(at_line(line) + "ghost atoms (Bq) are not supported");
    const auto dash = s.find('-');
    if (dash != std::string::npos && dash > 0) s = s.substr(0, dash);
    long z = 0;
    if (to_integer(s, &z)) {
        if (z >= 1 && z <= kMaxAtomicNumber) return static_cast<int>(z);
        throw std::runtime_error(at_line(line) + "atomic number " + s + " is not supported (elements have Z = 1 to 118)");
    }
    std::string letters;
    for (char c : s) {
        if (std::isalpha(static_cast<unsigned char>(c))) letters += c;
        else break;
    }
    if (letters.empty()) throw std::runtime_error(at_line(line) + "'" + label + "' is not an element symbol");
    if (lower(letters) == "x" || lower(letters) == "xx")
        throw std::runtime_error(at_line(line) + "dummy atoms (X) are not supported; give Cartesian coordinates of real atoms only");
    // the letters are the element symbol, in any case ("Cl", "CL", "cl"; the digits of "C1" are
    // a label number)
    const int zz = letters.size() <= 2 ? atomic_number(letters) : 0;
    if (zz == 0) throw std::runtime_error(at_line(line) + "unknown element '" + letters + "' (Z = 1 to 118)");
    return zz;
}

Atom make_atom(int Z, double x, double y, double z, double to_bohr) {
    Atom a;
    a.Z = Z;
    a.symbol = element_symbol(Z);
    a.position = {x * to_bohr, y * to_bohr, z * to_bohr};
    return a;
}

StructureFile parse_xyz(const std::vector<std::string>& L) {
    StructureFile f;
    f.format = "xyz";
    std::size_t i = 0;
    while (i < L.size() && blank(L[i])) ++i;
    if (i >= L.size()) throw std::runtime_error("the file is empty");
    long n = 0;
    {
        const auto fl = fields(L[i]);
        if (fl.empty() || !to_integer(fl[0], &n))
            throw std::runtime_error(at_line(i) + "an XYZ file starts with the number of atoms");
        if (n < 1) throw std::runtime_error(at_line(i) + "the number of atoms must be at least 1");
        if (n > 10000) throw std::runtime_error(at_line(i) + "too many atoms");
    }
    ++i;
    f.title = i < L.size() ? trim(L[i]) : "";
    ++i;
    const double to_bohr = 1.0 / kBohrInAngstrom;
    for (long k = 0; k < n; ++k, ++i) {
        if (i >= L.size())
            throw std::runtime_error("the file ends after " + std::to_string(k) + " of " + std::to_string(n) + " atoms");
        const auto fl = fields(L[i]);
        if (fl.empty())
            throw std::runtime_error(at_line(i) + "blank line: expected atom " + std::to_string(k + 1) + " of " +
                                     std::to_string(n));
        if (fl.size() < 4) throw std::runtime_error(at_line(i) + "expected 'symbol x y z'");
        double x, y, z;
        if (!to_number(fl[1], &x) || !to_number(fl[2], &y) || !to_number(fl[3], &z))
            throw std::runtime_error(at_line(i) + "the coordinates x y z must be numbers (angstrom)");
        f.atoms.push_back(make_atom(element_of_label(fl[0], i), x, y, z, to_bohr));
    }
    // further frames of a trajectory, or garbage
    while (i < L.size() && blank(L[i])) ++i;
    if (i < L.size()) {
        const auto fl = fields(L[i]);
        long m = 0;
        if (!fl.empty() && to_integer(fl[0], &m) && fl.size() == 1)
            f.warnings.push_back("The file contains several structures (frames); the first one is used.");
        else
            throw std::runtime_error(at_line(i) + "more lines than the " + std::to_string(n) +
                                     " atoms given in line 1 (correct the number of atoms)");
    }
    return f;
}

StructureFile parse_gaussian(const std::vector<std::string>& L) {
    StructureFile f;
    f.format = "gaussian";
    std::size_t i = 0;
    while (i < L.size() && (blank(L[i]) || trim(L[i])[0] == '%')) ++i;
    if (i >= L.size() || trim(L[i])[0] != '#')
        throw std::runtime_error(i < L.size() ? at_line(i) + "the route section (a line starting with #) is missing"
                                              : "the route section (a line starting with #) is missing");
    std::string route;
    while (i < L.size() && !blank(L[i])) route += " " + lower(trim(L[i++]));
    // geometry from a checkpoint file or a Z-matrix with variables cannot be read here
    if (route.find("geom=check") != std::string::npos || route.find("geom=allcheck") != std::string::npos ||
        route.find("geom(check") != std::string::npos || route.find("geom(allcheck") != std::string::npos)
        throw std::runtime_error("the route reads the geometry from a checkpoint file (Geom=Check); ΨLab needs the "
                                 "coordinates in the file");
    double to_bohr = 1.0 / kBohrInAngstrom;
    if (route.find("units=au") != std::string::npos || route.find("units(au") != std::string::npos ||
        route.find("units=bohr") != std::string::npos || route.find("units(bohr") != std::string::npos)
        to_bohr = 1.0;
    while (i < L.size() && blank(L[i])) ++i;
    std::string title;
    while (i < L.size() && !blank(L[i])) title += (title.empty() ? "" : " ") + trim(L[i++]);
    f.title = title;
    while (i < L.size() && blank(L[i])) ++i;
    if (i >= L.size()) throw std::runtime_error("the line with the charge and the multiplicity is missing");
    {
        const auto fl = fields(L[i]);
        long q = 0, m = 0;
        if (fl.size() < 2 || !to_integer(fl[0], &q) || !to_integer(fl[1], &m))
            throw std::runtime_error(at_line(i) + "expected the charge and the multiplicity, e.g. '0 1'");
        if (m < 1) throw std::runtime_error(at_line(i) + "the multiplicity must be at least 1");
        f.charge = static_cast<int>(q);
        f.multiplicity = static_cast<int>(m);
        f.has_charge = true;
        if (fl.size() > 2) f.warnings.push_back("Fragment charges and multiplicities are ignored; the first pair is used.");
    }
    ++i;
    for (; i < L.size() && !blank(L[i]); ++i) {
        const auto fl = fields(L[i]);
        // label x y z, or label freeze-code x y z
        std::vector<double> num;
        bool all_numbers = true;
        for (std::size_t k = 1; k < fl.size(); ++k) {
            double v;
            if (to_number(fl[k], &v)) num.push_back(v);
            else all_numbers = false;
        }
        if (!all_numbers || (num.size() != 3 && num.size() != 4)) {
            if (fl.size() == 1 || fl.size() == 3 || fl.size() == 5 || fl.size() == 7 || !all_numbers)
                throw std::runtime_error(at_line(i) + "this looks like a Z-matrix; ΨLab reads Cartesian coordinates only "
                                         "('symbol x y z')");
            throw std::runtime_error(at_line(i) + "expected 'symbol x y z'");
        }
        const int Z = element_of_label(fl[0], i);
        if (num.size() == 4) {
            if (num[0] != 0.0 && num[0] != -1.0 && num[0] != -2.0 && num[0] != -3.0)
                throw std::runtime_error(at_line(i) + "this looks like a Z-matrix; ΨLab reads Cartesian coordinates only "
                                         "('symbol x y z' or 'symbol freeze-code x y z')");
            num.erase(num.begin());
        }
        f.atoms.push_back(make_atom(Z, num[0], num[1], num[2], to_bohr));
    }
    if (f.atoms.empty()) throw std::runtime_error("the file contains no atoms after the charge and multiplicity line");
    return f;
}

std::string format_of(const std::string& hint, const std::vector<std::string>& L) {
    const std::string h = lower(hint);
    auto ends_with = [&](const char* ext) {
        const std::string e(ext);
        return h.size() >= e.size() && h.compare(h.size() - e.size(), e.size(), e) == 0;
    };
    if (h == "xyz" || ends_with(".xyz")) return "xyz";
    if (h == "gaussian" || h == "gjf" || h == "com" || ends_with(".gjf") || ends_with(".com") || ends_with(".gau"))
        return "gaussian";
    for (const auto& l : L) {
        const std::string t = trim(l);
        if (t.empty()) continue;
        if (t[0] == '%' || t[0] == '#') return "gaussian";
        const auto fl = fields(t);
        long n = 0;
        if (!fl.empty() && to_integer(fl[0], &n)) return "xyz";
        break;
    }
    if (!h.empty())
        throw std::runtime_error("unknown structure file type '" + hint + "' (ΨLab reads .xyz and Gaussian .gjf/.com files)");
    throw std::runtime_error("not a structure file: an XYZ file starts with the number of atoms, a Gaussian input with "
                             "% or # lines");
}

}  // namespace

StructureFile parse_structure_file(const std::string& text, const std::string& format_hint) {
    const std::vector<std::string> L = split_lines(text);
    const std::string fmt = format_of(format_hint, L);
    StructureFile f = fmt == "xyz" ? parse_xyz(L) : parse_gaussian(L);
    if (f.atoms.size() > 200) throw std::runtime_error("ΨLab accepts at most 200 atoms");
    return f;
}

}  // namespace macqc
