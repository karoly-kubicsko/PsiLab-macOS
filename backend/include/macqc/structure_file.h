// ΨLab 9.0: Cartesian structure files of custom molecules.
//
// XYZ (.xyz): line 1 the number of atoms N, line 2 a comment (ΨLab takes it as the name of the
// molecule; it may contain LaTeX such as "CO$_2$"), then N lines "symbol x y z" in angstrom (the
// symbol may also be an atomic number; further columns are ignored). A file with several frames
// (a trajectory) gives its first frame and a warning.
//
// Gaussian input (.gjf, .com): Link 0 lines (%...), the route section (# ... up to a blank line),
// the title section (up to a blank line), the line "charge multiplicity", then one line per atom
// "label x y z" or "label freeze-code x y z" (angstrom, or bohr with Units=AU in the route) up to
// a blank line. Labels such as "C1", "C(Fragment=1)" or "C-C_3" give the element; Z-matrices,
// ghost atoms (Bq, "C-Bq") and dummy atoms (X) are not accepted.
#pragma once

#include <string>
#include <vector>

#include "macqc/molecule.h"

namespace macqc {

struct StructureFile {
    std::string format;           // "xyz" or "gaussian"
    std::string title;            // the comment line (XYZ) or the title section (Gaussian input)
    std::vector<Atom> atoms;      // Z, symbol and position in bohr
    int charge = 0;               // Gaussian input; XYZ files have no charge (0)
    int multiplicity = 1;         // Gaussian input; XYZ files: 1
    bool has_charge = false;      // the file gave the charge and the multiplicity
    std::vector<std::string> warnings;
};

// Reads the text of a structure file. `format_hint` is "xyz", "gaussian", a file name whose
// extension decides (.xyz; .gjf, .com, .gau), or "" to recognise the format from the text.
// Throws std::runtime_error with an explanation (and the line number) if the text is not a
// valid structure.
StructureFile parse_structure_file(const std::string& text, const std::string& format_hint = "");

}  // namespace macqc
