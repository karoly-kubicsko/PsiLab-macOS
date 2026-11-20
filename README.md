# ΨLab — a native macOS quantum-chemistry app

**Current version: 16.0** — see [Versions](#versions).

ΨLab computes small molecules — closed and (15.0) open shells, with RHF, UHF and ROHF references —
with Hartree–Fock, Møller–Plesset perturbation theory (MP2–MP5), configuration interaction (CISD …
CISDTQ, spin-adapted Full CI), quadratic CI (QCISD … QCISDTQ) and coupled cluster (CCSD … CCSDTQ):

- **H₂, LiH, H₂O, H₂O₂, NH₃ and HeH⁺** (and custom molecules of H–Rn and Th–Lr) in 44 Gaussian
  basis sets from basissetexchange.org: the correlation-consistent sets cc-pVXZ, aug-cc-pVXZ,
  cc-pCVXZ and aug-cc-pCVXZ (X = D, T, Q, 5; spherical Gaussians up to h functions; 16.0: cc-pVXZ
  and cc-pCVXZ), the minimal STO-1G, STO-2G and STO-3G (15.0) and (14.0) the Pople sets
  6-31G … 6-31++G** (Cartesian d functions, 6D, as in Gaussian) and 6-311G … 6-311++G(3df,3pd) and
  the Karlsruhe sets def2-SVP … def2-QZVPPD for every element the Basis Set Exchange defines them
  for, with the **effective core potentials** of the def2 sets beyond Kr, with the methods **HF, CCSD, CCSD(T), CCSDT, CCSDTQ and Full CI** (all electrons correlated):
  RHF and CCSD(T) for all of them, CCSDT and CCSDTQ where they fit into memory, Full CI where it
  is feasible (H₂, HeH⁺, LiH). Before a calculation ΨLab shows its size (basis functions,
  two-electron integrals, amplitudes, determinants), memory and estimated time; it warns about
  long calculations and refuses impossible ones with an explanation (e.g. Full CI for H₂O:
  5.6·10¹¹ determinants already in aug-cc-pVDZ).
- **HeH⁺** also in a minimal STO-nG basis (n = 1, 2, 3), the Szabo–Ostlund calculation:
  Hartree–Fock identical to the book's program, Full CI with all states and their spin.

For these molecules it **optimises geometries** (redundant internal coordinates, analytic gradients
for HF, MP2–MP5, CISD … CISDTQ, QCISD … QCISDTQ, CCSD, CCSD(T) and Full CI, and (15.0) on UHF and
ROHF references also for CCSDT and CCSDTQ), finds **transition states** (Gaussian's Opt=TS:
partitioned RFO with a computed Hessian, CalcFC or CalcAll), runs **relaxed potential-energy scans**
of a bond, an angle or a dihedral angle (Gaussian's Opt=ModRedundant), computes **harmonic
frequencies, IR intensities and thermochemistry** (analytic Hessians and dipole derivatives for the
same methods on RHF, UHF and ROHF references; numerical Hessians from energies for closed-shell
CCSDT and CCSDTQ; point groups up to D₃ₕ) and runs nine job types,
from *energy* to *TS optimisation + frequencies + single-point energy*, each with a complete
output file that is written while the job runs.

It shows energies, orbital energies, SCF, coupled-cluster and Davidson iterations, the largest
cluster amplitudes, the leading CI determinants and natural-orbital occupations, the **molecular
orbitals** (canonical RHF orbitals and Full CI natural orbitals as rotatable isosurfaces or as
contour plots), a **3D molecule viewer** in the manner of PyMOL (ball and stick, sticks, spheres
or lines; rotate and zoom with the mouse, depth cueing, labels, bond lengths, click atoms to measure
distances, angles and dihedral angles) for the molecule, every optimisation step (played step by
step), every point of a relaxed scan and every normal mode (animated in 3D), thermochemical
tables and the output file, and scans potential-energy curves with the minima marked. SwiftUI
with Liquid Glass on macOS 26; runs on macOS 14 and later, on Apple Silicon and Intel Macs. The
same calculations run in the terminal on Linux and macOS
([`psilab`](#terminal-version-linux-and-macos)).

![App icon](app/icon/AppIcon-1024.png)

## Using the app

- **Download**: the repository's **Releases** page has every version, 1.0 to the current one
  (`PsiLab-<version>.app.zip`, built by CI from `main`); the newest build of a branch or pull
  request is under *Actions* ▸ *app* ▸ the run ▸ *Artifacts* ▸ **ΨLab-app** (kept 90 days; it
  contains `ΨLab.app.zip`). Unzip it; the app is not notarized: right-click `ΨLab.app` ▸ *Open*
  once, or run `xattr -dr com.apple.quarantine` on it.
- **Sidebar**: the **molecule** (HeH⁺, H₂, LiH, H₂O, H₂O₂, NH₃; 9.0: the kept custom molecules by
  name, and *Custom molecule…*, which opens an .xyz or Gaussian input file — the molecule is kept in
  the menu; its name, charge and coordinates in Å are shown, the name and charge editable; after an
  optimisation *Save in Molecule Menu* / *Save as XYZ…* / *Save as GJF…* (13.0: Gaussian input)
  keep the optimised coordinates with their level) and, for all of them except HeH⁺
  in STO-nG:
  - the **job** (a menu, in this order): 1 *Energy*, 2 *Geometry optimisation*, 3 *Frequencies*
    (at the input geometry; a warning if it is not a stationary point), 4 *Optimisation +
    frequencies* (frequencies at the optimised geometry), 5 *Optimisation + frequencies +
    single-point energy* (the single point with its own method and basis set, by default the
    same; the final energies add the thermal corrections), and (7.0) 6 *Relaxed scan*, 7 *TS
    optimisation*, 8 *TS optimisation + frequencies* (exactly one imaginary frequency required),
    9 *TS optimisation + frequencies + single-point energy*. Jobs 4, 5, 8 and 9 stop with a message
    if the optimisation does not converge;
  - the geometry — the *starting* geometry of an optimisation: the X–H bond length in Å
    (editable); H₂O the H–O–H angle; H₂O₂ R(O–O), the O–O–H angle and the H–O–O–H
    dihedral angle (C₂; 0° cis, 180° trans); NH₃ the H–N–H angle (C₃ᵥ; 120° planar). The reset
    button restores 0.74144 Å (H₂), 1.5949 Å (LiH), 0.9572 Å and 104.52° (H₂O), 1.4556 Å, 0.9670 Å,
    102.32° and 113.70° (H₂O₂), 1.012 Å and 106.7° (NH₃), 0.7743 Å (HeH⁺, Szabo–Ostlund's 1.4632 bohr); after an optimisation
    the optimised values (or the transition state) are shown with **Use optimised geometry**, which
    copies them into the inputs;
  - job 6, the **scanned coordinate**: suggested scans of the molecule (H₂O₂ torsion 0–180°, NH₃
    umbrella to the planar molecule, H₂O bend to 180°, bonds), or the type (bond, angle, dihedral)
    and its atoms, the start (the input value or a given one, Å or degrees), the number of steps and
    the step size; the sidebar shows the coordinate's value at the input geometry, the range and
    the equivalent Gaussian input line (e.g. `D 3 1 2 4 S 18 10`);
  - jobs 7–9, the **transition-state search**: the Hessian computed at the first step (CalcFC) or at
    every step (CalcAll), and the eigenvalue test (exactly one negative eigenvalue; off = Gaussian's
    NoEigenTest);
  - the **method** — HF, CCSD, CCSD(T), CCSDT, CCSDTQ or Full CI (a menu, in this order) — and
    the **basis set** (aug-cc-pVDZ … aug-cc-pV5Z, aug-cc-pCVDZ … aug-cc-pCV5Z; for HeH⁺ also
    STO-1G/2G/3G, which return to the Szabo–Ostlund calculation);
  - jobs 5 and 9: the method and basis set of the single point; frequency jobs: temperature (K)
    and pressure (Pa; buttons for 1 atm and 1 bar); optimisations, TS searches and the points of a
    scan: normal or tight convergence and the maximum number of steps;
  - *Size of the calculation*, updated as you type: basis functions per atom, two-electron
    integrals (kept in memory or integral-direct), coupled-cluster amplitudes and memory, Full CI
    determinants and memory, and for jobs each step — how it is computed (analytic, or numerical
    with the number of energies or gradients, in orange) and its time — and the estimated time of
    the whole job.
    Notes (e.g. H takes aug-cc-pVXZ in the core-valence sets) in grey, warnings in orange, and a
    red explanation when the calculation is not possible (then **Run** is disabled); the run bar and,
    before the first result, the detail area repeat the refusal or the long-calculation warning.

  For HeH⁺ in STO-nG: R in Å, the basis set (STO-1G/2G/3G or a correlation-consistent
  set, which switches to the general calculation above), Slater exponents ζ(He) and ζ(H), the
  exact-erf option; *Reset to Szabo–Ostlund defaults* restores R = 0.7743 Å (1.4632 bohr), ζ = 2.0925 / 1.24,
  STO-3G. Below: the appearance (System / Light / Dark).
- **Run** (⌘R) calculates on all processor cores. A job estimated to take longer than 30 s asks
  for confirmation first; while it runs, the sidebar shows one progress line (for example
  "Optimisation step 4 · E = -76.0418435426 E_h · max force 1.10e-04", "Hessian: displacement 3
  of 6 (S2 +)", "Lambda iteration 7 · max change 3.1e-06", SCF, coupled-cluster and Davidson
  iterations) and the elapsed time, and **Cancel** (⌘.) stops it (the output file says so).
- **Tabs** (six, so that they fit next to the sidebar at the minimum window width of 1000 pt):
  - *Summary*: energies, the job (final energy, optimised parameters, wavenumbers, ZPE, H, G, S),
    the molecule in the 3D viewer and its coordinates, basis set per atom, convergence, timings;
  - the **3D molecule viewer** (12.0; Summary, Optimisation, Frequencies, Relaxed scan), rendered
    with SceneKit in the manner of PyMOL: drag to turn the molecule about its centre (arcball),
    pinch or ⌥-scroll to zoom (a plain scroll scrolls the page); the display menu (palette button) chooses the representation —
    *ball and stick*, *sticks* (licorice), *spheres* (van der Waals radii, Bondi/Mantina) or
    *lines*, bonds split in the middle and coloured by their atoms (CPK) —, atom labels (element,
    or element and number), bond lengths in Å, depth cueing (the far side fades into the
    background) and the background (window, black, white); these settings are kept and shared by
    all viewers. **Click atoms** to pick up to four (pink shells, yellow dashes between them): the
    viewer gives the position of one, the distance of two (Å), the angle of three and the dihedral
    angle of four (degrees); a click beside the molecule clears the selection. *Reset view* returns
    to the starting orientation (the plane of a planar molecule, the heaviest atom on top),
    *Save image* writes a PNG;
  - *Methods*, with the pages *Hartree–Fock* (energy components and virial ratio, orbital
    energies, SCF iterations; for HeH⁺ also the integrals and the S/H/F/C/P matrices), *Coupled
    cluster* (energies, (T) correction, T1 diagnostic, amplitude counts and norms, largest
    amplitudes, convergence), *Full CI* (energies, ⟨S²⟩, leading determinants, natural-orbital
    occupations, Davidson iterations; for HeH⁺ the CI vector, singlet weights and all states) and
    (8.0) *Orbitals*: the canonical RHF orbitals (all occupied ones and up to ten virtual ones) and,
    after Full CI, its natural orbitals, as an energy-level diagram (occupied levels ↑↓, degenerate
    ones side by side, HOMO and LUMO marked) and a table (ε in Eₕ and eV, occupation); the chosen
    orbital as **isosurfaces** ψ = ±c (blue positive, red negative, with the molecule; drag to
    rotate, scroll to zoom; isovalue 0.002–0.5 a.u., default 0.05) or as a **contour plot** in the
    molecular plane or a plane of two principal axes (contours at ±c·2ᵏ, the nodal line dotted,
    the nuclei projected);
  - *Optimisation*: the **molecule of each step** in the 3D viewer (12.0; the steps turned onto
    each other, so that the molecule keeps its place and the atoms glide from step to step), with
    the bond lengths in Å and the starting geometry as a translucent molecule (*Starting
    geometry*), above the step's energy, ΔE, maximum force and geometry (15.0; beside them before); the buttons and the slider go to the
    first, previous, next or last step, **Play steps** shows the steps one after the other (0.7 s
    each, the last held 1.6 s, then from the start); then the start and optimised parameters with
    their change, whether H₂O kept C₂ᵥ, the gradient method, the steps table (energy, ΔE, maximum
    and RMS force and step, converged; a click on a step number shows that step), the energy and
    maximum-force charts (logarithmic, with the threshold; the shown step is marked) and, for the
    shown step, its Cartesian geometry, internal coordinates (value, force, step) and convergence
    table; **Use optimised geometry as input**. For a transition-state search also the number of
    negative Hessian eigenvalues of every step, the eigenvalues, the followed mode (eigenvalue,
    gradient component, overlap with the previous one, largest components) and the shifts λp, λn of
    the partitioned RFO step;
  - *Frequencies*, with the pages *Frequencies* (**the normal modes animated in the 3D viewer**
    (12.0): the atoms of the selected mode vibrate about the equilibrium geometry (translucent)
    along the mode, at a speed proportional to its wavenumber — keyframe animations played by
    SceneKit, r(t) = r₀ + a sin(2πt/T) d —, with *Play*/*Pause*, the amplitude (at 1 the largest
    atomic displacement is a quarter of the shortest bond) and arrows towards the turning points;
    the molecule can be turned and zoomed while it vibrates; a chooser lists all modes (number,
    symmetry, wavenumber, IR intensity); the table of modes with
    symmetry, wavenumber (imaginary in red), reduced mass, force constant and IR intensity; the IR
    stick spectrum, the selected mode in orange; how the Hessian was computed — analytic, or
    numerical with the displaced coordinates and the number of gradients or energies —, isotopes;
    7.0: the point group and, after a TS search, whether exactly one frequency is imaginary — its
    mode is selected and animated)
    and
    *Thermochemistry* (ZPE and the thermal corrections to U, H, A, G in Eₕ, kcal/mol and kJ/mol;
    E + ZPE, U, H, A, G; E, C_V and S by contribution in kcal/mol, cal/(mol K) and J/(mol K);
    moments of inertia, rotational constants and temperatures, vibrational temperatures);
  - *Scan*, with the pages *Relaxed scan* (7.0: the energy along the scanned coordinate in kcal/mol
    above the lowest point, with the parabolic minima and maxima; the optimised molecule at every
    point in the 3D viewer (12.0), played point by point; the table of points — target and measured value, energy, steps,
    maximum force; **Start a TS search from this point** copies its geometry into the inputs and
    selects job 8, Gaussian's route from a scan maximum to Opt=TS) and *Potential curve* (single-point
    energies whatever the job; a selectable step in Å — all X–H bonds stretched together at the
    other inputs — with the RHF minimum and that of the selected correlated method; long scans ask
    for confirmation and can be stopped);
  - *Output*: the output file, streamed while the job runs (it follows the end while *Follow
    output* is on), selectable, with a search field (Return: next match; arrows: next and
    previous); for HeH⁺ the FORTRAN-style report.

  *View* switches tabs and pages: Summary ⌘1, Hartree–Fock ⌘2, Coupled Cluster ⌘3, Full CI ⌘4,
  Orbitals ⌥⌘O, Optimisation ⌘5, Frequencies ⌘6, Thermochemistry ⌘7, Relaxed Scan ⌘0, Potential
  Curve ⌘8, Output ⌘9.
- **Menus**: *File ▸ Export Results as Text…* (⇧⌘E) saves, for the general calculations, the
  output file of the job — exactly the text of the *Output* tab and of `psilab` — and for HeH⁺ the summary, the
  FORTRAN-style report and, after a scan, the potential-curve table; 13.0: every exported text ends
  with the archive of the results (see below); *File ▸ Export Results as
  JSON…* (⌥⇧⌘E) saves the results as JSON; *File ▸ Export Potential Curve as Text…* saves the scan
  table with its minima; *Edit ▸ Copy Summary as Text* (⇧⌘C) copies the short summary.
  13.0: **File ▸ Open Output File…** (⇧⌘O) opens an output file — the app's text export or a
  `psilab` .out file of version 13.0 or later — or a JSON export of any version and shows all its
  results in every tab as after the calculation: summary, Hartree–Fock, coupled cluster and Full CI
  pages, the **orbitals** (computed again from the basis set and the MO coefficients of the file:
  isosurfaces, contour plots, level diagram), optimisation steps, frequencies and normal modes,
  thermochemistry, scans and the output text. A banner names the file; the sidebar keeps its inputs
  (*Run* calculates them). *File ▸ Save Optimised Geometry as XYZ… / as Gaussian Input…* write the
  optimised coordinates as an .xyz file or a Gaussian input file (.gjf): comment lines (`!`) with
  the name, the level and the energy and what the route cannot say; the route `#p method/basis`
  for the level of the optimisation — `HF`, `CCSD(Full)` or `CCSD(T,Full)` (all electrons
  correlated, as in ΨLab; Gaussian freezes the core by default), `HF` where Gaussian has no such
  method (CCSDT, CCSDTQ, Full CI); by name the sets Gaussian has built in, with the same default
  function types as ΨLab (14.0: the 6-31G family with 6D 7F, the 6-311G family, `Def2SVP`,
  `Def2TZVP`, `Def2TZVPP`, `Def2QZVP`, `Def2QZVPP`; aug-cc-pVXZ); aug-cc-pCVXZ (and 14.0:
  def2-mTZVP, def2-mTZVPP, def2-TZVPPD, def2-QZVPPD) as `Gen 5D 7F` with the functions ΨLab used
  written out after the coordinates (for aug-cc-pCVXZ H and He with aug-cc-pVXZ), or by name per
  element when they are not known —; the name as title (ASCII, without the characters Gaussian's
  title section must not contain: @ # ! – _ \ and the hyphen); charge and multiplicity; Cartesian
  coordinates in Å. ΨLab reads both back as custom molecules (13.0: Gaussian's `!` comments are
  skipped).
- **Build from source**: see [`app/SETUP.md`](app/SETUP.md) (Xcode 26, or XcodeGen with
  `xcode/project.yml`).

## Versions

### 16.0
**The changes asked for after 15.0** (the 15.0 release does not contain them): the
*Atom/Molecule* menu and the geometry parameters of the new molecules, the Dunning sets cc-pVXZ and
cc-pCVXZ before the Karlsruhe sets, the size of impossible calculations, a wider sidebar, Full CI at
the end of the method menu, en dashes in the output files and taller 3D viewers.

**The Dunning sets cc-pVXZ and cc-pCVXZ.** The basis-set menu (and `psilab`'s `basis` key) adds
cc-pVDZ, cc-pVTZ, cc-pVQZ, cc-pV5Z and the core-valence sets cc-pCVDZ, cc-pCVTZ, cc-pCVQZ, cc-pCV5Z;
the Dunning sets stand before the Karlsruhe sets in this order: cc-pVXZ, aug-cc-pVXZ, cc-pCVXZ,
aug-cc-pCVXZ (44 basis sets in all). The files come from the Basis Set Exchange (data version 0.12,
the version of the other files; its package reproduces aug-cc-pVTZ.gbs byte for byte) with every
element up to Kr it defines: cc-pVDZ, TZ and QZ H–Ar and Ca–Kr, cc-pV5Z H–Ar, Ca and Ga–Kr (the 5Z
sets of Sc–Zn contain i functions, beyond ΨLab's integrals), cc-pCVDZ, TZ and QZ Li–Ar and Ca,
cc-pCV5Z B–Ar and Ca. As in the aug-cc-pCVXZ sets, H and He take the valence functions (cc-pVXZ) in
the core-valence sets. The Gaussian input export writes cc-pVXZ by name and cc-pCVXZ as Gen with
ΨLab's functions.
Checked: SCF energies against PySCF with the same basis files agree to 1.6·10⁻¹⁰ E_h or better
(RHF: H₂O cc-pVDZ, TZ, QZ, cc-pCVDZ, TZ; N₂ cc-pV5Z, cc-pCVQZ; CO cc-pCV5Z, 290 functions; Cl₂
cc-pVDZ; CaH₂ cc-pCVDZ; HeH⁺ cc-pCVTZ; KrH⁺ cc-pVDZ; UHF: triplet O₂ cc-pVTZ).

**Taller 3D viewers.** The height of the molecule viewers follows the height of the tab: 78 % of the
visible area for the optimisation steps, the scan points and the normal modes (at least 560 and
580 pt, before 480 and 500 pt; at most 1000 pt) and 60 % in the Summary tab (at least 420 pt, before
280 pt; at most 800 pt). With the window filling the screen of a 16-inch MacBook Pro, about 630 pt
and 490 pt (1536 × 960 pt display, Dock shown).

**Atom/Molecule menu and size of the calculation.** The first section is *Atom/Molecule*: the atoms
H–Ar, then the molecules by their number of electrons (H₂⁺, H₂, HeH⁺, H₃⁺, H₅⁺, LiH, H₂O, CH₄, NH₃,
CO, N₂, NO, O₂, H₂CO, H₂O₂, F₂, CO₂, N₂O, NO₂, Cl₂), then the kept custom molecules. The molecules
added in 15.0 have geometry parameters like H₂O (`PresetGeometry.swift`): R for the diatomics and
CO₂, R(H–H) of the equilateral H₃⁺, R(N–O) and the angle of NO₂, R(N–N) and R(N–O) of N₂O, R(C–O),
R(C–H) and the angle H–C–H of H₂CO, R(C–H) of CH₄, and four distances of the C₂ᵥ H₅⁺; unchanged, the
backend's geometry is used, changed, the app sends the Cartesian coordinates; after an optimisation
*Use optimised geometry* copies the optimised parameters. The size of the calculation names another
basis set of the menu that gives the molecule exactly the same functions (e.g. LiH: def2-QZVP and
def2-QZVPP, 65 functions each, because the Basis Set Exchange defines the two sets identically for H
and Li; among H–Ar they differ for Na and Mg only; `backend/tools/basis_identity.py` generates the
table from the basis files). Checked: the basis functions per element, their sum and the number of
unique two-electron integrals agree with an independent count from the basis files for all 36 basis
sets of 15.0 and all 38 atoms and molecules. The estimated time and memory are shown also for
impossible calculations, readably (Full CI of Cl₂ in aug-cc-pCV5Z: 1.98·10⁶⁰ determinants, about
1.7·10⁵⁰ TB and 1.0·10⁷⁴ years; beyond 10³⁰⁰ "more than").

**Sidebar width, method menu, en dashes.** The sidebar is at least 340 pt wide, so the segmented
Hessian control of optimisations (Initial, CalcFC, CalcAll, RecalcFC) fits; the method menu ends
with Full CI (after coupled cluster). The output files and messages write names and ranges with en
dashes: Hartree–Fock, Møller–Plesset, Szabo–Ostlund, the bonds and internal coordinates R(O1–H2),
D(H3–O1–O2–H4), and element ranges such as H–Ar (columns are padded by characters, not bytes).

### 15.0
**Open shells, many more methods, analytic derivatives for them, new molecules, optimiser options,
spin-adapted Full CI and hardware acceleration.**

**Molecules and spin (tasks 0 and 4).** The molecule menu (and `psilab`'s `molecule` key) adds the
**atoms H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar** in their ground-state
multiplicities (²S, ¹S, ²S, ¹S, ²P, ³P, ⁴S, ³P, ²P, ¹S, …) and the molecules **H₂⁺, H₃⁺, H₅⁺, N₂, O₂,
F₂, Cl₂, NO, NO₂, CO, CO₂, N₂O, H₂CO, CH₄** (experimental rₑ for the diatomics, NO₂ 1.1934 Å and 134.1°, H₃⁺ an equilateral
triangle of 1.65 bohr, H₅⁺ the C₂ᵥ MP2/aug-cc-pVTZ structure; CO, CO₂, N₂O, H₂CO and CH₄ with the
geometries of the structure files in `examples/`, listed by formula only). **Charge and multiplicity** can be
set for every molecule: two steppers under the molecule menu of the app (with the numbers of α and
β electrons and an immediate check that they fit), `charge` and `multiplicity` in `psilab`.
**References**: RHF, **UHF** (Pople–Nesbet, DIIS on both Fock matrices, ⟨S²⟩ and Mulliken spin
populations) and **ROHF** (Roothaan's effective Fock operator with the canonicalisation of PySCF,
a pure spin state). The default is RHF for closed-shell singlets and UHF otherwise, as in
Gaussian (ROHF for open-shell Full CI); the app's *Reference* menu, `psilab`'s `reference` key or
Gaussian's prefixes (`UMP2`, `ROCCSD(T)`, `UHF`, `ROHF`) choose another.

**Methods (task 1)**: **MP2, MP3, MP4 (SDTQ), MP5, CISD, CISD(T), CISDT, CISDTQ, QCISD,
QCISD(T), QCISDT, QCISDTQ** next to HF, CCSD … CCSDTQ and Full CI, all electrons correlated, on
RHF, UHF and ROHF references. They run in a new **spin-orbital tensor engine**
(`so_engine.cpp`, `so_methods.cpp`): block-sparse tensors over the α/β spin blocks, contractions
through the blocked GEMM, a computation graph with common-subexpression elimination, and
**automatic differentiation** (reverse mode, forward tangents, forward-over-reverse). The
equations are generated (`backend/tools/so_equations.py`, 188 connected and disconnected diagrams
up to quadruples with a non-diagonal Fock matrix) and every term's sign and weight is fitted and
checked against brute force. Definitions: MPn is Rayleigh–Schrödinger perturbation theory with the
Møller–Plesset H₀ (non-canonical orbitals: f_oo + f_vv), evaluated order by order (MP5 needs the
third-order quadruples); CI-n is the eigenvalue problem in the excitations up to n; QCI-n keeps the
linear connected terms of CC-n plus (T₂T_m)_C (and ½T₂² in the quadruples equation), with E = ⟨0|F
T₁ + V T₂|0⟩ (Pople, Head-Gordon and Raghavachari, J. Chem. Phys. 87, 5968 (1987), for QCISD and
QCISD(T)); (T) is E = (1/36) Σ t₃ X + f (1/36) Σ t₃ Y with f = 1 for CCSD(T) and f = 2 for QCISD(T)
and CISD(T) (as in Psi4's fnocc). **Verified** (`macqc_so_tests`, references from
`reference/so_brute_force.py`, which evaluates every method from its definition with sparse matrices
over all determinants, independently of the generated equations): all 16 methods for H₂O (RHF),
OH and CH₂ (UHF) in STO-3G and LiH/6-31G agree within 1e-12 Eₕ; MP2, CISD, CCSD and CCSD(T) with
RHF and UHF agree with PySCF within 1e-11 Eₕ.

**Analytic gradients and Hessians.** The engine differentiates the Lagrangian L = E + Σ λ R
(Λ equations by reverse mode), builds relaxed one- and two-particle densities, and solves the
Z-vector equations; second derivatives follow the 2n + 2 rule (amplitude responses by forward
tangents, the Lagrangian's second derivatives by forward-over-reverse). Status:

| Method | RHF reference (closed shells) | UHF reference | ROHF reference |
|---|---|---|---|
| HF (the SCF) | analytic gradient and Hessian | analytic gradient and Hessian (coupled-perturbed UHF) | analytic gradient and Hessian |
| MP2–MP5, CISD … CISDTQ, QCISD … QCISDTQ | analytic gradient and Hessian | analytic gradient and Hessian | analytic gradient and Hessian |
| CCSD, CCSD(T) | analytic (10.0) | analytic gradient and Hessian | analytic gradient and Hessian |
| CCSDT, CCSDTQ | numerical (general-order engine of 4.0) | analytic gradient and Hessian | analytic gradient and Hessian |
| Full CI | analytic (10.0) | analytic gradient and Hessian (in the α orbitals) | analytic gradient and Hessian |

**Open-shell references** (`open_shell_correlated.cpp`). The energy is the Lagrangian L = ⟨Γ(t, λ),
H(C)⟩ + z·g(C), where g are the stationarity conditions of the SCF on its non-redundant rotations —
UHF: the occupied–virtual pairs of each spin; ROHF: the closed–open, closed–virtual and open–virtual
pairs of the common orbitals (non-canonical form, ε = Σ_σ F^σ N^σ) — and z their multipliers. z·g
is linear in the Hamiltonian, so its one- and two-particle density (Coulomb and exchange parts) joins
Γ; the Z-vector equations A z = −(Y − Yᵀ) use the UHF or ROHF orbital Hessian (applied with Fock
builds) and make the generalized Fock matrix Y of Γ_tot symmetric. Gradient: dE/dx = ⟨Γ_tot, H̃ˣ⟩ −
½⟨Y, S̃ˣ⟩ + V_nnˣ. Hessian (2n + 2 rule: first-order amplitudes tˣ and orbital responses Uˣ = −½S̃ˣ +
κˣ from the coupled-perturbed UHF or ROHF equations, zeroth-order λ and z): Eˣʸ = ⟨Γ_tot, Hˣʸ⟩ +
⟨Y[H^[x]], Uʸ⟩ + ⟨Y[H̃ʸ], Uˣ⟩ − ⟨Y, UˣᵀUʸ + S̃ˣUʸ + S̃ʸUˣ + UˣUʸ⟩ − ⟨W, Sˣʸ⟩ + ⟨dΓ[tˣ], H^[y]⟩ +
⟨dΓ[tʸ], H^[x]⟩ + tˣ L_tt tʸ + V_nnˣʸ, with H^[x] = H̃ˣ + H(Uˣ); the internal directions are
computed, translations and rotations follow from invariance, and three electric-field perturbations
give the dipole derivatives (IR intensities). **Full CI** on an open-shell reference uses the
orbitals of one set (ROHF, or the α orbitals of UHF; its energy is invariant to every orbital
rotation): spin-resolved densities of the CI vector, no Z-vector and no CPHF (the orbitals follow
the geometry by symmetric orthonormalisation), CI vector responses cˣ = −(H − E)⁻¹ Q H^[x] c kept
in the space of the requested spin by Löwdin's projector, and Eˣʸ = ⟨Γ, H^[xy]⟩ + 2⟨cˣ|H^[y]|c⟩.

Checked against 4-point differences (H₂O/6-31G, distorted): gradients of MP3, MP4, CISD, CISD(T),
QCISD, QCISD(T), CCSD(T), CISDT and QCISDT within 2.5e-10 Eₕ/bohr, Hessians within 2.0e-9 Eₕ/bohr²
and dipole derivatives within 5.4e-9 e. Open shells (`macqc_so_tests`): the orbital part alone (the
SCF energy as the "method") reproduces the UHF and ROHF SCF gradients and Hessians within 4e-13 and
4e-12; UHF gradients (MP3, QCISD, CISD(T)) and ROHF gradients (MP2, CCSD(T)) of OH agree with
differences of energies within 4.1e-10 Eₕ/bohr; Hessians against differences of the analytic
gradients: OH UCCSD 2.3e-9, OH ROCISD 2.1e-9, H₂O⁺ UMP3 6.2e-10, H₂O⁺ ROQCISD(T) 5.9e-10 (STO-3G),
H₂O⁺/6-31G UMP3, UMP4, UCISD, UQCISD(T), UCCSD(T), UCISDT, ROMP3, ROMP4, ROCISD, ROQCISD(T) and
ROCCSD(T) within 1.2e-9 Eₕ/bohr², dipole derivatives within 1.6e-9; open-shell Full CI: OH gradient
4.0e-10, Hessians of OH (UHF orbitals) 2.1e-9, H₂O⁺ (ROHF) 6.5e-10 and O₂ ³Σg⁻ (UHF) 3.7e-9, dipole
derivatives within 2e-10. Examples: H₂O QCISD(T)/6-31G *opt+freq*, analytic Hessian at the first
step and for the frequencies, 1669.16, 3602.37, 3765.72 cm⁻¹, 28 s on 4 cores; NO₂ ROMP2/STO-3G
*opt+freq* (analytic ROHF-based gradients and Hessian) 681.97, 1321.81, 2230.60 cm⁻¹; NO Full
CI/STO-3G (ROHF, 3300 CSFs) *opt+freq* with CalcFC 1709.59 cm⁻¹, 18 s.

**Optimiser (task 2).** ΨLab optimises in redundant internal coordinates with rational-function
(RFO) steps, a trust radius and BFGS updates (minima), and with Baker's partitioned RFO and Bofill
updates (transition states) — the Berny-type algorithm. New: **CalcFC, CalcAll and RecalcFC=N
for minimisations** (the app's *Hessian* control in *Optimisation*; `psilab`'s `opt_hessian =
calcfc | calcall | recalcfc=N`), and RecalcFC=N also for transition-state searches: the Hessian
of the method (analytic where available) at step 1 and every N-th step, updated (BFGS or Bofill) in
between; a computed Hessian of a minimisation has its negative eigenvalues replaced by their
magnitude. **GEDIIS** (Li and Frisch, J. Chem. Theory Comput. 2, 835 (2006), Gaussian's default
for minimisations) is available as `opt_algorithm = gediis` (app: *Step: GEDIIS*): once the RMS
force is below 0.01, the energy-represented DIIS interpolation of the last five points (the
coefficients minimise E(c) = Σ cᵢEᵢ − ½ Σ cᵢcⱼ (gᵢ − gⱼ)·(qᵢ − qⱼ) on the simplex), then the RFO
step from the interpolated point; it falls back to the RFO step when the result is longer than the
trust radius or not downhill in the quadratic model, and near convergence (RMS step < 2.5e-3). RFO
stays the default: on small molecules the two need the same number of steps or GEDIIS one fewer
(HF/6-31G, Lindh or analytic initial Hessian: H₂O₂ 11/11 and 10/9, NH₃ 8/7 and 7/7, distorted H₂O
9/9 and 8/8, trans-like H₂O₂ 13/13 and 12/11, RFO/GEDIIS); GEDIIS pays off for large, flexible
molecules with many soft modes, where the BFGS Hessian is poor.

**Basis sets (task 3).** **STO-1G, STO-2G and STO-3G for every molecule**: STO-2G and STO-3G from
the Basis Set Exchange (H–Xe), STO-1G (which the BSE does not have) generated by
`backend/tools/sto_ng.py` from the one-term least-squares fits of Hehre, Stewart and Pople (1s:
α = 0.270950 ζ², as in Szabo and Ostlund) with the ζ of STO-3G (H–Ca); the script verifies its
fits against Szabo and Ostlund's exponents and the BSE's STO-2G/STO-3G. **HeH⁺ takes every basis
set** in the general calculation; in STO-nG it remains the Szabo–Ostlund calculation unless the
app's *Szabo–Ostlund calculation* switch (or `szabo_ostlund = no`) is off.

**Spin-adapted Full CI (task 6).** The Full CI keeps the determinant (string) basis for σ = Hc —
the vectorised α/β-string algorithm (Knowles–Handy, Olsen) is the fastest way to apply H — and
projects every Davidson vector onto the requested total spin with Löwdin's projector built from
Ŝ² = Ŝ₋Ŝ₊ + M(M + 1). The solution therefore lies in the space of the **configuration state
functions** of that spin (their number, from the Weyl–Paldus formula, is printed), is spin-pure
for open shells, and cannot fall onto a lower state of another spin. Checked against PySCF: C ³P
−37.7162644, C ¹D −37.6586414 (6-31G), O₂ ³Σg⁻ −147.7440416 Eₕ (STO-3G). A CSF-based σ (GUGA)
would shorten the vectors by a factor of about 2–4 for low spin but makes each σ more expensive; the
truncated CI methods (CISD … CISDTQ) work with spin-orbital amplitudes, which are spin-adapted for
RHF references (the excitation spaces are closed under Ŝ²) and, like all UHF-based methods, not for
UHF references.

**Efficiency and acceleration (task 5).** All parts run on all cores (integrals, Fock builds,
transformations, contractions, Davidson). The matrix products — the core of coupled cluster,
CI, the transformations and the derivative code — can now also run on **NVIDIA GPUs (cuBLAS)** and
**AMD GPUs (hipBLAS, ROCm)**, and on an optimised CPU BLAS (OpenBLAS, MKL, BLIS) on Linux; all three
are loaded at run time (`accelerators.cpp`), so one binary runs everywhere and falls back to
the built-in blocked GEMM. Products of at least 2·10⁹ flop go to the GPU (`PSILAB_GPU_MIN_FLOPS`);
`PSILAB_GPU=off|cuda|hip` and `PSILAB_BLAS=off|<library>` override the choice; the output file
names the path in its header (*Matrices*). The calculations stay in double precision:
data-centre GPUs (A100/H100/B200, Instinct MI200/MI300) have full-rate FP64, consumer GPUs 1/32 to
1/64 of it. Apple GPUs have no FP64; on the Mac the products run in Accelerate on the matrix units
of Apple silicon. Not yet comparable with ORCA or Gaussian: the two-electron integrals (Rys/Obara–
Saika in double precision on the CPU, no density fitting) and the dense spin-orbital engine for
the new methods (no permutational packing of T₃/T₄: CISDTQ for LiH/6-31G takes about 2 minutes,
CCSDTQ longer, where the closed-shell CCSDTQ of 4.0 is used) are the limits.

**App (tasks 7 and 8).** The tab bar (*Summary … Output*) is again in the window's toolbar as in
12.0 (13.0 had wrapped the tab view in a stack for the opened-file banner, which moved the tabs
below a title bar; the banner is now inside each tab). The 3D viewers of the optimisation steps,
scan points (480 pt high) and normal modes (500 pt high) are larger and span the full width of their
panel: the data of the step or point and the chooser of the normal modes are below the viewer
instead of beside it. The method menu has sections (Hartree–Fock, Møller–Plesset, CI, QCI, coupled
cluster); the Gaussian input export writes `UMP2(Full)`, `QCISD(T,Full)` and so on.

**Terminal (`psilab`).** Every key above; `psilab input.inp > output.out` writes the output file
to the standard output (the summary goes to the standard error), `psilab input.inp output.out`
to a file. New keys: `charge`, `multiplicity`, `reference`, `opt_hessian`, `recalc_fc`,
`opt_algorithm`, `szabo_ostlund`; `molecule` takes the atoms and the new molecules.

### 14.0
**25 more basis sets**, 33 in all, for every molecule, method and job type: the Pople sets
**6-31G, 6-31G\*, 6-31+G, 6-31+G\*, 6-31+G\*\*, 6-31++G\*, 6-31++G\*\*** and **6-311G,
6-311G\*, 6-311+G, 6-311+G\*, 6-311+G\*\*, 6-311++G\*, 6-311++G\*\*, 6-311++G(2d,2p),
6-311++G(3df,3pd)**, and the Karlsruhe sets **def2-SVP, def2-mTZVP, def2-mTZVPP, def2-TZVP,
def2-TZVPP, def2-TZVPPD, def2-QZVP, def2-QZVPP, def2-QZVPPD** (Basis Set Exchange 0.12, every
element up to Kr for which it defines them; see [`basis/README.md`](basis/README.md)). The basis
menus of the app have one section per family (6-31G, 6-311G, def2, aug-cc-pVXZ, aug-cc-pCVXZ);
`psilab` also accepts Gaussian's spellings (`6-31G(d)`, `6-31+G(d,p)`, `def2SVP`).
**Cartesian d functions**: the 6-31G family has six Cartesian d functions per d shell (6D), as
Gaussian uses it by default (its manual lists 6-31G among the few sets with Cartesian d functions,
which holds for its polarised and diffuse versions; the Basis Set Exchange marks these shells
Cartesian), with
spherical f functions (7F; f occurs in 6-31G\* for Sc–Zn). The 6-311G family and the def2 sets are
spherical (5D, 7F), as in Gaussian. The convention matters: H₂O at the experimental geometry,
RHF/6-31G\* = −76.0105299763 Eₕ with 6D (19 functions), −76.0091323803 Eₕ with 5D (18 functions).
Integrals, derivative integrals, analytic gradients and Hessians, dipole derivatives, the orbital
plots and the output file handle Cartesian shells (each component x^i y^j z^k normalised, e.g.
⟨xx|yy⟩ = 1/3). The output names the function types (e.g. "6-31G\*: Cartesian d functions (6D),
spherical f functions (7F)"), and its appendix labels Cartesian functions by component (`O1 3dxy`).
**Verified against PySCF** with the same files (`reference/pyscf_basis_sets_reference.py`;
`macqc_basis_tests`): RHF energies of H₂O in all 25 sets and of molecules with Li, C, N, Cl, K, Ca,
Zn (6D 7F) and Br within 2.3e-12 Eₕ, CCSD(T) within 7e-11 Eₕ, analytic RHF and CCSD gradients within
3.1e-9 Eₕ/bohr, the RHF/6-31G\* Hessian within 7.9e-10 Eₕ/bohr², dipole moment and orbital values
within 3e-9; the analytic Hessians of RHF, CCSD, CCSD(T) and Full CI with Cartesian d functions agree
with 4-point differences of analytic gradients within 1.1e-9 Eₕ/bohr². Gaussian input files name
the sets Gaussian has built in and write out the others (`Gen 5D 7F`).
**All elements**: every set covers every element for which the Basis Set Exchange defines it — the
def2 sets H–Rn (def2-TZVPPD and def2-QZVPPD without Ce–Lu; def2-mTZVP and def2-mTZVPP also
Th–Lr), 6-311G and 6-311G\* also I (all-electron). Beyond Kr the def2 sets come with their
**effective core potentials** (ECP28 for Rb–Xe and the 4f-in-valence lanthanides Ce–Lu, ECP46 for
Cs–La, ECP60 for Hf–Rn; Stuttgart–Cologne potentials), which ΨLab 14.0 implements: semi-local ECP
integrals (the local part as three-centre overlaps, the projected parts by the expansion in
modified spherical Bessel functions with exact angular integrals and Gauss–Legendre radial
quadrature), part of the nuclear attraction of their atom, so that the analytic gradients and
Hessians of all methods include them; the electrons see the nuclear charge Z less the core
electrons. The BSE gives def2-mTZVP and def2-mTZVPP without ECPs; ΨLab attaches the ECPs their
valence bases belong to (def2-ECP for Rb–Rn, the small-core ECP60MWB for Th–Lr, as ORCA does). An
HI, CsF, AuH, PbO, LuF₃ or ThO₂ calculation agrees with PySCF to 1.7e-11 Eₕ (RHF), 4e-11 Eₕ
(CCSD(T)), gradients 3e-10 Eₕ/bohr and Hessians 3.4e-10 Eₕ/bohr²
(`reference/pyscf_ecp_reference.py`). **SCF guess**: the SCF starts from the superposition of
atomic densities (SAD: one SCF of each element in its basis functions and ECP, open shells spread
evenly over degenerate orbitals) instead of the core Hamiltonian, which failed for ionic
heavy-element molecules (CsF, LuF₃, ThO₂). Output, app and Gaussian input (`GenECP` with the
potentials written out) report the ECPs. See
[Pople and def2 basis sets](#pople-and-def2-basis-sets-version-140).

### 13.0
**Output files that can be opened again.** Every output file — the app's *Output* tab and text
export and `psilab`'s .out file — now ends with the basis set in full (exponents and contraction
coefficients), the molecular orbitals of the Orbitals tab with their AO coefficients and, before
the termination banner, an archive section with the complete results as JSON. *File ▸ Open Output
File…* (⇧⌘O) reads such a file (or a JSON export) and shows every tab as after the calculation,
the orbitals included. Optimised geometries can be saved as Gaussian input (*Save as GJF…*, *File ▸
Save Optimised Geometry as Gaussian Input…*) as well as .xyz. See
[Job types and the output file](#job-types-and-the-output-file-version-50).

### 12.0
**3D molecule viewer** in the manner of PyMOL (SceneKit), in place of the 2D drawings: the
molecule of the Summary, every optimisation step, every point of a relaxed scan and the normal
modes. Ball and stick, sticks, van der Waals spheres or lines, bonds coloured by their two atoms;
drag to turn (arcball), pinch or ⌥-scroll to zoom; depth cueing, atom labels, bond lengths, black,
white or window background, *Reset view* and *Save image* (PNG); click atoms to measure distances,
angles and dihedral angles. Optimisation steps and scan points glide into each other while they
are played; the starting geometry is a translucent molecule. The normal modes vibrate in 3D
(keyframe animations played by SceneKit, the molecule can be turned while it vibrates) with
arrows towards the turning points and a chooser of all modes.
**Faster analytic Hessians and coupled cluster**: the CCSD residual forms a few large matrix products instead of many small ones (every CCSD and CCSD(T) calculation); the analytic Hessians contract the derivative integrals with all perturbed densities in one pass over the integrals, batch the products of the orbital response, take one difference for the whole CCSD(T) density, and their response iterations need one residual per Jacobian product instead of two. Analytic Hessians including the energy and gradient at the geometry (aug-cc-pVDZ, 4 threads, the same Linux machine with the portable GEMM, identical results), 11.0 → 12.0: H₂O CCSD 14.1 → 12.1 s, CCSD(T) 23.3 → 19.7 s; HOF (C1) CCSD 103 → 78 s, CCSD(T) 212 → 160 s. Whole `psilab` jobs: H₂O CCSD(T)/aug-cc-pVDZ optimisation + frequencies 40.3 → 33.1 s, NH₃ CCSD/aug-cc-pVDZ frequencies 45.2 → 35.1 s, H₂O CCSD(T)/aug-cc-pVTZ energy 8.9 → 8.4 s. The agreement with differences of analytic gradients is unchanged (47 checks, ≤ 1.1e-8 Eₕ/bohr²).

### 11.0
**Enthalpy**: the thermochemistry tables (app, output file) give the definition of the enthalpy,
H = U + pV, instead of "H = U + kT"; for the ideal gas pV = RT per mole (k_BT per molecule), so the
values are unchanged. The *Starting geometry* (*Geometry*) section of the sidebar gives the bond
length in Å only (no longer also in bohr), and so does the Geometry section of HeH⁺ in the
Szabo–Ostlund STO-nG basis (R = 0.7743 Å, the book's 1.4632 bohr). The README gives the measured timings of the analytic
Full CI Hessian of LiH against the numerical route of earlier versions (6.0 s vs 4.1 s).

### 10.0
**Analytic Hessians of CCSD, CCSD(T) and Full CI.** Frequency jobs and transition-state searches
(CalcFC, CalcAll) with these methods compute exact second derivatives instead of differences of
analytic gradients, with analytic dipole derivatives (IR intensities) and the relaxed dipole
moment. The sidebar of a frequency or TS job has a *Hessian* section (*Analytic* by default,
*Numerical* for the central differences of analytic gradients of earlier versions; terminal:
`hessian = analytic | numerical`); the analytic Hessian falls back to the numerical one when it
would not fit into the memory limit. Minimisations (optimisations, the first point of a relaxed
scan) start from the analytic Hessian of the method — as before for HF, now also for CCSD,
CCSD(T) and Full CI — where it fits into the memory limit, otherwise from Lindh's model Hessian
(*Analytic initial Hessian* in the Optimisation section, on by default; terminal:
`initial_hessian = auto | analytic | model`). The analytic Hessians agree with 4-point differences of
analytic gradients within 1e-10–1e-8 Eₕ/bohr² (the accuracy of the differences), including
HOF without symmetry and NH₃ with degenerate orbitals; the central differences of earlier versions
were accurate to about 3e-5. See [Analytic Hessians](#analytic-hessians-version-100).

### 9.0
**Custom molecules** from structure files: *Custom molecule…* in the molecule menu (or *File ▸ Open
Structure File…*, ⌘O) reads an XYZ file (coordinates in Å; the comment line is the name of the
molecule, written in LaTeX, e.g. `CO$_2$`, `NH$_4^+$`, `\ce{H2SO4}`, shown as CO₂, NH₄⁺, H₂SO₄)
or a Gaussian input file (.gjf, .com: title, charge and multiplicity, Cartesian coordinates). The
app keeps every molecule it has read: the molecule menu lists them by name, like H₂ or H₂O. The
sidebar shows the name, formula, number of atoms and electrons, the editable name and charge, and the
coordinates in Å. After a converged geometry optimisation or transition-state search of any
molecule ΨLab offers to save the optimised coordinates — in the molecule menu, listed with the level
of the optimisation (e.g. "CO₂ · CCSD/aug-cc-pCVTZ"), or as an .xyz file whose comment line gives
the name, the level and the energy, or (13.0) as a Gaussian input file (.gjf). All methods and job types work for custom molecules (closed-shell
singlets); structures may contain any element (Z = 1 … 118), the basis sets cover H–Ar and Sc–Kr
(not K, Ca; heavier elements would need effective core potentials). See
[Custom molecules](#custom-molecules-version-90).

### 8.0
**Molecular orbitals** in the app (*Methods ▸ Orbitals*): the canonical RHF orbitals (all occupied
ones and up to ten virtual ones) and, after Full CI, the natural orbitals, as an energy-level
diagram and a table, and the chosen orbital as rotatable isosurfaces ψ = ±c or as a contour plot.
The backend computes the orbitals on grids and their isosurfaces (marching tetrahedra, vertices on
the exact isovalue); basis functions up to h, the orbitals and the Full CI natural orbitals agree
with PySCF. See [Molecular orbitals](#molecular-orbitals-version-80).

### 7.0
- **Relaxed potential-energy scans** (job 6, Gaussian's `Opt=ModRedundant` with `S`): a bond, an
  angle or a dihedral angle is frozen at start, start + step, … and everything else is optimised;
  the summary table gives the energies in kcal/mol and kJ/mol and the minima and maxima.
- **Transition states** (jobs 7–9, Gaussian's `Opt=TS`): partitioned RFO (Baker 1986) in redundant
  internal coordinates, the Hessian computed at the first step (CalcFC) or at every step
  (CalcAll), Bofill updates, trust radius 0.30, the eigenvalue test; then frequencies with exactly
  one imaginary mode and, in job 9, a single-point energy at the transition state.
- **H₂O₂ and NH₃** (nitrogen in all basis sets; improper torsions for the umbrella; point groups
  C₁ … D₃ₕ with symmetry species and symmetry numbers); **HeH⁺** in the correlation-consistent
  basis sets with every method and job type.
- The app: the scan and TS settings, the *Scan* tab with the relaxed scan, the TS steps and the
  transition-state check, drawings of non-planar molecules. See
  [Relaxed scans, transition states, H₂O₂ and NH₃](#relaxed-scans-transition-states-h2o2-and-nh3-version-70).

### 6.0
The molecule of each geometry-optimisation step is drawn in the *Optimisation* tab (bond lengths,
the H–O–H angle, the starting geometry dashed), with buttons, a slider and **Play steps**; every
normal mode is animated in the *Frequencies* tab (speed proportional to the wavenumber, adjustable
amplitude, arrows). The calculations and their results are those of 5.0; the output file and
`psilab --version` say 6.0.

### 5.0
Geometry optimisation, harmonic frequencies, IR intensities and thermochemistry, five job types,
analytic gradients (RHF, CCSD, CCSD(T), Full CI) and the analytic RHF Hessian, and a complete
output file written while the job runs — see
[Version 5.0 (next version)](#version-50-next-version) at the end of this file.

### 4.0
- **Coupled cluster**: CCSD, CCSD(T), CCSDT and CCSDTQ for H₂, LiH and H₂O, all electrons
  correlated, on the RHF orbitals; the method menu offers HF, CCSD, CCSD(T), CCSDT, CCSDTQ and
  Full CI. CCSD and CCSD(T): closed-shell, spin-adapted equations with the (vv|vv) ladder over AO
  integrals (H₂O CCSD(T) up to aug-cc-pCV5Z); CCSDT and CCSDTQ: a general-order spin-orbital
  engine whose 183 diagrams are generated and checked by `backend/tools/cc_diagrams.py`.
- Verified: CCSD and CCSD(T) agree with PySCF within 3e-10 Eₕ (H₂, LiH and H₂O, seven cases);
  CCSD, CCSDT and CCSDTQ agree with brute-force coupled cluster in the determinant space to 1e-12;
  CCSDTQ reproduces Full CI for LiH (four electrons) and every method reproduces Full CI for H₂.
- Estimates of amplitudes, memory and time; refusals above the memory limit (60 % of the RAM in
  the app) or three days; a *CC* tab; potential curves with any method; Accelerate in the app.
- Terminal version: `method = HF | CCSD | CCSD(T) | CCSDT | CCSDTQ | FCI`. The energies of all
  methods for H₂, LiH and H₂O/aug-cc-pVDZ are listed under
  [Coupled cluster (version 4.0)](#coupled-cluster-version-40).

### 3.0
- **H₂, LiH and H₂O** in eight correlation-consistent basis sets (aug-cc-pVXZ, aug-cc-pCVXZ,
  X = D … 5), read from Gaussian94 files from basissetexchange.org (`basis/`): a general Gaussian
  integral engine up to h functions, restricted Hartree–Fock (DIIS, integrals in memory or
  integral-direct, multithreaded) and a determinant-based Full CI (H₂ in all basis sets, LiH up to
  aug-cc-pCVTZ with 11 million determinants).
- Verified against PySCF (all 23 RHF energies and every Full CI energy identical to 10 decimals),
  against dense-matrix Full CI, and against published Hartree–Fock limits and exact energies.
- Size, memory and time estimates before running; warnings for long calculations; impossible
  ones are refused with an explanation (Full CI for H₂O, LiH from aug-cc-pVQZ on, LiH in
  aug-cc-pCV5Z which is not defined for Li).
- App: molecule, method and basis-set pickers, geometry inputs (bond length; H₂O bond length and
  angle), live size estimate, progress and cancellation of running calculations, new Summary,
  Hartree–Fock and Full CI views for the new molecules, a Report tab, potential curves for all
  molecules.
- Terminal version: `molecule`, `basis`, `angle`, `method = RHF | FCI`, `threads`, `memory`,
  `fci_memory`, `basis_dir` and `psilab --estimate`.

### 2.0
- **Terminal version for Linux** (`psilab`, also runs on macOS): reads a plain-text input file
  (`basis`, `R` in bohr or Å, `zeta_He`, `zeta_H`, `exact_erf`, `scan`, `scan_step`) and writes a
  plain-text output file; installs as `psilab input.inp output.txt`.
- **Bond length in Å and bohr**: the sidebar shows R (Å) and below it R (bohr), both editable;
  the summary, the minima and the exported text give both units.
- **Export Results as Text** (⇧⌘E): summary, full Hartree–Fock and Full CI report and the
  potential-curve table; JSON export moved to ⌥⇧⌘E.
- **App icon** and the new name **ΨLab** (formerly MacQC).
- Potential curve: FCI blue circles, HF orange diamonds; the minima are listed in a readout box
  that no longer overlaps the curves.

### 1.0
- Native macOS app (SwiftUI, Liquid Glass on macOS 26, runs on macOS 14+ on Apple Silicon and Intel).
- HeH⁺ in an STO-1G/2G/3G basis: restricted Hartree–Fock (a C++ port of the Szabo & Ostlund
  program, output identical to it) and Full CI with ⟨S²⟩ spin labels and singlet weights.
- Tabs Summary, Hartree–Fock, Full CI and Potential curve (scan R = 0.8 … 5.0 bohr with minima).
- Light, Dark and System appearance; Export Results as JSON; Copy Summary as Text.

## Terminal version (Linux and macOS)

The same calculation runs without the graphical interface as the terminal program `psilab`, on
Linux and macOS (Apple Silicon and Intel). It reads a plain-text input file and writes a
plain-text output file.

Build and install once (needs a C++17 compiler and CMake; on a Mac: `xcode-select --install` for
Apple's compiler and CMake, e.g. `brew install cmake`):

```bash
cmake -S . -B build && cmake --build build
sudo cmake --install build                       # installs /usr/local/bin/psilab
```

On macOS the matrix products go through Apple's Accelerate, as in the app, which is much faster
for coupled cluster, Full CI and the integral transformations than the portable blocked GEMM
(about 5 GFlop/s); `-DMACQC_USE_ACCELERATE=OFF` keeps the portable code. The output file names the
code path in its header (*Matrices*). Both builds are tested in CI (Linux and macOS).

Without administrator rights, install into your home folder instead with
`cmake --install build --prefix ~/.local` (then `~/.local/bin` must be on your `PATH`). After that
`psilab` runs from any folder:

```bash
psilab input.inp > output.out     # 15.0: the output file on the standard output
psilab input.inp output.txt       # read input.inp, write the results to output.txt
psilab input.inp -o               # output name input.out
psilab --example > input.inp      # start from a commented example input
```

The input file has one `key = value` per line; `#` starts a comment; keys are case-insensitive:

```text
basis     = STO-3G           # STO-1G, STO-2G or STO-3G
R         = 1.4632 bohr      # bond length; unit bohr (default) or angstrom, e.g. 0.7743 angstrom
zeta_He   = 2.0925           # Slater exponent of He
zeta_H    = 1.24             # Slater exponent of H
exact_erf = no               # no: the original program's erf; yes: the exact erf
scan      = yes              # potential curve, R = 0.8 ... 5.0 bohr
scan_step = 0.1              # bohr
```

Missing keys take the values shown (the Szabo–Ostlund case, without scan unless `scan = yes`).

For H₂, LiH and H₂O (restricted Hartree–Fock, see
[General molecules](#general-molecules-h2-lih-h2o-version-30)) the input names the
molecule and a basis set; `R` and `angle` default to the experimental equilibrium geometry:

```text
molecule = H2O               # H2, LiH or H2O (default HeH+)
basis    = aug-cc-pVTZ       # aug-cc-pVXZ or aug-cc-pCVXZ, X = D, T, Q, 5; 14.0: 6-31G ... 6-311++G(3df,3pd), def2-SVP ... def2-QZVPPD
R        = 0.9572 angstrom   # O-H (H2: H-H, LiH: Li-H); default 0.74144 / 1.5949 / 0.9572 Å
angle    = 104.52            # H-O-H in degrees (H2O only)
method   = RHF               # HF (RHF), CCSD, CCSD(T), CCSDT, CCSDTQ, or FCI (H2 and LiH up to aug-cc-pCVTZ)
threads  = 0                 # 0: all processor cores
memory   = 1024 MB           # two-electron integrals are kept in memory up to this size
fci_memory = 4 GB            # memory limit of Full CI and coupled cluster
job      = energy            # energy | opt | freq | opt+freq | opt+freq+sp (version 5.0)
sp_method = CCSD(T)          # opt+freq+sp: method of the final single point (default: method)
sp_basis = aug-cc-pVTZ       # opt+freq+sp: its basis set (default: basis)
temperature = 298.15 K       # frequency jobs: thermochemistry temperature
pressure = 101325 Pa         # and pressure (Pa, kPa, MPa, bar, atm or Torr)
opt_convergence = normal     # normal | tight (optimisation criteria, see below)
max_steps = 50               # optimisation steps
```

Version 7.0 adds H₂O₂ and NH₃, HeH⁺ in the correlation-consistent basis sets, relaxed scans and
transition-state searches:

```text
molecule = H2O2              # H2O2: O1, O2, H3 (on O1), H4 (on O2); NH3: N1, H2, H3, H4
R        = 0.9670 angstrom   # O-H (H2O2) or N-H (NH3); default 0.9670 / 1.012 Å
R_OO     = 1.4556 angstrom   # H2O2: O-O
angle    = 102.32            # degrees: O-O-H (H2O2) or H-N-H (NH3, 120 = planar); default 102.32 / 106.7
dihedral = 113.70            # H2O2: H-O-O-H in degrees (IUPAC sign; 0 cis, 180 trans)
job      = scan              # also: ts, ts+freq, ts+freq+sp (Opt=TS; 7, 8, 9)
scan     = D 3 1 2 4 S 18 10 # Gaussian ModRedundant: [B|A|D] N1 N2 [N3 [N4]] S nsteps stepsize (Å or degrees)
scan_start = 0               # first value (Å, bohr or degrees); default: the value of the input geometry
ts_hessian = calcfc          # ts jobs: calcfc (Hessian at the first step) or calcall (at every step)
eigen_test = yes             # ts jobs: stop unless the Hessian has exactly one negative eigenvalue
hessian  = analytic          # 10.0, freq and ts jobs with FCI, CCSD, CCSD(T): analytic or numerical
initial_hessian = auto       # 10.0, minimisations: auto/analytic (HF, FCI, CCSD, CCSD(T) where possible), model (Lindh)
```

Version 15.0 adds open shells, the atoms H–Ar, more molecules and methods, and optimiser options:

```text
molecule = O2                # also H2+, H3+, H5+, N2, F2, Cl2, NO, NO2, CO, CO2, N2O, H2CO, CH4, atoms H ... Ar
charge   = 0                 # every molecule (default: the molecule's own)
multiplicity = 3             # 2S+1 (default: the ground state, e.g. 3 for O2, 4 for N)
reference = UHF              # RHF | UHF | ROHF (default RHF for singlets, UHF otherwise)
method   = QCISD(T)          # HF, MP2 ... MP5, CISD, CISD(T), CISDT, CISDTQ, QCISD ... QCISDTQ, CCSD ... CCSDTQ, FCI; UMP2, ROMP2, ...
basis    = 6-31G*            # 15.0 also STO-1G, STO-2G, STO-3G for every molecule
job      = opt+freq
opt_hessian = recalcfc=3     # minimisations: auto | calcfc | calcall | recalcfc=N (Gaussian Opt=RecalcFC=N)
opt_algorithm = gediis       # minimisations: rfo (default) | gediis (hybrid GEDIIS/RFO)
```

`molecule = HeH+` with a Gaussian `basis` of the library (correlation-consistent; 14.0 also Pople
and def2, except the four 6-311++G sets, which do not define He) runs the general calculation (all
methods and job types); with STO-1G, STO-2G or STO-3G it remains the Szabo–Ostlund calculation.

The job keys and the output file are described under
[Job types and the output file](#job-types-and-the-output-file-version-50). For H₂, LiH and H₂O
`psilab` writes the output file while the job runs and flushes it after every line block, so
`tail -f water.out` follows SCF, coupled-cluster, Lambda, CPHF and optimisation iterations as they
happen; the terminal shows the final numbers (optimised energy and number of steps, wavenumbers,
G, H and ZPE, the final energy) when the job ends.

`psilab --estimate input.inp` prints the size of the calculation (basis functions, two-electron
integrals, memory) and a rough time estimate without computing; a calculation estimated to take
longer than 30 s prints a warning first. Calculations that cannot be done are refused with an
explanation (for example LiH in aug-cc-pCV5Z, which is not defined for Li, or Full CI for H₂O,
which would need 5.6·10¹¹ determinants already in aug-cc-pVDZ). The basis-set files are
found in `basis_dir` (input key), `$PSILAB_BASIS_DIR`, the installed copy
(`<prefix>/share/psilab/basis`, put there by `cmake --install`) or the source tree.
The output file contains the same text as the app's *File ▸ Export Results as Text*: for HeH⁺ the
summary, the full Hartree–Fock and Full CI report and, with `scan = yes`, the potential-curve table
(R in bohr and Å) with the HF and FCI minima; for H₂, LiH and H₂O the output file of the job
(version 5.0, see below). Mistakes in the input are reported with their line number;
the exit status is 0 on success, 1 if the calculation failed (e.g. no SCF convergence, an
optimisation that did not converge or a refused job) and 2 for input or file errors.

## The C++ backend

The calculations are done by a C++17 backend (standard library plus the vendored
`nlohmann/json` single header, MIT) that the app compiles directly and calls through a small C
API. The rest of this README documents the backend; its code and API keep the original working
name `macqc`.

### What it computes

For HeH⁺ (two electrons, one s-type STO-nG contracted function per atom, n = 1, 2, 3):

1. **Restricted Hartree–Fock** – a port of the Szabo & Ostlund (*Modern Quantum Chemistry*,
   Appendix B) FORTRAN program, via `reference/hartree_fock_sto3g.py`. The routine structure is
   kept (`calculate_integrals`, `collect_matrices`, `form_g`, `diag_2x2`, `scf`; primitive
   overlap/kinetic/nuclear-attraction/two-electron integrals; Boys `F0`; `derfother`), as are the
   details that make the printed output identical to the FORTRAN program: `PI = 3.1415926535898`,
   the Abramowitz–Stegun `erf` approximation (optionally `std::erf`), the analytic 2×2
   diagonalisation with its `E2 <= E1` swap rule, and the convergence test
   `sqrt(Σ(P − Pold)² / 4) < 1e-4` (max. 25 iterations, `P = 0` start).
2. **Full CI** in the converged HF MO basis: `h` and `(pq|rs)` are transformed to the MO basis, the
   4 determinants with M_s = 0 in 4 spin orbitals are coupled with Slater–Condon rules, the
   Hamiltonian is diagonalised with a built-in Jacobi solver, ⟨S²⟩ is evaluated for every
   eigenvector (singlet/triplet labels), and the ground state's singlet configuration weights
   `(1σ)²`, `(1σ)¹(2σ)¹`, `(2σ)²` are reported. The correlation energy is
   `E_FCI − E_HF(total)`, with `E_HF` the converged SCF energy as printed.

Reference values (R = 1.4632 bohr, STO-3G, ζ = 2.0925 / 1.24): E_HF,total = −2.860661992 Eh
(6 iterations), E_FCI = −2.880711848 Eh, E_corr = −0.020049857 Eh.

### General molecules: H2, LiH, H2O (version 3.0)

ΨLab 3.0 adds H₂, LiH and H₂O (closed shell) in the correlation-consistent basis sets aug-cc-pVXZ
and aug-cc-pCVXZ (X = D, T, Q, 5), up to h functions (l = 5) and up to 341 basis functions (H₂O in
aug-cc-pCV5Z), with restricted Hartree–Fock and Full CI, in the app and in `psilab`.

- **Basis sets**: the Gaussian94 files in [`basis/`](basis/README.md), downloaded from the Basis
  Set Exchange and read at run time (`basis.cpp`, G94 parser with `D` exponents, SP shells, scale
  factors). Spherical harmonics; primitive and contraction normalisation. aug-cc-pCVXZ is not
  defined for H, so aug-cc-pVXZ is used for H (and reported); aug-cc-pCV5Z does not exist for Li,
  so LiH/aug-cc-pCV5Z is refused.
- **Integrals** (`gto_integrals.cpp`): overlap and kinetic energy by Obara–Saika recurrences;
  nuclear attraction and electron repulsion by the Obara–Saika vertical recurrence per primitive
  quartet, contraction, then the Head-Gordon–Pople horizontal recurrence and the transformation to
  real solid harmonics (Helgaker, Jørgensen, Olsen, Eqs. 6.4.47–50). Boys function from a
  tabulated 7-term Taylor expansion with downward recursion (asymptotic form above T = 45),
  relative error ≈ 1e-15. Schwarz screening; 8-fold permutational symmetry.
- **RHF** (`rhf.cpp`): canonical orthogonalisation (overlap eigenvalues below 1e-7 are removed and
  reported), core-Hamiltonian guess, DIIS (8 vectors), convergence |ΔE| < 1e-10 Eₕ and
  max |X†(FPS − SPF)X| < 1e-7. The unique two-electron integrals are kept in memory when they fit
  into 1 GB (default), otherwise every Fock matrix is built integral-direct. Integrals and Fock
  builds run on all processor cores (`std::thread`). Householder/QL eigensolver.
- **Driver** (`molecular.cpp`): input validation, basis assignment, a size and time estimate
  before the calculation (warning above 30 s), the text report (geometry, basis per atom, SCF
  iterations, energy components, virial ratio, orbital energies).

- **Full CI** (`mo_integrals.cpp`, `ci.cpp`): all electrons correlated, determinants with
  M_S = 0, lowest singlet state. The two-electron integrals are transformed to the RHF orbitals
  integral-direct (shell quartets → (μν|rs) → (pq|rs), matrix products). The CI vector is a matrix
  C(I_α, I_β) over occupation strings; σ = HC is built from the one-spin string Hamiltonian F
  (σ = CF + FC) and, for the opposite-spin part, the (N−1)-electron resolution of the identity:
  D_jl(K_α,K_β) = ⟨K_α|a_j|J_α⟩⟨K_β|a_l|J_β⟩C(J_α,J_β), E = (ij|kl)·D as one matrix product
  (cost n⁶ for two electrons per spin, instead of n⁸), then the scatter back with a⁺_i a⁺_k.
  Davidson method with diagonal preconditioner; the trial vectors are kept symmetric (C = Cᵀ), which
  removes triplet components. Output: energy, correlation energy, ⟨S²⟩ (= |S₊Ψ|²), weight of the
  RHF determinant, leading determinants, natural-orbital occupation numbers.
- **Limits**: Full CI is refused, with an explanation, above 2·10⁷ determinants (160 MB per CI
  vector) or above the memory limit (`fciMemoryMB`, default 4 GB). Determinants per case:

  | Basis set | H₂ | LiH | H₂O |
  |---|---|---|---|
  | aug-cc-pVDZ | 324 | 246,016 | 5.6e11 — refused |
  | aug-cc-pVTZ | 2,116 | 5,503,716 | 2.4e15 — refused |
  | aug-cc-pVQZ | 8,464 | 62,015,625 — refused | 1.4e18 — refused |
  | aug-cc-pV5Z | 25,600 | 454,585,041 — refused | 2.5e20 — refused |
  | aug-cc-pCVDZ | 324 | 396,900 | 1.5e12 — refused |
  | aug-cc-pCVTZ | 2,116 | 11,029,041 | 9.3e15 — refused |
  | aug-cc-pCVQZ | 8,464 | 142,444,225 — refused | 6.8e18 — refused |
  | aug-cc-pCV5Z | 25,600 | – (basis not defined for Li) | 1.4e21 — refused |

  Before a calculation ΨLab estimates its memory and time and warns above 30 s.

**Verification** (`backend/tests/test_rhf.cpp`, `test_fci.cpp`, run by CTest):

- RHF total energies agree with [PySCF](https://pyscf.org) 2.15 — an independent program reading
  the same `.gbs` files with its own parser, at identical geometries — to better than 1e-9 Eₕ in
  all 23 cases (`reference/pyscf_rhf_reference.py` → `backend/tests/data/rhf_reference.json`;
  CI runs the cases up to 130 functions, the workflow input *full* runs all). Timings with 4
  threads on a 2.1 GHz x86-64 machine: H₂O aug-cc-pVQZ 2.9 s (integrals in memory), aug-cc-pCVQZ
  29 s, aug-cc-pV5Z 2.4 min, aug-cc-pCV5Z 3.9 min (integral-direct); LiH aug-cc-pV5Z 34 s; all
  other cases below 2 s.
- Integral invariants (traces and norms of S, T, V; norm of the ERI tensor; Σ(pp|qq), Σ(pq|pq))
  agree with PySCF to 1e-10, for a test basis with every l up to 5 and for H₂O/aug-cc-pV(D,T)Z;
  individual shell quartets up to (hh|hh) agree to 1e-13.
- Published Hartree–Fock limits from numerical (basis-set-free) HF calculations: H₂ at
  R = 1.4 bohr, −1.133629571 Eₕ (Kołos & Roothaan, Rev. Mod. Phys. 32, 219 (1960) give −1.13363);
  LiH at R = 3.015 bohr, −7.987352 Eₕ (Laaksonen, Pyykkö & Sundholm, Comput. Phys. Rep. 4, 313
  (1986)). The aug-cc-pVXZ energies decrease monotonically towards the limit from above:
  H₂ 4.8e-3, 6.0e-4, 1.6e-4, 1.9e-5 Eₕ and LiH 3.2e-3, 5.6e-4, 1.4e-4, 2.0e-5 Eₕ above it for
  X = D, T, Q, 5.

- Full CI energies agree with PySCF to every printed digit (10 decimals; `reference/pyscf_fci_reference.py`
  → `backend/tests/data/fci_reference.json`): H₂ in aug-cc-pV{D,T,Q,5}Z (reference: PySCF CCSD, which
  is exact for two electrons) and LiH in aug-cc-pVDZ, aug-cc-pCVDZ, aug-cc-pVTZ (PySCF's determinant
  Full CI, `fci.direct_spin0`). LiH/aug-cc-pCVTZ (11 million determinants) runs through the same
  code; it was not recomputed with PySCF, whose solver needed 2.2 h already for aug-cc-pVTZ (ΨLab:
  57 s). The string-driven σ algorithm is also compared with
  a dense Hamiltonian matrix built from spin-orbital Slater–Condon rules for 1, 2 and 5 electrons
  per spin (H₂, LiH and H₂O in truncated orbital spaces): agreement to 1e-12 Eₕ; ⟨S²⟩ = 0 to 1e-8.
- Exact non-relativistic energies: Full CI lies above them and approaches them with the basis set,
  e.g. H₂ at R = 1.4 bohr (exact −1.1744757142 Eₕ, Pachucki, Phys. Rev. A 82, 032509 (2010)):
  9.9e-3, 1.8e-3, 6.1e-4 Eₕ above it in aug-cc-pV{D,T,Q}Z. For LiH (exact ≈ −8.07055 Eₕ at
  3.015 bohr, explicitly correlated Gaussians) the core-valence set aug-cc-pCVDZ recovers most of the
  Li 1s correlation that aug-cc-pVDZ misses (1.9e-2 instead of 4.9e-2 Eₕ above the exact energy).

Full CI energies at the default geometries (E_h):

| Basis set | H₂ | LiH |
|---|---|---|
| aug-cc-pVDZ | −1.1646244795 | −8.0212084789 |
| aug-cc-pVTZ | −1.1726340136 | −8.0386964311 |
| aug-cc-pVQZ | −1.1738672462 | refused |
| aug-cc-pV5Z | −1.1742521870 | refused |
| aug-cc-pCVDZ | (= aug-cc-pVDZ) | −8.0518951298 |
| aug-cc-pCVTZ | (= aug-cc-pVTZ) | −8.0652388024 |

Timings with 4 threads (2.1 GHz x86-64; RHF, integral transformation and Davidson): H₂/aug-cc-pV5Z
21 s (1.6 GB), LiH/aug-cc-pVTZ 57 s (1.5 GB), LiH/aug-cc-pCVTZ 2.7 min (2.8 GB); all other Full CI
cases take a few seconds. The built-in estimates are 20–40 % above these times and within 5 % of
the memory.

RHF energies at the default geometries (E_h; identical to PySCF in every printed digit, with the same
threshold 1e-7 for removing near-linear dependencies — PySCF's own default, 1e-6, would remove one
function combination in LiH/aug-cc-pCVQZ, whose smallest overlap eigenvalue is 8.5e-7):

| Basis set | H₂ (R = 0.74144 Å) | LiH (R = 1.5949 Å) | H₂O (0.9572 Å, 104.52°) |
|---|---|---|---|
| aug-cc-pVDZ | -1.1287937367 (18) | -7.9841394069 (32) | -76.0414279605 (41) |
| aug-cc-pVTZ | -1.1330212852 (46) | -7.9867905250 (69) | -76.0606132999 (92) |
| aug-cc-pVQZ | -1.1334669066 (92) | -7.9872114378 (126) | -76.0660008181 (172) |
| aug-cc-pV5Z | -1.1336044053 (160) | -7.9873310001 (207) | -76.0673205023 (287) |
| aug-cc-pCVDZ | -1.1287937367 (18) | -7.9844061601 (36) | -76.0417253915 (45) |
| aug-cc-pCVTZ | -1.1330212852 (46) | -7.9869238626 (82) | -76.0608386487 (105) |
| aug-cc-pCVQZ | -1.1334669066 (92) | -7.9872421887 (155) | -76.0660945771 (201) |
| aug-cc-pCV5Z | -1.1336044053 (160) | not defined for Li | -76.0673328556 (341) |

(number of basis functions in parentheses; for H₂ the aug-cc-pCVXZ rows equal aug-cc-pVXZ, because
H takes the aug-cc-pVXZ functions)

### Coupled cluster (version 4.0)

ΨLab 4.0 adds CCSD, CCSD(T), CCSDT and CCSDTQ for H₂, LiH and H₂O on the canonical RHF orbitals,
all electrons correlated (no frozen core); the Fock matrix is used in full.

- **CCSD and CCSD(T)** (`cc_closed_shell.cpp`): closed-shell, spin-adapted CCSD equations
  (Stanton–Gauss intermediates in the form of PySCF's `rccsd`) written as matrix products. The
  integrals with at least one occupied index come from a transformation that keeps one occupied
  index in the ket (memory o·n·N(N+1)/2). The particle–particle ladder Σ_cd (ac|bd) τ_ij^cd runs in
  the AO basis: τ_ij is back-transformed with the virtual orbitals, split into its symmetric and
  antisymmetric parts, contracted with the AO integrals like an exchange matrix (o² "densities"
  at once) and transformed back. The AO integrals are kept in memory when they fit, otherwise
  recomputed with Schwarz screening, so no (vv|vv) integrals are stored. (T):
  E = Σ_ijk Σ_abc (4W_abc + W_bca + W_cab)(V_abc − V_cba)/(3D_ijk^abc), looping over i ≥ j ≥ k
  with matrix products (cost o³v⁴).
- **CCSDT and CCSDTQ** (`cc_general.cpp`): a general-order spin-orbital engine. Amplitudes t(A;I)
  are stored once per pair of sorted sets of spin orbitals, blocked by the number of α orbitals
  and addressed with the combinatorial number system (6.8·10⁷ for CCSDTQ of H₂O/aug-cc-pVDZ
  instead of 2.7·10¹¹ in a full spin-orbital array). The residuals ⟨μ|(H e^T)_C|0⟩ are the 183
  connected diagrams up to quadruples from `backend/tools/cc_diagrams.py`, which enumerates them,
  orders the cluster vertices of each for the lowest cost and determines every sign by fitting to
  a brute-force evaluation of ⟨μ|e^−T H e^T|0⟩ in a determinant space (exact to 1e-13). Each
  diagram is evaluated from the Hamiltonian vertex outwards, one cluster vertex at a time, as one
  matrix product per spin case (gather, product, signed scatter); diagrams that begin alike share
  their intermediates. The engine also runs CCSD (used for cross-checks).
- Both engines: Jacobi updates with DIIS, converged when the largest amplitude change is below
  1e-8 and |ΔE| below 1e-10 Eₕ. The report lists the energies, the (T) correction, the T1
  diagnostic, the amplitude norms per rank and the largest amplitudes.
- **Limits** (`cc_estimate.cpp`): memory and time are estimated before running; a calculation is
  refused above `fciMemoryMB` (60 % of the RAM in the app, at most 48 GB) or above three days.
  CCSDTQ for H₂O needs about 5.5 GB in aug-cc-pVDZ and is refused in larger basis sets; CCSDT
  for H₂O runs up to aug-cc-pVTZ (about 4 GB, estimated 77 min on 4 threads). For two electrons
  every method equals Full CI; for four electrons (LiH) CCSDTQ equals Full CI. ΨLab adds a note
  in both cases.

Verification:

- CCSD and CCSD(T) agree with PySCF (`reference/pyscf_cc_reference.py`: same basis files,
  geometries and linear-dependency threshold) within 3e-10 Eₕ for H₂ (aug-cc-pVDZ, aug-cc-pVTZ),
  LiH (aug-cc-pVDZ, aug-cc-pCVDZ) and H₂O (aug-cc-pVDZ, aug-cc-pCVDZ, aug-cc-pVTZ), with the AO
  integrals of the ladder in memory and recomputed.
- CCSD, CCSDT and CCSDTQ from both engines agree within 1e-12 Eₕ with brute-force coupled cluster
  in the determinant space (`backend/tests/gen_cc_reference.py`, no coupled-cluster equations
  coded): random integrals with a non-diagonal Fock matrix (6 electrons in 7 orbitals), H₂O/STO-3G
  and LiH/6-31G.
- Exactness: CCSD, CCSDT and CCSDTQ equal Full CI for H₂, and CCSDTQ equals Full CI for
  LiH/aug-cc-pVDZ (−8.0212084789 Eₕ, 246,016 determinants).

Energies at the default geometries (Eₕ, aug-cc-pVDZ):

| Method | H₂ | LiH | H₂O |
| --- | --- | --- | --- |
| RHF | −1.1287937367 | −7.9841394069 | −76.0414279605 |
| CCSD | −1.1646244795 | −8.0211854592 | −76.2707727838 |
| CCSD(T) | −1.1646244795 | −8.0212069372 | −76.2760053269 |
| CCSDT | −1.1646244795 | −8.0212097170 | −76.2762366785 |
| CCSDTQ | −1.1646244795 | −8.0212084789 | (about 5.5 GB) |
| Full CI | −1.1646244795 | −8.0212084789 | not possible |

Timings with 4 threads (2.1 GHz x86-64, portable matrix products; the app uses Accelerate):
H₂O CCSD(T) in aug-cc-pVDZ 0.5 s, aug-cc-pVTZ 6 s, aug-cc-pVQZ 51 s; H₂O CCSDT/aug-cc-pVDZ
2.6 min (20 iterations); LiH CCSDTQ/aug-cc-pVDZ 1.6 min.

### Analytic derivatives (version 5.0)

ΨLab 5.0 differentiates the energy with respect to the nuclear coordinates. The first stage is
the integral machinery and the analytic RHF gradient and Hessian (`gto_derivatives.cpp`,
`rhf_derivatives.cpp`).

- **Derivative integrals**: the derivative of a contracted spherical Gaussian with respect to its
  centre, d/dA_k [x^c exp(−a r²)] = 2a G(c + 1_k) − c_k G(c − 1_k), is written as "component
  shells" (angular momentum l + 1 with the contraction coefficients times 2a, and l − 1), each with
  an output transform whose rows are the 3(2l+1) derivatives of the spherical functions; second
  derivatives use l + 2, l and l − 2 components with 6(2l+1) rows. (14.0: the same for Cartesian
  shells, whose output transform is the normalisation of their (l+1)(l+2)/2 components.) The unchanged Obara–Saika /
  Head-Gordon–Pople engines evaluate them (shells up to l = 7 for second derivatives of h
  functions). Translational invariance gives the derivative with respect to the last centre:
  B for overlap and kinetic energy, the nucleus for nuclear attraction, D for electron repulsion
  (second derivatives from the AA, BB, CC, AB, AC, BC blocks).
- **RHF gradient**: dE/dx = Σ P h^x + Σ Γ (μν|λσ)^x − Σ W S^x + dV_nn/dx with the
  energy-weighted density W = 2 Σ_i ε_i C_i C_iᵀ, over the unique shell quartets (8-fold symmetry),
  multithreaded. When near-linearly dependent combinations are removed (canonical
  orthogonalisation), the kept orbital space rotates with the geometry; ΨLab adds the exact term
  2 Σ_{l removed, k kept} (e_lᵀ S^x e_k)(e_lᵀ F P S e_k) / (s_k (s_k − s_l)) (e, s: eigenvectors and
  eigenvalues of S), so the gradient is the derivative of the energy that is actually computed.
  PySCF omits this term.
- **RHF Hessian**: second-derivative integrals contracted with P, Γ and W, the explicit Fock
  derivatives F^(x) = h^x + G^x(P), and the coupled-perturbed HF equations for all 3N coordinates
  at once (Krylov subspace solver, orbital-energy preconditioner, one Fock build per coordinate and
  iteration). The same response gives the **dipole derivatives** (atomic polar tensors) for IR
  intensities. Without the second-order term of the orbital-space rotation, the analytic Hessian
  is exact only when no functions are removed; with removed functions, ΨLab differentiates the
  analytic gradients numerically instead (see the frequency section).

Verification (`backend/tests/test_gradients.cpp`, run by CTest; PySCF values from
`reference/pyscf_gradient_reference.py` → `backend/tests/data/gradient_reference.json`):

- Integral derivatives against 4-point central finite differences of the integrals, for a test
  basis with s … h shells on O, H and Li (62 functions): first derivatives of S, h = T + V
  (including the operator derivative), the dipole integrals and Σ Γ (μν|λσ) to 1e-9 relative
  (matrices of dS/dx, dh/dx and G^x(D) element by element); second derivatives to 1e-9 (S),
  1e-8 (h, rounding floor of the finite differences) and 3e-9 (ERI) relative.
- RHF gradient against finite differences of the energy, at distorted geometries without
  symmetry: H₂O (aug-cc-pVDZ, aug-cc-pCVDZ), LiH (aug-cc-pCVDZ), H₂ (aug-cc-pVTZ), all within
  2e-10 Eₕ/bohr. With two removed combinations (H₂O/aug-cc-pVDZ, threshold raised to 1e-2) the
  gradient including the rotation term agrees within 1.4e-10, while the term itself is
  3.2e-4 Eₕ/bohr.
- Against PySCF 2.15 (analytic gradient, analytic CPHF Hessian, dipole moment): energies within
  1e-9 Eₕ, gradients within 1.2e-10 Eₕ/bohr (H₂O aug-cc-pVDZ … aug-cc-pV5Z, LiH aug-cc-pVDZ and
  aug-cc-pCVTZ, H₂ aug-cc-pVTZ and aug-cc-pV5Z; with removed functions, ΨLab's gradient minus the
  rotation term), Hessians within 6e-10 Eₕ/bohr² (H₂O aug-cc-pVDZ and aug-cc-pVTZ, LiH
  aug-cc-pVDZ, H₂ aug-cc-pVTZ), dipole moments within 2.5e-9 e·bohr.
- RHF Hessian against finite differences of the analytic gradient (H₂O/aug-cc-pVDZ, distorted):
  1e-10 Eₕ/bohr²; dipole derivatives against finite differences of the dipole moment: 5e-10;
  translational sum rules of the Hessian and of the dipole derivatives: 1e-11.

Timings with 4 threads (2.1 GHz x86-64): H₂O/aug-cc-pVTZ gradient 1.2 s, Hessian 7 s;
H₂O/aug-cc-pVQZ gradient 10 s, Hessian 73 s (61 s derivative integrals, 12 s CPHF).

#### Correlated gradients (stage 2)

The gradients of CCSD, CCSD(T) and Full CI are analytic; CCSDT and CCSDTQ are differentiated
numerically (`gradient.cpp`, `correlated_gradient.cpp`, `cc_gradient.cpp`, `ci.cpp`).

- **Relaxed densities** (`correlated_gradient.cpp`): a method supplies its one- and
  two-particle response densities D and G in the RHF orbitals (E = Σ D h + Σ G (pq|rs) + V_nn,
  every Fock-matrix dependence written through h and (pq|rs)). ΨLab forms the orbital Lagrangian
  X_tp = 2 Σ h_tq D_qp + 4 Σ (tq|rs) G_pqrs, solves one Z-vector (CPHF) equation for the
  occupied–virtual orbital response, back-transforms G to the AO basis (stored once per
  8-fold-symmetric index class) and contracts everything with the derivative integrals of stage 1.
  The same relaxed density gives the relaxed dipole moment. The methods are invariant to
  occupied–occupied and virtual–virtual rotations, which the program checks (X_pq − X_qp below
  1e-11 in those blocks).
- **CCSD and CCSD(T)** (`cc_gradient.cpp`): closed-shell Λ equations (Stanton–Gauss
  intermediates in PySCF's formulation, DIIS), with the (vv|vv) ladder through the same
  AO-driven operator as the amplitude equations, so no (vv|vv) integrals are stored. The
  spin-traced densities are accumulated block by block into the packed two-particle density
  (the vvvv block slice by slice). CCSD(T): the Λ equations of the CCSD(T) Lagrangian get the (T)
  source terms, and the densities the (T) contributions, with triples denominators from the full
  Fock matrix (so the Lagrangian stays rotation-invariant); these terms are evaluated per ordered
  triple of occupied orbitals with v³ blocks (no o³v³ arrays).
- **Full CI** (`ci.cpp`): the two-particle density from the CI vector through (N−1)-electron
  intermediates A_xy(K) = ⟨K| a_y,β a_x,α |Ψ⟩ (opposite spins) and (N−2)-electron ones (same spin),
  each one matrix product. Full CI is invariant to all orbital rotations, so the Z-vector
  vanishes (checked).
- **Numerical gradients**: CCSDT and CCSDTQ, and correlated methods when near-linear
  dependencies were removed or the analytic gradient would not fit into the memory limit, use
  central differences of energies (step 0.005 bohr / rad) along internal coordinates: R for H₂
  and LiH (2 energies); for H₂O the symmetric stretch and the bend when the molecule is C₂ᵥ
  (4 energies, the antisymmetric derivative vanishes by symmetry), otherwise both O–H bonds and
  the angle (6 energies). The result says that it is numerical and why. The SCF of every displaced
  geometry starts from the orbitals of the reference geometry.

Verification (`backend/tests/test_correlated_gradients.cpp`; PySCF values from
`reference/pyscf_correlated_gradient_reference.py` →
`backend/tests/data/correlated_gradient_reference.json`):

- Against 4-point finite differences of the energy: CCSD and CCSD(T) for H₂O/aug-cc-pVDZ without
  symmetry within 1.7e-10 and 8.6e-11 Eₕ/bohr; LiH/aug-cc-pVDZ CCSD, CCSD(T) and Full CI within
  6e-11; H₂/aug-cc-pVTZ Full CI within 4e-10; in the *full* CI job also H₂O/aug-cc-pVTZ CCSD(T)
  within 2.7e-10 and LiH/aug-cc-pCVDZ Full CI within 1e-10.
- Against PySCF 2.15 (`pyscf.grad.ccsd`, `pyscf.grad.ccsd_t`, CASCI with all orbitals active):
  energies within 5e-12 Eₕ, gradients within 7.5e-10 Eₕ/bohr (CCSD, CCSD(T) of H₂O and LiH; CCSD of
  H₂) and 4.9e-9 (Full CI of LiH). For H₂/aug-cc-pVTZ, PySCF's CASCI gradient differs from PySCF's
  own CCSD gradient (exact for two electrons, energies equal to 1e-12) by 2.5e-7; ΨLab's Full CI
  gradient agrees with the CCSD one within 3e-10.
- The densities reproduce the energies within 3e-12 Eₕ; the blocked (T) terms agree with the
  full-array formulation within 2e-18; numerical gradients agree with the analytic ones within
  1e-5 Eₕ/bohr; the numerical CCSDT gradient of H₂ equals the analytic Full CI one within 8e-6.

Timings with 4 threads: H₂O CCSD(T) gradient in aug-cc-pVDZ 1.5 s, aug-cc-pVTZ 24 s; LiH Full CI
gradient in aug-cc-pVTZ 89 s.

### Analytic Hessians (version 10.0)

`correlated_hessian.cpp` (general part, Full CI), `cc_hessian.cpp` (CCSD, CCSD(T)),
`gto_derivatives.cpp` (first-derivative two-electron integrals of all coordinates in packed form),
`ci.cpp` (CI vector responses), `cc_closed_shell.cpp` and `cc_gradient.cpp` (residuals of the
amplitude and Λ equations, Lagrangian densities, the bilinear form of the (T) expression).

- **Framework**: the Hessian is the derivative of the relaxed-density gradient,
  E^xy = Σ D h^(xy) + Σ G g^(xy) + Σ D^(y) h^(x) + Σ G^(y) g^(x) + V_nn^xy, with D, G the
  Lagrangian densities of the method in the RHF orbitals and h^(x), g^(x) the MO integrals
  differentiated together with the orbitals (C(x) = C U(x)). U^x solves the coupled-perturbed HF
  equations of all 3N coordinates (U_oo = U_vv = −S/2: the methods are invariant to those
  rotations). The second-order response U^xy enters only through the orbital Lagrangian of the
  gradient: its symmetric part follows from orthonormality, its occupied–virtual part from the
  Z-vector of the gradient (Handy–Schaefer interchange, no second-order CPHF). The explicit
  second-derivative integrals are contracted with the relaxed densities (one-particle, two-particle
  with its separable Z-vector term, energy-weighted). With the RHF densities this reproduces the
  analytic RHF Hessian to 6e-12 Eₕ/bohr² — a test of the orbital part alone.
- **Internal directions**: the responses of the method are solved only along the 3N − 6 (3N − 5)
  Cartesian directions orthogonal to the translations and rotations; the other columns of the
  Hessian follow exactly from the invariance (H u = 0 for a translation, the rotated gradient for a
  rotation), with the gradient ⟨D, H^(x)⟩ + V_nn^x from the same quantities.
- **Full CI**: the CI vector responses (H − E) c^q = −(H^(q) − E^(q)) c by projected
  preconditioned conjugate gradients; the density responses are the bilinear densities
  [D(c + c^q) − D(c − c^q)] / 2 (exact, the densities are quadratic).
- **CCSD**: only the amplitude responses J t^q = −Ω(t; H^(q)) are solved: with the derivative of
  the Λ equations the multiplier responses drop out (2n + 2 rule),
  ⟨D^(q), H^(p)⟩ = L_t(H^(p))·t^q + L_t(H^(q))·t^p + t^p L_tt t^q. The Jacobian products are exact:
  the residual is a polynomial of degree 4 in the amplitudes, so the symmetric five-point formula
  has no truncation error; the iterations use a forward difference against the residual at t
  (12.0; one residual per product — 11.0: a central difference, two), corrected by exact
  residuals (defect correction), and the (vv|vv) ladder, linear in τ, is applied once per product. The pair
  terms t^p L_tt t^q are second derivatives of the scalar Lagrangian
  E_corr + 2 l₁·Ω₁ + (2 l₂ − l₂ᵀ)·Ω₂ (with the ladder part Z(λ)·τ).
- **CCSD(T)**: the (T) functional is homogeneous of degree one in the Hamiltonian but rational
  through the orbital-energy denominators, so the response term gains its curvature
  L_HH[H^(p), H^(q)]. It comes from the (T) Lagrangian density along H + a H^(q), evaluated in the
  semicanonical orbitals of the perturbed Fock matrix (the (T) functional is invariant to
  occupied–occupied and virtual–virtual rotations; degenerate orbitals need no special care), in
  one directional derivative together with the amplitude part. The (T) energy is a homogeneous
  quadratic function of the amplitudes: the (T) pair terms come from its bilinear form, all pairs
  in one pass over the triples.
- **Dipole derivatives**: the three components of a uniform field are further perturbations
  (CPHF with the dipole integrals, h~^F = r, no overlap or two-electron derivative); the mixed
  derivatives d²E/dq dF give dμ/dq (Full CI from the existing CI responses, CC with three more
  amplitude responses), completed with the translations and rotations; the relaxed dipole moment
  comes with them.

Verification (`backend/tests/test_correlated_hessian.cpp`, against 4-point differences, step 0.005
bohr, of analytic gradients and their relaxed dipoles; aug-cc-pVDZ):

| | Hessian (Eₕ/bohr²) | dipole derivatives (e) |
|---|---|---|
| Full CI: H₂, HeH⁺, H₃⁺ (C1), LiH | ≤ 1.4e-9 | ≤ 2.3e-10 |
| CCSD: LiH, H₂O, HOF (C1) | 1.2e-10, 9.3e-9, 1.6e-9 | 3.8e-10, 6.5e-8, – |
| CCSD(T): LiH, H₂O, HOF (C1), NH₃ (C₃ᵥ) | 1.2e-10, 1.1e-8, 1.5e-9, 4.8e-9 | 8.3e-10, 6.6e-8, –, – |

Without the (T) curvature term the CCSD(T) errors would be 6e-5 (LiH) and 2e-2 (H₂O). The
relaxed dipole moments agree with those of the analytic gradients within 2e-10 e bohr.

Timings of the frequency step (4 threads, including the energy and gradient at the geometry),
analytic vs the numerical route of earlier versions (2-point differences of analytic gradients
along the internal directions): Full CI LiH/aug-cc-pVDZ 6.0 s vs 4.1 s (2 gradients);
H₂O/aug-cc-pVDZ CCSD 13 s vs 8 s, CCSD(T) 21 s vs 11 s (6 gradients). The analytic Hessian also
solves 3 field responses for the IR intensities. For molecules this small the differences are
cheaper but about 10⁴ times less accurate; the cost of the analytic route grows with 3N − 3
responses instead of 2(3N − 6) gradients.

### Geometry optimisation (version 5.0)

`internal_coordinates.cpp`, `optimizer.cpp`: minimisation in redundant internal coordinates
(Pulay & Fogarasi, J. Chem. Phys. 96, 2856 (1992); Peng, Ayala, Schlegel & Frisch, J. Comput.
Chem. 17, 49 (1996)).

- **Coordinates**: bonds between atoms closer than 1.3 × the sum of their covalent radii (separate
  fragments are joined through their closest atoms), bends at every atom with two bonds, and two
  linear-bend components in perpendicular planes when an angle exceeds 175° (back to a bend below
  170°). Wilson B matrix (analytic for stretches and bends), G = BBᵀ with its generalised inverse,
  internal gradient g_q = G⁻ B g_x, projector P = G G⁻.
- **Step**: rational-function optimisation (lowest eigenvector of the augmented Hessian of the
  projected gradient and Hessian P H P + 1000 (1 − P)), scaled to the trust radius (0.3 at the
  start, adapted between 0.001 and 0.5 from the ratio of the actual to the predicted energy
  change); iterative back-transformation x ← x + Bᵀ G⁻ (q_target − q(x)) to 1e-10 bohr.
- **Hessian**: initial Hessian analytic (the RHF Hessian of stage 1 in internal coordinates,
  including the gradient term Σ g_q ∂²q/∂x², when the method is RHF and the basis has at most 100
  functions; 10.0: also the analytic CCSD, CCSD(T) and Full CI Hessians, without dipole
  derivatives, where they fit into the memory limit), otherwise Lindh's model Hessian (Chem. Phys. Lett. 241, 423 (1995)); Schlegel's is
  available. Negative eigenvalues of an analytic initial Hessian are replaced by their magnitude;
  if the start is a saddle point with no gradient along the negative-curvature direction (e.g. H₂O
  at exactly 180°), the first step follows that direction. BFGS updates (skipped, and reported,
  when the curvature condition fails).
- **Convergence** (as in Gaussian): maximum and RMS internal force (normal 4.5e-4 and 3.0e-4
  Eₕ/bohr, tight 1.5e-5 and 1.0e-5) and maximum and RMS displacement of the next step (normal
  1.8e-3 and 1.2e-3 bohr, tight 6.0e-5 and 4.0e-5), plus the energy change since the previous
  step (normal 1e-6, tight 1e-8 Eₕ); converged when all are met, or when both forces are 100 times
  below their thresholds. Each step reports which criteria are met. Default at most 50 steps; an
  optimisation that does not converge stops with a message, and no frequency or single-point step
  follows (jobs 4 and 5).
- **Symmetry**: H₂O is optimised in all three internal coordinates without constraints; the report
  says whether the C₂ᵥ start stays C₂ᵥ (|R(O–H1) − R(O–H2)|, 1e-14 bohr in practice). A change of the
  number of removed near-linear dependencies along the path is reported (the energy surface is
  then not continuous).

Verification (`backend/tests/test_optimizer.cpp`; references from
`reference/pyscf_optimization_reference.py` → `backend/tests/data/optimization_reference.json`):

- Tight optimisations from the default geometries against PySCF analytic gradients minimised by
  geomeTRIC (an independent optimiser): H₂O RHF/aug-cc-pVDZ and aug-cc-pVTZ, CCSD/aug-cc-pVDZ,
  CCSD(T)/aug-cc-pVDZ; LiH RHF, CCSD(T) and Full CI/aug-cc-pVDZ; H₂ Full CI/aug-cc-pVTZ. Bond
  lengths agree within 8e-6 Å, angles within 6e-4°, energies within 1e-10 Eₕ. (PySCF's own
  gradient scanner for CCSD(T) passes the CCSD Λ multipliers to the (T) gradient and returns the
  CCSD energy; the reference script uses a scanner with the CCSD(T) Λ equations instead — with
  PySCF's scanner the H₂O CCSD(T) optimum is 1.9e-3 Å and 0.12° off and 6.7e-6 Eₕ higher.)
- Difficult starts reach the same RHF minimum: an asymmetric start (1.05 Å, 0.90 Å, 120°) ends
  C₂ᵥ; 179°; exactly 180° (a saddle point, left along the bend, C₂ᵥ kept); 0.8 Å and 70°.
- Wilson B matrices (also linear bends) against finite differences to 1e-10; back-transformation
  residual below 1e-15; a CCSDT optimisation of H₂ with numerical gradients reaches the analytic
  Full CI optimum within 7e-6 Å; the step limit stops with a message.

### Harmonic frequencies and thermochemistry (version 5.0)

`hessian.cpp`, `frequencies.cpp`, `constants.cpp`.

- **Hessian**: RHF analytic (stage 1: second-derivative integrals and CPHF, with analytic dipole
  derivatives). 10.0: CCSD, CCSD(T) and Full CI analytic as well (see
  [Analytic Hessians](#analytic-hessians-version-100)); with `hessian = numerical`, as in
  versions 5–9, central differences of analytic gradients along internal
  coordinates (step 0.005 bohr / rad); for C₂ᵥ H₂O the symmetry-adapted S1 (A1, symmetric stretch),
  S2 (A1, bend) and S3 (B2, antisymmetric stretch), 6 gradients; for H₂ and LiH the bond, 2
  gradients. The Cartesian Hessian is Bᵀ H_q B + Σ g_q ∂²q/∂x², exact also away from a stationary
  point. The displaced geometries are superimposed on the reference (mass-weighted, Horn's
  quaternion method, i.e. the Eckart frame), so the relaxed dipole moments of the gradients give
  the dipole derivatives. CCSDT and CCSDTQ: second differences of energies along the same
  coordinates (mixed A1 × B2 terms vanish by symmetry; C₂ᵥ H₂O 9 energies, diatomics 2), dipole
  moments by finite fields μ = −dE/dF at the displaced geometries (10 and 4 more energies). With
  removed near-linear dependencies, RHF also differentiates its analytic gradients numerically.
  The output and the app say which was used. The SCF accepts a uniform electric field for this
  (and for the tests of the relaxed dipoles).
- **Normal modes**: masses of the most abundant isotopes (AME 2020: ¹H 1.007825031898 u,
  ⁷Li 7.016003434 u, ¹⁶O 15.994914619257 u); mass-weighted Hessian in the orthonormal complement of
  the translations and of the rotations about the principal axes through the centre of mass
  (Eckart conditions; 5 for linear molecules); harmonic wavenumbers (negative = imaginary),
  normalised Cartesian displacements, reduced masses, force constants (mdyn/Å), IR intensities
  (km/mol) from dμ/dQ, symmetry species for C₂ᵥ H₂O (A1, B2) and diatomics; the residual
  curvature of the translations and rotations is printed as a check.
- **Thermochemistry** (CODATA 2018): ideal-gas translation at the chosen pressure (Sackur–Tetrode),
  rigid rotor (linear for H₂ and LiH, nonlinear for H₂O; σ(H₂) = 2, σ(LiH) = 1, σ(H₂O) = 2), harmonic
  oscillators with the zero-point energy, closed-shell singlet (q_el = 1). Reported: zero-point
  energy; thermal corrections and totals of U, H = U + pV (pV = RT per mole of ideal gas), S (total and its electronic,
  translational, rotational and vibrational parts), A = U − TS, G = H − TS; C_V; moments of
  inertia, rotational constants (GHz, cm⁻¹), rotational and vibrational temperatures; in hartree,
  kcal/mol, kJ/mol and cal/(mol K) (also J/(mol K)). Imaginary frequencies and those below 10 cm⁻¹
  are left out of the vibrational sums, with a warning.

Verification (`backend/tests/test_frequencies.cpp`; `reference/pyscf_frequency_reference.py` →
`backend/tests/data/frequency_reference.json`):

- RHF wavenumbers, reduced masses and force constants against PySCF's analytic Hessian and
  harmonic analysis (given the same AME 2020 masses; by default PySCF uses standard atomic
  weights such as 1.008 for H, which shifts the wavenumbers of H₂O by 0.4 cm⁻¹) at the
  PySCF/geomeTRIC geometries of H₂O (aug-cc-pVDZ, aug-cc-pVTZ), LiH and H₂: within 2e-3 cm⁻¹,
  1e-6 u and 1e-5 mdyn/Å.
- Thermochemistry against a hand evaluation of the formulas with CODATA 2018 (in the reference
  script, independent of both programs): ZPE, U, H, G, A within 1e-9 Eₕ, entropies and C_V within
  5e-5 J/(mol K); against PySCF's `thermo` (older CODATA constants) within 3e-6 Eₕ and
  5e-3 J/(mol K).
- Numerical Hessians: finite differences of analytic RHF gradients reproduce the analytic
  wavenumbers within 0.05 cm⁻¹ and the IR intensities within 6e-4 km/mol; for CCSD(T) (H₂O) and
  Full CI (LiH) the gradient-based and the energy-based Hessians agree within 0.03 cm⁻¹, and the IR
  intensities from relaxed dipoles and from finite-field dipoles within 0.007 km/mol; relaxed
  dipole moments equal −dE/dF within 1e-9 e·bohr.
- At linear H₂O the degenerate bend is imaginary (twice) and left out of the thermochemistry with
  two warnings.

### Job types and the output file (version 5.0)

`jobs.cpp`, `output.cpp`. Five job types for H₂, LiH and H₂O (HeH⁺ keeps its 4.0 energy path):

| Job | What runs | Geometry of each step |
|---|---|---|
| 1 `energy` | single-point energy, as in 4.0 | input |
| 2 `opt` | geometry optimisation | input → optimised |
| 3 `freq` | Hessian, normal modes, IR intensities, thermochemistry; a warning if the input geometry is not stationary (largest gradient component above 4.5e-4 Eₕ/bohr) | input |
| 4 `opt+freq` | optimisation, then frequencies | frequencies at the optimised geometry |
| 5 `opt+freq+sp` | optimisation and frequencies with `method`/`basis`, then a single-point energy with `sp_method`/`sp_basis`; final U, H, G, A = single-point energy + thermal corrections of step 2 | optimised |

If the optimisation of job 4 or 5 does not converge within `max_steps`, the job stops after it
with an error termination that says why ("the frequency and single-point steps were not run:
they need a stationary point"). Every method can be used where its energy is possible: analytic
gradients for RHF, CCSD, CCSD(T) and Full CI; numerical gradients (central differences of
energies along internal coordinates) for CCSDT and CCSDTQ; analytic RHF Hessians; numerical
Hessians from analytic gradients (CCSD, CCSD(T), Full CI) or from energies (CCSDT, CCSDTQ). The
estimate lists each step with how it is computed, the number of energies, gradients or steps,
the time and the memory; jobs that would take more than three days are refused, and so are
CCSDTQ optimisations and frequencies of H₂O (no analytic CCSDTQ gradient: each step would need
CCSDTQ energies of hours each and more than 5 GB; the CCSDTQ energy itself remains possible with
enough memory, also as the single point of job 5). A job whose single point is impossible is
refused before anything runs.

The **output file** is written while the job runs and is identical in the app's *Output* view, in
*File ▸ Export Results as Text* and in `psilab`'s output file. In order:

1. header: program and version, start time, job type, threads; the input (keys as understood, for
   `psilab` also the input file as given);
2. molecule: electrons, charge, Cartesian coordinates in Å and bohr, bond lengths and angle,
   nuclear repulsion; basis set per element (primitives → contraction, functions), number of
   functions and shells, storage of the two-electron integrals;
3. estimates: size of Full CI or coupled cluster, each step of the job, total time, warnings or the
   refusal;
4. every energy calculation as it happens: SCF iterations (energy, change, orbital gradient,
   DIIS); CC iterations (energy, correlation energy, change, largest amplitude change) or Davidson
   iterations; Lambda iterations; CPHF / Z-vector iterations;
5. at the energy job's geometry or the optimised geometry: energy components (kinetic, nuclear
   attraction, electron repulsion, virial ratio), all orbital energies with HOMO, LUMO and gap,
   Mulliken and Löwdin charges, dipole moment (SCF density; relaxed density for correlated
   gradients), the coupled-cluster or Full CI report (T1 diagnostic, leading amplitudes or
   determinants, natural occupations);
6. optimisation, every step: the Cartesian geometry it uses, the energy and its change, the
   Cartesian gradient and how it was computed, the table of internal coordinates (name,
   definition, value in bohr/rad and Å/deg, force, step), the Hessian (initial or update), the
   RFO step (λ, length, trust radius, predicted change, ratio of actual to predicted change), and
   the convergence table (item, value, threshold, converged?); then "Stationary point found", a
   step summary, the optimised geometry and parameters (start in parentheses) and whether H₂O kept
   C₂ᵥ;
7. frequencies: the geometry used, how the Hessian and the dipole derivatives were obtained (and
   each displaced calculation of a numerical Hessian), isotopes and masses, the residual curvature
   of translations and rotations, the table of modes (symmetry, wavenumber, reduced mass, force
   constant, IR intensity) and the normalised Cartesian displacements;
8. thermochemistry: T, p, mass, symmetry number, moments of inertia, rotational constants and
   temperatures, vibrational temperatures; ZPE and the thermal corrections to U, H, A, G in Eₕ,
   kcal/mol and kJ/mol; E + ZPE, U, H, A, G; E, C_V and S of the electronic, translational,
   rotational and vibrational parts in kcal/mol, cal/(mol K) and J/(mol K), with partition
   functions;
9. single point of job 5 with its own SCF and correlated iterations;
10. summary (final geometry, final energy and its label, convergence, wavenumbers, thermochemistry);
11. (13.0) the **basis set in full** — for every atom (atoms with the same functions as an earlier
    one of their element are named once) each shell with the numbers of its spherical functions,
    the exponents (bohr⁻²) and the contraction coefficients of normalised primitives (every
    contracted function normalised) — and the **molecular orbitals of the Orbitals tab**: the RHF
    canonical orbitals (all occupied ones and ten virtual ones) and, after Full CI, its natural
    orbitals, five per block with number, HOMO/LUMO label, energy and occupation, and the AO
    coefficient of every basis function (labelled by atom, shell and m, real solid harmonics
    m = −l … l; 14.0: Cartesian functions by their component, e.g. `O1 3dxy`);
12. (13.0) the **archive**: the complete results — everything the app shows — as JSON (the result of
    the C API without its "output", plus "archive": format, program, time), ASCII with non-ASCII
    characters escaped, 100 columns per line between the lines `@@PSILAB-ARCHIVE-BEGIN 1` and
    `@@PSILAB-ARCHIVE-END`; joined without the line breaks it equals the JSON (tested). The app reads
    it with *File ▸ Open Output File…*; it is written for every job that was not cancelled (also
    for a refusal or an error, with what exists), not for the points of the app's potential curve,
    whose output is not kept;
13. "Normal termination" with the total time, or "Error termination" with the reason
    (refusal, non-converged optimisation, failed SCF) or "Calculation cancelled".

`psilab`'s output files of HeH⁺ in STO-nG (and the app's text export of them) end with the archive
too. [`examples/outputs/h2o-ccsd-opt-freq.out`](examples/outputs/h2o-ccsd-opt-freq.out) is an
example (H₂O CCSD/aug-cc-pVDZ optimisation + frequencies, `psilab` 13.0, from
[`h2o-ccsd-opt-freq.inp`](examples/outputs/h2o-ccsd-opt-freq.inp)).

Warnings (non-stationary geometry, removed near-linear dependencies, imaginary frequencies left out
of the thermochemistry, rises of the energy during the optimisation) appear where they occur and
again in the summary.

Verification (`backend/tests/test_jobs.cpp`, `backend/tests/test_psilab.cpp`):

- every job type with RHF for H₂ (aug-cc-pVTZ), LiH and H₂O (aug-cc-pVDZ), tight convergence,
  through the C API: optimised bond lengths and angle against PySCF + geomeTRIC within
  2e-5 Å and 2e-3°, energies within 1e-9 Eₕ; frequencies against PySCF's analytic Hessian at the
  same geometry within 0.01 cm⁻¹, at ΨLab's own optimised geometry within 0.05 cm⁻¹; ZPE, S and G
  against PySCF's `thermo`;
- job 5 (H₂O, RHF/aug-cc-pVDZ geometry and frequencies, CCSD(T)/aug-cc-pVTZ single point): the
  single point equals an energy job at the optimised geometry (8e-10 Eₕ; the energy job converges
  the amplitudes to 1e-8, the job to 1e-9), and E + ZPE,
  U, H, G, A equal the single-point energy plus the thermal corrections (1e-10 Eₕ);
- correlated methods (aug-cc-pVDZ): `opt+freq` with CCSD, CCSD(T), CCSDT, CCSDTQ and Full CI for
  H₂, where CCSD = CCSDT = CCSDTQ = Full CI exactly (two electrons): R_e, the wavenumber and the
  energy agree (analytic CCSD against analytic Full CI within 1e-6 Å and 0.05 cm⁻¹; numerical
  CCSDT and CCSDTQ within 5e-5 Å and 0.5 cm⁻¹); LiH `opt+freq` with CCSD, CCSD(T) and Full CI (in the
  *full* CI job also CCSDT and CCSDTQ, about 90 s per energy, with CCSDTQ = Full CI for four
  electrons); LiH CCSD(T) and Full CI and H₂O
  CCSD and CCSD(T) geometries against PySCF + geomeTRIC;
- refusals (CCSDTQ frequencies and optimisation of H₂O, Full CI of H₂O, an impossible single
  point in job 5), the counts of the numerical Hessians in the estimates (H₂O CCSDT: 9 energies +
  10 finite-field energies; CCSD(T): 6 gradients; diatomics: 2);
- the streamed output: read with `macqc_job_output` while the job runs, it is identical to the
  final `output`; progress text; cancellation during an optimisation step (the run returns
  "calculation cancelled" and the output says so);
- `psilab`: all new keys and their errors, the streamed text equals the output file.

### Relaxed scans, transition states, H2O2 and NH3 (version 7.0)

**Relaxed scans** (job `scan`) follow Gaussian's `Opt=ModRedundant`: the line
`N1 N2 [N3 [N4]] S nsteps stepsize` (two atoms a bond, three an angle, four a dihedral angle; the
letters B, A, D may precede them) gives nsteps + 1 points, the first at the input value or at
`scan_start`. At each point the coordinate is frozen and the others are optimised, starting from
the previous point's geometry, orbitals and BFGS-updated Hessian. Frozen coordinates use the
projector P̃ = P − PC(CᵀPC)⁻¹CᵀP (Peng, Ayala, Schlegel and Frisch, *J. Comput. Chem.* 17, 49
(1996)) for the RFO step and the convergence test, and a correction PC(CᵀPC)⁻¹(target − q) drives the
coordinate to its value; an angle frozen at 180° is imposed through its two linear-bend
components. The redundant internal coordinates now include torsions (IUPAC sign, Blondel–Karplus
B rows, differences reduced to (−π, π]) and improper torsions at atoms with three bonds (the
umbrella of NH₃). The summary table gives E, ΔE in kcal/mol and kJ/mol, the steps per point and the
other parameters along the scan, and minima and maxima from parabolas through three points.
Diatomics have nothing else to relax: their scan is a series of single points.

**Transition states** (jobs `ts`, `ts+freq`, `ts+freq+sp`) follow Gaussian's `Opt=TS`: the Hessian
is computed at the first point (CalcFC) or at every point (CalcAll) — analytic for RHF, central
differences of analytic gradients for CCSD, CCSD(T) and Full CI — and transformed to internal
coordinates as H_q = G⁻B(H_x − K)BᵀG⁻; between computed Hessians Bofill updates (Murtagh–Sargent and
Powell mixed with φ = (ξ·s)²/((ξ·ξ)(s·s))). The partitioned RFO step (Baker, *J. Comput. Chem.* 7,
385 (1986)) maximises along the followed mode (the lowest at the first step, then the one with the
largest overlap with the previous) and minimises along all others; the trust radius stays at 0.30
bohr/rad. The eigenvalue test stops the search unless the Hessian has exactly one negative
eigenvalue (`eigen_test = no`: Gaussian's NoEigenTest). `ts+freq` ends with an error unless the
frequencies give exactly one imaginary mode; `ts+freq+sp` adds a single point at the transition
state, and its U, H, G and A are that energy plus the thermal corrections.

**H₂O₂ and NH₃**: H₂O₂ (O1, O2, H3 on O1, H4 on O2; C₂ axis along z) from R(O–O), R(O–H), the
O–O–H angle and the H–O–O–H dihedral angle (defaults 1.4556 Å, 0.9670 Å, 102.32°, 113.70°); NH₃
(N1, H2–H4; C₃ᵥ) from R(N–H) and the H–N–H angle (1.012 Å, 106.7°; 120° gives the planar D₃ₕ
molecule). Nitrogen is in all eight basis files (BSE 0.12, the same data as the other elements).
Point groups C₁, Cₛ, Cᵢ, C₂, C₂ᵥ, C₂ₕ, C₃, C₃ᵥ and D₃ₕ are found from the symmetry operations (to
2·10⁻³ bohr, tightened if needed), the species of the normal modes from the characters, and the
symmetry number of the rotational partition function is the number of proper rotations.

**HeH⁺** with a correlation-consistent basis set runs the general calculation: every method and
job type except the TS search, which is refused for all diatomics. He has no core-valence sets;
aug-cc-pCVXZ falls back to aug-cc-pVXZ with a note.

Verification (all against independent PySCF/geomeTRIC calculations, in CI):

| Check | Result |
|---|---|
| H₂O RHF scans, R(O–H) frozen at 0.85 and 1.20 Å; angle 90° … 180° | ≤ 2.1·10⁻¹⁰ Eₕ, ≤ 9.4·10⁻⁷ Å, ≤ 5·10⁻⁵° |
| H₂O CCSD(T) scan point, R(O–H) frozen at 1.10 Å | 9.6·10⁻¹¹ Eₕ, 9.3·10⁻⁷ Å |
| H₂O₂ torsion scan 0–180° (RHF), NH₃ umbrella to the planar molecule | ≤ 10⁻⁹ Eₕ |
| RHF energies and gradients of H₂O₂ and NH₃ (aug-cc-pVDZ, aug-cc-pVTZ, aug-cc-pCVDZ) | ≤ 2·10⁻¹² Eₕ, ≤ 8·10⁻¹¹ Eₕ/bohr |
| Tight optimisations of H₂O₂ and NH₃ (RHF, CCSD(T)) vs geomeTRIC | ≤ 5·10⁻⁷ Å, 7·10⁻⁴°, 10⁻¹¹ Eₕ |
| Transition states vs geomeTRIC `transition=True`: trans and cis H₂O₂, planar NH₃ (RHF), NH₃ CCSD(T) | ≤ 3·10⁻¹⁰ Eₕ, ≤ 2·10⁻⁵ Å |
| Imaginary frequencies vs PySCF's analytic Hessian: 312.41i (Aᵤ, C₂ₕ), 626.74i (A₂, C₂ᵥ), 871.98i cm⁻¹ (A₂″, D₃ₕ) | ≤ 0.13 cm⁻¹ |

The CCSD(T)/aug-cc-pVDZ inversion barrier of NH₃ (planar transition state minus pyramidal minimum)
is 5.77 kcal/mol. The scripts are `reference/pyscf_scan_reference.py`,
`pyscf_molecules_reference.py` and `pyscf_ts_reference.py`; the test programs `macqc_scan_tests`,
`macqc_molecule_tests` and `macqc_ts_tests`.

### Molecular orbitals (version 8.0)

**Molecular orbitals**: the result contains the basis and the canonical RHF orbitals (all occupied
ones and up to ten virtual ones) and, after Full CI, the natural orbitals (eigenvectors of the Full
CI one-particle density, in AO coefficients C_RHF·U). The backend evaluates an orbital on a grid —
spherical shells as Σ_m M(m, c) x^a y^b z^c Σ_p c_p e^(−α_p r²), the same solid harmonics as the
integrals (14.0: Cartesian shells as N_abc x^a y^b z^c Σ_p c_p e^(−α_p r²)) — and extracts the isosurfaces ψ = ±c by marching tetrahedra (six tetrahedra around the
main diagonal of each cube, so neighbouring cubes split their faces alike; shared vertices). A
coarse grid (0.5 bohr) finds the region |ψ| ≥ c/2, a grid of at most 400 000 points covers it with a
margin of 1 bohr, and every vertex is then moved onto the exact isovalue by Newton steps along its
normal, with the normal recomputed from the gradient of the orbital.

The app (*Methods ▸ Orbitals*) shows the orbitals of the result as an energy-level diagram and a
table, and the chosen one as isosurfaces (SceneKit, the molecule as balls and sticks; seen from 30°
above the plane of the two longest principal axes, the third axis up) or as a contour plot at
±c·2ᵏ in the molecular plane or a plane of two principal axes (non-planar molecules start in a plane
that contains their unique axis; for NH₃ the σᵥ plane through N1 and H2).

Verification (against PySCF, in CI; `reference/pyscf_orbitals_reference.py`, `macqc_orbital_tests`).
The comparisons do not depend on how a program orders or signs its harmonics, or on the phases of
the orbitals:

| Check | Result |
|---|---|
| Basis functions l = 0 … 5 (O, Li, He, H aug-cc-pV5Z; N aug-cc-pCV5Z): Σ_m φ_lm(r) φ_lm(r′) | ≤ 1.1·10⁻¹⁵ (relative) |
| RHF orbitals of H₂O, NH₃, H₂O₂ at 24 points (up to the sign; NH₃ e pairs by Σψ²); density | ≤ 1.6·10⁻⁷; ≤ 6·10⁻¹⁰ (relative) |
| Full CI natural orbitals of LiH: occupations; density Σ n_k χ_k²; Cᵀ S C = 1 | 1.6·10⁻⁸; 7.5·10⁻⁹; 4·10⁻¹⁴ |
| Isosurfaces: |ψ| at the vertices of the H₂O HOMO after the Newton steps; a Gaussian sphere closed, oriented, volume | 1.5·10⁻⁷ (relative); 0.15 % |

### Layout

```
backend/include/macqc/   public headers (macqc_api.h is the C API)
backend/src/             linalg, integrals, scf, fci, report, calculation, api (HeH+, STO-nG);
                         molecule, basis, gto_integrals, ecp (14.0), rhf, molecular (H2, LiH, H2O);
                         gemm, mo_integrals, ci (Full CI for H2, LiH); job (cancel/progress);
                         cc_closed_shell (CCSD, (T)), cc_general (CCSDT, CCSDTQ), cc_estimate;
                         gto_derivatives (derivative integrals), rhf_derivatives (RHF gradient,
                         Hessian, CPHF, dipole derivatives); tensor (einsum), cc_gradient (CCSD
                         and CCSD(T) Lambda equations and densities), correlated_gradient
                         (relaxed densities, Z-vector), gradient (energies and gradients of
                         every method at any geometry); internal_coordinates, optimizer
                         (redundant internal coordinates, RFO, BFGS); hessian (analytic and
                         numerical Hessians, dipole derivatives), frequencies (normal modes,
                         thermochemistry), constants (CODATA 2018, AME 2020); jobs (job types,
                         relaxed scans, output file), symmetry (point groups, mode species),
                         orbitals (orbital sets, grids, isosurfaces)
backend/tools/           cc_diagrams.py: generates and verifies the coupled-cluster diagrams
basis/                   Gaussian94 basis-set files from basissetexchange.org (+ README)
backend/cli/main.cpp     developer command-line tool `macqc` (flags / JSON)
backend/cli/psilab*      terminal version `psilab` (text input file -> text output file)
examples/                example input files for `psilab` (heh.inp, h2o.inp)
backend/tests/           test programs (CTest; 14.0: test_basis_sets.cpp), reference data in tests/data/
reference/               the verified Python reference (HeH+); PySCF reference scripts (H2, LiH, H2O,
                         H2O2, NH3; scans, transition states, orbitals; 14.0 Pople and def2 sets)
third_party/nlohmann/    json.hpp (MIT)
CMakeLists.txt           library + CLI + tests (for testing; Xcode compiles the sources directly)
```

### Build, test and run

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

```bash
./build/macqc                     # text report (default: HeH+, STO-3G, print level 2)
./build/macqc --n 2 --R 1.6 --print 1
./build/macqc --json              # JSON instead of the text report
./build/macqc --input in.json     # read the API input JSON ("-" = stdin); flags override it
./build/macqc --help
```

Options: `--R --zeta1 --zeta2 --ZA --ZB --n --exact-erf --max-iter --conv --print --json --input`.
General molecules go through `--input`, e.g.
`echo '{"molecule":"H2O","basis":"aug-cc-pVTZ"}' | ./build/macqc --input -`.
The exit status is 0 on success and 1 (message on stderr) if the input is invalid or the SCF did
not converge.

The expected HF reports in `backend/tests/data/` are generated by the reference script
(`python3 backend/tests/gen_expected.py`, needs NumPy) and must never be edited by hand; the C++
text report is compared with them byte for byte (print levels 0, 1, 2; STO-1G/2G/3G; R = 2.5).

### Text report

The HF part reproduces the output of `reference/hartree_fock_sto3g.py` exactly (headings, spacing,
MATOUT layout, Fortran `D` numbers). The Full CI section is appended after it:
determinants, (level 2: MO `h`, CI Hamiltonian and MO integrals), all eigenvalues with ⟨S²⟩ and
spin label, the ground-state CI vector, singlet weights, FCI energies and the correlation energy.
Print level 0 is silent (empty report, use the JSON fields), level 1 is the summary, level 2 prints
every SCF iteration.

### C API

```c
char* macqc_run(const char* input_json);   // never throws; returns malloc'ed JSON
void  macqc_free(char* p);                 // release the string

// Jobs: cancellation and progress of long calculations (thread-safe).
macqc_job* macqc_job_create(void);
char* macqc_run_job(const char* input_json, macqc_job* job);   // like macqc_run
void  macqc_job_cancel(macqc_job* job);    // from another thread; the run then returns
                                           // "error": "calculation cancelled"
void  macqc_job_progress(const macqc_job* job, int* stage, int* iteration,
                         double* energy, double* residual);
      // stage 0 idle, 1 integrals, 2 SCF, 3 integral transformation, 4 Full CI, 5 done,
      // 6 coupled-cluster iterations (residual = largest amplitude change), 7 (T) correction
      // (iteration = percent done); 5.0: 8 gradient, 9 Hessian, 10 CPHF, 11 Lambda equations,
      // 12 geometry optimisation, 13 numerical Hessian
char* macqc_job_output(const macqc_job* job, long long offset, long long* total);
      // 5.0: the output file written so far, from byte `offset` on (free with macqc_free);
      // *total = its full length, so a reader polls with offset = the previous total
char* macqc_job_progress_text(const macqc_job* job);
      // 5.0: one line about the current step (free with macqc_free), e.g.
      // "Optimisation step 4 · E = -76.0418435426 E_h · max force 1.10e-04",
      // "Hessian: displacement 3 of 6 (S2 -)", "Orbital response iteration 7 · residual 3.1e-06"
void  macqc_job_free(macqc_job* job);      // after macqc_run_job has returned
```

The app polls `macqc_job_progress` four times a second for the progress line of the sidebar.
The output buffer of a job is a string guarded by a mutex inside the job handle: the calculation
appends each piece of the output file as it is produced and `macqc_job_output` copies the part a
reader has not seen yet, so the app's *Output* view follows the job without blocking it.

#### Input JSON

```json
{"R":1.4632,"zeta1":2.0925,"zeta2":1.24,"ZA":2,"ZB":1,"nGaussians":3,
 "exactErf":false,"maxIter":25,"convergence":1e-4,"printLevel":2}
```

All keys are optional (the values above are the defaults). Validation: `R`, `zeta1`, `zeta2`,
`ZA`, `ZB`, `convergence` > 0; `nGaussians` ∈ {1,2,3}; `maxIter` ≥ 1; `printLevel` ∈ {0,1,2}.
Invalid or malformed input yields an error JSON; the call never crashes.

General molecules are selected with `"molecule"` (`"HeH+"` is the default and keeps everything
above unchanged):

```json
{"molecule":"H2O","basis":"aug-cc-pVTZ","R":1.808846,"angle":104.52,"method":"RHF",
 "basisDirectory":"/path/to/basis","threads":0,"memoryMB":1024,"fciMemoryMB":4096,
 "scfMaxIter":100,"energyConvergence":1e-10,"gradientConvergence":1e-7,"estimateOnly":false,
 "printLevel":2}
```

`R` is in bohr (H₂O: O–H) and defaults to the experimental equilibrium bond length; `angle`
(degrees, H₂O) defaults to 104.52; `method` is `"RHF"` (also `"HF"`), `"CCSD"`, `"CCSD(T)"`,
`"CCSDT"`, `"CCSDTQ"` or `"FCI"` (the correlated methods follow RHF); `fciMemoryMB` is the
memory limit of Full CI and coupled cluster;
`basisDirectory` defaults to `$PSILAB_BASIS_DIR` or the compiled-in source folder. With
`"estimateOnly":true` only the molecule, the basis set, the size and time estimate, the warnings
and (for a Full CI or coupled-cluster calculation that is not possible) the refusal are
returned, without an error.

Version 5.0 adds the job keys (H₂, LiH, H₂O; HeH⁺ only runs the energy job):

```json
{"job":"opt+freq+sp","spMethod":"CCSD(T)","spBasis":"aug-cc-pVTZ",
 "temperature":298.15,"pressure":101325,"optConvergence":"normal","maxSteps":50}
```

`job` is `"energy"` (default), `"opt"`, `"freq"`, `"opt+freq"` or `"opt+freq+sp"` (other
spellings such as `"OPT FREQ"` or `"5"` are accepted); `spMethod`/`spBasis` default to `method`
and `basis`; `temperature` in K (0 < T ≤ 10⁵), `pressure` in Pa; `optConvergence` `"normal"` or
`"tight"`; `maxSteps` 1–1000. With `"estimateOnly":true` the result's `job.estimate` lists the
steps of the job (how each is computed, how many energies, gradients or steps, time, memory) and
the refusal, if any.

Version 7.0 adds the molecules `"H2O2"` and `"NH3"`, HeH⁺ in the correlation-consistent basis sets,
the scan and the transition-state jobs:

```json
{"molecule":"H2O2","R":1.82737,"ROO":2.75069,"angle":102.32,"dihedral":113.70,
 "job":"scan","scan":"D 3 1 2 4 S 18 10","scanStart":0,"tsHessian":"calcfc","eigenTest":true}
```

`ROO` (bohr) and `dihedral` (degrees, IUPAC sign) are those of H₂O₂; for H₂O₂ and NH₃ `R` is the
X–H bond (bohr) and `angle` the O–O–H or H–N–H angle. `job` may also be `"scan"` (6), `"ts"` (7,
also `"Opt=TS"`), `"ts+freq"` (8) or `"ts+freq+sp"` (9). `scan` is Gaussian's ModRedundant line
`[B|A|D] N1 N2 [N3 [N4]] S nsteps stepsize` (step in Å or degrees; without the letter the number
of atoms gives the type), `scanStart` the first value in Å or degrees (default: the input
geometry's). `tsHessian` is `"calcfc"` or `"calcall"`, `eigenTest` `true` or `false`. For
`"molecule":"HeH+"`, `"basis"` selects the calculation: a correlation-consistent set the general
one, `"STO-1G"`, `"STO-2G"` or `"STO-3G"` (or no basis) the Szabo–Ostlund one.

#### Output JSON

```text
{
  "input":  { ...the (defaulted) input echoed... },
  "error":  null | "message",
  "hf": {
    "converged", "iterations", "electronicEnergy", "nuclearRepulsion", "totalEnergy",
    "orbitalEnergies": [e1, e2],
    "C", "P", "S", "X", "H", "F", "G", "PS":   2x2 matrices as [[row0],[row1]],
    "integrals": {"S12","T11","T12","T22","V11A","V12A","V22A","V11B","V12B","V22B",
                  "V1111","V2111","V2121","V2211","V2221","V2222"},
    "scfHistory": [{"iteration","electronicEnergy","deltaP"}, ...]
  },
  "fci": {
    "totalEnergy", "electronicEnergy", "correlationEnergy",
    "states": [{"energy","S2","spin","vector"}, ...],      // ascending, total energies
    "groundStateVector": [c1..c4],
    "determinants": ["|1α 1β>", "|1α 2β>", "|2α 1β>", "|2α 2β>"],
    "singletWeights": [{"configuration":"(1sigma)^2","weight":…}, …]   // ground state
  },
  "report": "FORTRAN-style text report"
}
```

On invalid input `hf` and `fci` are `null`. If the SCF does not converge, `hf` is present with
`converged:false`, `fci` is `null` and `error` is set.

For general molecules the output is

```text
{
  "input": {...}, "error": null | "message",
  "molecule": {"name","electrons","coreElectrons" (14.0: in ECPs),"charge","R","angle","nuclearRepulsion",
               "atoms":[{"symbol","Z","position":[x,y,z]}]},                    // bohr
  "basis":    {"name","functions","shells","maxL",
               "cartesianD","functionTypes",                              // 14.0
               "atoms":[{"symbol","basis","contraction","primitives","functions","ecp","ecpElectrons"}]},
  "estimate": {"basisFunctions","uniqueIntegrals","integralMemoryMB","integralsStored",
               "threads","rhfSeconds","totalSeconds","longCalculation",
               "fci": {"orbitals","alphaElectrons","betaElectrons","alphaStrings","betaStrings",
                       "determinants","memoryMB","seconds","feasible","reason"} | null,
               "coupledCluster": {"method","engine","occupied","virtuals","amplitudes","memoryMB",
                       "ladderIntegralsStored","iterations","iterationSeconds","triplesSeconds",
                       "seconds","feasible","reason"} | null},
  "notes":    ["aug-cc-pCVDZ is not defined for H ..."],
  "warnings": ["Long calculation: ... Estimated time about 6 min on 8 threads."],
  "refusal":  null | "Full CI is not possible here: ...",
  "rhf": {"converged","iterations","totalEnergy","electronicEnergy","nuclearRepulsion",
          "kineticEnergy","nuclearAttractionEnergy","electronRepulsionEnergy","virialRatio",
          "orbitalEnergies","occupiedOrbitals","homo","lumo","smallestOverlapEigenvalue",
          "removedFunctions","integralsStored","uniqueIntegrals","threads",
          "integralSeconds","scfSeconds",
          "scfHistory":[{"iteration","totalEnergy","deltaE","gradient","diis"}]} | null,
  "fullCI": {"converged","iterations","totalEnergy","electronicEnergy","correlationEnergy","S2",
             "hfWeight","orbitals","alphaElectrons","betaElectrons","alphaStrings","betaStrings",
             "determinants","naturalOccupations","transformSeconds","ciSeconds",
             "leadingDeterminants":[{"coefficient","alpha":[MOs],"beta":[MOs],"excitation"}],
             "history":[{"iteration","totalEnergy","residual","subspace"}]} | null,
  "coupledCluster": {"method","engine","converged","iterations","occupied","virtuals",
             "referenceEnergy","correlationEnergy","ccEnergy","triplesCorrection" | null,
             "totalEnergy","t1Diagnostic","amplitudeCounts":[by rank],"amplitudeNorms":[by rank],
             "leadingAmplitudes":[{"value","rank","occupied":[MOs],"virtuals":[MOs],"spins"}],
             "history":[{"iteration","totalEnergy","correlationEnergy","deltaE","residual"}],
             "ladder","transformSeconds","ccSeconds","triplesSeconds"} | null,
  "hf": null, "fci": null,
  "report": "text report"
}
```

A refused calculation (e.g. LiH/aug-cc-pCV5Z, Full CI for H₂O, CCSDTQ for H₂O/aug-cc-pVTZ) has
`error` set and `rhf`, `fullCI` and `coupledCluster` `null`; for Full CI and coupled cluster the
reason is also in `refusal`. `correlationEnergy` and `ccEnergy` exclude (T); `totalEnergy`
includes it.

Version 5.0 adds, for H₂, LiH and H₂O (`null` where the job has no such step):

```text
  "job": {"type","title","termination":"normal"|"error"|"refused"|"estimate",
          "completed","message","seconds","warnings",
          "estimate": {"steps":[{"name","method","basis","how","analytic","evaluations",
                                 "evaluationKind":"steps"|"gradients"|"energies"|"Hessians",
                                 "seconds","memoryMB"}],
                       "totalSeconds","memoryMB","longCalculation","refusal","warnings","notes"},
          "initialGeometry","finalGeometry",           // {"atoms","R","R1","R2","angle"}, bohr, deg
          "finalEnergy","finalEnergyLabel",            // e.g. -76.33, "CCSD(T)/aug-cc-pVTZ"
          "singlePointMethod","singlePointBasis"},
  "optimization": {"converged","steps","message","initialHessian","gradient",
          "criteria":{"maxForce","rmsForce","maxDisplacement","rmsDisplacement","energy"},
          "initialEnergy","finalEnergy","initialGeometry","finalGeometry",
          "initialParameters","finalParameters":[{"name","value","angle"}],   // bohr or degrees
          "symmetricStart","symmetricEnd","symmetryDeviation","seconds",
          "history":[{"step","energy","deltaE","geometry",
                      "coordinates":[{"name","definition","value","gradient","step"}],
                      "convergence":[{"name","value","threshold","met","applicable"}],
                      "converged","hessianUpdate","trustRadius","stepNorm","stepScaled",
                      "predictedChange","ratio","warnings","seconds"}]},
  "frequencies": {"analytic","description","numericalReason","gradientEvaluations",
          "energyEvaluations","coordinates","intensities","intensityDescription","dipole",
          "linear","vibrations","masses","isotopes","geometry","energy","lowFrequencies",
          "modes":[{"wavenumber","imaginary","symmetry","reducedMass","forceConstant",
                    "irIntensity","displacement":[[x,y,z] per atom]}],"seconds"},
  "thermochemistry": {"temperature","pressure","symmetryNumber","linear","mass","moments",
          "rotationalConstantsGHz","rotationalConstantsCm","rotationalTemperatures",
          "vibrationalTemperatures","excludedModes","warnings",
          "terms":{"translation"|"rotation"|"vibration"|"electronic"|"total":
                   {"energy","entropy","cv","q"}},                 // J/mol, J/(mol K)
          "zeroPointEnergy","corrections":{"zeroPoint","thermalEnergy","enthalpy","gibbs",
          "helmholtz"},                                             // J/mol
          "electronicEnergy","totals":{...same keys...}},           // E_h
  "singlePoint": {"method","energy","converged","rhfEnergy","correlationEnergy","geometry",
          "warnings","seconds"},
  "output": "the complete output file"
```

`molecule`, `basis`, `rhf`, `fullCI`, `coupledCluster` and `report` then describe the final
energy calculation of the job (the optimised geometry; for `opt+freq+sp` the single point). In
`thermochemistry.totals` the energy is the job's final energy (for `opt+freq+sp` the single point)
plus the thermal corrections of the frequency step.

Version 7.0 adds:

```text
  "scan": {"coordinate":{"type","letter","atoms","definition","unit"},   // "angstrom" | "degrees"
          "start","stepSize","steps","inputValue","relaxed","converged","message","lowest",
          "highest","points":[{"index","value","measured","energy","relativeEnergy","converged",
          "steps","maxForce","geometry","parameters","optimizationSteps":[{"step","energy",
          "maxForce"}],"warnings","seconds"}],"extrema":[{"kind","index","value","energy"}],
          "seconds"},
  "optimization": {..., "transitionState","tsHessian","hessianEvaluations",
          "negativeEigenvalues","eigenTestFailed",
          "history":[{..., "ts":{"hessianEigenvalues","negativeEigenvalues","followedEigenvalue",
                      "followedGradient","followedOverlap","followedMode","lambdaP","lambdaN",
                      "hessianComputed"} | null}]},
  "frequencies": {..., "pointGroup"},                          // "C2", "C2v", "C2h", "C3v", "D3h", …
```

Geometries (`geometry`, `initialGeometry`, …) of H₂O₂ also give `ROO`, `R1`, `R2`, `angle1`,
`angle2` and `dihedral`, those of NH₃ `R1`–`R3`, `angle1`–`angle3` and `pyramidal`
(360° − ΣH–N–H).

Version 8.0 adds the molecular orbitals:

```text
  "orbitals": {"basis","functions","atoms":[{"symbol","Z","position","ecp" (14.0)}],
          "shells":[{"l","atom","center","exponents","coefficients","first","pure"}],
          "sets":[{"kind","title","orbitals":[{"number","energy","occupation","label"}],
                   "coefficients":[[AO coefficients] per orbital]}]}
```

`orbitals` (also for HeH⁺ in STO-nG, then with the two Szabo–Ostlund orbitals;
`null` for an estimate or a refused calculation) belongs to the final calculation of the job: the
shells' `coefficients` include the normalisation of the primitives and of the contraction, the
functions of a shell are the real solid harmonics m = −l … l from `first` on (14.0: with
`"pure": false` the (l+1)(l+2)/2 Cartesian functions xx, xy, xz, yy, yz, zz …, each normalised —
the d shells of the 6-31G family; a shell without `"pure"`, as written by 13.0, is spherical). The C API draws an
orbital from this object:

```c
char* macqc_orbital_surface(const char* orbitals_json, int set, int orbital, double isovalue,
                            int max_points, macqc_mesh* positive, macqc_mesh* negative);
char* macqc_orbital_plane(const char* orbitals_json, int set, int orbital, const double* origin,
                          const double* u, const double* v, int nu, int nv, float* values);
void macqc_mesh_free(macqc_mesh* mesh);
```

`macqc_orbital_surface` returns the isosurfaces ψ = +isovalue and ψ = −isovalue (`macqc_mesh`:
vertex positions in bohr, unit normals pointing out of the lobe, triangles counter-clockwise seen
from outside) on a grid of at most `max_points` points (0: 400 000); `macqc_orbital_plane` the
orbital at origin + i·u + j·v (i < nu fastest). Both return `NULL` or an error message (free it
with `macqc_free`).

### Adding to Xcode

How the backend is wired into the ΨLab Xcode target (the full step-by-step guide for the app is
[`app/SETUP.md`](app/SETUP.md); `xcode/project.yml` encodes the same settings for XcodeGen):

1. Copy or reference this repository's folders in the Xcode project (*File ▸ Add Files to …*,
   "Create folder references" or groups; do **not** tick "Copy items" if you want to keep a single
   source of truth).
2. In the app target's **Build Phases ▸ Compile Sources** add exactly these files:
   `backend/src/linalg.cpp`, `integrals.cpp`, `report.cpp`, `scf.cpp`, `fci.cpp`,
   `molecule.cpp`, `basis.cpp`, `gto_integrals.cpp`, `rhf.cpp`, `molecular.cpp`, `gemm.cpp`,
   `mo_integrals.cpp`, `ci.cpp`, `cc_general.cpp`, `cc_closed_shell.cpp`, `cc_estimate.cpp`,
   `gto_derivatives.cpp`, `rhf_derivatives.cpp`, `tensor.cpp`, `correlated_gradient.cpp`,
   `cc_gradient.cpp`, `gradient.cpp`, `internal_coordinates.cpp`, `optimizer.cpp`, `constants.cpp`,
   `frequencies.cpp`, `hessian.cpp`, `output.cpp`, `jobs.cpp`, `job.cpp`,
   `calculation.cpp`, `api.cpp`; add the folder `basis/`
   as a folder reference to *Copy Bundle Resources* (the app passes its path as `basisDirectory`). Do **not** add `backend/cli/main.cpp` or anything in
   `backend/tests/` (they have their own `main`). Headers (`backend/include/macqc/*.h` and
   `third_party/nlohmann/json.hpp`) need no compile phase but should be in the project.
3. **Build Settings ▸ Header Search Paths** (`HEADER_SEARCH_PATHS`): add
   `$(SRCROOT)/backend/include` and `$(SRCROOT)/third_party` (adjust the prefix to where the folders
   live relative to the `.xcodeproj`).
4. **Build Settings ▸ C++ Language Dialect** (`CLANG_CXX_LANGUAGE_STANDARD`): `C++17`
   (or `GNU++17`). C++ standard library: `libc++` (default).
5. Recommended: **Other C++ Flags** `-ffp-contract=off`. Apple clang contracts `a*b+c` into fused
   multiply-adds on Apple Silicon by default, which can change the last bits of the results
   relative to the reference output (the CMake build and CI use this flag).
6. **Build Settings ▸ Objective-C Bridging Header** (`SWIFT_OBJC_BRIDGING_HEADER`): point to a
   header such as `MacQC-Bridging-Header.h` containing

   ```c
   #include "macqc/macqc_api.h"
   ```

   The header is plain C (`extern "C"` guarded), so Swift imports `macqc_run`/`macqc_free`
   directly; no Objective-C++ wrapper is needed.
7. Call it from Swift:

   ```swift
   import Foundation

   func runMacQC(inputJSON: String) -> Data? {
       guard let raw = macqc_run(inputJSON) else { return nil }   // nil only if out of memory
       defer { macqc_free(raw) }
       return Data(bytes: raw, count: strlen(raw))                // UTF-8 JSON
   }

   // let data = runMacQC(inputJSON: #"{"R":1.4632,"nGaussians":3}"#)
   // let result = try JSONSerialization.jsonObject(with: data!)  // or Codable
   ```

   Run it off the main thread for scans; the library keeps no global state, so concurrent calls
   are safe.

### Custom molecules (version 9.0)

**Input.** JSON `"molecule": "custom"` with `"atoms"` (`[{"symbol": "C", "x": 0, "y": 0, "z": 1.16}, …]`,
Å), `"charge"` (0), `"multiplicity"` (1) and `"title"`; `psilab`: `molecule = custom` and
`geometry = file.xyz` (or .gjf/.com, relative to the input file; `charge`, `multiplicity`, `title`
override the file's; see `examples/co2.inp`). `macqc_parse_structure(text, hint)` reads XYZ files
(several frames: the first one, with a note) and Gaussian inputs (Link 0 and route lines, the title
section, charge and multiplicity, labels such as `C1`, `C(Fragment=1)`, freeze codes, `Units=AU`);
Z-matrices, ghost and dummy atoms are refused with the line number. Every element (Z = 1 … 118) can
be read; masses are those of the most abundant isotope (AME 2020 up to Kr).

**Basis sets.** The eight Gaussian94 files now contain every element the Basis Set Exchange defines
among H–Ar and Sc–Kr (see [`basis/README.md`](basis/README.md)); Sc–Zn stop at aug-cc-pVQZ (their 5Z
sets need i functions). K and Ca have no aug-cc-pVXZ sets; elements beyond Kr need effective core
potentials (aug-cc-pVXZ-PP), which ΨLab does not implement — such molecules are refused with this
explanation. RHF needs a closed-shell singlet: an odd number of electrons or a multiplicity other
than 1 is refused (change the charge).

**General code paths.** The built-in molecules have Z-matrix coordinates for their numerical
derivatives; a custom molecule (atoms in any order) uses the Cartesian directions orthogonal to the
three translations and the rotations about the centroid: 3N − 6 (linear: 3N − 5) directions, so a
numerical gradient takes 2(3N − 6) energies and a numerical Hessian 2(3N − 6) analytic gradients.
The Hessian is completed exactly in the rigid-motion space by invariance — H t = 0 for a translation
and H (n × (x − c)) = n × g for a rotation — so it is correct also away from a stationary point; the
atomic polar tensor is completed likewise (a translation changes μ by the charge times the shift, a
rotation by n × (μ − Q c)). Numerical Hessians from energies (CCSDT, CCSDTQ) are refused for custom
molecules of three or more atoms (use CCSD(T) for the Hessian). The point group is found from the
symmetry operations (C2 … C8, S3 … S8, reflections, inversion; candidate axes from atoms, pairs and
triangles of equivalent atoms and the principal axes): C1 … Oh and C∞v, D∞h, with the rotational
symmetry number of the thermochemistry (e.g. CH₄ 12, C₆H₆ 12, CO₂ 2); species of the normal modes for
C1, Cs, Ci, C2, C2v, C2h, C3, C3v, D3h, D2, D2h (Mulliken's axes), D3, D3d, Td and linear molecules
(Σ, Π with g/u). The potential curve of the app (X–H bonds stretched) is for the built-in molecules;
use the relaxed scan for custom molecules.

**Verification** (`macqc_custom_tests`; `reference/pyscf_custom_reference.py` in CI):

| Check | Result |
|---|---|
| RHF energies of 23 custom molecules against PySCF (CO₂, CH₄, C₂H₄, HCN, BeH₂, BH₃, LiF, NaH, SiH₄, PH₃, H₂S, HCl, HBr, Kr, Zn, ZnH₂, CuH, NH₄⁺, OH⁻, HOCl; HCl and SiH₄ in aug-cc-pCVDZ) | ≤ 2·10⁻¹⁰ Eₕ |
| CCSD(T)/aug-cc-pVDZ of HCl and CH₄ | ≤ 1.2·10⁻⁹ Eₕ |
| Analytic RHF gradient of HOCl | 8·10⁻¹⁰ Eₕ/bohr |
| Numerical gradient (6 energies) and Hessian (6 gradients, not at a minimum) of HOCl against the analytic ones | 1.4·10⁻⁵ Eₕ/bohr; 0.10 cm⁻¹, IR intensities 3·10⁻⁵ (relative) |
| NH₃ with the atoms in another order against the built-in NH₃ (energy, analytic frequencies, species, G) | 3·10⁻¹³ Eₕ, 10⁻⁹ cm⁻¹ |
| H₂O (H, O, H) CCSD/aug-cc-pVDZ opt+freq against the built-in H₂O | energy identical to 10⁻¹⁰ Eₕ; wavenumbers within 0.1 cm⁻¹, IR intensities within 10⁻⁴ km/mol |
| Point groups and symmetry numbers of 18 molecules (C1, Cs, C2v, C2, C3v, D3h, Td, D2h, D2d, D3d, C3h, D4h, D6h, Oh, D∞h, C∞v) | all correct |

### Pople and def2 basis sets (version 14.0)

ΨLab 14.0 adds 25 basis sets to the eight correlation-consistent ones; the library
(`basis_set_library()` in `basis.cpp`) has 33 sets in five families:

| Family | Sets | Function types | Elements (ΨLab) |
|---|---|---|---|
| Pople 6-31G | 6-31G, 6-31G\*, 6-31+G, 6-31+G\*, 6-31+G\*\*, 6-31++G\*, 6-31++G\*\* | Cartesian d (6D), spherical f (7F) | 6-31G, 6-31G\*: H–Kr; the diffuse sets H–Ar |
| Pople 6-311G | 6-311G, 6-311G\*, 6-311+G, 6-311+G\*, 6-311+G\*\*, 6-311++G\*, 6-311++G\*\*, 6-311++G(2d,2p), 6-311++G(3df,3pd) | spherical (5D, 7F) | 6-311G, 6-311G\*: H–Ca, Ga–Kr, I; 6-311+G…: H–Ca; 6-311++G…: H, Li–Ca; (3df,3pd): H, Li–Ar |
| Karlsruhe def2 | def2-SVP, def2-mTZVP, def2-mTZVPP, def2-TZVP, def2-TZVPP, def2-TZVPPD, def2-QZVP, def2-QZVPP, def2-QZVPPD | spherical (5D, 7F) | H–Rn (TZVPPD, QZVPPD: H–La, Hf–Rn; mTZVP, mTZVPP also Th–Lr), with ECPs from Rb on |
| Dunning | aug-cc-pVDZ … aug-cc-pV5Z | spherical | H–Ar, Sc–Kr (not K, Ca) |
| Dunning core-valence | aug-cc-pCVDZ … aug-cc-pCV5Z | spherical | Li–Ar (H, He: aug-cc-pVXZ) |

- **Names**: any case; Gaussian's notation for the Pople polarisation functions, `6-31G(d)` =
  6-31G\*, `6-31+G(d,p)` = 6-31+G\*\*, `6-311++G(d,p)` = 6-311++G\*\*; `def2SVP` = def2-SVP. The
  files use Psi4's names (`+` → p, `*` → s, `(`, `,`, `)` → `_`: `6-31pgss.gbs`,
  `6-311ppg_3df_3pd_.gbs`).
- **Cartesian d functions (6D)**: the six d functions x², xy, xz, y², yz, z² (each normalised; their
  span contains the five spherical d functions and the s-type combination r²) of every d shell of
  the 6-31G family, as Gaussian uses these sets by default and as the Basis Set Exchange marks them.
  `Shell::pure` is false for such a shell; its output transform (`cartesian_transform(l)`) is the
  diagonal of the factors √[(2l−1)!! / ((2i−1)!! (2j−1)!! (2k−1)!!)] that normalise x^i y^j z^k
  (the contraction coefficients carry the normalisation of x^l). Every integral (one-electron,
  electron repulsion, dipole), every derivative integral (first and second derivatives, so the
  analytic gradients and Hessians of all methods), the orbital values and the output go through the
  shell's transform, so nothing else changes. The f functions of 6-31G\* (Sc–Zn) remain spherical,
  as in Gaussian (6D 7F). The BSE marks the Na–Ar d functions of 6-311G\* Cartesian (they come from
  the McLean–Chandler sets); Gaussian, and ΨLab, use the whole 6-311G family with 5D.
- **Elements**: each file holds every element for which the BSE defines the set (the Pople sets up
  to Kr, 6-311G and 6-311G\* also I; the def2 sets with their ECPs, below). A molecule with an
  element outside a set is refused with the elements the set covers (e.g. "6-311++G\*\* is not
  available for He: ΨLab has it for H, Li-Ca"). def2-QZVP and def2-QZVPP are
  identical for H–Ne and Al–Ar (they differ for Na, Mg, K–Kr), so they give the same results for
  the built-in molecules.

Number of basis functions (6D for the 6-31G family):

| | H₂ | HeH⁺ | LiH | H₂O | NH₃ | H₂O₂ |
|---|---|---|---|---|---|---|
| 6-31G | 4 | 4 | 11 | 13 | 15 | 22 |
| 6-31G\* | 4 | 4 | 17 | 19 | 21 | 34 |
| 6-31+G\*\* | 10 | 10 | 24 | 29 | 34 | 48 |
| 6-31++G\*\* | 12 | 12 | 25 | 31 | 37 | 50 |
| 6-311G\* | 6 | 6 | 21 | 24 | 27 | 42 |
| 6-311++G\*\* | 14 | – | 29 | 36 | 43 | 58 |
| 6-311++G(3df,3pd) | 36 | – | 57 | 75 | 93 | 114 |
| def2-SVP | 10 | 10 | 14 | 24 | 29 | 38 |
| def2-TZVP | 12 | 12 | 20 | 43 | 49 | 74 |
| def2-TZVPPD | 34 | 34 | 39 | 74 | 88 | 114 |
| def2-QZVPPD | 66 | 66 | 68 | 132 | 162 | 198 |

(all 25 sets: [`basis/README.md`](basis/README.md)).

#### Effective core potentials (14.0)

Beyond Kr the def2 sets are valence basis sets for the Stuttgart–Cologne pseudopotentials (the
def2-ECP): 28 core electrons for Rb–Xe, 46 for Cs–La, 28 with the 4f shell in the valence for
Ce–Lu, 60 for Hf–Rn. The files hold them in Gaussian's ECP format (`I-ECP 3 28`, then the local
part U_L and U_l − U_L for l < L, each a sum of terms c·exp(−a r²)). For def2-mTZVP and
def2-mTZVPP, which the BSE gives without ECPs, ΨLab attaches the def2-ECP (Rb–Rn: their valence
bases are those of def2-TZVP, the same largest exponents) and the small-core ECP60MWB of the
actinides (Th–Lr: their (14s13p10d8f1g)/[10s9p5d4f1g] valence basis is the Cao–Dolg segmented set
of ECP60MWB, Turbomole's def-TZVP), as ORCA uses these sets; the files say so in their header.

- **Integrals** (`ecp.cpp`): ⟨a|U|b⟩ for Cartesian shells up to l = 9. The local part (a Gaussian in
  r about the ECP centre C) is a three-centre overlap (Obara–Saika). For a projected part, the
  angular projection of a shell function onto Y_lm on the sphere of radius r about C follows from
  the binomial expansion of (r − a)^c (a = A − C) and the expansion of exp(2α r a·r̂) in modified
  spherical Bessel functions; the angular integrals ∫ Y_λμ Y_lm x^p y^q z^s are exact (polynomials
  on the sphere), the radial integral of the two projections times exp(−ζ r²) r² is a 64-point
  Gauss–Legendre quadrature over ±9 widths about the peak of each Gaussian product (64 against 160
  points: 3.5e-14 relative). Scaled Bessel functions e^−x i_n(x): series, Miller's downward and
  upward recurrence by range.
- **Nuclear charges and derivatives**: the atom carries its ECP (`Atom::ecp`, set by `attach_ecps`
  when the basis is built); the electrons see Z less the core electrons (electron count, nuclear
  repulsion, dipoles, Mulliken charges). The ECP is part of the nuclear attraction of its atom, so
  every first and second derivative — with respect to the shell centres through the component
  shells, with respect to C by translational invariance — and hence every analytic gradient,
  Hessian, CPHF and response equation includes it.
- **Reports**: the molecule ("HI, 26 electrons (and 28 core electrons in ECPs)"), the basis block
  (`I def2-SVP (10s7p6d) -> [4s4p2d] 26 functions, ECP28`) and, in the appendix, the potential in
  full; JSON `molecule.coreElectrons`, `basis.atoms[].ecp`, `orbitals.atoms[].ecp`.

**Verification** (`macqc_basis_tests`; `reference/pyscf_ecp_reference.py` in CI; PySCF with the
same files, basis functions and ECPs):

| Check | Result |
|---|---|
| RHF energies: HI in def2-SVP, TZVP, QZVPP (h functions), TZVPPD, mTZVP; Xe/def2-QZVPP; CsF (ECP46); AuH (ECP60); PbO; ICl; LuF₃ (lanthanide ECP28, projectors up to g); ThO₂/def2-mTZVPP (ECP60MWB); HI/6-311G\* all-electron | ≤ 1.7·10⁻¹¹ Eₕ |
| CCSD(T)/def2-SVP of HI | 4·10⁻¹¹ Eₕ |
| Analytic RHF gradients (HI, CsF, PbO, LuF₃), CCSD gradient of HI | ≤ 3.3·10⁻¹⁰ Eₕ/bohr |
| RHF Hessians of HI and PbO | ≤ 3.4·10⁻¹⁰ Eₕ/bohr² |
| Analytic CCSD(T) Hessian of HI with ECP against 4-point differences of analytic gradients; dipole derivatives | 9.9·10⁻¹¹ Eₕ/bohr²; 3.8·10⁻¹⁰ e |
| AuH, PbO turned and shifted: energy unchanged | 4·10⁻¹³ Eₕ |
| e^−x i_n(x) against closed forms and the recurrence, every range | ≤ 3.3·10⁻¹⁴ |

**SCF guess (14.0)**: the SCF starts from the superposition of atomic densities (SAD) — for each
element (and ECP) an SCF of the neutral atom in its basis functions, filling the orbitals by
energy with the electrons of the last, partly filled group of degenerate orbitals spread evenly
over it, damped — instead of the core-Hamiltonian orbitals (`RHFOptions::guess = "core"`), which
left CsF, LuF₃ and ThO₂ unconverged after 300 iterations; with SAD they converge in 11–16.
Optimisations and finite differences still start from the orbitals of the previous geometry.

**Verification of the Pople and def2 sets up to Kr** (`macqc_basis_tests`; `reference/pyscf_basis_sets_reference.py` in CI; PySCF with
`cart=True` for the 6-31G family, the 6D 7F case of Zn from PySCF's Cartesian integrals transformed
to spherical f functions):

| Check | Result |
|---|---|
| RHF energies of H₂O in all 25 sets | ≤ 6.3·10⁻¹³ Eₕ |
| RHF energies of LiH, NH₃, CH₄, HCl (6D in 6-31G\*, 5D in 6-311G\*), KH (6-31G: Cartesian valence d), CaH₂, ZnH₂ (6-31G\*: 6D 7F; def2-TZVPP: g functions), HBr, HeH⁺ | ≤ 2.3·10⁻¹² Eₕ |
| CCSD(T) energies (H₂O 6-31G\* and def2-SVP, NH₃ 6-311++G\*\*) | ≤ 7·10⁻¹¹ Eₕ |
| Analytic RHF gradients (6-31+G\*\*, HCl 6-31G\*, def2-TZVP), analytic CCSD/6-31G\* gradient | ≤ 7.6·10⁻¹⁰, 3.1·10⁻⁹ Eₕ/bohr |
| RHF/6-31G\* Hessian, dipole moment, orbital energies, HOMO and LUMO values at points | 7.9·10⁻¹⁰ Eₕ/bohr², ≤ 3·10⁻⁹ |
| Analytic Hessians (RHF, CCSD/6-31+G\*, CCSD(T)/6-31G\*, Full CI LiH/6-31G\*) against 4-point differences of analytic gradients; dipole derivatives | ≤ 1.1·10⁻⁹ Eₕ/bohr²; ≤ 1.4·10⁻⁸ e |
| Cartesian d of O: ⟨xx\|xx⟩ = ⟨xy\|xy⟩ = 1, ⟨xx\|yy⟩ = 1/3, ⟨xx\|xy⟩ = 0; unit diagonal of S in all sets | ≤ 2·10⁻¹⁵ |
| Integral-direct against stored-integral RHF (6-31++G\*\*) | 10⁻¹³ Eₕ |

## App sources

The SwiftUI sources of ΨLab are in `app/MacQC/`, the app icon in `app/MacQC/Assets.xcassets`
(regenerated by `app/icon/make_icon.py`, master image `app/icon/AppIcon-1024.png`), and the
XcodeGen project spec in `xcode/project.yml`. See [`app/SETUP.md`](app/SETUP.md) for building.

## Continuous integration

- `.github/workflows/backend.yml` builds the C++ backend and runs its tests on `ubuntu-latest`
  (Python 3 + NumPy regenerate the expected reports and verify they are unchanged). A second job
  installs PySCF and recomputes the stored integral and RHF reference values; run manually with
  the input *full*, a third job checks every RHF case up to 341 basis functions. A macOS job
  (`macos-26`, Apple clang, Accelerate) builds the backend and `psilab`, runs all tests, installs
  `psilab` and runs an optimisation + frequencies + single-point job with it, with the basis-set
  files taken from the installed copy.
- `.github/workflows/app.yml` builds the universal ΨLab app on a `macos-26` runner (XcodeGen +
  Xcode), checks that the basis-set files are inside the bundle, runs a UI smoke test (launch, all
  tabs, potential-curve scans, menus, light and dark appearance; H₂O RHF, LiH Full CI, an H₂ curve,
  the refused Full CI of H₂O, a long-calculation warning, H₂O CCSD(T), the open method menu, LiH
  CCSDT and the refused CCSDTQ of H₂O/aug-cc-pVTZ; for 5.0 the open job menu, an H₂O RHF
  optimisation, LiH CCSD(T) frequencies and thermochemistry, H₂O optimisation + frequencies with
  a CCSD(T)/aug-cc-pVTZ single point, the Output tab while a job runs, Cancel and the refused
  CCSDTQ frequencies of H₂O; for 6.0 the molecule of the optimisation steps, still and while
  *Play steps* runs, and the animated normal modes of LiH and H₂O (two captures 0.3 s apart must
  differ; a few screenshots are printed into the log as small JPEGs); for 7.0 HeH⁺ Full CI/aug-cc-pVTZ
  optimisation + frequencies, the relaxed torsion scan of H₂O₂ (also played), the transition
  states of NH₃ and H₂O₂ (CalcAll), the scan settings of NH₃, the refused TS search of LiH; for 8.0
  the orbitals (isosurfaces of the H₂O HOMO and LUMO, a contour plot of NH₃, the natural orbitals of
  HeH⁺ Full CI, the Szabo–Ostlund orbitals of HeH⁺), selected with launch arguments
  `-psilabMolecule`, `-psilabBasis`, `-psilabMethod`, `-psilabJob`, `-psilabSPMethod`,
  `-psilabSPBasis`, `-psilabR` (Å), `-psilabAngle`, `-psilabROO` (Å), `-psilabDihedral`,
  `-psilabScan`, `-psilabScanStart`, `-psilabTSHessian`, `-psilabEigenTest`, `-psilabTab`,
  `-psilabOrbital`, `-psilabOrbitalSet`, `-psilabOrbitalView`, `-psilabIsovalue`,
  `-psilabAutoRun`, `-psilabPlaySteps`, (13.0) `-psilabOpenOutput` (opens an output file at launch,
  without a calculation; the app is sandboxed, so the path must lie in its container: the smoke test
  copies the example output of H₂O CCSD opt+freq there and shows its Summary, Optimisation,
  Frequencies and Orbitals pages)
  and the viewer settings `-viewer3D.*` (12.0); screenshots on the runner's 1024-pt-wide screen,
  where the window has its minimum width of 1000 pt) and publishes the **ΨLab-app** download.
  The first build of `main` with a new app version is published as a GitHub Release (tag
  `v<version>`, asset `PsiLab-<version>.app.zip`). A published release is not changed by later
  builds: a new release needs a new `MARKETING_VERSION`.
- `.github/workflows/releases.yml` lists the builds of versions 1.0 to 4.0, 6.0 and 7.0 (CI
  artifacts) and makes their releases hold them. For a tag at the version's commit it creates a
  missing release on the tag, or attaches the build and sets the notes of a release that lacks
  them. It never creates, moves or deletes a tag; a tag that is missing or at another commit is
  reported, and the run fails. It runs when the file changes on `main` or by hand.
- **Tags at older commits.** The workflow token (`GITHUB_TOKEN`) may create or move a tag at an
  older commit only while a *branch* points at a commit with the same workflow files
  (`.github/workflows/`); tags do not count. Otherwise GitHub refuses it with HTTP 403
  (*Resource not accessible by integration*; with `git push`: *refusing to allow a GitHub App to
  create or update workflow … without `workflows` permission*). GitHub does not document this;
  it is what the runs of November 2026 show: v1.0, v2.0 and v4.0 were created while branches
  pointed at commits with their workflow files, v6.0 (`9904c29`) was refused three times although
  the tag v5.0 has the same workflow files and was created at the first try once the branch
  `claude/release-6.0` pointed at `9904c29`, and v7.0 was moved to `c11b1d7`, the tip of
  `claude/fervent-thompson-voh6gk`. To set such a tag:
  1. a branch at the commit: `git push origin <commit>:refs/heads/<name>`, or on github.com open
     `…/tree/<commit>` and create it with the branch menu;
  2. the tag, either by a workflow run (`gh api -X POST repos/<owner>/<repo>/git/refs -f
     ref=refs/tags/v<version> -f sha=<commit>`; to move a tag, `-X PATCH
     …/git/refs/tags/v<version> -f sha=<commit> -F force=true`; the one-off workflow of commit
     `5033893`, run with *Run workflow* on its branch, did this for v6.0 and v7.0) or on
     github.com: *Releases* ▸ *Draft a new release*, tag `v<version>` (*Create new tag on
     publish*), that branch as target, no assets, without *Set as the latest release*;
  3. run `releases.yml`, which attaches the build and the notes. The branch can then be deleted.
- Until 8.0 every build of `main` replaced the release of its
  version: the builds of the 7.0 stages had replaced `v6.0`, and the orbital viewer, merged while
  numbered 7.0, had become `v7.0`. `v6.0` is again the build of 6.0, and `v7.0` is 6.0 with the
  scans, transition states, H₂O₂, NH₃ and HeH⁺ of stages 7.1–7.4 (commit `c11b1d7`), without
  the orbital viewer.

## Version 5.0 (next version)

ΨLab 5.0 is the next version of the application, after 4.0. It adds geometry optimisation,
harmonic frequencies, IR intensities and thermochemistry for H₂, LiH and H₂O, five job types and a
complete output file; HeH⁺ (Szabo–Ostlund) is unchanged (its Hartree–Fock text output is still
identical to `reference/hartree_fock_sto3g.py`).

**Job types** (app: the sidebar's *Job* menu; `psilab`: `job = …`; C API: `"job"`):
1 *energy*; 2 *opt* (geometry optimisation); 3 *freq* (frequencies and thermochemistry at the
input geometry, with a warning if it is not stationary); 4 *opt+freq* (frequencies at the
optimised geometry); 5 *opt+freq+sp* (then a single point with `sp_method` / `sp_basis`, final
U, H, G, A = single-point energy + thermal corrections). Jobs 4 and 5 stop with a message if the
optimisation does not converge. Every step prints the geometry it uses.

**What is analytic and what is numerical**

| Method | Gradient | Hessian | Dipole derivatives (IR) |
|---|---|---|---|
| HF (RHF) | analytic (derivative integrals, energy-weighted density) | analytic (second-derivative integrals, CPHF) | analytic (CPHF) |
| CCSD, CCSD(T) | analytic (Λ equations, relaxed densities, Z-vector) | central differences of analytic gradients along symmetry-adapted internal coordinates (H₂O 6, diatomics 2 gradients) | relaxed dipoles of those gradients |
| Full CI | analytic (relaxed densities, Z-vector) | as CCSD(T) | as CCSD(T) |
| CCSDT, CCSDTQ | numerical (central differences of energies along internal coordinates: H₂O 4, diatomics 2 energies) | second differences of energies (H₂O 9, diatomics 2) | finite fields μ = −dE/dF (H₂O 10, LiH 4 more energies) |

Correlated methods also differentiate numerically when near-linear dependencies were removed from
the basis or the analytic gradient would not fit into the memory limit, and RHF uses numerical
second derivatives of its analytic gradients then; the output file, the app and the JSON say
which was used and why. CCSDTQ optimisations and frequencies of H₂O are refused (hours per energy
and more than 5 GB); a single point with CCSDTQ in job 5 remains possible where the energy is.

**Verification** (each stage has its own test program, run by `ctest` and in CI; references from
PySCF 2.15 and geomeTRIC by the scripts in `reference/`):

| Stage | What | Result |
|---|---|---|
| 1 | first and second derivative integrals vs finite differences | 7e-10 relative (first); S 1e-11, h 2e-9, ERI 5e-9 (second) |
| 1 | RHF gradient vs finite differences / PySCF (aug-cc-pVDZ … aug-cc-pV5Z) | 2e-10 / 1.2e-10 Eₕ/bohr |
| 1 | RHF Hessian vs PySCF / finite differences; dipole derivatives | 6e-10 / 1e-10 Eₕ/bohr²; 5e-10 e |
| 2 | CCSD, CCSD(T), Full CI gradients vs finite differences | ≤ 1.7e-10 Eₕ/bohr |
| 2 | vs PySCF (`grad.ccsd`, `grad.ccsd_t`, CASCI) | ≤ 7.5e-10 (CC), 4.9e-9 (Full CI) |
| 3 | optimised geometries vs PySCF + geomeTRIC (8 cases) | ≤ 8e-6 Å, 6e-4°, 1e-10 Eₕ |
| 4 | RHF wavenumbers vs PySCF analytic Hessian | ≤ 2e-3 cm⁻¹ |
| 4 | thermochemistry vs the formulas evaluated by hand (CODATA 2018) / PySCF `thermo` | 1e-9 Eₕ, 5e-5 J/(mol K) / 3e-6 Eₕ |
| 4 | numerical vs analytic Hessians; relaxed dipoles vs −dE/dF | ≤ 0.05 cm⁻¹; 1e-9 e bohr |
| 5 | all job types through the C API vs PySCF; H₂: CCSD = CCSDT = CCSDTQ = Full CI; LiH: CCSDTQ = Full CI (*full* job) | 2e-5 Å, 0.01–0.05 cm⁻¹; 1.5e-9 Å (analytic), 7.4e-6 Å and 0.15 cm⁻¹ (numerical); 6e-12 Eₕ and 2.4e-3 cm⁻¹ |
| 5 | streamed output = final output file; cancellation; refusals; estimates | identical; < 1 ms during an optimisation step |

**Comparison with experiment.** Harmonic constants from optimisations and frequencies with the
largest feasible calculations (all electrons correlated):

| Quantity | ΨLab | Experiment |
|---|---|---|
| H₂ rₑ, Full CI/aug-cc-pV5Z | 0.741602 Å | 0.74144 Å |
| H₂ ωₑ, Full CI/aug-cc-pV5Z | 4403.42 cm⁻¹ | 4401.21 cm⁻¹ |
| LiH rₑ, Full CI/aug-cc-pCVTZ | 1.601196 Å | 1.5949 Å |
| LiH ωₑ, Full CI/aug-cc-pCVTZ | 1393.36 cm⁻¹ | 1405.65 cm⁻¹ |
| H₂O rₑ, θₑ, CCSD(T)/aug-cc-pVTZ | 0.958852 Å, 104.342° | 0.95782 Å, 104.48° |
| H₂O ω₁ (A₁), ω₂ (A₁), ω₃ (B₂), CCSD(T)/aug-cc-pVTZ | 3838.53, 1641.38, 3942.48 cm⁻¹ | 3832.2, 1648.5, 3942.5 cm⁻¹ |
| same, CCSD(T)/aug-cc-pCVTZ | 3813.59, 1645.52, 3920.97 cm⁻¹ (rₑ 0.960573 Å, θₑ 104.289°) | |
| S°(H₂O, 298.15 K, 1 bar), CCSD(T)/aug-cc-pVTZ | 188.69 J/(mol K) | 188.83 J/(mol K) |

Sources: H₂O harmonic frequencies from the NIST Computational Chemistry Comparison and Benchmark
Database (experimental data); H₂O equilibrium structure from Császár et al., J. Chem. Phys. 122,
214305 (2005) (rₑ = 0.95782 Å, θₑ = 104.48(5)°; the default geometry of ΨLab, 0.9572 Å and
104.52°, is the older, widely quoted value); S°(H₂O) from the NIST-JANAF tables (188.834) and the
CODATA key values for thermodynamics (188.835 ± 0.010); LiH and H₂ spectroscopic constants from
Huber & Herzberg, *Constants of Diatomic Molecules* (1979). The H₂ values (Full CI in a nearly
complete basis) differ by 2 cm⁻¹ and 2e-4 Å — what the Born–Oppenheimer approximation leaves
out. For LiH and H₂O the differences are those of the basis sets (aug-cc-pCVTZ and aug-cc-pVTZ are
not complete; for H₂O also the missing core-correlating functions in aug-cc-pVTZ) and, for
CCSD(T), of the method. The rigid-rotor, harmonic-oscillator entropy of H₂O is 0.14 J/(mol K)
below the tabulated value, the size of the anharmonic and rotation–vibration contributions it
leaves out.

**Terminal version on macOS**: `psilab` builds with CMake on macOS too, with Apple's Accelerate for
the matrix products; CI builds and tests it on Linux and on macOS (all tests, plus an installed
`psilab` running job 5).

**Not in this version**: analytic gradients of CCSDT and CCSDTQ (they need general-order Λ
equations; the gradients are numerical), analytic Hessians of the correlated methods (numerical
from analytic gradients), analytic correlated gradients with removed near-linear dependencies
(numerical then), transition-state searches, anharmonic frequencies, open-shell molecules.

**Usage**: choose the job in the sidebar's **Job** menu, the method and basis set (and for job 5
the single-point method and basis set), set the starting geometry, temperature, pressure and
convergence, and press **Run**; the *Optimisation*, *Frequencies* and *Output* tabs fill as the
job runs. In the terminal: `job = opt+freq` (see [Terminal version](#terminal-version-linux-and-macos));
`tail -f` follows the output file.
