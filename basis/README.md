# Basis sets

The Gaussian basis sets of ΨLab, in Gaussian94 format (`.gbs`), exactly as written by the
[Basis Set Exchange](https://www.basissetexchange.org) (BSE, data version 0.12), downloaded on
2026-11-12 with

```
https://www.basissetexchange.org/api/basis/<name>/format/gaussian94/?elements=<Z,...>
```

ΨLab reads these files at run time (`backend/src/basis.cpp`); they are not edited by hand. The
functions are used as spherical harmonics (5d, 7f, 9g, 11h), as the correlation-consistent sets are
designed — except (14.0) the d shells of the 6-31G family, which are six Cartesian functions (6D),
as Gaussian uses them (see [ΨLab 14.0](#ψlab-140-the-pople-6-31g-and-6-311g-families-and-the-def2-sets)). Contracted functions are normalised; the contraction coefficients refer to normalised
primitives (Gaussian94 convention).

| File | Basis set | H | He | Li | N | O | highest l |
|---|---|---|---|---|---|---|---|
| `aug-cc-pvdz.gbs`  | aug-cc-pVDZ  | [3s2p] | [3s2p] | [4s3p2d] | [4s3p2d] | [4s3p2d] | d |
| `aug-cc-pvtz.gbs`  | aug-cc-pVTZ  | [4s3p2d] | [4s3p2d] | [5s4p3d2f] | [5s4p3d2f] | [5s4p3d2f] | f |
| `aug-cc-pvqz.gbs`  | aug-cc-pVQZ  | [5s4p3d2f] | [5s4p3d2f] | [6s5p4d3f2g] | [6s5p4d3f2g] | [6s5p4d3f2g] | g |
| `aug-cc-pv5z.gbs`  | aug-cc-pV5Z  | [6s5p4d3f2g] | [6s5p4d3f2g] | [7s6p5d4f3g2h] | [7s6p5d4f3g2h] | [7s6p5d4f3g2h] | h |
| `aug-cc-pcvdz.gbs` | aug-cc-pCVDZ | – | – | [5s4p2d] | [5s4p2d] | [5s4p2d] | d |
| `aug-cc-pcvtz.gbs` | aug-cc-pCVTZ | – | – | [7s6p4d2f] | [7s6p4d2f] | [7s6p4d2f] | f |
| `aug-cc-pcvqz.gbs` | aug-cc-pCVQZ | – | – | [9s8p6d4f2g] | [9s8p6d4f2g] | [9s8p6d4f2g] | g |
| `aug-cc-pcv5z.gbs` | aug-cc-pCV5Z | – | – | – | [11s10p8d6f4g2h] | [11s10p8d6f4g2h] | h |

Nitrogen (for NH3) and helium (for HeH+ with correlation-consistent basis sets) were added in
ΨLab 7.0 with the offline copy of the same data, `pip install basis_set_exchange` (BSE data
version 0.12):

```
basis_set_exchange.get_basis(<name>, elements=[1, 2, 3, 7, 8], fmt="gaussian94", header=True)
```

With `elements=[1, 3, 8]` this reproduces the eight earlier files byte for byte, so the files
differ from them only by the inserted He and N blocks. The core-valence sets are not defined for
He (its 1s shell is the valence shell); like H, He then takes the aug-cc-pVXZ functions.

**ΨLab 9.0: more elements.** For custom molecules the eight files were extended with the same
offline copy of the BSE data (version 0.12):

```
basis_set_exchange.get_basis(<name>, elements=<Z, ...>, fmt="gaussian94", header=True)
```

with every element the BSE defines among H-Ar and Sc-Kr: aug-cc-pVDZ, TZ and QZ for H-Ar and
Sc-Kr, aug-cc-pV5Z for H-Ar and Ga-Kr (the 5Z sets of Sc-Zn contain i functions, l = 6, beyond
ΨLab's integrals), aug-cc-pCVDZ, TZ and QZ for Li-Ar, aug-cc-pCV5Z for B-Ar. The BSE has no
aug-cc-pVXZ sets for K and Ca; elements beyond Kr need effective core potentials (aug-cc-pVXZ-PP),
which ΨLab does not implement. The blocks of H, He, Li, N and O are unchanged byte for byte (checked
when the files were written). Molecules with other elements can be read and shown, and are refused
with this explanation when a calculation is asked for.

**Element coverage.**

- The core-valence sets aug-cc-pCVXZ are not defined for hydrogen, which has no core electrons.
  As is customary, ΨLab uses the matching aug-cc-pVXZ functions for H (e.g. aug-cc-pVDZ for H in
  an aug-cc-pCVDZ calculation) and says so in the output.
- aug-cc-pCV5Z is defined in the BSE for B–Ar only, not for Li. ΨLab therefore refuses LiH in
  aug-cc-pCV5Z with an explanation (use aug-cc-pCVQZ or aug-cc-pV5Z).

Number of basis functions per molecule (spherical):

| | aug-cc-pVDZ | TZ | QZ | 5Z | aug-cc-pCVDZ | CVTZ | CVQZ | CV5Z |
|---|---|---|---|---|---|---|---|---|
| H2  | 18 | 46 | 92 | 160 | 18 | 46 | 92 | 160 |
| HeH+ | 18 | 46 | 92 | 160 | 18 | 46 | 92 | 160 |
| LiH | 32 | 69 | 126 | 207 | 36 | 82 | 155 | – |
| H2O | 41 | 92 | 172 | 287 | 45 | 105 | 201 | 341 |
| NH3 | 50 | 115 | 218 | 367 | 54 | 128 | 247 | 421 |
| H2O2 | 64 | 138 | 252 | 414 | 72 | 164 | 310 | 522 |

## ΨLab 14.0: the Pople 6-31G and 6-311G families and the def2 sets

25 more files, written with the same offline copy of the BSE data (version 0.12) for every element
for which the BSE defines the set, with the effective core potentials (ECPs) of the def2 sets beyond
Kr:

```
basis_set_exchange.get_basis(<name>, fmt="gaussian94", header=True)
```

(ΨLab 14.0 first had these sets for H–Kr only; the blocks of H–Kr are unchanged byte for byte.)
def2-mTZVP and def2-mTZVPP are given by the BSE without ECPs for their elements beyond Kr; their
files add the def2-ECP of def2-SVP (Rb–Rn) and the small-core ECP60MWB ("Stuttgart RSC 1997 ECP"
of the BSE, Th–Lr) to the BSE's basis through `basis_set_exchange.write_formatted_basis_str`, and
their headers say so. These are the ECPs of their valence bases: for Rb–Rn they are derived from
def2-TZVP (the same largest exponents), for Th–Lr they are the (14s13p10d8f1g)/[10s9p5d4f1g]
reduction of the Cao–Dolg segmented basis of ECP60MWB (Turbomole's def-TZVP); ORCA uses these sets
with the same ECPs.

The file names follow Psi4's library (`+` → `p`, `*` → `s`, `(`, `,`, `)` → `_`):

| File | Basis set | Elements | Function types (ΨLab, as Gaussian) |
|---|---|---|---|
| `6-31g.gbs`, `6-31gs.gbs` | 6-31G, 6-31G\* | H–Kr | Cartesian d (6D), spherical f (7F) |
| `6-31pg.gbs`, `6-31pgs.gbs`, `6-31pgss.gbs` | 6-31+G, 6-31+G\*, 6-31+G\*\* | H–Ar | 6D |
| `6-31ppgs.gbs`, `6-31ppgss.gbs` | 6-31++G\*, 6-31++G\*\* | H–Ar | 6D |
| `6-311g.gbs`, `6-311gs.gbs` | 6-311G, 6-311G\* | H–Ca, Ga–Kr, I (all-electron) | spherical (5D, 7F) |
| `6-311pg.gbs`, `6-311pgs.gbs`, `6-311pgss.gbs` | 6-311+G, 6-311+G\*, 6-311+G\*\* | H–Ca | 5D |
| `6-311ppgs.gbs`, `6-311ppgss.gbs`, `6-311ppg_2d_2p_.gbs` | 6-311++G\*, 6-311++G\*\*, 6-311++G(2d,2p) | H, Li–Ca | 5D |
| `6-311ppg_3df_3pd_.gbs` | 6-311++G(3df,3pd) | H, Li–Ar | 5D 7F |
| `def2-svp.gbs`, `def2-tzvp.gbs`, `def2-tzvpp.gbs`, `def2-qzvp.gbs`, `def2-qzvpp.gbs` | def2-SVP, def2-TZVP, def2-TZVPP, def2-QZVP, def2-QZVPP | H–Rn, ECPs Rb–Rn | spherical (5D, 7F; up to g, def2-QZVPP up to h) |
| `def2-tzvppd.gbs`, `def2-qzvppd.gbs` | def2-TZVPPD, def2-QZVPPD | H–La, Hf–Rn, ECPs Rb–La, Hf–Rn | spherical |
| `def2-mtzvp.gbs`, `def2-mtzvpp.gbs` | def2-mTZVP, def2-mTZVPP | H–Rn, Th–Lr, ECPs Rb–Rn (def2-ECP), Th–Lr (ECP60MWB) | spherical |

**Function types.** The Gaussian94 format does not say whether a d shell has five spherical or six
Cartesian functions; ΨLab follows Gaussian's defaults, which are also the BSE's marks for the
6-31G family: the d shells of the 6-31G family are Cartesian (6D), every f shell and every shell of
the other sets is spherical. (The BSE marks the Na–Ar d functions of 6-311G\* Cartesian, from the
McLean–Chandler sets; Gaussian uses the 6-311G family with 5D throughout, and so does ΨLab.) For
H₂O, RHF/6-31G\* is −76.0105299763 Eₕ with 6D and −76.0091323803 Eₕ with 5D.

**The ECPs** (Gaussian's format: `I-ECP 3 28`, the local part U_L and U_l − U_L for l < L, terms
`2 exponent coefficient`, i.e. c·exp(−a r²)): 28 core electrons for Rb–Xe, 46 for Cs–La, 28 for the
lanthanides Ce–Lu (4f in the valence; local part l = 5, projectors up to g), 60 for Hf–Rn and
(ECP60MWB) Th–Lr. They are identical in all def2 files. ΨLab supports the Gaussian terms (n = 2),
which all of them are.

**def2-QZVP and def2-QZVPP** are identical for H–Ne and Al–Ar (the files differ for Na, Mg and
K–Kr). The Pople data are those of Gaussian 09 (6-31G family) and of the original BSE (6-311G
family); the def2 data those of Turbomole 7.3, def2-mTZVP(P) from S. Ehlert.

Contractions and numbers of basis functions of the built-in molecules (6D for the 6-31G family):

| Basis set | H | He | Li | N | O | H₂ | HeH⁺ | LiH | H₂O | NH₃ | H₂O₂ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 6-31G | [2s] | [2s] | [3s2p] | [3s2p] | [3s2p] | 4 | 4 | 11 | 13 | 15 | 22 |
| 6-31G\* | [2s] | [2s] | [3s2p1d] | [3s2p1d] | [3s2p1d] | 4 | 4 | 17 | 19 | 21 | 34 |
| 6-31+G | [2s] | [2s] | [4s3p] | [4s3p] | [4s3p] | 4 | 4 | 15 | 17 | 19 | 30 |
| 6-31+G\* | [2s] | [2s] | [4s3p1d] | [4s3p1d] | [4s3p1d] | 4 | 4 | 21 | 23 | 25 | 42 |
| 6-31+G\*\* | [2s1p] | [2s1p] | [4s3p1d] | [4s3p1d] | [4s3p1d] | 10 | 10 | 24 | 29 | 34 | 48 |
| 6-31++G\* | [3s] | [3s] | [4s3p1d] | [4s3p1d] | [4s3p1d] | 6 | 6 | 22 | 25 | 28 | 44 |
| 6-31++G\*\* | [3s1p] | [3s1p] | [4s3p1d] | [4s3p1d] | [4s3p1d] | 12 | 12 | 25 | 31 | 37 | 50 |
| 6-311G | [3s] | [3s] | [4s3p] | [4s3p] | [4s3p] | 6 | 6 | 16 | 19 | 22 | 32 |
| 6-311G\* | [3s] | [3s] | [4s3p1d] | [4s3p1d] | [4s3p1d] | 6 | 6 | 21 | 24 | 27 | 42 |
| 6-311+G | [3s] | [3s] | [5s4p] | [5s4p] | [5s4p] | 6 | 6 | 20 | 23 | 26 | 40 |
| 6-311+G\* | [3s] | [3s] | [5s4p1d] | [5s4p1d] | [5s4p1d] | 6 | 6 | 25 | 28 | 31 | 50 |
| 6-311+G\*\* | [3s1p] | [3s1p] | [5s4p1d] | [5s4p1d] | [5s4p1d] | 12 | 12 | 28 | 34 | 40 | 56 |
| 6-311++G\* | [4s] | – | [5s4p1d] | [5s4p1d] | [5s4p1d] | 8 | – | 26 | 30 | 34 | 52 |
| 6-311++G\*\* | [4s1p] | – | [5s4p1d] | [5s4p1d] | [5s4p1d] | 14 | – | 29 | 36 | 43 | 58 |
| 6-311++G(2d,2p) | [4s2p] | – | [5s4p2d] | [5s4p2d] | [5s4p2d] | 20 | – | 37 | 47 | 57 | 74 |
| 6-311++G(3df,3pd) | [4s3p1d] | – | [5s4p3d1f] | [5s4p3d1f] | [5s4p3d1f] | 36 | – | 57 | 75 | 93 | 114 |
| def2-SVP | [2s1p] | [2s1p] | [3s2p] | [3s2p1d] | [3s2p1d] | 10 | 10 | 14 | 24 | 29 | 38 |
| def2-mTZVP | [3s] | [3s] | [5s1p] | [5s3p1d] | [5s3p2d] | 6 | 6 | 11 | 30 | 28 | 54 |
| def2-mTZVPP | [2s1p] | [2s1p] | [5s2p] | [5s3p2d] | [5s3p2d] | 10 | 10 | 16 | 34 | 39 | 58 |
| def2-TZVP | [3s1p] | [3s1p] | [5s3p] | [5s3p2d1f] | [5s3p2d1f] | 12 | 12 | 20 | 43 | 49 | 74 |
| def2-TZVPP | [3s2p1d] | [3s2p1d] | [5s3p1d] | [5s3p2d1f] | [5s3p2d1f] | 28 | 28 | 33 | 59 | 73 | 90 |
| def2-TZVPPD | [3s3p1d] | [3s3p1d] | [5s4p1d] | [6s3p3d1f] | [6s4p3d1f] | 34 | 34 | 39 | 74 | 88 | 114 |
| def2-QZVP | [4s3p2d1f] | [4s3p2d1f] | [6s4p2d1f] | [7s4p3d2f1g] | [7s4p3d2f1g] | 60 | 60 | 65 | 117 | 147 | 174 |
| def2-QZVPP | [4s3p2d1f] | [4s3p2d1f] | [6s4p2d1f] | [7s4p3d2f1g] | [7s4p3d2f1g] | 60 | 60 | 65 | 117 | 147 | 174 |
| def2-QZVPPD | [4s4p2d1f] | [4s4p2d1f] | [6s4p2d1f] | [8s4p4d2f1g] | [8s5p4d2f1g] | 66 | 66 | 68 | 132 | 162 | 198 |

## References

Please cite the BSE and the original publications when you use these basis sets:

- B. P. Pritchard, D. Altarawy, B. Didier, T. D. Gibson, T. L. Windus, *A New Basis Set
  Exchange: An Open, Up-to-date Resource for the Molecular Sciences Community*,
  J. Chem. Inf. Model. **59**, 4814 (2019), doi:10.1021/acs.jcim.9b00725.
- D. Feller, J. Comput. Chem. **17**, 1571 (1996), doi:10.1002/(SICI)1096-987X(199610)17:13<1571::AID-JCC9>3.0.CO;2-P;
  K. L. Schuchardt et al., J. Chem. Inf. Model. **47**, 1045 (2007), doi:10.1021/ci600510j.
- H, O (cc-pVXZ): T. H. Dunning, Jr., J. Chem. Phys. **90**, 1007 (1989), doi:10.1063/1.456153.
- H, O (diffuse "aug" functions): R. A. Kendall, T. H. Dunning, Jr., R. J. Harrison,
  J. Chem. Phys. **96**, 6796 (1992), doi:10.1063/1.462569.
- Li (aug-cc-pVXZ): B. P. Prascher, D. E. Woon, K. A. Peterson, T. H. Dunning, Jr.,
  A. K. Wilson, Theor. Chem. Acc. **128**, 69 (2011), doi:10.1007/s00214-010-0764-0.
- O (core-valence functions, cc-pCVXZ): D. E. Woon, T. H. Dunning, Jr., J. Chem. Phys. **103**,
  4572 (1995), doi:10.1063/1.470645.
- Li (aug-cc-pCVXZ, X = D, T, Q): from the original BSE (Feller 1996; Schuchardt et al. 2007),
  diffuse functions from Prascher et al. (2011).
- 9.0, further elements: the references listed in the header and on the BSE page of each basis
  set (e.g. Al-Ar: D. E. Woon, T. H. Dunning, Jr., J. Chem. Phys. 98, 1358 (1993); Sc-Zn:
  N. B. Balabanov, K. A. Peterson, J. Chem. Phys. 123, 064107 (2005); Ga-Kr: A. K. Wilson,
  D. E. Woon, K. A. Peterson, T. H. Dunning, Jr., J. Chem. Phys. 110, 7667 (1999)).
- 14.0, 6-31G family: R. Ditchfield, W. J. Hehre, J. A. Pople, J. Chem. Phys. **54**, 724 (1971),
  doi:10.1063/1.1674902 (H); W. J. Hehre, R. Ditchfield, J. A. Pople, J. Chem. Phys. **56**, 2257
  (1972), doi:10.1063/1.1677527 (C–F); J. D. Dill, J. A. Pople, J. Chem. Phys. **62**, 2921 (1975),
  doi:10.1063/1.430801 (Li, B); J. S. Binkley, J. A. Pople, J. Chem. Phys. **66**, 879 (1977),
  doi:10.1063/1.433929 (Be); P. C. Hariharan, J. A. Pople, Theor. Chim. Acta **28**, 213 (1973),
  doi:10.1007/BF00533485 (d polarisation); M. M. Francl et al., J. Chem. Phys. **77**, 3654 (1982),
  doi:10.1063/1.444267, and M. S. Gordon et al., J. Am. Chem. Soc. **104**, 2797 (1982),
  doi:10.1021/ja00374a017 (Na–Ar); V. A. Rassolov et al., J. Chem. Phys. **109**, 1223 (1998),
  doi:10.1063/1.476673 (Sc–Zn), and J. Comput. Chem. **22**, 976 (2001), doi:10.1002/jcc.1058
  (K, Ca, Ga–Kr); diffuse functions: T. Clark et al., J. Comput. Chem. **4**, 294 (1983),
  doi:10.1002/jcc.540040303; G. W. Spitznagel et al., J. Comput. Chem. **8**, 1109 (1987),
  doi:10.1002/jcc.540080807; P. M. W. Gill et al., Chem. Phys. Lett. **197**, 499 (1992),
  doi:10.1016/0009-2614(92)85807-M.
- 14.0, 6-311G family: R. Krishnan, J. S. Binkley, R. Seeger, J. A. Pople, J. Chem. Phys. **72**,
  650 (1980), doi:10.1063/1.438955 (H–Ne); A. D. McLean, G. S. Chandler, J. Chem. Phys. **72**, 5639
  (1980), doi:10.1063/1.438980 (Na–Ar); J.-P. Blaudeau et al., J. Chem. Phys. **107**, 5016 (1997),
  doi:10.1063/1.474865 (K, Ca); L. A. Curtiss et al., J. Chem. Phys. **103**, 6104 (1995),
  doi:10.1063/1.470438 (Ga–Kr); multiple polarisation functions (2d,2p), (3df,3pd): M. J. Frisch,
  J. A. Pople, J. S. Binkley, J. Chem. Phys. **80**, 3265 (1984), doi:10.1063/1.447079; L. A. Curtiss
  et al., J. Chem. Phys. **109**, 7764 (1998), doi:10.1063/1.477422 (f of Ar).
- 14.0, def2: F. Weigend, R. Ahlrichs, Phys. Chem. Chem. Phys. **7**, 3297 (2005),
  doi:10.1039/b508541a; the diffuse def2-TZVPPD, def2-QZVPPD: D. Rappoport, F. Furche, J. Chem.
  Phys. **133**, 134105 (2010), doi:10.1063/1.3484283; def2-mTZVP: J. G. Brandenburg, C. Bannwarth,
  A. Hansen, S. Grimme, J. Chem. Phys. **148**, 064104 (2018), doi:10.1063/1.5012601 (B97-3c);
  def2-mTZVPP: S. Grimme, A. Hansen, S. Ehlert, J.-M. Mewes, J. Chem. Phys. **154**, 064103 (2021),
  doi:10.1063/5.0040021 (r²SCAN-3c).
- 14.0, all-electron iodine in 6-311G: M. N. Glukhovtsev, A. Pross, M. P. McGrath, L. Radom,
  J. Chem. Phys. **103**, 1878 (1995), doi:10.1063/1.469712.
- 14.0, the def2 ECPs (Stuttgart–Cologne pseudopotentials): Rb, Cs: G. Leininger et al., Chem. Phys.
  Lett. **255**, 274 (1996), doi:10.1016/0009-2614(96)00382-X; Sr, Ba: M. Kaupp et al., J. Chem.
  Phys. **94**, 1360 (1991), doi:10.1063/1.459993; Y–Cd, Hf–Hg: D. Andrae et al., Theor. Chim. Acta
  **77**, 123 (1990), doi:10.1007/BF01114537; In–Sb, Pb, Bi: B. Metz, H. Stoll, M. Dolg, J. Chem.
  Phys. **113**, 2563 (2000), doi:10.1063/1.1305880; Tl: B. Metz et al., Theor. Chem. Acc. **104**, 22
  (2000), doi:10.1007/s002149900101; Te–Xe, Po–Rn: K. A. Peterson et al., J. Chem. Phys. **119**,
  11113 (2003), doi:10.1063/1.1622924; La: M. Dolg et al., Theor. Chim. Acta **75**, 173 (1989),
  doi:10.1007/BF00528565, and **85**, 441 (1993), doi:10.1007/BF01112983; Ce–Yb: M. Dolg, H. Stoll,
  H. Preuss, J. Chem. Phys. **90**, 1730 (1989), doi:10.1063/1.456066; Lu: X. Cao, M. Dolg, J. Chem.
  Phys. **115**, 7348 (2001), doi:10.1063/1.1406535.
- 14.0, the actinide ECP60MWB (def2-mTZVP, def2-mTZVPP for Th–Lr): W. Küchle, M. Dolg, H. Stoll,
  H. Preuss, J. Chem. Phys. **100**, 7535 (1994), doi:10.1063/1.466847; the valence basis: X. Cao,
  M. Dolg, H. Stoll, J. Chem. Phys. **118**, 487 (2003), doi:10.1063/1.1521431; X. Cao, M. Dolg,
  J. Mol. Struct. THEOCHEM **673**, 203 (2004), doi:10.1016/j.theochem.2003.12.015.
