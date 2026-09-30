"""
degree_of_reduction.py
------------------------
Computes the REDUCTION DEGREE (gamma) of each metabolite, a MOLECULE-level
(not atom-level) quantity classical in microbial bioenergetics and
biochemical engineering (Roels, 1980, Biotechnol. Bioeng. 22:2457-2514;
Roels, 1983, "Energetics and Kinetics in Biotechnology", Elsevier). It
quantifies the number of electrons available per mole of compound,
relative to fully oxidised reference states: CO2 (carbon), H2O
(hydrogen/oxygen), NH3 (nitrogen), PO4(3-) (phosphorus), SO4(2-) (sulfur).

    gamma = 4*nC + nH - 2*nO - 3*nN + 5*nP + 6*nS - Q

where nC, nH, nO, nN, nP, nS are the numbers of atoms of each element in
the molecular formula, and Q is the net formal charge of the species.

DERIVATION / SELF-CONSISTENCY CHECK (why each coefficient has the sign and
magnitude it does): for every element, the coefficient equals the formal
oxidation state that element has in ITS OWN fully oxidised reference
compound, signed accordingly -- C is +4 in CO2, H is +1 in H2O, O is -2 in
H2O, N is -3 in NH3, S is +6 in SO4(2-), and P is +5 in PO4(3-).
gamma is therefore, equivalently, the total formal-oxidation-state sum the
molecule WOULD have if every atom were replaced by its reference form,
adjusted by -Q for the molecule's own net charge (needed so that
(de)protonation, an acid-base process, never changes gamma -- see below).

*** WHY THE PHOSPHORUS TERM MATTERS ***
The network contains many phosphorus-bearing metabolites (ATP, ADP,
NAD(P)(H), phospholipids, sugar phosphates...), so the +5*nP term is
required for gamma to be a correct ABSOLUTE, per-molecule quantity (e.g.
gamma of ATP considered on its own). It has, by contrast, no effect on any
ELECTRON-TRANSFER COUNT reported in this study: every donor/acceptor pair
used to compute a Delta-gamma is required, by construction, to have the
SAME number of phosphorus atoms on both sides (an explicit matching
criterion for ab initio pairs, matching.py; and simply
a chemical fact for every curated cofactor pair -- NAD+/NADH, NADP+/NADPH
and FAD/FADH2 do not gain or lose a phosphate group when they are
reduced). Because Delta-gamma = gamma_product - gamma_reactant and the
phosphorus contribution is IDENTICAL on both sides of every matched pair,
the +5*nP term cancels exactly and does not affect any electron count in
this study; it is nonetheless included so that gamma is correct as an
absolute quantity, not merely self-consistent by cancellation. See
Supplementary Validation and tests/test_validation.py for the reference
compounds used to validate this formula, and scripts/validate_reduction_degree.py
to reproduce that validation.

RATIONALE (molecule-level rather than atom-level):
A per-atom oxidation-number approach requires a unique structural
representation, including bond orders, ring/open-chain form and tautomeric
state. The reduction degree used here instead depends only on MOLECULAR
FORMULA and CHARGE, stable quantities already present in genome-scale
metabolic models (cobra `met.elements`, `met.charge`) and insensitive to
the particular structural representation chosen by a source database.
This provides the electron-content information required for robust
network-scale analysis without atom-resolved structures.

CHARGE TERM (-Q) -- WHY IT IS REQUIRED:
Roels' original formula (4C+H-2O-3N) was derived for neutral species
(fermentation substrates/biomass). For charged species -- essentially
all intracellular metabolites at physiological pH (deprotonated organic
acids, NAD+/NADH, phosphorylated intermediates...) -- a charge-adjustment
term is required, otherwise (de)protonation (an acid-base process, not a
redox process) is conflated with an actual change in electron content.
Verified here on two reference cases:
  - acetate (CH3COO-, charge -1) vs. acetic acid (CH3COOH, neutral): with
    NO charge term, gamma differs (7 vs. 8) even though (de)protonation is
    NOT a redox event -- gamma MUST be identical. WITH the charge term:
    gamma(acetate) = 7-(-1) = 8 = gamma(acetic acid). Correct.
  - NAD+ (C21H26N7O14P2, charge -1) -> NADH (C21H27N7O14P2, charge -2):
    WITHOUT the charge term, gamma(NADH)-gamma(NAD+) = 72-71 = +1 electron
    only (inconsistent with the universally established 2-electron
    transfer of this reaction, "NAD+ + 2e- + H+ -> NADH"). WITH the
    charge term: delta_gamma = 74-72 = +2 electrons. Correct.

SULFUR TERM (+6*nS):
Reference state SO4(2-) (S at oxidation state +6). This is the most
commonly cited standard coefficient (e.g. Rittmann & McCarty environmental
biotechnology; Heijnen's black-box biomass equations), though less
universally agreed in the literature than the C/H/O/N coefficients. Should
be re-examined if the analysis specifically targets sulfur metabolism.

PHOSPHORUS TERM (+5*nP):
Reference state PO4(3-) (P at oxidation state +5). P is
essentially always already at +5 in biological compounds, which means
this term does not, in practice, change a Delta-gamma between two
matched, phosphorus-count-conserving metabolites -- but it is required
for gamma itself (as an absolute, per-molecule quantity) to be correct.
"""

from typing import Dict, Optional

# Elements considered "safe" for the gamma calculation: either included in
# the formula (C, H, O, N, P, S) or simple ions that never change oxidation
# state in metabolism (Ca/K/Mg/Na/Cl remain simple fixed-charge ions). Any
# OTHER element (transition metals Fe/Cu/Mo/Co/Ni/Mn/Zn/Ag/Cd/Hg/W,
# redox-active metalloids such as Se or As, or "R"/"X" generic/undefined
# groups) makes gamma unreliable: such compounds are explicitly flagged
# as outside the supported chemical domain rather than assigned a potentially
# incorrect reduction degree.
SAFE_ELEMENTS = {"C", "H", "O", "N", "S", "P", "Ca", "Cl", "K", "Mg", "Na"}


def has_untracked_elements(elements: Dict[str, int]) -> bool:
    """True if the formula contains an element outside SAFE_ELEMENTS
    (transition metal, redox-active metalloid, generic formula group...)."""
    return any(el not in SAFE_ELEMENTS for el in (elements or {}))


def degree_of_reduction(elements: Dict[str, int], charge: int = 0) -> Optional[float]:
    """
    Computes gamma from an element-count dict (e.g. {'C': 6, 'H': 12, 'O': 6})
    and the net formal charge. Returns None if the formula is empty, or if
    it contains an element not covered by SAFE_ELEMENTS (returning nothing
    is preferable to returning a silently incorrect value).
    """
    if not elements:
        return None
    if has_untracked_elements(elements):
        return None
    c = elements.get("C", 0)
    h = elements.get("H", 0)
    o = elements.get("O", 0)
    n = elements.get("N", 0)
    p = elements.get("P", 0)
    s = elements.get("S", 0)
    charge = charge or 0
    return 4 * c + h - 2 * o - 3 * n + 5 * p + 6 * s - charge


def build_gamma_table(model) -> Dict[str, Optional[float]]:
    """Computes gamma for every metabolite in a cobrapy model."""
    table = {}
    for met in model.metabolites:
        try:
            table[met.id] = degree_of_reduction(met.elements, met.charge)
        except Exception:
            table[met.id] = None
    return table


if __name__ == "__main__":
    tests = {
        "glucose C6H12O6 (neutral)": ({"C": 6, "H": 12, "O": 6}, 0, 24),
        "ethanol C2H6O (neutral)": ({"C": 2, "H": 6, "O": 1}, 0, 12),
        "acetate CH3COO- (charge -1)": ({"C": 2, "H": 3, "O": 2}, -1, 8),
        "acetic acid CH3COOH (neutral)": ({"C": 2, "H": 4, "O": 2}, 0, 8),
        "phosphate PO4 3- (charge -3)": ({"O": 4, "P": 1}, -3, 0),
        "NAD+ (charge -1, 2 P)": ({"C": 21, "H": 26, "N": 7, "O": 14, "P": 2}, -1, 72),
        "NADH (charge -2, 2 P)": ({"C": 21, "H": 27, "N": 7, "O": 14, "P": 2}, -2, 74),
    }
    for name, (els, q, expected) in tests.items():
        g = degree_of_reduction(els, q)
        status = "OK" if g == expected else f"MISMATCH (expected {expected})"
        print(f"{name:35s} gamma={g:>4}  [{status}]")

    print()
    print("Delta gamma (NADH - NAD+) =",
          degree_of_reduction({"C": 21, "H": 27, "N": 7, "O": 14, "P": 2}, -2)
          - degree_of_reduction({"C": 21, "H": 26, "N": 7, "O": 14, "P": 2}, -1),
          "(expected: +2, a two-electron transfer; unaffected by the P term",
          "since NAD+ and NADH have the same number of phosphorus atoms)")
