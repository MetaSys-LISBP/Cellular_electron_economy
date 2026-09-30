"""
test_validation.py
-------------------
Regression tests for the publication-facing validation and core reduction-degree logic. Uses the bundled, local iML1515 model
(data/iML1515.xml) -- no network access required.

Run with:
    pytest tests/
"""
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pytest
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")


def test_52_reaction_validation_battery_perfect_accuracy(model, met_table):
    import importlib.util
    from pathlib import Path
    script = Path(__file__).resolve().parents[1] / "scripts" / "validation_battery.py"
    spec = importlib.util.spec_from_file_location("validation_battery", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    df = module.run_validation_battery(model, met_table)
    scorable = df[df["ground_truth"].isin(["redox", "non_redox"])]
    unscorable = df[df["ground_truth"].isin(["redox_unscorable", "non_redox_unscorable"])]
    assert len(df) == 52
    assert len(scorable) == 48
    assert len(unscorable) == 4
    assert df["correct"].all()


def test_degree_of_reduction_reference_values():
    from etn.degree_of_reduction import degree_of_reduction
    assert degree_of_reduction({"C": 6, "H": 12, "O": 6}, 0) == 24          # glucose
    assert degree_of_reduction({"C": 2, "H": 6, "O": 1}, 0) == 12           # ethanol
    # acetate and acetic acid must have identical gamma (protonation is not redox)
    g_acetate = degree_of_reduction({"C": 2, "H": 3, "O": 2}, -1)
    g_acetic_acid = degree_of_reduction({"C": 2, "H": 4, "O": 2}, 0)
    assert g_acetate == g_acetic_acid == 8
    # NAD+ -> NADH must be a 2-electron transfer once the charge term is included
    g_nad = degree_of_reduction({"C": 21, "H": 26, "N": 7, "O": 14, "P": 2}, -1)
    g_nadh = degree_of_reduction({"C": 21, "H": 27, "N": 7, "O": 14, "P": 2}, -2)
    assert g_nadh - g_nad == 2


def test_untracked_elements_return_none():
    from etn.degree_of_reduction import degree_of_reduction
    # a formula containing selenium must not silently return a value
    assert degree_of_reduction({"C": 10, "H": 17, "N": 3, "O": 6, "S": 1, "Se": 1}, -1) is None


def test_pseudo_and_exchange_reactions_excluded(model, met_table):
    from etn.matching import classify_reaction_molecular
    biomass_rxn = next(r for r in model.reactions if r.id.upper().startswith("BIOMASS"))
    result = classify_reaction_molecular(biomass_rxn, met_table)
    assert result["status"] == "pseudo_reaction"
    assert result["is_redox"] is False

    exchange_rxn = model.reactions.get_by_id("EX_glc__D_e")
    result = classify_reaction_molecular(exchange_rxn, met_table)
    assert result["status"] == "exchange_reaction"


def test_single_sided_flux_directionality(model, met_table):
    """Regression test for direction-aware donor/acceptor assignment when
    only ONE side of a redox reaction is formula-matched (e.g. GLUDy,
    glutamate dehydrogenase, whose deamination changes N count so the
    organic pair is not matched): reversing the flux must flip which
    FUNCTIONAL ROLE (donor vs acceptor) the resolved cofactor pair plays
    -- not just swap which member (reactant/product) of that pair is
    reported under an unchanged role label. GLUDy is written
    glu__L_c + nadp_c <=> akg_c + nadph_c + nh4_c (forward: glutamate
    oxidised, NADP+ reduced -- an ACCEPTOR pair). Under reverse flux, the
    true reaction is reductive amination, in which NADPH -- not NADP+ --
    donates electrons; the resolved acceptor pair must therefore be
    reported as a DONOR (identity = nadph_c) when flux is negative."""
    from etn.matching import classify_reaction_molecular
    from etn.electron_flow_analysis import signed_donor_acceptor_edges

    rxn = model.reactions.get_by_id("GLUDy")
    res = classify_reaction_molecular(rxn, met_table)
    assert len(res["pairs"]) == 1 and res["pairs"][0]["delta_gamma"] > 0, (
        "This test assumes GLUDy resolves only a single ACCEPTOR pair "
        "(nadp_c->nadph_c); if the matching algorithm changes and now "
        "also resolves the organic (glutamate/2-oxoglutarate) pair, this "
        "pair structure does not match the expected validation fixture.")

    class MockSolution:
        def __init__(self, fluxes):
            self.fluxes = fluxes

    edges_fwd = signed_donor_acceptor_edges(model, met_table, MockSolution({"GLUDy": 5.0}))
    row_fwd = edges_fwd[edges_fwd.reaction_id == "GLUDy"].iloc[0]
    assert row_fwd.donor == "(unresolved organic)"
    assert row_fwd.acceptor == "nadph_c"

    edges_rev = signed_donor_acceptor_edges(model, met_table, MockSolution({"GLUDy": -5.0}))
    row_rev = edges_rev[edges_rev.reaction_id == "GLUDy"].iloc[0]
    assert row_rev.donor == "nadph_c", (
        f"Expected NADPH (not NADP+) as the true donor under reverse flux, got {row_rev.donor!r}")
    assert row_rev.acceptor == "(unresolved organic)"

