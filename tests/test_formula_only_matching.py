"""Dependency-light tests for the formula-only matching core.

These tests do not require COBRApy and therefore exercise the central chemistry
logic even in minimal CI environments.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from etn.matching import (build_metabolite_table, build_matching_context,
                          classify_reaction_molecular, clear_matching_context_cache)
from etn.electron_flow_analysis import signed_donor_acceptor_edges


def _parse_formula(formula):
    out = {}
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        out[element] = int(count or 1)
    return out


class Collection(list):
    def get_by_id(self, identifier):
        return next(x for x in self if x.id == identifier)


class Met:
    def __init__(self, identifier, formula, charge=0):
        self.id = identifier
        self.name = identifier
        self.elements = _parse_formula(formula)
        self.charge = charge

    def __hash__(self):
        return hash(self.id)


class Rxn:
    def __init__(self, identifier, stoich, lb=-1000.0, ub=1000.0):
        self.id = identifier
        self.name = identifier
        self.metabolites = stoich
        self.lower_bound = lb
        self.upper_bound = ub
        self.subsystem = "test"
        self.model = None


class Model:
    def __init__(self, metabolites, reactions):
        self.metabolites = Collection(metabolites)
        self.reactions = Collection(reactions)
        for reaction in reactions:
            reaction.model = self


def test_transamination_is_formula_ambiguous_without_blacklist():
    val = Met("val__L_c", "C5H11NO2", 0)
    akg = Met("akg_c", "C5H4O5", -2)
    keto = Met("3mob_c", "C5H7O3", -1)
    glu = Met("glu__L_c", "C5H8NO4", -1)
    # Deliberately use an arbitrary reaction id: ambiguity detection must be
    # chemistry-based, not a lookup of VALTA.
    rxn = Rxn("synthetic_transamination", {akg: -1, val: -1, keto: 1, glu: 1})
    model = Model([val, akg, keto, glu], [rxn])
    clear_matching_context_cache()
    ctx = build_matching_context(model, use_cache=False)
    result = ctx.results_by_reaction[rxn.id]
    assert result["redox_status"] == "ambiguous"
    assert not result["is_redox"]


def test_empty_signature_o2_to_proton_is_not_matched():
    substrate = Met("s_c", "C9H9O4", -1)
    product = Met("p_c", "C9H8O6", -2)
    oxygen = Met("o2_c", "O2", 0)
    proton = Met("h_c", "H", 1)
    rxn = Rxn("synthetic_oxygenase", {substrate: -1, oxygen: -1, product: 1, proton: 1})
    model = Model([substrate, product, oxygen, proton], [rxn])
    clear_matching_context_cache()
    result = build_matching_context(model, use_cache=False).results_by_reaction[rxn.id]
    assert result["is_redox"]
    assert result["resolution_status"] == "donor_side_only"
    assert not any({pair["reactant"], pair["product"]} == {"o2_c", "h_c"}
                   for pair in result["pairs"])


def test_disproportionation_resolved_from_stoichiometric_occurrences():
    peroxide = Met("h2o2_c", "H2O2", 0)
    water = Met("h2o_c", "H2O", 0)
    oxygen = Met("o2_c", "O2", 0)
    rxn = Rxn("synthetic_catalase", {peroxide: -2, water: 2, oxygen: 1})
    model = Model([peroxide, water, oxygen], [rxn])
    clear_matching_context_cache()
    result = build_matching_context(model, use_cache=False).results_by_reaction[rxn.id]
    assert result["is_redox"]
    assert result["resolution_status"] == "fully_resolved"
    pairs = {(p["reactant"], p["product"], p["delta_gamma"]) for p in result["pairs"]}
    assert ("h2o2_c", "o2_c", -2.0) in pairs
    assert ("h2o2_c", "h2o_c", 2.0) in pairs


def test_reverse_flux_uses_reverse_acceptor_species():
    etoh = Met("etoh_c", "C2H6O", 0)
    acald = Met("acald_c", "C2H4O", 0)
    # Actual iML-like formulas/charges are unnecessary for this identity test;
    # the seed ids establish the NAD pair and formulas establish gamma change.
    nad = Met("nad_c", "C21H26N7O14P2", -1)
    nadh = Met("nadh_c", "C21H27N7O14P2", -2)
    rxn = Rxn("synthetic_adh", {etoh: -1, nad: -1, acald: 1, nadh: 1})
    model = Model([etoh, acald, nad, nadh], [rxn])
    table = build_metabolite_table(model)
    clear_matching_context_cache()

    class Solution:
        fluxes = {"synthetic_adh": -5.0}

    edges = signed_donor_acceptor_edges(model, table, Solution())
    row = edges.iloc[0]
    assert row.donor == "nadh_c"
    assert row.acceptor == "etoh_c"


def test_match_metabolite_pair_public_helper():
    from etn.matching import MetInfo, match_metabolite_pair
    mt = {
        "fum_c": MetInfo({"C": 4, "H": 2, "O": 4}, -2, 12.0, True),
        "succ_c": MetInfo({"C": 4, "H": 4, "O": 4}, -2, 14.0, True),
        "atp_c": MetInfo({"C": 10, "H": 12, "N": 5, "O": 13, "P": 3}, -4, 30.0, True),
        "adp_c": MetInfo({"C": 10, "H": 12, "N": 5, "O": 10, "P": 2}, -3, 30.0, True),
    }
    m = match_metabolite_pair("fum_c", "succ_c", mt)
    assert m is not None
    assert m["match_source"] == "formula"
    assert m["delta_gamma"] == 2.0
    assert m["electrons"] == 2.0
    assert match_metabolite_pair("atp_c", "adp_c", mt) is None
