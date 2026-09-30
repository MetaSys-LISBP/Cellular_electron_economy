"""Formula- and charge-based reconstruction of reaction-level electron transfer.

The algorithm uses only reaction stoichiometry, metabolite elemental formulae and
charges, plus a small seed list of universal redox-cofactor couples.  Matching is
performed globally within each reaction rather than greedily pair-by-pair.

Key design principles
---------------------
1. Exact transport of the same metabolite between compartments is removed first.
2. Universal seed redox couples are matched next, compartment-independently.
3. Remaining reactant/product occurrences are matched by maximum-cardinality,
   minimum-cost bipartite matching.  Formula-based candidates must conserve C,
   N, P and S and must share at least one element; this prevents chemically
   meaningless matches such as O2 <-> H+ without globally excluding small
   inorganic species (H+/H2, H2O/H2, formate/CO2, etc. remain possible).
4. Integer stoichiometric multiplicities are represented explicitly as repeated
   occurrences.  This resolves disproportionation reactions such as catalase and
   superoxide dismutase without reaction-specific overrides.
5. Formula-only 2x2 (and small n x n) reactions can be intrinsically ambiguous:
   a redox-looking assignment may compete with an equally stoichiometrically
   valid group-transfer assignment.  Such cases are not forced.  An internal
   redox-pair support library is inferred from reactions anchored by the universal
   seed cofactors.  This support is deliberately not recursively propagated.
   Ambiguous formula-only assignments are accepted only when their redox couples
   have this independent network support; otherwise the reaction is labelled ambiguous.

This preserves specificity without hard-coded reaction/pair blacklists.  The
ambiguity state is an explicit output because some mappings are not identifiable
from formula+charge alone.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Iterable, FrozenSet
from collections import defaultdict
import itertools
import math
import networkx as nx
import pandas as pd

from .degree_of_reduction import degree_of_reduction
from .reaction_electron_balance import compute_reaction_electron_balance

TOLERANCE = 1e-9
MAX_EXPLICIT_STOICH_OCCURRENCES = 8
MAX_AMBIGUITY_ENUMERATION = 7

# Canonical compartment-specific IDs are retained for API compatibility and reporting.
# Matching itself is compartment-independent and uses the corresponding base IDs.
KNOWN_REDOX_PAIRS = [
    ("nad_c", "nadh_c"), ("nadp_c", "nadph_c"),
    ("fad_c", "fadh2_c"), ("fmn_c", "fmnh2_c"),
    ("q8_c", "q8h2_c"), ("mqn8_c", "mql8_c"), ("2dmmq8_c", "2dmmql8_c"),
    ("o2_c", "h2o_c"), ("o2_c", "h2o2_c"),
    ("nh4_c", "no2_c"), ("no2_c", "no3_c"), ("no_c", "no2_c"),
    ("trdox_c", "trdrd_c"), ("gthox_c", "gthrd_c"),
    ("fdxo_2_2_c", "fdxr_2_2_c"), ("focytc_c", "ficytc_c"),
]


# Common compartment suffixes encountered in bacterial and eukaryotic GEMs.
# The E. coli c/p/e behaviour is unchanged; the additional labels enable the
# same compartment-independent seed-pair logic in mitochondrial and other
# compartmentalized reconstructions.
_COMPARTMENT_SUFFIXES = {"c", "p", "e", "m", "n", "r", "er", "erm", "lp", "v", "g", "ce", "mm", "x"}

def _metabolite_base_id(met_id: str) -> str:
    if "_" in met_id:
        head, tail = met_id.rsplit("_", 1)
        if head and tail in _COMPARTMENT_SUFFIXES:
            return head
    return met_id


_KNOWN_BASE_PAIRS: set[FrozenSet[str]] = {
    frozenset((_metabolite_base_id(a), _metabolite_base_id(b)))
    for a, b in KNOWN_REDOX_PAIRS
}


@dataclass(frozen=True)
class MetInfo:
    elements: Dict[str, int]
    charge: int
    gamma: Optional[float]
    ok: bool


@dataclass(frozen=True)
class _Occurrence:
    metabolite_id: str
    coefficient: float
    occurrence_index: int


@dataclass
class MatchingContext:
    """Network-wide context used to resolve formula-only ambiguities."""
    met_table: Dict[str, MetInfo]
    supported_redox_pairs: set[FrozenSet[str]]
    results_by_reaction: Dict[str, dict]


_CONTEXT_CACHE: Dict[int, MatchingContext] = {}


def clear_matching_context_cache() -> None:
    _CONTEXT_CACHE.clear()


def build_metabolite_table(model) -> Dict[str, MetInfo]:
    table: Dict[str, MetInfo] = {}
    for met in model.metabolites:
        try:
            els = met.elements
            g = degree_of_reduction(els, met.charge) if els else None
            table[met.id] = MetInfo(elements=els, charge=met.charge or 0,
                                    gamma=g, ok=g is not None)
        except Exception:
            table[met.id] = MetInfo(elements={}, charge=0, gamma=None, ok=False)
    return table


def _is_seed_pair(a: str, b: str) -> bool:
    return frozenset((_metabolite_base_id(a), _metabolite_base_id(b))) in _KNOWN_BASE_PAIRS


def _same_transport_species(a: str, b: str, met_table: Dict[str, MetInfo]) -> bool:
    ai, bi = met_table[a], met_table[b]
    return (_metabolite_base_id(a) == _metabolite_base_id(b)
            and ai.elements == bi.elements and ai.charge == bi.charge)


def _has_elemental_overlap(a: str, b: str, met_table: Dict[str, MetInfo],
                           *, include_hydrogen: bool = True) -> bool:
    ea, eb = met_table[a].elements, met_table[b].elements
    elements = set(ea) | set(eb)
    if not include_hydrogen:
        elements.discard("H")
    return any(min(ea.get(el, 0), eb.get(el, 0)) > 0 for el in elements)


def _strict_formula_candidate(a: str, b: str, met_table: Dict[str, MetInfo]) -> bool:
    """Redox-like formula correspondence.

    C/N/P/S are conserved within the pair; H/O/charge may change.  At least
    one atom type must be shared, preventing empty-signature cross-matches
    such as O2 -> H+ while retaining genuine small-inorganic chemistry.
    """
    if not _has_elemental_overlap(a, b, met_table, include_hydrogen=True):
        return False
    ea, eb = met_table[a].elements, met_table[b].elements
    return all(ea.get(el, 0) == eb.get(el, 0) for el in ("C", "N", "P", "S"))


def match_metabolite_pair(a: str, b: str, met_table: Dict[str, MetInfo]) -> Optional[dict]:
    """Match one reactant/product metabolite pair using the same local
    eligibility rules as the network matcher.

    This helper is intended for reference-pair validation and diagnostics,
    where the molecular correspondence is supplied explicitly. It does not
    resolve whole-reaction mapping ambiguity; that remains the responsibility
    of :func:`classify_reaction_molecular`.

    Returns ``None`` when the pair is not eligible, otherwise a dictionary
    containing the match source (``seed`` or ``formula``), ``delta_gamma``
    and electron count.
    """
    if a not in met_table or b not in met_table:
        return None
    ai, bi = met_table[a], met_table[b]
    if not ai.ok or not bi.ok or ai.gamma is None or bi.gamma is None:
        return None
    if _is_seed_pair(a, b):
        source = "seed"
    elif _strict_formula_candidate(a, b, met_table):
        source = "formula"
    else:
        return None
    delta = float(bi.gamma - ai.gamma)
    return {
        "reactant": a,
        "product": b,
        "match_source": source,
        "delta_gamma": delta,
        "electrons": abs(delta),
    }


def _formula_cost(a: str, b: str, met_table: Dict[str, MetInfo]) -> float:
    ai, bi = met_table[a], met_table[b]
    return (abs(ai.elements.get("H", 0) - bi.elements.get("H", 0))
            + abs(ai.elements.get("O", 0) - bi.elements.get("O", 0))
            + 0.1 * abs(ai.charge - bi.charge))


def _expand_occurrences(reaction, sign: int) -> List[_Occurrence]:
    """Expand small integer stoichiometric multiplicities into occurrences.

    This permits one chemical species to follow two redox branches in a
    disproportionation. Fractional coefficients (e.g. 0.5 O2) remain a single
    occurrence carrying their actual coefficient.
    """
    out: List[_Occurrence] = []
    for met, coef in reaction.metabolites.items():
        if coef * sign <= 0:
            continue
        magnitude = abs(float(coef))
        nearest = round(magnitude)
        if (abs(magnitude - nearest) <= TOLERANCE
                and 1 <= nearest <= MAX_EXPLICIT_STOICH_OCCURRENCES):
            unit_coef = 1.0 if coef > 0 else -1.0
            for i in range(int(nearest)):
                out.append(_Occurrence(met.id, unit_coef, i))
        else:
            out.append(_Occurrence(met.id, float(coef), 0))
    return out


def _maximum_cardinality_minimum_cost_matching(
    reactants: List[_Occurrence], products: List[_Occurrence],
    candidate_fn, met_table: Dict[str, MetInfo],
) -> List[Tuple[_Occurrence, _Occurrence]]:
    """Global bipartite matching: maximum cardinality, then minimum cost."""
    graph = nx.Graph()
    left = [("r", i) for i in range(len(reactants))]
    right = [("p", j) for j in range(len(products))]
    graph.add_nodes_from(left, bipartite=0)
    graph.add_nodes_from(right, bipartite=1)
    cardinality_bonus = 10000.0
    for i, r in enumerate(reactants):
        for j, p in enumerate(products):
            if candidate_fn(r.metabolite_id, p.metabolite_id):
                graph.add_edge(("r", i), ("p", j),
                               weight=cardinality_bonus - _formula_cost(
                                   r.metabolite_id, p.metabolite_id, met_table))
    matched = nx.algorithms.matching.max_weight_matching(
        graph, maxcardinality=True, weight="weight")
    result: List[Tuple[_Occurrence, _Occurrence]] = []
    for u, v in matched:
        if u[0] == "p":
            u, v = v, u
        if u[0] == "r" and v[0] == "p":
            result.append((reactants[u[1]], products[v[1]]))
    result.sort(key=lambda x: (x[0].metabolite_id, x[0].occurrence_index,
                               x[1].metabolite_id, x[1].occurrence_index))
    return result


def _remove_matched_occurrences(items: List[_Occurrence], matched: Iterable[_Occurrence]) -> List[_Occurrence]:
    ids = {id(x) for x in matched}
    return [x for x in items if id(x) not in ids]


def _seed_group_ratio(a: str, b: str, met_table: Dict[str, MetInfo]) -> tuple[int, int]:
    """Return the smallest occurrence ratio (n_a, n_b) for a seed pair.

    Most redox couples are 1:1.  Some universal couples change molecularity;
    oxidised/reduced glutathione is the canonical example (1 GSSG <-> 2 GSH).
    The ratio is inferred from the conserved heavy-element scaffold (C/N/P/S),
    so no pair-specific stoichiometric exception is required.  If the scaffold
    does not determine a ratio (e.g. O2/H2O), 1:1 is returned.
    """
    from fractions import Fraction
    ea, eb = met_table[a].elements, met_table[b].elements
    ratios = []
    for el in ("C", "N", "P", "S"):
        ca, cb = ea.get(el, 0), eb.get(el, 0)
        if ca == 0 and cb == 0:
            continue
        if ca == 0 or cb == 0:
            return (1, 1)
        ratios.append(Fraction(ca, cb))
    if not ratios or any(r != ratios[0] for r in ratios[1:]):
        return (1, 1)
    # n_a * count(a) == n_b * count(b), hence n_b/n_a = count(a)/count(b).
    r = ratios[0]
    return (r.denominator, r.numerator)


def _extract_nonunit_seed_groups(reactants, products, met_table):
    """Resolve seed pairs whose oxidised/reduced forms have non-1:1 molecularity.

    Returns synthetic seed matches plus the unconsumed occurrences.  Synthetic
    occurrences carry the grouped stoichiometric coefficient, allowing the
    normal gamma accounting to remain unchanged downstream.
    """
    grouped = []
    consumed_r, consumed_p = set(), set()
    # Work on actual metabolite-id pairs so compartments remain explicit.
    candidate_ids = sorted({(r.metabolite_id, p.metabolite_id)
                            for r in reactants for p in products
                            if _is_seed_pair(r.metabolite_id, p.metabolite_id)})
    for rid, pid in candidate_ids:
        nr, np = _seed_group_ratio(rid, pid, met_table)
        if (nr, np) == (1, 1):
            continue
        ravail = [r for r in reactants if r.metabolite_id == rid and id(r) not in consumed_r]
        pavail = [p for p in products if p.metabolite_id == pid and id(p) not in consumed_p]
        while len(ravail) >= nr and len(pavail) >= np:
            rs, ps = ravail[:nr], pavail[:np]
            rsynt = _Occurrence(rid, sum(x.coefficient for x in rs), rs[0].occurrence_index)
            psynt = _Occurrence(pid, sum(x.coefficient for x in ps), ps[0].occurrence_index)
            grouped.append((rsynt, psynt, "seed"))
            consumed_r.update(id(x) for x in rs)
            consumed_p.update(id(x) for x in ps)
            ravail, pavail = ravail[nr:], pavail[np:]
    reactants = [r for r in reactants if id(r) not in consumed_r]
    products = [p for p in products if id(p) not in consumed_p]
    return grouped, reactants, products


def _raw_match_occurrences(reaction, met_table: Dict[str, MetInfo]):
    reactants = _expand_occurrences(reaction, -1)
    products = _expand_occurrences(reaction, +1)

    # 0. Same metabolite transported between compartments.
    transport = _maximum_cardinality_minimum_cost_matching(
        reactants, products,
        lambda a, b: _same_transport_species(a, b, met_table), met_table)
    reactants = _remove_matched_occurrences(reactants, (r for r, _ in transport))
    products = _remove_matched_occurrences(products, (p for _, p in transport))

    # 1a. Seed pairs with non-unit molecularity (inferred from formula).
    grouped_seed_matches, reactants, products = _extract_nonunit_seed_groups(
        reactants, products, met_table)

    # 1b. Remaining universal redox seed pairs, compartment-independently.
    seed_matches = _maximum_cardinality_minimum_cost_matching(
        reactants, products, _is_seed_pair, met_table)
    reactants_after_seed = _remove_matched_occurrences(
        reactants, (r for r, _ in seed_matches))
    products_after_seed = _remove_matched_occurrences(
        products, (p for _, p in seed_matches))

    # 2. Global formula matching for all remaining occurrences.
    formula_matches = _maximum_cardinality_minimum_cost_matching(
        reactants_after_seed, products_after_seed,
        lambda a, b: _strict_formula_candidate(a, b, met_table), met_table)

    matches = (grouped_seed_matches
               + [(r, p, "seed") for r, p in seed_matches]
               + [(r, p, "formula") for r, p in formula_matches])
    return matches, reactants_after_seed, products_after_seed, formula_matches


def _aggregate_pair_contributions(matches, met_table: Dict[str, MetInfo]) -> List[dict]:
    aggregate: Dict[Tuple[str, str, str], float] = defaultdict(float)
    for r, p, source in matches:
        ri, pi = met_table[r.metabolite_id], met_table[p.metabolite_id]
        contribution = r.coefficient * ri.gamma + p.coefficient * pi.gamma
        if abs(contribution) > TOLERANCE:
            aggregate[(r.metabolite_id, p.metabolite_id, source)] += contribution
    rows = []
    for (r, p, source), contribution in sorted(aggregate.items()):
        if abs(contribution) > TOLERANCE:
            rows.append({"reactant": r, "product": p,
                         "match_source": source, "delta_gamma": contribution})
    return rows


def _pair_key(a: str, b: str) -> FrozenSet[str]:
    return frozenset((_metabolite_base_id(a), _metabolite_base_id(b)))


def _heavy_delta(a: str, b: str, met_table: Dict[str, MetInfo]) -> Tuple[int, int, int, int]:
    ea, eb = met_table[a].elements, met_table[b].elements
    return tuple(eb.get(el, 0) - ea.get(el, 0) for el in ("C", "N", "P", "S"))


def _has_plausible_group_transfer_alternative(reaction, met_table: Dict[str, MetInfo]) -> bool:
    """Detect a formula-level alternative assignment consistent with group transfer.

    Only small, square residual matching problems are enumerated. Alternative
    correspondences must retain at least one non-hydrogen atom in every pair;
    this rejects chemically nonsensical alternatives such as phosphonate->H2.
    """
    matches, residual_r, residual_p, strict_matches = _raw_match_occurrences(reaction, met_table)
    if (len(residual_r) != len(residual_p) or len(residual_r) < 2
            or len(residual_r) > MAX_AMBIGUITY_ENUMERATION):
        return False
    if len(strict_matches) != len(residual_r):
        return False

    from collections import Counter
    current = Counter((r.metabolite_id, p.metabolite_id) for r, p in strict_matches)
    for permutation in itertools.permutations(residual_p):
        alternative = list(zip(residual_r, permutation))
        alt_key = Counter((r.metabolite_id, p.metabolite_id) for r, p in alternative)
        if alt_key == current:
            continue
        if not all(_has_elemental_overlap(r.metabolite_id, p.metabolite_id,
                                          met_table, include_hydrogen=False)
                   for r, p in alternative):
            continue
        deltas = [_heavy_delta(r.metabolite_id, p.metabolite_id, met_table)
                  for r, p in alternative]
        if not any(any(v != 0 for v in d) for d in deltas):
            continue
        if all(sum(d[k] for d in deltas) == 0 for k in range(4)):
            return True
    return False


def _raw_classify_reaction(reaction, met_table: Dict[str, MetInfo]) -> dict:
    rid_upper = reaction.id.upper()
    if rid_upper.startswith("BIOMASS"):
        return {"reaction_id": reaction.id, "status": "pseudo_reaction",
                "missing_metabolites": "", "pairs": [], "is_redox": False,
                "redox_status": "not_applicable", "resolution_status": None}
    if reaction.id.startswith(("EX_", "DM_", "SK_")):
        return {"reaction_id": reaction.id, "status": "exchange_reaction",
                "missing_metabolites": "", "pairs": [], "is_redox": False,
                "redox_status": "not_applicable", "resolution_status": None}

    missing = [m.id for m in reaction.metabolites
               if not (met_table.get(m.id) and met_table[m.id].ok)]
    if missing:
        return {"reaction_id": reaction.id, "status": "missing_formula",
                "missing_metabolites": ";".join(missing), "pairs": [],
                "is_redox": False, "redox_status": "unscorable",
                "resolution_status": None}

    matches, _, _, _ = _raw_match_occurrences(reaction, met_table)
    pairs_info = _aggregate_pair_contributions(matches, met_table)
    balance = compute_reaction_electron_balance(pairs_info)
    return {
        "reaction_id": reaction.id,
        "status": "ok",
        "missing_metabolites": "",
        "pairs": pairs_info,
        "is_redox": bool(pairs_info),  # candidate status resolved by network-wide classification
        "redox_status": "candidate_redox" if pairs_info else "non_redox",
        "resolution_status": balance.resolution_status,
        "electrons_donated_resolved": balance.electrons_donated_resolved,
        "electrons_accepted_resolved": balance.electrons_accepted_resolved,
        "minimum_unresolved_donation": balance.minimum_unresolved_donation,
        "minimum_unresolved_acceptance": balance.minimum_unresolved_acceptance,
        "electrons_per_turnover_estimate": balance.electrons_per_turnover_estimate,
        "has_seed_redox_pair": any(p["match_source"] == "seed" for p in pairs_info),
    }


def _infer_supported_redox_pairs(model, raw_results: Dict[str, dict]) -> set[FrozenSet[str]]:
    """Infer formula-derived redox couples from independently seed-anchored reactions.

    Only balanced reactions containing at least one universal seed cofactor are
    used as evidence.  There is deliberately no recursive propagation: this
    prevents an erroneous formula-only assignment from becoming self-reinforcing
    elsewhere in the network.
    """
    support: set[FrozenSet[str]] = set()
    for res in raw_results.values():
        if not res.get("pairs") or not res.get("has_seed_redox_pair"):
            continue
        d, a = res["electrons_donated_resolved"], res["electrons_accepted_resolved"]
        if d <= TOLERANCE or abs(d - a) > TOLERANCE:
            continue
        for p in res["pairs"]:
            if p["match_source"] == "formula":
                support.add(_pair_key(p["reactant"], p["product"]))
    return support


def _finalise_result(reaction, raw: dict, met_table: Dict[str, MetInfo],
                     supported_pairs: set[FrozenSet[str]]) -> dict:
    result = dict(raw)
    pairs = result.get("pairs") or []
    if result.get("status") != "ok" or not pairs:
        return result

    d = result["electrons_donated_resolved"]
    a = result["electrons_accepted_resolved"]
    formula_keys = [_pair_key(p["reactant"], p["product"])
                    for p in pairs if p["match_source"] == "formula"]
    balanced = d > TOLERANCE and a > TOLERANCE and abs(d - a) <= TOLERANCE

    if result.get("has_seed_redox_pair"):
        confidence = "seed_anchored"
    elif balanced and formula_keys and any(k in supported_pairs for k in formula_keys):
        # One independently established pair is sufficient to anchor the
        # global assignment in a balanced reaction: once that molecular
        # correspondence is fixed, the competing cross-mapping is no longer
        # equally plausible.
        confidence = "network_anchored_ab_initio"
    elif balanced and _has_plausible_group_transfer_alternative(reaction, met_table):
        result.update({
            "status": "ambiguous_matching",
            "redox_status": "ambiguous",
            "is_redox": False,
            "resolution_status": "ambiguous",
            "matching_confidence": "ambiguous_formula_mapping",
        })
        return result
    elif balanced:
        confidence = "unique_ab_initio"
    else:
        # A unique one-sided/partially resolved gamma change is still direct
        # evidence of redox chemistry; the counterpart remains unresolved.
        confidence = "one_sided_ab_initio"

    result["redox_status"] = "redox"
    result["is_redox"] = True
    result["matching_confidence"] = confidence
    return result


def build_matching_context(model, met_table: Optional[Dict[str, MetInfo]] = None,
                           *, use_cache: bool = True) -> MatchingContext:
    key = id(model)
    if use_cache and key in _CONTEXT_CACHE:
        return _CONTEXT_CACHE[key]
    met_table = met_table or build_metabolite_table(model)
    raw_results = {rxn.id: _raw_classify_reaction(rxn, met_table)
                   for rxn in model.reactions}
    support = _infer_supported_redox_pairs(model, raw_results)
    final = {rxn.id: _finalise_result(rxn, raw_results[rxn.id], met_table, support)
             for rxn in model.reactions}
    context = MatchingContext(met_table=met_table,
                              supported_redox_pairs=support,
                              results_by_reaction=final)
    if use_cache:
        _CONTEXT_CACHE[key] = context
    return context


def classify_reaction_molecular(reaction, met_table: Dict[str, MetInfo],
                                context: Optional[MatchingContext] = None) -> dict:
    """Return the final reaction classification.

    When the reaction belongs to a model, the network-wide context is built once
    and cached so ambiguity resolution can use independently inferred redox-pair
    support.  For standalone reactions without a parent model, a conservative
    local result is returned: potentially ambiguous balanced formula-only cases
    are labelled ambiguous rather than forced.
    """
    if context is None:
        model = getattr(reaction, "model", None)
        if model is not None:
            context = build_matching_context(model, met_table)
    if context is not None and reaction.id in context.results_by_reaction:
        return dict(context.results_by_reaction[reaction.id])

    raw = _raw_classify_reaction(reaction, met_table)
    return _finalise_result(reaction, raw, met_table, set())


def classify_reactions_molecular(model) -> pd.DataFrame:
    context = build_matching_context(model)
    return pd.DataFrame([context.results_by_reaction[rxn.id] for rxn in model.reactions])
