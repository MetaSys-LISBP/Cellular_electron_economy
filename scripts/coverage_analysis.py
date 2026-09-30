"""
coverage_analysis.py
-------------------
Characterises exactly why each reaction of iML1515 is or is not
classified for electron transfer, broken down by precise exclusion
reason (transition metal, selenium, generic formula group, no
formula at all, transport/exchange/demand/sink, biomass), and reporting
coverage both as a fraction of ALL reactions and as a fraction of
CHEMICALLY ANALYSABLE reactions (i.e. excluding the structural
categories -- transport/exchange/demand/sink/biomass -- that were never
candidates for a redox call in the first place, since they either
involve a single metabolite by construction or are a modelling
artefact, not an enzymatic transformation).
"""
import os
import sys
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)

import pandas as pd
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from etn.data_acquisition import load_ecoli_model
from etn.degree_of_reduction import SAFE_ELEMENTS
from etn.matching import build_metabolite_table, classify_reaction_molecular

RESULTS_DIR = os.path.join(REPO_ROOT, "results")

TRANSITION_METALS = {"Fe", "Cu", "Mo", "Co", "Ni", "Mn", "Zn", "Ag", "Cd", "Hg", "W"}
METALLOIDS = {"Se", "As"}
GENERIC_FORMULA_SYMBOLS = {"X", "R"}


def classify_unscorable_reason(met) -> str:
    """For a metabolite that failed to score, identifies WHY: which
    specific element category caused it."""
    els = met.elements or {}
    if not els:
        return "no_formula"
    bad = [e for e in els if e not in SAFE_ELEMENTS]
    if any(e in TRANSITION_METALS for e in bad):
        return "transition_metal"
    if any(e in METALLOIDS for e in bad):
        return "metalloid"
    if any(e in GENERIC_FORMULA_SYMBOLS for e in bad):
        return "generic_placeholder"
    if bad:
        return f"other_element({','.join(sorted(bad))})"
    return "scorable"


def analyse_coverage(model, met_table):
    rows = []
    for rxn in model.reactions:
        rid_upper = rxn.id.upper()
        if rid_upper.startswith("BIOMASS"):
            rows.append({"reaction_id": rxn.id, "exclusion_category": "biomass_pseudoreaction"})
            continue
        if rxn.id.startswith("EX_"):
            rows.append({"reaction_id": rxn.id, "exclusion_category": "exchange"})
            continue
        if rxn.id.startswith("DM_"):
            rows.append({"reaction_id": rxn.id, "exclusion_category": "demand"})
            continue
        if rxn.id.startswith("SK_"):
            rows.append({"reaction_id": rxn.id, "exclusion_category": "sink"})
            continue

        res = classify_reaction_molecular(rxn, met_table)
        if res.get("redox_status") == "ambiguous":
            rows.append({"reaction_id": rxn.id, "exclusion_category": "ambiguous_formula_mapping"})
            continue
        if res["status"] == "ok":
            rows.append({"reaction_id": rxn.id,
                         "exclusion_category": "redox" if res["is_redox"] else "non_redox_scorable"})
            continue

        # status == "missing_formula": find the precise reason(s)
        reasons = set()
        for met in rxn.metabolites:
            r = classify_unscorable_reason(met)
            if r != "scorable":
                reasons.add(r)
        if not reasons:
            reasons = {"unknown"}
        rows.append({"reaction_id": rxn.id, "exclusion_category": "+".join(sorted(reasons))})

    df = pd.DataFrame(rows)
    return df


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    counts = df["exclusion_category"].value_counts().rename_axis("category").reset_index(name="n_reactions")
    counts["fraction_of_all_reactions"] = counts["n_reactions"] / len(df)

    # "chemically analysable" = everything except the structural categories
    structural = {"biomass_pseudoreaction", "exchange", "demand", "sink"}
    n_analysable = len(df[~df["exclusion_category"].isin(structural)])
    counts["fraction_of_analysable_reactions"] = counts.apply(
        lambda r: (r["n_reactions"] / n_analysable) if r["category"] not in structural else float("nan"), axis=1)
    return counts, n_analysable


if __name__ == "__main__":
    model = load_ecoli_model()
    met_table = build_metabolite_table(model)
    df = analyse_coverage(model, met_table)
    df.to_csv(f"{RESULTS_DIR}/coverage_per_reaction.csv", index=False)

    summary, n_analysable = summarise(df)
    summary.to_csv(f"{RESULTS_DIR}/coverage_summary.csv", index=False)

    print(f"Total reactions: {len(df)}")
    print(f"Chemically analysable reactions (excludes biomass/exchange/demand/sink): {n_analysable}")
    print()
    print(summary.to_string(index=False))

    n_redox = (df.exclusion_category == "redox").sum()
    n_nonredox = (df.exclusion_category == "non_redox_scorable").sum()
    print()
    print(f"Coverage (fraction of ALL reactions successfully classified redox/non-redox): "
          f"{(n_redox + n_nonredox) / len(df):.1%}")
    print(f"Coverage (fraction of CHEMICALLY ANALYSABLE reactions successfully classified): "
          f"{(n_redox + n_nonredox) / n_analysable:.1%}")
