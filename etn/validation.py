"""
validation.py
--------------
Publication-facing summary utilities for the reconstructed electron-transfer
network: matched-pair origin, ab initio redox-pair summaries, and redox-couple
usage statistics.  The 52-reaction categorical validation battery is defined
in ``scripts/validation_battery.py`` and is executed by
``scripts/run_validation.py``.
"""

from typing import List, Dict
import pandas as pd

from .matching import KNOWN_REDOX_PAIRS


def match_source_breakdown(classification_rows: List[dict]) -> pd.DataFrame:
    rows = []
    for res in classification_rows:
        if res.get("redox_status") == "ambiguous":
            continue
        for pair in res.get("pairs", []):
            rows.append({
                "reaction_id": res["reaction_id"],
                "reactant": pair["reactant"],
                "product": pair["product"],
                "source": "curated_cofactor_list" if pair.get("match_source") == "seed" else "ab_initio_formula_match",
                "electrons": abs(pair["delta_gamma"]),
            })
    df = pd.DataFrame(rows)
    return df


def ab_initio_redox_pairs(match_df: pd.DataFrame, met_name: Dict[str, str]) -> pd.DataFrame:
    """Deduplicated model-level metabolite pairs inferred by formula matching
    independently of the curated seed list. "Ab initio" identifies their
    algorithmic origin."""
    ab_initio = match_df[match_df.source == "ab_initio_formula_match"].copy()
    if len(ab_initio) == 0:
        return pd.DataFrame(columns=["reactant", "product", "reactant_name", "product_name",
                                      "n_reactions", "reactions", "mean_electrons"])
    grouped = ab_initio.groupby(["reactant", "product"]).agg(
        n_reactions=("reaction_id", "nunique"),
        reactions=("reaction_id", lambda s: ";".join(sorted(set(s)))),
        mean_electrons=("electrons", "mean"),
    ).reset_index()
    grouped["reactant_name"] = grouped["reactant"].map(met_name)
    grouped["product_name"] = grouped["product"].map(met_name)
    cols = ["reactant", "product", "reactant_name", "product_name",
            "n_reactions", "reactions", "mean_electrons"]
    return grouped[cols].sort_values("n_reactions", ascending=False).reset_index(drop=True)


def cofactor_usage_distribution(match_df: pd.DataFrame, met_name: Dict[str, str]) -> pd.DataFrame:
    """Number of distinct reactions in which each electron-carrier pair
    (curated or ab initio) is used."""
    grouped = match_df.groupby(["reactant", "product", "source"]).agg(
        n_reactions=("reaction_id", "nunique"),
        total_electrons=("electrons", "sum"),
    ).reset_index()
    grouped["pair_name"] = (grouped["reactant"].map(met_name).fillna(grouped["reactant"])
                             + " / " + grouped["product"].map(met_name).fillna(grouped["product"]))
    return grouped.sort_values("n_reactions", ascending=False).reset_index(drop=True)
