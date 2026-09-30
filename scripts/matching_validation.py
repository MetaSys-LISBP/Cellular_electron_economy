"""
matching_validation.py
-------------------------
Validates the specificity of ab initio reactant-product matching
(etn/matching.py): confirms representative non-redox group-transfer
reactions remain non-redox, explicitly forbidden chemically invalid
inorganic pairs are absent, genuine small-inorganic redox chemistry is
retained, and curated seed-list redox pairs remain intact.
"""
import os
import sys
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)

import pandas as pd
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from etn.data_acquisition import load_ecoli_model
from etn.matching import build_metabolite_table, classify_reaction_molecular

RESULTS_DIR = os.path.join(REPO_ROOT, "results")

def _base(mid):
    return mid[:-2] if len(mid) > 2 and mid[-2] == "_" and mid[-1] in {"c", "p", "e"} else mid

# Representative group-transfer reactions with a formula-compatible
# alternative correspondence elsewhere in the same reaction, spanning the
# transformation classes for which spurious ab initio cross-matches are
# possible in principle.
GROUP_TRANSFER_REACTIONS = {
    "VALTA": "amino-group transfer (transamination)",
    "NDPK1": "phosphoryl transfer",
    "NDPK8": "phosphoryl transfer",
    "ACOXT": "CoA transfer",
    "DHAPT": "phosphoryl transfer",
    "CRNBTCT": "CoA transfer",
}


def check_group_transfer_reactions(model, met_table) -> pd.DataFrame:
    rows = []
    for rxn_id, category in GROUP_TRANSFER_REACTIONS.items():
        res = classify_reaction_molecular(model.reactions.get_by_id(rxn_id), met_table)
        rows.append({"reaction_id": rxn_id, "category": category,
                      "classified_redox": res["is_redox"], "correct": not res["is_redox"]})
    return pd.DataFrame(rows)


def check_aminotransferases(model, met_table) -> pd.DataFrame:
    """Every reaction whose name identifies it as an aminotransferase or
    transaminase must be classified non-redox."""
    rows = []
    for rxn in model.reactions:
        name = (rxn.name or "").lower()
        if "aminotransferase" in name or "transaminase" in name:
            res = classify_reaction_molecular(rxn, met_table)
            rows.append({"reaction_id": rxn.id, "reaction_name": rxn.name,
                         "classified_redox": res["is_redox"], "correct": not res["is_redox"]})
    return pd.DataFrame(rows)


def check_invalid_inorganic_cross_matches(model, met_table) -> pd.DataFrame:
    """Screens every accepted matched pair for chemically empty O2<->H+
    cross-matches. The rule is tested as an outcome, not implemented as a
    hard-coded pair blacklist.
    """
    rows = []
    for rxn in model.reactions:
        res = classify_reaction_molecular(rxn, met_table)
        if res.get("redox_status") == "ambiguous":
            continue
        for p in res["pairs"]:
            bases = {_base(p["reactant"]), _base(p["product"])}
            if bases == {"o2", "h"}:
                rows.append({"reaction_id": rxn.id, "reactant": p["reactant"], "product": p["product"]})
    return pd.DataFrame(rows)


def check_seed_pairs_intact(model, met_table) -> pd.DataFrame:
    """Confirms every curated seed-list pair still resolves correctly
    where it is expected to (a representative reaction per pair)."""
    from etn.matching import KNOWN_REDOX_PAIRS
    representative = {
        ("nad_c", "nadh_c"): "ALCD2x", ("nadp_c", "nadph_c"): "G6PDH2r",
        ("o2_c", "h2o_c"): "CYTBDpp", ("no3_c", "no2_c"): "NO3R1pp",
    }
    rows = []
    for pair, rxn_id in representative.items():
        res = classify_reaction_molecular(model.reactions.get_by_id(rxn_id), met_table)
        found = pair in {(p["reactant"], p["product"]) for p in res["pairs"]}
        rows.append({"pair": pair, "reaction_id": rxn_id, "resolved": found})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    model = load_ecoli_model()
    met_table = build_metabolite_table(model)

    group_transfer = check_group_transfer_reactions(model, met_table)
    aminotransferases = check_aminotransferases(model, met_table)
    invalid_inorganic = check_invalid_inorganic_cross_matches(model, met_table)
    seed_pairs = check_seed_pairs_intact(model, met_table)

    group_transfer.to_csv(f"{RESULTS_DIR}/matching_validation_group_transfer.csv", index=False)
    aminotransferases.to_csv(f"{RESULTS_DIR}/matching_validation_aminotransferases.csv", index=False)
    invalid_inorganic.to_csv(f"{RESULTS_DIR}/matching_validation_trivial_species.csv", index=False)

    print(f"Group-transfer reactions correctly non-redox: {group_transfer['correct'].sum()}/{len(group_transfer)}")
    print(f"Aminotransferase reactions correctly non-redox: {aminotransferases['correct'].sum()}/{len(aminotransferases)}")
    print(f"Explicitly invalid inorganic cross-matches: {len(invalid_inorganic)}")
    print(f"Seed-list pairs resolved as expected: {seed_pairs['resolved'].sum()}/{len(seed_pairs)}")

    assert group_transfer["correct"].all(), "A group-transfer reaction was classified redox"
    assert aminotransferases["correct"].all(), "An aminotransferase reaction was classified redox"
    assert len(invalid_inorganic) == 0, "An explicitly invalid inorganic pair was matched"
    assert seed_pairs["resolved"].all(), "A curated seed-list pair failed to resolve"
    print("\nAll matching-specificity checks passed.")
