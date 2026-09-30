"""Reaction-centred summary of the molecular electron-transfer classification.

Electron counts and resolution status are taken from the central reaction-level
classification in :mod:`etn.matching` / :mod:`etn.reaction_electron_balance`.
Ambiguous formula-only mappings remain explicitly labelled and are never treated
as confident redox reactions by downstream topology or flux analyses.
"""

from typing import List
import pandas as pd


def build_electron_transfer_table(model, classification_rows: List[dict]) -> pd.DataFrame:
    """
    Parameters
    ----------
    model : cobrapy model (used to retrieve reaction/metabolite names)
    classification_rows : list of dicts, one per reaction, AS RETURNED BY
        matching.classify_reaction_molecular ('pairs'
        still a list of plain Python dicts, not JSON-serialised).
    """
    met_name = {m.id: m.name for m in model.metabolites}
    rxn_by_id = {r.id: r for r in model.reactions}

    def _uniq_join(items):
        seen = []
        for x in items:
            if x not in seen:
                seen.append(x)
        return ";".join(seen)

    rows = []
    for res in classification_rows:
        rxn = rxn_by_id[res["reaction_id"]]

        donor_names, acceptor_names = [], []
        donor_ids, acceptor_ids = [], []
        electrons_given = {}
        electrons_received = {}

        for pair in res["pairs"]:
            r_id, p_id, delta = pair["reactant"], pair["product"], pair["delta_gamma"]
            if delta < 0:
                donor_names.append(met_name.get(r_id, r_id))
                donor_ids.append(r_id)
                electrons_given["e"] = electrons_given.get("e", 0) + abs(delta)
            elif delta > 0:
                acceptor_names.append(met_name.get(p_id, p_id))
                acceptor_ids.append(p_id)
                electrons_received["e"] = electrons_received.get("e", 0) + delta

        total_electrons = res.get("electrons_per_turnover_estimate",
                                  max(sum(electrons_given.values()), sum(electrons_received.values())))
        status = res.get("redox_status",
                         res["status"] if res["status"] != "ok" else ("redox" if res["is_redox"] else "non_redox"))

        rows.append({
            "reaction_id": res["reaction_id"],
            "reaction_name": rxn.name,
            "reaction_equation": rxn.build_reaction_string(use_metabolite_names=True),
            "reversible": rxn.reversibility,
            "redox_status": status,
            "matching_confidence": res.get("matching_confidence", ""),
            "resolution_status": res.get("resolution_status", ""),
            "n_electrons_transferred": round(total_electrons, 2) if total_electrons else 0,
            "donor_metabolite_ids": _uniq_join(donor_ids),
            "donor_metabolite_names": _uniq_join(donor_names),
            "acceptor_metabolite_ids": _uniq_join(acceptor_ids),
            "acceptor_metabolite_names": _uniq_join(acceptor_names),
            "subsystem": rxn.subsystem,
        })

    columns = [
        "reaction_id", "reaction_name", "reaction_equation", "reversible",
        "redox_status", "matching_confidence", "resolution_status", "n_electrons_transferred",
        "donor_metabolite_ids", "donor_metabolite_names",
        "acceptor_metabolite_ids", "acceptor_metabolite_names",
        "subsystem",
    ]
    return pd.DataFrame(rows)[columns]
