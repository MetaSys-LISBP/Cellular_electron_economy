"""Sensitivity of the electron/metabolic network comparison to bookkeeping removal.

The primary comparison removes 18 explicitly documented compartment-specific
bookkeeping species from the conventional co-occurrence graph.  This script
repeats the comparison with no such removal, testing whether the conclusion
that the electron graph is topologically distinct depends on that choice.
"""
import json
import os
import sys

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)

import pandas as pd
from scipy.stats import spearmanr
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from etn.data_acquisition import load_ecoli_model
from etn.matching import classify_reactions_molecular
from etn.network_analysis import build_electron_graph, undirected_reaction_weighted_graph, deterministic_louvain_partition
from scripts.compare_networks import build_metabolic_network, jaccard_edge_similarity, partition_agreement

RESULTS_DIR = os.path.join(REPO_ROOT, "results")


def metrics(model, classification_rows, exclude_bookkeeping):
    ge = undirected_reaction_weighted_graph(build_electron_graph(model, classification_rows))
    gm = build_metabolic_network(model, exclude_bookkeeping=exclude_bookkeeping)
    common = set(ge.nodes()) & set(gm.nodes())
    ge = ge.subgraph(common).copy(); gm = gm.subgraph(common).copy()
    common_sorted = sorted(common)
    de = dict(ge.degree(weight="weight")); dm = dict(gm.degree(weight="weight"))
    rho, p = spearmanr([de[n] for n in common_sorted], [dm[n] for n in common_sorted])
    pe = deterministic_louvain_partition(ge, random_state=0)
    pm = deterministic_louvain_partition(gm, random_state=0)
    a = partition_agreement(pe, pm, common_sorted)
    return {
        "bookkeeping_removed": bool(exclude_bookkeeping),
        "n_common_nodes": len(common),
        "n_electron_edges": ge.number_of_edges(),
        "n_metabolic_edges": gm.number_of_edges(),
        "edge_jaccard": jaccard_edge_similarity(ge, gm, common),
        "weighted_degree_spearman_rho": float(rho),
        "weighted_degree_spearman_p": float(p),
        "adjusted_rand_index": float(a["adjusted_rand_index"]),
        "normalized_mutual_information": float(a["normalized_mutual_information"]),
        "n_electron_communities": len(set(pe.values())),
        "n_metabolic_communities": len(set(pm.values())),
    }


def main():
    model = load_ecoli_model()
    classification_rows = classify_reactions_molecular(model).to_dict("records")
    out = pd.DataFrame([
        metrics(model, classification_rows, True),
        metrics(model, classification_rows, False),
    ])
    out.to_csv(os.path.join(RESULTS_DIR, "network_comparison_sensitivity.csv"), index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
