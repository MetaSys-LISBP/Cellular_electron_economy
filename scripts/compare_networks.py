"""
compare_networks.py
-----------------------------------------
Quantitatively compares the FULL metabolic network
(a conventional metabolite co-occurrence network constructed from
non-structural reactions after explicit removal of documented bookkeeping
species) against the electron-transfer network reconstructed in
this study, restricted to the metabolites the two have in common. This
goes beyond the single Fisher-exact-test community/subsystem enrichment
already reported, by comparing the two networks' COMMUNITY STRUCTURE
directly against each other (not against a curated annotation), using
standard, simple, widely used metrics: Jaccard similarity of edge sets,
degree-distribution comparison, and Adjusted Rand Index / normalised
mutual information between their respective Louvain partitions.
"""
import os
import sys
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)

import numpy as np
import pandas as pd
import networkx as nx
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from etn.data_acquisition import load_ecoli_model
from etn.matching import build_metabolite_table, classify_reactions_molecular
from etn.network_analysis import build_electron_graph, undirected_reaction_weighted_graph, deterministic_louvain_partition

RESULTS_DIR = os.path.join(REPO_ROOT, "results")

# Exact compartment-specific bookkeeping species omitted from the conventional
# co-occurrence graph.  Redox carriers (NAD(H), NADP(H), FAD(H2), quinones) are
# intentionally retained because removing them would erase the chemistry that
# the electron-transfer network is designed to compare.
METABOLIC_NETWORK_BOOKKEEPING_IDS = frozenset({
    "h_c", "h_p", "h_e", "h2o_c", "h2o_p", "h2o_e",
    "atp_c", "adp_c", "amp_c", "pi_c", "pi_p", "ppi_c",
    "co2_c", "co2_p", "nh4_c", "o2_c", "o2_p", "coa_c",
})


def build_metabolic_network(model, exclude_bookkeeping=True):
    """
    Standard metabolite co-occurrence network: an (undirected) edge
    between two metabolites if they appear together in the same reaction
    (excluding BIOMASS/exchange/demand/sink reactions, for a fair
    comparison against the electron network, which excludes those too).
    `exclude_bookkeeping=True` additionally drops the exact compartment-specific
    bookkeeping species listed in ``METABOLIC_NETWORK_BOOKKEEPING_IDS`` (h_c, h2o_c, atp_c, adp_c, pi_c, co2_c...) that
    would otherwise dominate degree through their ubiquity in stoichiometric
    bookkeeping rather than a specific pathway relationship. The exact filter
    is documented here and in the Supplementary Information.
    """
    G = nx.Graph()
    for rxn in model.reactions:
        if rxn.id.upper().startswith("BIOMASS") or rxn.id.startswith(("EX_", "DM_", "SK_")):
            continue
        mets = [m.id for m in rxn.metabolites]
        if exclude_bookkeeping:
            mets = [m for m in mets if m not in METABOLIC_NETWORK_BOOKKEEPING_IDS]
        for i in range(len(mets)):
            for j in range(i + 1, len(mets)):
                if G.has_edge(mets[i], mets[j]):
                    G[mets[i]][mets[j]]["weight"] += 1
                else:
                    G.add_edge(mets[i], mets[j], weight=1)
    return G


def jaccard_edge_similarity(G1: nx.Graph, G2: nx.Graph, common_nodes: set) -> float:
    e1 = {frozenset(e) for e in G1.edges() if e[0] in common_nodes and e[1] in common_nodes}
    e2 = {frozenset(e) for e in G2.edges() if e[0] in common_nodes and e[1] in common_nodes}
    if not e1 and not e2:
        return float("nan")
    return len(e1 & e2) / len(e1 | e2)


def partition_agreement(part1: dict, part2: dict, common_nodes: list) -> dict:
    """Adjusted Rand Index and Normalised Mutual Information between two
    community partitions, restricted to nodes present in both."""
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    labels1 = [part1[n] for n in common_nodes]
    labels2 = [part2[n] for n in common_nodes]
    return {
        "adjusted_rand_index": adjusted_rand_score(labels1, labels2),
        "normalized_mutual_information": normalized_mutual_info_score(labels1, labels2),
        "n_common_nodes": len(common_nodes),
    }


def build_induced_common_subgraphs(model, classification_rows):
    """
    Builds the electron-transfer network and the full metabolic
    (co-occurrence) network, then returns both INDUCED to their common
    node set -- required before any comparison (edge Jaccard, degree
    correlation, community detection) is computed, since a node's
    community assignment depends on its complete neighbourhood: a
    partition of the full metabolic graph is not interchangeable with a
    partition of the graph restricted to only the nodes shared with the
    electron network. Shared by compare_networks.py and
    analyse_community_robustness.py so the induced-subgraph construction
    exists in exactly one place.

    Returns (G_electron_induced, G_metabolic_induced, common_nodes).
    """
    G_electron_multi = build_electron_graph(model, classification_rows)
    G_electron_full = undirected_reaction_weighted_graph(G_electron_multi)
    G_metabolic_full = build_metabolic_network(model)
    common_nodes = set(G_electron_full.nodes()) & set(G_metabolic_full.nodes())
    G_electron = G_electron_full.subgraph(common_nodes).copy()
    G_metabolic = G_metabolic_full.subgraph(common_nodes).copy()
    return G_electron, G_metabolic, common_nodes


if __name__ == "__main__":
    model = load_ecoli_model()
    met_table = build_metabolite_table(model)
    classification_rows = classify_reactions_molecular(model).to_dict("records")

    print("Building electron-transfer and metabolic networks, induced to their common node set...")
    G_electron, G_metabolic, common_nodes = build_induced_common_subgraphs(model, classification_rows)
    print(f"Common nodes: {len(common_nodes)}")
    print(f"Induced electron subgraph: {G_electron.number_of_nodes()} nodes, {G_electron.number_of_edges()} edges")
    print(f"Induced metabolic subgraph: {G_metabolic.number_of_nodes()} nodes, {G_metabolic.number_of_edges()} edges")

    jaccard = jaccard_edge_similarity(G_electron, G_metabolic, common_nodes)
    print(f"\nEdge Jaccard similarity (induced common-node subgraphs): {jaccard:.4f}")

    # degree distribution comparison, on the induced subgraphs
    deg_electron = dict(G_electron.degree(weight="weight"))
    deg_metabolic = dict(G_metabolic.degree(weight="weight"))
    deg_df = pd.DataFrame({
        "metabolite_id": sorted(common_nodes),
        "degree_electron_network": [deg_electron.get(n, 0) for n in sorted(common_nodes)],
        "degree_metabolic_network": [deg_metabolic.get(n, 0) for n in sorted(common_nodes)],
    })
    from scipy.stats import spearmanr
    rho, pval = spearmanr(deg_df["degree_electron_network"], deg_df["degree_metabolic_network"])
    print(f"Degree Spearman correlation (induced subgraphs): rho={rho:.3f}, P={pval:.2e}")
    deg_df.to_csv(f"{RESULTS_DIR}/degree_comparison_electron_vs_metabolic.csv", index=False)

    # community partition comparison -- Louvain run directly on each
    # INDUCED subgraph (not on the full graph with post-hoc filtering)
    print("\nDetecting communities on the induced common-node subgraphs (Louvain, random_state=0)...")
    part_electron = deterministic_louvain_partition(G_electron, random_state=0)
    part_metabolic = deterministic_louvain_partition(G_metabolic, random_state=0)

    common_list = sorted(common_nodes)
    agreement = partition_agreement(part_electron, part_metabolic, common_list)
    print(f"\nCommunity partition agreement (induced common-node subgraphs):")
    for k, v in agreement.items():
        print(f"  {k}: {v}")

    n_comm_electron = len(set(part_electron.values()))
    n_comm_metabolic = len(set(part_metabolic.values()))
    print(f"\nNumber of communities (induced common-node subgraphs): electron network = {n_comm_electron}, "
          f"metabolic network = {n_comm_metabolic}")

    import json
    summary = {
        "n_nodes_electron_induced": G_electron.number_of_nodes(),
        "n_edges_electron_induced": G_electron.number_of_edges(),
        "n_nodes_metabolic_induced": G_metabolic.number_of_nodes(),
        "n_edges_metabolic_induced": G_metabolic.number_of_edges(),
        "n_common_nodes": len(common_nodes),
        "jaccard_edge_similarity": jaccard,
        "degree_spearman_rho": rho,
        "degree_spearman_pvalue": pval,
        "n_communities_electron": n_comm_electron,
        "n_communities_metabolic": n_comm_metabolic,
        **agreement,
    }
    with open(f"{RESULTS_DIR}/network_comparison_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print("\nSummary written to results/network_comparison_summary.json")
