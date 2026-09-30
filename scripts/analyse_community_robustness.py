"""
analyse_community_robustness.py
----------------------------------
Evaluates the stability of Louvain community detection on the
electron-transfer network and on the induced common-node metabolic
network (compare_networks.py), by repeating detection across multiple
random seeds and reporting:
  - the number of communities found at each seed;
  - agreement (adjusted Rand index) between each seed's partition and a
    reference (seed-0) partition of the SAME network;
  - cross-network agreement (electron vs metabolic) at each seed.

Uses build_induced_common_subgraphs from compare_networks.py, so the
induced-subgraph construction exists in exactly one place.

Usage
-----
    python scripts/analyse_community_robustness.py [--seeds N]
"""
import os
import sys
import argparse
SCRIPTS_DIR = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(SCRIPTS_DIR, ".."))
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, SCRIPTS_DIR)

import pandas as pd
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from etn.data_acquisition import load_ecoli_model
from etn.matching import build_metabolite_table, classify_reactions_molecular
from etn.network_analysis import deterministic_louvain_partition
from compare_networks import build_induced_common_subgraphs, partition_agreement

RESULTS_DIR = os.path.join(REPO_ROOT, "results")


def run_robustness_analysis(G_electron, G_metabolic, common_nodes, seeds) -> pd.DataFrame:
    common_list = sorted(common_nodes)

    partitions_electron, partitions_metabolic = {}, {}
    rows = []
    for seed in seeds:
        part_e = deterministic_louvain_partition(G_electron, random_state=seed)
        part_m = deterministic_louvain_partition(G_metabolic, random_state=seed)
        partitions_electron[seed] = part_e
        partitions_metabolic[seed] = part_m

        cross = partition_agreement(part_e, part_m, common_list)
        rows.append({
            "seed": seed,
            "n_communities_electron": len(set(part_e.values())),
            "n_communities_metabolic": len(set(part_m.values())),
            "cross_network_ARI": cross["adjusted_rand_index"],
            "cross_network_NMI": cross["normalized_mutual_information"],
        })

    df = pd.DataFrame(rows)
    reference_seed = seeds[0]
    df["ARI_electron_vs_reference_seed"] = [
        partition_agreement(partitions_electron[reference_seed], partitions_electron[s], common_list)["adjusted_rand_index"]
        for s in seeds
    ]
    df["ARI_metabolic_vs_reference_seed"] = [
        partition_agreement(partitions_metabolic[reference_seed], partitions_metabolic[s], common_list)["adjusted_rand_index"]
        for s in seeds
    ]
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, default=10, help="number of random seeds to test (0..seeds-1)")
    args = parser.parse_args()

    model = load_ecoli_model()
    met_table = build_metabolite_table(model)
    classification_rows = classify_reactions_molecular(model).to_dict("records")
    G_electron, G_metabolic, common_nodes = build_induced_common_subgraphs(model, classification_rows)

    seeds = list(range(args.seeds))
    df = run_robustness_analysis(G_electron, G_metabolic, common_nodes, seeds)
    df.to_csv(f"{RESULTS_DIR}/community_robustness.csv", index=False)

    # Stability of the full static electron network used for the full-network community analyses in the supplementary documentation.
    from etn.network_analysis import build_electron_graph, undirected_reaction_weighted_graph
    from sklearn.metrics import adjusted_rand_score
    G_full = undirected_reaction_weighted_graph(build_electron_graph(model, classification_rows))
    full_parts = {seed: deterministic_louvain_partition(G_full, random_state=seed) for seed in seeds}
    nodes_full = sorted(G_full.nodes())
    ref = full_parts[seeds[0]]
    full_rows=[]
    for seed in seeds:
        pseed=full_parts[seed]
        full_rows.append({
            "seed":seed,
            "n_communities":len(set(pseed.values())),
            "ARI_vs_reference_seed":adjusted_rand_score([ref[n] for n in nodes_full],[pseed[n] for n in nodes_full]),
        })
    full_df=pd.DataFrame(full_rows)
    full_df.to_csv(f"{RESULTS_DIR}/community_robustness_full_electron.csv", index=False)

    print(df.to_string(index=False))
    print("\nFull electron-network community stability:")
    print(full_df.to_string(index=False))
    print()
    print(f"Electron network communities: {df.n_communities_electron.min()}-{df.n_communities_electron.max()} "
          f"across {args.seeds} seeds; ARI vs reference seed: "
          f"{df.ARI_electron_vs_reference_seed.min():.3f}-{df.ARI_electron_vs_reference_seed.max():.3f}")
    print(f"Metabolic network communities: {df.n_communities_metabolic.min()}-{df.n_communities_metabolic.max()} "
          f"across {args.seeds} seeds; ARI vs reference seed: "
          f"{df.ARI_metabolic_vs_reference_seed.min():.3f}-{df.ARI_metabolic_vs_reference_seed.max():.3f}")
