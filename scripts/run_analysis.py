"""
run_analysis.py
------------------
Runs the full analysis pipeline: network reconstruction, electron-
conservation diagnostics, validation, static topology (hubs,
connectivity, communities, subsystem enrichment, condition-independent
in/out degree), static path diagnostics, and the experimentally constrained
aerobic/anaerobic pFBA case study. Writes every result table to
<repo_root>/results/.

Usage
-----
    python scripts/run_analysis.py

Works whether or not the `etn` package has been installed (`pip install -e .`):
if not installed, the repository root is added to sys.path automatically.
"""
import os
import sys
import json
import pandas as pd
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)  # allow running without `pip install -e .`

from etn.data_acquisition import acquire_all, export_sbml
from etn.matching import build_metabolite_table, classify_reactions_molecular
from etn.reaction_electron_balance import compute_reaction_electron_balance, validate_electron_balance_invariants
from etn.electron_transfer_table import build_electron_transfer_table
from etn.network_analysis import (build_electron_graph, simple_weighted_graph, basic_stats, top_hubs,
                                   connected_components_report, detect_communities,
                                   build_bipartite_graph, subsystem_enrichment)
from etn.validation import (match_source_breakdown, ab_initio_redox_pairs,
                             cofactor_usage_distribution)
from etn.path_analysis import shortest_paths_report

RESULTS_DIR = os.path.join(REPO_ROOT, "results")
DATA_DIR = os.path.join(REPO_ROOT, "data")  # optional local cache, see README
os.makedirs(RESULTS_DIR, exist_ok=True)


def main():
    # ---------------- Part (i): network acquisition ----------------
    print("=" * 70); print("PART (i): network acquisition"); print("=" * 70)
    # include_structures=False (default): the primary pipeline uses the
    # molecular formula and charge already present in the bundled iML1515.xml.
    # Structure-enriched companion analyses can use a local ModelSEED cache.
    model, reactions_df, metabolites_df = acquire_all()
    reactions_df.to_csv(f"{RESULTS_DIR}/reactions_ecoli.csv", index=False)
    metabolites_df.to_csv(f"{RESULTS_DIR}/metabolites_ecoli.csv", index=False)
    export_sbml(model, f"{RESULTS_DIR}/ecoli_network.xml")

    # ---------------- Part (ii): molecular electron-transfer classification ----------------
    print(); print("=" * 70); print("PART (ii): electron-transfer network reconstruction"); print("=" * 70)
    met_table = build_metabolite_table(model)

    # Machine-readable metabolite electron-content table used by the
    # Supplementary Data package. Unsupported chemistry is retained with a
    # missing gamma value rather than assigned an artificial value.
    gamma_table = metabolites_df.copy()
    gamma_table["gamma"] = gamma_table["met_id"].map(
        {mid: info.gamma for mid, info in met_table.items()}
    )
    gamma_table["chemical_domain_supported"] = gamma_table["met_id"].map(
        {mid: bool(info.ok) for mid, info in met_table.items()}
    )
    gamma_table.to_csv(f"{RESULTS_DIR}/metabolite_degree_of_reduction.csv", index=False)

    classification_rows = classify_reactions_molecular(model).to_dict("records")

    electron_table = build_electron_transfer_table(model, classification_rows)
    electron_table.to_csv(f"{RESULTS_DIR}/electron_transfer_table.tsv", sep="\t", index=False)
    electron_table.loc[electron_table["redox_status"] == "ambiguous"].to_csv(
        f"{RESULTS_DIR}/ambiguous_reactions.csv", index=False
    )
    electron_table.loc[electron_table["redox_status"] == "unscorable"].to_csv(
        f"{RESULTS_DIR}/unscorable_reactions.csv", index=False
    )
    n_redox = (electron_table.redox_status == "redox").sum()
    print(f"{n_redox} electron-transfer reactions identified out of {len(electron_table)} total.")
    print(electron_table["redox_status"].value_counts().to_string())

    # ---------------- Electron-conservation diagnostics ----------------
    print(); print("=" * 70); print("ELECTRON-CONSERVATION DIAGNOSTICS"); print("=" * 70)
    diagnostics_rows = []
    for res in classification_rows:
        if not res["is_redox"]:
            continue
        balance = compute_reaction_electron_balance(res["pairs"])
        validate_electron_balance_invariants(balance, res["reaction_id"])  # raises on violation
        diagnostics_rows.append({
            "reaction_id": res["reaction_id"],
            "resolution_status": balance.resolution_status,
            "electrons_donated_resolved": balance.electrons_donated_resolved,
            "electrons_accepted_resolved": balance.electrons_accepted_resolved,
            "minimum_unresolved_donation": balance.minimum_unresolved_donation,
            "minimum_unresolved_acceptance": balance.minimum_unresolved_acceptance,
            "electrons_per_turnover_estimate": balance.electrons_per_turnover_estimate,
        })
    diagnostics = pd.DataFrame(diagnostics_rows)
    diagnostics.to_csv(f"{RESULTS_DIR}/electron_balance_diagnostics.csv", index=False)
    print("Electron-conservation invariants verified for all "
          f"{len(diagnostics)} redox reactions.")
    print(diagnostics["resolution_status"].value_counts().to_string())

    # ---------------- Validation ----------------
    print(); print("=" * 70); print("VALIDATION"); print("=" * 70)
    met_name = {m.id: m.name for m in model.metabolites}
    match_df = match_source_breakdown(classification_rows)
    match_df.to_csv(f"{RESULTS_DIR}/matched_pairs_all.csv", index=False)
    print(match_df["source"].value_counts().to_string())

    ab_initio = ab_initio_redox_pairs(match_df, met_name)
    ab_initio.to_csv(f"{RESULTS_DIR}/ab_initio_redox_pairs.csv", index=False)
    print(f"Distinct model-level redox pairs inferred ab initio (not supplied as seeds): {len(ab_initio)}")

    cofactor_dist = cofactor_usage_distribution(match_df, met_name)
    cofactor_dist.to_csv(f"{RESULTS_DIR}/cofactor_usage_distribution.csv", index=False)

    # ---------------- Network construction & topology ----------------
    print(); print("=" * 70); print("NETWORK TOPOLOGY"); print("=" * 70)
    # The model-direction-aware graph is the topology used for all network
    # analyses. The complete chemistry-resolved graph additionally retains the
    # single fully resolved reaction that is present in iML1515 but blocked by
    # zero model bounds (FHL), matching the complete reaction-edge representation
    # shown in Supplementary Fig. S2.
    G = build_electron_graph(model, classification_rows)
    G_complete = build_electron_graph(model, classification_rows, include_blocked_chemistry=True)
    edge_rows = []
    for donor, acceptor, key, attrs in G_complete.edges(keys=True, data=True):
        edge_rows.append({
            "donor": donor,
            "acceptor": acceptor,
            "reaction_id": attrs.get("reaction"),
            "reaction_name": attrs.get("reaction_name"),
            "subsystem": attrs.get("subsystem"),
            "electrons_per_turnover": attrs.get("electrons"),
            "direction": attrs.get("direction"),
            "edge_kind": attrs.get("kind"),
            "model_direction_allowed": bool(
                model.reactions.get_by_id(attrs.get("reaction")).upper_bound > 1e-9
                if attrs.get("direction") == "forward"
                else model.reactions.get_by_id(attrs.get("reaction")).lower_bound < -1e-9
            ),
        })
    pd.DataFrame(edge_rows).to_csv(f"{RESULTS_DIR}/electron_transfer_edges.csv", index=False)

    stats = basic_stats(G)
    stats["n_reaction_edges_complete"] = int(G_complete.number_of_edges())
    stats["n_reactions_complete"] = int(len({d["reaction"] for _, _, d in G_complete.edges(data=True)}))
    stats["n_blocked_reaction_edges_complete"] = int(G_complete.number_of_edges() - G.number_of_edges())
    with open(f"{RESULTS_DIR}/network_basic_stats.json", "w") as f:
        json.dump(stats, f, indent=2)
    print(json.dumps(stats, indent=2))

    hubs = top_hubs(G, 30)
    hubs.to_csv(f"{RESULTS_DIR}/top_hubs.csv", index=False)

    # Full (all-node) weighted-degree table used for the Supplementary topology
    # degree-rank diagnostic, showing the complete rather than top-N distribution.
    all_hubs = top_hubs(G, G.number_of_nodes())
    all_hubs.to_csv(f"{RESULTS_DIR}/all_weighted_degree.csv", index=False)

    components = connected_components_report(G)
    components.to_csv(f"{RESULTS_DIR}/connected_components.csv", index=False)
    print(f"Connected components: {len(components)}")

    # Static (condition-independent) weighted in/out degree per node --
    # electron-count-weighted by the reconstructed network itself, never
    # by any condition-specific flux. Distinct from the flux-weighted
    # donor/acceptor activity computed later for specific growth
    # conditions (etn/electron_flow_analysis.py).
    S = simple_weighted_graph(G)
    met_name_map = {m.id: m.name for m in model.metabolites}
    topo_degree = pd.DataFrame({
        "metabolite_id": list(S.nodes()),
    })
    topo_degree["name"] = topo_degree["metabolite_id"].map(met_name_map)
    d_in = dict(S.in_degree(weight="weight"))
    d_out = dict(S.out_degree(weight="weight"))
    topo_degree["D_in"] = topo_degree["metabolite_id"].map(d_in)
    topo_degree["D_out"] = topo_degree["metabolite_id"].map(d_out)
    topo_degree["D_total"] = topo_degree["D_in"] + topo_degree["D_out"]
    topo_degree = topo_degree.sort_values("D_total", ascending=False).reset_index(drop=True)
    topo_degree.to_csv(f"{RESULTS_DIR}/topological_in_out_degree.csv", index=False)

    communities = detect_communities(G)
    communities.to_csv(f"{RESULTS_DIR}/communities.csv", index=False)
    print(f"Louvain communities: {communities['community'].nunique()}")

    enrichment = subsystem_enrichment(model, classification_rows, communities)
    enrichment.to_csv(f"{RESULTS_DIR}/community_subsystem_enrichment.csv", index=False)
    n_sig = (enrichment.q_value < 0.05).sum()
    print(f"Significant community-subsystem enrichments (q<0.05): {n_sig} / {len(enrichment)}")

    # bipartite graph (basic stats only, saved as edge list for reproducibility)
    B = build_bipartite_graph(model, classification_rows)
    nx_edges = pd.DataFrame(list(B.edges()), columns=["node_1", "node_2"])
    nx_edges.to_csv(f"{RESULTS_DIR}/bipartite_graph_edges.csv", index=False)
    print(f"Bipartite graph: {B.number_of_nodes()} nodes "
          f"({sum(1 for n, d in B.nodes(data=True) if d['bipartite'] == 0)} reactions, "
          f"{sum(1 for n, d in B.nodes(data=True) if d['bipartite'] == 1)} metabolites), "
          f"{B.number_of_edges()} edges")

    # ---------------- Path analysis ----------------
    print(); print("=" * 70); print("PATH ANALYSIS"); print("=" * 70)
    paths = shortest_paths_report(G, model)
    paths.to_csv(f"{RESULTS_DIR}/shortest_paths.csv", index=False)
    print(paths[["label", "path_length_edges", "path_readable"]].to_string(index=False))

    # ---------------- Experimentally constrained aerobic/anaerobic case study ----------------
    print(); print("=" * 70); print("FLUX-WEIGHTED ELECTRON NETWORK (aerobic vs anaerobic glucose growth)"); print("=" * 70)
    from etn.fba_case_study import load_experimental_data, run_case_study
    exp_data = load_experimental_data()
    exp_data.to_csv(f"{RESULTS_DIR}/aerobic_anaerobic_constraints_annotated.csv", index=False)
    fba_result = run_case_study(model, met_table, exp_data)
    fba_result["state_summary"].to_csv(f"{RESULTS_DIR}/fba_condition_summary.csv", index=False)
    fba_result["subsystem_comparison"].to_csv(f"{RESULTS_DIR}/fba_subsystem_comparison.csv")
    for cond, table in fba_result["reaction_tables"].items():
        table.to_csv(f"{RESULTS_DIR}/fba_redox_reactions_{cond}.csv", index=False)
    summary = fba_result["state_summary"].set_index("condition")
    with open(f"{RESULTS_DIR}/fba_summary_stats.json", "w") as f:
        json.dump({
            "aerobic_cumulative_electron_transfer_activity": float(summary.at["aerobic", "total_cumulative_electron_transfer_activity"]),
            "anaerobic_cumulative_electron_transfer_activity": float(summary.at["anaerobic", "total_cumulative_electron_transfer_activity"]),
            "aerobic_predicted_growth_rate": float(summary.at["aerobic", "predicted_growth_rate"]),
            "anaerobic_predicted_growth_rate": float(summary.at["anaerobic", "predicted_growth_rate"]),
        }, f, indent=2)
    print(summary[["qGLC","qO2","predicted_growth_rate","total_cumulative_electron_transfer_activity"]].to_string())

    print(); print("All results written to", RESULTS_DIR)


if __name__ == "__main__":
    main()
