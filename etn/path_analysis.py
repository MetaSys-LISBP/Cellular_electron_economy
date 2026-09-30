"""Static shortest-path diagnostics for the reconstructed electron-transfer graph.

Paths are calculated on the reduced-pool electron-transfer graph and provide
structural diagnostics for the selected shortest paths. Condition-specific
electron routing is analysed separately from pFBA-weighted flow balances.
"""

from typing import List
import networkx as nx
import pandas as pd

from .network_analysis import simple_weighted_graph


SUBSTRATE_PRODUCT_PAIRS = [
    ("nadh_c", "h2o_c", "NADH -> H2O (aerobic respiratory relay)"),
    ("nadh_c", "etoh_c", "NADH -> ethanol (fermentative relay)"),
    ("g6p_c", "nadph_c", "glucose-6-phosphate -> NADPH"),
    ("succ_c", "h2o_c", "succinate -> H2O (respiratory relay)"),
]


def path_to_readable(S: nx.DiGraph, path: List[str]) -> str:
    names = nx.get_node_attributes(S, "name")
    return " -> ".join(names.get(n, n) for n in path)


def shortest_paths_report(G, model, pairs=None) -> pd.DataFrame:
    """Report one shortest structural electron-transfer path per source/target pair."""
    pairs = pairs or SUBSTRATE_PRODUCT_PAIRS
    S = simple_weighted_graph(G)
    rows = []
    for source, target, label in pairs:
        present_s, present_t = source in S, target in S
        path, length, min_weight = None, None, None
        if present_s and present_t:
            try:
                path = nx.shortest_path(S, source, target)
                length = len(path) - 1
                min_weight = min(S[u][v]["weight"] for u, v in zip(path[:-1], path[1:]))
            except nx.NetworkXNoPath:
                pass
        rows.append({
            "label": label,
            "source": source, "target": target,
            "source_in_network": present_s, "target_in_network": present_t,
            "path_length_edges": length,
            "path_readable": path_to_readable(S, path) if path else None,
            "min_edge_electrons": min_weight,
        })
    return pd.DataFrame(rows)


def all_paths_report(G, source: str, target: str, cutoff: int = 5) -> pd.DataFrame:
    """Report all simple structural paths up to ``cutoff`` edges."""
    S = simple_weighted_graph(G)
    if source not in S or target not in S:
        return pd.DataFrame(columns=["path", "length", "min_electrons", "reactions"])
    rows = []
    try:
        for path in nx.all_simple_paths(S, source, target, cutoff=cutoff):
            weights = [S[u][v]["weight"] for u, v in zip(path[:-1], path[1:])]
            reactions = []
            for u, v in zip(path[:-1], path[1:]):
                reactions.extend(S[u][v]["reactions"])
            rows.append({
                "path": path_to_readable(S, path),
                "length": len(path) - 1,
                "min_electrons": min(weights),
                "reactions": ";".join(sorted(set(reactions))),
            })
    except nx.NetworkXNoPath:
        pass
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values("length").reset_index(drop=True)
    return df
