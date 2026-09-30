"""
network_analysis.py
--------------------
Builds and analyses the E. coli electron-transfer network from the
molecular classification (matching.py), and provides
the basic topological characterisations of the reconstructed network: hubs
(centrality), connected components, communities (Louvain modularity),
paths, and a bipartite reaction-metabolite representation used to test
community/subsystem enrichment.

TWO GRAPH REPRESENTATIONS, TWO PURPOSES
1. Projected donor-acceptor graph (`build_electron_graph` /
   `simple_weighted_graph`): nodes = metabolites, directed edges are weighted
   by electrons transferred. This is used for directed hub/centrality analyses.
   Analyses requiring an undirected graph use
   `undirected_reaction_weighted_graph`, which collapses direction explicitly
   while counting each distinct reaction once per unordered metabolite pair.
2. Bipartite reaction-metabolite graph (`build_bipartite_graph`): nodes
   are BOTH reactions and metabolites, edges connect a reaction to each
   metabolite it involves. This is the standard representation for
   relating network structure to reaction-level annotation (EC number,
   pathway/subsystem) since, unlike the projected graph, it preserves
   reaction identity. It is used here specifically for the community-
   versus-subsystem enrichment analysis (`subsystem_enrichment`), which
   needs each metabolite's set of associated reactions/subsystems.
Both are simple, widely used constructs in the constraint-based/network
modelling community; neither requires anything beyond graph connectivity
(no elementary-flux-mode or full FBA machinery), consistent with the
"basic, broadly accepted methods" scope of this study.

GRAPH CONSTRUCTION (projected graph)
The PRIMARY electron-transfer graph represents actual electron-transfer
relationships -- electron donor pool -> electron acceptor pool -- and
only for reactions whose donor and acceptor pairs are FULLY RESOLVED and
BALANCED (electrons donated == electrons accepted; see
matching.classify_reaction_molecular's
'resolution_status'). For each such reaction, an edge d.reactant ->
a.product is drawn, weighted by electrons, annotated with the
originating reaction. When a reaction has more than one donor and/or
acceptor pair (e.g. NODOy: two donor pairs, NADPH and NO, feeding one
O2/H2O acceptor pair), the total is distributed across every donor x
acceptor combination by PROPORTIONAL allocation, weight(d_i, a_j) =
|delta_d_i| * |delta_a_j| / total_electrons -- the unique allocation
rule that (a) reduces to the single-pair case exactly, (b) always sums
back to the reaction's true total regardless of how many donor/acceptor
pairs are present (unlike a naive min(d_i, a_j) cross-product, which
only conserves the total by coincidence in the specific N-donors-to-
1-acceptor pattern that happens to dominate this network, and would
silently double-count a reaction with multiple pairs on BOTH sides).

Reactions that are only partially resolved, or that resolve only one
side of the transfer, do not contribute to the primary topology graph: a
"who transfers to whom" edge cannot be drawn without inventing the
missing molecular identity. Their resolved and minimum-unresolved electron
content is retained in the reaction-level output and can be used by the
condition-specific analysis. Ambiguous formula-only mappings are likewise
excluded from the confident topology.
"""

from typing import Dict, List, Optional
import networkx as nx
import pandas as pd


# ---------------------------------------------------------------------
# 1. Projected donor-acceptor graph
# ---------------------------------------------------------------------

def build_electron_graph(model, classification_rows: List[dict], include_blocked_chemistry: bool = False) -> nx.MultiDiGraph:
    """
    classification_rows: output of
    matching.classify_reaction_molecular (one entry
    per reaction, with 'pairs' as a list of plain Python dicts, and
    'resolution_status' as computed there).

    Only 'fully_resolved' reactions (balanced donor/acceptor pairs, both
    genuinely identified) contribute to this primary topology graph.
    Reaction directionality follows the model bounds. When
    ``include_blocked_chemistry=True``, a fully resolved reaction with both model
    bounds fixed to zero is retained once in its written forward direction. This
    is used only for the complete chemistry-resolved reaction-edge registry and
    the boundary-resolved Supplementary Fig. S2 representation; topology and
    flux analyses retain the default model-direction-aware graph. Nodes represent the reduced
    electron-bearing species before and after each transfer: a forward edge is
    donor_pair.reactant -> acceptor_pair.product, whereas the chemically
    reversed edge is acceptor_pair.product -> donor_pair.reactant. Thus a
    reversible reaction contributes the same reduced-pool relationship in
    opposite directions, independent of its arbitrary written orientation.
    """
    met_name = {m.id: m.name for m in model.metabolites}
    rxn_by_id = {r.id: r for r in model.reactions}

    G = nx.MultiDiGraph()

    def _add_node(mid):
        if mid not in G:
            G.add_node(mid, name=met_name.get(mid, mid))

    for res in classification_rows:
        if not res["is_redox"] or res.get("resolution_status") != "fully_resolved":
            continue
        rxn = rxn_by_id[res["reaction_id"]]
        donors = [p for p in res["pairs"] if p["delta_gamma"] < 0]
        acceptors = [p for p in res["pairs"] if p["delta_gamma"] > 0]
        total_electrons = sum(abs(p["delta_gamma"]) for p in donors)  # == sum for acceptors, fully_resolved

        forward_allowed = rxn.upper_bound > 1e-9
        reverse_allowed = rxn.lower_bound < -1e-9
        blocked_chemistry = (not forward_allowed and not reverse_allowed)
        if include_blocked_chemistry and blocked_chemistry:
            forward_allowed = True
        for d in donors:
            for a in acceptors:
                electrons = abs(d["delta_gamma"]) * abs(a["delta_gamma"]) / total_electrons
                if forward_allowed:
                    _add_node(d["reactant"]); _add_node(a["product"])
                    G.add_edge(d["reactant"], a["product"], reaction=rxn.id,
                               reaction_name=rxn.name, subsystem=rxn.subsystem,
                               electrons=electrons, kind="complete", direction="forward")
                if reverse_allowed:
                    _add_node(a["product"]); _add_node(d["reactant"])
                    G.add_edge(a["product"], d["reactant"], reaction=rxn.id,
                               reaction_name=rxn.name, subsystem=rxn.subsystem,
                               electrons=electrons, kind="complete", direction="reverse")

    return G


def build_unresolved_electron_registry(model, classification_rows: List[dict]) -> pd.DataFrame:
    """
    Companion registry (NOT part of the primary topology graph) for redox reactions excluded from build_electron_graph because they
    are not fully resolved and balanced: their resolved cofactor
    identity (donor or acceptor, whichever side matched) and electron
    content are recorded here, rather than silently dropped or forced
    into a same-couple self-loop. Used for transparency/diagnostics and
    by the condition-specific flux-weighting analysis (which explicitly
    marks the missing side as "(unresolved organic)"
    rather than omitting it -- see etn/electron_flow_analysis.py).
    """
    rxn_by_id = {r.id: r for r in model.reactions}
    rows = []
    for res in classification_rows:
        if not res["is_redox"] or res.get("resolution_status") == "fully_resolved":
            continue
        rxn = rxn_by_id[res["reaction_id"]]
        for p in res["pairs"]:
            role = "donor" if p["delta_gamma"] < 0 else "acceptor"
            rows.append({
                "reaction_id": rxn.id, "reaction_name": rxn.name, "subsystem": rxn.subsystem,
                "resolution_status": res["resolution_status"], "role": role,
                "reactant": p["reactant"], "product": p["product"],
                "electrons": abs(p["delta_gamma"]),
            })
    return pd.DataFrame(rows)


def simple_weighted_graph(G: nx.MultiDiGraph) -> nx.DiGraph:
    """Projects the MultiDiGraph to a simple DiGraph (parallel edges
    merged; weight = summed electrons; n_reactions = number of distinct
    reactions contributing to that edge)."""
    S = nx.DiGraph()
    S.add_nodes_from(G.nodes(data=True))
    for u, v, data in G.edges(data=True):
        if S.has_edge(u, v):
            S[u][v]["weight"] += data["electrons"]
            S[u][v]["n_reactions"] += 1
            S[u][v]["reactions"].append(data["reaction"])
        else:
            S.add_edge(u, v, weight=data["electrons"], n_reactions=1, reactions=[data["reaction"]])
    return S




def undirected_reaction_weighted_graph(G: nx.MultiDiGraph) -> nx.Graph:
    """Collapse the directed reaction-resolved ETN to an undirected graph.

    This projection is used only for analyses that require an undirected
    network (currently Louvain community detection and the direct comparison
    with a conventional metabolite co-occurrence graph).  Each distinct
    reaction contributes its electron-transfer weight once to an unordered
    metabolite pair.  Thus the two opposite edges generated for a reversible
    reaction are not double-counted, whereas genuinely distinct reactions
    connecting the same metabolite pair are additive.

    For a given reaction and unordered metabolite pair, forward and reverse
    directional contributions are first summed separately.  If both directions
    are present (the usual reversible-reaction case), they must agree within
    numerical tolerance and one contribution is retained.  This makes the
    directed-to-undirected reduction explicit and independent of NetworkX edge
    insertion/conversion behaviour.
    """
    if not isinstance(G, nx.MultiDiGraph):
        raise TypeError("undirected_reaction_weighted_graph expects a networkx.MultiDiGraph")

    U = nx.Graph()
    U.add_nodes_from(G.nodes(data=True))

    # (unordered node pair, reaction) -> directional totals
    grouped = {}
    for u, v, data in G.edges(data=True):
        if u == v:
            pair = (u, v)
        else:
            pair = tuple(sorted((u, v)))
        reaction = data.get("reaction")
        if reaction is None:
            raise ValueError("Every reaction-resolved ETN edge must carry a reaction identifier")
        key = (pair, reaction)
        direction_key = (u, v)
        grouped.setdefault(key, {})
        grouped[key][direction_key] = grouped[key].get(direction_key, 0.0) + float(data.get("electrons", 0.0))

    pair_contributions = {}
    for (pair, reaction), directional in grouped.items():
        vals = list(directional.values())
        if len(vals) > 1:
            reference = vals[0]
            if not all(abs(v - reference) <= 1e-9 for v in vals[1:]):
                raise ValueError(
                    f"Directional electron weights differ for reversible reaction {reaction} on pair {pair}: {vals}"
                )
            contribution = reference
        else:
            contribution = vals[0]
        pair_contributions.setdefault(pair, []).append((reaction, contribution))

    for (u, v), contributions in sorted(pair_contributions.items()):
        total = sum(w for _, w in contributions)
        reactions = sorted(r for r, _ in contributions)
        U.add_edge(u, v, weight=total, n_reactions=len(reactions), reactions=reactions)
    return U


def _explicit_undirected_weighted_projection(G: nx.Graph) -> nx.Graph:
    """Return an explicit undirected weighted projection of a generic graph.

    For directed input lacking reaction-level identity, reciprocal directional
    weights are summed.  Electron-network analyses should instead use
    :func:`undirected_reaction_weighted_graph`, which can avoid double-counting
    reversible reactions because reaction identity is still available.
    """
    if not G.is_directed():
        return G.copy()
    U = nx.Graph()
    U.add_nodes_from(G.nodes(data=True))
    for u, v, data in G.edges(data=True):
        w = float(data.get("weight", 1.0))
        if U.has_edge(u, v):
            U[u][v]["weight"] += w
        else:
            U.add_edge(u, v, weight=w)
    return U


def basic_stats(G: nx.MultiDiGraph) -> dict:
    S = simple_weighted_graph(G)
    components = list(nx.weakly_connected_components(S))
    return {
        "n_nodes": S.number_of_nodes(),
        "n_edges_simple": S.number_of_edges(),
        "n_edges_multi": G.number_of_edges(),
        "n_reactions": len({d["reaction"] for _, _, d in G.edges(data=True)}),
        "density": nx.density(S),
        "n_weakly_connected_components": len(components),
        "largest_component_size": len(max(components, key=len)),
    }


def top_hubs(G: nx.MultiDiGraph, top_n: int = 15) -> pd.DataFrame:
    """Ranks metabolites by total (electron-weighted) degree, betweenness
    centrality, and PageRank."""
    S = simple_weighted_graph(G)
    names = nx.get_node_attributes(S, "name")

    degree = dict(S.degree(weight="weight"))
    betweenness = nx.betweenness_centrality(S, weight=lambda u, v, d: 1 / d["weight"])
    pagerank = nx.pagerank(S, weight="weight")

    df = pd.DataFrame({"metabolite_id": list(S.nodes())})
    df["name"] = df["metabolite_id"].map(names)
    df["weighted_degree_electrons"] = df["metabolite_id"].map(degree)
    df["betweenness"] = df["metabolite_id"].map(betweenness)
    df["pagerank"] = df["metabolite_id"].map(pagerank)
    df["in_degree"] = df["metabolite_id"].map(dict(S.in_degree()))
    df["out_degree"] = df["metabolite_id"].map(dict(S.out_degree()))
    return df.sort_values("weighted_degree_electrons", ascending=False).head(top_n).reset_index(drop=True)


def connected_components_report(G: nx.MultiDiGraph) -> pd.DataFrame:
    """Lists every weakly connected component with its size and members
    (smallest components first, to inspect what isolates them)."""
    S = simple_weighted_graph(G)
    names = nx.get_node_attributes(S, "name")
    components = list(nx.weakly_connected_components(S))
    rows = []
    for i, comp in enumerate(sorted(components, key=len)):
        rows.append({
            "component_id": i,
            "size": len(comp),
            "members": "; ".join(sorted(names.get(m, m) for m in comp)),
        })
    return pd.DataFrame(rows)



def deterministic_louvain_partition(G, random_state: int = 0, resolution: float = 1.0) -> dict:
    """Louvain partition reproducible across Python processes.

    String-node hash randomisation can otherwise change iteration order inside
    python-louvain even with a fixed random_state. Nodes are therefore mapped
    to sorted integer IDs and edges inserted in sorted order before clustering.
    The partition is mapped back to the original node IDs.
    """
    try:
        import community as community_louvain
    except ModuleNotFoundError:
        community_louvain = None
    H0 = _explicit_undirected_weighted_projection(G)
    nodes = sorted(H0.nodes())
    to_int = {n: i for i, n in enumerate(nodes)}
    H = nx.Graph()
    H.add_nodes_from(range(len(nodes)))
    edge_rows = []
    for u, v, d in H0.edges(data=True):
        a, b = to_int[u], to_int[v]
        if a > b: a, b = b, a
        edge_rows.append((a, b, float(d.get("weight", 1.0))))
    for a, b, w in sorted(edge_rows):
        if H.has_edge(a, b): H[a][b]["weight"] += w
        else: H.add_edge(a, b, weight=w)
    if H.number_of_edges() == 0:
        part_int = {i: 0 for i in H.nodes()}
    else:
        if community_louvain is not None:
            part_int = community_louvain.best_partition(H, weight="weight", random_state=random_state, resolution=resolution)
        else:
            communities = nx.community.louvain_communities(H, weight="weight", resolution=resolution, seed=random_state)
            communities = sorted((sorted(c) for c in communities), key=lambda c: (c[0] if c else -1, len(c)))
            part_int = {node: cid for cid, comm in enumerate(communities) for node in comm}
    return {n: part_int[to_int[n]] for n in nodes}

def detect_communities(G: nx.MultiDiGraph, random_state: int = 0, resolution: float = 1.0) -> pd.DataFrame:
    """Louvain community detection with deterministic node/edge ordering.

    ``random_state`` controls the Louvain heuristic and ``resolution`` is the
    standard Louvain resolution parameter (1.0 for the reported partition), while
    :func:`deterministic_louvain_partition` removes dependence on Python string
    hash/iteration order by relabelling sorted node IDs before clustering.
    """
    U = undirected_reaction_weighted_graph(G)
    partition = deterministic_louvain_partition(U, random_state=random_state, resolution=resolution)
    names = nx.get_node_attributes(U, "name")
    df = pd.DataFrame({"metabolite_id": list(partition.keys()),
                       "community": list(partition.values())})
    df["name"] = df["metabolite_id"].map(names)
    return df.sort_values(["community", "metabolite_id"]).reset_index(drop=True)


def find_paths(G: nx.MultiDiGraph, source: str, target: str, cutoff: int = 6):
    """All simple paths between two metabolites (default max length 6)."""
    S = simple_weighted_graph(G)
    if source not in S or target not in S:
        return []
    try:
        return list(nx.all_simple_paths(S, source, target, cutoff=cutoff))
    except nx.NetworkXNoPath:
        return []


# ---------------------------------------------------------------------
# 2. Bipartite reaction-metabolite graph
# ---------------------------------------------------------------------

def build_bipartite_graph(model, classification_rows: List[dict]) -> nx.Graph:
    """
    Bipartite graph: one node type for reactions (bipartite=0), one for
    metabolites (bipartite=1). All confidently redox reactions with matched
    molecular pairs are included, including reactions that are not fully
    resolved enough to contribute to the primary projected topology.  This
    broader registry is used only for annotation/diagnostic purposes, not to
    add edges to the confident electron-transfer topology.
    """
    B = nx.Graph()
    for res in classification_rows:
        if not res["is_redox"]:
            continue
        rid = res["reaction_id"]
        B.add_node(rid, bipartite=0)
        for pair in res["pairs"]:
            for mid in (pair["reactant"], pair["product"]):
                B.add_node(mid, bipartite=1)
                B.add_edge(rid, mid)
    return B


# ---------------------------------------------------------------------
# 3. Community vs. subsystem enrichment
# ---------------------------------------------------------------------

def subsystem_enrichment(model, classification_rows: List[dict], communities_df: pd.DataFrame) -> pd.DataFrame:
    """
    Tests, for every (community, subsystem) pair, whether the metabolites
    of that community are significantly enriched in reactions annotated
    to that subsystem, using a right-tailed Fisher exact test (metabolite
    is counted once per distinct subsystem it touches through any
    redox-classified reaction it participates in) with Benjamini-Hochberg
    FDR correction across all tests. A standard, simple enrichment
    approach (as commonly used to compare data-driven modules to curated
    pathway annotations, e.g. gene-set/pathway enrichment analysis).
    """
    from scipy.stats import fisher_exact
    from statsmodels.stats.multitest import multipletests

    rxn_by_id = {r.id: r for r in model.reactions}
    met_subsystems: Dict[str, set] = {}
    for res in classification_rows:
        if not res["is_redox"]:
            continue
        subsystem = rxn_by_id[res["reaction_id"]].subsystem or "Unknown"
        for pair in res["pairs"]:
            for mid in (pair["reactant"], pair["product"]):
                met_subsystems.setdefault(mid, set()).add(subsystem)

    all_mets = set(communities_df["metabolite_id"])
    all_subsystems = sorted({s for sset in met_subsystems.values() for s in sset})
    comm_of = dict(zip(communities_df["metabolite_id"], communities_df["community"]))

    rows = []
    N = len(all_mets)
    for comm in sorted(communities_df["community"].unique()):
        comm_mets = {m for m, c in comm_of.items() if c == comm}
        n_comm = len(comm_mets)
        if n_comm < 2:
            continue
        for subsystem in all_subsystems:
            sub_mets = {m for m in all_mets if subsystem in met_subsystems.get(m, set())}
            n_sub = len(sub_mets)
            if n_sub < 2:
                continue
            overlap = len(comm_mets & sub_mets)
            table = [[overlap, n_comm - overlap],
                     [n_sub - overlap, N - n_comm - n_sub + overlap]]
            _, p = fisher_exact(table, alternative="greater")
            rows.append({
                "community": comm, "subsystem": subsystem,
                "community_size": n_comm, "subsystem_size": n_sub,
                "overlap": overlap, "p_value": p,
            })

    df = pd.DataFrame(rows)
    if len(df):
        df["q_value"] = multipletests(df["p_value"], method="fdr_bh")[1]
        df = df.sort_values("p_value").reset_index(drop=True)
    return df
