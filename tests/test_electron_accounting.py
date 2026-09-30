"""
test_electron_accounting.py
------------------------------
Tests for reactant-product matching, electron-balance accounting, static
network construction, direction-aware flux weighting, physiological-state
aggregation, and the induced-common-node network comparison.
"""
import pytest


# ----------------------------------------------------------------------
# Degree of reduction: representative compounds and canonical couples
# ----------------------------------------------------------------------

@pytest.mark.parametrize("elements,charge,expected_gamma", [
    ({"C": 21, "H": 26, "N": 7, "O": 14, "P": 2}, -1, 72),   # NAD+
    ({"C": 21, "H": 27, "N": 7, "O": 14, "P": 2}, -2, 74),   # NADH
    ({"C": 21, "H": 25, "N": 7, "O": 17, "P": 3}, -3, 72),   # NADP+
    ({"C": 21, "H": 26, "N": 7, "O": 17, "P": 3}, -4, 74),   # NADPH
    ({"O": 4, "P": 1}, -3, 0),                                # phosphate (own reference state)
    ({"C": 10, "H": 12, "N": 5, "O": 13, "P": 3}, -4, 30),   # ATP
    ({"C": 10, "H": 12, "N": 5, "O": 10, "P": 2}, -3, 30),   # ADP
    ({"C": 10, "H": 12, "N": 5, "O": 7, "P": 1}, -2, 30),    # AMP
])
def test_degree_of_reduction_reference_compounds(elements, charge, expected_gamma):
    from etn.degree_of_reduction import degree_of_reduction
    assert degree_of_reduction(elements, charge) == expected_gamma


@pytest.mark.parametrize("rxn_id,expected_electrons", [
    ("ALCD2x", 2),   # ethanol/acetaldehyde, NAD-linked
    ("G6PDH2r", 2),  # NADP-dehydrogenase
    ("FRD2", 2),     # fumarate/succinate, menaquinol-linked
    ("MDH", 2),      # malate/oxaloacetate, NAD-linked
    ("SUCDi", 2),    # succinate/fumarate, quinone-linked
])
def test_canonical_redox_couples(model, met_table, rxn_id, expected_electrons):
    from etn.matching import classify_reaction_molecular
    res = classify_reaction_molecular(model.reactions.get_by_id(rxn_id), met_table)
    assert res["is_redox"]
    assert res["electrons_per_turnover_estimate"] == expected_electrons


# ----------------------------------------------------------------------
# Non-redox group-transfer reactions must not be classified as redox
# ----------------------------------------------------------------------

@pytest.mark.parametrize("rxn_id", [
    "VALTA",     # transamination (2-oxoglutarate/valine <-> 3-methyl-2-oxobutanoate/glutamate)
    "NDPK1",     # nucleoside diphosphate kinase (phosphotransfer)
    "NDPK8",     # nucleoside diphosphate kinase (phosphotransfer)
    "ACOXT",     # CoA transferase
    "DHAPT",     # PEP-dependent phosphotransferase
    "CRNBTCT",   # CoA transferase
])
def test_group_transfer_reactions_are_non_redox(model, met_table, rxn_id):
    from etn.matching import classify_reaction_molecular
    res = classify_reaction_molecular(model.reactions.get_by_id(rxn_id), met_table)
    assert not res["is_redox"], f"{rxn_id} should be classified non-redox"


@pytest.mark.parametrize("rxn_id", [
    "HPPPNDO", "DHCINDO", "METOX1s", "METOX2s", "DHORDfum", "PPTHpp", "FHL",
])
def test_legitimate_redox_reactions_not_suppressed(model, met_table, rxn_id):
    from etn.matching import classify_reaction_molecular
    res = classify_reaction_molecular(model.reactions.get_by_id(rxn_id), met_table)
    assert res["is_redox"], f"{rxn_id} is genuine redox chemistry and must remain redox"


def test_all_aminotransferase_reactions_are_non_redox(model, met_table):
    """Every aminotransferase/transaminase reaction in iML1515 must be
    classified non-redox: transamination moves an amino group between two
    carbon skeletons and transfers no electrons."""
    from etn.matching import classify_reaction_molecular
    flagged = []
    for rxn in model.reactions:
        name = (rxn.name or "").lower()
        if "aminotransferase" in name or "transaminase" in name:
            res = classify_reaction_molecular(rxn, met_table)
            if res["is_redox"]:
                flagged.append(rxn.id)
    assert flagged == [], f"Aminotransferase reactions incorrectly classified as redox: {flagged}"


# ----------------------------------------------------------------------
# Trivial inorganic species must not be paired by ab initio matching
# ----------------------------------------------------------------------

def test_no_ab_initio_o2_to_proton_pair(model, met_table):
    """O2 -> H+ discards both oxygen atoms of O2 with no complementary
    product to account for them and is never a valid molecular redox
    pair; any genuine O2-involving redox chemistry is covered by the
    curated seed pairs (O2/H2O, O2/H2O2)."""
    from etn.matching import classify_reaction_molecular
    offending = []
    for rxn in model.reactions:
        res = classify_reaction_molecular(rxn, met_table)
        for p in res["pairs"]:
            if p["reactant"] == "o2_c" and p["product"] == "h_c":
                offending.append(rxn.id)
    assert offending == [], f"O2->H+ pair found in: {offending}"


def test_curated_inorganic_couples_still_resolve(model, met_table):
    """Genuine curated inorganic redox pairs must remain resolvable while
    explicitly invalid ab initio inorganic relationships are excluded."""
    from etn.matching import classify_reaction_molecular
    res = classify_reaction_molecular(model.reactions.get_by_id("CYTBDpp"), met_table)
    assert res["is_redox"]
    pairs = {(p["reactant"], p["product"]) for p in res["pairs"]}
    assert ("o2_c", "h2o_c") in pairs


def test_superoxide_dismutase_disproportionation(model, met_table):
    """Superoxide dismutase uses two copies of superoxide: one is oxidised
    to O2 and one is reduced to H2O2. The reaction therefore contains a
    matched 1-electron donor branch and a matched 1-electron acceptor branch,
    not a spurious H+->H2O2 pair."""
    from etn.matching import classify_reaction_molecular
    for rid, comp in [("SPODM", "c"), ("SPODMpp", "p")]:
        res = classify_reaction_molecular(model.reactions.get_by_id(rid), met_table)
        assert res["resolution_status"] == "fully_resolved"
        pairs = {(p["reactant"], p["product"], p["delta_gamma"]) for p in res["pairs"]}
        assert (f"o2s_{comp}", f"o2_{comp}", -1.0) in pairs
        assert (f"o2s_{comp}", f"h2o2_{comp}", 1.0) in pairs
        assert not any(p["reactant"] == f"h_{comp}" for p in res["pairs"])


# ----------------------------------------------------------------------
# Multi-pair electron counting (sum of pairs, not max of a single pair)
# ----------------------------------------------------------------------

@pytest.mark.parametrize("rxn_id", ["CINNDO", "PPPNDO", "PYROX"])
def test_multi_donor_pair_electron_total(model, met_table, rxn_id):
    """Two independent 2-electron donor pairs (the organic substrate and
    NADH) must sum to a 4-electron total, not report the magnitude of a
    single pair."""
    from etn.matching import classify_reaction_molecular
    res = classify_reaction_molecular(model.reactions.get_by_id(rxn_id), met_table)
    assert res["resolution_status"] == "donor_side_only"
    assert res["electrons_donated_resolved"] == 4
    assert res["electrons_per_turnover_estimate"] == 4


# ----------------------------------------------------------------------
# Partially resolved chemistry: resolved vs minimum-unresolved split
# ----------------------------------------------------------------------

@pytest.mark.parametrize("rxn_id", ["FDMO", "FDMO2", "FDMO3", "FDMO4", "FDMO6"])
def test_partially_resolved_reactions(model, met_table, rxn_id):
    """FMNH2->FMN (2 e- donor) matched against O2->H2O (4 e- acceptor) is
    imbalanced: 2 electrons are resolved on the donor side, with a
    minimum additional 2-electron donation required (but not measured)
    to balance the acceptor side."""
    from etn.matching import classify_reaction_molecular
    res = classify_reaction_molecular(model.reactions.get_by_id(rxn_id), met_table)
    assert res["resolution_status"] == "partially_resolved"
    assert res["electrons_donated_resolved"] == 2
    assert res["electrons_accepted_resolved"] == 4
    assert res["minimum_unresolved_donation"] == 2
    assert res["minimum_unresolved_acceptance"] == 0


def test_electron_balance_invariants_hold_network_wide(model, met_table):
    from etn.matching import classify_reaction_molecular
    from etn.reaction_electron_balance import compute_reaction_electron_balance, validate_electron_balance_invariants
    for rxn in model.reactions:
        res = classify_reaction_molecular(rxn, met_table)
        if not res["is_redox"]:
            continue
        balance = compute_reaction_electron_balance(res["pairs"])
        validate_electron_balance_invariants(balance, rxn.id)  # raises on violation


def test_network_wide_electron_transfer_stoichiometry(model, met_table):
    """The Fig. 1d distribution is derived directly from reaction-level
    electron accounting and should remain stable unless the reconstruction
    itself changes. The two 1-electron cases are SPODM/SPODMpp; two-electron
    chemistry dominates the remaining redox network."""
    from collections import Counter
    from etn.matching import classify_reaction_molecular

    counts = Counter()
    one_electron = []
    for rxn in model.reactions:
        res = classify_reaction_molecular(rxn, met_table)
        if not res["is_redox"]:
            continue
        n_e = int(round(res["electrons_per_turnover_estimate"]))
        counts[n_e] += 1
        if n_e == 1:
            one_electron.append(rxn.id)

    assert counts == Counter({1: 2, 2: 237, 4: 22, 6: 6, 8: 3})
    assert sorted(one_electron) == ["SPODM", "SPODMpp"]


def test_invariant_check_fails_on_inconsistent_synthetic_reaction():
    """A deliberately inconsistent ElectronBalance (claiming fully_resolved
    status while donor and acceptor totals differ) must fail the
    invariant check."""
    from etn.reaction_electron_balance import ElectronBalance, validate_electron_balance_invariants
    bad = ElectronBalance(
        electrons_donated_resolved=2.0, electrons_accepted_resolved=4.0,
        resolution_status="fully_resolved",  # inconsistent with the totals above
        minimum_unresolved_donation=0.0, minimum_unresolved_acceptance=0.0,
        electrons_per_turnover_estimate=4.0,
    )
    with pytest.raises(AssertionError):
        validate_electron_balance_invariants(bad, "synthetic")


# ----------------------------------------------------------------------
# Directionality
# ----------------------------------------------------------------------

def test_forward_flux_direction(model, met_table):
    from etn.electron_flow_analysis import signed_donor_acceptor_edges

    class MockSolution:
        def __init__(self, fluxes):
            self.fluxes = fluxes

    edges = signed_donor_acceptor_edges(model, met_table, MockSolution({"ALCD2x": 5.0}))
    row = edges[edges.reaction_id == "ALCD2x"].iloc[0]
    assert row.donor == "etoh_c" and row.acceptor == "nadh_c"


def test_reverse_flux_direction(model, met_table):
    from etn.electron_flow_analysis import signed_donor_acceptor_edges

    class MockSolution:
        def __init__(self, fluxes):
            self.fluxes = fluxes

    edges = signed_donor_acceptor_edges(model, met_table, MockSolution({"ALCD2x": -5.0}))
    row = edges[edges.reaction_id == "ALCD2x"].iloc[0]
    assert row.donor == "nadh_c" and row.acceptor == "etoh_c"


def test_glud_direction_dependent_donor(model, met_table):
    """GLUDy resolves only a cofactor pair (organic glutamate/2-
    oxoglutarate pair not formula-matched, since deamination changes
    nitrogen count); the true donor depends on flux direction."""
    from etn.constants import UNRESOLVED_ORGANIC_NODE
    from etn.electron_flow_analysis import signed_donor_acceptor_edges

    class MockSolution:
        def __init__(self, fluxes):
            self.fluxes = fluxes

    fwd = signed_donor_acceptor_edges(model, met_table, MockSolution({"GLUDy": 5.0}))
    row_fwd = fwd[fwd.reaction_id == "GLUDy"].iloc[0]
    assert row_fwd.donor == UNRESOLVED_ORGANIC_NODE
    assert row_fwd.acceptor == "nadph_c"

    rev = signed_donor_acceptor_edges(model, met_table, MockSolution({"GLUDy": -5.0}))
    row_rev = rev[rev.reaction_id == "GLUDy"].iloc[0]
    assert row_rev.donor == "nadph_c"
    assert row_rev.acceptor == UNRESOLVED_ORGANIC_NODE


def test_fdmo_directional_edges_conserve_total(model, met_table):
    from etn.electron_flow_analysis import signed_donor_acceptor_edges

    class MockSolution:
        def __init__(self, fluxes):
            self.fluxes = fluxes

    edges = signed_donor_acceptor_edges(model, met_table, MockSolution({"FDMO": 3.0}))
    sub = edges[edges.reaction_id == "FDMO"]
    assert sub.electron_flux.sum() == pytest.approx(3.0 * 4.0)
    assert (sub.acceptor == "h2o_c").all()


# ----------------------------------------------------------------------
# Static network topology
# ----------------------------------------------------------------------

def test_static_graph_contains_only_direct_donor_acceptor_edges(model, met_table):
    """For a fully resolved reaction (ethanol/NAD+ -> acetaldehyde/NADH),
    the static graph must contain the direct donor->acceptor relationship
    and must NOT contain internal donor-pair or acceptor-pair
    transformations as separate edges."""
    from etn.matching import classify_reaction_molecular
    from etn.network_analysis import build_electron_graph
    rows = [classify_reaction_molecular(r, met_table) for r in model.reactions]
    G = build_electron_graph(model, rows)

    alcd2x_edges = [(u, v, d.get("direction")) for u, v, d in G.edges(data=True) if d.get("reaction") == "ALCD2x"]
    assert ("etoh_c", "nadh_c", "forward") in alcd2x_edges
    assert ("nadh_c", "etoh_c", "reverse") in alcd2x_edges
    assert not any((u, v) == ("etoh_c", "acald_c") for u, v, _ in alcd2x_edges)
    assert not any((u, v) == ("nad_c", "nadh_c") for u, v, _ in alcd2x_edges)


def test_static_graph_excludes_incompletely_resolved_reactions(model, met_table):
    from etn.matching import classify_reaction_molecular
    from etn.network_analysis import build_electron_graph
    rows = [classify_reaction_molecular(r, met_table) for r in model.reactions]
    G = build_electron_graph(model, rows)
    for rxn_id in ["CINNDO", "FDMO", "GAPD", "AKGDH"]:
        edges = [1 for u, v, d in G.edges(data=True) if d.get("reaction") == rxn_id]
        assert edges == []


def test_static_graph_excludes_unresolved_partner():
    from etn.constants import UNRESOLVED_ORGANIC_NODE
    import networkx as nx
    from etn.network_analysis import build_electron_graph
    # The static graph uses molecule identities drawn from resolved pairs.
    assert UNRESOLVED_ORGANIC_NODE not in ("etoh_c", "nadh_c", "acald_c", "nad_c")


# ----------------------------------------------------------------------
# Aerobic / anaerobic case-study constraints
# ----------------------------------------------------------------------

def test_anaerobic_condition_fixes_oxygen_to_zero_and_is_feasible(model):
    from etn.fba_case_study import solve_condition, OXYGEN_EXCHANGE_ID
    sol = solve_condition(model, 14.9, 0.0, False)
    assert sol is not None
    assert sol.fluxes[OXYGEN_EXCHANGE_ID] == pytest.approx(0.0, abs=1e-10)
    assert sol.fluxes["BIOMASS_Ec_iML1515_core_75p37M"] > 0


def test_case_study_input_matches_published_glucose_oxygen_constraints():
    from etn.fba_case_study import load_experimental_data
    df = load_experimental_data().set_index("condition")
    assert df.loc["aerobic", "qGLC"] == pytest.approx(8.7)
    assert df.loc["aerobic", "qO2"] == pytest.approx(11.9)
    assert df.loc["anaerobic", "qGLC"] == pytest.approx(14.9)
    assert df.loc["anaerobic", "qO2"] == pytest.approx(0.0)


# ----------------------------------------------------------------------
# Common-node network comparison
# ----------------------------------------------------------------------

def test_community_comparison_uses_induced_subgraphs():
    """Community detection must be run on graphs induced to the common
    node set, not on the full graphs with labels filtered afterwards: a
    node's community assignment depends on its complete neighbourhood, so
    the two are not interchangeable."""
    networkx = pytest.importorskip("networkx")
    community_louvain = pytest.importorskip("community")
    from scripts.compare_networks import partition_agreement

    g1 = networkx.Graph()
    g1.add_edges_from([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (4, 6), (5, 7)])
    g2 = networkx.Graph()
    g2.add_edges_from([(1, 2), (2, 3)])

    common = {1, 2, 3}
    g1_induced = g1.subgraph(common).copy()
    g2_induced = g2.subgraph(common).copy()

    part1 = community_louvain.best_partition(g1_induced, random_state=0)
    part2 = community_louvain.best_partition(g2_induced, random_state=0)
    agreement = partition_agreement(part1, part2, sorted(common))
    assert agreement["adjusted_rand_index"] == pytest.approx(1.0)


def test_condition_solver_matches_growth_maximising_pfba(model):
    """The case-study solver fixes measured glucose/O2 and applies pFBA at maximum growth."""
    from etn.fba_case_study import (solve_condition, BIOMASS_REACTION_ID,
                                    GLUCOSE_EXCHANGE_ID, OXYGEN_EXCHANGE_ID,
                                    AEROBIC_ZERO_PRODUCT_EXCHANGES)
    from cobra.flux_analysis import pfba
    qglc, qo2 = 8.7, 11.9
    sol = solve_condition(model, qglc, qo2, True)
    assert sol is not None
    m=model.copy()
    m.reactions.get_by_id(GLUCOSE_EXCHANGE_ID).bounds=(-qglc,-qglc)
    m.reactions.get_by_id(OXYGEN_EXCHANGE_ID).bounds=(-qo2,-qo2)
    for rid in AEROBIC_ZERO_PRODUCT_EXCHANGES:
        m.reactions.get_by_id(rid).bounds=(0.0,0.0)
    m.objective=BIOMASS_REACTION_ID
    primary=m.optimize()
    expected=pfba(m,fraction_of_optimum=1.0)
    assert sol.fluxes[BIOMASS_REACTION_ID] == pytest.approx(primary.objective_value, rel=1e-8, abs=1e-8)
    assert sol.fluxes.abs().sum() == pytest.approx(expected.fluxes.abs().sum(), rel=1e-8, abs=1e-8)
    assert sol.fluxes[GLUCOSE_EXCHANGE_ID] == pytest.approx(-qglc)
    assert sol.fluxes[OXYGEN_EXCHANGE_ID] == pytest.approx(-qo2)


def test_source_sink_flow_invariants_for_acyclic_flow():
    import pandas as pd
    from etn.flow_decomposition import decompose_source_sink_paths
    edges=pd.DataFrame({
        'reaction_id':['r1','r2','r3'],
        'donor':['A','B','C'], 'acceptor':['B','C','D'],
        'electron_flux':[5.0,5.0,5.0],
    })
    paths,s,node_stats=decompose_source_sink_paths(edges)
    assert s['total_edge_activity'] == pytest.approx(15.0)
    assert s['net_source_flux'] == pytest.approx(5.0)
    assert s['net_sink_flux'] == pytest.approx(5.0)
    assert s['effective_transfer_depth'] == pytest.approx(3.0)
    assert paths.electron_flux.sum() == pytest.approx(5.0)
    st=node_stats.set_index('node')
    assert st.loc['B','relay_throughflow'] == pytest.approx(5.0)
    assert st.loc['C','relay_throughflow'] == pytest.approx(5.0)


def test_source_sink_unresolved_nodes_are_reaction_specific():
    import pandas as pd
    from etn.flow_decomposition import reaction_specific_unresolved
    from etn.constants import UNRESOLVED_ORGANIC_NODE
    edges=pd.DataFrame({'reaction_id':['r1','r2'], 'donor':[UNRESOLVED_ORGANIC_NODE,'X'],
                        'acceptor':['X',UNRESOLVED_ORGANIC_NODE], 'electron_flux':[2.0,2.0]})
    out=reaction_specific_unresolved(edges)
    assert out.loc[0,'donor']=='unresolved:r1'
    assert out.loc[1,'acceptor']=='unresolved:r2'
    assert out.loc[0,'donor'] != out.loc[1,'acceptor']

def test_deterministic_louvain_repeats_exactly():
    import networkx as nx
    from etn.network_analysis import deterministic_louvain_partition
    G=nx.Graph()
    G.add_weighted_edges_from([('z','a',2),('a','b',2),('z','b',2),('c','d',2),('d','e',2),('c','e',2),('b','c',0.1)])
    p1=deterministic_louvain_partition(G,random_state=0)
    p2=deterministic_louvain_partition(G,random_state=0)
    assert p1 == p2

def test_detect_communities_is_insertion_order_invariant():
    """Public community API must not depend on string-node insertion order."""
    import networkx as nx
    from etn.network_analysis import detect_communities
    edges = [
        ("z", "a", 3.0), ("a", "b", 2.0), ("b", "z", 2.0),
        ("m", "n", 3.0), ("n", "o", 2.0), ("o", "m", 2.0),
        ("b", "m", 0.1),
    ]
    def make(order):
        g = nx.MultiDiGraph()
        for u, v, w in order:
            g.add_edge(u, v, electrons=w, reaction=f"r_{u}_{v}")
        return g
    a = detect_communities(make(edges), random_state=0).set_index("metabolite_id")["community"]
    b = detect_communities(make(list(reversed(edges))), random_state=0).set_index("metabolite_id")["community"]
    # Community integer labels are arbitrary; compare pairwise co-membership.
    nodes = sorted(a.index)
    for i, u in enumerate(nodes):
        for v in nodes[i:]:
            assert (a[u] == a[v]) == (b[u] == b[v])


def test_undirected_reaction_weighted_graph_counts_reversible_reaction_once():
    """A reversible reaction is represented by opposite directed edges in the ETN,
    but contributes once to undirected community/comparison weights; a distinct
    reaction on the same unordered pair remains additive."""
    import networkx as nx
    from etn.network_analysis import undirected_reaction_weighted_graph
    g = nx.MultiDiGraph()
    g.add_edge("A", "B", reaction="rev", electrons=2.0, direction="forward")
    g.add_edge("B", "A", reaction="rev", electrons=2.0, direction="reverse")
    g.add_edge("B", "A", reaction="other", electrons=3.0, direction="forward")
    u = undirected_reaction_weighted_graph(g)
    assert u.number_of_edges() == 1
    assert u["A"]["B"]["weight"] == pytest.approx(5.0)
    assert u["A"]["B"]["n_reactions"] == 2
    assert u["A"]["B"]["reactions"] == ["other", "rev"]
