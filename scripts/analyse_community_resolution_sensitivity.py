"""Louvain resolution and enrichment-definition sensitivity diagnostics."""
import os
import sys
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")); sys.path.insert(0, REPO_ROOT)

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import adjusted_rand_score
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from etn.data_acquisition import load_ecoli_model
from etn.matching import classify_reactions_molecular
from etn.network_analysis import (
    build_electron_graph, undirected_reaction_weighted_graph, deterministic_louvain_partition,
    detect_communities, subsystem_enrichment,
)

RESULTS = os.path.join(REPO_ROOT, "results")
RESOLUTIONS = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]


def enrichment_from_subset(model, rows, communities, fully_resolved_only):
    rxn_by = {r.id: r for r in model.reactions}
    met_subsystems = {}
    for res in rows:
        if not res["is_redox"]:
            continue
        if fully_resolved_only and res.get("resolution_status") != "fully_resolved":
            continue
        subsystem = rxn_by[res["reaction_id"]].subsystem or "Unknown"
        for pair in res["pairs"]:
            for mid in (pair["reactant"], pair["product"]):
                met_subsystems.setdefault(mid, set()).add(subsystem)
    all_mets = set(communities.metabolite_id)
    all_subsystems = sorted({s for ss in met_subsystems.values() for s in ss})
    comm_of = dict(zip(communities.metabolite_id, communities.community))
    out=[]; N=len(all_mets)
    for comm in sorted(communities.community.unique()):
        cm={m for m,c in comm_of.items() if c==comm}; nc=len(cm)
        if nc<2: continue
        for sub in all_subsystems:
            sm={m for m in all_mets if sub in met_subsystems.get(m,set())}; ns=len(sm)
            if ns<2: continue
            ov=len(cm&sm)
            _,p=fisher_exact([[ov,nc-ov],[ns-ov,N-nc-ns+ov]],alternative="greater")
            out.append({"community":comm,"subsystem":sub,"community_size":nc,"subsystem_size":ns,"overlap":ov,"p_value":p})
    df=pd.DataFrame(out)
    if len(df): df["q_value"]=multipletests(df.p_value,method="fdr_bh")[1]
    return df


def main():
    model=load_ecoli_model(); rows=classify_reactions_molecular(model).to_dict("records")
    multi=build_electron_graph(model,rows); undirected=undirected_reaction_weighted_graph(multi)
    ref=deterministic_louvain_partition(undirected,random_state=0,resolution=1.0)
    resolution_rows=[]
    for resolution in RESOLUTIONS:
        p=deterministic_louvain_partition(undirected,random_state=0,resolution=resolution)
        common=sorted(ref)
        resolution_rows.append({
            "resolution":resolution,
            "n_communities":len(set(p.values())),
            "ARI_vs_resolution_1":adjusted_rand_score([ref[n] for n in common],[p[n] for n in common]),
        })
    pd.DataFrame(resolution_rows).to_csv(os.path.join(RESULTS,"community_resolution_sensitivity.csv"),index=False)

    communities=detect_communities(multi,random_state=0,resolution=1.0)
    primary=subsystem_enrichment(model,rows,communities)
    strict=enrichment_from_subset(model,rows,communities,True)
    summary=pd.DataFrame([
        {"annotation_definition":"all confidently redox matched-pair reactions (primary)","n_tests":len(primary),"n_q_lt_0.05":int((primary.q_value<0.05).sum())},
        {"annotation_definition":"fully resolved edge-generating reactions only (sensitivity)","n_tests":len(strict),"n_q_lt_0.05":int((strict.q_value<0.05).sum())},
    ])
    summary.to_csv(os.path.join(RESULTS,"community_enrichment_sensitivity_summary.csv"),index=False)
    strict.sort_values("p_value").to_csv(os.path.join(RESULTS,"community_subsystem_enrichment_fully_resolved_sensitivity.csv"),index=False)
    print(pd.DataFrame(resolution_rows).to_string(index=False))
    print("\n",summary.to_string(index=False))
    if len(strict): print("\nStrict significant tests:\n",strict[strict.q_value<0.05].sort_values('q_value').to_string(index=False))

if __name__ == "__main__": main()
