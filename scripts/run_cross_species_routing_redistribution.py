#!/usr/bin/env python3
"""Quantify whole-network electron-routing redistribution across paired electron-disposal-capacity states.

For each of the five publication cross-species glucose comparisons, compute the
Jensen-Shannon distance and total-variation distance between normalized
reaction-resolved donor->acceptor electron-flux distributions in the
high-respiratory-capacity and low-respiratory-capacity states. This complements Fig. 4b without
changing its main-panel metrics.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
DATA = ROOT / "data" / "cross_species"
OUT = ROOT / "results" / "publication"
OUT.mkdir(parents=True, exist_ok=True)

from etn.data_acquisition import load_ecoli_model
from etn.fba_case_study import solve_condition
from etn.cross_species import (
    apply_yeast_anaerobic, chemical_model, etn_metrics,
    load_cross_species_models, normalize_yeast_for_etn, solve_growth_pfba,
)
from etn.matching import clear_matching_context_cache


def edge_distribution(edges: pd.DataFrame) -> pd.Series:
    # Reaction ID + donor + acceptor identifies a resolved electron-transfer edge.
    e = edges.copy()
    e = e[np.isfinite(e["electron_flux"]) & (e["electron_flux"].abs() > 1e-12)].copy()
    e["key"] = e["reaction_id"].astype(str) + "|" + e["donor"].astype(str) + "|" + e["acceptor"].astype(str)
    s = e.groupby("key", sort=True)["electron_flux"].apply(lambda x: float(np.abs(x).sum()))
    total = float(s.sum())
    return s / total if total > 0 else s


def compare(org: str, hi_edges: pd.DataFrame, lo_edges: pd.DataFrame, depth_hi: float, depth_lo: float):
    p = edge_distribution(hi_edges)
    q = edge_distribution(lo_edges)
    keys = p.index.union(q.index)
    pv = p.reindex(keys, fill_value=0.0).to_numpy(float)
    qv = q.reindex(keys, fill_value=0.0).to_numpy(float)
    jsd = float(jensenshannon(pv, qv, base=2.0))
    tvd = float(0.5 * np.abs(pv - qv).sum())
    changed = int(np.sum(np.abs(pv-qv) > 1e-4))
    return {
        "organism": org,
        "jensen_shannon_distance_bits": jsd,
        "total_variation_distance": tvd,
        "n_union_edges": int(len(keys)),
        "n_edges_changed_gt_1e-4_fraction": changed,
        "transfer_depth_change_pct": float(100 * (depth_lo / depth_hi - 1)),
    }


def metrics_edges(chem, sol, q):
    m,e,_ = etn_metrics(chem, sol, q)
    return m,e


def main():
    clear_matching_context_cache()
    ecoli = load_ecoli_model()
    others = load_cross_species_models(DATA)
    yeast, bac, salm, pichia = (others[k] for k in ["S. cerevisiae","B. subtilis","S. enterica","K. phaffii"])
    chems = {
        "E. coli": chemical_model(ecoli),
        "S. cerevisiae": chemical_model(normalize_yeast_for_etn(yeast), yeast=True),
        "B. subtilis": chemical_model(bac),
        "S. enterica": chemical_model(salm),
        "K. phaffii": chemical_model(pichia),
    }
    pairs = {}

    # E. coli: manuscript physiological states.
    hi = solve_condition(ecoli, 8.7, 11.9, True); lo = solve_condition(ecoli, 14.9, 0.0, False)
    pairs["E. coli"] = (*metrics_edges(chems["E. coli"], hi, 8.7), *metrics_edges(chems["E. coli"], lo, 14.9))

    # S. cerevisiae.
    ya = yeast.copy(); ya.objective="r_2111"; ya.reactions.get_by_id("r_1714").bounds=(-10,-10); ya.reactions.get_by_id("r_1992").lower_bound=-1000
    hi = solve_growth_pfba(ya,"r_2111")
    yn = yeast.copy(); apply_yeast_anaerobic(yn, DATA/"aminoAcid_Bjorkeroth2020.tsv"); yn.objective="r_2111"; yn.reactions.get_by_id("r_1714").bounds=(-10,-10)
    lo = solve_growth_pfba(yn,"r_2111")
    pairs["S. cerevisiae"] = (*metrics_edges(chems["S. cerevisiae"], hi, 10), *metrics_edges(chems["S. cerevisiae"], lo, 10))

    # B. subtilis.
    bs=[]
    for o2 in [(-1000,1000),(0,1000)]:
        mm=bac.copy(); mm.objective="BiomassRepsolRed"; mm.reactions.get_by_id("EX_glc__D_e").bounds=(-10,-10); mm.reactions.get_by_id("EX_nh4_e").bounds=(-10,1000); mm.reactions.get_by_id("EX_no3_e").bounds=(0,1000); mm.reactions.get_by_id("EX_o2_e").bounds=o2
        bs.append(solve_growth_pfba(mm,"BiomassRepsolRed"))
    pairs["B. subtilis"] = (*metrics_edges(chems["B. subtilis"], bs[0],10), *metrics_edges(chems["B. subtilis"],bs[1],10))

    # S. enterica.
    ss=[]
    for cap in [20,0]:
        mm=salm.copy(); mm.objective="BIOMASS_iRR1083_1"; mm.reactions.get_by_id("EX_glc__D_e").bounds=(-10,-10); mm.reactions.get_by_id("EX_no3_e").bounds=(0,1000); mm.reactions.get_by_id("EX_o2_e").bounds=(-cap,1000) if cap else (0,1000)
        ss.append(solve_growth_pfba(mm,"BIOMASS_iRR1083_1"))
    pairs["S. enterica"] = (*metrics_edges(chems["S. enterica"],ss[0],10), *metrics_edges(chems["S. enterica"],ss[1],10))

    # K. phaffii: high-respiratory-capacity vs publication low-oxygen endpoint (cap 0.5).
    ps=[]
    for cap in [20,0.5]:
        mm=pichia.copy(); mm.objective="Ex_biomass"; mm.reactions.get_by_id("Ex_glc_D").bounds=(-10,-10); mm.reactions.get_by_id("Ex_o2").bounds=(-cap,0)
        ps.append(solve_growth_pfba(mm,"Ex_biomass"))
    pairs["K. phaffii"] = (*metrics_edges(chems["K. phaffii"],ps[0],10), *metrics_edges(chems["K. phaffii"],ps[1],10))

    rows=[]
    detail=[]
    for org,(mhi,ehi,mlo,elo) in pairs.items():
        rows.append(compare(org,ehi,elo,mhi["effective_transfer_depth"],mlo["effective_transfer_depth"]))
        for state,edges in [("high_respiratory_capacity",ehi),("low_respiratory_capacity",elo)]:
            d=edge_distribution(edges)
            for key,val in d.items():
                rid,donor,acceptor=key.split("|",2)
                detail.append({"organism":org,"state":state,"reaction_id":rid,"donor":donor,"acceptor":acceptor,"normalized_edge_electron_flux":float(val)})
    df=pd.DataFrame(rows).sort_values("jensen_shannon_distance_bits",ascending=False)
    df.to_csv(OUT/"cross_species_routing_redistribution.csv",index=False)
    pd.DataFrame(detail).to_csv(OUT/"cross_species_routing_edge_distributions.csv",index=False)
    summary={
        "metric_definition":"Jensen-Shannon distance (base 2) and total variation distance between normalized reaction-resolved donor-acceptor electron-flux distributions in paired high-respiratory-capacity and low-respiratory-capacity states.",
        "n_organisms":int(len(df)),
        "jsd_min":float(df.jensen_shannon_distance_bits.min()),
        "jsd_max":float(df.jensen_shannon_distance_bits.max()),
        "tvd_min":float(df.total_variation_distance.min()),
        "tvd_max":float(df.total_variation_distance.max()),
    }
    (OUT/"cross_species_routing_redistribution_statistics.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(df.to_string(index=False))
    print(json.dumps(summary,indent=2))

if __name__ == "__main__": main()
