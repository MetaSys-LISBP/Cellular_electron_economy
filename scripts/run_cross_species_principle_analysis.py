#!/usr/bin/env python3
"""Test the carbon-entry/electron-disposal organizing principle across species.

This publication analysis applies paired high- and low-electron-disposal-capacity states to
all six manuscript carbon sources across E. coli, B. subtilis, S. enterica, K. phaffii and S. cerevisiae. Every organism–substrate combination is attempted; combinations lacking a canonical substrate exchange or positive-growth state are retained explicitly in the feasibility table.
"""
from __future__ import annotations
import json, sys, subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUT = ROOT / "results" / "publication"
OUT.mkdir(parents=True, exist_ok=True)

from etn.data_acquisition import load_ecoli_model
from etn.cross_species import load_bacillus_xlsx, load_standard_gem_xlsx, load_cross_species_models, normalize_yeast_for_etn, apply_yeast_anaerobic, chemical_model, etn_metrics, solve_growth_pfba
from cobra.io import read_sbml_model
from etn.degree_of_reduction import degree_of_reduction
from etn.matching import build_metabolite_table, build_matching_context, clear_matching_context_cache

SUBSTRATE_CHEMISTRY = {
    "glucose": (6.0, 24.0),
    "fructose": (6.0, 24.0),
    "galactose": (6.0, 24.0),
    "xylose": (5.0, 20.0),
    "gluconate": (6.0, 22.0),
    "glycerol": (3.0, 14.0),
}

CONFIG = {
    "E. coli": {
        "objective": "BIOMASS_Ec_iML1515_core_75p37M", "o2": "EX_o2_e", "low_o2": 0.0,
        "subs": {"glucose":"EX_glc__D_e","fructose":"EX_fru_e","xylose":"EX_xyl__D_e","galactose":"EX_gal_e","gluconate":"EX_glcn_e","glycerol":"EX_glyc_e"},
        "close_acceptors": ["EX_no3_e","EX_tmao_e","EX_dmso_e","EX_fum_e"],
    },
    "B. subtilis": {
        "objective": "BiomassRepsolRed", "o2": "EX_o2_e", "low_o2": 0.0,
        "subs": {"glucose":"EX_glc__D_e","fructose":"EX_fru_e","xylose":"EX_xyl__D_e","galactose":"EX_gal_e","gluconate":"EX_glcn__D_e","glycerol":"EX_glyc_e"},
        "close_acceptors": ["EX_no3_e","EX_fum_e"],
    },
    "S. enterica": {
        "objective": "BIOMASS_iRR1083_1", "o2": "EX_o2_e", "low_o2": 0.0,
        "subs": {"glucose":"EX_glc__D_e","fructose":"EX_fru_e","xylose":"EX_xyl__D_e","galactose":"EX_gal_e","gluconate":"EX_glcn_e","glycerol":"EX_glyc_e"},
        "close_acceptors": ["EX_no3_e","EX_tmao_e","EX_dmso_e","EX_fum_e"],
    },
    "K. phaffii": {
        "objective": "Ex_biomass", "o2": "Ex_o2", "low_o2": 0.5,
        "subs": {"glucose":"Ex_glc_D","fructose":"Ex_fru","galactose":"Ex_gal","xylose":"Ex_xyl_D","gluconate":None,"glycerol":"Ex_glyc"},
        "close_acceptors": ["Ex_no3","Ex_fum"],
    },
    "S. cerevisiae": {
        "objective": "r_2111", "o2": "r_1992", "low_o2": 0.0,
        "subs": {"glucose":"r_1714","fructose":"r_1709","galactose":"r_1710","xylose":"r_1718","gluconate":None,"glycerol":"r_1808"},
        "close_acceptors": [],
    },
}

def load_model(org):
    data = ROOT / "data" / "cross_species"
    if org == "E. coli":
        return load_ecoli_model()
    if org == "B. subtilis":
        return load_bacillus_xlsx(data / "iBB1018.xlsx")
    if org == "S. enterica":
        return read_sbml_model(str(data / "iYS1720.xml.gz"))
    if org == "K. phaffii":
        return load_standard_gem_xlsx(data / "iMT1026_model.xlsx", "iMT1026")
    if org == "S. cerevisiae":
        return load_standard_gem_xlsx(data / "yeast-GEM-v9.1.0.xlsx", "yeastGEM_v9.1.0")
    raise KeyError(org)

def close_exchange(model, rid, pichia=False):
    if rid not in model.reactions: return
    r = model.reactions.get_by_id(rid)
    r.bounds = (0.0, 0.0) if pichia else (0.0, max(0.0, float(r.upper_bound)))

def set_uptake(model, rid, value=10.0, pichia=False):
    r = model.reactions.get_by_id(rid)
    r.bounds = (-float(value), 0.0 if pichia else max(1000.0, float(r.upper_bound)))

def set_o2(model, rid, cap, pichia=False):
    r = model.reactions.get_by_id(rid)
    r.bounds = ((-float(cap), 0.0) if cap else (0.0, 0.0)) if pichia else ((-float(cap), max(1000.0,float(r.upper_bound))) if cap else (0.0, max(1000.0,float(r.upper_bound))))

def substrate_metadata(model, ex):
    rxn = model.reactions.get_by_id(ex)
    mets = list(rxn.metabolites)
    met = next((m for m in mets if m.compartment in {"e","extracellular"}), mets[0])
    try:
        gamma = float(degree_of_reduction(met.elements, met.charge))
    except Exception:
        gamma = np.nan
    try: carbon = float(met.elements.get("C", np.nan))
    except Exception: carbon = np.nan
    return met.formula, met.charge, carbon, gamma

def prepare_state(base, org, sub, state):
    cfg = CONFIG[org]
    ex = cfg["subs"].get(sub)
    if not ex:
        raise KeyError(f"canonical {sub} exchange unavailable in reconstruction")
    if ex not in base.reactions:
        raise KeyError(f"configured exchange {ex} absent from reconstruction")

    if org == "S. cerevisiae" and state == "low-disposal-capacity":
        m = base.copy()
        apply_yeast_anaerobic(m, ROOT / "data" / "cross_species" / "aminoAcid_Bjorkeroth2020.tsv")
    else:
        m = base.copy()
    m.objective = cfg["objective"]

    pichia = org == "K. phaffii"
    # Close every tested carbon-source exchange first so that the selected
    # carbon-entry route is unambiguous. Missing exchanges are simply absent.
    for sx in cfg["subs"].values():
        if sx:
            close_exchange(m, sx, pichia)
    for ax in cfg["close_acceptors"]:
        close_exchange(m, ax, pichia)

    if org == "B. subtilis" and "EX_nh4_e" in m.reactions:
        m.reactions.get_by_id("EX_nh4_e").bounds = (-10.0, 1000.0)

    # Fixed uptake makes the per-substrate comparison exact rather than relying
    # on the growth optimum to saturate a one-sided uptake bound.
    r = m.reactions.get_by_id(ex)
    r.bounds = (-10.0, -10.0)

    cap = 20.0 if state == "high-disposal-capacity" else cfg["low_o2"]
    if org == "S. cerevisiae":
        o2 = m.reactions.get_by_id(cfg["o2"])
        o2.bounds = (-cap, max(0.0, float(o2.upper_bound))) if cap else (0.0, max(0.0, float(o2.upper_bound)))
    else:
        set_o2(m, cfg["o2"], cap, pichia)
    return m, cap


def run_one(org, only_substrate=None):
    base = load_model(org); cfg = CONFIG[org]; rows = []; feasibility = []
    clear_matching_context_cache()
    chem_source = normalize_yeast_for_etn(base) if org == "S. cerevisiae" else base
    chem = chemical_model(chem_source, yeast=(org == "S. cerevisiae"))
    mt = build_metabolite_table(chem); build_matching_context(chem, mt)

    substrates = list(SUBSTRATE_CHEMISTRY)
    if only_substrate is not None:
        substrates = [only_substrate]
    for sub in substrates:
        ex = cfg["subs"].get(sub)
        if not ex or ex not in base.reactions:
            for state in ("low-disposal-capacity", "high-disposal-capacity"):
                feasibility.append({"organism":org,"substrate":sub,"state":state,"status":"unavailable","growth":np.nan,"reason":"canonical substrate exchange absent from reconstruction"})
            continue
        for state in ("low-disposal-capacity", "high-disposal-capacity"):
            try:
                m, cap = prepare_state(base, org, sub, state)
                sol = solve_growth_pfba(m, cfg["objective"])
                growth = float(sol.fluxes[cfg["objective"]])
                metrics, _, _ = etn_metrics(chem, sol, 10.0)
                actual_o2 = max(0.0, -float(sol.fluxes[cfg["o2"]]))
                co2_candidates = ["EX_co2_e", "Ex_co2", "r_1672"]
                co2_id = next((x for x in co2_candidates if x in m.reactions), None)
                co2 = float(sol.fluxes[co2_id]) if co2_id else np.nan
                formula, charge, carbon, gamma = substrate_metadata(base, ex)
                # Canonical chemistry values are used if model metabolite formula
                # metadata are incomplete or inconsistent across reconstructions.
                carbon_ref, gamma_ref = SUBSTRATE_CHEMISTRY[sub]
                if not np.isfinite(carbon): carbon = carbon_ref
                if not np.isfinite(gamma): gamma = gamma_ref
                row = {"organism":org,"substrate":sub,"exchange":ex,"state":state,"oxygen_cap":cap,"actual_o2":actual_o2,"growth":growth,"co2":co2,"formula":formula,"charge":charge,"C":carbon,"gamma":gamma,**metrics}
                row["net_per_C"] = metrics["net_per_glucose"] / carbon
                rows.append(row)
                feasibility.append({"organism":org,"substrate":sub,"state":state,"status":"feasible","growth":growth,"reason":""})
            except Exception as exc:
                feasibility.append({"organism":org,"substrate":sub,"state":state,"status":"infeasible_or_failed","growth":np.nan,"reason":str(exc)[:300]})
    tag = org.replace(" ","_").replace(".","") + ("_"+only_substrate if only_substrate else "")
    pd.DataFrame(rows).to_csv(OUT/f"cross_species_carbon_acceptor_states_{tag}.csv", index=False)
    pd.DataFrame(feasibility).to_csv(OUT/f"cross_species_carbon_acceptor_feasibility_{tag}.csv", index=False)
    print(org, "complete", len(rows), flush=True)

def merge():
    state_files=[x for x in OUT.glob("cross_species_carbon_acceptor_states_*.csv") if x.name != "cross_species_carbon_acceptor_states.csv"]
    feas_files=[x for x in OUT.glob("cross_species_carbon_acceptor_feasibility_*.csv") if x.name != "cross_species_carbon_acceptor_feasibility.csv"]
    state_frames=[]
    for x in state_files:
        try:
            df=pd.read_csv(x)
        except pd.errors.EmptyDataError:
            continue
        if not df.empty:
            state_frames.append(df)
    if not state_frames:
        raise RuntimeError("No cross-species carbon/acceptor state files were produced.")
    raw=pd.concat(state_frames,ignore_index=True)
    # Deduplicate by biological state in case a single-organism and per-substrate run coexist.
    raw=raw.drop_duplicates(["organism","substrate","state"],keep="last").sort_values(["organism","substrate","state"])
    raw.to_csv(OUT/"cross_species_carbon_acceptor_states.csv",index=False)
    feas=pd.concat([pd.read_csv(x) for x in feas_files],ignore_index=True)
    feas=feas.drop_duplicates(["organism","substrate","state"],keep="last").sort_values(["organism","substrate","state"])
    feas.to_csv(OUT/"cross_species_carbon_acceptor_feasibility.csv",index=False)
    pairs=[]
    for (org,sub),g in raw.groupby(["organism","substrate"]):
        if set(g.state)!={"low-disposal-capacity","high-disposal-capacity"}: continue
        lo=g[g.state=="low-disposal-capacity"].iloc[0]; hi=g[g.state=="high-disposal-capacity"].iloc[0]
        pairs.append({"organism":org,"substrate":sub,"low_net_per_substrate":lo.net_per_glucose,"high_net_per_substrate":hi.net_per_glucose,"electron_flux_expansion":hi.net_per_glucose-lo.net_per_glucose,"low_co2":lo.co2/10.0,"high_co2":hi.co2/10.0,"additional_carbon_oxidation":(hi.co2-lo.co2)/10.0,"low_growth":lo.growth,"high_growth":hi.growth,"low_transfer_depth":lo.effective_transfer_depth,"high_transfer_depth":hi.effective_transfer_depth,"low_actual_o2":lo.actual_o2,"high_actual_o2":hi.actual_o2,"C":lo.C,"gamma":lo.gamma})
    non_ecoli=pd.DataFrame(pairs)
    non_ecoli.to_csv(OUT/"cross_species_carbon_acceptor_pairs_non_ecoli.csv",index=False)

    # The publication-facing comparison is now the systematic five-organism screen itself.
    combined = non_ecoli.copy().sort_values(["organism", "substrate"])
    combined["low_e_per_C"] = combined["low_net_per_substrate"] / combined["C"]
    combined["high_e_per_C"] = combined["high_net_per_substrate"] / combined["C"]
    combined["expansion_e_per_C"] = combined["electron_flux_expansion"] / combined["C"]
    combined["low_fraction_substrate_gamma"] = combined["low_net_per_substrate"] / combined["gamma"]
    combined["high_fraction_substrate_gamma"] = combined["high_net_per_substrate"] / combined["gamma"]
    combined["expansion_fraction_substrate_gamma"] = combined["electron_flux_expansion"] / combined["gamma"]
    combined.to_csv(OUT/"cross_species_carbon_acceptor_pairs.csv", index=False)

    stats={
        "n_pairs":int(len(combined)),
        "n_candidate_pairs":30,
        "n_fully_feasible_pairs":int(len(combined)),
        "organisms":sorted(combined.organism.unique().tolist()),
        "all_pairs_expand":bool((combined.electron_flux_expansion>1e-8).all()),
        "low_baseline_min":float(combined.low_net_per_substrate.min()),
        "low_baseline_max":float(combined.low_net_per_substrate.max()),
        "expansion_min":float(combined.electron_flux_expansion.min()),
        "expansion_max":float(combined.electron_flux_expansion.max()),
    }
    per={}
    norm_summary=[]
    for org,g in combined.groupby("organism"):
        per[org]={
            "n":int(len(g)),
            "baseline_min":float(g.low_net_per_substrate.min()),
            "baseline_max":float(g.low_net_per_substrate.max()),
            "baseline_range":float(g.low_net_per_substrate.max()-g.low_net_per_substrate.min()),
            "baseline_per_C_min":float(g.low_e_per_C.min()),
            "baseline_per_C_max":float(g.low_e_per_C.max()),
            "baseline_per_C_range_factor":float(g.low_e_per_C.max()/g.low_e_per_C.min()),
            "baseline_fraction_gamma_min":float(g.low_fraction_substrate_gamma.min()),
            "baseline_fraction_gamma_max":float(g.low_fraction_substrate_gamma.max()),
            "baseline_fraction_gamma_range_factor":float(g.low_fraction_substrate_gamma.max()/g.low_fraction_substrate_gamma.min()),
            "all_expand":bool((g.electron_flux_expansion>1e-8).all()),
            "all_expand_per_C":bool((g.expansion_e_per_C>1e-8).all()),
            "all_expand_fraction_gamma":bool((g.expansion_fraction_substrate_gamma>1e-8).all()),
            "expansion_min":float(g.electron_flux_expansion.min()),
            "expansion_max":float(g.electron_flux_expansion.max())
        }
        norm_summary.append({"organism":org, **per[org]})
    stats["per_organism"]=per
    pd.DataFrame(norm_summary).to_csv(OUT/"cross_species_carbon_acceptor_normalization_summary.csv",index=False)
    # CO2 relationships are reported within individual models only; pooled values across
    # organisms are not biologically comparable because biomass and carbon-output terms differ.
    co2_per={}
    for org,g in non_ecoli.groupby("organism"):
        if len(g)>=3 and np.std(g.additional_carbon_oxidation)>1e-10 and np.std(g.electron_flux_expansion)>1e-10:
            r,_=pearsonr(g.additional_carbon_oxidation,g.electron_flux_expansion)
            rho,_=spearmanr(g.additional_carbon_oxidation,g.electron_flux_expansion)
            co2_per[org]={"n":int(len(g)),"pearson_r":float(r),"spearman_rho":float(rho)}
    stats["within_model_expansion_vs_CO2"]=co2_per
    (OUT/"cross_species_carbon_acceptor_statistics.json").write_text(json.dumps(stats,indent=2)+"\n")
    print(json.dumps(stats,indent=2))

def _cleanup_intermediate_pair_files():
    """Remove per-pair scratch CSVs after the canonical merged outputs exist."""
    for pattern in ("cross_species_carbon_acceptor_states_*.csv", "cross_species_carbon_acceptor_feasibility_*.csv"):
        for path in OUT.glob(pattern):
            if path.name in {"cross_species_carbon_acceptor_states.csv", "cross_species_carbon_acceptor_feasibility.csv"}:
                continue
            path.unlink(missing_ok=True)


def main():
    if "--merge" in sys.argv:
        merge(); return
    if "--organism" in sys.argv:
        org=sys.argv[sys.argv.index("--organism")+1]
        sub=sys.argv[sys.argv.index("--substrate")+1] if "--substrate" in sys.argv else None
        run_one(org, sub); return

    # Run one fresh interpreter per organism. This reuses the static chemistry
    # reconstruction across six substrates but prevents cross-organism solver/model
    # state from accumulating.
    script = Path(__file__).resolve()
    for org in ("E. coli", "B. subtilis", "S. enterica", "K. phaffii", "S. cerevisiae"):
        subprocess.run([sys.executable, str(script), "--organism", org], cwd=ROOT, check=True)
    merge()
    _cleanup_intermediate_pair_files()

if __name__=="__main__": main()
