#!/usr/bin/env python3
"""Supplementary carbon-entry coverage screen beyond the six paired-state substrates.

Tests acetate (acetyl-CoA entry), pyruvate (pyruvate node) and succinate (TCA entry)
across the same five reconstructions under the same low- and high-O2 framework used
for the main carbon-entry screen.  The purpose is coverage, not to force a balanced
factorial: aerobic-only states are retained explicitly and paired feasibility is reported.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
OUT=ROOT/'results'/'publication'; OUT.mkdir(parents=True,exist_ok=True)

from scripts.run_cross_species_principle_analysis import (
    load_model, CONFIG, close_exchange, set_o2, substrate_metadata
)
from etn.cross_species import normalize_yeast_for_etn, apply_yeast_anaerobic, chemical_model, etn_metrics, solve_growth_pfba
from etn.matching import build_metabolite_table, build_matching_context, clear_matching_context_cache

NEW_CARBON={
 'acetate': {'C':2.0,'gamma':8.0,'entry':'acetyl-CoA / two-carbon oxidation'},
 'pyruvate': {'C':3.0,'gamma':10.0,'entry':'pyruvate node'},
 'succinate': {'C':4.0,'gamma':14.0,'entry':'TCA-cycle intermediate'},
}
NEW_EX={
 'E. coli': {'acetate':'EX_ac_e','pyruvate':'EX_pyr_e','succinate':'EX_succ_e'},
 'B. subtilis': {'acetate':'EX_ac_e','pyruvate':'EX_pyr_e','succinate':'EX_succ_e'},
 'S. enterica': {'acetate':'EX_ac_e','pyruvate':'EX_pyr_e','succinate':'EX_succ_e'},
 'K. phaffii': {'acetate':'Ex_ac','pyruvate':'Ex_pyr','succinate':'Ex_succ'},
 'S. cerevisiae': {'acetate':'r_1634','pyruvate':'r_2033','succinate':'r_2056'},
}

# All nine explicitly tested carbon routes are closed before opening the selected one.
def all_tested_exchanges(org):
    return [x for x in CONFIG[org]['subs'].values() if x] + list(NEW_EX[org].values())


def prepare(base,org,sub,state):
    cfg=CONFIG[org]; ex=NEW_EX[org][sub]
    if ex not in base.reactions: raise KeyError(f'{ex} absent')
    if org=='S. cerevisiae' and state=='low-respiratory-capacity':
        m=base.copy(); apply_yeast_anaerobic(m, ROOT/'data'/'cross_species'/'aminoAcid_Bjorkeroth2020.tsv')
    else: m=base.copy()
    m.objective=cfg['objective']; pichia=(org=='K. phaffii')
    for sx in all_tested_exchanges(org):
        if sx in m.reactions: close_exchange(m,sx,pichia)
    for ax in cfg['close_acceptors']:
        close_exchange(m,ax,pichia)
    if org=='B. subtilis' and 'EX_nh4_e' in m.reactions:
        m.reactions.get_by_id('EX_nh4_e').bounds=(-10.0,1000.0)
    m.reactions.get_by_id(ex).bounds=(-10.0,-10.0)
    cap=20.0 if state=='high-respiratory-capacity' else cfg['low_o2']
    if org=='S. cerevisiae':
        o2=m.reactions.get_by_id(cfg['o2'])
        o2.bounds=(-cap,max(0.0,float(o2.upper_bound))) if cap else (0.0,max(0.0,float(o2.upper_bound)))
    else:
        set_o2(m,cfg['o2'],cap,pichia)
    return m,cap


def run_org(org):
    base=load_model(org); cfg=CONFIG[org]
    clear_matching_context_cache()
    chem_src=normalize_yeast_for_etn(base) if org=='S. cerevisiae' else base
    chem=chemical_model(chem_src,yeast=(org=='S. cerevisiae'))
    mt=build_metabolite_table(chem); build_matching_context(chem,mt)
    rows=[]; feas=[]
    for sub,meta in NEW_CARBON.items():
        ex=NEW_EX[org][sub]
        if ex not in base.reactions:
            for state in ('low-respiratory-capacity','high-respiratory-capacity'):
                feas.append({'organism':org,'substrate':sub,'entry_point':meta['entry'],'state':state,'status':'unavailable','growth':np.nan,'reason':'exchange absent'})
            continue
        for state in ('low-respiratory-capacity','high-respiratory-capacity'):
            try:
                m,cap=prepare(base,org,sub,state)
                sol=solve_growth_pfba(m,cfg['objective'])
                growth=float(sol.fluxes[cfg['objective']])
                metrics,_,_=etn_metrics(chem,sol,10.0)
                actual_o2=max(0.0,-float(sol.fluxes[cfg['o2']]))
                co2id=next((x for x in ['EX_co2_e','Ex_co2','r_1672'] if x in m.reactions),None)
                co2=float(sol.fluxes[co2id]) if co2id else np.nan
                try: formula,charge,C,gamma=substrate_metadata(base,ex)
                except Exception: formula,charge,C,gamma='',np.nan,np.nan,np.nan
                if not np.isfinite(C): C=meta['C']
                if not np.isfinite(gamma): gamma=meta['gamma']
                rows.append({'organism':org,'substrate':sub,'entry_point':meta['entry'],'exchange':ex,'state':state,'oxygen_cap':cap,'actual_o2':actual_o2,'growth':growth,'co2':co2,'formula':formula,'charge':charge,'C':C,'gamma':gamma,'net_e_per_C':metrics['net_per_glucose']/C,'fraction_substrate_gamma':metrics['net_per_glucose']/gamma,**metrics})
                feas.append({'organism':org,'substrate':sub,'entry_point':meta['entry'],'state':state,'status':'feasible','growth':growth,'reason':''})
            except Exception as exc:
                feas.append({'organism':org,'substrate':sub,'entry_point':meta['entry'],'state':state,'status':'infeasible_or_failed','growth':np.nan,'reason':str(exc)[:300]})
    tag=org.replace(' ','_').replace('.','')
    pd.DataFrame(rows).to_csv(OUT/f'extended_carbon_entry_states_{tag}.csv',index=False)
    pd.DataFrame(feas).to_csv(OUT/f'extended_carbon_entry_feasibility_{tag}.csv',index=False)
    print(org,'states',len(rows),flush=True)


def merge():
    states=[]; feas=[]
    for org in CONFIG:
        tag=org.replace(' ','_').replace('.','')
        ps=OUT/f'extended_carbon_entry_states_{tag}.csv'; pf=OUT/f'extended_carbon_entry_feasibility_{tag}.csv'
        if ps.exists():
            try:
                d=pd.read_csv(ps)
                if not d.empty: states.append(d)
            except pd.errors.EmptyDataError: pass
        if pf.exists(): feas.append(pd.read_csv(pf))
    st=pd.concat(states,ignore_index=True) if states else pd.DataFrame()
    ff=pd.concat(feas,ignore_index=True) if feas else pd.DataFrame()
    st.to_csv(OUT/'extended_carbon_entry_states.csv',index=False); ff.to_csv(OUT/'extended_carbon_entry_feasibility.csv',index=False)
    pairs=[]
    if not st.empty:
        for (org,sub),g in st.groupby(['organism','substrate']):
            lo=g[g.state=='low-respiratory-capacity']; hi=g[g.state=='high-respiratory-capacity']
            row={'organism':org,'substrate':sub,'entry_point':g.entry_point.iloc[0],'low_feasible':bool(len(lo)),'high_feasible':bool(len(hi))}
            if len(lo):
                row.update(low_growth=float(lo.growth.iloc[0]),low_net_e_per_substrate=float(lo.net_per_glucose.iloc[0]))
            if len(hi):
                row.update(high_growth=float(hi.growth.iloc[0]),high_net_e_per_substrate=float(hi.net_per_glucose.iloc[0]))
            if len(lo) and len(hi):
                row['electron_flux_expansion']=float(hi.net_per_glucose.iloc[0]-lo.net_per_glucose.iloc[0])
                row['growth_ratio_low_to_high']=float(lo.growth.iloc[0]/hi.growth.iloc[0]) if hi.growth.iloc[0] else np.nan
            pairs.append(row)
    pp=pd.DataFrame(pairs); pp.to_csv(OUT/'extended_carbon_entry_summary.csv',index=False)
    stats={
      'candidate_organism_substrate_pairs':15,
      'high_disposal_feasible_pairs':int((pp.high_feasible==True).sum()) if len(pp) else 0,
      'paired_low_high_feasible_pairs':int(((pp.low_feasible==True)&(pp.high_feasible==True)).sum()) if len(pp) else 0,
      'high_feasible_by_substrate':{s:int(((pp.substrate==s)&(pp.high_feasible==True)).sum()) for s in NEW_CARBON} if len(pp) else {},
      'paired_feasible_by_substrate':{s:int(((pp.substrate==s)&(pp.low_feasible==True)&(pp.high_feasible==True)).sum()) for s in NEW_CARBON} if len(pp) else {},
    }
    (OUT/'extended_carbon_entry_statistics.json').write_text(json.dumps(stats,indent=2)+'\n')
    print(json.dumps(stats,indent=2))

if __name__=='__main__':
    if '--merge' in sys.argv: merge()
    elif '--organism' in sys.argv: run_org(sys.argv[sys.argv.index('--organism')+1])
    else:
        for org in CONFIG: run_org(org)
        merge()
