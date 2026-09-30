#!/usr/bin/env python3
"""Extended capability-based terminal electron-acceptor screen.

The screen includes O2, nitrate, TMAO, DMSO, fumarate, nitrite, Salmonella
tetrathionate and thiosulfate.  All candidate exchanges are closed in the
baseline.  A candidate is called functional only if the model is feasible *and*
actual uptake is > tolerance.  Carbon-bearing fumarate is retained but excluded
from clean CO2-based mechanistic correlations.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
OUT=ROOT/'results'/'publication'; OUT.mkdir(parents=True,exist_ok=True)
DATA=ROOT/'data'/'cross_species'

from scripts.run_systematic_acceptor_screen import load_model, close_exchange, open_uptake
from etn.cross_species import normalize_yeast_for_etn, chemical_model, etn_metrics, solve_growth_pfba, apply_yeast_anaerobic
from etn.matching import build_metabolite_table, build_matching_context, clear_matching_context_cache

# e_per_acceptor is the stoichiometric electron demand for the represented reduction.
CONFIG={
 'E. coli': {'objective':'BIOMASS_Ec_iML1515_core_75p37M','glc':'EX_glc__D_e','o2':'EX_o2_e',
   'acceptors':{
    'oxygen':('EX_o2_e',4.0,False),'nitrate':('EX_no3_e',2.0,False),'nitrite':('EX_no2_e',6.0,False),
    'TMAO':('EX_tmao_e',2.0,False),'DMSO':('EX_dmso_e',2.0,False),'fumarate':('EX_fum_e',2.0,True),
    'thiosulfate':('EX_tsul_e',2.0,False)}},
 'B. subtilis': {'objective':'BiomassRepsolRed','glc':'EX_glc__D_e','o2':'EX_o2_e',
   'acceptors':{'oxygen':('EX_o2_e',4.0,False),'nitrate':('EX_no3_e',2.0,False),'nitrite':('EX_no2_e',6.0,False),'fumarate':('EX_fum_e',2.0,True)}},
 'S. enterica': {'objective':'BIOMASS_iRR1083_1','glc':'EX_glc__D_e','o2':'EX_o2_e',
   'acceptors':{
    'oxygen':('EX_o2_e',4.0,False),'nitrate':('EX_no3_e',2.0,False),'nitrite':('EX_no2_e',6.0,False),
    'TMAO':('EX_tmao_e',2.0,False),'DMSO':('EX_dmso_e',2.0,False),'fumarate':('EX_fum_e',2.0,True),
    'tetrathionate':('EX_tet_e',2.0,False),'thiosulfate':('EX_tsul_e',2.0,False)}},
 'K. phaffii': {'objective':'Ex_biomass','glc':'Ex_glc_D','o2':'Ex_o2',
   'acceptors':{'oxygen':('Ex_o2',4.0,False),'nitrate':('Ex_no3',2.0,False),'fumarate':('Ex_fum',2.0,True)}},
 'S. cerevisiae': {'objective':'r_2111','glc':'r_1714','o2':'r_1992',
   'acceptors':{'oxygen':('r_1992',4.0,False),'fumarate':('r_1798',2.0,True)}},
}
CAP=20.0; TOL=1e-8


def prepare(base,org,selected=None):
    cfg=CONFIG[org]
    # Use model-provided anaerobic physiology for yeast in non-O2/baseline states.
    if org=='S. cerevisiae' and selected!='oxygen':
        m=base.copy(); apply_yeast_anaerobic(m,DATA/'aminoAcid_Bjorkeroth2020.tsv')
    else: m=base.copy()
    m.objective=cfg['objective']; m.reactions.get_by_id(cfg['glc']).bounds=(-10.0,-10.0)
    for _,(rid,_,_) in cfg['acceptors'].items():
        if rid in m.reactions: close_exchange(m,rid)
    if org=='B. subtilis' and 'EX_nh4_e' in m.reactions:
        m.reactions.get_by_id('EX_nh4_e').bounds=(-10.0,1000.0)
    if selected is not None:
        rid=cfg['acceptors'][selected][0]
        if rid not in m.reactions: raise KeyError(f'{rid} absent')
        open_uptake(m,rid,CAP)
    return m


def get_co2(m,sol):
    rid=next((x for x in ['EX_co2_e','Ex_co2','r_1672'] if x in m.reactions),None)
    return float(sol.fluxes[rid])/10.0 if rid else np.nan


def run_org(org):
    base=load_model(org); cfg=CONFIG[org]
    clear_matching_context_cache()
    chem_src=normalize_yeast_for_etn(base) if org=='S. cerevisiae' else base
    chem=chemical_model(chem_src,yeast=(org=='S. cerevisiae'))
    mt=build_metabolite_table(chem); build_matching_context(chem,mt)
    rows=[]
    # Explicit no-candidate-acceptor baseline.
    try:
        m=prepare(base,org,None); sol=solve_growth_pfba(m,cfg['objective']); met,_,_=etn_metrics(chem,sol,10.0)
        baseline={'growth':float(sol.fluxes[cfg['objective']]),'co2_per_glucose':get_co2(m,sol),**met}
        rows.append({'organism':org,'acceptor':'none','exchange':'','status':'baseline','carbon_bearing':False,'capacity':0.0,'actual_acceptor_uptake':0.0,'terminal_e_capacity':0.0,**baseline})
    except Exception as exc:
        baseline=None; rows.append({'organism':org,'acceptor':'none','exchange':'','status':'baseline_failed','carbon_bearing':False,'capacity':0.0,'reason':str(exc)[:250]})
    for acc,(rid,e_per,carbon_bearing) in cfg['acceptors'].items():
        if rid not in base.reactions:
            rows.append({'organism':org,'acceptor':acc,'exchange':rid,'status':'unavailable','carbon_bearing':carbon_bearing,'capacity':CAP,'reason':'exchange absent from reconstruction'})
            continue
        try:
            m=prepare(base,org,acc); sol=solve_growth_pfba(m,cfg['objective']); met,_,_=etn_metrics(chem,sol,10.0)
            uptake=max(0.0,-float(sol.fluxes[rid])); status='functional' if uptake>TOL else 'available_but_unused'
            row={'organism':org,'acceptor':acc,'exchange':rid,'status':status,'carbon_bearing':carbon_bearing,'capacity':CAP,
                 'actual_acceptor_uptake':uptake,'terminal_e_capacity':uptake*e_per,'growth':float(sol.fluxes[cfg['objective']]),'co2_per_glucose':get_co2(m,sol),**met}
            if baseline is not None:
                row['electron_flux_expansion']=met['net_per_glucose']-baseline['net_per_glucose']
                row['additional_carbon_oxidation']=row['co2_per_glucose']-baseline['co2_per_glucose']
                row['growth_increment']=row['growth']-baseline['growth']
            rows.append(row)
        except Exception as exc:
            rows.append({'organism':org,'acceptor':acc,'exchange':rid,'status':'infeasible_or_failed','carbon_bearing':carbon_bearing,'capacity':CAP,'reason':str(exc)[:250]})
    tag=org.replace(' ','_').replace('.','')
    pd.DataFrame(rows).to_csv(OUT/f'extended_acceptor_screen_{tag}.csv',index=False)
    print(org,{k:int(v) for k,v in pd.DataFrame(rows).status.value_counts().items()},flush=True)


def merge():
    frames=[]
    for org in CONFIG:
        tag=org.replace(' ','_').replace('.',''); p=OUT/f'extended_acceptor_screen_{tag}.csv'
        if p.exists(): frames.append(pd.read_csv(p))
    df=pd.concat(frames,ignore_index=True); df.to_csv(OUT/'extended_acceptor_screen.csv',index=False)
    fun=df[df.status=='functional'].copy(); fun.to_csv(OUT/'extended_acceptor_functional_states.csv',index=False)
    clean=fun[~fun.carbon_bearing.astype(bool)].copy()
    stats={
      'n_organisms':len(CONFIG),
      'candidate_capabilities_attempted':int(sum(len(c['acceptors']) for c in CONFIG.values())),
      'functional_capabilities_total':int(len(fun)),
      'functional_noncarbon_capabilities':int(len(clean)),
      'functional_by_organism':{org:sorted(g.acceptor.tolist()) for org,g in fun.groupby('organism')},
      'nonfunctional_or_unavailable_by_organism':{org:sorted(g.acceptor.tolist()) for org,g in df[(df.acceptor!='none') & (df.status!='functional')].groupby('organism')},
    }
    # Cross-state relation: all values are per fixed 10 mmol glucose, so pooled comparison is interpretable.
    for x in ['terminal_e_capacity','additional_carbon_oxidation']:
        t=clean[[x,'electron_flux_expansion']].dropna()
        if len(t)>=3 and t[x].std()>1e-12 and t.electron_flux_expansion.std()>1e-12:
            stats[f'pooled_expansion_vs_{x}']={'n':int(len(t)),'pearson_r':float(pearsonr(t[x],t.electron_flux_expansion)[0]),'spearman_rho':float(spearmanr(t[x],t.electron_flux_expansion)[0])}
    within={}
    for org,g in clean.groupby('organism'):
        if len(g)>=3:
            t=g[['additional_carbon_oxidation','electron_flux_expansion']].dropna()
            if len(t)>=3 and t.additional_carbon_oxidation.std()>1e-12 and t.electron_flux_expansion.std()>1e-12:
                within[org]={'n':int(len(t)),'pearson_r':float(pearsonr(t.additional_carbon_oxidation,t.electron_flux_expansion)[0]),'spearman_rho':float(spearmanr(t.additional_carbon_oxidation,t.electron_flux_expansion)[0])}
    stats['within_organism_expansion_vs_additional_carbon_oxidation']=within
    # Normalize each multi-acceptor architecture to its own maximum response to test
    # whether acceptor chemistry changes amplitude but preserves the oxidation->flux relation.
    scaled=[]
    for org,g in clean.groupby('organism'):
        if len(g)>=3 and g.electron_flux_expansion.max()>0 and g.additional_carbon_oxidation.max()>0:
            q=g.copy()
            q['expansion_scaled']=q.electron_flux_expansion/q.electron_flux_expansion.max()
            q['oxidation_scaled']=q.additional_carbon_oxidation/q.additional_carbon_oxidation.max()
            scaled.append(q)
    if scaled:
        sc=pd.concat(scaled,ignore_index=True)
        sc.to_csv(OUT/'extended_acceptor_scaled_multiarchitecture.csv',index=False)
        stats['scaled_multiarchitecture_expansion_vs_oxidation']={
            'n':int(len(sc)),
            'organisms':sorted(sc.organism.unique().tolist()),
            'pearson_r':float(pearsonr(sc.oxidation_scaled,sc.expansion_scaled)[0]),
            'spearman_rho':float(spearmanr(sc.oxidation_scaled,sc.expansion_scaled)[0])
        }
    (OUT/'extended_acceptor_screen_statistics.json').write_text(json.dumps(stats,indent=2)+'\n')
    print(json.dumps(stats,indent=2))

if __name__=='__main__':
    if '--merge' in sys.argv: merge()
    elif '--organism' in sys.argv: run_org(sys.argv[sys.argv.index('--organism')+1])
    else:
        for org in CONFIG: run_org(org)
        merge()
