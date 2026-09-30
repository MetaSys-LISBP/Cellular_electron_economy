#!/usr/bin/env python3
"""Systematically test encoded terminal electron acceptors on a common glucose background.

The screen is capability-based rather than a forced balanced factorial: every curated
respiratory acceptor exchange encoded by each reconstruction is attempted with oxygen
closed, alongside an oxygen series and a no-terminal-acceptor reference.  Results are
reported even when the organism cannot sustain positive growth.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from cobra.io import read_sbml_model

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
OUT=ROOT/'results'/'publication'; OUT.mkdir(parents=True,exist_ok=True)
DATA=ROOT/'data'/'cross_species'
from etn.data_acquisition import load_ecoli_model
from etn.cross_species import load_bacillus_xlsx, load_standard_gem_xlsx, normalize_yeast_for_etn, chemical_model, etn_metrics, solve_growth_pfba, apply_yeast_anaerobic
from etn.matching import build_metabolite_table, build_matching_context, clear_matching_context_cache

CAPS=[0.0,20.0]
CONFIG={
 'E. coli': {'objective':'BIOMASS_Ec_iML1515_core_75p37M','glc':'EX_glc__D_e','o2':'EX_o2_e','acceptors':{'oxygen':('EX_o2_e',4.0),'nitrate':('EX_no3_e',2.0),'TMAO':('EX_tmao_e',2.0),'DMSO':('EX_dmso_e',2.0),'fumarate':('EX_fum_e',2.0)},'pichia':False},
 'B. subtilis': {'objective':'BiomassRepsolRed','glc':'EX_glc__D_e','o2':'EX_o2_e','acceptors':{'oxygen':('EX_o2_e',4.0),'nitrate':('EX_no3_e',2.0),'fumarate':('EX_fum_e',2.0)},'pichia':False},
 'S. enterica': {'objective':'BIOMASS_iRR1083_1','glc':'EX_glc__D_e','o2':'EX_o2_e','acceptors':{'oxygen':('EX_o2_e',4.0),'nitrate':('EX_no3_e',2.0),'TMAO':('EX_tmao_e',2.0),'DMSO':('EX_dmso_e',2.0),'fumarate':('EX_fum_e',2.0)},'pichia':False},
 'K. phaffii': {'objective':'Ex_biomass','glc':'Ex_glc_D','o2':'Ex_o2','acceptors':{'oxygen':('Ex_o2',4.0),'nitrate':('Ex_no3',2.0),'fumarate':('Ex_fum',2.0)},'pichia':True},
 'S. cerevisiae': {'objective':'r_2111','glc':'r_1714','o2':'r_1992','acceptors':{'oxygen':('r_1992',4.0),'fumarate':('r_1798',2.0)},'pichia':False},
}

def load_model(org):
 if org=='E. coli': return load_ecoli_model()
 if org=='B. subtilis': return load_bacillus_xlsx(DATA/'iBB1018.xlsx')
 if org=='S. enterica': return read_sbml_model(str(DATA/'iYS1720.xml.gz'))
 if org=='K. phaffii': return load_standard_gem_xlsx(DATA/'iMT1026_model.xlsx','iMT1026')
 if org=='S. cerevisiae': return load_standard_gem_xlsx(DATA/'yeast-GEM-v9.1.0.xlsx','yeastGEM_v9.1.0')
 raise KeyError(org)

def close_exchange(m,rid):
 if rid not in m.reactions: return
 r=m.reactions.get_by_id(rid)
 # All curated exchanges use negative flux for uptake in these reconstructions.
 r.bounds=(0.0,max(0.0,float(r.upper_bound)))

def open_uptake(m,rid,cap):
 r=m.reactions.get_by_id(rid)
 if cap<=0: r.bounds=(0.0,max(0.0,float(r.upper_bound)))
 else: r.bounds=(-float(cap),max(0.0,float(r.upper_bound)))

def prepare(base,org,acceptor,cap):
 cfg=CONFIG[org]
 # Use the model-provided anaerobic edits for yeast whenever oxygen is absent.
 if org=='S. cerevisiae' and acceptor!='oxygen':
  m=base.copy(); apply_yeast_anaerobic(m, DATA/'aminoAcid_Bjorkeroth2020.tsv')
 else: m=base.copy()
 m.objective=cfg['objective']
 # Fixed glucose input.
 g=m.reactions.get_by_id(cfg['glc']); g.bounds=(-10.0,-10.0)
 # Close every curated terminal acceptor, then open only the selected one.
 for _name,(rid,_e) in cfg['acceptors'].items(): close_exchange(m,rid)
 # B. subtilis publication medium retains defined ammonium so nitrate is not needed as N source.
 if org=='B. subtilis' and 'EX_nh4_e' in m.reactions: m.reactions.get_by_id('EX_nh4_e').bounds=(-10.0,1000.0)
 rid,_=cfg['acceptors'][acceptor]
 open_uptake(m,rid,cap)
 return m

def run_org(org):
 base=load_model(org); cfg=CONFIG[org]
 clear_matching_context_cache()
 chem_src=normalize_yeast_for_etn(base) if org=='S. cerevisiae' else base
 chem=chemical_model(chem_src,yeast=(org=='S. cerevisiae'))
 mt=build_metabolite_table(chem); build_matching_context(chem,mt)
 rows=[]
 for acc,(rid,e_per) in cfg['acceptors'].items():
  if rid not in base.reactions:
   rows.append({'organism':org,'acceptor':acc,'exchange':rid,'capacity':np.nan,'status':'unavailable','reason':'exchange absent from reconstruction'})
   continue
  for cap in CAPS:
   # cap=0 is retained for each acceptor so feasibility is explicit; duplicates are harmless.
   try:
    m=prepare(base,org,acc,cap)
    sol=solve_growth_pfba(m,cfg['objective'])
    metrics,_,_=etn_metrics(chem,sol,10.0)
    flux=float(sol.fluxes[rid]); uptake=max(0.0,-flux)
    co2_id=next((x for x in ['EX_co2_e','Ex_co2','r_1672'] if x in m.reactions),None)
    co2=float(sol.fluxes[co2_id]) if co2_id else np.nan
    rows.append({'organism':org,'acceptor':acc,'exchange':rid,'capacity':cap,'status':'feasible','reason':'','growth':float(sol.fluxes[cfg['objective']]),'actual_acceptor_uptake':uptake,'theoretical_terminal_e_capacity':uptake*e_per,'co2_flux':co2,'co2_per_glucose':co2/10.0,**metrics})
   except Exception as exc:
    rows.append({'organism':org,'acceptor':acc,'exchange':rid,'capacity':cap,'status':'infeasible_or_failed','reason':str(exc)[:250]})
 df=pd.DataFrame(rows)
 tag=org.replace(' ','_').replace('.','')
 df.to_csv(OUT/f'systematic_acceptor_screen_{tag}.csv',index=False)
 print(org,df.groupby(['acceptor','status']).size().to_dict(),flush=True)

def merge():
 frames=[]
 for org in CONFIG:
  tag=org.replace(' ','_').replace('.','')
  p=OUT/f'systematic_acceptor_screen_{tag}.csv'
  if p.exists(): frames.append(pd.read_csv(p))
 if not frames: raise RuntimeError('no acceptor screen outputs')
 df=pd.concat(frames,ignore_index=True)
 df.to_csv(OUT/'systematic_acceptor_screen.csv',index=False)
 # One baseline per organism: feasible no-acceptor cap=0 where available. Use the median
 # across duplicate cap=0 rows, which should be numerically identical.
 summary=[]
 for org,gorg in df.groupby('organism'):
  base=gorg[(gorg.capacity==0)&(gorg.status=='feasible')]
  base_net=float(base.net_per_glucose.median()) if len(base) else np.nan
  base_co2=float(base.co2_per_glucose.median()) if len(base) else np.nan
  for acc,g in gorg[(gorg.capacity>0)&(gorg.status=='feasible')].groupby('acceptor'):
   # strongest tested feasible state for this acceptor
   r=g.sort_values('capacity').iloc[-1]
   summary.append({'organism':org,'acceptor':acc,'max_tested_capacity':r.capacity,'actual_acceptor_uptake':r.actual_acceptor_uptake,'theoretical_terminal_e_capacity':r.theoretical_terminal_e_capacity,'net_e_per_glucose':r.net_per_glucose,'baseline_net_e_per_glucose':base_net,'electron_flux_expansion':r.net_per_glucose-base_net if np.isfinite(base_net) else np.nan,'co2_per_glucose':r.co2_per_glucose,'baseline_co2_per_glucose':base_co2,'additional_carbon_oxidation':r.co2_per_glucose-base_co2 if np.isfinite(base_co2) else np.nan,'growth':r.growth})
 summ=pd.DataFrame(summary)
 summ.to_csv(OUT/'systematic_acceptor_screen_summary.csv',index=False)
 # Pooled mechanistic correlations exclude fumarate because the acceptor itself is carbon-bearing
 # and therefore confounds CO2-based carbon-oxidation accounting.
 test=summ[(summ.acceptor!='fumarate') & summ.electron_flux_expansion.notna()]
 stats={'n_organisms_attempted':len(CONFIG),'n_acceptor_capabilities_feasible':int(len(summ)),'feasible_acceptors_by_organism':{org:sorted(g.acceptor.unique().tolist()) for org,g in summ.groupby('organism')}}
 for x in ['theoretical_terminal_e_capacity','additional_carbon_oxidation']:
  t=test[[x,'electron_flux_expansion']].dropna()
  if len(t)>=3 and t[x].std()>1e-12 and t.electron_flux_expansion.std()>1e-12:
   stats[f'expansion_vs_{x}']={'n':int(len(t)),'pearson_r':float(pearsonr(t[x],t.electron_flux_expansion)[0]),'spearman_rho':float(spearmanr(t[x],t.electron_flux_expansion)[0])}
 (OUT/'systematic_acceptor_screen_statistics.json').write_text(json.dumps(stats,indent=2)+'\n')
 print(json.dumps(stats,indent=2))

if __name__=='__main__':
 if '--merge' in sys.argv: merge()
 elif '--organism' in sys.argv: run_org(sys.argv[sys.argv.index('--organism')+1])
 else:
  for org in CONFIG: run_org(org)
  merge()
