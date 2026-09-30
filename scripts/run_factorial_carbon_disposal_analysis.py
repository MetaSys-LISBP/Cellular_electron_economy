from pathlib import Path
import pandas as pd, numpy as np, json, sys, multiprocessing as mp
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/publication'; OUT.mkdir(parents=True,exist_ok=True)
subs={'glucose':'EX_glc__D_e','fructose':'EX_fru_e','galactose':'EX_gal_e','xylose':'EX_xyl__D_e','gluconate':'EX_glcn_e','glycerol':'EX_glyc_e'}
caps=[0,1,2,5,10,15,20]

def worker(args):
  sub,ex=args
  import sys, pandas as pd, numpy as np
  import cobra
  from cobra.flux_analysis import pfba
  model=cobra.io.read_sbml_model(str(ROOT/'data/iML1515.xml'))
  met=pd.read_csv(ROOT/'results/metabolite_degree_of_reduction.csv').set_index('met_id'); gamma=met['gamma'].to_dict()
  pairs=pd.read_csv(ROOT/'results/matched_pairs_all.csv'); table=pd.read_csv(ROOT/'results/electron_transfer_table.tsv',sep='\t').set_index('reaction_id')
  specs={}
  for rid,g in pairs.groupby('reaction_id'):
    if rid not in table.index or table.loc[rid,'redox_status']!='redox': continue
    donors=[]; accs=[]
    for _,r in g.iterrows():
      ga=gamma.get(r.reactant); gb=gamma.get(r['product'])
      if ga is None or gb is None or pd.isna(ga) or pd.isna(gb): continue
      delta=float(gb-ga); e=float(r.electrons); d=(r.reactant,r['product'],e)
      if delta<-1e-9: donors.append(d)
      elif delta>1e-9: accs.append(d)
    D=sum(x[2] for x in donors); A=sum(x[2] for x in accs)
    total=float(table.loc[rid,'n_electrons_transferred']) if pd.notna(table.loc[rid,'n_electrons_transferred']) else max(D,A)
    if D<A-1e-9: donors.append((f'unresolved:{rid}',f'unresolved:{rid}',A-D))
    elif A<D-1e-9: accs.append((f'unresolved:{rid}',f'unresolved:{rid}',D-A))
    if donors and accs and total>0: specs[rid]=(donors,accs,total)
  def metrics(sol):
    out={}; inn={}; cum=0.0
    for rid,(donors,accs,total) in specs.items():
      flux=float(sol.fluxes.get(rid,0));
      if abs(flux)<1e-6: continue
      fwd=flux>0
      for d in donors:
       for a in accs:
        ef=abs(flux)*(d[2]*a[2]/total); cum+=ef
        donor=d[0] if fwd else a[1]; accept=a[1] if fwd else d[0]
        out[donor]=out.get(donor,0)+ef; inn[accept]=inn.get(accept,0)+ef
    nodes=set(out)|set(inn); net=sum(max(out.get(n,0)-inn.get(n,0),0) for n in nodes)
    return net,cum,cum/net
  # set non-state-specific bounds once
  allsubs=list(subs.values()); testedacc=['EX_o2_e','EX_no3_e','EX_tmao_e','EX_dmso_e','EX_fum_e']
  for rid in allsubs: model.reactions.get_by_id(rid).bounds=(0,1000)
  for rid in testedacc: model.reactions.get_by_id(rid).bounds=(0,1000)
  model.reactions.get_by_id(ex).bounds=(-10,-10); model.objective='BIOMASS_Ec_iML1515_core_75p37M'
  met_ex=next(mm for mm in model.reactions.get_by_id(ex).metabolites if mm.compartment=='e'); C=int(met_ex.elements.get('C',0)); gam=float(gamma.get(met_ex.id,np.nan))
  rows=[]
  for cap in caps:
    model.reactions.get_by_id('EX_o2_e').bounds=(-cap,1000) if cap else (0,1000)
    sol=pfba(model)
    net,cum,depth=metrics(sol)
    rows.append(dict(substrate=sub,oxygen_cap=cap,growth=float(sol.fluxes['BIOMASS_Ec_iML1515_core_75p37M']),actual_o2=max(0,-float(sol.fluxes['EX_o2_e'])),co2=float(sol.fluxes['EX_co2_e']),net_e_flux=net,net_per_substrate=net/10,cumulative_per_substrate=cum/10,effective_transfer_depth=depth,C=C,gamma=gam,net_per_C=net/10/C,net_fraction_gamma=net/10/gam))
  return rows

if __name__=='__main__':
  with mp.Pool(processes=6) as pool:
    chunks=pool.map(worker, list(subs.items()))
  rows=[r for c in chunks for r in c]; df=pd.DataFrame(rows); df.to_csv(OUT/'factorial_carbon_oxygen_states.csv',index=False)
  scaled_rows=[]
  for sub,g in df.groupby('substrate'):
    g=g.sort_values('oxygen_cap'); b=float(g[g.oxygen_cap==0].net_per_substrate.iloc[0]); full=float(g[g.oxygen_cap==20].net_per_substrate.iloc[0]-b)
    for _,r in g.iterrows(): scaled_rows.append(dict(substrate=sub,oxygen_cap=r.oxygen_cap,scaled_expansion=(r.net_per_substrate-b)/full,baseline=b,full_expansion=full))
  scaled=pd.DataFrame(scaled_rows); scaled.to_csv(OUT/'factorial_scaled_expansion_curves.csv',index=False)
  order=list(subs); non=[1,2,5,10,15,20]
  M=np.array([[float(df[(df.substrate==s)&(df.oxygen_cap==c)].net_per_substrate.iloc[0]-df[(df.substrate==s)&(df.oxygen_cap==0)].net_per_substrate.iloc[0]) for c in non] for s in order])
  sv=np.linalg.svd(M,compute_uv=False); rank1=float(sv[0]**2/(sv**2).sum()); corrs=[]
  for i,s1 in enumerate(order):
    for s2 in order[i+1:]:
      a=scaled[(scaled.substrate==s1)&(scaled.oxygen_cap>0)].scaled_expansion.values; b=scaled[(scaled.substrate==s2)&(scaled.oxygen_cap>0)].scaled_expansion.values; corrs.append(float(np.corrcoef(a,b)[0,1]))
  def anova(col):
    P=df.pivot(index='substrate',columns='oxygen_cap',values=col).loc[order,caps].values; grand=P.mean(); ssS=len(caps)*((P.mean(1)-grand)**2).sum(); ssO=len(order)*((P.mean(0)-grand)**2).sum(); ssT=((P-grand)**2).sum(); ssI=ssT-ssS-ssO; return {'substrate_pct':100*ssS/ssT,'oxygen_pct':100*ssO/ssT,'interaction_pct':100*ssI/ssT}
  stats={'rank1_expansion_variance_fraction':rank1,'rank1_expansion_variance_pct':100*rank1,'scaled_curve_min_pairwise_r':min(corrs),'scaled_curve_mean_pairwise_r':float(np.mean(corrs)),'raw':anova('net_per_substrate'),'per_carbon':anova('net_per_C'),'per_gamma':anova('net_fraction_gamma')}
  (OUT/'factorial_carbon_oxygen_statistics.json').write_text(json.dumps(stats,indent=2)+'\n'); print(json.dumps(stats,indent=2))
