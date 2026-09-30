#!/usr/bin/env python3
from pathlib import Path
import pandas as pd, numpy as np, json
ROOT=Path(__file__).resolve().parents[1]; EXP=ROOT/'data/experimental'; OUT=ROOT/'results/publication'
s=pd.read_csv(EXP/'steinsiek2014/steinsiek2014_byproducts_means.csv')
rows=[]
for strain,g in s.groupby('strain'):
 a0=float(g.loc[g.aerobiosis_percent==0,'acetate_mean'].iloc[0]); a150=float(g.loc[g.aerobiosis_percent==150,'acetate_mean'].iloc[0]); rows.append((strain,a0,a150,100*(1-a150/a0)))
pd.DataFrame(rows,columns=['strain','acetate_0pct','acetate_150pct','suppression_pct']).to_csv(OUT/'steinsiek2014_architecture_summary.csv',index=False)
a=pd.read_csv(EXP/'anand2022/anand2022_replicate_phenotypes.csv'); rank={'1H':1,'3H':2,'2H':3,'4H':4}; z=a.groupby(['architecture','state']).agg(growth_mean=('growth','mean'),growth_sd=('growth','std'),acetate_carbon_fraction_mean=('acetate_carbon_fraction','mean'),acetate_carbon_fraction_sd=('acetate_carbon_fraction','std')).reset_index(); z['aerobicity_rank']=z.architecture.map(rank); z.to_csv(OUT/'anand2022_architecture_summary.csv',index=False)
print('Respiratory-architecture support tables generated.')
