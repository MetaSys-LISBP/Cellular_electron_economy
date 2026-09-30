#!/usr/bin/env python3
"""Descriptive decomposition of cross-species low-respiratory-capacity baselines.

This analysis is intentionally supplementary: the 5x6 screen is incomplete
because some organism/substrate states are biologically/model-infeasible.
We therefore use least-squares additive models only to ask how much of the
observed feasible-state variation is captured by organism and substrate terms;
we do not assign inferential P values or interpret a saturated interaction.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
IN=ROOT/'results/publication/cross_species_carbon_acceptor_pairs.csv'
OUT=ROOT/'results/publication/cross_species_baseline_decomposition.json'

def design(df, include_org=True, include_sub=True):
    cols=[np.ones(len(df))]
    if include_org:
        for v in sorted(df.organism.unique())[1:]: cols.append((df.organism==v).astype(float).to_numpy())
    if include_sub:
        for v in sorted(df.substrate.unique())[1:]: cols.append((df.substrate==v).astype(float).to_numpy())
    return np.column_stack(cols)

def fit_r2(df,ycol,io=True,isub=True):
    y=df[ycol].to_numpy(float); X=design(df,io,isub); b=np.linalg.lstsq(X,y,rcond=None)[0]; yh=X@b
    rss=float(np.sum((y-yh)**2)); tss=float(np.sum((y-y.mean())**2)); return {'r2':1-rss/tss if tss else 1.,'rss':rss,'rank':int(np.linalg.matrix_rank(X))}

def main():
    df=pd.read_csv(IN)
    out={'n_feasible_pairs':int(len(df)),'note':'Descriptive least-squares models on the incomplete feasible-state matrix; no inferential P values.'}
    for y in ['low_net_per_substrate','low_e_per_C','low_fraction_substrate_gamma']:
        full=fit_r2(df,y,True,True); org=fit_r2(df,y,True,False); sub=fit_r2(df,y,False,True); intercept=fit_r2(df,y,False,False)
        out[y]={'organism_only':org,'substrate_only':sub,'additive_organism_plus_substrate':full,
                'incremental_R2_of_substrate_given_organism':full['r2']-org['r2'],
                'incremental_R2_of_organism_given_substrate':full['r2']-sub['r2']}
    OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__': main()
