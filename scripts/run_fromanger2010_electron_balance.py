#!/usr/bin/env python3
"""Whole-cell electron-equivalent analysis of Fromanger et al. 2010 C. shehatae data."""
from pathlib import Path
import json
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
IN=ROOT/'data/experimental/fromanger2010/fromanger2010_carbon_partition.csv'
OUT=ROOT/'results/publication'; OUT.mkdir(parents=True,exist_ok=True)

SUBSTRATE_GAMMA_PER_C=4.0
GAMMA_PER_C={'ethanol':6.0,'xylitol':4.4,'arabitol':4.4,'ribitol':4.4,'glycerol':14.0/3.0,'biomass':4.34}
ACID_MID=3.25; ACID_LOW=3.0; ACID_HIGH=3.5

def calculate(row, acid_gamma):
    retained=0.0
    for key,gamma in GAMMA_PER_C.items():
        retained += float(row[f'{key}_Cmol_per_CmolS'])*gamma
    retained += float(row['organic_acids_Cmol_per_CmolS'])*acid_gamma
    frac=(SUBSTRATE_GAMMA_PER_C-retained)/SUBSTRATE_GAMMA_PER_C
    return retained, frac

df=pd.read_csv(IN)
out=[]
for _,r in df.iterrows():
    retained, frac=calculate(r,ACID_MID)
    _, frac_a=calculate(r,ACID_LOW)
    _, frac_b=calculate(r,ACID_HIGH)
    lo=min(frac_a,frac_b); hi=max(frac_a,frac_b)
    d=r.to_dict()
    d.update({
      'substrate_gamma_e_per_C':SUBSTRATE_GAMMA_PER_C,
      'retained_e_equiv_per_Cmol_substrate':retained,
      'electron_to_O2_fraction':frac,
      'electron_to_O2_pct':100*frac,
      'electron_to_O2_pct_acid_gamma_low':100*lo,
      'electron_to_O2_pct_acid_gamma_high':100*hi,
      'acid_gamma_per_C_midpoint':ACID_MID,
    })
    if pd.notna(r.get('CO2_Cmol_per_CmolS')) and pd.notna(r.get('RQ_mol_CO2_per_mol_O2')):
        rqfrac=float(r['CO2_Cmol_per_CmolS'])/float(r['RQ_mol_CO2_per_mol_O2'])
        d['aerobic_O2_fraction_from_CO2_RQ']=rqfrac
        d['aerobic_O2_pct_from_CO2_RQ']=100*rqfrac
        d['aerobic_balance_minus_RQ_percentage_points']=100*(frac-rqfrac)
    else:
        d['aerobic_O2_fraction_from_CO2_RQ']=float('nan')
        d['aerobic_O2_pct_from_CO2_RQ']=float('nan')
        d['aerobic_balance_minus_RQ_percentage_points']=float('nan')
    out.append(d)
res=pd.DataFrame(out)
res.to_csv(OUT/'fromanger2010_electron_balance.csv',index=False)

lim=res[res.state_type=='oxygen_limited']
summary={
 'source':'Fromanger et al. 2010, DOI 10.1007/s10295-009-0688-7',
 'observable':'whole-cell fraction of substrate electron equivalents delivered to O2 inferred by degree-of-reduction balance',
 'glucose_oxygen_limited_min_pct':float(lim[lim.substrate=='glucose'].electron_to_O2_pct.min()),
 'glucose_oxygen_limited_max_pct':float(lim[lim.substrate=='glucose'].electron_to_O2_pct.max()),
 'xylose_oxygen_limited_min_pct':float(lim[lim.substrate=='xylose'].electron_to_O2_pct.min()),
 'xylose_oxygen_limited_max_pct':float(lim[lim.substrate=='xylose'].electron_to_O2_pct.max()),
 'glucose_aerobic_pct':float(res[(res.substrate=='glucose')&(res.state_type=='aerobic_reference')].electron_to_O2_pct.iloc[0]),
 'xylose_aerobic_pct':float(res[(res.substrate=='xylose')&(res.state_type=='aerobic_reference')].electron_to_O2_pct.iloc[0]),
 'max_organic_acid_gamma_sensitivity_percentage_points':float((res.electron_to_O2_pct_acid_gamma_high-res.electron_to_O2_pct_acid_gamma_low).abs().max()),
 'aerobic_crosscheck_max_abs_difference_percentage_points':float(res.aerobic_balance_minus_RQ_percentage_points.abs().max()),
 'reported_redox_balance_closure_limit_oxygen_limited_pct':10.0,
 'reported_redox_balance_closure_limit_aerobic_pct':6.0,
}
(OUT/'fromanger2010_electron_balance_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(res[['substrate','condition','specific_OUR_mmol_O2_gDCW_h','electron_to_O2_pct','aerobic_O2_pct_from_CO2_RQ']].to_string(index=False))
print(json.dumps(summary,indent=2))
