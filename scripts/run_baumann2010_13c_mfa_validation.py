#!/usr/bin/env python3
"""Reproduce the Baumann et al. K. phaffii oxygen-gradient electron-generation summary."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
IN=ROOT/'data'/'experimental'/'baumann2010'/'baumann2010_condition_metrics.csv'
OUT=ROOT/'results'/'publication'; OUT.mkdir(parents=True,exist_ok=True)

df=pd.read_csv(IN)
df['recalculated_source_per_glucose']=df['electron_source_mmol_e_gDCW_h']/df['glucose_uptake_mmol_gDCW_h']
df['observable']='electron generation calculated from 13C-MFA flux distributions'
df['balanced_source_to_sink_flux']=False
cols=['strain','oxygen_inlet_percent','glucose_uptake_mmol_gDCW_h','electron_source_mmol_e_gDCW_h','electron_source_per_glucose','propagated_sd_e_per_glucose','recalculated_source_per_glucose','observable','balanced_source_to_sink_flux']
df[cols].to_csv(OUT/'baumann2010_13c_mfa_electron_source_metrics.csv',index=False)
rows=[]
for strain,g in df.groupby('strain',sort=False):
    hi=g.loc[g.oxygen_inlet_percent.idxmax()]; lo=g.loc[g.oxygen_inlet_percent.idxmin()]
    rows.append({
        'strain':strain,'high_oxygen_percent':hi.oxygen_inlet_percent,'low_oxygen_percent':lo.oxygen_inlet_percent,
        'source_per_glucose_high':hi.electron_source_per_glucose,'source_per_glucose_low':lo.electron_source_per_glucose,
        'source_compression_percent':100*(hi.electron_source_per_glucose-lo.electron_source_per_glucose)/hi.electron_source_per_glucose,
        'glucose_uptake_high':hi.glucose_uptake_mmol_gDCW_h,'glucose_uptake_low':lo.glucose_uptake_mmol_gDCW_h,
        'glucose_uptake_increase_percent':100*(lo.glucose_uptake_mmol_gDCW_h-hi.glucose_uptake_mmol_gDCW_h)/hi.glucose_uptake_mmol_gDCW_h,
        'absolute_source_change_percent':100*(lo.electron_source_mmol_e_gDCW_h-hi.electron_source_mmol_e_gDCW_h)/hi.electron_source_mmol_e_gDCW_h,
    })
changes=pd.DataFrame(rows)
changes.to_csv(OUT/'baumann2010_13c_mfa_headline_changes.csv',index=False)
summary={
    'source':'Baumann et al. BMC Systems Biology 4, 141 (2010)',
    'doi':'10.1186/1752-0509-4-141','observable':'electron generation','not_balanced_source_to_sink':True,
}
for _,r in changes.iterrows():
    key='control' if r['strain']=='Control' else 'fab'
    summary[f'{key}_compression_percent']=float(r.source_compression_percent)
    summary[f'{key}_glucose_uptake_increase_percent']=float(r.glucose_uptake_increase_percent)
    summary[f'{key}_absolute_source_change_percent']=float(r.absolute_source_change_percent)
(OUT/'baumann2010_13c_mfa_validation_checks.json').write_text(json.dumps(summary,indent=2)+'\n')
assert np.allclose(df.electron_source_per_glucose,df.recalculated_source_per_glucose,rtol=0,atol=5e-3)
print('Wrote Baumann oxygen-gradient validation outputs')
