#!/usr/bin/env python3
"""Predictive transfer of an oxygen-response shape from glucose to maltose in Weusthuis et al. (1994)."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression

ROOT=Path(__file__).resolve().parents[1]
IN=ROOT/'data'/'experimental'/'weusthuis1994'/'s_cerevisiae_oxygen_series.csv'
OUT=ROOT/'results'/'publication'; OUT.mkdir(parents=True,exist_ok=True)

def analyse(censored_glucose_value_mM=0.0, write_detail=False):
    d=pd.read_csv(IN).copy()
    cens=(d['substrate'].eq('Glucose') & d['ethanol_out_mM'].isna() & d['ethanol_out_censored_lt_mM'].notna())
    d.loc[cens,'ethanol_out_mM']=float(censored_glucose_value_mM)
    ethanol_feed_mM=12.0
    d['q_ethanol']=d['dilution_h']*(d['ethanol_out_mM']-ethanol_feed_mM)/d['dry_weight_g_L']
    glc=d[d['substrate'].eq('Glucose') & d['ethanol_out_mM'].notna()].sort_values('oxygen_in_mmol_L_h').copy()
    mal=d[d['substrate'].eq('Maltose') & d['ethanol_out_mM'].notna()].sort_values('oxygen_in_mmol_L_h').copy()
    xg=glc['oxygen_in_mmol_L_h'].to_numpy(float); yg=glc['q_ethanol'].to_numpy(float)
    iso=IsotonicRegression(increasing=False,out_of_bounds='clip').fit(xg,yg)
    glc['q_ethanol_isotonic']=iso.predict(xg)
    glc_low=float(glc.loc[np.isclose(glc['oxygen_in_mmol_L_h'],0.0),'q_ethanol'].iloc[0])
    glc_high=float(iso.predict([30.8])[0])
    def f(x):
        yy=iso.predict(np.asarray(x,dtype=float))
        return (yy-glc_low)/(glc_high-glc_low)
    mal_zero=mal[np.isclose(mal['oxygen_in_mmol_L_h'],0.0)]
    mal_low=float(mal_zero['q_ethanol'].mean()); mal_low_sd=float(mal_zero['q_ethanol'].std(ddof=1))
    mal_high=float(mal.loc[np.isclose(mal['oxygen_in_mmol_L_h'],30.8),'q_ethanol'].iloc[0])
    held_x=np.array([5.4,10.8,18.7,22.3,24.8])
    held=mal[mal['oxygen_in_mmol_L_h'].isin(held_x)].sort_values('oxygen_in_mmol_L_h').copy()
    held['q_pred']=mal_low+(mal_high-mal_low)*f(held['oxygen_in_mmol_L_h'].to_numpy(float))
    obs=held['q_ethanol'].to_numpy(float); pred=held['q_pred'].to_numpy(float)
    r=float(np.corrcoef(obs,pred)[0,1])
    r2=float(1-np.sum((obs-pred)**2)/np.sum((obs-obs.mean())**2))
    rmse=float(np.sqrt(np.mean((obs-pred)**2)))
    curve_x=np.r_[np.arange(0.0,31.0,2.0),30.8]
    curve=pd.DataFrame({'oxygen_in_mmol_L_h':curve_x})
    curve['q_pred_maltose']=mal_low+(mal_high-mal_low)*f(curve_x)
    if write_detail:
        glc_out=glc[['oxygen_in_mmol_L_h','ethanol_out_mM','ethanol_out_censored_lt_mM','dry_weight_g_L','dilution_h','q_ethanol','q_ethanol_isotonic']].copy()
        glc_out['censored_glucose_imputation_mM']=float(censored_glucose_value_mM)
        glc_out.to_csv(OUT/'weusthuis1994_glucose_response.csv',index=False)
        anchors=pd.DataFrame([
            {'substrate':'Maltose','role':'low_anchor','oxygen_in_mmol_L_h':0.0,'q_ethanol':mal_low,'q_ethanol_sd':mal_low_sd},
            {'substrate':'Maltose','role':'high_anchor','oxygen_in_mmol_L_h':30.8,'q_ethanol':mal_high,'q_ethanol_sd':np.nan},
        ])
        anchors.to_csv(OUT/'weusthuis1994_maltose_anchors.csv',index=False)
        held[['oxygen_in_mmol_L_h','q_ethanol','q_pred']].to_csv(OUT/'weusthuis1994_maltose_heldout_predictions.csv',index=False)
        curve.to_csv(OUT/'weusthuis1994_prediction_curve.csv',index=False)
        all_rows=[]
        for _,rr in anchors.iterrows():
            all_rows.append({'series':'maltose','role':rr.role,'oxygen_in_mmol_L_h':rr.oxygen_in_mmol_L_h,'observed_q_ethanol':rr.q_ethanol,'predicted_q_ethanol':rr.q_ethanol,'sd_q_ethanol':rr.q_ethanol_sd})
        for _,rr in held.iterrows():
            all_rows.append({'series':'maltose','role':'held_out','oxygen_in_mmol_L_h':rr.oxygen_in_mmol_L_h,'observed_q_ethanol':rr.q_ethanol,'predicted_q_ethanol':rr.q_pred,'sd_q_ethanol':np.nan})
        pd.DataFrame(all_rows).to_csv(OUT/'weusthuis1994_predictive_transfer_points.csv',index=False)
    return {'r':r,'predictive_r2':r2,'rmse_mmol_gDW_h':rmse,'n_heldout':len(held),'glucose_censored_imputation_mM':float(censored_glucose_value_mM),'maltose_low_anchor_q':mal_low,'maltose_high_anchor_q':mal_high}

primary=analyse(0.0,write_detail=True)
sens=[]
for val in [0.0,0.25,0.5,0.75,1.0]:
    s=analyse(val,write_detail=False); sens.append(s)
pd.DataFrame(sens).to_csv(OUT/'weusthuis1994_censoring_sensitivity.csv',index=False)
(OUT/'weusthuis1994_predictive_transfer_statistics.json').write_text(json.dumps(primary,indent=2)+'\n')
assert round(primary['r'],3)==0.988
assert round(primary['predictive_r2'],3)==0.956
assert round(primary['rmse_mmol_gDW_h'],2)==0.58
print(json.dumps(primary,indent=2))
