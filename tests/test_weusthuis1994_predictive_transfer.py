from pathlib import Path
import json
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def test_weusthuis_prediction_contract():
    s=json.loads((ROOT/'results/publication/weusthuis1994_predictive_transfer_statistics.json').read_text())
    assert s['n_heldout']==5
    assert round(s['r'],3)==0.988
    assert round(s['predictive_r2'],3)==0.956
    assert round(s['rmse_mmol_gDW_h'],2)==0.58
    h=pd.read_csv(ROOT/'results/publication/weusthuis1994_maltose_heldout_predictions.csv')
    assert list(h['oxygen_in_mmol_L_h'])==[5.4,10.8,18.7,22.3,24.8]
