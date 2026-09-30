from pathlib import Path
import json
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]

def test_fromanger_endpoint_balance():
    df=pd.read_csv(ROOT/'results/publication/fromanger2010_electron_balance.csv')
    for substrate in ('glucose','xylose'):
        g=df[(df.substrate==substrate)&(df.state_type=='oxygen_limited')].sort_values('specific_OUR_mmol_O2_gDCW_h')
        assert len(g)==3
        assert g.electron_to_O2_pct.iloc[-1] > g.electron_to_O2_pct.iloc[0]
        a=df[(df.substrate==substrate)&(df.state_type=='aerobic_reference')].iloc[0]
        assert a.electron_to_O2_pct > g.electron_to_O2_pct.max()
    s=json.loads((ROOT/'results/publication/fromanger2010_electron_balance_summary.json').read_text())
    assert abs(s['glucose_oxygen_limited_min_pct']-2.0186666667)<1e-8
    assert abs(s['glucose_oxygen_limited_max_pct']-15.91325)<1e-8
    assert abs(s['xylose_oxygen_limited_min_pct']-7.0800833333)<1e-8
    assert abs(s['xylose_oxygen_limited_max_pct']-23.3851666667)<1e-8
    assert s['max_organic_acid_gamma_sensitivity_percentage_points'] < 0.4
    assert s['aerobic_crosscheck_max_abs_difference_percentage_points'] < 1.3
