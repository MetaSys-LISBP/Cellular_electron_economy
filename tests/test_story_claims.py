from pathlib import Path
import json, pandas as pd, numpy as np
ROOT=Path(__file__).resolve().parents[1]
def test_factorial_headline_values():
 s=json.loads((ROOT/'results/publication/factorial_carbon_oxygen_statistics.json').read_text())
 assert 99.4 < s['rank1_expansion_variance_pct'] < 99.6
 assert 0.95 < s['scaled_curve_min_pairwise_r'] < 0.97
 assert 4.5 < s['raw']['interaction_pct'] < 4.8
def test_ecom4la_headline_values():
 d=pd.read_csv(ROOT/'results/publication/ecom4la_state_constrained_summary.csv')
 assert np.isclose(d.iloc[0].e_per_glucose,10.33,atol=.05)
 assert np.isclose(d.iloc[1].e_per_glucose,4.13,atol=.05)
