from pathlib import Path
import json
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
PUB=ROOT/'results'/'publication'

def test_extended_carbon_entry_screen():
    f=pd.read_csv(PUB/'extended_carbon_entry_feasibility.csv')
    assert len(f)==30  # 5 organisms x 3 additional substrates x 2 electron-disposal-capacity states
    s=pd.read_csv(PUB/'extended_carbon_entry_summary.csv')
    assert set(s.substrate)=={'acetate','pyruvate','succinate'}
    paired=s[s.low_feasible.astype(bool) & s.high_feasible.astype(bool)]
    assert len(paired)==1
    r=paired.iloc[0]
    assert r.organism=='E. coli' and r.substrate=='pyruvate'
    assert r.growth_ratio_low_to_high < 0.02

def test_extended_acceptor_screen():
    st=json.loads((PUB/'extended_acceptor_screen_statistics.json').read_text())
    assert st['candidate_capabilities_attempted']==24
    assert st['functional_capabilities_total']==21
    assert st['functional_noncarbon_capabilities']==17
    sc=st['scaled_multiarchitecture_expansion_vs_oxidation']
    assert sc['n']==15
    assert set(sc['organisms'])=={'E. coli','B. subtilis','S. enterica'}
    assert sc['pearson_r'] > 0.99
    assert sc['spearman_rho'] > 0.95

def test_systematic_carbon_screen_is_complete():
    f=pd.read_csv(PUB/'cross_species_carbon_acceptor_feasibility.csv')
    assert len(f)==60  # 30 combinations x 2 states
    pairs=pd.read_csv(PUB/'cross_species_carbon_acceptor_pairs.csv')
    assert len(pairs)==22
    assert 'S. cerevisiae' in set(pairs.organism)
    assert (pairs.high_net_per_substrate > pairs.low_net_per_substrate).all()
