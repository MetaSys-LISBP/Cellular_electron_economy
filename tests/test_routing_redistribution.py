import json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

def test_cross_species_routing_redistribution_release_output():
    p = ROOT / 'results' / 'publication' / 'cross_species_routing_redistribution.csv'
    df = pd.read_csv(p)
    assert len(df) == 5
    assert set(df['organism']) == {'E. coli','S. enterica','B. subtilis','S. cerevisiae','K. phaffii'}
    assert df['jensen_shannon_distance_bits'].between(0.70, 0.90).all()
    assert df['total_variation_distance'].between(0.65, 0.90).all()
    # The robustness metric confirms extensive rerouting but does not replace
    # the directional transfer-depth/NADH metrics used in main Fig. 4b.
    assert df['transfer_depth_change_pct'].min() < -20
    assert df['transfer_depth_change_pct'].max() > 5

def test_cross_species_routing_statistics_consistent():
    df = pd.read_csv(ROOT / 'results' / 'publication' / 'cross_species_routing_redistribution.csv')
    s = json.loads((ROOT / 'results' / 'publication' / 'cross_species_routing_redistribution_statistics.json').read_text())
    assert s['n_organisms'] == 5
    assert abs(s['jsd_min'] - df['jensen_shannon_distance_bits'].min()) < 1e-12
    assert abs(s['jsd_max'] - df['jensen_shannon_distance_bits'].max()) < 1e-12
