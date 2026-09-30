from pathlib import Path
import json
import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1]

def test_toya_source_totals_and_labels():
    t=pd.read_csv(ROOT/'results/publication/toya2012_electron_source_totals.csv').set_index('condition')
    assert 'WT anaerobic' in t.index
    assert t.loc['WT anaerobic','source_e_per_glucose'] == pytest.approx(4.270, abs=1e-9)
    assert t.loc['WT nitrate','source_e_per_glucose'] == pytest.approx(6.384, abs=1e-9)
    assert t.loc['ArcA nitrate','source_e_per_glucose'] == pytest.approx(8.458, abs=1e-9)
    assert t.loc['WT nitrate','conservative_low'] > t.loc['WT anaerobic','conservative_high']

def test_toya_decomposition_sums():
    t=pd.read_csv(ROOT/'results/publication/toya2012_electron_source_totals.csv').set_index('condition')
    d=pd.read_csv(ROOT/'results/publication/toya2012_electron_source_decomposition.csv')
    sums=d.groupby('condition').e_per_glucose.sum()
    for c,v in sums.items(): assert v == pytest.approx(t.loc[c,'source_e_per_glucose'], abs=1e-12)

def test_denby_terminal_delivery_and_matched_model():
    d=json.loads((ROOT/'results/publication/denby2015_summary.json').read_text())
    assert d['measured_terminal_delivery_e_per_glucose'] == pytest.approx(4.58, abs=1e-12)
    assert d['measured_terminal_delivery_e_sd_per_glucose'] == pytest.approx(0.14, abs=1e-12)
    assert d['reported_source_generation_e_per_glucose'] == pytest.approx(5.44, abs=1e-12)
    assert d['model_fermentation_net_e_per_glucose'] == pytest.approx(4.103284122741754, rel=1e-9)
    assert d['model_TMAO_net_e_per_glucose'] == pytest.approx(6.0425007113670315, rel=1e-9)

def test_perrenoud_is_explicitly_qualitative():
    p=pd.read_csv(ROOT/'results/publication/perrenoud2005_acceptor_comparison.csv').set_index('experimental_condition')
    assert p.loc['nitrate','reported_respiratory_cyclic_TCA_flux_mmol_g_h'] == pytest.approx(0.2)
    assert p.loc['DMSO','reported_respiratory_cyclic_TCA_flux_mmol_g_h'] == pytest.approx(0.0)
    assert p.loc['nitrate','external_acceptor_mM'] == 40.0
    assert p.loc['nitrate','model_standardized_uptake_cap'] == 10.0
