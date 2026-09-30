from pathlib import Path
import csv
import re

ROOT=Path(__file__).resolve().parents[1]

def test_publication_source_data_contract():
    d=ROOT/'source_data'/'publication'
    required=[
        'Fig2a_condition_summary.csv','Fig2b_carrier_relay.csv','Fig2e_Gonzalez_13C_MFA_electron_flux_metrics.csv',
        'Fig3ab_acceptor_titration.csv','Fig3f_factorial_carbon_oxygen_states.csv',
        'Fig4a_cross_species_flux_compression.csv','Fig4e_Steinsiek2014_predictive_transfer.csv',
        'Fig4f_Weusthuis1994_predictive_transfer.csv',
        'SuppFigS2_iML1515_ETN_nodes.csv','SuppFigS2_iML1515_ETN_edges.csv',
        'SuppFigS1_classification_summary.csv','SuppFigS6_ECOM4LA_state_constrained_summary.csv',
        'SuppFigS10_Toya2012_source_totals.csv','SuppFigS14_Baumann2010_oxygen_gradient.csv',
        'SuppFigS15_extended_acceptor_functional_states.csv','SuppFigS17_extended_carbon_entry_summary.csv',
        'SuppFigS18_Fromanger2010_electron_balance.csv','SuppFigS19_Steinsiek2014_predictive_transfer.csv',
        'SuppFigS20_Weusthuis1994_glucose_response.csv',
    ]
    missing=[x for x in required if not (d/x).exists()]
    assert not missing, missing
    supp_numbers={int(m.group(1)) for p in d.glob('SuppFigS*.csv') if (m:=re.match(r'SuppFigS(\d+)_',p.name))}
    assert supp_numbers == set(range(1,21)), supp_numbers
    assert not any(d.glob('Fig5*'))
    with (ROOT/'FIGURE_SOURCE_DATA_MANIFEST.csv').open(newline='') as fh:
        manifest_files = {row['file'] for row in csv.DictReader(fh)}
    actual_files = {p.name for p in d.glob('*.csv')}
    assert actual_files == manifest_files

def test_publication_figure_contract():
    d=ROOT/'figures'/'publication'
    stems=['Fig2_glucose_decoupling_and_disposal_capacity',
           'Fig3_carbon_source_disposal_capacity_ecoli','Fig4_cross_species_organizing_principles']
    for stem in stems:
        assert (d/f'{stem}.png').exists(), stem
        assert (d/f'{stem}.pdf').exists(), stem
    assert not (d/'Fig1_electron_economy.png').exists()
    assert not (d/'Fig1_electron_economy.pdf').exists()
    names=[p.name for p in d.glob('*.png')]
    supp_numbers={int(m.group(1)) for n in names if (m:=re.match(r'SuppFigS(\d+)_',n))}
    assert supp_numbers == set(range(1,21)), supp_numbers
    for suffix in ('png','pdf','svg'):
        assert (d/f'SuppFigS2_iML1515_ETN.{suffix}').exists()
    assert (d/'SuppFigS14_Baumann2010_oxygen_gradient.png').exists()
    assert (d/'SuppFigS15_extended_terminal_acceptor_generality.png').exists()
    assert (d/'SuppFigS17_extended_carbon_entry_coverage.png').exists()
    assert (d/'SuppFigS18_fromanger2010_electron_balance.png').exists()
    assert (d/'SuppFigS19_Steinsiek2014_predictive_transfer.png').exists()
    assert (d/'SuppFigS20_Weusthuis1994_predictive_transfer.png').exists()

def test_complete_etn_edge_contract():
    import csv, json
    edge_path = ROOT/'results'/'electron_transfer_edges.csv'
    with edge_path.open(newline='') as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 280
    fhl = [r for r in rows if r['reaction_id'] == 'FHL']
    assert len(fhl) == 1
    assert fhl[0]['direction'] == 'forward'
    assert fhl[0]['model_direction_allowed'] == 'False'
    stats = json.loads((ROOT/'results'/'network_basic_stats.json').read_text())
    assert stats['n_edges_simple'] == 263
    assert stats['n_edges_multi'] == 279
    assert stats['n_reaction_edges_complete'] == 280
    assert stats['n_blocked_reaction_edges_complete'] == 1
