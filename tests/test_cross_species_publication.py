from pathlib import Path
import pandas as pd
from etn.matching import _metabolite_base_id

ROOT = Path(__file__).resolve().parents[1]


def test_compartment_stripping_handles_eukaryotic_suffixes():
    assert _metabolite_base_id("nadh_c") == "nadh"
    assert _metabolite_base_id("nadh_m") == "nadh"
    assert _metabolite_base_id("nadph_er") == "nadph"
    assert _metabolite_base_id("q6h2_m") == "q6h2"


def test_cross_species_publication_endpoint_convergence():
    path = ROOT / "results" / "cross_species" / "electron_flux_compression_summary.csv"
    df = pd.read_csv(path)
    assert set(df.organism) == {"E. coli", "S. cerevisiae", "B. subtilis", "S. enterica", "K. phaffii"}
    assert (df.low_oxygen_net_e_per_glucose < df.high_respiratory_capacity_net_e_per_glucose).all()
    assert df.low_oxygen_net_e_per_glucose.between(3.9, 4.3).all()
    assert 1.5 < 100 * df.low_oxygen_net_e_per_glucose.std(ddof=1) / df.low_oxygen_net_e_per_glucose.mean() < 3.0



def test_cross_species_routing_responses_are_not_universal():
    path = ROOT / "results" / "cross_species" / "electron_flux_compression_summary.csv"
    df = pd.read_csv(path).set_index("organism")
    assert {"transfer_depth_change_pct", "nadh_relay_change_pct"}.issubset(df.columns)
    # E. coli retains/slightly increases NADH relay, whereas the other four models reduce it strongly.
    assert df.loc["E. coli", "nadh_relay_change_pct"] > 0
    assert (df.drop(index="E. coli").nadh_relay_change_pct < -40).all()
    # Effective transfer depth decreases in bacteria such as E. coli/S. enterica but rises in both yeasts.
    assert df.loc["E. coli", "transfer_depth_change_pct"] < -15
    assert df.loc["S. enterica", "transfer_depth_change_pct"] < -15
    assert df.loc["S. cerevisiae", "transfer_depth_change_pct"] > 0
    assert df.loc["K. phaffii", "transfer_depth_change_pct"] > 0
