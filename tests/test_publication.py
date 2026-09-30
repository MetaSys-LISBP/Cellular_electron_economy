from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def test_cross_species_electron_flux_compression():
    df = pd.read_csv(ROOT / "results/cross_species/electron_flux_compression_summary.csv")
    assert len(df) == 5
    assert (df.low_oxygen_net_e_per_glucose < df.high_respiratory_capacity_net_e_per_glucose).all()
    assert abs(df.low_oxygen_net_e_per_glucose.mean() - 4.10) < 0.01


def test_acceptor_carbon_oxidation_relationship():
    s = json.loads((ROOT / "results/publication/acceptor_carbon_oxidation_statistics.json").read_text())
    assert s["net_electron_flux_vs_CO2_pearson_r"] > 0.996
    assert s["net_electron_flux_vs_CO2_linear_R2"] > 0.993


def test_substrate_core_is_dominant_not_rigid_endpoint():
    core = pd.read_csv(ROOT / "results/publication/substrate_core_redox_contribution.csv")
    assert core.fraction_of_net.min() > 0.948
    mfa = pd.read_csv(ROOT / "results/publication/gonzalez_13c_mfa_electron_flux_metrics.csv").set_index("condition")
    # Experimental anaerobic xylose is higher than the pFBA core/endpoint, as described in the manuscript.
    assert mfa.loc["Anaerobic xylose", "net_e_per_substrate"] > 4.3


def test_gonzalez_13c_mfa_compression_and_flux_balance():
    df = pd.read_csv(ROOT / "results/publication/gonzalez_13c_mfa_electron_flux_metrics.csv").set_index("condition")
    assert df.loc["Anaerobic glucose", "net_e_per_substrate"] < df.loc["Aerobic glucose", "net_e_per_substrate"]
    assert df.loc["Anaerobic xylose", "net_e_per_substrate"] < df.loc["Aerobic xylose", "net_e_per_substrate"]
    assert df.loc["Anaerobic glucose", "q_substrate"] > df.loc["Aerobic glucose", "q_substrate"]
    assert df.loc["Anaerobic glucose", "absolute_net_e_flux"] < df.loc["Aerobic glucose", "absolute_net_e_flux"]
    assert df.source_sink_imbalance_e_per_substrate.abs().max() < 2e-5


def test_effective_transfer_depth_panel_values():
    df = pd.read_csv(ROOT / "results/electron_path_summary_by_condition.csv").set_index("condition")
    assert abs(df.loc["aerobic", "effective_transfer_depth"] - 2.7540685) < 1e-6
    assert abs(df.loc["anaerobic", "effective_transfer_depth"] - 2.1183880) < 1e-6
    assert df.loc["anaerobic", "effective_transfer_depth"] < df.loc["aerobic", "effective_transfer_depth"]


def test_gonzalez_uptake_sd_available_for_decoupling_figure():
    df = pd.read_csv(ROOT / "results/publication/gonzalez_13c_mfa_electron_flux_metrics.csv").set_index("condition")
    expected = {"Aerobic glucose": 0.5, "Anaerobic glucose": 1.0, "Aerobic xylose": 0.5, "Anaerobic xylose": 1.1}
    for condition, sd in expected.items():
        assert abs(df.loc[condition, "q_sd"] - sd) < 1e-12


def test_gonzalez_normalized_oxygen_expansion():
    df = pd.read_csv(ROOT / "results/publication/gonzalez_13c_mfa_electron_flux_metrics.csv").set_index("condition")
    for substrate in ("glucose", "xylose"):
        aer = df.loc[f"Aerobic {substrate}"]
        anox = df.loc[f"Anaerobic {substrate}"]
        assert aer.net_e_per_C > anox.net_e_per_C
        assert aer.net_e_fraction_of_substrate_gamma > anox.net_e_fraction_of_substrate_gamma
    # Similar oxygen-dependent expansion after normalizing the two different substrates.
    dg = df.loc["Aerobic glucose", "net_e_per_C"] - df.loc["Anaerobic glucose", "net_e_per_C"]
    dx = df.loc["Aerobic xylose", "net_e_per_C"] - df.loc["Anaerobic xylose", "net_e_per_C"]
    assert abs(dg - 0.6929156261) < 1e-8
    assert abs(dx - 0.6574732110) < 1e-8


def test_cross_species_baseline_normalization_robustness():
    df = pd.read_csv(ROOT / "results/publication/cross_species_carbon_acceptor_pairs.csv")
    assert len(df) == 22
    assert (df.expansion_e_per_C > 0).all()
    assert (df.expansion_fraction_substrate_gamma > 0).all()
    spans = df.groupby("organism").agg(per_C_min=("low_e_per_C","min"), per_C_max=("low_e_per_C","max"), gamma_min=("low_fraction_substrate_gamma","min"), gamma_max=("low_fraction_substrate_gamma","max"))
    # E. coli, S. enterica and K. phaffii retain appreciable normalized substrate dependence.
    for org in ("E. coli", "S. enterica", "K. phaffii"):
        assert spans.loc[org,"per_C_max"] / spans.loc[org,"per_C_min"] > 1.25
        assert spans.loc[org,"gamma_max"] / spans.loc[org,"gamma_min"] > 1.25
    # In B. subtilis, the raw substrate difference is largely explained by substrate size/reducing content.
    assert spans.loc["B. subtilis","per_C_max"] / spans.loc["B. subtilis","per_C_min"] < 1.05
    assert spans.loc["B. subtilis","gamma_max"] / spans.loc["B. subtilis","gamma_min"] < 1.05
    # The three feasible S. cerevisiae hexose states also largely collapse after normalization.
    assert spans.loc["S. cerevisiae","per_C_max"] / spans.loc["S. cerevisiae","per_C_min"] < 1.05
    assert spans.loc["S. cerevisiae","gamma_max"] / spans.loc["S. cerevisiae","gamma_min"] < 1.05
