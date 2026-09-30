from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def test_jouhten_condition_means_and_sample_sd():
    df = pd.read_csv(ROOT / "results/publication/jouhten2008_13c_mfa_electron_source_metrics.csv").set_index("oxygen_inlet_percent")
    expected_mean = {20.9: 10.28, 2.8: 10.10, 1.0: 6.61, 0.5: 4.96, 0.0: 3.64}
    expected_sd = {20.9: 0.2545584412, 2.8: 0.3959797975, 1.0: 0.0707106781, 0.5: 0.4242640687, 0.0: 0.0}
    for oxygen, value in expected_mean.items():
        assert abs(df.loc[oxygen, "mean_electron_source_per_glucose"] - value) < 1e-10
        assert abs(df.loc[oxygen, "sd_electron_source_per_glucose"] - expected_sd[oxygen]) < 1e-9


def test_jouhten_gradient_compresses_while_glucose_uptake_increases():
    df = pd.read_csv(ROOT / "results/publication/jouhten2008_13c_mfa_electron_source_metrics.csv").set_index("oxygen_inlet_percent")
    assert df.loc[0.0, "mean_electron_source_per_glucose"] < df.loc[20.9, "mean_electron_source_per_glucose"]
    assert df.loc[0.0, "mean_glucose_uptake_mmol_gCDW_h"] > df.loc[20.9, "mean_glucose_uptake_mmol_gCDW_h"]
    compression = 100 * (1 - df.loc[0.0, "mean_electron_source_per_glucose"] / df.loc[20.9, "mean_electron_source_per_glucose"])
    assert abs(compression - 64.5914396887) < 1e-8
