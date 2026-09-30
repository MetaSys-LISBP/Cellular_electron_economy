import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_steinsiek_predictive_transfer_contract():
    out = ROOT / "results/publication"
    stats = json.loads((out / "steinsiek2014_predictive_transfer_statistics.json").read_text())
    points = pd.read_csv(out / "steinsiek2014_predictive_transfer_points.csv")
    held = points[points.role == "held_out_prediction"]
    assert len(held) == 12
    assert round(stats["predictive_r2"], 3) == 0.897
    assert round(stats["rmse_mmol_gDW_h"], 2) == 0.48
    assert stats["uncertainty_propagation"]["status"] == "not_recomputed"


def test_steinsiek_robustness_outputs_present():
    out = ROOT / "results/publication"
    for name in [
        "steinsiek2014_predictive_transfer_source_choice.csv",
        "steinsiek2014_predictive_transfer_leave_one_architecture_out.csv",
        "steinsiek2014_predictive_transfer_source_level_omission.csv",
    ]:
        assert (out / name).exists(), name
