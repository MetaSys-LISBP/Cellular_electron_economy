#!/usr/bin/env python3
"""Calculate electron generation from Jouhten et al. (2008) S. cerevisiae 13C-MFA flux distributions.

The published Figure 2 reports two replicate net-flux distributions normalized to
100 glucose for five inlet-O2 conditions. Additional file 4 defines reaction
stoichiometry. Only source-generating redox reactions are required because the
original MFA formulation did not constrain NADH/NADPH balances; the resulting
observable is electron generation, not complete source-to-sink electron flux.
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "jouhten2008"
OUT = ROOT / "results" / "publication"
OUT.mkdir(parents=True, exist_ok=True)

EXPECTED_O2 = [20.9, 2.8, 1.0, 0.5, 0.0]


def main() -> None:
    flux = pd.read_csv(DATA / "jouhten2008_source_reaction_fluxes.csv")
    sto = pd.read_csv(DATA / "jouhten2008_reaction_electron_stoichiometry.csv")
    meta = pd.read_csv(DATA / "jouhten2008_condition_metadata.csv")

    if set(flux["reaction_id"]) != set(sto["reaction_id"]):
        raise AssertionError("Jouhten flux and electron-stoichiometry reaction sets differ")
    if set(flux["replicate"]) != {"I", "II"}:
        raise AssertionError("Expected exactly replicates I and II")

    merged = flux.merge(sto[["reaction_id", "electron_source_coefficient"]], on="reaction_id", validate="many_to_one")
    merged["electron_source_per_100_glucose"] = merged["flux_per_100_glucose"] * merged["electron_source_coefficient"]
    merged["electron_source_per_glucose"] = merged["electron_source_per_100_glucose"] / 100.0
    merged.to_csv(OUT / "jouhten2008_13c_mfa_reaction_electron_accounting.csv", index=False)

    rep = (merged.groupby(["oxygen_inlet_percent", "replicate"], as_index=False)
                  ["electron_source_per_glucose"].sum())
    rep = rep.merge(meta, on=["oxygen_inlet_percent", "replicate"], how="left", validate="one_to_one")
    rep.to_csv(OUT / "jouhten2008_13c_mfa_replicate_metrics.csv", index=False)

    rows = []
    for o2 in EXPECTED_O2:
        s = rep[rep["oxygen_inlet_percent"] == o2].sort_values("replicate")
        if len(s) != 2:
            raise AssertionError(f"Expected two Jouhten replicates at {o2}% O2")
        vals = s["electron_source_per_glucose"].to_numpy(float)
        qvals = s["glucose_uptake_mmol_gCDW_h"].to_numpy(float)
        rows.append({
            "oxygen_inlet_percent": o2,
            "replicate_I_electron_source_per_glucose": float(s.loc[s.replicate == "I", "electron_source_per_glucose"].iloc[0]),
            "replicate_II_electron_source_per_glucose": float(s.loc[s.replicate == "II", "electron_source_per_glucose"].iloc[0]),
            "mean_electron_source_per_glucose": float(vals.mean()),
            "sd_electron_source_per_glucose": float(vals.std(ddof=1)),
            "mean_glucose_uptake_mmol_gCDW_h": float(qvals.mean()),
            "sd_between_replicate_glucose_uptake_mmol_gCDW_h": float(qvals.std(ddof=1)),
        })
    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUT / "jouhten2008_13c_mfa_electron_source_metrics.csv", index=False)

    expected = {20.9: 10.28, 2.8: 10.10, 1.0: 6.61, 0.5: 4.96, 0.0: 3.64}
    for o2, val in expected.items():
        got = float(metrics.loc[metrics.oxygen_inlet_percent == o2, "mean_electron_source_per_glucose"].iloc[0])
        if not np.isclose(got, val, atol=1e-12):
            raise AssertionError(f"Jouhten {o2}% O2 mean: expected {val}, got {got}")

    checks = {
        "n_conditions": int(len(metrics)),
        "replicates_per_condition": 2,
        "aerobic_mean_e_per_glucose": float(metrics.loc[metrics.oxygen_inlet_percent == 20.9, "mean_electron_source_per_glucose"].iloc[0]),
        "aerobic_sd_e_per_glucose": float(metrics.loc[metrics.oxygen_inlet_percent == 20.9, "sd_electron_source_per_glucose"].iloc[0]),
        "anoxic_mean_e_per_glucose": float(metrics.loc[metrics.oxygen_inlet_percent == 0.0, "mean_electron_source_per_glucose"].iloc[0]),
        "anoxic_sd_e_per_glucose": float(metrics.loc[metrics.oxygen_inlet_percent == 0.0, "sd_electron_source_per_glucose"].iloc[0]),
        "compression_percent": float(100.0 * (1.0 - metrics.loc[metrics.oxygen_inlet_percent == 0.0, "mean_electron_source_per_glucose"].iloc[0] / metrics.loc[metrics.oxygen_inlet_percent == 20.9, "mean_electron_source_per_glucose"].iloc[0])),
    }
    (OUT / "jouhten2008_13c_mfa_validation_checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(metrics.to_string(index=False))
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
