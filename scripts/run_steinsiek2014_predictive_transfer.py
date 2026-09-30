#!/usr/bin/env python3
"""Predict acetate responses across respiratory architectures.

The response shape is inferred from one architecture after normalization
between its 0% and 150% aerobiosis endpoints. Target architectures contribute
only their two endpoints; the four intermediate states are held out.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/experimental/steinsiek2014/steinsiek2014_byproducts_means.csv"
OUT = ROOT / "results/publication"

LEVELS = [0, 20, 50, 80, 100, 150]
HELD_OUT = [20, 50, 80, 100]
PARENT = "TBE029"
TARGETS = ["TBE031_bo_only", "TBE032_bdII_only", "TBE042_bdI_only"]


def metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    observed = np.asarray(observed, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    ss_res = float(np.sum((observed - predicted) ** 2))
    ss_tot = float(np.sum((observed - observed.mean()) ** 2))
    return {
        "n": int(observed.size),
        "pearson_r": float(np.corrcoef(observed, predicted)[0, 1]),
        "predictive_r2": float(1.0 - ss_res / ss_tot),
        "rmse_mmol_gDW_h": float(np.sqrt(np.mean((observed - predicted) ** 2))),
    }


def normalized_response(series: pd.Series, high: int = 150) -> pd.Series:
    denominator = float(series.loc[high] - series.loc[0])
    if abs(denominator) < 1e-12:
        raise ValueError("Identical endpoint values cannot define a response shape")
    return (series - float(series.loc[0])) / denominator


def predict_from_shape(target: pd.Series, shape: pd.Series, high: int = 150) -> pd.Series:
    return float(target.loc[0]) + (float(target.loc[high]) - float(target.loc[0])) * shape


def two_way_decomposition(matrix: np.ndarray) -> dict[str, float]:
    grand = matrix.mean()
    row_effect = matrix.mean(axis=1, keepdims=True) - grand
    col_effect = matrix.mean(axis=0, keepdims=True) - grand
    residual = matrix - grand - row_effect - col_effect
    ss_total = float(np.sum((matrix - grand) ** 2))
    return {
        "architecture_pct": 100.0 * float(matrix.shape[1] * np.sum(row_effect ** 2)) / ss_total,
        "aerobiosis_pct": 100.0 * float(matrix.shape[0] * np.sum(col_effect ** 2)) / ss_total,
        "interaction_pct": 100.0 * float(np.sum(residual ** 2)) / ss_total,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(INPUT)
    pivot = data.pivot(index="strain", columns="aerobiosis_percent", values="acetate_mean").loc[:, LEVELS]

    source_shape = normalized_response(pivot.loc[PARENT])
    point_rows: list[dict[str, object]] = []
    for strain in TARGETS:
        pred = predict_from_shape(pivot.loc[strain], source_shape)
        for level in LEVELS:
            point_rows.append({
                "source_architecture": PARENT,
                "target_architecture": strain,
                "aerobiosis_percent": level,
                "observed_acetate_mmol_gDW_h": float(pivot.loc[strain, level]),
                "predicted_acetate_mmol_gDW_h": float(pred.loc[level]),
                "role": "scaling_endpoint" if level in (0, 150) else "held_out_prediction",
            })
    points = pd.DataFrame(point_rows)
    held = points[points.role == "held_out_prediction"]
    primary = metrics(held.observed_acetate_mmol_gDW_h, held.predicted_acetate_mmol_gDW_h)

    # Alternative high endpoint: 100% aerobiosis.
    shape_100 = normalized_response(pivot.loc[PARENT], high=100)
    alt_obs, alt_pred = [], []
    for strain in TARGETS:
        pred = predict_from_shape(pivot.loc[strain], shape_100, high=100)
        for level in [20, 50, 80]:
            alt_obs.append(float(pivot.loc[strain, level]))
            alt_pred.append(float(pred.loc[level]))
    endpoint_100 = metrics(np.asarray(alt_obs), np.asarray(alt_pred))

    # Use each architecture in turn as the source response and predict all others.
    source_rows = []
    for source in pivot.index:
        shape = normalized_response(pivot.loc[source])
        obs, pred = [], []
        for target in [s for s in pivot.index if s != source]:
            target_pred = predict_from_shape(pivot.loc[target], shape)
            obs.extend(pivot.loc[target, HELD_OUT].astype(float).tolist())
            pred.extend(target_pred.loc[HELD_OUT].astype(float).tolist())
        source_rows.append({"source_architecture": source, **metrics(np.asarray(obs), np.asarray(pred))})
    source_choice = pd.DataFrame(source_rows)

    # Leave-one-architecture-out: predict each target from the mean normalized
    # response of the other three architectures.
    loo_rows = []
    for target in pivot.index:
        donors = [s for s in pivot.index if s != target]
        mean_shape = pd.concat([normalized_response(pivot.loc[s]) for s in donors], axis=1).mean(axis=1)
        pred = predict_from_shape(pivot.loc[target], mean_shape)
        for level in HELD_OUT:
            loo_rows.append({
                "target_architecture": target,
                "aerobiosis_percent": level,
                "observed_acetate_mmol_gDW_h": float(pivot.loc[target, level]),
                "predicted_acetate_mmol_gDW_h": float(pred.loc[level]),
            })
    loo = pd.DataFrame(loo_rows)
    loo_metrics = metrics(loo.observed_acetate_mmol_gDW_h, loo.predicted_acetate_mmol_gDW_h)

    # Omit each intermediate point from the parental response and replace it by
    # linear interpolation between the nearest retained levels.
    omission_rows = []
    for omitted in HELD_OUT:
        kept = [x for x in LEVELS if x != omitted]
        shape = source_shape.copy()
        shape.loc[omitted] = float(np.interp(omitted, kept, source_shape.loc[kept]))
        obs, pred = [], []
        for target in TARGETS:
            target_pred = predict_from_shape(pivot.loc[target], shape)
            obs.extend(pivot.loc[target, HELD_OUT].astype(float).tolist())
            pred.extend(target_pred.loc[HELD_OUT].astype(float).tolist())
        omission_rows.append({"omitted_source_level_percent": omitted, **metrics(np.asarray(obs), np.asarray(pred))})
    omission = pd.DataFrame(omission_rows)

    # Endpoint-independent summaries of shared response structure.
    baseline_subtracted = pivot.loc[:, [20, 50, 80, 100, 150]].to_numpy(float)
    baseline_subtracted = pivot.loc[:, [0]].to_numpy(float) - baseline_subtracted
    singular_values = np.linalg.svd(baseline_subtracted, full_matrices=False, compute_uv=False)
    rank_one_pct = 100.0 * float(singular_values[0] ** 2 / np.sum(singular_values ** 2))
    variance = two_way_decomposition(pivot.to_numpy(float))

    statistics = {
        **primary,
        "source_architecture": PARENT,
        "target_architectures": TARGETS,
        "source_endpoint_percent": [0, 150],
        "held_out_aerobiosis_percent": HELD_OUT,
        "alternative_100_percent_endpoint": endpoint_100,
        "leave_one_architecture_out": loo_metrics,
        "rank_one_baseline_subtracted_variance_pct": rank_one_pct,
        "two_way_variance_pct": variance,
        "uncertainty_propagation": {
            "status": "not_recomputed",
            "reason": "The supplied Steinsiek input contains means but not pointwise acetate uncertainties.",
        },
    }

    points.to_csv(OUT / "steinsiek2014_predictive_transfer_points.csv", index=False)
    source_choice.to_csv(OUT / "steinsiek2014_predictive_transfer_source_choice.csv", index=False)
    loo.to_csv(OUT / "steinsiek2014_predictive_transfer_leave_one_architecture_out.csv", index=False)
    omission.to_csv(OUT / "steinsiek2014_predictive_transfer_source_level_omission.csv", index=False)
    (OUT / "steinsiek2014_predictive_transfer_statistics.json").write_text(
        json.dumps(statistics, indent=2) + "\n", encoding="utf-8"
    )

    assert round(primary["predictive_r2"], 3) == 0.897
    assert round(primary["rmse_mmol_gDW_h"], 2) == 0.48


if __name__ == "__main__":
    main()
