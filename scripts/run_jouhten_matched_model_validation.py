#!/usr/bin/env python3
"""Matched genome-scale / 13C-MFA validation for Jouhten et al. (2008).

The experimental observable is electron generation calculated from the published
13C-MFA central-carbon flux maps. The genome-scale model is constrained only by
measured extracellular physiology: mean glucose uptake from the two labelled
replicates and the published mean oxygen-uptake rate for each inlet-O2 condition.
Internal 13C fluxes are never used to constrain the genome-scale model.

Published average specific oxygen uptake rates (mmol O2 g biomass^-1 h^-1):
20.9, 2.8, 1.0, 0.5, 0.0% inlet O2 -> 2.7, 2.5, 1.7, 1.2, 0.0.
Source: Jouhten et al., BMC Systems Biology 2008, 2:60,
doi:10.1186/1752-0509-2-60.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUT = ROOT / "results" / "publication"
OUT.mkdir(parents=True, exist_ok=True)

from etn.cross_species import (
    load_standard_gem_xlsx,
    normalize_yeast_for_etn,
    apply_yeast_anaerobic,
    chemical_model,
    etn_metrics,
    solve_growth_pfba,
)
from etn.matching import build_metabolite_table, build_matching_context, clear_matching_context_cache

OUR = {20.9: 2.7, 2.8: 2.5, 1.0: 1.7, 0.5: 1.2, 0.0: 0.0}
OBJECTIVE = "r_2111"
GLC_EX = "r_1714"
O2_EX = "r_1992"


def main() -> None:
    exp = pd.read_csv(OUT / "jouhten2008_13c_mfa_electron_source_metrics.csv")
    expected = [20.9, 2.8, 1.0, 0.5, 0.0]
    assert exp["oxygen_inlet_percent"].tolist() == expected

    model_path = ROOT / "data" / "cross_species" / "yeast-GEM-v9.1.0.xlsx"
    base = load_standard_gem_xlsx(model_path, "yeastGEM_v9.1.0")
    base.objective = OBJECTIVE

    clear_matching_context_cache()
    chem_source = normalize_yeast_for_etn(base)
    chem = chemical_model(chem_source, yeast=True)
    mt = build_metabolite_table(chem)
    build_matching_context(chem, mt)

    rows = []
    for _, r in exp.iterrows():
        inlet = float(r.oxygen_inlet_percent)
        qglc = float(r.mean_glucose_uptake_mmol_gCDW_h)
        our = float(OUR[inlet])

        m = base.copy()
        m.objective = OBJECTIVE
        if inlet == 0.0:
            apply_yeast_anaerobic(m, ROOT / "data" / "cross_species" / "aminoAcid_Bjorkeroth2020.tsv")

        # Fix measured mean glucose uptake. Other common carbon-source exchanges
        # remain at their reconstruction defaults (non-uptake unless supplied).
        glc = m.reactions.get_by_id(GLC_EX)
        glc.bounds = (-qglc, -qglc)

        o2 = m.reactions.get_by_id(O2_EX)
        if our > 0:
            o2.bounds = (-our, max(0.0, float(o2.upper_bound)))
        else:
            o2.bounds = (0.0, max(0.0, float(o2.upper_bound)))

        sol = solve_growth_pfba(m, OBJECTIVE)
        metrics, _, _ = etn_metrics(chem, sol, qglc)
        predicted = float(metrics["net_per_glucose"])
        observed = float(r.mean_electron_source_per_glucose)
        obs_sd = float(r.sd_electron_source_per_glucose)
        growth = float(sol.fluxes[OBJECTIVE])
        actual_o2 = max(0.0, -float(sol.fluxes[O2_EX]))

        rows.append({
            "oxygen_inlet_percent": inlet,
            "published_mean_OUR_mmol_gCDW_h": our,
            "mean_glucose_uptake_mmol_gCDW_h": qglc,
            "experimental_source_generation_e_per_glucose": observed,
            "experimental_sd_e_per_glucose": obs_sd,
            "model_source_generation_e_per_glucose": predicted,
            "absolute_error_e_per_glucose": abs(predicted - observed),
            "relative_error_percent": 100.0 * (predicted - observed) / observed,
            "predicted_growth_h-1": growth,
            "actual_model_O2_uptake_mmol_gCDW_h": actual_o2,
        })

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "jouhten2008_matched_model_experiment.csv", index=False)

    y = res.experimental_source_generation_e_per_glucose.to_numpy(float)
    yhat = res.model_source_generation_e_per_glucose.to_numpy(float)
    pr, pp = pearsonr(y, yhat)
    rmse = float(np.sqrt(np.mean((yhat-y)**2)))
    nrmse_range = float(rmse / (y.max()-y.min()))
    mape = float(np.mean(np.abs((yhat-y)/y))*100.0)
    # Regression of predicted versus experimental, not a fitted correction.
    slope, intercept = np.polyfit(y, yhat, 1)
    ss_res = float(np.sum((yhat-y)**2))
    ss_tot = float(np.sum((y-y.mean())**2))
    r2_identity = float(1.0 - ss_res/ss_tot)

    stats = {
        "n_conditions": int(len(res)),
        "pearson_r": float(pr),
        "pearson_p": float(pp),
        "rmse_e_per_glucose": rmse,
        "nrmse_experimental_range": nrmse_range,
        "mean_absolute_percent_error": mape,
        "predicted_vs_experimental_regression_slope": float(slope),
        "predicted_vs_experimental_regression_intercept": float(intercept),
        "identity_R2": r2_identity,
        "model_constraint_note": "Model constrained by measured mean glucose uptake and published mean OUR. The internal 13C fluxes provide an independent experimental comparison.",
    }
    (OUT / "jouhten2008_matched_model_experiment_statistics.json").write_text(json.dumps(stats, indent=2)+"\n")

    print(res.to_string(index=False))
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
