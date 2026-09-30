#!/usr/bin/env python3
"""Independent experimental validation of non-O2 terminal-electron-acceptor effects.

Reproduces three publication-facing analyses from curated primary-source transcriptions:
- Toya et al. 2012: 13C-constrained central fluxes, WT anaerobic (O2−) vs WT nitrate and ArcA nitrate;
- Denby et al. 2015: anaerobic MG1655 fermentation -> TMAO respiration;
- Perrenoud & Sauer 2005: qualitative nitrate-versus-DMSO respiratory TCA comparison.

Important observable distinctions are preserved:
Toya -> electron generation calculated from reported central-carbon fluxes, not complete source-to-sink flux.
Denby -> measured TMAO terminal delivery and a published redox/source balance; model-matched
         source-to-sink flux is reported separately.
Perrenoud -> qualitative routing comparison only; extracellular concentration is never
             equated to model uptake capacity.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUT = ROOT / "results" / "publication"
OUT.mkdir(parents=True, exist_ok=True)
DATA = ROOT / "data" / "experimental"

from etn.data_acquisition import load_ecoli_model
from etn.cross_species import chemical_model, etn_metrics, solve_growth_pfba
from scripts.run_mechanistic_analysis import BIOMASS, GLUCOSE, close_tested_acceptors, solve_acceptor


def run_toya() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    src = pd.read_csv(DATA / "toya2012" / "toya2012_table_SII_redox_subset.csv")
    sinks = pd.read_csv(DATA / "toya2012" / "toya2012_table_SII_reduced_product_sinks.csv")
    # Supplementary Table SII uses the O2− notation for the WT anaerobic state.
    # The main article identifies the corresponding comparison as anaerobic versus nitrate respiration.
    specs = [
        ("WT anaerobic", "wt_no_o2"),
        ("WT nitrate", "wt_nitrate"),
        ("ArcA nitrate", "arca_nitrate"),
    ]
    totals, details, sink_rows = [], [], []
    for label, key in specs:
        e = src["electrons_per_flux"].astype(float).to_numpy()
        flux = src[key].astype(float).to_numpy()
        lo = src[key + "_low"].astype(float).to_numpy()
        hi = src[key + "_high"].astype(float).to_numpy()
        # Only the oxidative direction contributes to electron generation. A negative SDH flux in the
        # no-O2 state is fumarate reduction and is therefore an electron sink.
        contrib = np.maximum(flux, 0.0) * e / 100.0
        contrib_lo = np.maximum(lo, 0.0) * e / 100.0
        contrib_hi = np.maximum(hi, 0.0) * e / 100.0
        totals.append({
            "condition": label,
            "source_e_per_glucose": float(contrib.sum()),
            "conservative_low": float(contrib_lo.sum()),
            "conservative_high": float(contrib_hi.sum()),
        })
        for cls in ["Glycolysis", "Oxidative PPP", "PDH", "Oxidative TCA"]:
            mask = src["reaction_class"].eq(cls).to_numpy()
            details.append({
                "condition": label,
                "source_class": cls,
                "e_per_glucose": float(contrib[mask].sum()),
            })
        sink_value = float(
            (sinks[key].astype(float) * sinks["sink_electrons_per_flux"].astype(float) / 100.0).sum()
        )
        sink_rows.append({"condition": label, "explicit_fermentation_sink_e_per_glucose": sink_value})

    totals = pd.DataFrame(totals)
    details = pd.DataFrame(details)
    sink_df = pd.DataFrame(sink_rows)
    totals.to_csv(OUT / "toya2012_electron_source_totals.csv", index=False)
    details.to_csv(OUT / "toya2012_electron_source_decomposition.csv", index=False)
    sink_df.to_csv(OUT / "toya2012_explicit_fermentation_sinks.csv", index=False)

    wt = totals.set_index("condition")
    stats = {
        "wt_nitrate_vs_no_o2_change_pct": 100.0 * (wt.loc["WT nitrate", "source_e_per_glucose"] / wt.loc["WT anaerobic", "source_e_per_glucose"] - 1.0),
        "conservative_minimum_change_pct": 100.0 * (wt.loc["WT nitrate", "conservative_low"] / wt.loc["WT anaerobic", "conservative_high"] - 1.0),
        "arca_vs_wt_nitrate_change_pct": 100.0 * (wt.loc["ArcA nitrate", "source_e_per_glucose"] / wt.loc["WT nitrate", "source_e_per_glucose"] - 1.0),
    }
    (OUT / "toya2012_statistics.json").write_text(json.dumps(stats, indent=2) + "\n")
    return totals, details, sink_df


def _denby_model_states(base, chem) -> pd.DataFrame:
    # Denby steady state: 45.8 +/- 1.4 mM TMA from 20 mM feed glucose = 2.29 +/- 0.07
    # TMAO turnovers per glucose. With model glucose fixed at 10, matched TMAO uptake is 22.9.
    specs = [("fermentation", 0.0), ("tmao_matched", 22.9)]
    rows = []
    for label, cap in specs:
        model = base.copy()
        model.objective = BIOMASS
        model.reactions.get_by_id(GLUCOSE).bounds = (-10.0, -10.0)
        close_tested_acceptors(model)
        if cap:
            model.reactions.get_by_id("EX_tmao_e").bounds = (-cap, 1000.0)
        sol = solve_growth_pfba(model, BIOMASS)
        metrics, _, stats = etn_metrics(chem, sol, 10.0)
        row = {"condition": label, "tmao_uptake_cap": cap, "growth": float(sol.fluxes[BIOMASS]), **metrics}
        for name, rid in [
            ("acetate", "EX_ac_e"), ("formate", "EX_for_e"), ("ethanol", "EX_etoh_e"),
            ("succinate", "EX_succ_e"), ("co2", "EX_co2_e"), ("TMA", "EX_tma_e"),
        ]:
            row[name] = float(sol.fluxes[rid])
        tma_stats = stats[stats["node"].astype(str).str.contains("tma_p", case=False, regex=False)]
        row["tma_terminal_delivery_e_flux"] = float(tma_stats["inflow"].sum()) if len(tma_stats) else 0.0
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "denby2015_model_matched_states.csv", index=False)
    return out


def run_denby(base, chem) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    steady = pd.read_csv(DATA / "denby2015" / "denby2015_table1_steady_states.csv").set_index("variable")
    redox = pd.read_csv(DATA / "denby2015" / "denby2015_fig1_redox_balance.csv").set_index("quantity")
    glucose = float(redox.loc["feed_glucose_mM", "value"])
    tma = float(steady.loc["TMA", "tmao_state_mean"])
    tma_sd = float(steady.loc["TMA", "tmao_state_sd"])
    glycolysis_re = float(redox.loc["glycolysis_RE", "value"])
    pdh_re = float(redox.loc["PDH_RE", "value"])
    orotate_re = float(redox.loc["orotate_RE", "value"])
    succ_re = float(redox.loc["succinate_sink_RE", "value"])
    ethanol_re = float(redox.loc["ethanol_sink_RE", "value"])
    tmao_re = float(redox.loc["TMAO_sink_RE", "value"])
    direct = {
        "measured_TMA_per_glucose": tma / glucose,
        "measured_TMA_per_glucose_sd": tma_sd / glucose,
        "measured_terminal_delivery_e_per_glucose": 2.0 * tma / glucose,
        "measured_terminal_delivery_e_sd_per_glucose": 2.0 * tma_sd / glucose,
        # Denby expresses RE as two-electron reducing equivalents; convert to electrons.
        "reported_source_generation_e_per_glucose": 2.0 * (glycolysis_re + pdh_re + orotate_re) / glucose,
        "reported_internal_sink_e_per_glucose": 2.0 * (succ_re + ethanol_re) / glucose,
        "reported_TMAO_sink_e_per_glucose": 2.0 * tmao_re / glucose,
    }
    (OUT / "denby2015_experimental_electron_balance.json").write_text(json.dumps(direct, indent=2) + "\n")

    model = _denby_model_states(base, chem)
    m = model.set_index("condition")
    products = []
    for metabolite in ["acetate", "formate", "ethanol", "succinate"]:
        e0 = float(steady.loc[metabolite, "fermentation_mean"])
        e1 = float(steady.loc[metabolite, "tmao_state_mean"])
        m0 = float(m.loc["fermentation", metabolite]) / 10.0
        m1 = float(m.loc["tmao_matched", metabolite]) / 10.0
        products.append({
            "metabolite": metabolite,
            "experimental_fermentation_mM": e0,
            "experimental_TMAO_mM": e1,
            "experimental_fold": e1 / e0,
            "model_fermentation_per_glucose": m0,
            "model_TMAO_per_glucose": m1,
            "model_fold": m1 / m0 if abs(m0) > 1e-12 else np.nan,
        })
    products = pd.DataFrame(products)
    products.to_csv(OUT / "denby2015_product_redistribution.csv", index=False)
    rho = float(spearmanr(products["experimental_fold"], products["model_fold"]).statistic)
    summary = {
        **direct,
        "model_fermentation_net_e_per_glucose": float(m.loc["fermentation", "net_per_glucose"]),
        "model_TMAO_net_e_per_glucose": float(m.loc["tmao_matched", "net_per_glucose"]),
        "model_TMAO_expansion_pct": 100.0 * (float(m.loc["tmao_matched", "net_per_glucose"]) / float(m.loc["fermentation", "net_per_glucose"]) - 1.0),
        "model_TMAO_CO2_flux": float(m.loc["tmao_matched", "co2"]),
        "product_fold_rank_spearman_rho": rho,
    }
    (OUT / "denby2015_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary, model, products


def run_perrenoud(base, chem) -> pd.DataFrame:
    exp = pd.read_csv(DATA / "perrenoud2005" / "perrenoud2005_acceptor_routing_summary.csv")
    # Standardized model comparison only: equal uptake cap of 10, never matched to 40 mM medium concentration.
    states = {
        "strict anaerobic": solve_acceptor(base, chem, "oxygen", 0.0),
        "nitrate": solve_acceptor(base, chem, "nitrate", 10.0),
        "DMSO": solve_acceptor(base, chem, "dmso", 10.0),
    }
    rows = []
    for r in exp.itertuples(index=False):
        s = states[r.condition]
        rows.append({
            "experimental_condition": r.condition,
            "external_acceptor_mM": float(r.external_acceptor_mM),
            "reported_respiratory_cyclic_TCA_flux_mmol_g_h": float(r.reported_respiratory_cyclic_TCA_flux_mmol_g_h),
            "interpretation": r.interpretation,
            "model_standardized_uptake_cap": 0.0 if r.condition == "strict anaerobic" else 10.0,
            "model_actual_acceptor_uptake": float(s["actual"]),
            "model_net_e_per_glucose": float(s["net_per_glucose"]),
            "model_cumulative_e_per_glucose": float(s["cumulative_per_glucose"]),
            "model_CO2_flux": float(s["co2"]),
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "perrenoud2005_acceptor_comparison.csv", index=False)
    return out


def main() -> None:
    totals, details, sinks = run_toya()
    base = load_ecoli_model()
    chem = chemical_model(base)
    denby_summary, denby_model, denby_products = run_denby(base, chem)
    perrenoud = run_perrenoud(base, chem)
    publication_summary = {
        "toya": totals.to_dict(orient="records"),
        "denby": denby_summary,
        "perrenoud": perrenoud.to_dict(orient="records"),
    }
    (OUT / "alternative_acceptor_experimental_validation_summary.json").write_text(
        json.dumps(publication_summary, indent=2) + "\n"
    )
    print(totals.to_string(index=False))
    print(json.dumps(denby_summary, indent=2))
    print(perrenoud.to_string(index=False))


if __name__ == "__main__":
    main()
