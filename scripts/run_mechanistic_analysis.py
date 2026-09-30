#!/usr/bin/env python3
"""Reproduce the publication-facing terminal-acceptor and carbon-substrate analyses.

This script tests two mechanistic questions introduced in the publication-facing
manuscript:

1. whether terminal electron acceptors permit expansion of net electron flux;
2. whether the low-electron-flux state depends on carbon-entry chemistry.

All flux states are growth-maximising pFBA solutions of iML1515.  Electron
electron fluxes are calculated with the same chemistry-driven ETN machinery used by the
primary publication workflow.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUT = ROOT / "results" / "publication"
OUT.mkdir(parents=True, exist_ok=True)

from etn.data_acquisition import load_ecoli_model
from etn.cross_species import chemical_model, etn_metrics, solve_growth_pfba
from etn.degree_of_reduction import degree_of_reduction
from etn.matching import build_metabolite_table, build_matching_context

BIOMASS = "BIOMASS_Ec_iML1515_core_75p37M"
GLUCOSE = "EX_glc__D_e"
ACCEPTOR_EXCHANGES = {
    "oxygen": "EX_o2_e",
    "nitrate": "EX_no3_e",
    "tmao": "EX_tmao_e",
    "dmso": "EX_dmso_e",
    "fumarate": "EX_fum_e",
}

# Main-text Fig. 3a uses a shared electron-equivalent capacity axis.
# Values are nominal electron-accepting capacities in e- mmol gDW^-1 h^-1.
# Every point is solved explicitly (no interpolation).
FIG3A_ELECTRON_CAPACITY_GRID = [0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0]
ACCEPTOR_ELECTRONS_PER_MMOL = {
    "oxygen": 4.0,
    "nitrate": 2.0,
    "tmao": 2.0,
}
SUBSTRATE_EXCHANGES = {
    "glucose": "EX_glc__D_e",
    "fructose": "EX_fru_e",
    "xylose": "EX_xyl__D_e",
    "galactose": "EX_gal_e",
    "gluconate": "EX_glcn_e",
    "glycerol": "EX_glyc_e",
}
CORE_SOURCE_REACTIONS = {
    "glucose": ["GAPD"],
    "fructose": ["GAPD"],
    "xylose": ["GAPD"],
    "galactose": ["GAPD"],
    "gluconate": ["GAPD", "GND"],
    "glycerol": ["G3PD2", "GAPD"],
}


def close_tested_acceptors(model) -> None:
    for rid in ACCEPTOR_EXCHANGES.values():
        model.reactions.get_by_id(rid).bounds = (0.0, 1000.0)


def close_tested_substrates(model) -> None:
    for rid in SUBSTRATE_EXCHANGES.values():
        model.reactions.get_by_id(rid).bounds = (0.0, 1000.0)


def solve_acceptor(base, chem, acceptor: str, cap: float):
    model = base.copy()
    model.objective = BIOMASS
    model.reactions.get_by_id(GLUCOSE).bounds = (-10.0, -10.0)
    close_tested_acceptors(model)
    rid = ACCEPTOR_EXCHANGES[acceptor]
    model.reactions.get_by_id(rid).bounds = (-float(cap), 1000.0) if cap else (0.0, 1000.0)
    sol = solve_growth_pfba(model, BIOMASS)
    metrics, _, _ = etn_metrics(chem, sol, 10.0)
    return {
        "acceptor": acceptor,
        "cap": float(cap),
        "actual": max(0.0, -float(sol.fluxes[rid])),
        "growth": float(sol.fluxes[BIOMASS]),
        "co2": float(sol.fluxes["EX_co2_e"]),
        **metrics,
    }


def reaction_electron_activity(chem, solution, reaction_id: str) -> tuple[float, float, float]:
    """Return abs flux, e/turnover and cumulative electron activity for one reaction."""
    mt = build_metabolite_table(chem)
    ctx = build_matching_context(chem, mt)
    result = ctx.results_by_reaction[reaction_id]
    ne = float(result.get("electrons_per_turnover_estimate") or 0.0)
    flux = abs(float(solution.fluxes[reaction_id]))
    return flux, ne, flux * ne


def run_acceptor_panel(base, chem):
    caps = {
        "oxygen": [0, 2, 5, 10, 20],
        "nitrate": [0, 5, 10, 20, 50],
        "tmao": [0, 10, 20, 50, 100],
    }
    rows = [solve_acceptor(base, chem, acc, cap) for acc, values in caps.items() for cap in values]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "acceptor_titration.csv", index=False)

    endpoint_specs = [
        ("fermentation", "oxygen", 0),
        ("oxygen", "oxygen", 20),
        ("nitrate", "nitrate", 10),
        ("tmao", "tmao", 10),
        ("dmso", "dmso", 10),
        ("fumarate", "fumarate", 10),
    ]
    endpoints = []
    for label, acc, cap in endpoint_specs:
        row = solve_acceptor(base, chem, acc, cap)
        endpoints.append({
            "acceptor": label,
            "cap": cap,
            "actual_uptake": row["actual"],
            "growth": row["growth"],
            "cumulative_e_flux": row["cumulative_e_flux"],
            "net_e_flux": row["net_e_flux"],
            "effective_transfer_depth": row["effective_transfer_depth"],
            "cumulative_per_glucose": row["cumulative_per_glucose"],
            "net_per_glucose": row["net_per_glucose"],
        })
    pd.DataFrame(endpoints).to_csv(OUT / "acceptor_endpoints.csv", index=False)

    r, p = pearsonr(df["net_per_glucose"], df["co2"])
    rho, ps = spearmanr(df["net_per_glucose"], df["co2"])
    rc, pc = pearsonr(df["cumulative_per_glucose"], df["co2"])
    fit = np.polyfit(df["co2"], df["net_per_glucose"], 1)
    pred = np.polyval(fit, df["co2"])
    r2 = 1.0 - float(((df["net_per_glucose"] - pred) ** 2).sum()) / float(((df["net_per_glucose"] - df["net_per_glucose"].mean()) ** 2).sum())
    stats = {
        "net_electron_flux_vs_CO2_pearson_r": float(r),
        "net_electron_flux_vs_CO2_pearson_p": float(p),
        "net_electron_flux_vs_CO2_spearman_rho": float(rho),
        "net_electron_flux_vs_CO2_spearman_p": float(ps),
        "net_electron_flux_vs_CO2_linear_R2": float(r2),
        "cumulative_electron_flux_vs_CO2_pearson_r": float(rc),
        "cumulative_electron_flux_vs_CO2_pearson_p": float(pc),
    }
    (OUT / "acceptor_carbon_oxidation_statistics.json").write_text(json.dumps(stats, indent=2) + "\n")
    return df, stats


def run_fig3_acceptor_electron_capacity_panel(base, chem):
    """Run the simulations used in main-text Fig. 3a-b.

    Oxygen, nitrate and TMAO are compared on the same nominal
    electron-accepting-capacity grid. The electron-equivalent x values are
    converted to acceptor-specific molar uptake caps before each pFBA solve:
    4 e-/mmol for O2 and 2 e-/mmol for nitrate or TMAO.

    The approved grid is 0, 5, 10, 20, 30, 40, 50, 60 and 80
    e- mmol gDW^-1 h^-1. No point is interpolated.
    """
    rows = []
    for acceptor in ["oxygen", "nitrate", "tmao"]:
        e_per_mmol = ACCEPTOR_ELECTRONS_PER_MMOL[acceptor]
        for nominal_e_capacity in FIG3A_ELECTRON_CAPACITY_GRID:
            molar_cap = float(nominal_e_capacity) / e_per_mmol
            solved = solve_acceptor(base, chem, acceptor, molar_cap)
            actual_uptake = float(solved["actual"])
            rows.append({
                "acceptor": acceptor,
                "nominal_e_capacity": float(nominal_e_capacity),
                "molar_cap": molar_cap,
                "actual_uptake": actual_uptake,
                "actual_e_capacity": actual_uptake * e_per_mmol,
                "net_per_glucose": float(solved["net_per_glucose"]),
                "co2": float(solved["co2"]),
                "growth": float(solved["growth"]),
            })

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "acceptor_electron_capacity_titration.csv", index=False)

    r, p = pearsonr(df["net_per_glucose"], df["co2"])
    rho, ps = spearmanr(df["net_per_glucose"], df["co2"])
    fit = np.polyfit(df["co2"], df["net_per_glucose"], 1)
    pred = np.polyval(fit, df["co2"])
    r2 = 1.0 - float(((df["net_per_glucose"] - pred) ** 2).sum()) / float(((df["net_per_glucose"] - df["net_per_glucose"].mean()) ** 2).sum())
    stats = {
        "n_exact_simulations": int(len(df)),
        "electron_capacity_grid": FIG3A_ELECTRON_CAPACITY_GRID,
        "acceptor_electrons_per_mmol": ACCEPTOR_ELECTRONS_PER_MMOL,
        "net_electron_flux_vs_CO2_pearson_r": float(r),
        "net_electron_flux_vs_CO2_pearson_p": float(p),
        "net_electron_flux_vs_CO2_spearman_rho": float(rho),
        "net_electron_flux_vs_CO2_spearman_p": float(ps),
        "net_electron_flux_vs_CO2_linear_R2": float(r2),
    }
    (OUT / "acceptor_electron_capacity_statistics.json").write_text(json.dumps(stats, indent=2) + "\n")
    return df, stats


def run_substrate_panel(base, chem):
    rows = []
    anoxic_solutions = {}
    for substrate, exchange in SUBSTRATE_EXCHANGES.items():
        for state, oxygen_cap in [("anoxic", 0.0), ("high-respiratory-capacity", 20.0)]:
            model = base.copy()
            model.objective = BIOMASS
            close_tested_substrates(model)
            close_tested_acceptors(model)
            model.reactions.get_by_id(exchange).bounds = (-10.0, -10.0)
            model.reactions.get_by_id("EX_o2_e").bounds = (-oxygen_cap, 1000.0) if oxygen_cap else (0.0, 1000.0)
            sol = solve_growth_pfba(model, BIOMASS)
            metrics, _, _ = etn_metrics(chem, sol, 10.0)

            ex_rxn = model.reactions.get_by_id(exchange)
            extracellular = next(m for m in ex_rxn.metabolites if m.compartment == "e")
            formula = extracellular.formula
            gamma = degree_of_reduction(extracellular.elements, extracellular.charge)
            carbon = int(extracellular.elements.get("C", 0))
            rows.append({
                "substrate": substrate,
                "exchange": exchange,
                "formula": formula,
                "charge": extracellular.charge,
                "C": carbon,
                "gamma": gamma,
                "state": state,
                "growth": float(sol.fluxes[BIOMASS]),
                "actual_o2": max(0.0, -float(sol.fluxes["EX_o2_e"])),
                **metrics,
                "net_per_C": metrics["net_per_glucose"] / carbon,
                "net_fraction_substrate_gamma": metrics["net_per_glucose"] / gamma,
                "cumulative_per_C": metrics["cumulative_per_glucose"] / carbon,
            })
            if state == "anoxic":
                anoxic_solutions[substrate] = sol

    raw = pd.DataFrame(rows)
    raw.to_csv(OUT / "substrate_flux_states.csv", index=False)

    anox = raw[raw.state == "anoxic"].copy()
    oxy = raw[raw.state == "high-respiratory-capacity"][["substrate", "net_per_glucose", "cumulative_per_glucose", "effective_transfer_depth"]].copy()
    summary = anox.merge(oxy, on="substrate", suffixes=("_anoxic", "_oxygen"))
    summary["compression_pct"] = 100.0 * (1.0 - summary["net_per_glucose_anoxic"] / summary["net_per_glucose_oxygen"])
    summary["expansion_factor"] = summary["net_per_glucose_oxygen"] / summary["net_per_glucose_anoxic"]
    summary.to_csv(OUT / "substrate_flux_compression.csv", index=False)

    source_rows = []
    for substrate, reaction_ids in CORE_SOURCE_REACTIONS.items():
        sol = anoxic_solutions[substrate]
        details, total = [], 0.0
        for rid in reaction_ids:
            flux, ne, activity = reaction_electron_activity(chem, sol, rid)
            total += activity
            details.append(f"{rid}:{flux:.4f}x{ne:g}e={activity:.4f}")
        total_net = float(summary.loc[summary.substrate == substrate, "net_e_flux"].iloc[0])
        source_rows.append({
            "substrate": substrate,
            "selected_reactions": ";".join(reaction_ids),
            "details": " | ".join(details),
            "selected_e_flux": total,
            "selected_e_per_substrate": total / 10.0,
            "net_e_flux": total_net,
            "net_e_per_substrate": total_net / 10.0,
            "fraction_of_net": total / total_net,
        })
    source = pd.DataFrame(source_rows)
    source["pred_error"] = source["selected_e_per_substrate"] - source["net_e_per_substrate"]
    source["abs_error"] = source["pred_error"].abs()
    source.to_csv(OUT / "substrate_core_redox_contribution.csv", index=False)

    r, p = pearsonr(source["selected_e_per_substrate"], source["net_e_per_substrate"])
    stats = {
        "core_vs_total_pearson_r": float(r),
        "core_vs_total_pearson_p": float(p),
        "mean_fraction_of_net": float(source["fraction_of_net"].mean()),
        "min_fraction_of_net": float(source["fraction_of_net"].min()),
        "max_fraction_of_net": float(source["fraction_of_net"].max()),
        "mean_absolute_error_e_per_substrate": float(source["abs_error"].mean()),
    }
    (OUT / "substrate_core_redox_statistics.json").write_text(json.dumps(stats, indent=2) + "\n")
    return raw, summary, source, stats


def main():
    model = load_ecoli_model()
    chem = chemical_model(model)
    # Build once; the matching context is cached for all subsequent flux states.
    mt = build_metabolite_table(chem)
    build_matching_context(chem, mt)
    _, acceptor_stats = run_acceptor_panel(model, chem)
    _, fig3_acceptor_stats = run_fig3_acceptor_electron_capacity_panel(model, chem)
    _, _, _, substrate_stats = run_substrate_panel(model, chem)
    print(json.dumps({
        "acceptor_molar_titration": acceptor_stats,
        "fig3_acceptor_electron_capacity": fig3_acceptor_stats,
        "substrate": substrate_stats,
    }, indent=2))


if __name__ == "__main__":
    main()
