#!/usr/bin/env python3
"""Reproduce the publication cross-species ETN and electron-flux analyses."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from cobra.flux_analysis import pfba
from cobra.io import read_sbml_model

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
DATA = ROOT / "data" / "cross_species"
OUT = ROOT / "results" / "cross_species"
OUT.mkdir(parents=True, exist_ok=True)

from etn.data_acquisition import load_ecoli_model
from etn.fba_case_study import solve_condition
from etn.cross_species import (
    apply_yeast_anaerobic,
    chemical_model,
    classify_and_qc,
    etn_metrics,
    load_cross_species_models,
    normalize_yeast_for_etn,
    relay_by_base,
    sha256,
    solve_growth_pfba,
)
from etn.matching import build_metabolite_table, build_matching_context, clear_matching_context_cache


def reaction_activity(chem_model, solution, reaction_id: str) -> float:
    if reaction_id not in chem_model.reactions or reaction_id not in solution.fluxes.index:
        return 0.0
    mt = build_metabolite_table(chem_model)
    ctx = build_matching_context(chem_model, mt)
    row = ctx.results_by_reaction[reaction_id]
    n_e = float(row.get("electrons_per_turnover_estimate") or 0.0)
    return abs(float(solution.fluxes[reaction_id])) * n_e


def append_metrics(phys, carriers, organism, condition, solution, chem, glucose, oxygen, objective):
    metrics, _, stats = etn_metrics(chem, solution, glucose)
    phys.append({
        "organism": organism,
        "condition": condition,
        "glucose_uptake": float(glucose),
        "oxygen_uptake": float(oxygen),
        "growth": float(solution.fluxes[objective]),
        **metrics,
    })
    for row in relay_by_base(stats):
        carriers.append({"organism": organism, "condition": condition, **row})
    return metrics




def write_publication_summary(phys_df: pd.DataFrame, carrier_df: pd.DataFrame) -> None:
    """Write paired cross-species compression and routing-response summary tables."""
    paired = phys_df[phys_df.organism.isin(["E. coli", "S. cerevisiae", "B. subtilis", "S. enterica", "K. phaffii"])].copy()
    low_names = {"E. coli": "anaerobic", "S. cerevisiae": "anaerobic", "B. subtilis": "anaerobic",
                 "S. enterica": "anoxic", "K. phaffii": "severe oxygen limitation"}
    high_names = {"E. coli": "aerobic", "S. cerevisiae": "aerobic", "B. subtilis": "aerobic",
                  "S. enterica": "high-respiratory-capacity", "K. phaffii": "high-respiratory-capacity"}
    rows = []
    for org in high_names:
        hi = paired[(paired.organism == org) & (paired.condition == high_names[org])].iloc[0]
        lo = paired[(paired.organism == org) & (paired.condition == low_names[org])].iloc[0]
        hi_nadh = carrier_df[(carrier_df.organism == org) & (carrier_df.condition == high_names[org]) & (carrier_df.carrier.str.lower() == "nadh")]
        lo_nadh = carrier_df[(carrier_df.organism == org) & (carrier_df.condition == low_names[org]) & (carrier_df.carrier.str.lower() == "nadh")]
        if len(hi_nadh) != 1 or len(lo_nadh) != 1:
            raise AssertionError(f"Expected one NADH-relay value per paired state for {org}")
        hi_nadh_value = float(hi_nadh.relay.iloc[0])
        lo_nadh_value = float(lo_nadh.relay.iloc[0])
        rows.append({
            "organism": org,
            "high_respiratory_capacity_net_e_per_glucose": hi.net_per_glucose,
            "low_oxygen_net_e_per_glucose": lo.net_per_glucose,
            "change_pct": 100 * (lo.net_per_glucose / hi.net_per_glucose - 1),
            "high_respiratory_capacity_cumulative_e_per_glucose": hi.cumulative_per_glucose,
            "low_oxygen_cumulative_e_per_glucose": lo.cumulative_per_glucose,
            "high_respiratory_capacity_transfer_depth": hi.effective_transfer_depth,
            "low_oxygen_transfer_depth": lo.effective_transfer_depth,
            "transfer_depth_change_pct": 100 * (lo.effective_transfer_depth / hi.effective_transfer_depth - 1),
            "high_respiratory_capacity_nadh_relay": hi_nadh_value,
            "low_oxygen_nadh_relay": lo_nadh_value,
            "nadh_relay_change_pct": 100 * (lo_nadh_value / hi_nadh_value - 1),
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "electron_flux_compression_summary.csv", index=False)
    low = summary.low_oxygen_net_e_per_glucose
    stats = {
        "low_oxygen_mean": float(low.mean()), "low_oxygen_sd": float(low.std(ddof=1)),
        "low_oxygen_cv_pct": float(100 * low.std(ddof=1) / low.mean()),
        "low_oxygen_min": float(low.min()), "low_oxygen_max": float(low.max()),
    }
    (OUT / "electron_flux_compression_stats.json").write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps(stats, indent=2))

def main():
    if "--summary-only" in sys.argv:
        phys_df = pd.read_csv(OUT / "cross_species_physiology.csv")
        carrier_df = pd.read_csv(OUT / "cross_species_carrier_relay.csv")
        write_publication_summary(phys_df, carrier_df)
        return

    ecoli = load_ecoli_model()
    others = load_cross_species_models(DATA)
    yeast = others["S. cerevisiae"]
    bac = others["B. subtilis"]
    salmonella = others["S. enterica"]
    pichia = others["K. phaffii"]
    ijo = others["iJO1366"]

    # -------- Static reconstruction/QC --------
    # Reconstruct each chemical ETN exactly once and reuse it for all flux states.
    clear_matching_context_cache()
    y_norm = normalize_yeast_for_etn(yeast)
    static = [
        ("E. coli", ecoli, chemical_model(ecoli)),
        ("S. cerevisiae", yeast, chemical_model(y_norm, yeast=True)),
        ("B. subtilis", bac, chemical_model(bac)),
        ("S. enterica", salmonella, chemical_model(salmonella)),
        ("K. phaffii", pichia, chemical_model(pichia)),
    ]
    static_chem = {organism: chem for organism, _full, chem in static}
    qcs, hubs = [], []
    qc_path, hubs_path = OUT / "cross_species_qc.csv", OUT / "cross_species_hubs.csv"
    if qc_path.exists() and hubs_path.exists() and "--force-qc" not in sys.argv:
        qcs = pd.read_csv(qc_path).to_dict("records")
    else:
        for organism, full, chem in static:
            qc, _, _, top = classify_and_qc(full, chem, organism)
            qcs.append(qc)
            for rank, (metabolite_id, name, degree) in enumerate(top, 1):
                hubs.append({
                    "organism": organism, "rank": rank, "metabolite_id": metabolite_id,
                    "name": name, "weighted_degree": degree,
                })
        pd.DataFrame(qcs).to_csv(qc_path, index=False)
        pd.DataFrame(hubs).to_csv(hubs_path, index=False)

    # -------- Dynamic states --------
    phys, carriers = [], []
    low_solutions = {}
    low_chems = {}

    echem = static_chem["E. coli"]
    for cond, qg, qo, zeros in [("aerobic", 8.7, 11.9, True), ("anaerobic", 14.9, 0.0, False)]:
        sol = solve_condition(ecoli, qg, qo, zeros)
        append_metrics(phys, carriers, "E. coli", cond, sol, echem, qg, qo,
                       "BIOMASS_Ec_iML1515_core_75p37M")
        if cond == "anaerobic":
            low_solutions["E. coli"] = sol; low_chems["E. coli"] = echem

    # Yeast aerobic
    ya = yeast.copy(); ya.objective = "r_2111"
    ya.reactions.get_by_id("r_1714").bounds = (-10, -10)
    ya.reactions.get_by_id("r_1992").lower_bound = -1000
    sol = solve_growth_pfba(ya, "r_2111")
    ychem_a = static_chem["S. cerevisiae"]
    append_metrics(phys, carriers, "S. cerevisiae", "aerobic", sol, ychem_a, 10,
                   max(0.0, -float(sol.fluxes["r_1992"])), "r_2111")

    # Yeast official anaerobic configuration
    yn = yeast.copy(); apply_yeast_anaerobic(yn, DATA / "aminoAcid_Bjorkeroth2020.tsv")
    yn.objective = "r_2111"; yn.reactions.get_by_id("r_1714").bounds = (-10, -10)
    sol = solve_growth_pfba(yn, "r_2111")
    ychem_n = static_chem["S. cerevisiae"]
    append_metrics(phys, carriers, "S. cerevisiae", "anaerobic", sol, ychem_n, 10, 0.0, "r_2111")
    low_solutions["S. cerevisiae"] = sol; low_chems["S. cerevisiae"] = ychem_n

    # B. subtilis, nitrate closed
    bchem = static_chem["B. subtilis"]
    for cond, o2bounds in [("aerobic", (-1000, 1000)), ("anaerobic", (0, 1000))]:
        bm = bac.copy(); bm.objective = "BiomassRepsolRed"
        bm.reactions.get_by_id("EX_glc__D_e").bounds = (-10, -10)
        bm.reactions.get_by_id("EX_nh4_e").bounds = (-10, 1000)
        bm.reactions.get_by_id("EX_no3_e").bounds = (0, 1000)
        bm.reactions.get_by_id("EX_o2_e").bounds = o2bounds
        sol = solve_growth_pfba(bm, "BiomassRepsolRed")
        append_metrics(phys, carriers, "B. subtilis", cond, sol, bchem, 10,
                       max(0.0, -float(sol.fluxes["EX_o2_e"])), "BiomassRepsolRed")
        if cond == "anaerobic":
            low_solutions["B. subtilis"] = sol; low_chems["B. subtilis"] = bchem

    # Salmonella, nitrate closed
    schem = static_chem["S. enterica"]
    salmonella_titration = []
    for cap in [20, 15, 10, 5, 2, 1, 0]:
        sm = salmonella.copy(); sm.objective = "BIOMASS_iRR1083_1"
        sm.reactions.get_by_id("EX_glc__D_e").bounds = (-10, -10)
        sm.reactions.get_by_id("EX_no3_e").bounds = (0, 1000)
        sm.reactions.get_by_id("EX_o2_e").bounds = (-cap, 1000) if cap else (0, 1000)
        sol = solve_growth_pfba(sm, "BIOMASS_iRR1083_1")
        metrics, _, _ = etn_metrics(schem, sol, 10)
        row = {
            "organism": "S. enterica", "oxygen_cap": cap,
            "growth": float(sol.fluxes["BIOMASS_iRR1083_1"]),
            "actual_oxygen_uptake": max(0.0, -float(sol.fluxes["EX_o2_e"])), **metrics,
        }
        salmonella_titration.append(row)
        if cap in (20, 0):
            cond = "high-respiratory-capacity" if cap == 20 else "anoxic"
            append_metrics(phys, carriers, "S. enterica", cond, sol, schem, 10,
                           row["actual_oxygen_uptake"], "BIOMASS_iRR1083_1")
            if cap == 0:
                low_solutions["S. enterica"] = sol; low_chems["S. enterica"] = schem
    pd.DataFrame(salmonella_titration).to_csv(OUT / "salmonella_oxygen_titration.csv", index=False)

    # K. phaffii oxygen titration; 0.5 is the physiological low-oxygen endpoint.
    pchem = static_chem["K. phaffii"]
    pichia_titration = []
    for cap in [20, 10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0]:
        pm = pichia.copy(); pm.objective = "Ex_biomass"
        pm.reactions.get_by_id("Ex_glc_D").bounds = (-10, -10)
        pm.reactions.get_by_id("Ex_o2").bounds = (-cap, 0) if cap else (0, 0)
        sol = solve_growth_pfba(pm, "Ex_biomass")
        metrics, _, _ = etn_metrics(pchem, sol, 10)
        row = {
            "organism": "K. phaffii", "oxygen_cap": cap,
            "growth": float(sol.fluxes["Ex_biomass"]),
            "actual_oxygen_uptake": max(0.0, -float(sol.fluxes["Ex_o2"])), **metrics,
        }
        pichia_titration.append(row)
        if cap in (20, 0.5):
            cond = "high-respiratory-capacity" if cap == 20 else "severe oxygen limitation"
            append_metrics(phys, carriers, "K. phaffii", cond, sol, pchem, 10,
                           row["actual_oxygen_uptake"], "Ex_biomass")
            if cap == 0.5:
                low_solutions["K. phaffii"] = sol; low_chems["K. phaffii"] = pchem
    pd.DataFrame(pichia_titration).to_csv(OUT / "pichia_oxygen_titration.csv", index=False)


    phys_df = pd.DataFrame(phys)
    phys_df.to_csv(OUT / "cross_species_physiology.csv", index=False)
    carrier_df = pd.DataFrame(carriers)
    carrier_df.to_csv(OUT / "cross_species_carrier_relay.csv", index=False)

    # -------- Central-carbon source attribution in low-O2/anoxic states --------
    source_rows = []
    source_rxns = {
        "E. coli": ["GAPD"],
        "S. enterica": ["GAPD"],
        "B. subtilis": ["GAPD"],
        "S. cerevisiae": ["r_0486", "r_0466", "r_0889"],  # GAPD + two oxidative PPP steps
    }
    for organism, rxns in source_rxns.items():
        sol, chem = low_solutions[organism], low_chems[organism]
        selected_activity = sum(reaction_activity(chem, sol, rid) for rid in rxns)
        row = phys_df[(phys_df.organism == organism) & phys_df.condition.str.contains("anaerobic|anoxic", regex=True)].iloc[0]
        source_rows.append({
            "organism": organism,
            "reaction_ids": ";".join(rxns),
            "selected_source_e_flux": selected_activity,
            "selected_source_e_per_glucose": selected_activity / float(row.glucose_uptake),
            "total_net_e_flux": float(row.net_e_flux),
            "fraction_of_net_source_flux": selected_activity / float(row.net_e_flux),
        })
    pd.DataFrame(source_rows).to_csv(OUT / "low_oxygen_source_attribution.csv", index=False)

    # -------- Independent E. coli model robustness --------
    robust = []
    # standardised iML1515
    for model_name, base, objective in [
        ("iML1515", ecoli, "BIOMASS_Ec_iML1515_core_75p37M"),
        ("iJO1366", ijo, "BIOMASS_Ec_iJO1366_core_53p95M"),
    ]:
        chem = chemical_model(base)
        for cond, cap in [("aerobic", 20), ("anaerobic", 0)]:
            mm = base.copy(); mm.objective = objective
            mm.reactions.get_by_id("EX_glc__D_e").bounds = (-10, -10)
            mm.reactions.get_by_id("EX_o2_e").bounds = (-cap, 1000) if cap else (0, 1000)
            # For a comparable fermentative test, do not impose E. coli manuscript product non-detections here.
            sol = solve_growth_pfba(mm, objective)
            metrics, _, _ = etn_metrics(chem, sol, 10)
            robust.append({
                "model": model_name, "condition": cond,
                "growth": float(sol.fluxes[objective]),
                "oxygen_uptake": max(0.0, -float(sol.fluxes["EX_o2_e"])), **metrics,
            })
    pd.DataFrame(robust).to_csv(OUT / "ecoli_model_robustness.csv", index=False)

    # -------- Input provenance/checksums --------
    meta = pd.read_csv(DATA / "model_metadata.tsv", sep="\t")
    checksums = {}
    for line in (DATA / "checksums.sha256").read_text().splitlines():
        if line.strip():
            digest, name = line.split(None, 1); checksums[name.strip()] = digest
    meta["sha256"] = meta["bundled_file"].map(lambda x: checksums.get(Path(str(x)).name, sha256(ROOT / "data" / "iML1515.xml") if str(x) == "../iML1515.xml" else ""))
    meta.to_csv(OUT / "cross_species_model_metadata.csv", index=False)

    # Publication-facing convergence and routing-response summary.
    write_publication_summary(phys_df, carrier_df)


if __name__ == "__main__":
    main()
