"""Flux-variability analysis for the aerobic/anaerobic physiological case study.

The primary paper reports one parsimonious (pFBA) completion of experimentally
measured glucose/oxygen constraints.  This script asks which conclusions are
fixed across all growth-optimal FBA solutions.  It uses the exact same condition
configuration as ``etn.fba_case_study`` and performs FVA at 100% of maximal
growth.
"""
import argparse
import os
import sys

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)

import pandas as pd
from cobra.flux_analysis import flux_variability_analysis
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from etn.data_acquisition import load_ecoli_model
from etn.fba_case_study import (
    BIOMASS_REACTION_ID,
    configure_condition_model,
    load_experimental_data,
    solve_row,
)

RESULTS_DIR = os.path.join(REPO_ROOT, "results")

# Reactions needed to test the main physiological interpretation.  The table
# deliberately mixes measured/constrained exchanges, pFBA-predicted product
# exchanges, fermentative reductions, and respiratory/transhydrogenase fluxes.
KEY_REACTIONS = [
    BIOMASS_REACTION_ID,
    "EX_glc__D_e", "EX_o2_e", "EX_ac_e", "EX_for_e", "EX_etoh_e", "EX_succ_e",
    "ALCD2x", "ACALD", "THD2pp", "NADH16pp", "CYTBO3_4pp", "CYTBDpp", "CYTBD2pp",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--all-reactions", action="store_true",
                        help="also export full-model FVA tables (slower)")
    args = parser.parse_args()
    os.makedirs(RESULTS_DIR, exist_ok=True)
    model = load_ecoli_model()
    data = load_experimental_data()
    rows = []

    for r in data.itertuples(index=False):
        m = configure_condition_model(
            model, float(r.qGLC), float(r.qO2),
            bool(int(r.aerobic_zero_fermentation_products)),
        )
        optimum = m.slim_optimize(error_value=None)
        if optimum is None:
            raise RuntimeError(f"{r.condition}: infeasible growth optimization")

        # FVA of the reactions required for the principal robustness claims is
        # fast enough for the default publication workflow.  Full-model FVA is
        # available explicitly with --all-reactions.
        fva = flux_variability_analysis(
            m, reaction_list=[m.reactions.get_by_id(rid) for rid in KEY_REACTIONS],
            fraction_of_optimum=1.0, processes=1
        )
        fva.index.name = "reaction_id"
        fva = fva.reset_index()
        if args.all_reactions:
            fva_all = flux_variability_analysis(
                m, reaction_list=m.reactions, fraction_of_optimum=1.0, processes=1
            )
            fva_all.index.name = "reaction_id"
            fva_all = fva_all.reset_index()
            fva_all.insert(0, "condition", r.condition)
            fva_all.insert(1, "maximal_growth_rate", float(optimum))
            fva_all.to_csv(os.path.join(RESULTS_DIR, f"fva_all_reactions_{r.condition}.csv"), index=False)

        pfba_sol = solve_row(model, r)
        if pfba_sol is None:
            raise RuntimeError(f"{r.condition}: pFBA failed")

        for rid in KEY_REACTIONS:
            if rid not in fva.reaction_id.values:
                continue
            z = fva.loc[fva.reaction_id.eq(rid)].iloc[0]
            rxn = m.reactions.get_by_id(rid)
            tol = 1e-7
            p = float(pfba_sol.fluxes[rid])
            lo, hi = float(z.minimum), float(z.maximum)
            # Normalize solver-scale numerical noise around zero for the
            # publication-facing table while retaining a conservative tolerance.
            if abs(p) < tol: p = 0.0
            if abs(lo) < tol: lo = 0.0
            if abs(hi) < tol: hi = 0.0
            span = max(0.0, hi - lo)
            fixed = span <= tol
            sign_invariant = (lo > tol) or (hi < -tol) or (abs(lo) <= tol and abs(hi) <= tol)
            rows.append({
                "condition": r.condition,
                "reaction_id": rid,
                "reaction_name": rxn.name,
                "pFBA_flux": p,
                "FVA_min": lo,
                "FVA_max": hi,
                "FVA_span": span,
                "flux_fixed_at_growth_optimum": bool(fixed),
                "sign_or_zero_invariant": bool(sign_invariant),
            })

    key = pd.DataFrame(rows)
    key.to_csv(os.path.join(RESULTS_DIR, "fva_key_reactions_by_condition.csv"), index=False)

    # Compact flux-space robustness summary.  This does not assert
    # invariance of every electron-transfer edge: in particular NADH16pp is
    # allowed to vary aerobically across alternative optimal solutions.
    def row(cond, rid):
        return key[(key.condition == cond) & (key.reaction_id == rid)].iloc[0]

    summary = pd.DataFrame([
        {
            "finding": "terminal_oxygen_use",
            "aerobic_min": row("aerobic", "CYTBO3_4pp").FVA_min,
            "aerobic_max": row("aerobic", "CYTBO3_4pp").FVA_max,
            "anaerobic_min": row("anaerobic", "CYTBO3_4pp").FVA_min,
            "anaerobic_max": row("anaerobic", "CYTBO3_4pp").FVA_max,
            "interpretation": "Bo3 terminal oxidase use is fixed aerobically and zero anaerobically at maximal growth under these constraints.",
        },
        {
            "finding": "ethanol_branch_ALCD2x",
            "aerobic_min": row("aerobic", "ALCD2x").FVA_min,
            "aerobic_max": row("aerobic", "ALCD2x").FVA_max,
            "anaerobic_min": row("anaerobic", "ALCD2x").FVA_min,
            "anaerobic_max": row("anaerobic", "ALCD2x").FVA_max,
            "interpretation": "The ethanol-forming reduction is zero aerobically and fixed in the anaerobic growth optimum.",
        },
        {
            "finding": "respiratory_NADH16pp_allocation",
            "aerobic_min": row("aerobic", "NADH16pp").FVA_min,
            "aerobic_max": row("aerobic", "NADH16pp").FVA_max,
            "anaerobic_min": row("anaerobic", "NADH16pp").FVA_min,
            "anaerobic_max": row("anaerobic", "NADH16pp").FVA_max,
            "interpretation": "Exact aerobic NADH-dehydrogenase allocation is not unique; local pFBA edge magnitudes should not be interpreted as experimentally fixed.",
        },
    ])
    summary.to_csv(os.path.join(RESULTS_DIR, "fva_robustness_summary.csv"), index=False)
    print(key.to_string(index=False))
    print("\nCompact robustness summary:\n", summary.to_string(index=False))


if __name__ == "__main__":
    main()
