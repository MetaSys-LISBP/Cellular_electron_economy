"""Experimentally constrained aerobic-versus-anaerobic glucose case study.

The physiological example uses wild-type *Escherichia coli* K-12 MG1655
cultivated in M9 minimal medium with glucose as the sole carbon source under
aerobic or fermentative (anaerobic) conditions, as characterised by Chen et
al. (Metab. Eng. 13, 38-48, 2011; DOI: 10.1016/j.ymben.2010.11.004).

The experimentally measured glucose uptake rates are 8.7 and 14.9 mmol
gDW^-1 h^-1 in aerobic and anaerobic cultures, respectively. Oxygen uptake is
11.9 mmol gDW^-1 h^-1 aerobically and zero anaerobically. These exact values
are also reproduced in a later pFBA application of the Chen dataset (Jamialahmadi et al., *Molecular BioSystems* 12, 3459–3466, 2016;
DOI: 10.1039/C6MB00532B).

For the aerobic state, Chen et al. report that lactate, succinate, formate and
ethanol were not detected. Their exchange reactions are therefore fixed to
zero. Acetate secretion is deliberately *not* constrained: Chen et al. showed
that genome-scale FBA constrained by measured glucose and oxygen uptake could
predict aerobic product secretion. For the anaerobic state, no organic
product exchange is fixed, allowing pFBA to determine the fermentative product
mixture consistent with glucose uptake and the absence of oxygen.

For each condition, biomass growth is first maximised, then parsimonious FBA
(pFBA) selects the minimum-total-flux solution at 100% of the growth optimum.
The resulting flux state is therefore a model-based completion of the measured
external constraints, not a direct measurement of every intracellular flux.
"""

import os
from typing import Dict, List
import pandas as pd

from .matching import classify_reaction_molecular

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DEFAULT_DATA_PATH = os.path.join(REPO_ROOT, "data", "aerobic_anaerobic_glucose_constraints.tsv")

BIOMASS_REACTION_ID = "BIOMASS_Ec_iML1515_core_75p37M"
GLUCOSE_EXCHANGE_ID = "EX_glc__D_e"
OXYGEN_EXCHANGE_ID = "EX_o2_e"
ATP_MAINTENANCE_ID = "ATPM"

# Chen et al. explicitly report no detectable lactate, succinate, formate or
# ethanol in the aerobic culture. Both lactate stereoisomer exchanges present
# in iML1515 are constrained for completeness.
AEROBIC_ZERO_PRODUCT_EXCHANGES = (
    "EX_lac__D_e", "EX_lac__L_e", "EX_succ_e", "EX_for_e", "EX_etoh_e"
)


def load_experimental_data(path: str = None) -> pd.DataFrame:
    """Load the two experimentally constrained physiological states."""
    df = pd.read_csv(path or DEFAULT_DATA_PATH, sep="\t")
    required = {"condition", "qGLC", "qO2", "aerobic_zero_fermentation_products"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    df["condition"] = df["condition"].str.lower()
    if set(df["condition"]) != {"aerobic", "anaerobic"}:
        raise ValueError("Expected exactly the conditions 'aerobic' and 'anaerobic'.")
    return df


def configure_condition_model(model, qGLC: float, qO2: float,
                              aerobic_zero_fermentation_products: bool = False):
    """Return a model copy carrying exactly the published condition constraints.

    This helper is shared by pFBA and flux-space robustness analyses so that
    alternative-optimum calculations cannot silently diverge from the primary
    physiological workflow.
    """
    m = model.copy()
    m.reactions.get_by_id(GLUCOSE_EXCHANGE_ID).bounds = (-float(qGLC), -float(qGLC))
    m.reactions.get_by_id(OXYGEN_EXCHANGE_ID).bounds = (-float(qO2), -float(qO2))
    if bool(aerobic_zero_fermentation_products):
        for rid in AEROBIC_ZERO_PRODUCT_EXCHANGES:
            m.reactions.get_by_id(rid).bounds = (0.0, 0.0)
    m.objective = BIOMASS_REACTION_ID
    return m


def solve_condition(model, qGLC: float, qO2: float,
                    aerobic_zero_fermentation_products: bool = False):
    """Return a growth-maximising pFBA solution for one measured state."""
    from cobra.flux_analysis import pfba

    m = configure_condition_model(
        model, qGLC, qO2, aerobic_zero_fermentation_products
    )
    primary = m.optimize()
    if primary.status != "optimal":
        return None
    sol = pfba(m, fraction_of_optimum=1.0)
    return sol if sol.status == "optimal" else None


def solve_row(model, row):
    return solve_condition(
        model, row.qGLC, row.qO2,
        bool(int(row.aerobic_zero_fermentation_products)),
    )


def electron_flux_table(model, met_table, solution, min_flux: float = 1e-6) -> pd.DataFrame:
    """Cumulative redox activity for all active redox reactions in one state."""
    rows = []
    for rxn in model.reactions:
        flux = solution.fluxes.get(rxn.id, 0.0)
        if abs(flux) < min_flux:
            continue
        res = classify_reaction_molecular(rxn, met_table)
        if not res["is_redox"]:
            continue
        electrons_per_turnover = res["electrons_per_turnover_estimate"]
        rows.append({
            "reaction_id": rxn.id,
            "reaction_name": rxn.name,
            "flux": flux,
            "resolution_status": res["resolution_status"],
            "electrons_per_turnover": electrons_per_turnover,
            "minimum_unresolved_donation": res["minimum_unresolved_donation"],
            "minimum_unresolved_acceptance": res["minimum_unresolved_acceptance"],
            "electron_flux": abs(flux) * electrons_per_turnover,
            "subsystem": rxn.subsystem,
        })
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values("electron_flux", ascending=False).reset_index(drop=True)
    return df


def run_case_study(model, met_table, data: pd.DataFrame) -> Dict[str, object]:
    """Solve both states and summarize cumulative activity by subsystem."""
    state_rows: List[dict] = []
    subsystem_series = {}
    reaction_tables = {}

    for i, row in data.iterrows():
        solution = solve_row(model, row)
        if solution is None:
            state_rows.append({"state": i, "condition": row.condition, "status": "infeasible"})
            continue
        et = electron_flux_table(model, met_table, solution)
        reaction_tables[row.condition] = et
        subsystem_series[row.condition] = et.groupby("subsystem")["electron_flux"].sum()
        state_rows.append({
            "state": i,
            "condition": row.condition,
            "status": "optimal",
            "qGLC": float(row.qGLC),
            "qO2": float(row.qO2),
            "predicted_growth_rate": float(solution.fluxes[BIOMASS_REACTION_ID]),
            "ATPM_flux": float(solution.fluxes[ATP_MAINTENANCE_ID]),
            "predicted_acetate_exchange": float(solution.fluxes["EX_ac_e"]),
            "predicted_formate_exchange": float(solution.fluxes["EX_for_e"]),
            "predicted_ethanol_exchange": float(solution.fluxes["EX_etoh_e"]),
            "predicted_succinate_exchange": float(solution.fluxes["EX_succ_e"]),
            "total_cumulative_electron_transfer_activity": float(et["electron_flux"].sum()),
        })

    summary = pd.DataFrame(state_rows)
    if (summary.status != "optimal").any():
        raise RuntimeError("At least one physiological state is infeasible.")

    all_subsystems = sorted(set().union(*(s.index for s in subsystem_series.values())))
    sub = pd.DataFrame(index=all_subsystems)
    for cond in ("aerobic", "anaerobic"):
        sub[cond] = subsystem_series[cond].reindex(all_subsystems).fillna(0.0)
    sub["delta_anaerobic_minus_aerobic"] = sub["anaerobic"] - sub["aerobic"]
    sub = sub.sort_values("delta_anaerobic_minus_aerobic", key=abs, ascending=False)

    return {
        "state_summary": summary,
        "subsystem_comparison": sub,
        "reaction_tables": reaction_tables,
    }

