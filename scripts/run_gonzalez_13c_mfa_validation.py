#!/usr/bin/env python3
"""Calculate electron fluxes from Gonzalez et al. (2017) integrated 13C-MFA flux distributions.

The input tables in ``data/experimental/gonzalez2017`` were extracted from the
supplementary material of Gonzalez, Long & Antoniewicz, Metabolic Engineering
39, 9–18 (2017), DOI: 10.1016/j.ymben.2016.11.003.

The 13C-MFA model explicitly contains NADH, NADPH and FADH2.  We therefore
calculate conserved electron flux directly from the fitted reaction fluxes:

- two electrons per NADH, NADPH or FADH2 equivalent;
- carrier-to-carrier conversions are relay rather than new electron source;
- direct formate -> CO2 + H2 is a two-electron transfer;
- the exact condition-specific biomass NAD(P)H stoichiometry is included.

Fluxes in the source workbook are normalized to 100 mol substrate, so derived
values are divided by 100 to report electron equivalents per substrate.
"""
from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "gonzalez2017"
OUT = ROOT / "results" / "publication"
OUT.mkdir(parents=True, exist_ok=True)

REDUCED = ("NADH", "NADPH", "FADH2")
SUBSTRATE_CHEMISTRY = {
    # C atoms and total degree of reduction (electron equivalents per molecule)
    # for the neutral substrates used in Gonzalez et al.
    "glucose": {"C": 6.0, "gamma": 24.0},
    "xylose": {"C": 5.0, "gamma": 20.0},
}

CONDITION_KEYS = {
    "Aerobic glucose": "aerobic_glucose",
    "Anaerobic glucose": "anaerobic_glucose",
    "Aerobic xylose": "aerobic_xylose",
    "Anaerobic xylose": "anaerobic_xylose",
}


def split_equation(equation: str) -> tuple[str, str]:
    for arrow in ("<=>", "->", "=>"):
        if arrow in equation:
            left, right = equation.split(arrow, 1)
            return left.strip(), right.strip()
    raise ValueError(f"No reaction arrow in {equation!r}")


def species_coefficient(side: str, species: str) -> float:
    """Return total coefficient of a named species in a simple reaction side."""
    total = 0.0
    for raw in side.split("+"):
        term = raw.strip()
        # Remove atom-mapping parentheses present in the model workbook if any.
        term = re.sub(r"\s*\([^)]*\)\s*$", "", term).strip()
        match = re.match(r"^(?:(\d+(?:\.\d+)?)\s+)?(.+)$", term)
        if not match:
            continue
        coef = float(match.group(1) or 1.0)
        name = match.group(2).strip()
        if name == species:
            total += coef
    return total


def reduced_counts(equation: str) -> tuple[float, float]:
    """Reduced-cofactor equivalents consumed and produced in written direction."""
    left, right = split_equation(equation)
    consumed = sum(species_coefficient(left, s) for s in REDUCED)
    produced = sum(species_coefficient(right, s) for s in REDUCED)
    return consumed, produced


def biomass_equations() -> dict[str, str]:
    df = pd.read_csv(DATA / "condition_specific_biomass_reactions.csv")
    return dict(zip(df["condition"], df["biomass_reaction"]))


def contribution(equation: str, flux: float, reaction_no: int) -> dict[str, float]:
    """Electron accounting for one fitted net flux value."""
    if not np.isfinite(flux):
        return {"source": 0.0, "sink": 0.0, "cumulative": 0.0}
    consumed, produced = reduced_counts(equation)
    if flux < 0:
        consumed, produced = produced, consumed
    magnitude = abs(float(flux))
    source = 2.0 * magnitude * max(produced - consumed, 0.0)
    sink = 2.0 * magnitude * max(consumed - produced, 0.0)
    cumulative = 2.0 * magnitude * max(produced, consumed)

    # Direct cofactor-independent transfer of the two formate electrons to H2.
    if int(reaction_no) == 39:
        source += 2.0 * magnitude
        sink += 2.0 * magnitude
        cumulative += 2.0 * magnitude
    return {"source": source, "sink": sink, "cumulative": cumulative}


def interval_extrema(equation: str, lb: float, ub: float, reaction_no: int, field: str) -> tuple[float, float]:
    """Conservative marginal interval for one reaction contribution."""
    candidates = [contribution(equation, float(lb), reaction_no)[field], contribution(equation, float(ub), reaction_no)[field]]
    if lb <= 0.0 <= ub:
        candidates.append(contribution(equation, 0.0, reaction_no)[field])
    return min(candidates), max(candidates)


def numeric(value):
    try:
        x = float(value)
        return x if np.isfinite(x) else np.nan
    except Exception:
        return np.nan


def main():
    fluxes = pd.read_csv(DATA / "reaction_fluxes.csv")
    meta = pd.read_csv(DATA / "condition_metadata.csv")
    biomass = biomass_equations()

    accounting_rows = []
    metrics = []

    for display, key in CONDITION_KEYS.items():
        best_col, lb_col, ub_col = f"{key}_best", f"{key}_lb95", f"{key}_ub95"
        source_best = sink_best = cumulative_best = 0.0
        source_min = source_max = 0.0
        cumulative_min = cumulative_max = 0.0

        for _, row in fluxes.iterrows():
            reaction_no = int(row["reaction_no"])
            equation = str(row["reaction"])
            if reaction_no == 82:
                equation = biomass[key]
            best = numeric(row[best_col])
            lb, ub = numeric(row[lb_col]), numeric(row[ub_col])
            if np.isnan(best):
                continue
            c = contribution(equation, best, reaction_no)
            source_best += c["source"]
            sink_best += c["sink"]
            cumulative_best += c["cumulative"]
            if not np.isnan(lb) and not np.isnan(ub):
                slo, shi = interval_extrema(equation, lb, ub, reaction_no, "source")
                clo, chi = interval_extrema(equation, lb, ub, reaction_no, "cumulative")
                source_min += slo; source_max += shi
                cumulative_min += clo; cumulative_max += chi

            accounting_rows.append({
                "condition": display,
                "reaction_no": reaction_no,
                "reaction": equation,
                "flux_best_per_100_substrate": best,
                "electron_source_per_100_substrate": c["source"],
                "electron_sink_per_100_substrate": c["sink"],
                "cumulative_electron_flux_per_100_substrate": c["cumulative"],
            })

        m = meta.loc[meta.condition == display].iloc[0]
        q = float(m.substrate_uptake_mmol_gDW_h)
        net_per_substrate = source_best / 100.0
        cumulative_per_substrate = cumulative_best / 100.0
        chemistry = SUBSTRATE_CHEMISTRY[str(m.substrate)]
        c_atoms = chemistry["C"]
        gamma_substrate = chemistry["gamma"]
        metrics.append({
            "condition": display,
            "substrate": m.substrate,
            "state": m.state,
            "q_substrate": q,
            "q_sd": float(m.substrate_uptake_sd),
            "substrate_C": c_atoms,
            "substrate_gamma": gamma_substrate,
            "net_e_per_substrate": net_per_substrate,
            "net_e_lb_envelope": source_min / 100.0,
            "net_e_ub_envelope": source_max / 100.0,
            "net_e_per_C": net_per_substrate / c_atoms,
            "net_e_lb_per_C": (source_min / 100.0) / c_atoms,
            "net_e_ub_per_C": (source_max / 100.0) / c_atoms,
            "net_e_fraction_of_substrate_gamma": net_per_substrate / gamma_substrate,
            "net_e_fraction_gamma_lb": (source_min / 100.0) / gamma_substrate,
            "net_e_fraction_gamma_ub": (source_max / 100.0) / gamma_substrate,
            "cumulative_e_per_substrate": cumulative_per_substrate,
            "cumulative_lb_envelope": cumulative_min / 100.0,
            "cumulative_ub_envelope": cumulative_max / 100.0,
            "effective_transfer_depth": cumulative_best / source_best,
            "absolute_net_e_flux": q * net_per_substrate,
            "absolute_cumulative_e_flux": q * cumulative_per_substrate,
            "source_sink_imbalance_e_per_substrate": (source_best - sink_best) / 100.0,
        })

    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(OUT / "gonzalez_13c_mfa_electron_flux_metrics.csv", index=False)
    accounting = pd.DataFrame(accounting_rows)
    accounting.to_csv(OUT / "gonzalez_13c_mfa_reaction_electron_accounting.csv", index=False)

    # Two-carbon-source × two-oxygen-state normalization robustness.  The paired
    # aerobic-minus-anaerobic expansion is expressed per substrate, per carbon
    # atom and relative to the substrate degree of reduction.
    norm_rows = []
    for substrate in ("glucose", "xylose"):
        g = metrics_df.loc[metrics_df.substrate.eq(substrate)].set_index("state")
        aer, anox = g.loc["aerobic"], g.loc["anaerobic"]
        norm_rows.append({
            "substrate": substrate,
            "substrate_C": float(aer.substrate_C),
            "substrate_gamma": float(aer.substrate_gamma),
            "aerobic_net_e_per_substrate": float(aer.net_e_per_substrate),
            "anaerobic_net_e_per_substrate": float(anox.net_e_per_substrate),
            "expansion_e_per_substrate": float(aer.net_e_per_substrate - anox.net_e_per_substrate),
            "aerobic_net_e_per_C": float(aer.net_e_per_C),
            "anaerobic_net_e_per_C": float(anox.net_e_per_C),
            "expansion_e_per_C": float(aer.net_e_per_C - anox.net_e_per_C),
            "aerobic_fraction_substrate_gamma": float(aer.net_e_fraction_of_substrate_gamma),
            "anaerobic_fraction_substrate_gamma": float(anox.net_e_fraction_of_substrate_gamma),
            "expansion_fraction_substrate_gamma": float(aer.net_e_fraction_of_substrate_gamma - anox.net_e_fraction_of_substrate_gamma),
            "anaerobic_low_per_C": float(anox.net_e_lb_per_C),
            "anaerobic_high_per_C": float(anox.net_e_ub_per_C),
            "anaerobic_low_fraction_gamma": float(anox.net_e_fraction_gamma_lb),
            "anaerobic_high_fraction_gamma": float(anox.net_e_fraction_gamma_ub),
        })
    pd.DataFrame(norm_rows).to_csv(OUT / "gonzalez_13c_mfa_normalization_summary.csv", index=False)

    # Focused anaerobic source decomposition used in the manuscript/SI.
    selected = {
        "Anaerobic glucose": {6: "GAPD", 9: "G6PDH", 10: "GND", 22: "PDH", 25: "ICDH", 48: "Serine biosynthesis", 39: "Formate → H2", 63: "Histidine biosynthesis"},
        "Anaerobic xylose": {6: "GAPD", 39: "Formate → H2", 22: "PDH", 33: "NADP malic enzyme", 25: "ICDH", 48: "Serine biosynthesis", 63: "Histidine biosynthesis"},
    }
    source_rows = []
    for cond, mapping in selected.items():
        total = float(metrics_df.loc[metrics_df.condition == cond, "net_e_per_substrate"].iloc[0])
        for reaction_no, label in mapping.items():
            row = accounting[(accounting.condition == cond) & (accounting.reaction_no == reaction_no)].iloc[0]
            per_sub = float(row.electron_source_per_100_substrate) / 100.0
            source_rows.append({
                "condition": cond,
                "reaction_no": reaction_no,
                "source": label,
                "electron_source_per_substrate": per_sub,
                "fraction_of_net_source": per_sub / total,
            })
    pd.DataFrame(source_rows).to_csv(OUT / "gonzalez_anaerobic_source_decomposition.csv", index=False)

    # Publication-facing validation checks.
    checks = {
        "max_abs_source_sink_imbalance_e_per_substrate": float(metrics_df.source_sink_imbalance_e_per_substrate.abs().max()),
        "aerobic_glucose_implied_O2_uptake": None,
        "aerobic_xylose_implied_O2_uptake": None,
    }
    # Reactions 67 and 68 consume one reduced cofactor per 1/2 O2 in the 13C-MFA model.
    for cond, label in [("Aerobic glucose", "aerobic_glucose_implied_O2_uptake"), ("Aerobic xylose", "aerobic_xylose_implied_O2_uptake")]:
        sub = accounting[(accounting.condition == cond) & (accounting.reaction_no.isin([67, 68]))]
        # Electron sink/4 = mol O2 equivalent.
        checks[label] = float(sub.electron_sink_per_100_substrate.sum() / 4.0 / 100.0 * float(meta.loc[meta.condition == cond, "substrate_uptake_mmol_gDW_h"].iloc[0]))
    (OUT / "gonzalez_13c_mfa_validation_checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(metrics_df.to_string(index=False))
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
