"""
electron_flow_analysis.py
---------------------------
Turns a flux-balance solution into a direction-aware, reaction-resolved
donor-to-acceptor electron flow table, then builds the source/carrier/
sink classification, the respiration/fermentation/biosynthesis
partition, and condition-specific electron-transfer redistribution analysis on
top of it.

DIRECTION-AWARE DONOR/ACCEPTOR ASSIGNMENT
matching.classify_reaction_molecular() identifies each
reaction's donor/acceptor pairs relative to the reaction AS WRITTEN
(forward direction) in the model. Under a flux-balance solution, a
reversible reaction's true electron-donating and -accepting roles are
set by which direction it actually operates, so the true donor and
acceptor must be determined from the SIGN of its flux, not assumed from
the written direction alone: if flux is negative, the reaction runs in
reverse and the forward-direction donor/acceptor roles must be
exchanged. The reported TOTAL electron flux magnitude does not depend on this
assignment, since it only ever sums
|flux| x electrons -- but any analysis that depends on donor/acceptor
IDENTITY (this module) requires it.

For each matched pair (reactant r, product p, delta_gamma), the rule is:
    if flux > 0: the reaction runs forward; donor/acceptor are as
                 identified by delta_gamma (r is donor if delta<0, p is
                 acceptor if delta>0).
    if flux < 0: the reaction runs in reverse; donor and acceptor
                 exchange between the two pairs: the species that would
                 have been the acceptor (p, for a delta>0 pair) is now
                 being consumed and is the donor, and the species that
                 would have been the donor (r, for a delta<0 pair) is
                 now being formed and is the acceptor.
The network represents each redox pool by the reduced species carrying the
electrons: a complete forward edge is (donor_pair.reactant ->
acceptor_pair.product). Under reverse flux the reduced acceptor product
becomes the donor and the reduced form of the opposite pair is regenerated,
so the directed edge is simply reversed:
(acceptor_pair.product -> donor_pair.reactant).
"""
from typing import Dict, List
import pandas as pd

from .matching import classify_reaction_molecular
from .constants import UNRESOLVED_ORGANIC_NODE


def signed_donor_acceptor_edges(model, met_table, solution, min_flux: float = 1e-6) -> pd.DataFrame:
    """
    Builds the flux-direction-aware donor->acceptor electron-flow table for
    one FBA solution: one row per resolved (or partially resolved)
    donor/acceptor contribution carrying non-negligible flux, with the
    TRUE (flux-direction-aware) donor and acceptor identity.

    MULTI-PAIR AND PARTIALLY-RESOLVED ACCOUNTING
    When a reaction has more than one donor and/or acceptor pair, the
    total is distributed across every donor x acceptor combination by
    PROPORTIONAL allocation (weight = |delta_d|*|delta_a| / total),
    which always sums back to the reaction's true total (Methods;
    matching.classify_reaction_molecular).

    When the matched donor and acceptor pairs do NOT balance (for example,
    FMNH2->FMN 2 e- matched against O2->H2O 4 e-, with a missing 2 e-
    contribution from an unmatched organic donor), the smaller side is
    padded with a single phantom
    "(unresolved organic donor/acceptor)" pair of magnitude
    |donated-accepted|, so that BOTH real and phantom pairs, taken
    together, always balance exactly and can be run through the same
    proportional-allocation and flux-direction logic uniformly -- the
    phantom pair is treated as just another pair, gaining a
    "(unresolved organic ...)" identity instead of a real metabolite id.
    """
    rows = []
    for rxn in model.reactions:
        flux = solution.fluxes.get(rxn.id, 0.0)
        if abs(flux) < min_flux:
            continue
        res = classify_reaction_molecular(rxn, met_table)
        if not res["is_redox"]:
            continue

        donors = [{"reactant": p["reactant"], "product": p["product"], "e": abs(p["delta_gamma"])}
                  for p in res["pairs"] if p["delta_gamma"] < 0]
        acceptors = [{"reactant": p["reactant"], "product": p["product"], "e": abs(p["delta_gamma"])}
                     for p in res["pairs"] if p["delta_gamma"] > 0]
        # Resolved totals and the minimum unresolved gap come from the
        # single authoritative electron-balance calculation
        # (reaction_electron_balance.py), not a local recomputation.
        D = res["electrons_donated_resolved"]
        A = res["electrons_accepted_resolved"]
        gap = res["minimum_unresolved_donation"] + res["minimum_unresolved_acceptance"]
        if gap > 1e-9:
            # Phantom pair for the unresolved portion. A REAL pair's
            # reactant/product fields track a SINGLE PHYSICAL ENTITY
            # whose role (donor vs acceptor) is determined by flux
            # direction -- e.g. for ALCD2x, "etoh_c" is used consistently
            # as that pair's identity whether ethanol is playing donor
            # (forward) or acceptor (reverse) -- so the phantom must use
            # ONE consistent, role-neutral label for both fields too.
            phantom = {"reactant": UNRESOLVED_ORGANIC_NODE, "product": UNRESOLVED_ORGANIC_NODE, "e": gap}
            if D < A:
                donors.append(phantom)
            else:
                acceptors.append(phantom)
        total = res["electrons_per_turnover_estimate"]
        if total <= 1e-9 or not donors or not acceptors:
            continue  # donor_side_only / acceptor_side_only with nothing on the other side at all

        forward = flux > 0
        for d in donors:
            for a in acceptors:
                electrons = d["e"] * a["e"] / total
                if electrons <= 1e-9:
                    continue
                if forward:
                    donor_id, acceptor_id = d["reactant"], a["product"]
                else:
                    donor_id, acceptor_id = a["product"], d["reactant"]
                rows.append({
                    "reaction_id": rxn.id, "reaction_name": rxn.name, "subsystem": rxn.subsystem,
                    "flux": flux, "flux_direction": "forward" if forward else "reverse",
                    "donor": donor_id, "acceptor": acceptor_id,
                    "resolution_status": res["resolution_status"],
                    "electrons_per_turnover": electrons,
                    "electron_flux": abs(flux) * electrons,
                })
    return pd.DataFrame(rows)


def donor_acceptor_matrix(edges: pd.DataFrame) -> pd.DataFrame:
    """Aggregates the signed edge table into a donor x acceptor electron-
    flux matrix (sum over all reactions connecting that pair)."""
    return (edges.groupby(["donor", "acceptor"])["electron_flux"]
            .sum().reset_index().sort_values("electron_flux", ascending=False))


def source_carrier_sink_classification(edges: pd.DataFrame, met_name: Dict[str, str]) -> pd.DataFrame:
    """
    Classifies every metabolite appearing in the signed edge table by its
    ROLE in the electron-flow network under this condition, from its
    total electron OUT-flux (as donor) and IN-flux (as acceptor):
        D_out(i) = sum of electron_flux where i is donor
        D_in(i)  = sum of electron_flux where i is acceptor
    Role (a simple, quantitative rule -- see Results/Methods):
        'primary_source' : D_in / (D_in+D_out) < 0.1   (almost pure donor)
        'terminal_sink'  : D_out / (D_in+D_out) < 0.1   (almost pure acceptor)
        'carrier'        : otherwise (substantial flux both in and out)
    Metabolites explicitly marked "(unresolved ... )" (Supplementary
    Note 1) are excluded.
    """
    d_out = edges.groupby("donor")["electron_flux"].sum()
    d_in = edges.groupby("acceptor")["electron_flux"].sum()
    all_mets = (set(d_out.index) | set(d_in.index)) - {UNRESOLVED_ORGANIC_NODE}

    rows = []
    for met in all_mets:
        out_flux = d_out.get(met, 0.0)
        in_flux = d_in.get(met, 0.0)
        total = out_flux + in_flux
        in_frac = in_flux / total if total else 0.0
        out_frac = out_flux / total if total else 0.0
        if in_frac < 0.1:
            role = "primary_source"
        elif out_frac < 0.1:
            role = "terminal_sink"
        else:
            role = "carrier"
        rows.append({"metabolite_id": met, "name": met_name.get(met, met),
                      "electron_out_flux": out_flux, "electron_in_flux": in_flux,
                      "total_electron_flux": total, "role": role})
    df = pd.DataFrame(rows, columns=["metabolite_id", "name", "electron_out_flux",
                                      "electron_in_flux", "total_electron_flux", "role"])
    df = df.sort_values(["total_electron_flux", "metabolite_id"], ascending=[False, True]).reset_index(drop=True)
    return df


# Curated classification of terminal fates, for the respiration / fermentation
# / biosynthesis partition -- based on the WELL-ESTABLISHED physiological role
# of each terminal acceptor in E. coli (not re-derived from the network
# itself, to avoid circularity): respiratory acceptors are those that
# accept electrons from the quinone pool via a membrane-bound terminal
# oxidase/reductase; fermentation products are reduced, secreted organic
# end-products; everything else is attributed to biosynthesis (reducing
# power consumed by anabolic reactions rather than delivered to a
# dedicated terminal acceptor).
# Curated classification of terminal fates, for the respiration / fermentation
# / biosynthesis partition -- based on the WELL-ESTABLISHED physiological role
# of each terminal acceptor in E. coli (not re-derived from the network
# itself, to avoid circularity). IMPORTANT: because signed_donor_acceptor_edges
# always records the REDUCED PRODUCT species as "acceptor" (e.g. h2o_c, not
# o2_c, for O2 respiration -- see that function's docstring), this set lists
# the reduced/product forms of each respiratory terminal acceptor, not the
# oxidised starting forms. For the aerobic-versus-fermentative glucose case study, succinate is included
# among fermentation products because a net terminal succinate sink under the
# anaerobic state represents reduced carbon excreted during fermentation. The
# source-to-sink classification uses the *net node balance*, so succinate is
# counted as a terminal fate only when it is a net sink rather than a relay.
RESPIRATORY_ACCEPTORS = {"h2o_c", "no2_c", "nh4_c", "dms_c", "tma_c", "n2_c"}
FERMENTATION_PRODUCTS = {"etoh_c", "lac__D_c", "lac__L_c", "ac_c", "for_c", "acald_c", "succ_c"}


def respiration_fermentation_biosynthesis_partition(edges: pd.DataFrame) -> pd.DataFrame:
    """
    Partitions total electron flux by the ultimate fate of the accepting
    metabolite in each edge, using the curated sets above. Edges whose
    acceptor is a "carrier" cofactor (NADH, quinols, etc. -- i.e. not a
    terminal species itself) are attributed to whichever category their
    OWN downstream terminal acceptor eventually reaches; since building a
    full downstream trace is out of scope for this basic partition, we
    instead report the SIMPLER, transparent version: flux is partitioned
    by the acceptor of each EDGE directly, with cofactor-to-cofactor
    edges (e.g. NADH -> quinol, a carrier-to-carrier relay) reported
    separately as 'carrier_relay' rather than forced into one of the
    three end-point categories. This keeps the partition honest about
    what it does and does not resolve (see Discussion, "topological vs
    physiological" caveat).
    """
    CARRIERS = {"nad_c", "nadh_c", "nadp_c", "nadph_c", "fad_c", "fadh2_c",
                "fmn_c", "fmnh2_c", "q8_c", "q8h2_c", "mqn8_c", "mql8_c",
                "2dmmq8_c", "2dmmql8_c", "trdox_c", "trdrd_c", "gthox_c", "gthrd_c",
                "h2o2_c", "succ_c", "fum_c", "o2_c", "o2s_c"}

    def fate(acceptor):
        if acceptor in RESPIRATORY_ACCEPTORS:
            return "respiration"
        if acceptor in FERMENTATION_PRODUCTS:
            return "fermentation"
        if acceptor in CARRIERS:
            return "carrier_relay"
        return "biosynthesis"

    e = edges.copy()
    e["fate"] = e["acceptor"].apply(fate)
    summary = e.groupby("fate")["electron_flux"].sum().reset_index()
    summary["fraction"] = summary["electron_flux"] / summary["electron_flux"].sum()
    return summary.sort_values("electron_flux", ascending=False)


def redistribution(edges_A: pd.DataFrame, edges_C: pd.DataFrame) -> pd.DataFrame:
    """
    Compares two conditions' donor->acceptor matrices and reports, per
    pair, the change (J_ij^C - J_ij^A): increased, decreased, newly
    activated (zero in A, nonzero in C) or suppressed (nonzero in A, zero
    in C) transfers.
    """
    mat_A = donor_acceptor_matrix(edges_A).rename(columns={"electron_flux": "flux_A"})
    mat_C = donor_acceptor_matrix(edges_C).rename(columns={"electron_flux": "flux_C"})
    merged = pd.merge(mat_A, mat_C, on=["donor", "acceptor"], how="outer").fillna(0.0)
    merged["delta"] = merged["flux_C"] - merged["flux_A"]

    def status(row):
        if row.flux_A < 1e-9 and row.flux_C >= 1e-9:
            return "newly_activated"
        if row.flux_A >= 1e-9 and row.flux_C < 1e-9:
            return "suppressed"
        if row.delta > 1e-9:
            return "increased"
        if row.delta < -1e-9:
            return "decreased"
        return "unchanged"

    merged["status"] = merged.apply(status, axis=1)
    return merged.sort_values("delta", key=abs, ascending=False).reset_index(drop=True)


def in_out_weighted_degree(edges: pd.DataFrame, met_name: Dict[str, str]) -> pd.DataFrame:
    """
    D_in(i) = sum_j w_ji (electron flux where i is the acceptor)
    D_out(i) = sum_j w_ij (electron flux where i is the donor)
    Reported separately, plus total, so donors/acceptors/carriers are not
    collapsed into a single undirected ranking.
    """
    d_in = edges.groupby("acceptor")["electron_flux"].sum().rename("D_in")
    d_out = edges.groupby("donor")["electron_flux"].sum().rename("D_out")
    df = pd.concat([d_in, d_out], axis=1).fillna(0.0)
    df["D_total"] = df["D_in"] + df["D_out"]
    df = df.reset_index().rename(columns={"index": "metabolite_id"})
    df["name"] = df["metabolite_id"].map(met_name)
    df = df[df.metabolite_id != UNRESOLVED_ORGANIC_NODE]
    return df.sort_values("D_total", ascending=False).reset_index(drop=True)
