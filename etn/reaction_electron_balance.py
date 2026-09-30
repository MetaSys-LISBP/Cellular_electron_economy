"""
reaction_electron_balance.py
-------------------------------
The single, authoritative implementation of reaction-level electron
accounting. Every module that needs a reaction's electron count calls
compute_reaction_electron_balance() defined here; no other module
re-implements this logic.

A reaction's matched donor/acceptor pairs (matching.py)
are combined into:

    electrons_donated_resolved     sum of |delta gamma| over donor pairs
                                    (delta gamma < 0), each backed by a
                                    specific, identified reactant/product
    electrons_accepted_resolved    sum of |delta gamma| over acceptor
                                    pairs (delta gamma > 0), likewise
                                    backed by specific identified species
    resolution_status              one of:
        fully_resolved          - donor and acceptor sides both present
                                   and equal within tolerance: every
                                   electron donated is accounted for by an
                                   identified acceptor and vice versa
        partially_resolved       - donor and acceptor sides both present
                                   but unequal: at least one side's total
                                   is under-identified
        donor_side_only          - only donor pairs identified, no
                                   acceptor pair matched at all
        acceptor_side_only       - only acceptor pairs identified, no
                                   donor pair matched at all
    minimum_unresolved_donation    the SMALLEST additional donor
                                    contribution that would be needed to
                                    balance the identified acceptor total
                                    (0 for fully_resolved). This is an
                                    INFERRED LOWER BOUND, not a measured
                                    quantity: a reaction could in
                                    principle carry additional, exactly
                                    cancelling unresolved donor AND
                                    acceptor contributions beyond this
                                    minimum, which no formula-level
                                    accounting can detect. The name
                                    deliberately avoids implying the true
                                    unresolved donation is known.
    minimum_unresolved_acceptance  symmetric definition for the acceptor
                                    side.
    electrons_per_turnover_estimate
                                    max(electrons_donated_resolved,
                                    electrons_accepted_resolved): the best
                                    available point estimate of total
                                    electron-transfer magnitude, exact for
                                    fully_resolved reactions (where both
                                    sides agree) and a lower bound for
                                    every other status.

For a fully_resolved reaction, electrons_donated_resolved ==
electrons_accepted_resolved is a required invariant, checked by
validate_electron_balance_invariants() below; violation raises
AssertionError rather than being silently accepted.
"""
from dataclasses import dataclass
from typing import List, Optional


@dataclass(frozen=True)
class ElectronBalance:
    electrons_donated_resolved: float
    electrons_accepted_resolved: float
    resolution_status: Optional[str]
    minimum_unresolved_donation: float
    minimum_unresolved_acceptance: float
    electrons_per_turnover_estimate: float


def compute_reaction_electron_balance(pairs: List[dict], tolerance: float = 1e-9) -> ElectronBalance:
    """
    pairs: the 'pairs' list returned by
    matching.classify_reaction_molecular (each a dict
    with at least 'delta_gamma'). An empty list (non-redox reaction)
    returns a zeroed, status=None result.
    """
    donated = sum(abs(p["delta_gamma"]) for p in pairs if p["delta_gamma"] < 0)
    accepted = sum(abs(p["delta_gamma"]) for p in pairs if p["delta_gamma"] > 0)
    has_donor = any(p["delta_gamma"] < 0 for p in pairs)
    has_acceptor = any(p["delta_gamma"] > 0 for p in pairs)

    if not has_donor and not has_acceptor:
        return ElectronBalance(0.0, 0.0, None, 0.0, 0.0, 0.0)

    if has_donor and has_acceptor:
        if abs(donated - accepted) <= tolerance:
            status = "fully_resolved"
        else:
            status = "partially_resolved"
    elif has_donor:
        status = "donor_side_only"
    else:
        status = "acceptor_side_only"

    min_unresolved_donation = max(0.0, accepted - donated)
    min_unresolved_acceptance = max(0.0, donated - accepted)
    estimate = max(donated, accepted)

    return ElectronBalance(
        electrons_donated_resolved=donated,
        electrons_accepted_resolved=accepted,
        resolution_status=status,
        minimum_unresolved_donation=min_unresolved_donation,
        minimum_unresolved_acceptance=min_unresolved_acceptance,
        electrons_per_turnover_estimate=estimate,
    )


def validate_electron_balance_invariants(balance: ElectronBalance, reaction_id: str = "") -> None:
    """
    Raises AssertionError if a conservation invariant that should always
    hold is violated. Intended to be called for every redox-classified
    reaction across the network (see scripts/run_validation.py).
    """
    resolved_plus_min_unresolved_donor = balance.electrons_donated_resolved + balance.minimum_unresolved_donation
    resolved_plus_min_unresolved_acceptor = balance.electrons_accepted_resolved + balance.minimum_unresolved_acceptance
    assert abs(resolved_plus_min_unresolved_donor - balance.electrons_per_turnover_estimate) <= 1e-6, (
        f"{reaction_id}: donor accounting inconsistent with the turnover estimate")
    assert abs(resolved_plus_min_unresolved_acceptor - balance.electrons_per_turnover_estimate) <= 1e-6, (
        f"{reaction_id}: acceptor accounting inconsistent with the turnover estimate")

    if balance.resolution_status == "fully_resolved":
        assert abs(balance.electrons_donated_resolved - balance.electrons_accepted_resolved) <= 1e-6, (
            f"{reaction_id}: fully_resolved reaction must have donated==accepted")
        assert balance.minimum_unresolved_donation == 0.0 and balance.minimum_unresolved_acceptance == 0.0, (
            f"{reaction_id}: fully_resolved reaction must have zero minimum-unresolved contributions")
    elif balance.resolution_status in ("donor_side_only", "acceptor_side_only"):
        assert (balance.minimum_unresolved_donation == 0.0) != (balance.minimum_unresolved_acceptance == 0.0) or \
            balance.electrons_per_turnover_estimate == 0.0, (
            f"{reaction_id}: one-sided reaction should have exactly one nonzero minimum-unresolved side")
