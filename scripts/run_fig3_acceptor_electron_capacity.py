#!/usr/bin/env python3
"""Run only the exact acceptor simulations used in Fig. 3A–B.

The x axis is nominal electron-accepting capacity, not molar acceptor uptake.
O2, nitrate and TMAO are each simulated at 0, 5, 10, 20, 30, 40, 50, 60
and 80 e- mmol gDW^-1 h^-1. Every plotted point is solved explicitly.
"""
from __future__ import annotations

import json

from run_mechanistic_analysis import (
    load_ecoli_model,
    chemical_model,
    build_metabolite_table,
    build_matching_context,
    run_fig3_acceptor_electron_capacity_panel,
)


def main() -> None:
    model = load_ecoli_model()
    chem = chemical_model(model)
    # Build/caches the matching context once, then reuse it across all 27 states.
    mt = build_metabolite_table(chem)
    build_matching_context(chem, mt)
    _, stats = run_fig3_acceptor_electron_capacity_panel(model, chem)
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
