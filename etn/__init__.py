"""
etn — Electron Transfer Network reconstruction and genome-scale electron flux
analysis.

Reconstructs, validates, and topologically characterises the electron-transfer
network of a metabolic reconstruction (E. coli iML1515), using
each metabolite's molecular formula and charge (the "reduction degree" of
Roels, 1980) to identify electron donors and acceptors without requiring a
manually curated redox annotation for every reaction.

Typical usage
-------------
>>> from etn.data_acquisition import load_ecoli_model
>>> from etn.matching import build_metabolite_table, classify_reactions_molecular
>>> from etn.electron_transfer_table import build_electron_transfer_table
>>>
>>> model = load_ecoli_model()
>>> met_table = build_metabolite_table(model)
>>> classification = classify_reactions_molecular(model)

See scripts/run_analysis.py for the full reconstruction and analysis
pipeline, and README.md for an overview of the package layout and
methodology.
"""

__version__ = "1.0.0"

from .degree_of_reduction import degree_of_reduction, build_gamma_table
from .matching import (
    build_metabolite_table,
    classify_reaction_molecular,
    classify_reactions_molecular,
    match_metabolite_pair,
    KNOWN_REDOX_PAIRS,
)
from .electron_transfer_table import build_electron_transfer_table
from .network_analysis import (
    build_electron_graph,
    simple_weighted_graph,
    basic_stats,
    top_hubs,
    connected_components_report,
    detect_communities,
    build_bipartite_graph,
    subsystem_enrichment,
)

__all__ = [
    "degree_of_reduction", "build_gamma_table",
    "build_metabolite_table", "classify_reaction_molecular", "classify_reactions_molecular",
    "match_metabolite_pair", "KNOWN_REDOX_PAIRS", "build_electron_transfer_table",
    "build_electron_graph", "simple_weighted_graph", "basic_stats", "top_hubs",
    "connected_components_report", "detect_communities", "build_bipartite_graph",
    "subsystem_enrichment",
]
