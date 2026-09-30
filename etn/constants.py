"""
constants.py
-------------
Shared constants used across the package: the unresolved-partner
identifier, numerical tolerance, and the biomass/exchange
reaction prefixes used to exclude non-enzymatic pseudo-reactions from
classification.
"""

# Identity for the unresolved electron-donating or electron-accepting partner of
# a reaction whose organic side could not be formula-matched (e.g. a
# reaction where only a cofactor pair, such as NAD+/NADH, is resolved,
# leaving the true organic substrate/product unidentified). Used only in
# the condition-specific, flux-weighted electron-activity representation
# (etn/electron_flow_analysis.py) -- never in the static electron-transfer
# graph (etn/network_analysis.py), which contains only fully chemically
# resolved donor/acceptor relationships.
UNRESOLVED_ORGANIC_NODE = "(unresolved organic)"

# Numerical tolerance for floating-point equality checks on reduction
# degree / electron-balance quantities.
TOLERANCE = 1e-9

# Standard BiGG-namespace prefixes for structural pseudo-reactions
# (biomass objective, exchange/demand/sink) that are not enzymatic
# transformations and are excluded from electron-transfer classification
# a priori.
EXCLUDED_REACTION_PREFIXES = ("EX_", "DM_", "SK_")
BIOMASS_PREFIX = "BIOMASS"
