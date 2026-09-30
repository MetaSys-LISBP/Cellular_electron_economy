# Results manifest

This manifest groups the principal outputs supporting **Organizing principles of the cellular electron economy**. It is an analysis-level guide rather than an exhaustive list of generated files. Exact panel-to-file mappings are provided in `FIGURE_SOURCE_DATA_MANIFEST.csv`.

## 1. Electron-transfer-network reconstruction and validation

Principal outputs:

- `results/electron_transfer_table.tsv` — reaction-level electron-transfer assignments and stoichiometries;
- `results/electron_transfer_edges.csv` — complete chemistry-resolved directed donor–acceptor reaction edges, including the chemically resolved but model-blocked FHL edge;
- `results/matched_pairs_all.csv` — reactant–product matching results;
- `results/ambiguous_reactions.csv` — reactions retained as chemically ambiguous;
- `results/coverage_summary.csv` and `results/coverage_per_reaction.csv` — reconstruction coverage;
- `results/electron_balance_diagnostics.csv` — electron-balance diagnostics;
- `results/metabolite_degree_of_reduction.csv` — calculated metabolite degrees of reduction;
- `results/network_basic_stats.json` — model-direction-aware topology statistics plus the complete reaction-edge count;
- `results/iML1515_ETN_boundary_nodes.xml` — SBML representation with duplicated boundary source/sink nodes.

Related figures and Source Data:

- `figures/publication/SuppFigS1_reconstruction_details.*`
- `figures/publication/SuppFigS2_iML1515_ETN.*`
- `source_data/publication/SuppFigS1_*`
- `source_data/publication/SuppFigS2_*`

## 2. *E. coli* aerobic/anaerobic electron economy

Principal outputs:

- `results/fba_condition_summary.csv`
- `results/fba_redox_reactions_aerobic.csv`
- `results/fba_redox_reactions_anaerobic.csv`
- `results/electron_path_summary_by_condition.csv`
- `results/electron_carrier_relay_by_condition.csv`
- `results/electron_sink_delivery_by_condition.csv`
- `results/electron_terminal_fates_interpreted_by_condition.csv`
- `results/fva_robustness_summary.csv`
- `results/publication/gonzalez_13c_mfa_electron_flux_metrics.csv`
- `results/publication/ecom4la_state_constrained_summary.csv`
- `results/publication/anand2022_architecture_summary.csv`

Related figures and Source Data:

- main Fig. 2;
- Supplementary Figs. S3–S7;
- `source_data/publication/Fig2*`;
- `source_data/publication/SuppFigS3_*` to `SuppFigS7_*`.

## 3. Carbon source and disposal capacity

Principal outputs:

- `results/publication/acceptor_titration.csv`
- `results/publication/acceptor_electron_capacity_titration.csv`
- `results/publication/acceptor_carbon_oxidation_statistics.json`
- `results/publication/acceptor_endpoints.csv`
- `results/publication/substrate_flux_compression.csv`
- `results/publication/substrate_core_redox_contribution.csv`
- `results/publication/factorial_carbon_oxygen_states.csv`
- `results/publication/factorial_scaled_expansion_curves.csv`
- `results/publication/factorial_carbon_oxygen_statistics.json`
- `results/publication/denby2015_*`
- `results/publication/toya2012_*`
- `results/publication/perrenoud2005_*`

Related figures and Source Data:

- main Fig. 3;
- Supplementary Figs. S5 and S8–S12;
- `source_data/publication/Fig3*`;
- corresponding `source_data/publication/SuppFigS5_*` and `SuppFigS8_*` to `SuppFigS12_*` files.

## 4. Cross-species organization and robustness

Principal outputs:

- `results/cross_species/cross_species_physiology.csv`
- `results/cross_species/cross_species_carrier_relay.csv`
- `results/cross_species/electron_flux_compression_summary.csv`
- `results/publication/cross_species_carbon_acceptor_feasibility.csv`
- `results/publication/cross_species_carbon_acceptor_pairs.csv`
- `results/publication/cross_species_carbon_acceptor_normalization_summary.csv`
- `results/publication/cross_species_routing_redistribution.csv`
- `results/publication/extended_acceptor_screen.csv`
- `results/publication/extended_carbon_entry_summary.csv`

## 5. Independent electron-generation and endpoint analyses

- Jouhten *S. cerevisiae*: `results/publication/jouhten2008_*`;
- Baumann *K. phaffii*: `results/publication/baumann2010_13c_mfa_*`;
- Fromanger *C. shehatae*: `results/publication/fromanger2010_electron_balance.csv` and its summary JSON.

Electron generation is calculated from electron-generating reactions in the published metabolic flux distributions. It is kept distinct from balanced net electron flux. The Fromanger analysis instead calculates a whole-cell electron-equivalent partition from reported carbon yields and product balances.

## 6. Predictive-transfer analyses

- respiratory architectures: `results/publication/steinsiek2014_*`;
- carbon sources: `results/publication/weusthuis1994_*`.

These outputs support main Fig. 4 and Supplementary Figs. S19–S20. Cross-species routing and carbon-source analyses additionally support Supplementary Figs. S13 and S15–S18.

## Publication-facing artifacts

- `figures/publication/` — generated computational main and supplementary figures;
- `source_data/publication/` — panel-level figure Source Data;
- `FIGURE_SOURCE_DATA_MANIFEST.csv` — exact figure-to-file mapping;
- `INPUT_PROVENANCE_SHA256.csv` — input provenance and checksums.
