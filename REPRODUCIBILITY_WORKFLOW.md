# Reproducibility workflow

Run all commands from the repository root using the reference environment in `requirements.txt`. The workflow calculates electron-transfer-network and electron-economy quantities from the bundled metabolic reconstructions, experimental metabolic flux distributions and model-derived flux distributions.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate  # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The reference environment uses Python 3.13.5. Input integrity can be checked against `INPUT_PROVENANCE_SHA256.csv` and the dataset-specific checksum files.

## Complete workflow

```bash
python scripts/run_all.py
```

This command performs the stages below in order and generates:

- reconstruction, validation and network outputs in `results/`;
- publication analyses in `results/publication/`;
- cross-species outputs in `results/cross_species/`;
- panel-level Source Data in `source_data/publication/`;
- the supplementary Excel table;
- computational main Figs. 2–4 and Supplementary Figs. S1–S20 in `figures/publication/`.

The complete workflow is intentionally comprehensive. The staged commands below expose progress and allow individual analysis groups to be rerun.

## Stage 1 — Chemical reconstruction and validation

```bash
python scripts/validate_reduction_degree.py
python scripts/matching_validation.py
python scripts/run_validation.py
python scripts/validation_battery.py
python scripts/run_analysis.py
python scripts/coverage_analysis.py
```

Principal outputs include:

- `results/metabolite_degree_of_reduction.csv`
- `results/electron_transfer_table.tsv`
- `results/electron_transfer_edges.csv`
- `results/matched_pairs_all.csv`
- `results/ambiguous_reactions.csv`
- `results/coverage_summary.csv`
- `results/electron_balance_diagnostics.csv`
- `results/iML1515_ETN_boundary_nodes.xml`

## Stage 2 — *E. coli* electron flux and network analyses

```bash
python scripts/run_electron_flow_analyses.py
python scripts/run_electron_path_decomposition.py
python scripts/compare_networks.py
python scripts/analyse_community_resolution_sensitivity.py
python scripts/analyse_community_robustness.py
python scripts/analyse_network_comparison_sensitivity.py
python scripts/analyse_flux_robustness.py
```

This stage calculates net electron flux, reaction-level electron fluxes, carrier relay, effective transfer depth, terminal electron delivery, topology and robustness for the aerobic and anaerobic *E. coli* states.

## Stage 3 — Experimental calculations and respiratory perturbations

```bash
python scripts/run_gonzalez_13c_mfa_validation.py
python scripts/run_jouhten_13c_mfa_validation.py
python scripts/run_jouhten_matched_model_validation.py
python scripts/run_acceptor_experimental_validation.py
python scripts/run_mechanistic_analysis.py
python data/experimental/portnoy2010_ecom4la/run_ecom4la_state_constrained.py
python scripts/run_respiratory_architecture_support.py
python scripts/run_factorial_carbon_disposal_analysis.py
python scripts/run_fromanger2010_electron_balance.py
python scripts/run_baumann2010_13c_mfa_validation.py
python scripts/run_steinsiek2014_predictive_transfer.py
python scripts/run_weusthuis1994_predictive_transfer.py
```

The Gonzalez analysis calculates balanced electron-flux quantities from published <sup>13</sup>C-MFA reaction-flux distributions. The Jouhten, Toya and Baumann analyses calculate electron generation from the electron-generating reactions resolved by the corresponding published metabolic flux distributions. These observables are not treated as interchangeable.

## Stage 4 — Cross-species and coverage analyses

```bash
python scripts/run_cross_species_analysis.py
python scripts/run_cross_species_principle_analysis.py --organism "B. subtilis"
python scripts/run_cross_species_principle_analysis.py --organism "S. enterica"
python scripts/run_cross_species_principle_analysis.py --organism "K. phaffii"
python scripts/run_cross_species_principle_analysis.py --merge
python scripts/run_extended_acceptor_screen.py
python scripts/run_extended_carbon_entry_screen.py
python scripts/run_cross_species_routing_redistribution.py
```

The staged organism commands generate organism-specific outputs before the merged cross-species tables are assembled. Infeasible organism–substrate states are retained explicitly.

## Stage 5 — Publication tables, Source Data and figures

```bash
python scripts/make_excel_table.py
python scripts/make_source_data.py
python scripts/make_figures.py
```

`make_source_data.py` writes the panel-level files listed in `FIGURE_SOURCE_DATA_MANIFEST.csv`. `make_figures.py` calls `scripts/main_figures.py` for main Figs. 2–4 and generates Supplementary Figs. S1–S20. Supplementary Fig. S2 is generated from the reconstructed network tables by `scripts/make_ecoli_etn_network_figure.py`.

Important figure-specific analyses include:

- Fig. 2 and Supplementary Figs. S3–S6 and S11: *E. coli* decoupling, Gonzalez, ECOM4LA and respiratory-architecture analyses;
- Fig. 3 and Supplementary Figs. S7–S12: acceptors, carbon sources and the factorial response surface;
- Fig. 4 and Supplementary Figs. S13–S20: cross-species analyses, electron generation, whole-cell accounting and predictive transfer.

## Tests

Run the tests separately after the analysis workflow:

```bash
pytest -q
```

`scripts/run_all.py` does not run `pytest`. The tests cover chemical accounting, matching, reconstruction coverage, calculated experimental observables, cross-species analyses, figure contracts, predictive-transfer analyses and manuscript-facing numerical claims.

## Output inventories

- `RESULTS_MANIFEST.md` groups the principal results by scientific analysis.
- `FIGURE_SOURCE_DATA_MANIFEST.csv` maps figures to Source Data files.
- `source_data/publication/README.md` explains figure-level observable and uncertainty boundaries.
- `INPUT_PROVENANCE_SHA256.csv` records the input files and checksums.

Generated outputs are bundled with the repository, so they can be inspected without rerunning the complete workflow.
