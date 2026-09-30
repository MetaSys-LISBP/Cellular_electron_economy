# Cellular electron economy

This repository accompanies the manuscript **Organizing principles of the cellular electron economy** (DOI: [XXXXXX](https://doi.org/XXXXXX)). It contains the code, metabolic reconstructions, curated experimental inputs and generated outputs used to reconstruct electron-transfer networks, calculate genome-scale electron fluxes and reproduce the quantitative analyses, figures and Source Data reported in the manuscript and Supplementary Information.

## What the framework calculates

The workflow has three connected levels:

1. **Electron transfer network reconstruction.** Reaction chemistry, metabolite formulae and charge are used to identify direction-aware electron donors and acceptors. Exact transport events are excluded, universal redox couples provide independent anchors, and unresolved or ambiguous reactions are retained explicitly in the reconstruction diagnostics.
2. **Genome-scale electron flux analysis.** A metabolic reaction-flux distribution is combined with reaction-specific electron-transfer stoichiometry to calculate reaction-resolved electron fluxes across the reconstructed network.
3. **Cellular electron economy.** The calculated fluxes are summarized as net electron flux, carrier relay flux, effective transfer depth and terminal electron delivery. Together, these quantities describe the magnitude and routing of electron flow without introducing an additional metabolic degree of freedom or conservation law.

Depending on flux map coverage, the code calculates either net source-to-sink electron flux from sufficiently resolved metabolic flux distributions or electron generation from the available electron-generating reactions.

The principal *Escherichia coli* electron transfer network outputs are `results/electron_transfer_table.tsv`, `results/electron_transfer_edges.csv` and `results/iML1515_ETN_boundary_nodes.xml`. Supplementary Fig. S2 and its node and edge tables provide the corresponding visual and machine-readable network representations. The complete chemistry-resolved edge registry contains 280 directed reaction edges; the model-direction-aware topology used for network statistics contains 279 reaction-labelled edges collapsed to 263 directed simple edges. The single difference is FHL, which is chemically resolved but blocked by zero model bounds and is retained in the complete network representation.

## Scientific analyses

The analyses follow the organization of the manuscript.

### 1. Formalizing the electron economy

The repository reconstructs the *E. coli* electron-transfer network from reaction chemistry and evaluates reconstruction coverage, electron balance, donor–acceptor matching, independent chemical anchors, ambiguity, topology and robustness. These analyses support Fig. 1 and Supplementary Figs. S1–S2.

Principal entry points include `scripts/validate_reduction_degree.py`, `scripts/matching_validation.py`, `scripts/run_validation.py`, `scripts/validation_battery.py`, `scripts/run_analysis.py` and `scripts/coverage_analysis.py`.

### 2. Electron and carbon flows decouple

Aerobic and anaerobic glucose-grown *E. coli* states test whether carbon uptake determines electron flow. The analysis calculates net electron flux, reaction-level rerouting, carrier relay, effective transfer depth and terminal electron delivery. Independent Gonzalez <sup>13</sup>C-MFA flux maps test the same decoupling on glucose and xylose. ECOM4LA and respiratory-chain perturbations distinguish oxygen availability from respiratory disposal capacity. These analyses support Fig. 2 and Supplementary Figs. S3–S7.

Principal entry points include `scripts/run_electron_flow_analyses.py`, `scripts/run_electron_path_decomposition.py`, `scripts/run_gonzalez_13c_mfa_validation.py`, `data/experimental/portnoy2010_ecom4la/run_ecom4la_state_constrained.py` and `scripts/run_respiratory_architecture_support.py`.

### 3. A carbon source–disposal relation

Alternative terminal acceptors, six carbon sources and a six-substrate × seven-O₂-capacity factorial analysis separate the roles of carbon source and disposal capacity. Carbon source sets the baseline and maximal electron flow, whereas disposal capacity determines how much of that range is realized. Gonzalez, Denby, Toya and Perrenoud datasets provide complementary experimental checks while retaining the observable resolved by each dataset. These analyses support Fig. 3 and Supplementary Figs. S5 and S8–S12.

Principal entry points include `scripts/run_acceptor_experimental_validation.py`, `scripts/run_mechanistic_analysis.py` and `scripts/run_factorial_carbon_disposal_analysis.py`.

### 4. Shared principles across microorganisms

Five genome-scale reconstructions test the response to oxygen restriction and the separation between carbon-source-dependent baselines/maximum electron flux and disposal-capacity-dependent increases. Network-wide routing remains organism specific. Jouhten and Baumann <sup>13</sup>C-MFA datasets provide independent electron-generation analyses; Fromanger whole-cell balances provide complementary electron-equivalent accounting; and the Steinsiek and Weusthuis datasets test predictive transfer across respiratory architectures and carbon sources. Extended acceptor, carbon-source and routing analyses test robustness and generality. These analyses support Fig. 4 and Supplementary Figs. S13–S20.

Principal entry points include `scripts/run_cross_species_analysis.py`, `scripts/run_cross_species_principle_analysis.py`, `scripts/run_cross_species_routing_redistribution.py`, `scripts/run_jouhten_13c_mfa_validation.py`, `scripts/run_baumann2010_13c_mfa_validation.py`, `scripts/run_fromanger2010_electron_balance.py`, `scripts/run_steinsiek2014_predictive_transfer.py` and `scripts/run_weusthuis1994_predictive_transfer.py`.

## Repository structure

| Path | Contents |
| --- | --- |
| `etn/` | Core reconstruction, matching, electron accounting and network-analysis code |
| `data/` | Metabolic models, physiological constraints and static figure assets |
| `data/cross_species/` | Cross-species reconstructions, metadata and condition definitions |
| `data/experimental/` | Curated experimental inputs transcribed from published studies, with provenance |
| `scripts/` | Reconstruction, analysis, validation, Source Data and figure entry points |
| `results/` | Reconstruction, validation and calculated analysis outputs |
| `results/publication/` | Publication-facing numerical results |
| `results/cross_species/` | Cross-species physiological, routing and quality-control outputs |
| `source_data/publication/` | Panel-level Source Data for main and supplementary figures |
| `figures/publication/` | Computational main and supplementary figures in PNG and PDF formats; Supplementary Fig. S2 is also supplied as SVG |
| `tests/` | Unit, numerical, workflow and publication-contract tests |
| `REPRODUCIBILITY_WORKFLOW.md` | Staged reproduction commands and expected outputs |
| `RESULTS_MANIFEST.md` | Analysis-to-output inventory |
| `FIGURE_SOURCE_DATA_MANIFEST.csv` | Exact figure-to-Source-Data mapping |
| `INPUT_PROVENANCE_SHA256.csv` | Input-file provenance and SHA256 checksums |

## Installation

The reference environment uses Python 3.13.5. From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate  # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The `pyproject.toml` file also describes the reusable `etn` package and its minimum dependency versions. The pinned `requirements.txt` file is the reference environment for complete reproduction.

## Reproduction

### Complete workflow

Run the complete reconstruction, analysis and figure workflow from the repository root:

```bash
python scripts/run_all.py
```

The complete workflow generates reconstruction and validation outputs, publication analyses, cross-species results, Source Data, the supplementary Excel table, and all computational main and supplementary figures. It is intentionally comprehensive and may take substantially longer than an individual analysis.

Run the test suite separately after the workflow completes:

```bash
pytest -q
```

`run_all.py` does not invoke the test suite. The staged commands in `REPRODUCIBILITY_WORKFLOW.md` reproduce the same workflow while exposing progress by analysis group. Individual analyses can also be rerun using the entry points listed above.

### Generated outputs

- `results/`: reconstruction, validation and analysis-level outputs.
- `results/publication/`: numerical tables used by publication analyses.
- `source_data/publication/`: panel-level Source Data generated by `scripts/make_source_data.py`.
- `figures/publication/`: computational figures generated by `scripts/make_figures.py`; main-figure assembly for Figs. 2–4 is implemented in `scripts/main_figures.py`.
- `FIGURE_SOURCE_DATA_MANIFEST.csv`: machine-readable mapping from each figure to its Source Data files.

Existing generated outputs are included so that results can be inspected without rerunning the complete workflow.

## Experimental and model inputs

The input hierarchy and observable boundaries are documented in `data/README.md`. Briefly:

| Dataset | System | Observable used here | Analysis role |
| --- | --- | --- | --- |
| Chen et al. | *E. coli* | Glucose and O₂ uptake constraints | Aerobic/anaerobic model states |
| Gonzalez et al. | *E. coli* | <sup>13</sup>C-MFA reaction fluxes | Net electron-flux validation on glucose and xylose |
| Portnoy et al. | *E. coli* ECOM4LA | Physiological measurements and state-constrained model outputs | Oxygen availability versus respiratory disposal capacity |
| Steinsiek et al. | *E. coli* respiratory architectures | Acetate and ethanol production | Architecture controls and predictive transfer |
| Anand et al. | Synthetic respiratory chains | Replicated extracellular phenotypes | Respiratory-architecture support |
| Denby et al. | *E. coli* with TMAO | Reducing-equivalent balance and extracellular products | Non-O₂ terminal-acceptor validation |
| Toya et al. | *E. coli* with nitrate | <sup>13</sup>C-derived redox reaction fluxes | Electron-generation analysis |
| Perrenoud and Sauer | *E. coli* with nitrate or DMSO | Acceptor-specific central-carbon routing | Acceptor-specific routing check |
| Jouhten et al. | *S. cerevisiae* oxygen series | <sup>13</sup>C-MFA reaction fluxes | Electron-generation analysis and matched-model comparison |
| Baumann et al. | *K. phaffii* oxygen series | <sup>13</sup>C-MFA reaction fluxes | Electron-generation analysis |
| Fromanger et al. | *C. shehatae* | Whole-cell carbon and product balances | Electron-equivalent partition |
| Weusthuis et al. | *S. cerevisiae* glucose/maltose series | Ethanol production across oxygen supply | Carbon-source predictive transfer |

Dataset-specific README files document provenance, transformations, observable boundaries and primary publications. Input integrity is recorded in `INPUT_PROVENANCE_SHA256.csv` and dataset-specific checksum files.

## Figures and Source Data

`scripts/make_figures.py` generates computational main Figs. 2–4 and Supplementary Figs. S1–S20. Supplementary Fig. S2 is generated by `scripts/make_ecoli_etn_network_figure.py` from the reconstructed network tables. Main Figs. 2–4 are assembled by `scripts/main_figures.py`.

Panel-level numerical inputs are stored in `source_data/publication/` and mapped in `FIGURE_SOURCE_DATA_MANIFEST.csv`. Figure 1 is conceptual and therefore has no numerical Source Data. `source_data/publication/README.md` explains the observable boundaries and uncertainty conventions.

## Tests and validation

The test suite checks electron accounting, donor–acceptor matching, reconstruction coverage, calculations from experimental flux distributions, cross-species analyses, predictive-transfer analyses, figure contracts and manuscript-facing numerical claims. Run:

```bash
pytest -q
```

## Citation

Please cite the accompanying manuscript.

> Organizing principles of the cellular electron economy.
> Pierre Millard and Jean-Charles Portais
> bioRxiv preprint, 2026, doi: [XXXXXX](https://doi.org/XXXXXX)

## Licence

The repository code is distributed under the MIT License; see `LICENSE`. Individual bundled source materials remain subject to the terms identified in their dataset documentation.

## Author

Pierre Millard  
Toulouse Biotechnology Institute, INRAE  
millard@insa-toulouse.fr
