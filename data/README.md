# Input data

This directory contains the metabolic reconstructions, physiological constraints, experimental metabolic flux distributions and static figure assets used by the repository. Published studies provide metabolic flux distributions or physiological measurements; electron fluxes and electron-generation quantities are calculated within this repository.

## Directory structure

| Path | Contents |
| --- | --- |
| `iML1515.xml` | *Escherichia coli* K-12 MG1655 genome-scale metabolic reconstruction used for the principal analyses |
| `aerobic_anaerobic_glucose_constraints.tsv` | Measured glucose and O₂ constraints for the principal aerobic/anaerobic *E. coli* comparison |
| `cross_species/` | Additional organism reconstructions, metadata, condition definitions and yeast anaerobic-model adjustments |
| `experimental/` | Curated experimental inputs transcribed from published studies |
| `figure_assets/` | Static conceptual or network assets used during figure assembly |

Input-file integrity is recorded in the repository-level `INPUT_PROVENANCE_SHA256.csv` and in dataset-specific checksum files where applicable.

## Principal *E. coli* reconstruction and constraints

`iML1515.xml` is the *E. coli* K-12 MG1655 genome-scale metabolic reconstruction used for electron-transfer-network reconstruction and most *E. coli* model analyses.

`aerobic_anaerobic_glucose_constraints.tsv` contains external flux constraints for wild-type *E. coli* K-12 MG1655 grown in M9 minimal medium with glucose as the sole carbon source. The biological measurements are from Chen X, Alonso AP, Allen DK, Reed JL & Shachar-Hill Y, *Metabolic Engineering* 13, 38–48 (2011), DOI: 10.1016/j.ymben.2010.11.004. The exact glucose and oxygen uptake constraints used here are 8.7 and 11.9 mmol gDW⁻¹ h⁻¹ aerobically and 14.9 and 0 mmol gDW⁻¹ h⁻¹ anaerobically.

Only glucose and oxygen uptake are used as quantitative measured constraints. Aerobic lactate, succinate, formate and ethanol exchanges are fixed to zero from the reported non-detections. Growth and all other unconstrained exchange rates are model outputs.

## Cross-species reconstructions

`cross_species/` contains the additional reconstructions and condition definitions used for *Bacillus subtilis*, *Salmonella enterica*, *Komagataella phaffii* and *Saccharomyces cerevisiae*. The directory also contains an independent *E. coli* reconstruction used for robustness analysis.

See `cross_species/README.md` and `cross_species/model_metadata.tsv` for model identities, source information and organism-specific configuration. `cross_species/condition_config.tsv` records the condition definitions used in the standardized comparisons. Infeasible model–substrate combinations are retained in the generated outputs rather than silently omitted.

## Experimental inputs

| Directory | System and input | Quantity calculated or analysis supported |
| --- | --- | --- |
| `experimental/gonzalez2017/` | *E. coli* glucose/xylose × aerobic/anaerobic <sup>13</sup>C-MFA reaction-flux distributions | Net electron flux and reaction-level source decomposition |
| `experimental/portnoy2010_ecom4la/` | *E. coli* WT and ECOM4LA physiological measurements and published <sup>13</sup>C checks | State-constrained electron flux and respiratory disposal-capacity analysis |
| `experimental/steinsiek2014/` | Acetate and ethanol production across *E. coli* respiratory architectures | Architecture controls and predictive transfer |
| `experimental/anand2022/` | Replicated synthetic respiratory-chain phenotypes | Respiratory-architecture support |
| `experimental/denby2015/` | TMAO perturbation, product measurements and reducing-equivalent balance | Terminal electron delivery and matched-model analysis |
| `experimental/toya2012/` | *E. coli* anaerobic/nitrate <sup>13</sup>C-derived reaction-flux distributions | Electron generation and pathway decomposition |
| `experimental/perrenoud2005/` | Nitrate-versus-DMSO central-carbon routing | Acceptor-specific routing comparison |
| `experimental/jouhten2008/` | *S. cerevisiae* oxygen-gradient <sup>13</sup>C-MFA reaction-flux distributions | Electron generation and matched-model comparison |
| `experimental/baumann2010/` | *K. phaffii* oxygen-gradient <sup>13</sup>C-MFA reaction-flux distributions | Electron generation |
| `experimental/fromanger2010/` | *C. shehatae* whole-cell carbon yields and oxygen uptake | Electron-equivalent partition |
| `experimental/weusthuis1994/` | *S. cerevisiae* glucose/maltose oxygen series | Carbon-source predictive transfer |

Each dataset README identifies the primary publication, files used, transformations performed and observable boundary. The experimental values are not represented as direct measurements of electron flux. Electron quantities are calculated by the repository from the supplied metabolic flux distributions, extracellular measurements or whole-cell balances, as appropriate for each dataset.

## Observable boundaries

- Gonzalez flux distributions support balanced source-to-sink electron-flux calculations within the supplied formulation.
- Jouhten, Toya and Baumann flux distributions resolve electron generation from reported central-carbon reactions but do not provide complete source-to-sink electron balances.
- Denby and Fromanger provide endpoint or whole-cell electron-equivalent balances rather than reaction-resolved genome-scale electron flux.
- Steinsiek, Anand and Weusthuis provide extracellular physiological observables used for architecture controls or predictive-transfer analyses.
- ECOM4LA measured constraints are completed with a state-constrained genome-scale model; measured and model-derived quantities remain identified separately.
