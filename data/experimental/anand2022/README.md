# Anand et al. (2022) synthetic respiratory-chain phenotypes

Primary source: Anand A et al. **Laboratory evolution of synthetic electron transport system variants reveals a larger metabolic respiratory system and its plasticity.** *Nature Communications* 13, 3682 (2022).

`anand2022_replicate_phenotypes.csv` contains the replicated growth, glucose-uptake and acetate-secretion values used for the synthetic unbranched electron-transport-system architectures before and after adaptive evolution.

The repository calculates the fraction of consumed glucose carbon secreted as acetate from the reported extracellular rates. These data provide a physiological respiratory-architecture control; they are not treated as direct electron-flux measurements.

`scripts/run_respiratory_architecture_support.py` calculates the condition means and sample standard deviations stored in `results/publication/anand2022_architecture_summary.csv`. The panel-level values supporting Supplementary Fig. S7 are copied to `source_data/publication/SuppFigS7_Anand2022.csv`.
