# Portnoy et al. (2010) ECOM4LA respiratory-disposal analysis

Primary source: Portnoy VA et al. **Deletion of genes encoding cytochrome oxidases and quinol monooxygenase blocks the aerobic–anaerobic shift in Escherichia coli K-12 MG1655.** *Applied and Environmental Microbiology* 76, 6529–6540 (2010).

Portnoy et al. compared oxic wild-type *E. coli* MG1655 with ECOM4LA, an evolved strain lacking the major terminal cytochrome oxidase systems together with `ygiN`. Oxygen is present in both conditions, but oxygen utilization is nearly abolished in ECOM4LA. This provides a direct test of oxygen availability versus respiratory disposal capacity.

`run_ecom4la_state_constrained.py` fixes the reported glucose, oxygen, D-lactate and acetate rates in iML1515, maximizes growth, applies pFBA and calculates electron fluxes from the resulting balanced model-derived metabolic flux distributions using the reconstructed electron-transfer stoichiometries.

The directory contains:

- measured physiological constraints and direct endpoint electron-delivery calculations;
- state-constrained electron-flux, sink and path outputs;
- published ^13C-derived oxidative checks;
- strict and relaxed product-secretion robustness calculations;
- edge and illustrative-path tables for the wild-type and ECOM4LA states.

Measured quantities, direct endpoint calculations and genome-scale model-derived electron-flux quantities remain explicitly distinct. `results/publication/ecom4la_state_constrained_summary.csv` contains the principal state summary, and the corresponding figure-facing files support Supplementary Fig. S5.
