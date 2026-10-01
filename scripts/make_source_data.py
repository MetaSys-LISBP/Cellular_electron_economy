#!/usr/bin/env python3
"""Generate machine-readable Source Data for the manuscript and Supplementary Information.

"""
from pathlib import Path
import shutil, json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "source_data" / "publication"
OUT.mkdir(parents=True, exist_ok=True)
for p in OUT.iterdir():
    if p.is_file(): p.unlink()

def cp(name, rel):
    src=ROOT/rel
    if not src.exists(): raise FileNotFoundError(src)
    shutil.copy2(src, OUT/name)

# Main quantitative figures (Fig. 1 is conceptual and has no numerical Source Data).
main = {
# Fig. 2 — glucose-focused discovery/causality
"Fig2a_condition_summary.csv":"results/fba_condition_summary.csv",
"Fig2a_electron_path_summary.csv":"results/electron_path_summary_by_condition.csv",
"Fig2b_carrier_relay.csv":"results/electron_carrier_relay_by_condition.csv",
"Fig2c_edge_changes.csv":"results/redistribution_aerobic_vs_anaerobic.csv",
"Fig2d_terminal_fates.csv":"results/electron_terminal_fates_interpreted_by_condition.csv",
"Fig2e_Gonzalez_13C_MFA_electron_flux_metrics.csv":"results/publication/gonzalez_13c_mfa_electron_flux_metrics.csv",
# Fig. 3 — E. coli carbon-entry and electron-disposal organization
"Fig3ab_acceptor_titration.csv":"results/publication/acceptor_electron_capacity_titration.csv",
"Fig3c_substrate_flux_compression.csv":"results/publication/substrate_flux_compression.csv",
"Fig3d_Gonzalez_13C_MFA_electron_flux_metrics.csv":"results/publication/gonzalez_13c_mfa_electron_flux_metrics.csv",
"Fig3e_substrate_core_redox_contribution.csv":"results/publication/substrate_core_redox_contribution.csv",
"Fig3f_factorial_carbon_oxygen_states.csv":"results/publication/factorial_carbon_oxygen_states.csv",
"Fig3f_inset_scaled_expansion_curves.csv":"results/publication/factorial_scaled_expansion_curves.csv",
# Fig. 4 — cross-species generalization
"Fig4a_cross_species_flux_compression.csv":"results/cross_species/electron_flux_compression_summary.csv",
"Fig4a_cross_species_physiology.csv":"results/cross_species/cross_species_physiology.csv",
"Fig4b_cross_species_carrier_relay.csv":"results/cross_species/cross_species_carrier_relay.csv",
"Fig4c_Jouhten2008_condition_mean_SD.csv":"results/publication/jouhten2008_13c_mfa_electron_source_metrics.csv",
"Fig4c_Jouhten2008_replicate_metrics.csv":"results/publication/jouhten2008_13c_mfa_replicate_metrics.csv",
"Fig4c_Jouhten2008_matched_model_experiment.csv":"results/publication/jouhten2008_matched_model_experiment.csv",
"Fig4d_cross_species_carbon_acceptor_pairs.csv":"results/publication/cross_species_carbon_acceptor_pairs.csv",
"Fig4d_cross_species_carbon_acceptor_feasibility.csv":"results/publication/cross_species_carbon_acceptor_feasibility.csv",
"Fig4e_Steinsiek2014_predictive_transfer.csv":"results/publication/steinsiek2014_predictive_transfer_points.csv",
"Fig4f_Weusthuis1994_predictive_transfer.csv":"results/publication/weusthuis1994_predictive_transfer_points.csv",
"Fig4f_Weusthuis1994_prediction_curve.csv":"results/publication/weusthuis1994_prediction_curve.csv",
}
for a,b in main.items(): cp(a,b)

# Factorial statistics used in the Supplementary factorial diagnostics.
fstats=json.loads((ROOT/'results/publication/factorial_carbon_oxygen_statistics.json').read_text())
pd.DataFrame([
    ('rank1_expansion_variance_pct',fstats['rank1_expansion_variance_pct'],'percent'),
    ('scaled_curve_min_pairwise_r',fstats['scaled_curve_min_pairwise_r'],'correlation'),
    ('scaled_curve_mean_pairwise_r',fstats['scaled_curve_mean_pairwise_r'],'correlation'),
    ('raw_substrate_pct',fstats['raw']['substrate_pct'],'percent variance'),
    ('raw_disposal_pct',fstats['raw']['oxygen_pct'],'percent variance'),
    ('raw_interaction_pct',fstats['raw']['interaction_pct'],'percent variance'),
    ('per_carbon_interaction_pct',fstats['per_carbon']['interaction_pct'],'percent variance'),
    ('per_gamma_interaction_pct',fstats['per_gamma']['interaction_pct'],'percent variance'),
],columns=['metric','value','unit']).to_csv(OUT/'SuppFigS12_factorial_statistics.csv',index=False)

# Supplementary quantitative panels.
supp = {
"SuppFigS2_iML1515_ETN_nodes.csv":"data/figure_assets/SuppFigS2_iML1515_ETN_nodes.csv",
"SuppFigS2_iML1515_ETN_edges.csv":"data/figure_assets/SuppFigS2_iML1515_ETN_edges.csv",
"SuppFigS3_electron_path_summary.csv":"results/electron_path_summary_by_condition.csv",
"SuppFigS3_physiological_state_summary.csv":"results/physiological_state_summary.csv",
"SuppFigS3_source_carrier_sink_aerobic.csv":"results/source_carrier_sink_aerobic.csv",
"SuppFigS3_source_carrier_sink_anaerobic.csv":"results/source_carrier_sink_anaerobic.csv",
"SuppFigS3_carrier_relay.csv":"results/electron_carrier_relay_by_condition.csv",
"SuppFigS13_cross_species_carrier_relay.csv":"results/cross_species/cross_species_carrier_relay.csv",
"SuppFigS8_acceptor_endpoints.csv":"results/publication/acceptor_endpoints.csv",
"SuppFigS4_Gonzalez_13C_MFA_electron_flux_metrics.csv":"results/publication/gonzalez_13c_mfa_electron_flux_metrics.csv",
"SuppFigS4_substrate_flux_compression.csv":"results/publication/substrate_flux_compression.csv",
"SuppFigS4_Gonzalez_normalization_summary.csv":"results/publication/gonzalez_13c_mfa_normalization_summary.csv",
"SuppFigS11_Gonzalez_anaerobic_source_decomposition.csv":"results/publication/gonzalez_anaerobic_source_decomposition.csv",
"SuppFigS7_Denby2015_model_matched_states.csv":"results/publication/denby2015_model_matched_states.csv",
"SuppFigS7_Denby2015_product_redistribution.csv":"results/publication/denby2015_product_redistribution.csv",
"SuppFigS10_Toya2012_source_totals.csv":"results/publication/toya2012_electron_source_totals.csv",
"SuppFigS10_Toya2012_source_decomposition.csv":"results/publication/toya2012_electron_source_decomposition.csv",
"SuppFigS9_Perrenoud2005_acceptor_routing.csv":"results/publication/perrenoud2005_acceptor_comparison.csv",
"SuppFigS16_cross_species_carbon_acceptor_pairs.csv":"results/publication/cross_species_carbon_acceptor_pairs.csv",
"SuppFigS16_cross_species_carbon_acceptor_feasibility.csv":"results/publication/cross_species_carbon_acceptor_feasibility.csv",
"SuppFigS16_Ecoli_core_redox_contribution.csv":"results/publication/substrate_core_redox_contribution.csv",
"SuppFigS16_normalization_summary.csv":"results/publication/cross_species_carbon_acceptor_normalization_summary.csv",
"SuppFigS18_Fromanger2010_electron_balance.csv":"results/publication/fromanger2010_electron_balance.csv",
"SuppFigS5_ECOM4LA_measured_endpoint_delivery.csv":"results/publication/ecom4la_measured_endpoint_electron_delivery.csv",
"SuppFigS5_ECOM4LA_state_constrained_summary.csv":"results/publication/ecom4la_state_constrained_summary.csv",
"SuppFigS5_ECOM4LA_published_13C_checks.csv":"results/publication/ecom4la_published_13C_checks.csv",
"SuppFigS5_ECOM4LA_constraint_robustness.csv":"results/publication/constraint_robustness.csv",
"SuppFigS5_substrate_flux_compression.csv":"results/publication/substrate_flux_compression.csv",
"SuppFigS6_Steinsiek2014.csv":"data/experimental/steinsiek2014/steinsiek2014_byproducts_means.csv",
"SuppFigS6_Anand2022.csv":"data/experimental/anand2022/anand2022_replicate_phenotypes.csv",
"SuppFigS12_factorial_states.csv":"results/publication/factorial_carbon_oxygen_states.csv",
"SuppFigS12_scaled_curves.csv":"results/publication/factorial_scaled_expansion_curves.csv",
"SuppFigS15_extended_acceptor_functional_states.csv":"results/publication/extended_acceptor_functional_states.csv",
"SuppFigS15_extended_acceptor_scaled_multiarchitecture.csv":"results/publication/extended_acceptor_scaled_multiarchitecture.csv",
"SuppFigS17_extended_carbon_entry_summary.csv":"results/publication/extended_carbon_entry_summary.csv",
"SuppFigS17_extended_carbon_entry_feasibility.csv":"results/publication/extended_carbon_entry_feasibility.csv",
"SuppFigS14_Baumann2010_oxygen_gradient.csv":"results/publication/baumann2010_13c_mfa_electron_source_metrics.csv",
"SuppFigS19_Steinsiek2014_predictive_transfer.csv":"results/publication/steinsiek2014_predictive_transfer_points.csv",
"SuppFigS19_Steinsiek2014_source_choice.csv":"results/publication/steinsiek2014_predictive_transfer_source_choice.csv",
"SuppFigS19_Steinsiek2014_source_level_omission.csv":"results/publication/steinsiek2014_predictive_transfer_source_level_omission.csv",
"SuppFigS20_Weusthuis1994_glucose_response.csv":"results/publication/weusthuis1994_glucose_response.csv",
"SuppFigS20_Weusthuis1994_maltose_predictions.csv":"results/publication/weusthuis1994_maltose_heldout_predictions.csv",
"SuppFigS20_Weusthuis1994_censoring_sensitivity.csv":"results/publication/weusthuis1994_censoring_sensitivity.csv",
}
for a,b in supp.items(): cp(a,b)


# Detailed reconstruction source data (Supplementary Fig. S1).
et = pd.read_csv(ROOT / "results/electron_transfer_table.tsv", sep="\t")
pairs = pd.read_csv(ROOT / "results/matched_pairs_all.csv")
hubs = pd.read_csv(ROOT / "results/top_hubs.csv").head(8).copy()
status_order = ["non_redox", "not_applicable", "redox", "unscorable", "ambiguous"]
labels = {"non_redox":"Non-redox","not_applicable":"Structural exclusion","redox":"Redox","unscorable":"Chemically unscorable","ambiguous":"Ambiguous"}
counts = et["redox_status"].value_counts().reindex(status_order, fill_value=0).astype(int)
pd.DataFrame({"category":[labels[x] for x in status_order],"redox_status":status_order,"n_reactions":counts.values}).to_csv(OUT / "SuppFigS1_classification_summary.csv", index=False)
stoich=(pd.to_numeric(et.loc[et.redox_status.eq("redox"),"n_electrons_transferred"],errors="raise").astype(int).value_counts().sort_index())
pd.DataFrame({"electrons_per_turnover":stoich.index.astype(int),"n_redox_reactions":stoich.values.astype(int)}).to_csv(OUT / "SuppFigS1_redox_stoichiometry.csv", index=False)
pd.DataFrame([
    ("Seed matching","curated_cofactor_list",int(pairs.source.eq("curated_cofactor_list").sum())),
    ("Formula / stoichiometry","ab_initio_formula_match",int(pairs.source.eq("ab_initio_formula_match").sum())),
],columns=["origin","source_code","n_matched_pair_instances"]).to_csv(OUT / "SuppFigS1_pair_origin_summary.csv", index=False)
hubs[["metabolite_id","name","weighted_degree_electrons"]].to_csv(OUT / "SuppFigS1_electron_weighted_hubs.csv", index=False)

# Denby observables used only in Supplementary Fig. S7.
denby=json.loads((ROOT/'results/publication/denby2015_summary.json').read_text())
pd.DataFrame([
    ("Measured TMAO terminal delivery", denby["measured_terminal_delivery_e_per_glucose"], denby["measured_terminal_delivery_e_sd_per_glucose"], "experimental"),
    ("Reported source generation", denby["reported_source_generation_e_per_glucose"], float("nan"), "experimental redox balance"),
    ("Model fermentation net flux", denby["model_fermentation_net_e_per_glucose"], float("nan"), "iML1515/ETN"),
    ("Model TMAO net flux", denby["model_TMAO_net_e_per_glucose"], float("nan"), "iML1515/ETN at measured TMAO/glucose ratio"),
], columns=["observable","e_per_glucose","sd_e_per_glucose","provenance"]).to_csv(OUT/'SuppFigS7_Denby2015_TMAO_observables.csv',index=False)

assert counts.to_dict() == {"non_redox":1783,"not_applicable":339,"redox":270,"unscorable":314,"ambiguous":6}
assert stoich.to_dict() == {1:2,2:237,4:22,6:6,8:3}

readme = """# Figure Source Data

This directory contains panel-level numerical Source Data generated by `scripts/make_source_data.py` from the repository analysis outputs. Exact mappings are listed in the repository-level `FIGURE_SOURCE_DATA_MANIFEST.csv`.

Figure 1 is conceptual and therefore has no numerical Source Data. Reconstruction diagnostics are reported in Supplementary Fig. S1.

## Main figures

- Fig. 2: `Fig2*` files — *E. coli* aerobic/anaerobic electron economy, routing and Gonzalez validation.
- Fig. 3: `Fig3*` files — terminal acceptors, carbon sources and the carbon source × disposal capacity response surface. Factorial statistics are supplied with Supplementary Fig. S12.
- Fig. 4a–d: `Fig4a_*` to `Fig4d_*` — cross-species analyses and Jouhten electron generation.
- Fig. 4e: `Fig4e_Steinsiek2014_predictive_transfer.csv` — respiratory-architecture predictive transfer.
- Fig. 4f: `Fig4f_Weusthuis1994_predictive_transfer.csv` and `Fig4f_Weusthuis1994_prediction_curve.csv` — carbon-source predictive transfer.

## Supplementary figures

Supplementary Fig. S2 presents the reconstructed iML1515 electron-transfer network and is supported by node and edge tables. Other files are named `SuppFigS#_*` according to their corresponding supplementary figure. Fromanger whole-cell electron-equivalent balances support Supplementary Fig. S18; respiratory-architecture transfer is documented in Supplementary Fig. S19; and carbon-source transfer is documented in Supplementary Fig. S20.

## Observable boundaries

- Gonzalez intervals are conservative envelopes obtained by independently propagating the reported marginal reaction-flux limits; they are not joint confidence intervals.
- Jouhten, Toya and Baumann files report electron generation calculated from electron-generating reactions in the published metabolic flux distributions. They do not report balanced genome-scale source-to-sink electron flux.
- Fromanger files report a whole-cell electron-equivalent partition calculated from published carbon yields and degree-of-reduction balance.
- ECOM4LA files distinguish measured physiological inputs and endpoint quantities from electron-flux quantities calculated using state-constrained model flux distributions.
- Steinsiek and Weusthuis analyses test response transfer using measured acetate or ethanol production as physiological observables.
- Model-only states are deterministic unless an uncertainty source is explicitly identified.

## Uncertainty conventions

- `s.d.` denotes reported or calculated sample standard deviation where biological replicates are available.
- Gonzalez uncertainty envelopes retain the published marginal reaction-flux limits.
- Toya and Baumann uncertainty values are calculated by independent propagation of reported marginal uncertainties because covariance matrices are unavailable.
- Deterministic genome-scale model calculations do not receive artificial error bars.
"""
(OUT/'README.md').write_text(readme)

# Machine-readable map from every figure-facing CSV to its repository path.
def natural_key(path):
    import re
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r'(\d+)', path.name)]

manifest_rows=[]
for path in sorted(OUT.glob('*.csv'), key=natural_key):
    stem=path.name.split('_',1)[0]
    panel_or_figure=(stem.replace('SuppFig','SuppFig') if stem.startswith('SuppFig') else stem[:4])
    manifest_rows.append((panel_or_figure,path.name,f"source_data/publication/{path.name}"))
pd.DataFrame(manifest_rows,columns=['panel_or_figure','file','repository_path']).to_csv(
    ROOT/'FIGURE_SOURCE_DATA_MANIFEST.csv',index=False
)
print(f"Generated {len(list(OUT.glob('*.csv')))} source-data CSVs in {OUT}")
