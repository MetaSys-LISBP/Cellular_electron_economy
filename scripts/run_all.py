"""Run the complete analysis workflow.

This convenience runner is intentionally comprehensive and may take longer than
interactive front-end timeouts. The staged commands in
REPRODUCIBILITY_WORKFLOW.md provide the same reproducible workflow.
"""
from pathlib import Path
import subprocess, sys
ROOT=Path(__file__).resolve().parents[1]

def run(rel,*args):
    cmd=[sys.executable,str(ROOT/rel),*args]
    print('\n===',' '.join([rel,*args]),'===',flush=True)
    subprocess.run(cmd,cwd=ROOT,check=True)

def main():
    for s in [
        'scripts/validate_reduction_degree.py','scripts/matching_validation.py','scripts/run_validation.py','scripts/validation_battery.py','scripts/run_analysis.py','scripts/coverage_analysis.py',
        'scripts/run_electron_flow_analyses.py','scripts/run_electron_path_decomposition.py','scripts/compare_networks.py','scripts/analyse_community_resolution_sensitivity.py','scripts/analyse_community_robustness.py','scripts/analyse_network_comparison_sensitivity.py','scripts/analyse_flux_robustness.py',
        'scripts/run_cross_species_analysis.py','scripts/run_cross_species_routing_redistribution.py',
        'scripts/run_gonzalez_13c_mfa_validation.py','scripts/run_jouhten_13c_mfa_validation.py','scripts/run_jouhten_matched_model_validation.py',
        'scripts/run_acceptor_experimental_validation.py','scripts/run_mechanistic_analysis.py',
        'data/experimental/portnoy2010_ecom4la/run_ecom4la_state_constrained.py','scripts/run_respiratory_architecture_support.py',
        'scripts/run_factorial_carbon_disposal_analysis.py','scripts/run_extended_acceptor_screen.py','scripts/run_extended_carbon_entry_screen.py']:
        run(s)
    # Staged cross-species principle route: equivalent canonical outputs, clearer progress.
    for org in ['B. subtilis','S. enterica','K. phaffii']:
        run('scripts/run_cross_species_principle_analysis.py','--organism',org)
    run('scripts/run_cross_species_principle_analysis.py','--merge')
    run('scripts/run_fromanger2010_electron_balance.py')
    run('scripts/run_baumann2010_13c_mfa_validation.py')
    run('scripts/run_steinsiek2014_predictive_transfer.py')
    run('scripts/run_weusthuis1994_predictive_transfer.py')
    for s in ['scripts/make_excel_table.py','scripts/make_source_data.py','scripts/make_figures.py']:
        run(s)
    print('\nComplete analysis workflow finished successfully.',flush=True)
if __name__=='__main__': main()
