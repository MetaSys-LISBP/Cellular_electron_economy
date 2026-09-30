"""Condition-specific electron-flow analyses for aerobic and anaerobic glucose growth."""
import os, sys
REPO_ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),'..')); sys.path.insert(0,REPO_ROOT)
import pandas as pd
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from etn.data_acquisition import load_ecoli_model
from etn.matching import build_metabolite_table
from etn.fba_case_study import load_experimental_data, solve_row
from etn.electron_flow_analysis import (
    signed_donor_acceptor_edges, donor_acceptor_matrix,
    source_carrier_sink_classification, respiration_fermentation_biosynthesis_partition,
    redistribution, in_out_weighted_degree,
)

RESULTS_DIR=os.path.join(REPO_ROOT,'results')


def main():
    model=load_ecoli_model(); mt=build_metabolite_table(model)
    met_name={m.id:m.name for m in model.metabolites}
    data=load_experimental_data()
    edges={}; solutions={}
    for _,row in data.iterrows():
        sol=solve_row(model,row)
        if sol is None: raise RuntimeError(f"pFBA infeasible for {row.condition}")
        solutions[row.condition]=sol
        edges[row.condition]=signed_donor_acceptor_edges(model,mt,sol)

    for cond in ('aerobic','anaerobic'):
        mat=donor_acceptor_matrix(edges[cond])
        mat.to_csv(f'{RESULTS_DIR}/donor_acceptor_matrix_{cond}.csv',index=False)
        scs=source_carrier_sink_classification(edges[cond],met_name)
        scs.to_csv(f'{RESULTS_DIR}/source_carrier_sink_{cond}.csv',index=False)
        part=respiration_fermentation_biosynthesis_partition(edges[cond])
        part.to_csv(f'{RESULTS_DIR}/partition_{cond}.csv',index=False)
        deg=in_out_weighted_degree(edges[cond],met_name)
        deg.to_csv(f'{RESULTS_DIR}/in_out_degree_{cond}.csv',index=False)
        print(f'\n=== {cond}: top donor->acceptor transfers ===')
        print(mat.head(12).to_string(index=False))

    redist=redistribution(edges['aerobic'],edges['anaerobic'])
    redist.to_csv(f'{RESULTS_DIR}/redistribution_aerobic_vs_anaerobic.csv',index=False)
    print('\n=== Aerobic -> anaerobic redistribution (top 20 by |delta|) ===')
    print(redist.head(20).to_string(index=False))

    # Exchange and growth predictions for transparent auditing of the pFBA states.
    rows=[]
    for _,r in data.iterrows():
        sol=solutions[r.condition]
        rows.append({
            'condition':r.condition,'measured_qGLC':float(r.qGLC),'measured_qO2':float(r.qO2),
            'predicted_growth_rate':float(sol.fluxes['BIOMASS_Ec_iML1515_core_75p37M']),
            'predicted_acetate_exchange':float(sol.fluxes['EX_ac_e']),
            'predicted_formate_exchange':float(sol.fluxes['EX_for_e']),
            'predicted_ethanol_exchange':float(sol.fluxes['EX_etoh_e']),
            'predicted_succinate_exchange':float(sol.fluxes['EX_succ_e']),
            'predicted_lactate_D_exchange':float(sol.fluxes['EX_lac__D_e']),
            'predicted_CO2_exchange':float(sol.fluxes['EX_co2_e']),
        })
    pd.DataFrame(rows).to_csv(f'{RESULTS_DIR}/physiological_state_summary.csv',index=False)
    print('\nAll condition-specific electron-flow results written to',RESULTS_DIR)

if __name__=='__main__': main()
