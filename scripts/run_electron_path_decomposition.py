"""Source-to-sink electron routing for aerobic and anaerobic pFBA states.

The reported net source flux, terminal sink delivery, relay throughflow and
effective transfer depth are invariant to the non-unique path decomposition.
A deterministic path decomposition is exported only for traceability.
"""
import os, sys
REPO_ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),'..')); sys.path.insert(0,REPO_ROOT)
import pandas as pd
from etn.data_acquisition import load_ecoli_model
from etn.matching import build_metabolite_table
from etn.fba_case_study import load_experimental_data, solve_row
from etn.electron_flow_analysis import signed_donor_acceptor_edges, RESPIRATORY_ACCEPTORS, FERMENTATION_PRODUCTS
from etn.flow_decomposition import decompose_source_sink_paths

RESULTS=os.path.join(REPO_ROOT,'results')
CARRIERS={'nadh_c':'NADH','nadph_c':'NADPH','fadh2_c':'FADH2','q8h2_c':'Ubiquinol-8',
          'mql8_c':'Menaquinol-8','2dmmql8_c':'2-Demethylmenaquinol-8'}
FATES=['respiration','fermentation','biosynthesis_or_other','unresolved_sink']
INTERPRETED_FATES=['respiration','fermentation_associated','biosynthesis_or_other','unresolved_other']

def sink_fate(s):
    if s in RESPIRATORY_ACCEPTORS: return 'respiration'
    if s in FERMENTATION_PRODUCTS: return 'fermentation'
    if str(s).startswith('unresolved:'): return 'unresolved_sink'
    return 'biosynthesis_or_other'

def main():
    model=load_ecoli_model(); mt=build_metabolite_table(model); data=load_experimental_data()
    all_paths=[]; summaries=[]; nodes_all=[]; sink_rows=[]; source_rows=[]; carrier_rows=[]
    for _,row in data.iterrows():
        cond=row.condition; sol=solve_row(model,row)
        if sol is None: raise RuntimeError(f'pFBA infeasible for {cond}')
        edges=signed_donor_acceptor_edges(model,mt,sol)
        paths,summary,node_stats=decompose_source_sink_paths(edges)
        paths['condition']=cond; all_paths.append(paths)
        summaries.append(dict(condition=cond,**summary))
        node_stats['condition']=cond; nodes_all.append(node_stats)
        for r in node_stats.itertuples():
            if r.balance < -1e-9:
                sink_rows.append({'condition':cond,'sink':r.node,'sink_fate':sink_fate(r.node),
                                  'electron_flux':-r.balance})
            elif r.balance > 1e-9:
                source_rows.append({'condition':cond,'source':r.node,
                                    'source_type':'unresolved_reaction_source' if str(r.node).startswith('unresolved:') else 'resolved_metabolite_source',
                                    'electron_flux':r.balance})
        nd=node_stats.set_index('node')
        for cid,cname in CARRIERS.items():
            carrier_rows.append({'condition':cond,'carrier':cid,'carrier_name':cname,
                                 'relay_throughflow':float(nd.loc[cid,'relay_throughflow']) if cid in nd.index else 0.0})

    pd.concat(all_paths,ignore_index=True).to_csv(os.path.join(RESULTS,'electron_source_sink_paths_by_condition.csv'),index=False)
    summ=pd.DataFrame(summaries); summ.to_csv(os.path.join(RESULTS,'electron_path_summary_by_condition.csv'),index=False)
    pd.concat(nodes_all,ignore_index=True).to_csv(os.path.join(RESULTS,'electron_node_flow_statistics_by_condition.csv'),index=False)
    sinks=pd.DataFrame(sink_rows); sinks.to_csv(os.path.join(RESULTS,'electron_sink_delivery_by_condition.csv'),index=False)
    pd.DataFrame(source_rows).to_csv(os.path.join(RESULTS,'electron_source_injection_by_condition.csv'),index=False)
    carriers=pd.DataFrame(carrier_rows); carriers.to_csv(os.path.join(RESULTS,'electron_carrier_relay_by_condition.csv'),index=False)

    sf=sinks.groupby(['condition','sink_fate'],as_index=False).electron_flux.sum()
    grid=pd.MultiIndex.from_product([['aerobic','anaerobic'],FATES],names=['condition','sink_fate']).to_frame(index=False)
    sf=grid.merge(sf,on=['condition','sink_fate'],how='left').fillna({'electron_flux':0.0})
    sf['total_terminal_flux']=sf.groupby('condition').electron_flux.transform('sum')
    sf['fraction']=sf['electron_flux']/sf['total_terminal_flux']
    sf.to_csv(os.path.join(RESULTS,'electron_terminal_fates_by_condition.csv'),index=False)

    # Interpretation layer for the aerobic/anaerobic case study. The reverse
    # ACALD reaction consumes NADH to reduce acetyl-CoA toward acetaldehyde,
    # but group transfer prevents a molecularly resolved organic acceptor. Its
    # reaction-specific unresolved sink is therefore reported transparently as
    # fermentation-associated rather than mixed with unrelated unresolved
    # biosynthetic sinks. This does not alter reconstruction or edge inference.
    interp=sinks.copy()
    def interpreted_fate(row):
        if row.sink_fate == 'respiration': return 'respiration'
        if row.sink_fate == 'fermentation': return 'fermentation_associated'
        if row.sink_fate == 'biosynthesis_or_other': return 'biosynthesis_or_other'
        if str(row.sink) == 'unresolved:ACALD': return 'fermentation_associated'
        return 'unresolved_other'
    interp['interpreted_fate']=interp.apply(interpreted_fate,axis=1)
    it=interp.groupby(['condition','interpreted_fate'],as_index=False).electron_flux.sum()
    igrid=pd.MultiIndex.from_product([['aerobic','anaerobic'],INTERPRETED_FATES],names=['condition','interpreted_fate']).to_frame(index=False)
    it=igrid.merge(it,on=['condition','interpreted_fate'],how='left').fillna({'electron_flux':0.0})
    it['total_terminal_flux']=it.groupby('condition').electron_flux.transform('sum')
    it['fraction']=it['electron_flux']/it['total_terminal_flux']
    it.to_csv(os.path.join(RESULTS,'electron_terminal_fates_interpreted_by_condition.csv'),index=False)

    print('\n=== Source-to-sink invariant summary ==='); print(summ.to_string(index=False))
    print('\n=== Terminal electron fates ==='); print(sf.to_string(index=False))
    print('\n=== Interpreted terminal electron fates ==='); print(it.to_string(index=False))
    print('\n=== Carrier relay throughflow ==='); print(carriers.to_string(index=False))

if __name__=='__main__': main()
