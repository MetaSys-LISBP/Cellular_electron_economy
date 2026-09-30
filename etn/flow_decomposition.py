"""Source-to-sink analysis of condition-specific electron-transfer flow.

The input is the reaction-resolved directional edge table returned by
``signed_donor_acceptor_edges``.  Edge weights are electron-transfer activity
(electron equivalents x mmol gDW^-1 h^-1).  Node balance (outflow - inflow)
defines net electron sources (>0), terminal sinks (<0), and balanced relay
nodes (=0) within the reaction-resolved electron-transfer layer.

Two classes of quantities are distinguished:

* decomposition-invariant quantities: total cumulative edge activity, net
  source/sink electron flux, effective transfer depth (activity / net source flux),
  source injection, terminal sink delivery, and node throughflow;
* illustrative path quantities: one deterministic source-to-sink path
  decomposition obtained by shortest-residual-path routing. Individual paths
  report the selected deterministic decomposition.

Generic unresolved partners are made reaction-specific before analysis so
unresolved chemistry from unrelated reactions cannot cancel artificially.
"""
from __future__ import annotations
from collections import defaultdict, deque
from typing import Dict, Iterable, List, Tuple
import math
import pandas as pd

from .constants import UNRESOLVED_ORGANIC_NODE


def reaction_specific_unresolved(edges: pd.DataFrame) -> pd.DataFrame:
    """Replace the generic unresolved node by a reaction-specific identity."""
    e = edges.copy()
    if e.empty:
        return e
    for side in ("donor", "acceptor"):
        mask = e[side].eq(UNRESOLVED_ORGANIC_NODE)
        e.loc[mask, side] = e.loc[mask, "reaction_id"].map(lambda r: f"unresolved:{r}")
    return e


def aggregate_directed_flow(edges: pd.DataFrame, min_flow: float = 1e-12) -> Dict[Tuple[str, str], float]:
    """Aggregate reaction-resolved edges into a positive donor->acceptor map."""
    if edges.empty:
        return {}
    g = edges.groupby(["donor", "acceptor"], as_index=False)["electron_flux"].sum()
    return {(str(r.donor), str(r.acceptor)): float(r.electron_flux)
            for r in g.itertuples() if float(r.electron_flux) > min_flow and r.donor != r.acceptor}


def node_flow_statistics(flow: Dict[Tuple[str, str], float]) -> pd.DataFrame:
    """Return inflow, outflow, net balance and relay throughflow for all nodes.

    ``balance = outflow - inflow``.  Positive balance is net source injection;
    negative balance is net sink delivery.  ``relay_throughflow=min(in,out)``
    is decomposition-invariant and measures how much electron flow is relayed
    by a node rather than injected/removed there.
    """
    inflow = defaultdict(float); outflow = defaultdict(float)
    for (u, v), f in flow.items():
        outflow[u] += f; inflow[v] += f
    nodes = sorted(set(inflow) | set(outflow))
    rows=[]
    for n in nodes:
        i=float(inflow[n]); o=float(outflow[n]); b=o-i
        rows.append({"node":n,"inflow":i,"outflow":o,"balance":b,
                     "relay_throughflow":min(i,o)})
    return pd.DataFrame(rows)


def _shortest_residual_path(residual: Dict[Tuple[str, str], float], source: str,
                            sinks: set[str], tol: float) -> List[str] | None:
    adj=defaultdict(list)
    for (u,v),f in residual.items():
        if f > tol: adj[u].append(v)
    for u in adj: adj[u].sort()
    q=deque([source]); parent={source:None}; target=None
    while q:
        u=q.popleft()
        if u in sinks and u != source:
            target=u; break
        for v in adj.get(u,[]):
            if v not in parent:
                parent[v]=u; q.append(v)
    if target is None: return None
    path=[]; x=target
    while x is not None:
        path.append(x); x=parent[x]
    return list(reversed(path))


def decompose_source_sink_paths(edges: pd.DataFrame, tol: float = 1e-9) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    """Return one deterministic path decomposition, invariant summary, node stats."""
    e=reaction_specific_unresolved(edges)
    original=aggregate_directed_flow(e,min_flow=tol*1e-3)
    stats=node_flow_statistics(original)
    supply={r.node:float(r.balance) for r in stats.itertuples() if r.balance > tol}
    demand={r.node:float(-r.balance) for r in stats.itertuples() if r.balance < -tol}
    total_supply=sum(supply.values()); total_demand=sum(demand.values())
    if abs(total_supply-total_demand) > max(tol,1e-8*max(total_supply,total_demand,1.0)):
        raise ValueError(f"Reaction-resolved electron flow is imbalanced: {total_supply} vs {total_demand}")

    residual=dict(original); rows=[]; pid=0
    for source in sorted(supply):
        while supply[source] > tol:
            sinks={n for n,d in demand.items() if d > tol}
            path=_shortest_residual_path(residual,source,sinks,tol)
            if path is None:
                raise ValueError(f"Cannot route residual source flow from {source}")
            sink=path[-1]
            amount=min(supply[source],demand[sink],*(residual[(u,v)] for u,v in zip(path[:-1],path[1:])))
            if amount <= tol: raise RuntimeError("Non-positive path amount")
            for u,v in zip(path[:-1],path[1:]):
                residual[(u,v)]-=amount
                if residual[(u,v)] < tol: residual[(u,v)]=0.0
            supply[source]-=amount; demand[sink]-=amount
            rows.append({"path_id":pid,"source":source,"sink":sink,"path":" -> ".join(path),
                         "n_transfer_steps":len(path)-1,"electron_flux":amount})
            pid+=1

    paths=pd.DataFrame(rows)
    total_activity=float(sum(original.values()))
    net_flux=float(total_supply)
    # This ratio is independent of how paths/cycles are decomposed and has a
    # direct interpretation as cumulative transfer events per net e-equivalent.
    depth=total_activity/net_flux if net_flux > 0 else float('nan')
    # Report residual edge activity as a diagnostic of the selected path and
    # cycle decomposition.
    residual_activity=float(sum(max(0.0,f) for f in residual.values()))
    summary={
        "total_edge_activity":total_activity,
        "net_source_flux":net_flux,
        "net_sink_flux":float(total_demand),
        "effective_transfer_depth":depth,
        "illustrative_n_paths":int(len(paths)),
        "illustrative_residual_activity":residual_activity,
    }
    return paths,summary,stats
