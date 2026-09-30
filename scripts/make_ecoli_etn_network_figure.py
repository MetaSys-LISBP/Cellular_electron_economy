#!/usr/bin/env python3
"""Generate the boundary-resolved E. coli iML1515 electron-transfer network.

The layout uses automated source and sink classification, compact boundary
columns, collision avoidance and outward label placement.

Inputs
------
1. electron_transfer_table.tsv
2. iML1515.xml (SBML model)

Outputs
-------
PNG, PDF, SVG, plus node/edge/classification tables.
"""

from __future__ import annotations

import argparse
import math
import re
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch


SPECIAL_LABELS = {
    "nadh": "NADH",
    "nadph": "NADPH",
    "fadh2": "FADH₂",
    "q8h2": "Ubiquinol-8",
    "mql8": "Menaquinol-8",
    "2dmmql8": "DMK-8H₂",
    "h2o": "H₂O",
    "h2o2": "H₂O₂",
    "h2s": "H₂S",
    "nh4": "NH₄⁺",
    "h2": "H₂",
    "etoh": "Ethanol",
    "succ": "Succinate",
    "g6p": "G6P",
    "3pg": "3PG",
    "mal__L": "Malate",
    "asp__L": "Aspartate",
    "glyc3p": "Glycerol-3-P",
    "glyclt": "Glycolate",
    "mlthf": "5,10-CH₂-THF",
    "thf": "THF",
    "5mthf": "5mTHF",
    "e4p": "E4P",
    "ipdp": "IPDP",
    "dmpp": "DMPP",
    "dhor__S": "Dihydroorotate",
    "23dhmb": "23dhmb",
    "23dhmp": "23dhmp",
    "23ddhb": "23ddhb",
    "thdp": "THDP",
    "hom__L": "Homoserine",
    "skm": "Shikimate",
    "3c2hmp": "3c2hmp",
    "2oph": "2-Octaprenylphenol",
    "pdx5p": "pdx5p",
    "imp": "IMP",
    "histd": "Histidinol",
    "pppg9": "pppg9",
    "dscl": "dscl",
    "for": "Formate",
    "lac__D": "D-Lactate",
    "lac__L": "L-Lactate",
    "glc__D": "D-Glucose",
    "acald": "Acetaldehyde",
}


def split_field(value):
    if pd.isna(value):
        return []
    return [x.strip() for x in str(value).split(";") if x.strip()]


def parse_bool(value):
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def chemical_id(met_id):
    return re.sub(r"_(c|p|e)$", "", str(met_id))


def name_without_compartment(name):
    if name is None or pd.isna(name):
        return ""
    return str(name).replace(" H2O H2O", "")


def short_label(chem_id, candidate_name=None):
    if chem_id in SPECIAL_LABELS:
        return SPECIAL_LABELS[chem_id]
    nm = name_without_compartment(candidate_name)
    if nm and len(nm) <= 20:
        return nm
    lab = chem_id.replace("__", "-")
    return lab if len(lab) <= 20 else lab[:18] + "…"


def reconstruct_reaction_edges(table):
    resolved = table[
        table["redox_status"].astype(str).str.lower().eq("redox")
        & table["resolution_status"].astype(str).eq("fully_resolved")
    ].copy()

    rows = []
    names = {}

    for r in resolved.itertuples(index=False):
        donors = split_field(r.donor_metabolite_ids)
        acceptors = split_field(r.acceptor_metabolite_ids)
        donor_names = split_field(r.donor_metabolite_names)
        acceptor_names = split_field(r.acceptor_metabolite_names)

        for i, met in enumerate(donors):
            names.setdefault(met, donor_names[i] if i < len(donor_names) else met)
        for i, met in enumerate(acceptors):
            names.setdefault(met, acceptor_names[i] if i < len(acceptor_names) else met)

        if not donors or not acceptors:
            continue

        n_pairs = len(donors) * len(acceptors)
        edge_e = float(r.n_electrons_transferred) / n_pairs if n_pairs else np.nan

        for donor in donors:
            for acceptor in acceptors:
                rows.append(
                    {
                        "donor_original": donor,
                        "acceptor_original": acceptor,
                        "reaction_id": r.reaction_id,
                        "reaction_name": r.reaction_name,
                        "direction": "forward",
                        "electron_weight": edge_e,
                    }
                )
                if parse_bool(r.reversible):
                    rows.append(
                        {
                            "donor_original": acceptor,
                            "acceptor_original": donor,
                            "reaction_id": r.reaction_id,
                            "reaction_name": r.reaction_name,
                            "direction": "reverse",
                            "electron_weight": edge_e,
                        }
                    )
    return pd.DataFrame(rows), names


def find_exchangeable_chemical_ids(sbml_path):
    tree = ET.parse(sbml_path)
    root = tree.getroot()
    core_ns = root.tag.split("}")[0].strip("{")
    ns = {"s": core_ns}

    exchangeable = set()
    for rxn in root.findall(".//s:listOfReactions/s:reaction", ns):
        rid = rxn.attrib.get("id", "")
        sbo = rxn.attrib.get("sboTerm", "")
        if not (rid.startswith("R_EX_") or rid.startswith("EX_") or sbo == "SBO:0000627"):
            continue
        species_refs = (
            rxn.findall("./s:listOfReactants/s:speciesReference", ns)
            + rxn.findall("./s:listOfProducts/s:speciesReference", ns)
        )
        for ref in species_refs:
            sid = ref.attrib.get("species", "")
            if sid.startswith("M_"):
                sid = sid[2:]
            if sid.endswith("_e"):
                exchangeable.add(chemical_id(sid))
    return exchangeable


def representative_names(original_names):
    by_chem = defaultdict(list)
    for met, name in original_names.items():
        by_chem[chemical_id(met)].append(name)
    result = {}
    for chem, names in by_chem.items():
        names = [name_without_compartment(n) for n in names if n]
        result[chem] = min(names, key=len) if names else chem
    return result


def classify_and_remap(edges_original, original_names, exchangeable):
    edges = edges_original.copy()
    edges["donor_chem"] = edges["donor_original"].map(chemical_id)
    edges["acceptor_chem"] = edges["acceptor_original"].map(chemical_id)

    donor_count = Counter(edges["donor_chem"])
    acceptor_count = Counter(edges["acceptor_chem"])
    all_chems = sorted(set(donor_count) | set(acceptor_count))
    rep_names = representative_names(original_names)

    rows = []
    display_map = {}
    for chem in all_chems:
        is_boundary = chem in exchangeable
        has_out = donor_count[chem] > 0
        has_in = acceptor_count[chem] > 0

        if not is_boundary:
            node_id = f"{chem}::carrier"
            display_map[(chem, "donor")] = node_id
            display_map[(chem, "acceptor")] = node_id
            rows.append(
                {
                    "node_id": node_id,
                    "chemical_id": chem,
                    "name": rep_names.get(chem, chem),
                    "role": "carrier",
                    "boundary_accessible": False,
                    "donor_edge_count": donor_count[chem],
                    "acceptor_edge_count": acceptor_count[chem],
                }
            )
            continue

        if has_out:
            source_id = f"{chem}::source"
            display_map[(chem, "donor")] = source_id
            rows.append(
                {
                    "node_id": source_id,
                    "chemical_id": chem,
                    "name": rep_names.get(chem, chem),
                    "role": "source",
                    "boundary_accessible": True,
                    "donor_edge_count": donor_count[chem],
                    "acceptor_edge_count": 0,
                }
            )

        if has_in:
            sink_id = f"{chem}::sink"
            display_map[(chem, "acceptor")] = sink_id
            rows.append(
                {
                    "node_id": sink_id,
                    "chemical_id": chem,
                    "name": rep_names.get(chem, chem),
                    "role": "terminal_sink",
                    "boundary_accessible": True,
                    "donor_edge_count": 0,
                    "acceptor_edge_count": acceptor_count[chem],
                }
            )

    nodes = pd.DataFrame(rows)
    remapped = []
    for r in edges.itertuples(index=False):
        remapped.append(
            {
                "source_node": display_map[(r.donor_chem, "donor")],
                "target_node": display_map[(r.acceptor_chem, "acceptor")],
                "donor_original": r.donor_original,
                "acceptor_original": r.acceptor_original,
                "donor_chemical_id": r.donor_chem,
                "acceptor_chemical_id": r.acceptor_chem,
                "reaction_id": r.reaction_id,
                "reaction_name": r.reaction_name,
                "direction": r.direction,
                "electron_weight": r.electron_weight,
            }
        )
    display_edges = pd.DataFrame(remapped)

    in_degree = Counter(display_edges["target_node"])
    out_degree = Counter(display_edges["source_node"])
    nodes["in_degree"] = nodes["node_id"].map(lambda n: in_degree[n])
    nodes["out_degree"] = nodes["node_id"].map(lambda n: out_degree[n])
    nodes["degree"] = nodes["in_degree"] + nodes["out_degree"]
    nodes["label"] = nodes.apply(lambda r: short_label(r["chemical_id"], r["name"]), axis=1)
    return nodes, display_edges


def node_size(degree):
    # Binned connectivity scale; deliberately compressed to preserve readability.
    if degree >= 100:
        return 1500
    if degree >= 50:
        return 1200
    if degree >= 25:
        return 900
    if degree >= 12:
        return 680
    if degree >= 6:
        return 510
    if degree >= 3:
        return 390
    return 290


def robust_rescale(values, half_range):
    values = np.asarray(values, dtype=float)
    if len(values) == 0:
        return values
    lo, hi = np.quantile(values, [0.01, 0.99])
    if hi <= lo:
        return np.zeros_like(values)
    values = np.clip(values, lo, hi)
    return half_range * (2 * (values - lo) / (hi - lo) - 1)


def spread_two_columns(ids, spring, x_inner, x_outer, y_low=-0.93, y_high=0.93):
    """Two compact boundary columns, vertically staggered for label readability."""
    ids = sorted(ids, key=lambda n: spring[n][1], reverse=True)
    cols = [ids[::2], ids[1::2]]
    pos = {}
    nmax = max(len(cols[0]), len(cols[1]), 1)
    half_step = (y_high - y_low) / max(2 * nmax, 1)
    for j, (x, col) in enumerate(zip((x_inner, x_outer), cols)):
        if not col:
            continue
        top = y_high - (half_step if j else 0.0)
        bottom = y_low + (half_step if j else 0.0)
        ys = np.linspace(top, bottom, len(col))
        for n, y in zip(col, ys):
            pos[n] = np.array([x, y], dtype=float)
    return pos


def marker_radius_data(size, x_span=2.52, y_span=2.08, fig_w=20.0, fig_h=15.0):
    """Approximate scatter-marker radius in data units from size in pt^2."""
    r_pt = math.sqrt(size / math.pi)
    r_in = r_pt / 72.0
    rx = r_in * x_span / fig_w
    ry = r_in * y_span / fig_h
    return 1.15 * max(rx, ry)


def collision_pack(ids, pos, radii, bounds, n_iter=2500):
    """Deterministic pairwise disk packing for the intracellular cloud."""
    ids = list(ids)
    xmin, xmax, ymin, ymax = bounds
    for _ in range(n_iter):
        disp = {n: np.zeros(2, dtype=float) for n in ids}
        total = 0.0
        for i, a in enumerate(ids):
            xa, ya = pos[a]
            for b in ids[i + 1:]:
                xb, yb = pos[b]
                dx, dy = xb - xa, yb - ya
                dist = math.hypot(dx, dy)
                target = 1.18 * (radii[a] + radii[b])
                if dist < target:
                    overlap = target - dist
                    total += overlap
                    if dist < 1e-12:
                        # deterministic fallback direction
                        dx, dy, dist = 0.613, 0.790, 1.0
                    ux, uy = dx / dist, dy / dist
                    push = 0.55 * overlap
                    disp[a] -= np.array([ux, uy]) * push
                    disp[b] += np.array([ux, uy]) * push
        for n in ids:
            pos[n] += disp[n]
            r = radii[n]
            pos[n][0] = np.clip(pos[n][0], xmin + r, xmax - r)
            pos[n][1] = np.clip(pos[n][1], ymin + r, ymax - r)
        if total < 1e-5:
            break


def make_layout(nodes, edges, seed=31):
    graph = nx.Graph()
    graph.add_nodes_from(nodes['node_id'])
    agg = edges.groupby(['source_node', 'target_node']).size().reset_index(name='weight')
    for r in agg.itertuples(index=False):
        graph.add_edge(r.source_node, r.target_node, weight=r.weight)

    spring = nx.spring_layout(
        graph,
        seed=seed,
        weight='weight',
        k=1.72 / math.sqrt(max(graph.number_of_nodes(), 1)),
        iterations=1800,
        scale=1.0,
    )

    source_ids = nodes.loc[nodes.role == 'source', 'node_id'].tolist()
    carrier_ids = nodes.loc[nodes.role == 'carrier', 'node_id'].tolist()
    sink_ids = nodes.loc[nodes.role == 'terminal_sink', 'node_id'].tolist()

    pos = {}
    # Boundary nodes are now much closer to the intracellular cloud.
    pos.update(spread_two_columns(source_ids, spring, -0.675, -0.740))
    pos.update(spread_two_columns(sink_ids, spring, 0.675, 0.740))

    # Preserve an organic force-directed structure for intracellular nodes.
    cx = robust_rescale([spring[n][0] for n in carrier_ids], 0.585)
    cy = robust_rescale([spring[n][1] for n in carrier_ids], 0.925)
    for n, x, y in zip(carrier_ids, cx, cy):
        pos[n] = np.array([x, y], dtype=float)

    radii = {
        r.node_id: marker_radius_data(node_size(r.degree))
        for r in nodes.itertuples(index=False)
    }
    collision_pack(
        carrier_ids,
        pos,
        radii,
        bounds=(-0.625, 0.625, -0.955, 0.955),
        n_iter=3000,
    )
    return pos, radii


def draw_all_edges(ax, edges, pos):
    pair_counts = edges.groupby(['source_node', 'target_node']).size().to_dict()
    pair_seen = defaultdict(int)
    for e in edges.itertuples(index=False):
        u, v = e.source_node, e.target_node
        count = pair_counts[(u, v)]
        idx = pair_seen[(u, v)]
        pair_seen[(u, v)] += 1
        parallel = 0.0 if count == 1 else np.linspace(-0.11, 0.11, count)[idx]
        dy = pos[v][1] - pos[u][1]
        base = 0.016 * (1 if dy >= 0 else -1)
        ax.add_patch(FancyArrowPatch(
            posA=pos[u], posB=pos[v],
            arrowstyle='-|>', mutation_scale=5.8,
            linewidth=0.55, color='black', alpha=0.28,
            shrinkA=7, shrinkB=7,
            connectionstyle=f'arc3,rad={base + parallel}',
            zorder=1,
        ))


def rectangles_overlap(a, b, pad=0.002):
    ax0, ax1, ay0, ay1 = a
    bx0, bx1, by0, by1 = b
    return not (ax1 + pad < bx0 or bx1 + pad < ax0 or ay1 + pad < by0 or by1 + pad < ay0)


def label_bbox(x, y, text, fontsize, ha, va):
    """Approximate text rectangle in data coordinates for greedy placement."""
    w = 0.0060 * len(text) * (fontsize / 6.0)
    h = 0.0210 * (fontsize / 6.0)
    if ha == 'left':
        x0, x1 = x, x + w
    elif ha == 'right':
        x0, x1 = x - w, x
    else:
        x0, x1 = x - w / 2, x + w / 2
    if va == 'bottom':
        y0, y1 = y, y + h
    elif va == 'top':
        y0, y1 = y - h, y
    else:
        y0, y1 = y - h / 2, y + h / 2
    return (x0, x1, y0, y1)


def carrier_label_positions(nodes, pos, radii):
    """Greedy 8-direction label placement avoiding nodes and earlier labels."""
    carrier = nodes[nodes.role == 'carrier'].copy()
    carrier = carrier.sort_values(['degree', 'label'], ascending=[False, True])

    # Node rectangles used as obstacles.
    node_boxes = []
    for r in nodes.itertuples(index=False):
        x, y = pos[r.node_id]
        rr = radii[r.node_id] * 1.12
        node_boxes.append((r.node_id, (x-rr, x+rr, y-rr, y+rr)))

    placed = []
    result = {}
    dirs = [
        (1,0,'left','center'), (-1,0,'right','center'),
        (0,1,'center','bottom'), (0,-1,'center','top'),
        (0.707,0.707,'left','bottom'), (-0.707,0.707,'right','bottom'),
        (0.707,-0.707,'left','top'), (-0.707,-0.707,'right','top'),
    ]

    for r in carrier.itertuples(index=False):
        if r.degree >= 50:
            fs, weight = 8.8, 'bold'
        elif r.degree >= 15:
            fs, weight = 7.3, 'bold'
        elif r.degree >= 6:
            fs, weight = 5.9, 'normal'
        else:
            fs, weight = 5.0, 'normal'

        x, y = pos[r.node_id]
        offset = radii[r.node_id] + 0.008
        radial = np.array([x, y], dtype=float)
        if np.linalg.norm(radial) > 1e-9:
            radial /= np.linalg.norm(radial)

        best = None
        best_score = float('inf')
        for dx, dy, ha, va in dirs:
            tx, ty = x + dx*offset, y + dy*offset
            box = label_bbox(tx, ty, r.label, fs, ha, va)
            score = 0.0

            # Hard penalty for overlapping any node other than itself.
            for nid, nb in node_boxes:
                if nid != r.node_id and rectangles_overlap(box, nb, pad=0.001):
                    score += 60.0
            # Strong penalty for overlapping labels already placed.
            for pb in placed:
                if rectangles_overlap(box, pb, pad=0.002):
                    score += 100.0

            # Prefer labels extending away from network centre.
            score -= 4.0 * (dx*radial[0] + dy*radial[1])
            # Mild penalty for leaving central permitted region.
            if box[0] < -0.70 or box[1] > 0.70 or box[2] < -0.985 or box[3] > 0.985:
                score += 15.0

            if score < best_score:
                best_score = score
                best = (tx, ty, ha, va, fs, weight, box)

        tx, ty, ha, va, fs, weight, box = best
        placed.append(box)
        result[r.node_id] = (tx, ty, ha, va, fs, weight)
    return result


def make_figure(table_path, model_path, output_prefix, seed=31):
    table = pd.read_csv(table_path, sep='\t')
    raw_edges, original_names = reconstruct_reaction_edges(table)
    exchangeable = find_exchangeable_chemical_ids(model_path)
    nodes, edges = classify_and_remap(raw_edges, original_names, exchangeable)
    pos, radii = make_layout(nodes, edges, seed=seed)

    role_color = {
        'source': '#4E79A7',
        'carrier': '#59A14F',
        'terminal_sink': '#E15759',
    }

    fig, ax = plt.subplots(figsize=(20, 15), dpi=240)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    draw_all_edges(ax, edges, pos)

    for role in ('source', 'carrier', 'terminal_sink'):
        sub = nodes[nodes.role == role]
        ids = sub.node_id.tolist()
        ax.scatter(
            [pos[n][0] for n in ids],
            [pos[n][1] for n in ids],
            s=[node_size(d) for d in sub.degree],
            c=role_color[role],
            edgecolors='black', linewidths=0.75, zorder=3,
        )

    # Boundary labels: unboxed and placed outside the columns.
    for r in nodes[nodes.role == 'source'].itertuples(index=False):
        x, y = pos[r.node_id]
        fs = 5.1 if r.degree < 6 else 5.8
        ax.text(x-0.014, y, r.label, ha='right', va='center', fontsize=fs, zorder=4)
    for r in nodes[nodes.role == 'terminal_sink'].itertuples(index=False):
        x, y = pos[r.node_id]
        fs = 5.1 if r.degree < 6 else 5.8
        ax.text(x+0.014, y, r.label, ha='left', va='center', fontsize=fs, zorder=4)

    # Central labels: greedy obstacle-aware placement, no white boxes.
    cpos = carrier_label_positions(nodes, pos, radii)
    for r in nodes[nodes.role == 'carrier'].itertuples(index=False):
        tx, ty, ha, va, fs, weight = cpos[r.node_id]
        ax.text(tx, ty, r.label, ha=ha, va=va, fontsize=fs,
                fontweight=weight, zorder=4)

    ax.text(-0.705, 1.040, 'Sources', ha='center', va='bottom',
            fontsize=17, fontweight='bold')
    ax.text(0.0, 1.040, 'Intracellular carriers / intermediates',
            ha='center', va='bottom', fontsize=17, fontweight='bold')
    ax.text(0.705, 1.040, 'Terminal sinks', ha='center', va='bottom',
            fontsize=17, fontweight='bold')

    n_source = int((nodes.role == 'source').sum())
    n_carrier = int((nodes.role == 'carrier').sum())
    n_sink = int((nodes.role == 'terminal_sink').sum())
    n_dual = (nodes[nodes.boundary_accessible]
              .groupby('chemical_id')['role'].nunique().gt(1).sum())

    ax.text(-1.17, 1.155, 'E. coli electron-transfer network',
            ha='left', va='bottom', fontsize=22, fontweight='bold')
    ax.text(
        -1.17, 1.112,
        f'Complete resolved ETN: {len(edges)} directed reaction edges; '
        f'{n_source} source nodes, {n_carrier} intracellular nodes and '
        f'{n_sink} sink nodes. {n_dual} boundary metabolites are represented '
        f'as separate source and sink copies.',
        ha='left', va='bottom', fontsize=10.3,
    )

    legend = [
        Line2D([0],[0], marker='o', linestyle='None', label='Source',
               markerfacecolor=role_color['source'], markeredgecolor='black', markersize=9),
        Line2D([0],[0], marker='o', linestyle='None',
               label='Intracellular carrier / intermediate',
               markerfacecolor=role_color['carrier'], markeredgecolor='black', markersize=9),
        Line2D([0],[0], marker='o', linestyle='None', label='Terminal sink',
               markerfacecolor=role_color['terminal_sink'], markeredgecolor='black', markersize=9),
    ]
    ax.legend(handles=legend, loc='lower center', bbox_to_anchor=(0.5,-0.047),
              ncol=3, frameon=False, fontsize=10, columnspacing=1.7,
              handletextpad=0.5)

    ax.set_xlim(-1.19, 1.19)
    ax.set_ylim(-1.055, 1.205)
    ax.axis('off')

    output_prefix = Path(output_prefix)
    png = output_prefix.with_suffix('.png')
    pdf = output_prefix.with_suffix('.pdf')
    svg = output_prefix.with_suffix('.svg')
    node_csv = output_prefix.with_name(output_prefix.name + '_nodes.csv')
    edge_csv = output_prefix.with_name(output_prefix.name + '_edges.csv')

    fig.savefig(png, bbox_inches='tight', facecolor='white', dpi=300)
    fig.savefig(pdf, bbox_inches='tight', facecolor='white')
    fig.savefig(svg, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    nodes.to_csv(node_csv, index=False)
    edges.to_csv(edge_csv, index=False)

    print(f'Nodes: {len(nodes)} | Edges: {len(edges)}')
    print(nodes.role.value_counts().to_string())
    for p in (png, pdf, svg, node_csv, edge_csv):
        print('Saved:', p)
    return {'png':png, 'pdf':pdf, 'svg':svg, 'nodes':node_csv, 'edges':edge_csv}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('table', type=Path)
    parser.add_argument('model', type=Path)
    parser.add_argument('--output-prefix', type=Path,
                        default=Path('SuppFigS2_iML1515_ETN'))
    parser.add_argument('--seed', type=int, default=31)
    args = parser.parse_args()
    make_figure(args.table, args.model, args.output_prefix, seed=args.seed)


if __name__ == '__main__':
    main()
