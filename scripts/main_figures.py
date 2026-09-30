#!/usr/bin/env python3
"""Generate main-text Figures 2–4.

This module defines the visual appearance of manuscript Figures 2, 3 and 4 at
183-mm double-column publication width, including typography, colors, closed axes, panel
spacing, limits, annotations and legends.

Both PNG (600 dpi) and vector PDF files are written to figures/publication/.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

MM = 1 / 25.4
FIG_W = 183 * MM

STATE_COLORS = {"high": "#1f77b4", "low": "#ff7f0e"}
ACCEPTOR_COLORS = {"oxygen": "#1f77b4", "nitrate": "#ff7f0e", "tmao": "#2ca02c"}
SUBSTRATE_COLORS = {
    "glucose": "#1f77b4",
    "fructose": "#ff7f0e",
    "galactose": "#2ca02c",
    "xylose": "#d62728",
    "gluconate": "#9467bd",
    "glycerol": "#8c564b",
}
ORGANISM_COLORS = {
    "E. coli": "#1f77b4",
    "B. subtilis": "#ff7f0e",
    "S. enterica": "#2ca02c",
    "K. phaffii": "#d62728",
    "S. cerevisiae": "#9467bd",
}
ORGANISM_MARKERS = {
    "E. coli": "o",
    "B. subtilis": "s",
    "S. enterica": "^",
    "K. phaffii": "D",
    "S. cerevisiae": "P",
}


def _set_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 6.5,
        "axes.labelsize": 6.6,
        "xtick.labelsize": 5.8,
        "ytick.labelsize": 5.8,
        "legend.fontsize": 5.3,
        "axes.linewidth": 0.55,
        "lines.linewidth": 0.9,
    })


def _panel_letter(ax, letter: str, x: float = -0.08, y: float = 1.03) -> None:
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=7.4, fontweight="bold",
            va="bottom", ha="left")


def _close_axes(ax) -> None:
    for side in ["top", "right", "bottom", "left"]:
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(0.55)
    ax.tick_params(direction="out", width=0.55)


def _save(fig, out: Path, stem: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{stem}.png", dpi=600, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(out / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def generate_figure2(root: Path, out: Path) -> None:
    _set_style()
    condition_summary = pd.read_csv(root / "results/fba_condition_summary.csv").set_index("condition")
    flow_summary = pd.read_csv(root / "results/electron_path_summary_by_condition.csv").set_index("condition")
    relay = pd.read_csv(root / "results/electron_carrier_relay_by_condition.csv")
    edges = pd.read_csv(root / "results/redistribution_aerobic_vs_anaerobic.csv")
    terminal = pd.read_csv(root / "results/electron_terminal_fates_interpreted_by_condition.csv")
    mfa = pd.read_csv(root / "results/publication/gonzalez_13c_mfa_electron_flux_metrics.csv")

    # Publication panels a-e.
    fig = plt.figure(figsize=(FIG_W, 116 * MM))
    outer = fig.add_gridspec(2, 1, height_ratios=[0.98, 1.0], hspace=0.72)
    top = outer[0].subgridspec(1, 2, width_ratios=[1.08, 0.92], wspace=0.90)
    lower = outer[1].subgridspec(1, 3, width_ratios=[1.05, 0.72, 1.02], wspace=0.72)
    fig.subplots_adjust(left=0.068, right=0.993, top=0.968, bottom=0.165)

    # 2a
    ax = fig.add_subplot(top[0])
    q_a = float(condition_summary.loc["aerobic", "qGLC"])
    q_n = float(condition_summary.loc["anaerobic", "qGLC"])
    net_a = float(flow_summary.loc["aerobic", "net_source_flux"])
    net_n = float(flow_summary.loc["anaerobic", "net_source_flux"])
    aerobic_vals = np.array([q_a, net_a, net_a / q_a])
    anaerobic_vals = np.array([q_n, net_n, net_n / q_n])
    relative = 100 * anaerobic_vals / aerobic_vals
    labels = ["Glucose uptake", "Absolute electron flux", "Electron flux / glucose"]
    ypos = np.arange(3)
    formats = ["{:.1f}", "{:.1f}", "{:.2f}"]
    for y, r in zip(ypos, relative):
        ax.plot([100, r], [y, y], lw=0.9, color="black", zorder=1)
    ax.scatter(np.repeat(100, 3), ypos, s=14, color=STATE_COLORS["high"], label="Aerobic", zorder=3)
    ax.scatter(relative, ypos, s=14, color=STATE_COLORS["low"], label="Anaerobic", zorder=3)
    ax.axvline(100, ls="--", lw=0.6, color="black")
    ax.set_xlim(40, 200)
    ax.set_xticks([40, 80, 120, 160, 200])
    ax.set_yticks(ypos, labels)
    ax.invert_yaxis()
    ax.set_ylim(2.40, -0.55)
    ax.set_xlabel("Anaerobic relative to aerobic (%)")
    _panel_letter(ax, "a")
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.55, 1.03), ncol=2,
              handletextpad=0.35, columnspacing=0.8, borderaxespad=0.0)
    for i, (y, a_val, n_val, r, fmt) in enumerate(zip(ypos, aerobic_vals, anaerobic_vals, relative, formats)):
        if i == 1:
            ax.annotate(fmt.format(n_val), xy=(r, y), xytext=(-5, 6), textcoords="offset points",
                        ha="right", va="bottom", fontsize=5.8)
            ax.annotate(fmt.format(a_val), xy=(100, y), xytext=(5, 6), textcoords="offset points",
                        ha="left", va="bottom", fontsize=5.8)
        elif i == 2:
            ax.annotate(fmt.format(a_val), xy=(100, y), xytext=(0, 8), textcoords="offset points",
                        ha="center", va="bottom", fontsize=5.8)
            ax.annotate(fmt.format(n_val), xy=(r, y), xytext=(0, 10), textcoords="offset points",
                        ha="center", va="bottom", fontsize=5.8)
        else:
            ax.annotate(fmt.format(a_val), xy=(100, y), xytext=(0, 6), textcoords="offset points",
                        ha="center", va="bottom", fontsize=5.8)
            ax.annotate(fmt.format(n_val), xy=(r, y), xytext=(0, 6), textcoords="offset points",
                        ha="center", va="bottom", fontsize=5.8)
    _close_axes(ax)

    # 2b
    ax = fig.add_subplot(top[1])
    picks = [("nadh_c", "NADH"), ("nadph_c", "NADPH"), ("fadh2_c", "FADH$_2$"), ("q8h2_c", "Ubiquinol-8")]
    names, aa, bb = [], [], []
    for cid, name in picks:
        gg = relay[relay["carrier"] == cid]
        if len(gg):
            names.append(name)
            aa.append(float(gg[gg["condition"] == "aerobic"]["relay_throughflow"].iloc[0]))
            bb.append(float(gg[gg["condition"] == "anaerobic"]["relay_throughflow"].iloc[0]))
    x = np.arange(len(names)); w = 0.36
    ax.bar(x - w/2, aa, width=w, color=STATE_COLORS["high"], label="Aerobic")
    ax.bar(x + w/2, bb, width=w, color=STATE_COLORS["low"], label="Anaerobic")
    ax.set_xticks(x, names, rotation=28, ha="right")
    ax.set_ylabel("Carrier relay flux\n(e$^{-}$ mmol gDW$^{-1}$ h$^{-1}$)")
    _panel_letter(ax, "b")
    ax.legend(frameon=False, loc="upper right", handletextpad=0.35, borderaxespad=0.15)
    _close_axes(ax)

    # 2c
    ax = fig.add_subplot(lower[0])
    resolved = edges[
        ~edges["donor"].str.contains("unresolved", case=False, na=False)
        & ~edges["acceptor"].str.contains("unresolved", case=False, na=False)
    ].copy()
    resolved["abs_delta"] = resolved["delta"].abs()
    sel = resolved.sort_values("abs_delta", ascending=False).head(8).copy()
    pretty = {
        ("q8h2_c", "h2o_c"): "Ubiquinol-8 → H$_2$O",
        ("nadh_c", "q8h2_c"): "NADH → Ubiquinol-8",
        ("nadh_c", "etoh_c"): "NADH → Ethanol",
        ("g6p_c", "nadph_c"): "G6P → NADPH",
        ("nadh_c", "nadph_c"): "NADH → NADPH",
        ("nadh_c", "fadh2_c"): "NADH → FADH$_2$",
        ("3pg_c", "nadh_c"): "3PG → NADH",
        ("mal__L_c", "nadh_c"): "Malate → NADH",
    }
    sel["label"] = [pretty.get((d, a), f"{d} → {a}") for d, a in zip(sel["donor"], sel["acceptor"])]
    sel = sel.sort_values("delta")
    yy = np.arange(len(sel))
    ax.barh(yy, sel["delta"].values, color="#7f7f7f")
    ax.axvline(0, lw=0.55, color="black")
    ax.set_yticks(yy, sel["label"].values)
    ax.set_xlabel("Anaerobic − aerobic edge flux\n(e$^{-}$ mmol gDW$^{-1}$ h$^{-1}$)")
    _panel_letter(ax, "c")
    _close_axes(ax)

    # 2d
    ax = fig.add_subplot(lower[1])
    fate_order = ["respiration", "fermentation_associated", "biosynthesis_or_other", "unresolved_other"]
    fate_labels = ["Respiration", "Fermentation-associated", "Biosynthesis / other", "Other unresolved"]
    fate_colors = ["#4c78a8", "#f58518", "#54a24b", "#b279a2"]
    bottom = np.zeros(2)
    for fate, label, color in zip(fate_order, fate_labels, fate_colors):
        vals = []
        for cond in ["aerobic", "anaerobic"]:
            hit = terminal[(terminal["condition"] == cond) & (terminal["interpreted_fate"] == fate)]
            vals.append(100 * float(hit["fraction"].iloc[0]) if len(hit) else 0)
        ax.bar([0, 1], vals, bottom=bottom, label=label, width=0.62, color=color)
        bottom += np.asarray(vals)
    ax.set_xticks([0, 1], ["Aerobic", "Anaerobic"])
    ax.set_ylim(0, 100)
    ax.set_ylabel("Terminal electron delivery (%)")
    _panel_letter(ax, "d")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, -0.17), ncol=1,
              fontsize=4.7, handletextpad=0.35, labelspacing=0.28, borderaxespad=0.0)
    _close_axes(ax)

    # 2e
    ax = fig.add_subplot(lower[2])
    g = mfa[mfa["substrate"] == "glucose"].set_index("state")
    a = g.loc["aerobic"]; n = g.loc["anaerobic"]
    for row, label, color in [(a, "Aerobic", STATE_COLORS["high"]), (n, "Anaerobic", STATE_COLORS["low"])]:
        x0 = float(row.q_substrate)
        y0 = float(row.absolute_net_e_flux)
        x_sd = float(row.q_sd)
        y_sd = float(row.net_e_per_substrate) * x_sd
        ax.errorbar(x0, y0, xerr=x_sd, yerr=y_sd, fmt="o", ms=3.6, color=color,
                    ecolor="black", capsize=1.8, lw=0.8)
        ax.annotate(label, (x0, y0), xytext=((4, 4) if label == "Aerobic" else (4, -9)),
                    textcoords="offset points", fontsize=5.5, color=color)
    ax.set_xlabel("Measured glucose uptake\n(mmol gDW$^{-1}$ h$^{-1}$)")
    ax.set_ylabel("Net electron flux\n(e$^{-}$ mmol gDW$^{-1}$ h$^{-1}$)")
    ax.set_xlim(7.3, 14.7)
    ax.set_xticks([8, 10, 12, 14])
    ax.set_ylim(40, 90)
    _panel_letter(ax, "e")
    _close_axes(ax)

    _save(fig, out, "Fig2_glucose_decoupling_and_disposal_capacity")


def generate_figure3(root: Path, out: Path) -> None:
    _set_style()
    sim = pd.read_csv(root / "results/publication/acceptor_electron_capacity_titration.csv")
    subs = pd.read_csv(root / "results/publication/substrate_flux_compression.csv").set_index("substrate")
    mfa = pd.read_csv(root / "results/publication/gonzalez_13c_mfa_electron_flux_metrics.csv")
    core = pd.read_csv(root / "results/publication/substrate_core_redox_contribution.csv")
    factorial = pd.read_csv(root / "results/publication/factorial_carbon_oxygen_states.csv")
    scaled = pd.read_csv(root / "results/publication/factorial_scaled_expansion_curves.csv")

    fig = plt.figure(figsize=(FIG_W, 104 * MM))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.0, 1.28], height_ratios=[1.0, 1.0])
    fig.subplots_adjust(left=0.065, right=0.993, top=0.97, bottom=0.11, wspace=0.56, hspace=0.62)
    order = ["glucose", "fructose", "galactose", "xylose", "gluconate", "glycerol"]
    label_map = {"oxygen": "O$_2$", "nitrate": "Nitrate", "tmao": "TMAO"}

    # 3a
    ax = fig.add_subplot(gs[0, 0])
    for acceptor in ["oxygen", "nitrate", "tmao"]:
        g = sim[sim["acceptor"] == acceptor].sort_values("nominal_e_capacity")
        ax.plot(g["nominal_e_capacity"], g["net_per_glucose"], marker="o", lw=0.95, ms=3.2,
                label=label_map[acceptor], color=ACCEPTOR_COLORS[acceptor])
    ax.set_xlim(-2, 82)
    ax.set_xticks(np.arange(0, 81, 10))
    ax.set_ylim(0, 16)
    ax.set_xlabel("Nominal electron-accepting capacity\n(e$^{-}$ mmol gDW$^{-1}$ h$^{-1}$)")
    ax.set_ylabel("Net electron flux\n(e$^{-}$ per glucose)")
    _panel_letter(ax, "a")
    ax.legend(frameon=False, loc="upper left", handletextpad=0.35, borderaxespad=0.15)
    _close_axes(ax)

    # 3b
    ax = fig.add_subplot(gs[0, 1])
    x_all = sim["co2"].to_numpy(float)
    y_all = sim["net_per_glucose"].to_numpy(float)
    r = np.corrcoef(x_all, y_all)[0, 1]
    coef = np.polyfit(x_all, y_all, 1)
    xx = np.linspace(x_all.min(), x_all.max(), 200)
    for acceptor in ["oxygen", "nitrate", "tmao"]:
        g = sim[sim["acceptor"] == acceptor]
        ax.scatter(g["co2"], g["net_per_glucose"], s=13, label=label_map[acceptor], color=ACCEPTOR_COLORS[acceptor])
    ax.plot(xx, np.polyval(coef, xx), lw=0.75, color="black")
    ax.set_ylim(0, 16)
    ax.set_xlabel("Net CO$_2$ production\n(mmol gDW$^{-1}$ h$^{-1}$)")
    ax.set_ylabel("Net electron flux\n(e$^{-}$ per glucose)")
    _panel_letter(ax, "b")
    ax.text(0.04, 0.76, f"Pearson r = {r:.3f}", transform=ax.transAxes, va="top", fontsize=5.5)
    ax.legend(frameon=False, loc="upper left", ncol=3, handletextpad=0.35, columnspacing=0.8,
              borderaxespad=0.15)
    _close_axes(ax)

    # 3c
    ax = fig.add_subplot(gs[0, 2])
    d = subs.loc[order]
    x = np.arange(len(order)); w = 0.36
    ax.bar(x - w/2, d["net_per_glucose_anoxic"], w, label="Zero respiratory capacity", color=STATE_COLORS["low"])
    ax.bar(x + w/2, d["net_per_glucose_oxygen"], w, label="High respiratory capacity", color=STATE_COLORS["high"])
    ax.set_xticks(x, [s.capitalize() for s in order], rotation=20, ha="right")
    ax.set_ylabel("Net electron flux\n(e$^{-}$ per substrate)")
    ax.set_ylim(0, 15)
    _panel_letter(ax, "c")
    ax.legend(frameon=False, loc="upper left", ncol=2, handlelength=1.4, columnspacing=0.8,
              handletextpad=0.35, borderaxespad=0.15)
    _close_axes(ax)

    # 3d
    ax = fig.add_subplot(gs[1, 0])
    states = ["aerobic", "anaerobic"]
    x = np.arange(2); w = 0.36
    for j, substrate in enumerate(["glucose", "xylose"]):
        g = mfa[mfa["substrate"] == substrate].set_index("state").loc[states]
        y = g["net_e_per_substrate"].to_numpy(float)
        lo = y - g["net_e_lb_envelope"].to_numpy(float)
        hi = g["net_e_ub_envelope"].to_numpy(float) - y
        xpos = x + (j - 0.5) * w
        ax.bar(xpos, y, width=w, label=substrate.capitalize(), color=SUBSTRATE_COLORS[substrate])
        ax.errorbar(xpos, y, yerr=np.vstack([lo, hi]), fmt="none", ecolor="black",
                    elinewidth=0.7, capsize=1.8, capthick=0.7, zorder=3)
    ax.set_xticks(x, ["Aerobic", "Anaerobic"])
    ax.set_ylabel(r"$^{13}$C-MFA net electron flux" + "\n" + r"(e$^{-}$ per substrate)")
    ax.set_ylim(bottom=0)
    _panel_letter(ax, "d")
    ax.legend(frameon=False, loc="upper right", handletextpad=0.35, borderaxespad=0.15)
    _close_axes(ax)

    # 3e
    ax = fig.add_subplot(gs[1, 1])
    core_by = core.set_index("substrate").loc[order]
    components = {"GAPD": [], "GND": [], "G3PD2": [], "Other": []}
    for substrate, row in core_by.iterrows():
        vals = {"GAPD": 0., "GND": 0., "G3PD2": 0.}
        scale = float(row.selected_e_per_substrate) / float(row.selected_e_flux) if float(row.selected_e_flux) else 0
        for part in str(row.details).split(" | "):
            rxn, rest = part.split(":", 1)
            if rxn in vals:
                vals[rxn] = float(rest.split("=")[-1]) * scale
        selected = sum(vals.values())
        for k in ["GAPD", "GND", "G3PD2"]:
            components[k].append(vals[k])
        components["Other"].append(max(0, float(row.net_e_per_substrate) - selected))
    x = np.arange(len(order)); bottom = np.zeros(len(order))
    comp_colors = {"GAPD": "#4c78a8", "GND": "#f58518", "G3PD2": "#54a24b", "Other": "#bab0ac"}
    for k in ["GAPD", "GND", "G3PD2", "Other"]:
        v = np.asarray(components[k])
        ax.bar(x, v, bottom=bottom, label=k, color=comp_colors[k])
        bottom += v
    ax.set_xticks(x, [s.capitalize() for s in order], rotation=24, ha="right")
    ax.set_ylabel("Anoxic net electron flux\n(e$^{-}$ per substrate)")
    ax.set_ylim(bottom=0)
    _panel_letter(ax, "e")
    ax.legend(frameon=False, ncol=1, fontsize=5.0, handletextpad=0.35,
              loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    _close_axes(ax)

    # 3f
    ax = fig.add_subplot(gs[1, 2])
    for sub in order:
        g = factorial[factorial["substrate"] == sub].sort_values("oxygen_cap")
        ax.plot(g["oxygen_cap"], g["net_per_substrate"], marker="o", ms=3.0, lw=0.9,
                label=sub.capitalize(), color=SUBSTRATE_COLORS[sub])
    ax.set_xlabel("O$_2$-uptake capacity\n(mmol gDW$^{-1}$ h$^{-1}$)")
    ax.set_ylabel("Net electron flux\n(e$^{-}$ per substrate)")
    ax.set_ylim(bottom=0)
    ax.set_xlim(0, 20)
    ax.set_xticks([0, 5, 10, 15, 20])
    _panel_letter(ax, "f")
    ax.legend(frameon=False, ncol=2, loc="upper left", fontsize=5.0,
              handletextpad=0.35, columnspacing=0.75, borderaxespad=0.15)
    ax.text(0.08, 0.66, r"$J_e(s,d)\ \approx\ B_s + A_s\,f(d)$", transform=ax.transAxes, fontsize=5.8,
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="none", alpha=0.90))
    inset = ax.inset_axes([0.64, 0.10, 0.30, 0.30])
    for sub in order:
        g = scaled[scaled["substrate"] == sub].sort_values("oxygen_cap")
        inset.plot(g["oxygen_cap"], g["scaled_expansion"], marker="o", ms=1.9, lw=0.75, color=SUBSTRATE_COLORS[sub])
    inset.set_xlim(0, 20)
    inset.set_ylim(0, 1.05)
    inset.set_xticks([0, 10, 20])
    inset.set_yticks([0, 0.5, 1.0])
    inset.set_xlabel("O$_2$ cap.", fontsize=4.8, labelpad=1)
    inset.set_ylabel("Scaled\nexpansion", fontsize=4.8, labelpad=1)
    inset.tick_params(labelsize=4.4)
    _close_axes(ax)
    _close_axes(inset)

    _save(fig, out, "Fig3_carbon_source_disposal_capacity_ecoli")


def _connected_clusters(g: pd.DataFrame, xtol: float = 0.12, ytol: float = 0.12):
    idx = list(g.index)
    adj = {i: set() for i in idx}
    for a_idx in range(len(idx)):
        for b_idx in range(a_idx + 1, len(idx)):
            ia, ib = idx[a_idx], idx[b_idx]
            if (abs(g.loc[ia, "low_net_per_substrate"] - g.loc[ib, "low_net_per_substrate"]) <= xtol and
                abs(g.loc[ia, "high_net_per_substrate"] - g.loc[ib, "high_net_per_substrate"]) <= ytol):
                adj[ia].add(ib); adj[ib].add(ia)
    seen, clusters = set(), []
    for i in idx:
        if i in seen:
            continue
        stack, cluster = [i], []
        while stack:
            q = stack.pop()
            if q in seen:
                continue
            seen.add(q); cluster.append(q); stack.extend(adj[q] - seen)
        clusters.append(cluster)
    return clusters


def generate_figure4(root: Path, out: Path) -> None:
    _set_style()
    comp = pd.read_csv(root / "results/cross_species/electron_flux_compression_summary.csv")
    jouhten = pd.read_csv(root / "results/publication/jouhten2008_matched_model_experiment.csv")
    pairs = pd.read_csv(root / "results/publication/cross_species_carbon_acceptor_pairs.csv")
    fromanger = pd.read_csv(root / "results/publication/fromanger2010_electron_balance.csv")
    steinsiek_points = pd.read_csv(root / "results/publication/steinsiek2014_predictive_transfer_points.csv")
    steinsiek_stats = __import__("json").loads((root / "results/publication/steinsiek2014_predictive_transfer_statistics.json").read_text())

    fig = plt.figure(figsize=(FIG_W, 118 * MM))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.96])
    fig.subplots_adjust(left=0.065, right=0.992, top=0.97, bottom=0.115, wspace=0.58, hspace=0.70)

    org_order = ["E. coli", "S. cerevisiae", "B. subtilis", "S. enterica", "K. phaffii"]
    comp_plot = comp.set_index("organism").loc[org_order].reset_index()

    # 4a
    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(len(comp_plot)); w = 0.36
    ax.bar(x - w/2, comp_plot["high_respiratory_capacity_net_e_per_glucose"], width=w, label="High disposal capacity", color=STATE_COLORS["high"])
    ax.bar(x + w/2, comp_plot["low_oxygen_net_e_per_glucose"], width=w, label="Low disposal capacity", color=STATE_COLORS["low"])
    ax.set_xticks(x, comp_plot["organism"], rotation=25, ha="right")
    ax.set_ylabel("Net electron flux\n(e$^{-}$ per glucose)")
    ax.set_ylim(0, 14)
    _panel_letter(ax, "a")
    ax.legend(frameon=False, loc="upper left", ncol=2, handlelength=1.4, columnspacing=0.8)
    _close_axes(ax)

    # 4b
    ax = fig.add_subplot(gs[0, 1])
    ax.axhline(0, lw=0.55, color="black")
    ax.axvline(0, lw=0.55, color="black")
    label_offsets = {
        "E. coli": (3, 3),
        "S. cerevisiae": (-4, 4),
        "B. subtilis": (-6, -1),
        "S. enterica": (3, 4),
        "K. phaffii": (3, 5),
    }
    alignments = {
        "E. coli": ("left", "bottom"),
        "S. cerevisiae": ("right", "bottom"),
        "B. subtilis": ("right", "center"),
        "S. enterica": ("left", "bottom"),
        "K. phaffii": ("left", "bottom"),
    }
    for _, rr in comp_plot.iterrows():
        org = rr.organism
        x0 = rr["transfer_depth_change_pct"]
        y0 = rr["nadh_relay_change_pct"]
        ax.scatter(x0, y0, s=18, color=ORGANISM_COLORS[org], marker=ORGANISM_MARKERS[org], zorder=3)
        dx, dy = label_offsets[org]
        ha, va = alignments[org]
        ax.annotate(org, (x0, y0), xytext=(dx, dy), textcoords="offset points",
                    fontsize=5.2, color=ORGANISM_COLORS[org], ha=ha, va=va)
    xvals = comp_plot["transfer_depth_change_pct"].to_numpy(float)
    yvals = comp_plot["nadh_relay_change_pct"].to_numpy(float)
    ax.set_xlim(min(-22, xvals.min()-1.5), max(6.5, xvals.max()+2.8))
    ax.set_ylim(min(-63, yvals.min()-4.0), max(4.5, yvals.max()+4.0))
    ax.set_xlabel("Change in effective transfer depth (%)")
    ax.set_ylabel("Change in NADH relay flux (%)")
    _panel_letter(ax, "b")
    _close_axes(ax)

    # 4c
    ax = fig.add_subplot(gs[0, 2])
    j = jouhten.sort_values("published_mean_OUR_mmol_gCDW_h")
    ax.errorbar(j["published_mean_OUR_mmol_gCDW_h"], j["experimental_source_generation_e_per_glucose"],
                yerr=j["experimental_sd_e_per_glucose"], marker="o", ms=3.2, lw=0.9, capsize=1.8,
                label="$^{13}$C-MFA", color="#1f77b4")
    ax.plot(j["published_mean_OUR_mmol_gCDW_h"], j["model_source_generation_e_per_glucose"],
            marker="s", ms=3.0, lw=0.9, label="Genome-scale model", color="#7f7f7f")
    ax.set_xlabel("Mean O$_2$ uptake\n(mmol gCDW$^{-1}$ h$^{-1}$)")
    ax.set_ylabel("Electron generation\n(e$^{-}$ per glucose)")
    _panel_letter(ax, "c")
    ax.legend(frameon=False, loc="lower right", handlelength=1.4)
    ax.text(0.04, 0.95, "r = 0.987\nRMSE = 0.86", transform=ax.transAxes, va="top", fontsize=5.2)
    _close_axes(ax)

    # 4d
    ax = fig.add_subplot(gs[1, 0])
    pairs_plot = pairs.copy()
    jitter_templates = {
        2: [(-0.07, -0.07), (0.07, 0.07)],
        3: [(-0.09, -0.06), (0.09, -0.06), (0.0, 0.09)],
        4: [(-0.09, -0.08), (0.09, -0.08), (-0.09, 0.08), (0.09, 0.08)],
    }
    pairs_plot["plot_x"] = pairs_plot["low_net_per_substrate"]
    pairs_plot["plot_y"] = pairs_plot["high_net_per_substrate"]
    d_org_order = ["E. coli", "B. subtilis", "S. enterica", "K. phaffii", "S. cerevisiae"]
    for org in d_org_order:
        g = pairs_plot[pairs_plot["organism"] == org]
        for cluster in _connected_clusters(g):
            if len(cluster) > 1:
                cluster = sorted(cluster, key=lambda i: pairs_plot.loc[i, "substrate"])
                for i, (dx, dy) in zip(cluster, jitter_templates[len(cluster)]):
                    pairs_plot.loc[i, "plot_x"] += dx
                    pairs_plot.loc[i, "plot_y"] += dy
    for org in d_org_order:
        g = pairs_plot[pairs_plot["organism"] == org]
        ax.scatter(g["plot_x"], g["plot_y"], s=18, marker=ORGANISM_MARKERS[org], label=org,
                   zorder=3, color=ORGANISM_COLORS[org])
    lim = [1.7, 13.5]
    ax.plot(lim, lim, ls="--", lw=0.65, color="black")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("Electron flux at low disposal capacity\n(e$^{-}$ per substrate)")
    ax.set_ylabel("Electron flux at high disposal capacity\n(e$^{-}$ per substrate)")
    _panel_letter(ax, "d")
    ax.text(0.98, 0.34, "22 feasible pairs / 30 attempted", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=5.0)
    ax.legend(frameon=False, fontsize=5.0, ncol=2, loc="lower right",
              handletextpad=0.4, columnspacing=0.8)
    _close_axes(ax)

    # 4e: transfer of the parental acetate-response shape across
    # respiratory architectures.
    ax = fig.add_subplot(gs[1, 1])
    colors = {
        "TBE031_bo_only": "#ff7f0e",
        "TBE032_bdII_only": "#2ca02c",
        "TBE042_bdI_only": "#9467bd",
    }
    labels = {
        "TBE031_bo_only": "bo-only",
        "TBE032_bdII_only": "bd-II-only",
        "TBE042_bdI_only": "bd-I-only",
    }
    for target, color in colors.items():
        g = steinsiek_points[steinsiek_points.target_architecture == target].sort_values("aerobiosis_percent")
        ax.plot(g.aerobiosis_percent, g.predicted_acetate_mmol_gDW_h, color=color, lw=0.9)
        anchors = g[g.role == "scaling_endpoint"]
        held = g[g.role == "held_out_prediction"]
        ax.scatter(anchors.aerobiosis_percent, anchors.observed_acetate_mmol_gDW_h,
                   marker="s", s=18, color=color, zorder=3)
        ax.scatter(held.aerobiosis_percent, held.observed_acetate_mmol_gDW_h,
                   marker="o", s=18, facecolor="white", edgecolor=color,
                   linewidth=0.9, zorder=3)
        ax.text(112, float(g.loc[g.aerobiosis_percent == 100, "predicted_acetate_mmol_gDW_h"].iloc[0]),
                labels[target], color=color, fontsize=5.0, va="center")
    ax.set_xlabel("Aerobiosis (%)")
    ax.set_ylabel("Acetate formation\n(mmol gDW$^{-1}$ h$^{-1}$)")
    ax.set_xlim(-5, 155); ax.set_ylim(0, 9)
    _panel_letter(ax, "e")
    ax.text(0.04, 0.08, f"$R^2_{{pred}}$ = {steinsiek_stats['predictive_r2']:.3f}\nRMSE = {steinsiek_stats['rmse_mmol_gDW_h']:.2f}",
            transform=ax.transAxes, fontsize=5.2)
    _close_axes(ax)

    # 4f
    ax = fig.add_subplot(gs[1, 2])
    from sklearn.isotonic import IsotonicRegression

    weusthuis = pd.read_csv(root / "data/experimental/weusthuis1994/s_cerevisiae_oxygen_series.csv")
    # Primary treatment of reported non-detections (<1 mM): zero culture ethanol.
    censored_glucose = (
        (weusthuis["substrate"] == "Glucose")
        & weusthuis["ethanol_out_mM"].isna()
        & weusthuis["ethanol_out_censored_lt_mM"].notna()
    )
    weusthuis.loc[censored_glucose, "ethanol_out_mM"] = 0.0
    ethanol_feed_mM = 12.0
    weusthuis["q_ethanol"] = (
        weusthuis["dilution_h"]
        * (weusthuis["ethanol_out_mM"] - ethanol_feed_mM)
        / weusthuis["dry_weight_g_L"]
    )
    glc = weusthuis[(weusthuis["substrate"] == "Glucose") & weusthuis["ethanol_out_mM"].notna()].sort_values("oxygen_in_mmol_L_h")
    mal = weusthuis[(weusthuis["substrate"] == "Maltose") & weusthuis["ethanol_out_mM"].notna()].sort_values("oxygen_in_mmol_L_h")

    # Infer a monotonic response to disposal capacity from glucose only.
    xg = glc["oxygen_in_mmol_L_h"].to_numpy(float)
    yg = glc["q_ethanol"].to_numpy(float)
    iso = IsotonicRegression(increasing=False, out_of_bounds="clip")
    iso.fit(xg, yg)
    glc_low = float(glc.loc[np.isclose(glc["oxygen_in_mmol_L_h"], 0.0), "q_ethanol"].iloc[0])
    glc_high = float(iso.predict([30.8])[0])

    def glucose_response_shape(x):
        yy = iso.predict(np.asarray(x, dtype=float))
        return (yy - glc_low) / (glc_high - glc_low)

    # Maltose contributes only the low- and high-disposal-capacity anchors.
    mal_zero = mal[np.isclose(mal["oxygen_in_mmol_L_h"], 0.0)]
    mal_low = float(mal_zero["q_ethanol"].mean())
    mal_low_sd = float(mal_zero["q_ethanol"].std(ddof=1))
    mal_high = float(mal.loc[np.isclose(mal["oxygen_in_mmol_L_h"], 30.8), "q_ethanol"].iloc[0])

    heldout_x = np.array([5.4, 10.8, 18.7, 22.3, 24.8])
    heldout = mal[mal["oxygen_in_mmol_L_h"].isin(heldout_x)].sort_values("oxygen_in_mmol_L_h").copy()
    heldout["q_pred"] = mal_low + (mal_high - mal_low) * glucose_response_shape(heldout["oxygen_in_mmol_L_h"].to_numpy(float))
    obs = heldout["q_ethanol"].to_numpy(float)
    pred = heldout["q_pred"].to_numpy(float)
    r = float(np.corrcoef(obs, pred)[0, 1])
    r2_pred = float(1.0 - np.sum((obs - pred) ** 2) / np.sum((obs - obs.mean()) ** 2))
    rmse = float(np.sqrt(np.mean((obs - pred) ** 2)))

    # Show the transferred prediction at equal 2-unit O2-feed increments, plus the 30.8 endpoint.
    x_curve = np.r_[np.arange(0.0, 31.0, 2.0), 30.8]
    y_curve = mal_low + (mal_high - mal_low) * glucose_response_shape(x_curve)
    ax.plot(x_curve, y_curve, lw=0.9, color="black", label="Prediction", zorder=2)

    exp_color = ORGANISM_COLORS["S. cerevisiae"]
    ax.errorbar([0.0], [mal_low], yerr=[mal_low_sd], fmt="s", ms=3.4, color=exp_color,
                ecolor="black", capsize=1.8, elinewidth=0.65, label="Anchors", zorder=3)
    ax.plot([30.8], [mal_high], linestyle="none", marker="s", ms=3.4, color=exp_color, zorder=3)
    ax.plot(heldout["oxygen_in_mmol_L_h"], heldout["q_ethanol"], linestyle="none", marker="o",
            ms=3.8, markerfacecolor="white", markeredgecolor=exp_color, markeredgewidth=0.9,
            label="Held-out states", zorder=3)

    ax.set_xlabel("O$_2$ feed\n(mmol l$^{-1}$ h$^{-1}$)")
    ax.set_ylabel("Net ethanol production\n(mmol gDW$^{-1}$ h$^{-1}$)")
    ax.set_xlim(-0.8, 32.4)
    ax.set_ylim(-0.7, 11.9)
    ax.set_xticks([0, 10, 20, 30])
    ax.set_yticks([0, 2, 4, 6, 8, 10])
    _panel_letter(ax, "f")
    ax.legend(frameon=False, loc="lower left", handlelength=1.4, handletextpad=0.4,
              borderaxespad=0.15)
    ax.text(0.98, 0.95, f"r = {r:.3f}\nRMSE = {rmse:.2f}\n$R^2_{{pred}}$ = {r2_pred:.3f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.2)
    _close_axes(ax)

    _save(fig, out, "Fig4_cross_species_organizing_principles")


def generate_all(root: Path | None = None, out: Path | None = None) -> None:
    root = Path(root) if root else Path(__file__).resolve().parents[1]
    out = Path(out) if out else root / "figures" / "publication"
    generate_figure2(root, out)
    generate_figure3(root, out)
    generate_figure4(root, out)


if __name__ == "__main__":
    generate_all()
