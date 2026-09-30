from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "publication" / "cross_species_carbon_acceptor_pairs.csv"
FEAS = ROOT / "results" / "publication" / "cross_species_carbon_acceptor_feasibility.csv"


def test_cross_species_principle_publication_result():
    df = pd.read_csv(RESULT)
    assert set(df["organism"]) == {"E. coli", "B. subtilis", "S. enterica", "K. phaffii", "S. cerevisiae"}
    assert len(df) == 22
    assert (df["high_net_per_substrate"] > df["low_net_per_substrate"]).all()
    # Basal carbon-entry dependence is architecture-specific: it remains strong in
    # E. coli, S. enterica and K. phaffii, whereas the tested B. subtilis and
    # S. cerevisiae hexose states largely collapse after stoichiometric normalization.
    for org in ("E. coli", "S. enterica", "K. phaffii"):
        group=df[df.organism==org]
        assert group["low_net_per_substrate"].max() - group["low_net_per_substrate"].min() > 0.5
    se_glcn = df[(df.organism == "S. enterica") & (df.substrate == "gluconate")].iloc[0]
    assert abs(se_glcn.low_net_per_substrate - 2.1385686678) < 1e-6


def test_cross_species_principle_feasibility_is_explicit():
    df = pd.read_csv(FEAS)
    # Infeasible combinations are retained in the audit rather than silently dropped.
    assert ((df.organism == "B. subtilis") & (df.substrate == "gluconate") & (df.state == "low-disposal-capacity") & (df.status != "feasible")).any()
    assert ((df.organism == "K. phaffii") & (df.substrate == "galactose") & (df.status != "feasible")).any()
