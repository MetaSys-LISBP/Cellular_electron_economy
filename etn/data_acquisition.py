"""
data_acquisition.py
--------------------
Retrieves the complete E. coli genome-scale metabolic network (reactions,
metabolites) and associates a SMILES/InChI structure with each metabolite
where available. Also exports the full network in SBML format.

Sources:
  - Metabolic network (reactions, stoichiometry, reversibility, GPR,
    subsystems): iML1515 (2,712 reactions, 1,877 metabolites, 1,516 genes;
    Monk et al. 2017, Nat Biotechnol 35:904-908), the most recent and most
    complete BiGG genome-scale reconstruction of E. coli K-12 MG1655. The
    SBML file is bundled locally at `data/iML1515.xml` (~9.2 MB) so that
    every analysis in this repository runs fully OFFLINE, with no
    dependency on the BiGG Models web service at runtime. This is the
    single model used throughout this study; no other reconstruction
    is used for any reported result.
  - Molecular structures (SMILES, InChI, InChIKey): the ModelSEED
    Biochemistry Database (Henry et al. 2010, Nat Biotechnol 28:977-982;
    Seaver et al. 2021, Nucleic Acids Res 49:D575-D588;
    https://github.com/ModelSEED/ModelSEEDDatabase), cross-referenced via
    the `seed.compound` annotation already present on every metabolite of
    BiGG models. The full InChI (ModelSEED provides only SMILES and
    InChIKey, not the full InChI string) is then recomputed from the
    SMILES with RDKit. Structural annotations are optional and are not
    required by the formula-and-charge electron-transfer reconstruction.
    Retrieval requires network access on first call (see data/README.md to
    cache it locally); all reported reconstruction results use the bundled
    genome-scale model and its formula/charge annotations.

Coverage achieved (measured empirically on iML1515): ~90% of metabolites
carry a usable SMILES via this route (~81% once metabolites lacking any
`seed.compound` cross-reference are also mis-typed values are filtered
out). The remainder (mostly ill-defined species / generic macromolecules
/ lipids with unspecified chain length) are explicitly flagged
(smiles=None) rather than silently dropped.
"""

import os
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

MODELSEED_COMPOUNDS_URL = (
    "https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/"
    "master/Biochemistry/compounds.tsv"
)

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
LOCAL_MODEL_PATH = os.path.join(REPO_ROOT, "data", "iML1515.xml")


def load_ecoli_model(local_path: str = None):
    """
    Loads the iML1515 E. coli genome-scale model from the local, bundled
    SBML file (data/iML1515.xml by default) via cobrapy -- fully offline,
    no BiGG Models web access required. `local_path` can point to an
    alternative SBML file if desired, but this repository's results were all
    produced with the
    bundled file.
    """
    import cobra

    path = local_path or LOCAL_MODEL_PATH
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"iML1515 SBML file not found at {path}. This repository expects "
            f"the model to be bundled at data/iML1515.xml -- see data/README.md "
            f"for how to (re)download it if it is missing."
        )
    model = cobra.io.read_sbml_model(path)
    print(f"[data_acquisition] Loaded model: {model.id} "
          f"({len(model.reactions)} reactions, {len(model.metabolites)} metabolites, "
          f"{len(model.genes)} genes) from {path}")
    return model


def load_modelseed_structures(local_path: str = None) -> pd.DataFrame:
    """
    Loads the ModelSEED structure table (id, name, formula, charge,
    smiles, inchikey, ...). Uses a local file if provided AND it exists,
    otherwise downloads it from GitHub (raw.githubusercontent.com).
    """
    if local_path and os.path.isfile(local_path):
        source = local_path
    else:
        source = MODELSEED_COMPOUNDS_URL
    print(f"[data_acquisition] Loading ModelSEED structures from: {source}")
    df = pd.read_csv(source, sep="\t", dtype=str)
    return df


def _get_seed_ids(met) -> list:
    """Retrieves the ModelSEED ID(s) (`seed.compound`) of a cobra metabolite."""
    ann = met.annotation.get("seed.compound", [])
    if isinstance(ann, str):
        return [ann]
    return list(ann) if ann else []


def build_metabolites_table(model, modelseed_df: pd.DataFrame) -> pd.DataFrame:
    """
    Builds the model's metabolite table, enriched with SMILES / InChI /
    InChIKey by cross-referencing ModelSEED (seed.compound annotation).
    """
    smiles_map = dict(zip(modelseed_df["id"], modelseed_df["smiles"]))
    inchikey_map = dict(zip(modelseed_df["id"], modelseed_df["inchikey"]))

    rows = []
    for met in model.metabolites:
        seed_ids = _get_seed_ids(met)
        smiles, inchikey, seed_used = None, None, None
        for sid in seed_ids:
            s = smiles_map.get(sid)
            if s and isinstance(s, str) and s.lower() != "null":
                smiles, seed_used = s, sid
                inchikey = inchikey_map.get(sid)
                break

        inchi = None
        if smiles:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    inchi = Chem.MolToInchi(mol)
                    if not inchikey or inchikey.lower() == "null":
                        inchikey = Chem.InchiToInchiKey(inchi)
            except Exception:
                pass

        rows.append({
            "met_id": met.id,
            "name": met.name,
            "formula": met.formula,
            "charge": met.charge,
            "compartment": met.compartment,
            "seed_compound_id": seed_used,
            "smiles": smiles,
            "inchi": inchi,
            "inchikey": inchikey,
            "has_structure": smiles is not None,
        })

    return pd.DataFrame(rows)


def build_reactions_table(model) -> pd.DataFrame:
    """
    Builds the model's reaction table: id, name, equation (in met_id),
    flux bounds, reversibility, gene-protein-reaction rule (GPR),
    subsystem, annotations (EC, KEGG, BiGG...).
    """
    rows = []
    for rxn in model.reactions:
        rows.append({
            "reaction_id": rxn.id,
            "name": rxn.name,
            "equation": rxn.build_reaction_string(use_metabolite_names=False),
            "equation_readable": rxn.build_reaction_string(use_metabolite_names=True),
            "lower_bound": rxn.lower_bound,
            "upper_bound": rxn.upper_bound,
            "reversible": rxn.reversibility,
            "subsystem": rxn.subsystem,
            "gene_reaction_rule": rxn.gene_reaction_rule,
            "ec_number": rxn.annotation.get("ec-code"),
            "kegg_reaction": rxn.annotation.get("kegg.reaction"),
            "metanetx_reaction": rxn.annotation.get("metanetx.reaction"),
        })
    return pd.DataFrame(rows)


def _canonicalise_sbml_group_member_order(path: str) -> None:
    """Sort SBML Groups members by ``idRef`` for byte-stable exports.

    COBRApy stores ``Group.members`` as sets, so the biological content is
    deterministic but the serialized member order can vary between Python
    processes. Sorting only the self-closing ``groups:member`` lines inside
    each ``groups:listOfMembers`` block preserves SBML semantics while making
    repeated exports byte-reproducible.
    """
    import re
    from pathlib import Path

    sbml_path = Path(path)
    text = sbml_path.read_text(encoding="utf-8")
    block_re = re.compile(
        r"(?P<open>^[ \t]*<groups:listOfMembers>\n)"
        r"(?P<body>(?:^[ \t]*<groups:member\b[^\n]*/>\n)+)"
        r"(?P<close>^[ \t]*</groups:listOfMembers>)",
        flags=re.MULTILINE,
    )

    def sort_block(match):
        lines = match.group("body").splitlines(keepends=True)

        def member_key(line):
            id_match = re.search(r'groups:idRef="([^"]+)"', line)
            return (id_match.group(1) if id_match else line.strip(), line)

        return match.group("open") + "".join(sorted(lines, key=member_key)) + match.group("close")

    canonical = block_re.sub(sort_block, text)
    if canonical != text:
        sbml_path.write_text(canonical, encoding="utf-8")


def export_sbml(model, path: str):
    """Export the full metabolic network in byte-stable SBML (FBC v2)."""
    from cobra.io import write_sbml_model

    write_sbml_model(model, path)
    _canonicalise_sbml_group_member_order(path)
    print(f"[data_acquisition] SBML network written to: {path}")


def acquire_all(local_modelseed_path: str = None, local_model_path: str = None,
                 include_structures: bool = False):
    """
    Main entry point for network + (optionally) structure acquisition.
    Returns (model, reactions_df, metabolites_df).

    include_structures=False (the default) skips ModelSEED SMILES/InChI
    acquisition entirely: those annotations are used only for optional
    structure-enriched analyses, not by the formula/charge-based reconstruction
    that produces the reported results in this study. This keeps the primary pipeline
    (scripts/run_analysis.py) reproducible fully offline from the
    bundled data/iML1515.xml alone: data/modelseed_compounds.tsv is not
    bundled (it is a large, optional cache -- see data/README.md), so
    calling load_modelseed_structures() unconditionally would silently
    attempt a network download (pd.read_csv on a GitHub raw URL) when
    that local file is absent. Pass include_structures=True explicitly
    (and ideally a cached local_modelseed_path) if you specifically want
    the SMILES-enriched metabolite table for structure-enriched analyses.
    """
    model = load_ecoli_model(local_path=local_model_path)
    if include_structures:
        modelseed_df = load_modelseed_structures(local_modelseed_path)
        metabolites_df = build_metabolites_table(model, modelseed_df)
        n_ok = metabolites_df["has_structure"].sum()
        n_tot = len(metabolites_df)
        print(f"[data_acquisition] Structures (SMILES/InChI) found for "
              f"{n_ok}/{n_tot} metabolites ({100*n_ok/n_tot:.1f}%)")
    else:
        metabolites_df = pd.DataFrame(
            [{"met_id": m.id, "name": m.name, "formula": m.formula, "charge": m.charge,
              "compartment": m.compartment} for m in model.metabolites])
    reactions_df = build_reactions_table(model)

    return model, reactions_df, metabolites_df


if __name__ == "__main__":
    model, reactions_df, metabolites_df = acquire_all()
    reactions_df.to_csv("reactions_ecoli.csv", index=False)
    metabolites_df.to_csv("metabolites_ecoli.csv", index=False)
    export_sbml(model, "ecoli_network.xml")
    print(f"[data_acquisition] Wrote: reactions_ecoli.csv ({len(reactions_df)} rows), "
          f"metabolites_ecoli.csv ({len(metabolites_df)} rows)")
