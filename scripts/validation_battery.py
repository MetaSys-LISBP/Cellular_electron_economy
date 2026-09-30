"""
validation_battery.py
------------------------
52-reaction categorical validation battery designed to
cover the full diversity of redox and non-redox chemistry present in the
network, and reporting standard classification-performance metrics
(accuracy, precision, recall, specificity, F1) rather than a single
"N/N correct" headline number.

IMPORTANT, EXPLICIT LIMITATION:
this is a HAND-CURATED reference set, built
and labelled by the same team that built the classification method. It
is NOT a blind, independently sourced gold standard, and should not be
presented as one. It is deliberately built to provide broad and adversarial coverage of
NAD(P)-linked
dehydrogenases, flavoprotein reductases, quinone-linked dehydrogenases,
an oxidase, transhydrogenases, a ribonucleotide reductase, isomerases,
kinases, a hydratase, transporters, a group-transfer/non-redox
phosphatase, and -- importantly -- three reactions that ARE genuinely
redox but involve an iron centre this method cannot score, included to
quantify the false-negative rate of the metal exclusion honestly rather
than omit it).

Ground-truth categories:
  'redox'              -- genuinely redox, chemically scorable by this
                           method (no metal/metalloid/generic group)
  'non_redox'          -- genuinely non-redox, chemically scorable
  'redox_unscorable'   -- genuinely redox, but involves an element this
                           method deliberately does not score (e.g. Fe);
                           the CORRECT behaviour is status=missing_formula,
                           not a redox/non-redox guess
  'non_redox_unscorable' -- genuinely non-redox, but also involves an
                           unscorable element (metal chelation etc.)

Confusion-matrix metrics (accuracy, precision, recall, specificity, F1)
are computed ONLY over the chemically SCORABLE subset ('redox' and
'non_redox'), since the 'unscorable' categories are, by design, meant to
return "cannot classify" rather than a binary call -- see the coverage
analysis (scripts/coverage_analysis.py) for how the unscorable fraction
is characterised separately.
"""
import os
import sys
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)

import pandas as pd
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from etn.data_acquisition import load_ecoli_model
from etn.matching import build_metabolite_table, classify_reaction_molecular

RESULTS_DIR = os.path.join(REPO_ROOT, "results")

# (reaction_id, ground_truth_category, expected_electrons_or_None, category_label)
VALIDATION_SET = [
    # -- isomerases (non-redox) --
    ("PGI", "non_redox", None, "isomerase"),
    ("TPI", "non_redox", None, "isomerase"),
    ("RPI", "non_redox", None, "isomerase"),
    # -- kinases (non-redox) --
    ("PFK", "non_redox", None, "kinase"),
    ("PYK", "non_redox", None, "kinase"),
    ("HEX1", "non_redox", None, "kinase"),
    ("ADK1", "non_redox", None, "kinase"),
    # -- hydratase (non-redox, one bond from a redox pair) --
    ("FUM", "non_redox", None, "hydratase (near redox pair)"),
    # -- transporters (non-redox) --
    ("PIuabcpp", "non_redox", None, "transporter"),
    ("NO3t7pp", "non_redox", None, "transporter (antiport of a redox pair)"),
    ("NH4tpp", "non_redox", None, "transporter"),
    ("GLCptspp", "non_redox", None, "transporter (PTS, group translocation)"),
    # -- group transfer, non-redox --
    ("PPA", "non_redox", None, "phosphate group transfer (pyrophosphatase)"),
    # -- NAD(P)-dependent dehydrogenases (redox, 2e-) --
    ("ALCD2x", "redox", 2, "NAD-dehydrogenase"),
    ("GAPD", "redox", 2, "NAD-dehydrogenase + phosphoryl transfer"),
    ("MDH", "redox", 2, "NAD-dehydrogenase"),
    ("G6PDH2r", "redox", 2, "NADP-dehydrogenase"),
    ("ICDHyr", "redox", 2, "NADP-dehydrogenase"),
    ("LDH_D", "redox", 2, "NAD-dehydrogenase"),
    ("HACD1", "redox", 2, "NAD-dehydrogenase (beta-oxidation)"),
    ("ALCD19", "redox", 2, "NAD-dehydrogenase"),
    ("GLUDy", "redox", 2, "NADP-dehydrogenase"),
    # -- decarboxylating dehydrogenase complexes (redox + decarboxylation) --
    ("AKGDH", "redox", 2, "decarboxylating dehydrogenase"),
    ("PDH", "redox", 2, "decarboxylating dehydrogenase + CoA transfer"),
    # -- quinone-linked dehydrogenases (redox, 2e-) --
    ("SUCDi", "redox", 2, "quinone-dehydrogenase"),
    ("MDH2", "redox", 2, "quinone-dehydrogenase"),
    ("MDH3", "redox", 2, "menaquinone-dehydrogenase"),
    ("GLYCTO2", "redox", 2, "quinone-dehydrogenase"),
    ("POX", "redox", 2, "quinone-linked oxidative decarboxylation"),
    ("DHORD2", "redox", 2, "quinone-dehydrogenase"),
    ("DHORD5", "redox", 2, "menaquinone-dehydrogenase"),
    ("FDH4pp", "redox", 2, "quinone-dehydrogenase (formate/CO2)"),
    ("FDH5pp", "redox", 2, "menaquinone-dehydrogenase (formate/CO2)"),
    ("FRD2", "redox", 2, "menaquinol-linked reductase"),
    ("FRD3", "redox", 2, "menaquinol-linked reductase"),
    # -- flavoprotein reductases (redox, 2e-, non-NAD terminal step) --
    ("FADRx", "redox", 2, "flavoprotein reductase"),
    ("FMNRx2", "redox", 2, "flavoprotein reductase"),
    ("GTHOr", "redox", 2, "glutathione reductase"),
    # -- respiratory chain (redox, 2e-) --
    ("NADH16pp", "redox", 2, "respiratory complex"),
    ("CYTBO3_4pp", "redox", 2, "terminal oxidase"),
    ("CYTBDpp", "redox", 2, "terminal oxidase"),
    ("CYTBD2pp", "redox", 2, "terminal oxidase"),
    # -- anaerobic respiration --
    ("NO3R1pp", "redox", 2, "anaerobic respiration"),
    ("NTRIR2x", "redox", 6, "anaerobic respiration, 6-electron"),
    # -- transhydrogenases (redox, cofactor-to-cofactor) --
    ("THD2pp", "redox", 2, "transhydrogenase"),
    ("NADTRHD", "redox", 2, "transhydrogenase"),
    # -- O2-linked oxidase --
    ("ASPO6", "redox", 2, "amino-acid oxidase"),
    # -- ribonucleotide reductase: genuinely redox, but thioredoxin is
    #    represented with a generic "X" formula group in this model
    #    (not a real small-molecule formula), so it is correctly
    #    unscorable for the same reason as a metal; see the coverage
    #    breakdown (Methods) --
    ("RNDR1", "redox_unscorable", None, "ribonucleotide reductase (generic-formula cofactor, thioredoxin)"),
    # -- other dehydrogenase --
    ("SHCHD2", "redox", 2, "NAD-dehydrogenase (tetrapyrrole biosynthesis)"),
    # -- metal-centred redox: genuinely redox, but NOT scorable by this
    #    method (Fe not in SAFE_ELEMENTS) -- correct behaviour is
    #    status=missing_formula, quantified separately, not guessed --
    ("FEROpp", "redox_unscorable", None, "Fe(II)/Fe(III) redox (metal-excluded)"),
    ("FE3Ri", "redox_unscorable", None, "Fe(II)/Fe(III) redox (metal-excluded)"),
    # -- metal-involving, non-redox (chelation) --
    ("SHCHF", "non_redox_unscorable", None, "metal chelation (metal-excluded)"),
]


def run_validation_battery(model, met_table) -> pd.DataFrame:
    rows = []
    for rid, truth, expected_e, category in VALIDATION_SET:
        rxn = model.reactions.get_by_id(rid)
        res = classify_reaction_molecular(rxn, met_table)
        donors = [p["reactant"] for p in res["pairs"] if p["delta_gamma"] < 0]
        acceptors = [p["product"] for p in res["pairs"] if p["delta_gamma"] > 0]
        n_electrons = res.get("electrons_per_turnover_estimate", 0.0)

        if truth in ("redox_unscorable", "non_redox_unscorable"):
            predicted_ok = res["status"] == "missing_formula"
            failure_mode = "" if predicted_ok else f"expected missing_formula, got status={res['status']}"
        else:
            expected_redox = (truth == "redox")
            predicted_ok = (res["status"] == "ok") and (res["is_redox"] == expected_redox) and (
                expected_e is None or n_electrons == expected_e)
            if not predicted_ok:
                if res["status"] != "ok":
                    failure_mode = f"unexpectedly unscored (status={res['status']})"
                elif res["is_redox"] != expected_redox:
                    failure_mode = "wrong redox/non-redox call"
                else:
                    failure_mode = f"electron count mismatch (got {n_electrons}, expected {expected_e})"
            else:
                failure_mode = ""

        rows.append({
            "reaction_id": rid,
            "category": category,
            "ground_truth": truth,
            "predicted_status": res["status"],
            "predicted_is_redox": res["is_redox"],
            "predicted_electrons": n_electrons,
            "donor(s)": ";".join(donors),
            "acceptor(s)": ";".join(acceptors),
            "correct": predicted_ok,
            "failure_mode": failure_mode,
        })
    return pd.DataFrame(rows)


def confusion_matrix_metrics(df: pd.DataFrame) -> dict:
    """Computed ONLY on the chemically scorable subset (ground_truth in
    {'redox', 'non_redox'}); see module docstring."""
    scorable = df[df.ground_truth.isin(["redox", "non_redox"])].copy()
    tp = ((scorable.ground_truth == "redox") & (scorable.predicted_is_redox) & (scorable.predicted_status == "ok")).sum()
    fn = ((scorable.ground_truth == "redox") & ~((scorable.predicted_is_redox) & (scorable.predicted_status == "ok"))).sum()
    tn = ((scorable.ground_truth == "non_redox") & (~scorable.predicted_is_redox) & (scorable.predicted_status == "ok")).sum()
    fp = ((scorable.ground_truth == "non_redox") & (scorable.predicted_is_redox) & (scorable.predicted_status == "ok")).sum()
    n = len(scorable)
    accuracy = (tp + tn) / n if n else float("nan")
    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    recall = tp / (tp + fn) if (tp + fn) else float("nan")
    specificity = tn / (tn + fp) if (tn + fp) else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else float("nan")

    unscorable = df[df.ground_truth.isin(["redox_unscorable", "non_redox_unscorable"])]
    unscorable_correct = unscorable["correct"].sum()

    return {
        "n_scorable": int(n), "TP": int(tp), "FN": int(fn), "TN": int(tn), "FP": int(fp),
        "accuracy": accuracy, "precision": precision, "recall_sensitivity": recall,
        "specificity": specificity, "F1": f1,
        "n_unscorable_cases": int(len(unscorable)),
        "n_unscorable_correctly_flagged": int(unscorable_correct),
    }


if __name__ == "__main__":
    model = load_ecoli_model()
    met_table = build_metabolite_table(model)
    df = run_validation_battery(model, met_table)
    df.to_csv(f"{RESULTS_DIR}/validation_battery_table.csv", index=False)

    print(f"Expanded validation set: {len(df)} reactions across "
          f"{df.category.nunique()} categories")
    print()
    wrong = df[~df.correct]
    if len(wrong):
        print("INCORRECT predictions:")
        print(wrong[["reaction_id", "category", "ground_truth", "predicted_status",
                      "predicted_is_redox", "failure_mode"]].to_string(index=False))
    else:
        print("All predictions correct on both the scorable and unscorable subsets.")
    print()

    # This validation script is part of the executable QC suite: any wrong
    # classification or expected electron count must fail the command rather
    # than merely print a warning.
    assert wrong.empty, (
        f"52-reaction validation battery failed for {len(wrong)} reaction(s): "
        + ", ".join(wrong.reaction_id.astype(str).tolist())
    )

    metrics = confusion_matrix_metrics(df)
    for k, v in metrics.items():
        print(f"  {k}: {v:.4f}" if isinstance(v, float) else f"  {k}: {v}")

    import json
    with open(f"{RESULTS_DIR}/validation_battery_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2, default=float)
