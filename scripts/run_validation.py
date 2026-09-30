"""
run_validation.py
--------------------
Runs the complete validation suite in sequence:
  1. Reduction-degree formula (reference compounds, network-wide
     consistency check) -- validate_reduction_degree.py
  2. 52-reaction categorical validation battery spanning the network's
     chemical diversity -- validation_battery.py
  3. Classification coverage (fraction of reactions successfully
     classified, and the precise breakdown of exclusion reasons) --
     coverage_analysis.py
  4. Matching specificity (group-transfer reactions, aminotransferase
     screen, invalid-inorganic-pair screen, seed-pair integrity) --
     matching_validation.py

Each step is also independently runnable as its own script. Writes every
validation result table to results/.
"""
import os
import runpy

SCRIPTS_DIR = os.path.dirname(__file__)

VALIDATION_STEPS = [
    "validate_reduction_degree.py",
    "validation_battery.py",
    "coverage_analysis.py",
    "matching_validation.py",
]

if __name__ == "__main__":
    for step in VALIDATION_STEPS:
        print("\n" + "=" * 78)
        print(f"{step}")
        print("=" * 78)
        runpy.run_path(os.path.join(SCRIPTS_DIR, step), run_name="__main__")
    print("\n" + "=" * 78)
    print("Validation suite complete.")
    print("=" * 78)
