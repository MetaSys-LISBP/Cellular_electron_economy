"""
conftest.py
------------
Shared pytest fixtures (session-scoped, loaded once for the whole test
suite): the bundled iML1515 model, and its metabolite table (formula,
charge, gamma). No network access required.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pytest
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")


@pytest.fixture(scope="session")
def model():
    pytest.importorskip("cobra", reason="COBRApy is required for model-based tests")
    from etn.data_acquisition import load_ecoli_model
    return load_ecoli_model()


@pytest.fixture(scope="session")
def met_table(model):
    from etn.matching import build_metabolite_table
    return build_metabolite_table(model)
