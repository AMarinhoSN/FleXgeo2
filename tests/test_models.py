from __future__ import annotations

import subprocess
import sys
import typing
from dataclasses import fields

import pandas as pd
import pytest

from flexgeo2 import models
from flexgeo2.models import (
    AnalysisResult,
    DistanceResult,
    ResidueClusteringResult,
    ResidueRangeClusteringResult,
)


@pytest.mark.parametrize(
    "result_class",
    [AnalysisResult, DistanceResult, ResidueClusteringResult, ResidueRangeClusteringResult],
)
def test_result_tables_are_typed_as_data_frames(result_class: type) -> None:
    # pandas is imported for type checkers only, so resolve the annotations with it.
    hints = typing.get_type_hints(result_class, globalns={**vars(models), "pd": pd})
    tables = [field.name for field in fields(result_class) if field.name.endswith("_df")]

    assert tables
    assert {name: hints[name] for name in tables} == dict.fromkeys(tables, pd.DataFrame)


def test_importing_flexgeo2_does_not_import_pandas() -> None:
    # Keeps `import flexgeo2` and `flexgeo2 --help` fast; pandas loads when a run needs it.
    code = "import sys, flexgeo2; print('pandas' in sys.modules)"
    completed = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )

    assert completed.stdout.strip() == "False"
