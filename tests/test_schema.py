import os
import sys

import pandas as pd
import pytest

path = os.path.normpath(os.path.dirname(__file__))

sys.path.insert(0, os.path.join(path, os.pardir, "bin"))
from check_samplesheet_v2 import schema  # noqa: E402


@pytest.fixture
def samples():
    data = pd.DataFrame([["s1", "bam1.bam"], ["s2", "bam1.bam"]])
    data.columns = ["sample", "bam"]
    return data


def test_schema(samples):
    rows = samples.to_dict("records")
    schema.validate(rows)
