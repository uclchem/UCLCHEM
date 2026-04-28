import pytest
import numpy as np
from uclchem.utils import get_dtype

get_dtypes_data = [
    ("fp64", np.float64),
    ("fp32", np.float32),
    ("fp16", np.float16),
    (np.float64, np.float64),
    (np.dtype(np.float64), np.dtype(np.float64)),
]


@pytest.mark.parametrize("input, expected", get_dtypes_data)
def test_get_dtype(input, expected):
    assert get_dtype(input) == expected
