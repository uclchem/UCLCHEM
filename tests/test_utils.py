import numpy as np
import pytest

from uclchem.utils import get_dtype, would_overflow

get_dtypes_data = [
    ("fp64", np.float64),
    ("fp32", np.float32),
    ("fp16", np.float16),
    (np.float64, np.float64),
    (np.dtype(np.float64), np.dtype(np.float64)),
]


@pytest.mark.parametrize(("input", "expected"), get_dtypes_data)
def test_get_dtype(input, expected):
    assert get_dtype(input) == expected


overflow_data = [
    (0, np.float16, False),
    (0, np.float32, False),
    (0, np.float64, False),
    (np.finfo(np.float64).max, np.float64, False),
    (np.finfo(np.float64).max, np.float32, True),
    (np.finfo(np.float32).max, np.float16, True),
    (-np.finfo(np.float64).max, np.float64, False),
    (-np.finfo(np.float64).max, np.float32, True),
    (-np.finfo(np.float32).max, np.float16, True),
]


@pytest.mark.parametrize(("value", "dtype", "expected"), overflow_data)
def test_would_overflow(value, dtype, expected):
    assert would_overflow(value, dtype) == expected
