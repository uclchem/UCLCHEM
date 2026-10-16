import numpy as np
import pytest

from uclchem.utils import get_dtype, get_fuzzy_unknown_keys

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


test_get_fuzzy_unknown_keys_data = [
    (["abcd"], ["abce", "d"], {("abcd", "abce")}),
    (["abcd"], ["abcd", "d"], set()),
    (["abcdE"], ["abcde", "d"], {("abcdE", "abcde")}),
]


@pytest.mark.parametrize(
    ("keys", "expected_keys", "expected_unknown"), test_get_fuzzy_unknown_keys_data
)
def test_get_fuzzy_unknown_keys(keys, expected_keys, expected_unknown):
    unknown_keys = get_fuzzy_unknown_keys(keys, expected_keys)
    assert unknown_keys == expected_unknown
