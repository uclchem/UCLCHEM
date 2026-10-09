import numpy as np
import pytest

from uclchem.utils import (
    get_dtype,
    get_lowercase_copy,
    remove_keys_with_none,
    would_overflow,
)

get_dtypes_data = [
    ("fp64", np.float64),
    ("fp32", np.float32),
    ("fp16", np.float16),
    (np.float64, np.float64),
    (np.dtype(np.float64), np.dtype(np.float64)),
]


@pytest.mark.parametrize(("input_dtype", "expected"), get_dtypes_data)
def test_get_dtype(input_dtype, expected):
    assert get_dtype(input_dtype) == expected


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


get_lowercase_copy_data = [
    ({}, {}),
    ({"a": "b"}, {"a": "b"}),
    ({"a": "A"}, {"a": "A"}),
    ({"A": "A"}, {"a": "A"}),
    ({"A": "A", "b": "B"}, {"a": "A", "b": "B"}),
    ({"A": "A", "B": "B"}, {"a": "A", "b": "B"}),
]


@pytest.mark.parametrize(("dct", "expected_dct"), get_lowercase_copy_data)
def test_get_lowercase_copy(dct, expected_dct):
    assert get_lowercase_copy(dct) == expected_dct


remove_keys_with_none_data = [
    ({}, {}),
    ({"a": None}, {}),
    ({"a": "None"}, {"a": "None"}),
    ({"a": ""}, {"a": ""}),
    ({"a": None, "b": "B"}, {"b": "B"}),
]


@pytest.mark.parametrize(("dct", "expected_dct"), remove_keys_with_none_data)
def test_remove_keys_with_none(dct, expected_dct):
    assert remove_keys_with_none(dct) == expected_dct
