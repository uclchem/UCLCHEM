import numpy as np
import pytest

from uclchem.makerates.utils import (
    array_to_string,
    pad_to_max_length,
)


def test_array_to_string_1d_int():
    arr = np.array([1, 2, 3, 4])
    result = array_to_string("arr1", arr, value_type="int", parameter=True)
    assert "integer, parameter :: arr1 (4)=(/1,2,3,4/)" in result.replace("\n", "")


def test_array_to_string_1d_float():
    arr = np.array([1.0, 2.0, 3.0])
    result = array_to_string("arr2", arr, value_type="float", parameter=True)
    assert (
        "real(dp), parameter :: arr2 (3)=(/1.0000e+00_dp,2.0000e+00_dp,3.0000e+00_dp/)"
        in result.replace("\n", "")
    )


def test_array_to_string_2d_int():
    arr = np.array([[1, 2, 3], [4, 5, 6]])
    result = array_to_string("arr3", arr, value_type="int", parameter=True)
    expected = """integer, parameter :: arr3(2,3) = RESHAPE((/ 1,4,2,5,3,6 /), (/ 2,3 /))
"""
    assert result == expected


def test_array_to_string_2d_ones():
    arr = np.ones((5, 7), dtype=int)
    result = array_to_string("arr_ones", arr, value_type="int", parameter=True)
    expected = """integer, parameter :: arr_ones(5,7) = RESHAPE((/ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1&
    &,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 /), (/ 5,7 /))
"""
    assert result == expected


def test_array_to_string_2d_float():
    arr = np.array([[1.0, 2.0], [3.0, 4.0]])
    result = array_to_string("arr4", arr, value_type="float", parameter=True)
    expected = """real(dp), parameter :: arr4(2,2) = RESHAPE((/ 1.0000e+00_dp,3.0000e+00_dp&
    &,2.0000e+00_dp,4.0000e+00_dp /), (/ 2,2 /))
"""
    assert result == expected


def test_array_to_string_2d_string():
    arr = np.array([["A", "B"], ["C", "D"]])
    result = array_to_string("arr5", arr, value_type="string", parameter=True)
    expected = """character(LEN=1), parameter :: arr5(2,2) = RESHAPE((/ "A","C","B","D" /), (/ 2&
    &,2 /))
"""
    assert result == expected


pad_to_max_length_data = [
    (["string"], ["string"]),
    (["short", "very long"], ["short    ", "very long"]),
    (["", "short"], ["     ", "short"]),
    (["", " "], [" ", " "]),
]


@pytest.mark.parametrize(("strings", "expected_strings"), pad_to_max_length_data)
def test_pad_to_max_length(strings, expected_strings):
    padded_strings = pad_to_max_length(strings)
    assert all(len(string) == len(padded_strings[0]) for string in padded_strings)
    assert padded_strings == expected_strings
