import numpy as np
import pytest

from uclchem.makerates.io_functions import array_to_string


def test_array_to_string_1d_int():
    arr = np.array([1, 2, 3, 4])
    result = array_to_string("arr1", arr, type="int", parameter=True)
    assert "INTEGER(dp), PARAMETER :: arr1 (4)=(/1,2,3,4/)" in result.replace("\n", "")


def test_array_to_string_1d_float():
    arr = np.array([1.0, 2.0, 3.0])
    result = array_to_string("arr2", arr, type="float", parameter=True)
    assert (
        "REAL(dp), PARAMETER :: arr2 (3)=(/1.0000e+00,2.0000e+00,3.0000e+00/)"
        in result.replace("\n", "")
    )


def test_array_to_string_2d_int():
    arr = np.array([[1, 2, 3], [4, 5, 6]])
    result = array_to_string("arr3", arr, type="int", parameter=True)
    expected = """INTEGER(dp), PARAMETER :: arr3(2,3) = RESHAPE((/ 1,2,3,4,5,6 /), (/ 2&
    &, 3 /))
"""
    assert result == expected


def test_array_to_string_2d_ones():
    arr = np.ones((5, 7), dtype=int)
    result = array_to_string("arr_ones", arr, type="int", parameter=True)
    expected = """INTEGER(dp), PARAMETER :: arr_ones(5,7) = RESHAPE((/ 1,1,1,1,1,1,1,1,1,1&
    &,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 /), (/ 5, 7 /))
"""
    assert result == expected


def test_array_to_string_2d_float():
    arr = np.array([[1.0, 2.0], [3.0, 4.0]])
    result = array_to_string("arr4", arr, type="float", parameter=True)
    expected = """REAL(dp), PARAMETER :: arr4(2,2) = RESHAPE((/ 1.0000e+00,2.0000e+00&
    &,3.0000e+00,4.0000e+00 /), (/ 2, 2 /))
"""
    assert result == expected


def test_array_to_string_2d_string():
    arr = np.array([["A", "B"], ["C", "D"]])
    result = array_to_string("arr5", arr, type="string", parameter=True)
    expected = """CHARACTER(Len=1), PARAMETER :: arr5(2,2) = RESHAPE((/ "A","B","C","D" /)&
    &, (/ 2, 2 /))
"""
    assert result == expected
