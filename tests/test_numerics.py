import numpy as np
import pytest
import uclchemwrap

is_equal_data = [
    (0.3, 0.3, None, True),
    (0.3, 0.1 + 0.2, 1e-10, True),
    (0.3, 0.1 + 0.2, 0, False),
]


@pytest.mark.parametrize(("a", "b", "atol", "expected"), is_equal_data)
def test_is_equal(a, b, atol, expected):
    assert bool(uclchemwrap.numerics.is_equal(a, b, atol=atol)) == expected


def test_evaluate_polynomial():
    rng = np.random.default_rng()
    for size in range(1, 10):
        coeffs = rng.random(size)
        x = rng.random()
        y_horner = uclchemwrap.numerics.evaluate_polynomial(coeffs, x)

        assert y_horner == np.polynomial.Polynomial(coeffs)(x)


def test_integrate_trapezoid():
    rng = np.random.default_rng()
    for size in range(2, 10):
        x, y = np.split(rng.random((2, size)), 2)
        x = np.sort(x)
        integral = uclchemwrap.numerics.integrate_trapezoid(x, y)

        assert np.isclose(integral, np.trapezoid(y, x))


def test_logspace():
    rng = np.random.default_rng()
    for size in range(2, 10):
        start, stop = rng.random(2)
        x = uclchemwrap.numerics.logspace(start, stop, size)

        assert np.allclose(x, np.logspace(start, stop, num=size))


def test_pair_insertion_sort():
    rng = np.random.default_rng()
    for size in range(2, 10):
        x = rng.random(size)
        x_original = x.copy()
        uclchemwrap.numerics.pair_insertion_sort(x)

        x_sort = np.sort(x_original)
        assert np.all(x == x_sort)


def test_pair_insertion_sort_with_perm():
    rng = np.random.default_rng()
    for size in range(2, 10):
        x = rng.random(size)
        x_original = x.copy()

        perms = uclchemwrap.numerics.pair_insertion_sort_with_perm(x)

        x_argsort = np.argsort(x_original)
        x_sort = np.sort(x_original)

        assert np.all(x == x_sort)
        assert np.all(perms == x_argsort + 1)
