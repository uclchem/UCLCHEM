import numpy as np
import uclchemwrap


def test_evaluate_polynomial():
    for size in range(1, 10):
        coeffs = np.random.random(size)
        x = np.random.random()
        y_horner = uclchemwrap.numerics.evaluate_polynomial(coeffs, x)

        assert y_horner == np.polynomial.Polynomial(coeffs)(x)


def test_integrate_trapezoid():
    for size in range(2, 10):
        x, y = np.split(np.random.random((2, size)), 2)
        x = np.sort(x)
        integral = uclchemwrap.numerics.integrate_trapezoid(x, y)

        assert np.isclose(integral, np.trapezoid(y, x))


def test_logspace():
    for size in range(2, 10):
        start, stop = np.random.random(2)
        x = uclchemwrap.numerics.logspace(start, stop, size)

        assert np.allclose(x, np.logspace(start, stop, num=size))
