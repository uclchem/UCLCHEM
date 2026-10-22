import pytest

import uclchem


def test_ode_conservation(tolerance: float = 1e-15):
    result = uclchem.tests.test_ode_conservation()
    assert all(discrep < tolerance for discrep in result.values()), (
        f"ODE not conserved with total rate of change {result}"
    )


def main():
    # Run the tests using pytest
    pytest.main(["-v", __file__])


if __name__ == "__main__":
    main()
