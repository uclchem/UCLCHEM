import pytest

import uclchem


def test_ode_conservation():
    result = uclchem.tests.test_ode_conservation()

    for key, value in result.items():
        print(f"{key} {value:.2e}")
    assert (
        result["H"] < 1e-15
    ), f"H not conserved with total rate of change {result['H']:.2e}"


def main():

    # Run the tests using pytest
    pytest.main(["-v", __file__])


if __name__ == "__main__":
    main()
