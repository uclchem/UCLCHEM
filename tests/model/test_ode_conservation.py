import uclchem


def test_ode_conservation():
    result = uclchem.tests.test_ode_conservation()
    assert all(discrep < 1e-15 for discrep in result.values()), (
        f"ODE not conserved with total rate of change {result}"
    )
