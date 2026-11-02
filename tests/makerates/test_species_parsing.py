import tempfile
from collections import Counter

import pandas as pd
import pytest

from uclchem.makerates.species import Species
from uclchem.utils import get_species_table

test_sodium_norway_problem_data = [
    pd.Series({"NAME": "NA", "MASS": 14}),  # Usually taken as NaN by pandas
    pd.Series(
        {"NAME": "NAN", "MASS": 14}
    ),  # Not taken as NaN by pandas (because we also capitalize the "a")
]


@pytest.mark.parametrize("series", test_sodium_norway_problem_data)
def test_sodium_norway_problem(series):
    series = add_dummy_values(series)
    df = pd.DataFrame([series])
    with tempfile.NamedTemporaryFile("w", dir=".") as file:
        df.to_csv(file, index=False)
        file.flush()
        read_df = get_species_table(file=file.name)
        assert all(read_df.notna())
    assert df.equals(read_df)


test_data = [
    (
        pd.Series({"NAME": "H", "MASS": 1}),
        Counter("H"),
        1,
        False,
    ),
    (pd.Series({"NAME": "#H", "MASS": 1}), Counter("H"), 1, False),
    (pd.Series({"NAME": "#H2", "MASS": 2}), Counter("HH"), 2, True),
    (pd.Series({"NAME": "CH4", "MASS": 16}), Counter({"C": 1, "H": 4}), 16, False),
    (pd.Series({"NAME": "(C)60", "MASS": 12 * 60}), Counter({"C": 60}), 12 * 60, False),
]


def add_dummy_values(series: pd.Series) -> pd.Series:
    new_columns = {
        "BINDING_ENERGY": 0,
        "SOLID_FRACTION": 0,
        "MONO_FRACTION": 0,
        "VOLCANO_FRACTION": 0,
        "ENTHALPY": 0,
        "DESORPTION_PREF": 0,
        "DIFFUSION_BARRIER": 0,
        "DIFFUSION_PREF": 0,
        "IX": 0,
        "IY": 0,
        "IZ": 0,
        "SYMMETRY_NUMBER": 0,
    }
    dct = series.to_dict()
    dct.update(new_columns)
    return pd.Series(data=dct)


@pytest.mark.parametrize(
    "series, expected_constituents, expected_mass, expected_linear", test_data
)
def test_species_class_parsing(
    series, expected_constituents, expected_mass, expected_linear
):
    spec = Species(add_dummy_values(series))

    assert series["NAME"] == spec.name

    constituents = spec.find_constituents()
    assert constituents == expected_constituents

    assert spec.mass == expected_mass
    assert spec.is_linear() == expected_linear


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
