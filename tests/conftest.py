from copy import deepcopy
from functools import cache
from pathlib import Path
from types import MappingProxyType

import pytest

import uclchem


@pytest.fixture
def uclchem_root() -> Path:
    return Path(__file__).parent.parent.resolve()


DEFAULT_PARAMETER_DICTIONARY = MappingProxyType(
    {
        "initialTemp": 10,
        "initialDens": 1e4,
        "finalTime": 1e5,
    }
)


@pytest.fixture
@cache
def basic_cloud() -> uclchem.model.Cloud:
    return uclchem.model.Cloud(param_dict=DEFAULT_PARAMETER_DICTIONARY)


@pytest.fixture
def basic_cloud_copy(basic_cloud) -> uclchem.model.Cloud:
    return deepcopy(basic_cloud)


@pytest.fixture()
def fresh_cloud() -> uclchem.model.Cloud:
    return uclchem.model.Cloud(param_dict=DEFAULT_PARAMETER_DICTIONARY)
