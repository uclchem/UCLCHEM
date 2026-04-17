import logging
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

import uclchem


@pytest.fixture(scope="module")
def common_output_directory(request):
    """Create temporary directory for test outputs."""
    temp_dir = Path(tempfile.mkdtemp())
    yield temp_dir
    shutil.rmtree(temp_dir, ignore_errors=True)


def test_configure_suppress_logging(caplog):
    uclchem.utils.configure_logging(level="DEBUG", stream=None)
    uclchem_logger = logging.getLogger("uclchem")

    logging.getLogger().handlers = (
        uclchem_logger.handlers
    )  # Set root logger same handler as UCLCHEM logger
    uclchem_logger.propagate = True

    with caplog.at_level("DEBUG", logger="uclchem"):
        uclchem_logger.critical("Testing")
    assert not caplog.text


def test_configure_stdout_logging(caplog):
    uclchem.utils.configure_logging(level="DEBUG", stream=sys.stdout)
    uclchem_logger = logging.getLogger("uclchem")
    uclchem_logger.propagate = True  # Have to set to True for caplog handler to find it

    with caplog.at_level("DEBUG", logger="uclchem"):
        uclchem_logger.critical("Testing")
    assert caplog.text, "Expected output, but did not get anything in caplog"


def test_configure_file_logging(caplog):
    with (
        caplog.at_level("DEBUG", logger="uclchem"),
        tempfile.NamedTemporaryFile(mode="w+") as file,
    ):
        uclchem.utils.configure_logging(level="DEBUG", stream=file.name)
        uclchem_logger = logging.getLogger("uclchem")
        uclchem_logger.propagate = (
            True  # Have to set to True for caplog handler to find it
        )
        logging.getLogger().handlers = (
            uclchem_logger.handlers
        )  # Set root logger same handler as UCLCHEM logger

        uclchem_logger.critical("Testing")
        lines = file.readlines()
    assert caplog.text, "Expected output, but did not get anything in caplog"
    assert lines


def test_configure_logging_levels():
    level_names_mapping = logging.getLevelNamesMapping()
    level_names_mapping.pop("NOTSET")
    level_names_mapping.pop("WARN")

    for message_level in level_names_mapping:
        for logger_level in level_names_mapping:
            with tempfile.NamedTemporaryFile(mode="w+") as file:
                uclchem.utils.configure_logging(level=logger_level, stream=file.name)
                uclchem_logger = logging.getLogger("uclchem")
                log_method = getattr(uclchem_logger, message_level.lower())
                uclchem_logger.propagate = True

                logging.getLogger().handlers = uclchem_logger.handlers

                log_method("Testing")

                lines = file.readlines()
            assert bool(lines) == bool(
                level_names_mapping[message_level] >= level_names_mapping[logger_level]
            ), (
                f"Logging mismatch between message level {message_level} and logger level {logger_level}"
            )


def test_model_logging(common_output_directory):
    output_file = common_output_directory / "test_model_logging_debug.dat"
    with tempfile.NamedTemporaryFile(mode="w+") as file:
        uclchem.utils.configure_logging(level="DEBUG", stream=file.name)
        model = uclchem.model.Cloud(param_dict={"outputFile": output_file})
        model.check_error(only_error=True)

        lines = file.readlines()
    assert lines

    output_file = common_output_directory / "test_model_logging_warn.dat"
    with tempfile.NamedTemporaryFile(mode="w+") as file:
        uclchem.utils.configure_logging(level="WARNING", stream=file.name)
        model = uclchem.model.Cloud(param_dict={"outputFile": output_file})
        model.check_error(only_error=True)

        lines = file.readlines()

        lines = file.readlines()
    assert not lines
