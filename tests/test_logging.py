import logging
import sys

import uclchem


def test_configure_suppress_logging(caplog):
    uclchem.utils.configure_logging(level="DEBUG", stream=None)
    uclchem_logger = logging.getLogger("uclchem")

    with caplog.at_level("DEBUG", logger="uclchem"):
        uclchem_logger.critical("Testing")

    assert not caplog.text


def test_configure_stdout_logging(caplog, capsys):
    uclchem.utils.configure_logging(level="DEBUG", stream=sys.stdout)
    uclchem_logger = logging.getLogger("uclchem")
    uclchem_logger.propagate = True  # Have to set to True for caplog handler to find it

    with caplog.at_level("DEBUG", logger="uclchem"):
        uclchem_logger.critical("Testing")

    captured = capsys.readouterr()
    assert "Testing\n" in captured.out
    assert caplog.text, "Expected output, but did not get anything in caplog"


def test_configure_file_logging(caplog, tmp_path):

    tmp_logfile = tmp_path / "uclchem.log"
    with (
        caplog.at_level("DEBUG", logger="uclchem"),
        tmp_logfile.open(mode="w+") as file,
    ):
        uclchem.utils.configure_logging(level="DEBUG", stream=file.name)
        uclchem_logger = logging.getLogger("uclchem")
        uclchem_logger.propagate = (
            True  # Have to set to True for caplog handler to find it
        )

        uclchem_logger.critical("Testing")
        lines = file.readlines()
    assert caplog.text, "Expected output, but did not get anything in caplog"
    assert lines


def test_configure_logging_levels(tmp_path):
    level_names_mapping = logging.getLevelNamesMapping()
    level_names_mapping.pop("NOTSET")
    level_names_mapping.pop("WARN")

    for message_level in level_names_mapping:
        for logger_level in level_names_mapping:
            with (tmp_path / f"{message_level}_{logger_level}.txt").open(
                mode="w+"
            ) as file:
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


def test_model_logging(tmp_path):
    output_file = tmp_path / "test_model_logging_debug.dat"
    with output_file.open(mode="w+") as file:
        uclchem.utils.configure_logging(level="DEBUG", stream=file)
        model = uclchem.model.Cloud(param_dict={"outputFile": output_file})
        model.check_error(only_error=True)

        lines = file.readlines()
    assert lines

    output_file = tmp_path / "test_model_logging_warn.dat"
    with output_file.open(mode="w+") as file:
        uclchem.utils.configure_logging(level="WARNING", stream=file.name)
        model = uclchem.model.Cloud(param_dict={"outputFile": output_file})
        model.check_error(only_error=True)

        lines = file.readlines()
    assert not lines
