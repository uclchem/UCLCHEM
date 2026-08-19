import itertools
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

UCLCHEM_ROOT_DIR = Path(__file__).parent.parent
MAKERATES_DIR = UCLCHEM_ROOT_DIR / "Makerates"
PYTHON_ROOT_DIR = UCLCHEM_ROOT_DIR / "src" / "uclchem"
FORTRAN_ROOT_DIR = UCLCHEM_ROOT_DIR / "src" / "fortran_src"


def test_makerates_invalid_log_level_raises():
    result = subprocess.run(
        [
            sys.executable,
            "MakeRates.py",
            "--verbosity-file",
            "INVALID_LOG_LEVEL",
        ],
        capture_output=True,
        text=True,
        cwd=MAKERATES_DIR,
    )

    assert "invalid choice" in result.stderr


VERBOSITIES = ("DEBUG", "INFO", "WARNING")
verbosity_combinations = tuple(itertools.product(VERBOSITIES, VERBOSITIES))


@pytest.mark.parametrize(("verbosity_stdout", "verbosity_file"), verbosity_combinations)
def test_makerates_verbosity(verbosity_stdout, verbosity_file):
    result = subprocess.run(
        [
            sys.executable,
            "MakeRates.py",
            "--verbosity-stdout",
            verbosity_stdout,
            "--verbosity-file",
            verbosity_file,
        ],
        check=True,
        capture_output=True,
        text=True,
        cwd=MAKERATES_DIR,
    )

    def _verify_only_allowed_log_levels(lines: list[str], level: str) -> bool:
        verbosity_index = VERBOSITIES.index(level)
        verbosity_good = False
        for line in lines:
            if level in line:
                verbosity_good = True
            if "Configured logging: " in line:
                # Gives both logging levels in this line,
                # so possibly also prints the level lower than what is allowed for this
                # logger if the other logger is lower.
                continue
            not_allowed_level_found = any(
                verb in line for verb in VERBOSITIES[:verbosity_index]
            )
            if not_allowed_level_found:
                return False
        return verbosity_good

    logfile = MAKERATES_DIR / "makerates.log"
    lines = logfile.read_text().split("\n")

    assert _verify_only_allowed_log_levels(lines, verbosity_file)
    assert _verify_only_allowed_log_levels(result.stdout.split("\n"), verbosity_stdout)


def test_makerates_dry_run_doesnt_write(tmp_path: Path):
    written_files = [
        PYTHON_ROOT_DIR / "species.csv",
        PYTHON_ROOT_DIR / "reactions.csv",
        FORTRAN_ROOT_DIR / "network.f90",
        FORTRAN_ROOT_DIR / "f2py_constants.f90",
        FORTRAN_ROOT_DIR / "odes.f90",
    ]

    tmp_filepaths = [tmp_path / file.name for file in written_files]

    try:
        for filepath, tmp_filepath in zip(written_files, tmp_filepaths, strict=True):
            assert filepath.exists()
            shutil.move(filepath, tmp_filepath)

        _result = subprocess.run(
            [sys.executable, "MakeRates.py", "--dry-run"],
            check=True,
            capture_output=True,
            text=True,
            cwd=MAKERATES_DIR,
        )

        for filepath in written_files:
            assert not filepath.exists()
    finally:
        for filepath, tmp_filepath in zip(written_files, tmp_filepaths, strict=True):
            if filepath.exists():
                filepath.unlink()
            shutil.move(tmp_filepath, filepath)
