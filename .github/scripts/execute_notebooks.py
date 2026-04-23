#!/usr/bin/env python3
"""Execute jupyter notebooks.

Convert .py sources (using jupytext), copy ipynb to executed_notebooks,
execute them in place, and continue on errors.

Exit code: 0 even if some notebooks failed (mirrors previous behavior).
Writes a simple log file executed_notebooks/run.log.
"""

import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path.cwd()
SRC_DIR = ROOT / "notebooks"
OUT_DIR = ROOT / "executed_notebooks"
LOG = OUT_DIR / "run.log"

OUT_DIR.mkdir(exist_ok=True, parents=True)

with LOG.open("a") as lf:
    lf.write(f"Run started: {datetime.utcnow().isoformat()}Z\n")

# Convert .py to .ipynb if any
py_files = [str(path) for path in SRC_DIR.glob("*.py")]
if py_files:
    try:
        subprocess.check_call(
            [sys.executable, "-m", "jupytext", "--to", "ipynb"] + py_files
        )
        with LOG.open("a") as lf:
            lf.write(f"Converted {len(py_files)} .py -> .ipynb\n")
    except Exception as e:
        with LOG.open("a") as lf:
            lf.write(f"jupytext conversion failed: {e}\n")

# Allow override of execution timeout via EXEC_TIMEOUT env var (seconds)
EXEC_TIMEOUT = os.environ.get("EXEC_TIMEOUT")
if EXEC_TIMEOUT:
    try:
        timeout_val = int(EXEC_TIMEOUT)
    except Exception:
        timeout_val = 7200
else:
    timeout_val = 7200

with LOG.open("a") as lf:
    lf.write(f"Using EXEC_TIMEOUT={timeout_val}\n")

# Copy ipynb sources into executed_notebooks
ipynbs = [str(path) for path in SRC_DIR.glob("*.ipynb")]
if ipynbs:
    for p in ipynbs:
        try:
            shutil.copy(p, OUT_DIR)
        except Exception as e:
            with LOG.open("a") as lf:
                lf.write(f"Failed to copy {p}: {e}\n")
else:
    with LOG.open("a") as lf:
        lf.write("No ipynb sources found in notebooks/\n")

# Execute notebooks in executed_notebooks
exec_ipynbs = [str(path) for path in OUT_DIR.glob("*.ipynb")]
for nb in exec_ipynbs:
    with LOG.open("a") as lf:
        lf.write(f"Executing: {nb}\n")
    try:
        # Use nbconvert as in previous workflow for parity
        subprocess.check_call(
            [
                sys.executable,
                "-m",
                "jupyter",
                "nbconvert",
                "--to",
                "notebook",
                "--inplace",
                f"--ExecutePreprocessor.timeout={timeout_val}",
                "--execute",
                nb,
            ]
        )
        with LOG.open("a") as lf:
            lf.write(f"Executed successfully: {nb}\n")
    except subprocess.CalledProcessError as e:
        with LOG.open("a") as lf:
            lf.write(f"Notebook {nb} failed: returncode={e.returncode}\n")
    except Exception as e:
        with LOG.open("a") as lf:
            lf.write(f"Unexpected error executing {nb}: {e}\n")

with LOG.open("a") as lf:
    lf.write(f"Run finished: {datetime.utcnow().isoformat()}Z\n")

# Print log to stdout for GitHub Actions logs
with LOG.open() as lf:
    sys.stdout.write(lf.read())

# Always exit 0 to mirror previous behavior (continue on notebook failures)
sys.exit(0)
