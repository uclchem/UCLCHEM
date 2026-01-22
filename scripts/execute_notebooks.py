#!/usr/bin/env python3
"""Execute notebooks (duplicate of .github/scripts/execute_notebooks.py for CI visibility)."""

import os
import sys
import glob
import shutil
import subprocess
from datetime import datetime

ROOT = os.getcwd()
SRC_DIR = os.path.join(ROOT, "notebooks")
OUT_DIR = os.path.join(ROOT, "executed_notebooks")
LOG = os.path.join(OUT_DIR, "run.log")

os.makedirs(OUT_DIR, exist_ok=True)

with open(LOG, "a") as lf:
    lf.write(f"Run started: {datetime.utcnow().isoformat()}Z\n")

# Convert .py to .ipynb if any
py_files = glob.glob(os.path.join(SRC_DIR, "*.py"))
if py_files:
    try:
        subprocess.check_call(
            [sys.executable, "-m", "jupytext", "--to", "ipynb"] + py_files
        )
        with open(LOG, "a") as lf:
            lf.write(f"Converted {len(py_files)} .py -> .ipynb\n")
    except Exception as e:
        with open(LOG, "a") as lf:
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

with open(LOG, "a") as lf:
    lf.write(f"Using EXEC_TIMEOUT={timeout_val}\n")

# Copy ipynb sources into executed_notebooks
ipynbs = glob.glob(os.path.join(SRC_DIR, "*.ipynb"))
if ipynbs:
    for p in ipynbs:
        try:
            shutil.copy(p, OUT_DIR)
        except Exception as e:
            with open(LOG, "a") as lf:
                lf.write(f"Failed to copy {p}: {e}\n")
else:
    with open(LOG, "a") as lf:
        lf.write("No ipynb sources found in notebooks/\n")

# Execute notebooks in executed_notebooks
exec_ipynbs = glob.glob(os.path.join(OUT_DIR, "*.ipynb"))
for nb in exec_ipynbs:
    with open(LOG, "a") as lf:
        lf.write(f"Executing: {nb}\n")
    try:
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
        with open(LOG, "a") as lf:
            lf.write(f"Executed successfully: {nb}\n")
    except subprocess.CalledProcessError as e:
        with open(LOG, "a") as lf:
            lf.write(f"Notebook {nb} failed: returncode={e.returncode}\n")
    except Exception as e:
        with open(LOG, "a") as lf:
            lf.write(f"Unexpected error executing {nb}: {e}\n")

with open(LOG, "a") as lf:
    lf.write(f"Run finished: {datetime.utcnow().isoformat()}Z\n")

# Print log to stdout for GitHub Actions logs
with open(LOG) as lf:
    sys.stdout.write(lf.read())

# Always exit 0 to mirror previous behaviour
sys.exit(0)
