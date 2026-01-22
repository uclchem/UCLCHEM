#!/usr/bin/env python3
"""Write executed_notebooks/meta.json based on environment variables and installed package version."""

import json
import os
from datetime import datetime

REF = os.environ.get("REF") or os.environ.get("INPUT_REF") or ""
COMMIT_SHA = os.environ.get("COMMIT_SHA") or os.environ.get("GITHUB_SHA") or ""
BUILD_DATE = os.environ.get("BUILD_DATE") or ""
GITHUB_ACTOR = os.environ.get("GITHUB_ACTOR") or os.environ.get("GITHUB_ACTOR") or ""

# Discover installed package version
try:
    import importlib.metadata as _im

    PKG_VER = _im.version("uclchem")
except Exception:
    PKG_VER = "0.0.0"


meta = {
    "ref": REF,
    "uclchem_version": PKG_VER,
    "commit_sha": COMMIT_SHA,
    "built_at": BUILD_DATE,
    "built_at_ts": datetime.utcnow().isoformat() + "Z",
    "built_by": GITHUB_ACTOR,
}

os.makedirs("executed_notebooks", exist_ok=True)
with open("executed_notebooks/meta.json", "w") as fh:
    json.dump(meta, fh, indent=2)

print(json.dumps(meta, indent=2))
