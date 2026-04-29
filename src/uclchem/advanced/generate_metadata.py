"""Utility to regenerate the ``fortran_parameters`` section of ``fortran_metadata.yaml``.

Parses all Fortran source files in ``src/fortran_src/`` and extracts module-scope
``PARAMETER`` declarations, then writes the result back into the YAML file.  The
``internal_parameters`` and ``file_path_parameters`` sections are left untouched as
they require manual curation.

Usage::

    uclchem-generate-metadata            # update YAML in-place
    uclchem-generate-metadata --dry-run  # print diff, do not write
    uclchem-generate-metadata --check    # exit 1 if YAML would change (CI use)
"""

from __future__ import annotations

import argparse
import difflib
import re
import sys
from pathlib import Path

import yaml

from uclchem.advanced.worker_state import _MODULE_NAMES
from uclchem.utils import UCLCHEM_ROOT_DIR

_FORTRAN_SRC = UCLCHEM_ROOT_DIR.parent / "fortran_src"
_METADATA_PATH = Path(__file__).parent / "fortran_metadata.yaml"

# Fortran types that can have PARAMETER attribute
_TYPE_RE = re.compile(
    r"^\s*(?:INTEGER|REAL|LOGICAL|CHARACTER|COMPLEX|DOUBLE\s+PRECISION)"
    r"[^:]*,\s*PARAMETER\s*::\s*(.+)",
    re.IGNORECASE,
)

# Lines that increase nesting depth (we only want module-scope PARAMETERs)
_NEST_OPEN = re.compile(r"^\s*(?:SUBROUTINE|FUNCTION|CONTAINS)\b", re.IGNORECASE)
_NEST_CLOSE = re.compile(r"^\s*END\s+(?:SUBROUTINE|FUNCTION)\b", re.IGNORECASE)

# Fortran MODULE declaration
_MODULE_RE = re.compile(r"^\s*MODULE\s+(\w+)\s*$", re.IGNORECASE)


def _strip_comment(line: str) -> str:
    """Remove Fortran inline comment (everything from ``!`` onward).

    Parameters
    ----------
    line : str
        line to remove comment of

    Returns
    -------
    line : str
        Line with comments stripped, respecting character literals.
    """
    # Respect character literals by scanning manually
    in_str = False
    quote = ""
    for i, ch in enumerate(line):
        if in_str:
            if ch == quote:
                in_str = False
        elif ch in {"'", '"'}:
            in_str = True
            quote = ch
        elif ch == "!":
            return line[:i]
    return line


def _extract_param_names(rhs: str) -> list[str]:
    """Extract variable names from the RHS of a ``PARAMETER ::`` declaration.

    Handles comma-separated names with optional array dimensions and initializers::

    Parameters
    ----------
    rhs : str
        right hand side of ``PARAMETER ::`` declaration.

    Returns
    -------
    result : list[str]
        List of parameter names in lowercase.

    Examples
    --------
    >>> param_names = _extract_param_names("a = 1.0, b(10) = (/.../)")
    >>> print(param_names)
    ['a', 'b']
    """
    names: list[str] = []
    # Split on commas that are not inside parentheses
    depth = 0
    current: list[str] = []
    for ch in rhs:
        if ch == "(":
            depth += 1
            current.append(ch)
        elif ch == ")":
            depth -= 1
            current.append(ch)
        elif ch == "," and depth == 0:
            names.append("".join(current).strip())
            current = []
        else:
            current.append(ch)
    if current:
        names.append("".join(current).strip())

    result: list[str] = []
    for tok in names:
        # Take the part before '(' (array dim) or '=' (initializer)
        name = re.split(r"[=(]", tok)[0].strip()
        if re.match(r"^\w+$", name):
            result.append(name.lower())
    return result


def parse_fortran_parameters(src_dir: str | Path) -> dict[str, list[str]]:
    """Parse all ``.f90`` files in ``src_dir`` and return module-scope PARAMETERs.

    Handles Fortran continuation lines (ending with ``&`` and starting next line with ``&``).

    Parameters
    ----------
    src_dir : str | Path
        path to fortran source directory.

    Returns
    -------
    result : dict[str, list[str]]
        Mapping of f2py module name (lowercase)
        to sorted list of PARAMETER names.
    """
    src_dir = Path(src_dir)
    known_modules = set(_MODULE_NAMES)
    result: dict[str, list[str]] = {}

    for f90 in sorted(src_dir.glob("*.f90")):
        module_name: str | None = None
        params: list[str] = []
        depth = 0  # nesting level; 0 = module scope
        continuation = ""  # accumulated continuation lines

        with Path(f90).open(encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                line = _strip_comment(raw).rstrip()

                # Handle Fortran continuation: lines ending with & continue on next line
                if continuation:
                    # Previous line ended with &, prepend it
                    line = continuation + line.lstrip("&").lstrip()
                    continuation = ""

                if line.endswith("&"):
                    # This line continues on the next; accumulate and skip processing
                    continuation = line[:-1].rstrip()
                    continue

                # Detect MODULE declaration (must be depth 0, i.e. file scope)
                if module_name is None:
                    m = _MODULE_RE.match(line)
                    if m:
                        candidate = m.group(1).lower()
                        if candidate in known_modules:
                            module_name = candidate

                if module_name is None:
                    continue

                # Track nesting so we only grab module-scope PARAMETERs
                if _NEST_CLOSE.match(line):
                    depth = max(0, depth - 1)
                elif _NEST_OPEN.match(line):
                    depth += 1

                if depth > 0:
                    continue

                m = _TYPE_RE.match(line)
                if m:
                    params.extend(_extract_param_names(m.group(1)))

        if module_name and params:
            result[module_name] = sorted(set(params))

    return result


def _load_yaml(path: str | Path) -> dict:
    """Load a yaml file to a dictionary.

    Parameters
    ----------
    path : str | Path
        Path to yaml file.

    Returns
    -------
    dict
        loaded dictionary.
    """
    with Path(path).open() as f:
        return yaml.safe_load(f) or {}


def _dump_yaml(data: dict) -> str:
    """Dump yaml.

    Parameters
    ----------
    data : dict
        Data to dump.

    Returns
    -------
    str
        Dumped dictionary.
    """
    return yaml.dump(data, default_flow_style=False, sort_keys=False, allow_unicode=True)


def _merge(existing: dict, detected: dict[str, list[str]]) -> dict:
    """Merge ``detected`` into the ``fortran_parameters`` section of ``existing``.

    The ``global`` key and any other hand-maintained keys not present in
    *detected* are left untouched.  Auto-detected module keys are replaced.

    Parameters
    ----------
    existing : dict
        Existing dictionary
    detected : dict[str, list[str]]
        dictionary with key ``fortran_parameters`` to
        merge into ``existing``.

    Returns
    -------
    merged : dict
        New merged dictionary.
    """
    merged = dict(existing)
    fp: dict = dict(merged.get("fortran_parameters", {}))

    fp.update(detected)

    merged["fortran_parameters"] = fp
    return merged


def main(argv: list[str] | None = None) -> None:
    """Entry point for ``uclchem-generate-metadata``."""
    parser = argparse.ArgumentParser(
        description="Regenerate the fortran_parameters section of fortran_metadata.yaml."
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--dry-run",
        action="store_true",
        help="Print a unified diff of the changes without writing.",
    )
    mode.add_argument(
        "--check",
        action="store_true",
        help="Exit with status 1 if the YAML would change (useful in CI).",
    )
    args = parser.parse_args(argv)

    detected = parse_fortran_parameters(_FORTRAN_SRC)

    existing = _load_yaml(_METADATA_PATH)
    merged = _merge(existing, detected)

    old_text = _dump_yaml(existing)
    new_text = _dump_yaml(merged)

    if old_text == new_text:
        print("fortran_metadata.yaml is already up to date.")
        return

    if args.dry_run or args.check:
        diff = difflib.unified_diff(
            old_text.splitlines(keepends=True),
            new_text.splitlines(keepends=True),
            fromfile="fortran_metadata.yaml (current)",
            tofile="fortran_metadata.yaml (updated)",
        )
        sys.stdout.writelines(diff)
        if args.check:
            sys.exit(1)
        return

    with Path(_METADATA_PATH).open("w") as f:
        f.write(new_text)
    print(f"Updated {_METADATA_PATH}")
    for mod, names in sorted(detected.items()):
        print(f"  {mod}: {names}")


if __name__ == "__main__":
    main()
