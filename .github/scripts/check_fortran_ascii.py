#!/usr/bin/env python3
"""Check Fortran source files for non-ASCII characters.

This script scans Fortran files (.f90, .f, .f95, .f03, .f08, .f15, .f18) for
non-ASCII characters and reports any violations. It can be used as a pre-commit hook.

Usage:
    python check_fortran_ascii.py file1.f90 file2.f90 ...
    python check_fortran_ascii.py --scan-directory src/fortran_src/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def check_file_for_non_ascii(filepath: str | Path) -> list[tuple[int, str]] | None:
    """Check a single file for non-ASCII characters.

    Parameters
    ----------
    filepath : str | Path
        Path to the file to check

    Returns
    -------
    list[tuple[int, str]] | None
        List of tuples (line_number, line_content) with non-ASCII characters
    """
    violations = []

    try:
        with Path(filepath).open("rb") as f:
            for line_num, line_bytes in enumerate(f, 1):
                try:
                    line_bytes.decode("ascii")
                except UnicodeDecodeError:
                    try:
                        line_text = line_bytes.decode("utf-8", errors="replace")
                    except Exception:
                        line_text = repr(line_bytes)
                    violations.append((line_num, line_text.rstrip("\n")))
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return None

    return violations


def main() -> int:
    """Check multiple files for non-ascii characters.

    Returns
    -------
        int: success flag (0 if no violations, 1 if violations found)
    """
    parser = argparse.ArgumentParser(
        description="Check Fortran source files for non-ASCII characters"
    )
    parser.add_argument("files", nargs="*", help="Fortran source files to check")
    parser.add_argument(
        "--scan-directory", help="Scan all Fortran files in a directory recursively"
    )

    args = parser.parse_args()

    files_to_check = []

    # Gather files from command line arguments
    if args.files:
        files_to_check.extend(args.files)

    # Gather files from directory scan
    if args.scan_directory:
        scan_path = Path(args.scan_directory)
        fortran_extensions = {".f90", ".f", ".f95", ".f03", ".f08", ".f15", ".f18"}
        for ext in fortran_extensions:
            files_to_check.extend(str(p) for p in scan_path.glob(f"**/*{ext}"))

    if not files_to_check:
        print("No files to check", file=sys.stderr)
        return 1

    found_violations = False

    for filepath in files_to_check:
        path = Path(filepath)
        if not path.exists():
            print(f"File not found: {filepath}", file=sys.stderr)
            found_violations = True
            continue

        violations = check_file_for_non_ascii(path)

        if violations is None:
            found_violations = True
            continue

        if violations:
            found_violations = True
            print(f"\n{filepath}:")
            for line_num, line_content in violations:
                print(f"  Line {line_num}: {line_content}")

    if found_violations:
        print("\nNon-ASCII characters found in Fortran source files", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
