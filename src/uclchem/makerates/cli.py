"""Command-line interface for makerates.

Usage::

    uclchem-makerates config.yaml                    # Generate network
    uclchem-makerates config.yaml --output-dir ./out # Specify output directory
    uclchem-makerates --help                          # Show help
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from uclchem.makerates._output_resolver import (
    ProjectRootError,
    _is_valid_project_root,
    _stored_project_root,
)
from uclchem.makerates.makerates import run_makerates


def main(argv: list[str] | None = None) -> None:
    """Entry point for ``uclchem-makerates`` CLI."""
    parser = argparse.ArgumentParser(
        description="Generate a chemical reaction network from a makerates configuration.",
        prog="uclchem-makerates",
    )
    parser.add_argument(
        "config",
        type=str,
        help="Path to the YAML configuration file for makerates.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help=(
            "Output directory for generated files. "
            "If the directory contains src/uclchem and src/fortran_src subdirectories, "
            "files are written to those subdirectories (project-root layout). "
            "Otherwise all files are written to the single directory (flat layout). "
            "If not specified, uses the project root recorded at install time."
        ),
    )
    parser.add_argument(
        "--no-write",
        action="store_true",
        help="Load and validate configuration without writing output files.",
    )

    args = parser.parse_args(argv)

    # Determine the output_directory to forward to run_makerates.
    # If the user did not supply --output-dir, check the stored project root here
    # so we can emit a clear error rather than a confusing one from deep inside
    # the call stack.
    output_directory = args.output_dir
    if output_directory is None and not args.no_write:
        stored_root = _stored_project_root()
        if stored_root is None:
            print(
                "Error: No output directory was specified and no project root "
                "was recorded at install time. "
                "Pass --output-dir <path> pointing to your UCLCHEM project root.",
                file=sys.stderr,
            )
            sys.exit(1)
        if not _is_valid_project_root(stored_root):
            print(
                f"Error: The project root recorded at install time ('{stored_root}') "
                "no longer has the expected src/uclchem and src/fortran_src structure. "
                "This typically means the package was installed from PyPI (no local "
                "source tree) or the project directory was moved. "
                "Pass --output-dir <project-root> to specify where to write files.",
                file=sys.stderr,
            )
            sys.exit(1)
        # Pass the stored root as an explicit path so run_makerates resolves it
        # via Tier 1 of resolve_output_dirs (project-root layout detection).
        output_directory = str(stored_root)

    try:
        config_path = Path(args.config)
        network = run_makerates(
            configuration=config_path,
            write_files=not args.no_write,
            output_directory=output_directory,
        )
        print(
            f"Network generated successfully with {len(network.get_species_list())} species "
            f"and {len(network.get_reaction_list())} reactions"
        )
    except (ProjectRootError, Exception) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
