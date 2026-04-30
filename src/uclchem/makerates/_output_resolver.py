"""Resolve (output_dir, fortran_src_dir) for makerates write_outputs.

Resolution priority (highest to lowest):

1. Explicit directory passed by the caller.
   - If it has src/uclchem + src/fortran_src children  → project-root layout
   - Otherwise                                          → flat layout (all files in one dir)

2. Stored project root from _project_root.py (local source install).
   - Validated at runtime; if invalid, raises ProjectRootError.

3. Legacy relative paths ../src/uclchem and ../src/fortran_src
   (caller explicitly requests this tier by passing use_legacy_relative=True).

4. None provided, no stored root, not legacy → raises ProjectRootError.

"""

from __future__ import annotations

from pathlib import Path


class ProjectRootError(RuntimeError):
    """Raised when makerates cannot determine where to write output files."""


def _stored_project_root() -> Path | None:
    """Return the project root recorded at install time, or None if unavailable.

    Returns
    -------
    Path | None
        project root, or None if it could not be found.

    """
    try:
        from uclchem._project_root import (  # noqa: PLC0415 # type: ignore[import]
            _PROJECT_ROOT,
        )

        return Path(_PROJECT_ROOT)
    except ImportError:
        return None


def _is_valid_project_root(root: Path) -> bool:
    """Determine whether ``root`` has the expected ``src/`` layout.

    Parameters
    ----------
    root : Path
        Path to project root

    Returns
    -------
    bool
        whether ``root`` has the expected project structure.

    """
    return (root / "src" / "uclchem").is_dir() and (root / "src" / "fortran_src").is_dir()


def resolve_output_dirs(
    explicit_dir: str | Path | None,
    use_legacy_relative: bool = False,
) -> tuple[Path, Path]:
    """Return ``(output_dir, fortran_src_dir)`` following the documented priority.

    Parameters
    ----------
    explicit_dir : str | Path | None
        Directory explicitly supplied by the user / CLI.
    use_legacy_relative : bool
        If True, fall through to the legacy ``../src/``
        relative-path tier instead of raising.  Pass True only from
        ``run_makerates()`` when called with no output_directory and no
        CLI involvement (programmatic / legacy Makerates/ usage). Default = False.

    Returns
    -------
    tuple[Path, Path]
        Tuple ``(output_dir, fortran_src_dir)`` both as resolved ``Path``s.

    Raises
    ------
    ProjectRootError
        When no valid output location can be determined.

    """
    # --- Tier 1: explicit directory ------------------------------------------------
    if explicit_dir is not None:
        root = Path(explicit_dir)
        if (root / "src" / "uclchem").is_dir() and (
            root / "src" / "fortran_src"
        ).is_dir():
            # Project-root layout
            return root / "src" / "uclchem", root / "src" / "fortran_src"
        else:
            # Flat layout — all files go into the single directory
            return root, root

    # --- Tier 2: stored project root -----------------------------------------------
    stored_root = _stored_project_root()
    if stored_root is not None:
        if _is_valid_project_root(stored_root):
            return stored_root / "src" / "uclchem", stored_root / "src" / "fortran_src"
        else:
            msg = (
                f"UCLCHEM was installed with project root '{stored_root}', "
                "but that directory no longer has the expected src/uclchem and "
                "src/fortran_src structure. This happens when the project was "
                "moved or when the package was installed from PyPI (where no "
                "local source tree exists). "
                "Pass --output-dir <project-root> to specify where to write files."
            )
            raise ProjectRootError(msg)

    # --- Tier 3: legacy relative paths (programmatic / Makerates/ usage) -----------
    if use_legacy_relative:
        return Path("../src/uclchem"), Path("../src/fortran_src")

    # --- Tier 4: nothing worked ----------------------------------------------------
    msg = (
        "makerates could not determine an output directory. "
        "No --output-dir was given, no stored project root is available, "
        "and legacy relative-path mode was not requested. "
        "Run uclchem-makerates with --output-dir pointing to your project root, "
        "or re-install from source."
    )
    raise ProjectRootError(msg)
