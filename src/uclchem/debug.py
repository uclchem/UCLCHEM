"""Functions to help debugging UCLCHEM."""

from pathlib import Path

from uclchemwrap import uclchemwrap as wrap


def get_f2py_signature(*, write: bool = False) -> str:
    """Get the signature of the UCLCHEM fortran code.

    Parameters
    ----------
    write : bool
        Write to disk. Defaults to False.

    Returns
    -------
    signature : str
        Signature of the UCLCHEM fortran code from the f2py wrapper

    """
    signature = wrap.__doc__
    if write:
        Path("signature_file.txt").write_text(signature)
    return signature
