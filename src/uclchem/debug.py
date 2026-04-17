"""Functions to help debugging UCLCHEM."""

from pathlib import Path

from uclchemwrap import uclchemwrap as wrap


def get_f2py_signature(write: bool = False) -> str:
    """Get the signature of the UCLCHEM fortran code.

    Args:
        write (bool): Write to disk. Defaults to False.

    Returns:
        signature (str): Signature of the UCLCHEM fortran code from the f2py wrapper

    """
    signature = wrap.__doc__
    if write:
        with Path("signature_file.txt").open("w") as fh:
            fh.write(signature)
    return signature
