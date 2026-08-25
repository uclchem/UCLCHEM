"""Utilities for MakeRates."""

import re
from collections.abc import Iterable
from typing import Any

import numpy as np


def normalize_species_name(name: str) -> str:
    """Normalize a species name to a canonical form.

    Empty strings are preserved as empty strings. Other falsy values (like None)
    are converted to "NAN". Grain prefixes (#/@) are preserved as-is.
    A chemical isomer prefix — a single alphabetic character followed by a hyphen
    (e.g. 'o-', 'p-', 'a-', 'l-') — is lowercased so that input is case-insensitive.
    The chemical formula part is uppercased. All other names are simply uppercased.

    Parameters
    ----------
    name : str
        Raw species name string to normalize.

    Returns
    -------
    str
        Normalized species name

    Examples
    --------
    >>> normalize_species_name('H2O')
    'H2O'
    >>> normalize_species_name('o-H2')
    'o-H2'
    >>> # Elements are uppercased
    >>> normalize_species_name('o-h2')
    'o-H2'
    >>> # Prefix is lowercased
    >>> normalize_species_name('O-H2')
    'o-H2'
    >>> # Phase indications are kept
    >>> normalize_species_name('#o-H2')
    '#o-H2'

    >>> # Negative ion: len==2, not a prefix
    >>> normalize_species_name('C-')
    'C-'
    >>> # Same for electrons
    >>> normalize_species_name('E-')
    'E-'

    >>> # Empty string is kept
    >>> normalize_species_name('')
    ''
    >>> # Falsy non-string value is turned into 'NAN'
    >>> normalize_species_name(None)
    'NAN'

    """
    # Preserve empty strings; convert other falsy values to "NAN"
    if name == "":  # ruff: ignore[compare-to-empty-string]
        return ""
    if not name:
        return "NAN"
    grain_prefix = ""
    rest = name
    if rest[0] in {"#", "@"}:
        grain_prefix = rest[0]
        rest = rest[1:]
    # A chemical prefix is exactly one alpha char + '-' with more formula after it.
    if len(rest) > 2 and rest[1] == "-" and rest[0].isalpha():  # ruff: ignore[magic-value-comparison]
        return grain_prefix + rest[0].lower() + "-" + rest[2:].upper()
    return grain_prefix + rest.upper()


def find_number_of_consecutive_digits(string: str, start: int) -> int:
    """Determine the number of consecutive digits in a string, starting from some index `start`.

    Parameters
    ----------
    string : str
        the string
    start : int
        the starting index

    Returns
    -------
    num_digits : int
        the number of consecutive digits in the string starting from `start`.

    Examples
    --------
    >>> find_number_of_consecutive_digits("Hello123", 0)
    0
    >>> find_number_of_consecutive_digits("Hello123", 5)
    3
    >>> find_number_of_consecutive_digits("Hello123", 6)
    2
    >>> find_number_of_consecutive_digits("He1llo23", 2)
    1

    """
    num_digits = 0
    while start + num_digits < len(string) and string[start + num_digits].isdigit():
        num_digits += 1
    return num_digits


def is_number(s: Any) -> bool:
    """Try to convert input to a float, if it succeeds, return True.

    Parameters
    ----------
    s : Any
        Input object to check

    Returns
    -------
    bool
        True if a number, False if not.

    """
    try:
        float(s)
        return True
    except ValueError:
        return False


def sanitize_input_float(row: list[Any], index: int, default: Any = 0.0) -> float:
    """Sanitize the input. If the index is out of bounds of the row or the value.

    from the row cannot be turned into a float, use the `default` value.
    Otherwise, just gets the value from the row.

    Parameters
    ----------
    row : list[Any]
        list of objects
    index : int
        index within list to use
    default : Any
        default value to use. Default = 0.0.

    Returns
    -------
    float
        sanitized value.

    """
    output = default
    if len(row) > index and is_number(row[index]):
        output = float(row[index])
    return output


def get_default_coolants() -> list[dict[str, str]]:
    """Get the default coolant configuration for UCLCHEM.

    Returns
    -------
    list[dict[str, str]]
        List of coolant dictionaries with 'file' and 'name' keys.

    """
    return [
        {"file": "ly-a.dat", "name": "H"},
        {"file": "12c+_nometa.dat", "name": "C+"},
        {"file": "16o.dat", "name": "O"},
        {"file": "12c.dat", "name": "C"},
        {"file": "12co.dat", "name": "CO"},
        {"file": "p-h2.dat", "name": "p-H2"},
        {"file": "o-h2.dat", "name": "o-H2"},
    ]


def check_reaction(reaction_row: list[Any], keep_list: list[str]) -> bool:
    """Check a row parsed from a reaction file and checks it only contains acceptable things.

    It checks if all species in the reaction are present,
    and adds the temperature range is none is specified.

    Parameters
    ----------
    reaction_row : list[Any]
        List parsed from a reaction file
        and formatted to be able to called Reaction(reaction_row)
    keep_list : list[str]
        list of species strings that are
        acceptable in the reactant or product bits of row

    Returns
    -------
    bool
        Whether the row contains acceptable entries.

    Raises
    ------
    ValueError
        If custom desorb or freeze reactions contain species not in the
        species list.

    """
    # Convert empty strings in species slots to "NAN" for placeholder slots
    # Row entries are heterogeneous (str | float); a numeric 0.0 must NOT be
    # treated as empty, so we compare to "" explicitly rather than using falsiness.
    for i in range(7):
        if reaction_row[i] == "":  # ruff: ignore[compare-to-empty-string] heterogeneous-row
            reaction_row[i] = "NAN"

    if all(normalize_species_name(x) in keep_list for x in reaction_row[0:7]):
        if reaction_row[10] == "":  # ruff: ignore[compare-to-empty-string] heterogeneous-row
            reaction_row[10] = 0.0
            reaction_row[11] = 10000.0
        if len(reaction_row) >= 13 and reaction_row[12] == "":  # ruff: ignore[magic-value-comparison, compare-to-empty-string]
            reaction_row[12] = 0.0
        if len(reaction_row) >= 14 and reaction_row[13] == "":  # ruff: ignore[magic-value-comparison, compare-to-empty-string]
            reaction_row[13] = False
        return True
    if reaction_row[1] in {"DESORB", "FREEZE"}:
        reac_error = "Desorb or freeze reaction in custom input contains species not in species list"
        reac_error += f"\nReaction was {reaction_row}"
        raise ValueError(reac_error)
    return False


FORTRAN_LINE_LENGTH = 80


def truncate_line(input_string: str, line_length: int = FORTRAN_LINE_LENGTH) -> str:
    """Take a string and add line endings at regular intervals.

    Keeps us from overshooting fortran's line limits and, frankly,
    makes for nicer odes.f90 even if human readability isn't very important.

    Parameters
    ----------
    input_string : str
        Line of code to be truncated
    line_length : int
        rough line length. Default = :data:`FORTRAN_LINE_LENGTH`.

    Returns
    -------
    result : str
        Code string with line endings at regular intervals

    """
    result = ""
    i = 0
    # we only want to split at operators to make it look nice
    splits = ["*", "+", ",", '"']
    while len(input_string[i:]) > line_length:
        j = i + line_length
        if "\n" in input_string[i:j]:
            j = input_string[i:j].index("\n") + i + 1
            result += input_string[i:j]
        else:
            while (
                input_string[j] not in splits
                or input_string[j - 1 : j + 1].lower() == "e+"
            ):
                j -= 1
            result += input_string[i:j] + "&\n    &"
        i = j
    result += input_string[i:]
    return result


def replace_value_with_name(
    string: str, value: int | float, replace_string: str, *, truncate: bool = True
) -> str:
    """Replace all instances of `value` with a string `replace_string`.

    Uses func:`array_to_string` to determine how `value` would be formatted as a string.

    Parameters
    ----------
    string : str
        string to reformat
    value : int | float
        value to replace
    replace_string : str
        string to put instead of ``value``
    truncate : bool
        Whether to truncate the line using func:`truncate_line`.
        Default = True.

    Returns
    -------
    replaced_string : str
        string with ``value`` replaced.

    Raises
    ------
    TypeError
        If ``value`` is not an instance of ``int`` or ``float``.

    Examples
    --------
    >>> replace_value_with_name("(/0,1,2/)", 2, "X")
    '(/0,1,X/)'

    >>> replace_value_with_name("(/0.0000e+00_dp,1.0000e+00_dp,2.0000e+00_dp/)", 2.0, "X")
    '(/0.0000e+00_dp,1.0000e+00_dp,X/)'

    >>> # Replaces all instances of 'value'
    >>> replace_value_with_name(
    ...     "(/0.0000e+00_dp,1.0000e+00_dp,2.0000e+00_dp,1.0000e+00_dp/)",
    ...     1.0,
    ...     "X",
    ... )
    '(/0.0000e+00_dp,X,2.0000e+00_dp,X/)'

    """
    if string.endswith("\n"):
        suffix = "\n"
    else:
        suffix = ""
    string = "".join([line.strip() for line in string.split("\n")]).replace("&", "")
    if isinstance(value, int):
        value_string = str(value)
    elif isinstance(value, float):
        # Use array_to_string to find how 'value' would be formatted.
        array_string = array_to_string("", [value], value_type="float")
        value_string = array_string.split("/")[1]
    else:
        msg = f"replace_value_with_name is not supported for type {type(value)}. Supported types: 'float', 'int'."
        raise TypeError(msg)
    # Use regex to replace only complete tokens surrounded by array delimiters (,  (/  /)
    # or list-style delimiters ([ ] space).  A plain str.replace() can corrupt adjacent
    # values if line-continuation stripping ever places two numbers adjacent to each other.
    replaced_string = re.sub(
        r"(?<=[,\s(/\[])" + re.escape(value_string) + r"(?=[,\s)/\]])",
        replace_string,
        string,
    )
    if truncate:
        replaced_string = truncate_line(replaced_string)
    replaced_string += suffix
    return replaced_string


def array_to_string(
    name: str,
    array: list | np.ndarray,
    value_type: str = "int",
    length_name: str | None = None,
    *,
    parameter: bool = True,
) -> str:
    """Write an array to fortran source code.

    Parameters
    ----------
    name : str
        Variable name of array in Fortran
    array : list | np.ndarray
        List of values of array
    value_type : str
        The array's type. Must be one of "int", "float", "string" or "logical".
        Defaults to "int".
    length_name : str | None
        Name to put in the ``dimension`` statement.
        If None, simply but the length of the array (or shape for 2D arrays).
        Default = None.
    parameter : bool
        Whether the array is a Fortran parameter (constant).
        Defaults to True.

    Returns
    -------
    out_string : str
        String containing the Fortran code to declare this array.

    Raises
    ------
    ValueError
        Raises an error if type isn't "int", "float", "string" or "logical".
    ValueError
        If the shape of `array` is 2-dimensional, but `length_name` does not contain
        a comma.

    """
    if value_type not in {"int", "float", "string", "logical"}:
        msg = f"value_type '{value_type}' is not a valid type for array_to_string. "
        msg += "Must be one of 'int', 'float', 'string' or 'logical'."
        raise ValueError(msg)

    # Check for 2D array
    arr = np.array(array)
    if arr.ndim == 2:  # ruff: ignore[magic-value-comparison]
        if length_name is None:
            shape_name = arr.shape
        else:
            if "," not in length_name:
                msg = f"length_name '{length_name}' should contain a comma to indicate a 2D array"
                raise ValueError(msg)
            shape_name: list[str] = [i.strip() for i in length_name.split(",")]  # type: ignore[no-redef]
        shape_string: str = ",".join(str(s) for s in shape_name)  # type: ignore[no-redef]
        flat = arr.flatten(order="F")
        if value_type == "int":
            dtype = "integer"
            values = ",".join(str(int(v)) for v in flat)
        elif value_type == "float":
            dtype = "real(dp)"
            values = ",".join(f"{float(v):.4e}_dp" for v in flat)
        elif value_type == "string":
            padded_strings = pad_to_max_length([str(v) for v in flat])
            string_length = len(padded_strings[0])
            dtype = f"character(LEN={string_length})"
            values = ",".join('"' + v + '"' for v in padded_strings)
        elif value_type == "logical":
            dtype = "logical"
            values = ",".join(".true." if v else ".false." for v in flat)
        param_str = ", parameter" if parameter else ""
        out_string = f"{dtype}{param_str} :: {name}({shape_string}) = RESHAPE((/ {values} /), (/ {shape_string} /))\n"
    else:
        length_name: str = str(len(arr)) if length_name is None else length_name  # type: ignore[no-redef]
        if parameter:
            out_string = ", parameter :: " + name + f" ({length_name})=(/"
        else:
            out_string = " :: " + name + f" ({length_name})=(/"
        if value_type == "int":
            out_string = "integer" + out_string
            values = ",".join(str(value) for value in arr)
        elif value_type == "float":
            out_string = "real(dp)" + out_string
            values = ",".join(f"{value:.4e}_dp" for value in arr)
        elif value_type == "string":
            padded_strings = pad_to_max_length(arr)
            string_length = len(padded_strings[0])
            out_string = f"character(LEN={string_length:.0f})" + out_string
            values = ",".join(f'"{value}"' for value in padded_strings)
        elif value_type == "logical":
            out_string = "logical" + out_string
            values = ",".join(f"{'.true.' if value else '.false.'}" for value in arr)
        out_string += values + "/)\n"
    out_string = truncate_line(out_string)
    return out_string


def separate_common_terms(string: str, term_to_separate: str) -> str:
    """Separate out common terms in a string of an equation.

    Parameters
    ----------
    string : str
        String to clean up.
    term_to_separate : str
        Term to separate out of the parentheses.

    Returns
    -------
    string : str
        The string that evaluates to the same value, but with the common
        term taken out.

    Raises
    ------
    ValueError
        If a term `string` does not contain ``f"*{term_to_separate}"``.

    Notes
    -----
    Currently only supports multiplications, meaning that
    only strings of the form ``"A*X + B*X*Y - ..."`` are supported.
    Common terms, but with divisions are not supported.

    Examples
    --------
    >>> separate_common_terms("A*B + C*B", "B")
    '(A + C)*B'
    >>> separate_common_terms("A*B + C*B*D", "B")
    '(A + C*D)*B'

    >>> # Only takes the first instance of `term_to_separate`
    >>> separate_common_terms("A*B*B + C*B", "B")
    '(A*B + C)*B'
    >>> # Even if there are multiple that COULD be separated (potential enhancement)
    >>> separate_common_terms("A*B*B - C*B*B", "B")
    '(A*B - C*B)*B'

    >>> # Does not care about whitespace
    >>> separate_common_terms("A *    B", "B")
    '(A)*B'
    >>> # And keeps whitespace that was between the terms
    >>> separate_common_terms("A*B   +  C*B", "B")
    '(A   +  C)*B'

    """
    if "+" in string or "-" in string:
        split_string = re.split(r"(\s*[-+]\s*)", string)
        plus_and_minus = split_string[1::2]
        split_terms = split_string[::2]
    else:
        # If it is just a single term
        plus_and_minus = []
        split_terms = [string]

    mult_regexp = re.compile(r"\s*\*\s*" + re.escape(term_to_separate))
    div_regexp = re.compile(r"\s*/\s*" + re.escape(term_to_separate))

    for term_idx, term in enumerate(split_terms):
        mult_m = mult_regexp.search(term)
        div_m = div_regexp.search(term)
        if mult_m is None:
            msg = f"Multiplication by '{term_to_separate}' was not found in term '{term}'"
            raise ValueError(msg)
        if div_m is not None:
            msg = f"Division by '{term_to_separate}' found in term '{term}'"
            raise ValueError(msg)
        split_terms[term_idx] = re.sub(mult_regexp, "", term, count=1)

    string = (
        "("
        + "".join(x + y for x, y in zip(split_terms[:-1], plus_and_minus, strict=True))
        + f"{split_terms[-1]})*{term_to_separate}"
    )
    return string


def pad_to_max_length(strings: Iterable[str]) -> list[str]:
    """Pad some strings to the maximum length of a string in that iterable.

    Parameters
    ----------
    strings : Iterable[str]
        Iterable of strings to pad

    Returns
    -------
    list[str]
        List of padded strings

    """
    max_length = max(len(string) for string in strings)
    return [string.ljust(max_length) for string in strings]
