# Auto-generated from model.json - DO NOT EDIT
# Regenerate with: uv run model-types <path-to-model.json> -o src/model/types.py

from pathlib import Path
from typing import TypedDict, Literal

class InstanceInputs(TypedDict, total=False):
    # Query chemical formula, e.g. "NaCl" or "LiFePO4"
    composition: str
    # Number of nearest neighbours to return (1-1000)
    k: int  # default: 100
    # Only return compositions containing all of these elements, e.g. ["Li", "O"]
    must_contain: list[str]  # default: []
    # Only return compositions containing none of these elements
    must_exclude: list[str]  # default: []
    # Only return compositions found in an experimental database (needs an index built with database lookup tables; ignored otherwise)
    experimental_only: bool  # default: False
    # Only return compositions with an associated crystal structure (needs an index built with database lookup tables; ignored otherwise)
    structures_only: bool  # default: False

class DatasetInputs(TypedDict, total=False):
    # CSV file with a 'composition' column of chemical formulas
    file: tuple[Path, str]
    # Number of nearest neighbours per query (1-1000)
    k: int  # default: 10
    # Only return compositions containing all of these elements
    must_contain: list[str]  # default: []
    # Only return compositions containing none of these elements
    must_exclude: list[str]  # default: []
    # Only return compositions found in an experimental database (needs an index built with database lookup tables)
    experimental_only: bool  # default: False
    # Only return compositions with an associated crystal structure (needs an index built with database lookup tables)
    structures_only: bool  # default: False