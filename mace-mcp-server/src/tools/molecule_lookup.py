"""
Molecule name to SMILES lookup via PubChem.

Wraps PubChemPy to resolve a human-readable molecule name (e.g. "aspirin",
"ethanol") into its canonical SMILES representation.
"""

from __future__ import annotations

import logging
from typing import Any, Dict

logger = logging.getLogger("mace_mcp.molecule_lookup")


def molecule_name_to_smiles(name: str) -> Dict[str, Any]:
    """Convert a molecule name to its canonical SMILES string via PubChem.

    Parameters
    ----------
    name : str
        The common or IUPAC name of the molecule (e.g. "aspirin", "water",
        "ethanol").

    Returns
    -------
    dict
        On success::

            {
                "status": "completed",
                "molecule_name": "aspirin",
                "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"
            }

        On failure::

            {"status": "failed", "error": "..."}
    """
    if not name or not name.strip():
        return {"status": "failed", "error": "Molecule name must not be empty."}

    name = name.strip()

    try:
        import pubchempy
    except ImportError:
        return {
            "status": "failed",
            "error": (
                "pubchempy is not installed. Install it with: pip install pubchempy"
            ),
        }

    try:
        compounds = pubchempy.get_compounds(str(name), "name")
    except Exception as e:
        return {
            "status": "failed",
            "error": f"PubChem lookup failed for '{name}': {e}",
        }

    if not compounds:
        return {
            "status": "failed",
            "error": (
                f"No compounds found in PubChem for '{name}'. "
                "Check the spelling or try a different name."
            ),
        }

    smiles = compounds[0].canonical_smiles
    logger.info("Resolved '%s' -> SMILES: %s", name, smiles)

    return {
        "status": "completed",
        "molecule_name": name,
        "smiles": smiles,
    }
