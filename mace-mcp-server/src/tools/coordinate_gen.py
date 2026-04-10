"""
SMILES to 3D coordinate file generation.

Uses RDKit for 3D embedding and UFF force-field optimisation, then writes
the structure to an XYZ file via ASE.
"""

from __future__ import annotations

import logging
import os
from typing import Any, Dict

logger = logging.getLogger("mace_mcp.coordinate_gen")


def smiles_to_coordinate_file(
    smiles: str,
    output_file: str = "molecule.xyz",
    random_seed: int = 2025,
) -> Dict[str, Any]:
    """Convert a SMILES string to a 3D coordinate file (XYZ format).

    The function generates a 3D molecular structure from a SMILES string
    using RDKit's distance-geometry embedding and UFF force-field
    optimisation, then writes the result as an XYZ file via ASE.

    Parameters
    ----------
    smiles : str
        SMILES string representation of the molecule
        (e.g. ``"CCO"`` for ethanol).
    output_file : str, optional
        Path to save the output XYZ coordinate file.  Defaults to
        ``"molecule.xyz"`` in the current directory.
    random_seed : int, optional
        Random seed for RDKit 3D structure generation (default 2025).

    Returns
    -------
    dict
        On success::

            {
                "status": "completed",
                "artifact": "coordinate_file",
                "format": "xyz",
                "path": "/absolute/path/to/molecule.xyz",
                "smiles": "CCO",
                "n_atoms": 9
            }

        On failure::

            {"status": "failed", "error": "..."}
    """
    if not smiles or not smiles.strip():
        return {"status": "failed", "error": "SMILES string must not be empty."}

    smiles = smiles.strip()

    # --- RDKit imports ---
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return {
            "status": "failed",
            "error": (
                "RDKit is not installed. Install it with: pip install rdkit-pypi"
            ),
        }

    # --- ASE import ---
    try:
        from ase import Atoms
        from ase.io import write as ase_write
    except ImportError:
        return {
            "status": "failed",
            "error": ("ASE is not installed. Install it with: pip install ase"),
        }

    # --- Generate 3D structure ---
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "status": "failed",
            "error": f"Invalid SMILES string: '{smiles}'.",
        }

    # Add explicit hydrogens
    mol = Chem.AddHs(mol)

    # Embed in 3D
    embed_result = AllChem.EmbedMolecule(mol, randomSeed=random_seed)
    if embed_result != 0:
        return {
            "status": "failed",
            "error": (
                f"Failed to generate 3D coordinates for SMILES '{smiles}'. "
                "The molecule may be too complex or the SMILES may be invalid."
            ),
        }

    # UFF force-field optimisation
    opt_result = AllChem.UFFOptimizeMolecule(mol)
    if opt_result != 0:
        logger.warning(
            "UFF optimisation did not converge for '%s'; using best geometry obtained.",
            smiles,
        )

    # --- Extract atomic data ---
    conf = mol.GetConformer()
    numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    positions = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]

    # --- Write XYZ file via ASE ---
    atoms = Atoms(numbers=numbers, positions=positions)

    try:
        ase_write(output_file, atoms)
    except Exception as e:
        return {
            "status": "failed",
            "error": f"Failed to write coordinate file '{output_file}': {e}",
        }

    abs_path = os.path.abspath(output_file)
    logger.info(
        "Generated XYZ file: %s (%d atoms) from SMILES: %s",
        abs_path,
        len(numbers),
        smiles,
    )

    return {
        "status": "completed",
        "artifact": "coordinate_file",
        "format": "xyz",
        "path": abs_path,
        "smiles": smiles,
        "n_atoms": len(numbers),
    }
