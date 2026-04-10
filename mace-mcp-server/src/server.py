#!/usr/bin/env python3
"""
MACE MCP Server
===============

A Model Context Protocol (MCP) server that exposes MACE machine-learning
interatomic potential calculations as tools.  Any MCP-compatible client
(Claude Desktop, an OpenAI tool-calling agent, the MCP Inspector, etc.)
can connect and run molecular simulations without writing chemistry code.

The server wraps the three tools from the ALCF AI Science Training Series
(``tools.py``):

- ``molecule_name_to_smiles``   -- resolve a molecule name to SMILES via PubChem
- ``smiles_to_coordinate_file`` -- generate 3D coordinates from SMILES (RDKit + ASE)
- ``run_mace_calculation``      -- compute energies with the MACE-MP potential

Together they form an end-to-end pipeline:
**molecule name -> SMILES -> 3D structure -> energy / optimised geometry**.

Usage
-----
The server communicates over **stdio** and is designed to be launched as a
subprocess by an MCP client::

    python src/server.py

It can also be tested interactively with the MCP Inspector::

    npx @modelcontextprotocol/inspector python src/server.py
"""

from __future__ import annotations

import json
import logging
from typing import Optional

from mcp.server.fastmcp import FastMCP

from tools.molecule_lookup import molecule_name_to_smiles as _molecule_name_to_smiles
from tools.coordinate_gen import smiles_to_coordinate_file as _smiles_to_coordinate_file
from tools.mace_calc import run_mace_calculation as _run_mace_calculation

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("mace_mcp.server")

# ---------------------------------------------------------------------------
# MCP Server
# ---------------------------------------------------------------------------

mcp = FastMCP("MACE Chemistry Server")


# ---- Tool 1: Molecule name -> SMILES --------------------------------------


@mcp.tool()
def molecule_name_to_smiles(name: str) -> str:
    """Convert a molecule name to its canonical SMILES string via PubChem.

    Use this tool when you have a common or IUPAC molecule name (e.g.
    "aspirin", "caffeine", "water") and need the SMILES representation
    for downstream processing.

    Parameters
    ----------
    name : str
        The molecule name to look up (e.g. "aspirin", "ethanol", "benzene").
    """
    logger.info("molecule_name_to_smiles: name=%s", name)
    result = _molecule_name_to_smiles(name)
    return json.dumps(result, indent=2)


# ---- Tool 2: SMILES -> coordinate file ------------------------------------


@mcp.tool()
def smiles_to_coordinate_file(
    smiles: str,
    output_file: str = "molecule.xyz",
    random_seed: int = 2025,
) -> str:
    """Convert a SMILES string to a 3D coordinate file (XYZ format).

    Generates a 3D molecular structure from a SMILES string using RDKit's
    distance-geometry embedding and UFF force-field pre-optimisation.
    The result is written as an XYZ file via ASE.

    Use this after ``molecule_name_to_smiles`` to create a structure file
    that can be passed to ``run_mace_calculation``.

    Parameters
    ----------
    smiles : str
        SMILES string (e.g. "CCO" for ethanol, "O" for water).
    output_file : str
        Path for the output XYZ file (default "molecule.xyz").
    random_seed : int
        Random seed for reproducible 3D embedding (default 2025).
    """
    logger.info("smiles_to_coordinate_file: smiles=%s output=%s", smiles, output_file)
    result = _smiles_to_coordinate_file(
        smiles=smiles,
        output_file=output_file,
        random_seed=random_seed,
    )
    return json.dumps(result, indent=2)


# ---- Tool 3: MACE calculation ---------------------------------------------


@mcp.tool()
def run_mace_calculation(
    input_file: str,
    mace_model_name: str = "small",
    device: str = "cpu",
    optimize: bool = False,
    fmax: float = 0.05,
    max_steps: int = 200,
) -> str:
    """Run a MACE machine-learning potential energy calculation.

    Computes the potential energy of a molecular structure using the
    MACE-MP foundation model.  Optionally performs a BFGS geometry
    optimisation to find the energy minimum.

    Use this after ``smiles_to_coordinate_file`` to compute the energy
    of a generated structure.

    Parameters
    ----------
    input_file : str
        Path to a structure file readable by ASE (typically an XYZ file
        produced by ``smiles_to_coordinate_file``).
    mace_model_name : str
        MACE-MP model variant: "small", "medium", or "large" (default "small").
    device : str
        Compute device: "cpu" or "cuda" (default "cpu").
    optimize : bool
        If True, run BFGS geometry optimisation (default False).
    fmax : float
        Force convergence threshold in eV/A for optimisation (default 0.05).
    max_steps : int
        Maximum optimisation steps (default 200).
    """
    logger.info(
        "run_mace_calculation: file=%s model=%s device=%s optimize=%s",
        input_file,
        mace_model_name,
        device,
        optimize,
    )
    result = _run_mace_calculation(
        input_file=input_file,
        mace_model_name=mace_model_name,
        device=device,
        optimize=optimize,
        fmax=fmax,
        max_steps=max_steps,
    )
    return json.dumps(result, indent=2)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    logger.info("Starting MACE MCP Server v0.1.0")
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
