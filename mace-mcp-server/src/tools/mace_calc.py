"""
MACE machine-learning interatomic potential calculations.

Provides single-point energy evaluation and optional geometry optimisation
using the MACE foundation model (``mace_mp``) via ASE.
"""

from __future__ import annotations

import logging
import os
from typing import Any, Dict, Literal

logger = logging.getLogger("mace_mcp.mace_calc")


def run_mace_calculation(
    input_file: str,
    mace_model_name: str = "small",
    device: Literal["cpu", "cuda"] = "cpu",
    optimize: bool = False,
    fmax: float = 0.05,
    max_steps: int = 200,
) -> Dict[str, Any]:
    """Run a MACE energy calculation, optionally with geometry optimisation.

    Uses the MACE-MP foundation model to compute the potential energy of
    a molecular structure.  If ``optimize=True``, a BFGS geometry
    optimisation is performed before reporting the final energy.

    Parameters
    ----------
    input_file : str
        Path to a structure file readable by ASE (typically ``.xyz``).
    mace_model_name : str, optional
        Name of the MACE-MP model to load (default ``"small"``).
        Common choices: ``"small"``, ``"medium"``, ``"large"``.
    device : {"cpu", "cuda"}, optional
        Device to run the calculation on (default ``"cpu"``).
    optimize : bool, optional
        If ``True``, run a BFGS geometry optimisation (default ``False``).
    fmax : float, optional
        Maximum force convergence criterion in eV/A for optimisation
        (default 0.05).
    max_steps : int, optional
        Maximum number of BFGS optimisation steps (default 200).

    Returns
    -------
    dict
        On success (single-point)::

            {
                "status": "completed",
                "mode": "single_point",
                "input_file": "molecule.xyz",
                "mace_model": "small",
                "device": "cpu",
                "energy_eV": -123.456
            }

        On success (geometry optimisation)::

            {
                "status": "completed",
                "mode": "geometry_optimization",
                "converged": true,
                "input_file": "molecule.xyz",
                "mace_model": "small",
                "device": "cpu",
                "energy_eV": -123.789,
                "final_positions": [[x, y, z], ...],
                "final_cell": [[ax, ay, az], ...],
                "fmax": 0.05,
                "max_steps": 200
            }

        On failure::

            {"status": "failed", "error": "..."}
    """
    # --- Validate input file ---
    if not input_file or not input_file.strip():
        return {"status": "failed", "error": "Input file path must not be empty."}

    input_file = input_file.strip()

    if not os.path.isfile(input_file):
        return {
            "status": "failed",
            "error": f"Input structure file '{input_file}' does not exist.",
        }

    # --- ASE import ---
    try:
        from ase.io import read as ase_read
        from ase.optimize import BFGS
    except ImportError:
        return {
            "status": "failed",
            "error": "ASE is not installed. Install it with: pip install ase",
        }

    # --- MACE import ---
    try:
        from mace.calculators import mace_mp
    except ImportError:
        return {
            "status": "failed",
            "error": ("MACE is not installed. Install it with: pip install mace-torch"),
        }

    # --- Normalise device ---
    dev = device.lower() if device else "cpu"
    if dev not in ("cpu", "cuda"):
        dev = "cpu"

    # --- Read structure ---
    try:
        atoms = ase_read(input_file)
    except Exception as e:
        return {
            "status": "failed",
            "error": f"Could not read '{input_file}' with ASE: {e}",
        }

    # --- Create MACE calculator ---
    try:
        calc = mace_mp(model=mace_model_name, device=dev)
    except Exception as e:
        return {
            "status": "failed",
            "error": (
                f"Could not load MACE model '{mace_model_name}'. Original error: {e}"
            ),
        }

    atoms.calc = calc

    # --- Single-point calculation ---
    if not optimize:
        try:
            energy = float(atoms.get_potential_energy())
        except Exception as e:
            return {
                "status": "failed",
                "error": f"MACE single-point calculation failed: {e}",
            }

        logger.info(
            "Single-point energy: %.6f eV (file=%s, model=%s, device=%s)",
            energy,
            input_file,
            mace_model_name,
            dev,
        )

        return {
            "status": "completed",
            "message": "MACE single-point energy computed.",
            "mode": "single_point",
            "input_file": input_file,
            "mace_model": mace_model_name,
            "device": dev,
            "energy_eV": energy,
        }

    # --- Geometry optimisation ---
    try:
        opt = BFGS(atoms)
        opt.run(fmax=fmax, steps=max_steps)
        converged = True
    except Exception as e:
        return {
            "status": "failed",
            "error": f"Geometry optimisation failed: {e}",
        }

    try:
        final_energy = float(atoms.get_potential_energy())
    except Exception as e:
        return {
            "status": "failed",
            "error": f"Could not retrieve energy after optimisation: {e}",
        }

    logger.info(
        "Geometry optimisation completed: %.6f eV "
        "(converged=%s, file=%s, model=%s, device=%s)",
        final_energy,
        converged,
        input_file,
        mace_model_name,
        dev,
    )

    return {
        "status": "completed",
        "message": "MACE geometry optimisation completed.",
        "mode": "geometry_optimization",
        "converged": converged,
        "input_file": input_file,
        "mace_model": mace_model_name,
        "device": dev,
        "energy_eV": final_energy,
        "final_positions": atoms.get_positions().tolist(),
        "final_cell": atoms.get_cell().tolist(),
        "fmax": fmax,
        "max_steps": max_steps,
    }
