#!/usr/bin/env python3
"""
Example MCP client for the MACE MCP Server.

Demonstrates the complete end-to-end workflow: **molecule name -> SMILES ->
3D coordinates -> MACE energy calculation**.

No LLM is involved -- this script calls MCP tools directly to verify
the server works end-to-end.

Workflow
--------
1. Connect to the MACE MCP server over stdio.
2. Discover available tools.
3. Demo 1: Look up "ethanol" -> generate XYZ -> MACE single-point energy.
4. Demo 2: Look up "aspirin" -> generate XYZ -> MACE geometry optimisation.
5. Demo 3: Direct SMILES input ("O" for water) -> XYZ -> MACE energy.

Usage
-----
    python src/example_client.py
"""

from __future__ import annotations

import asyncio
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict

from mcp.client.session import ClientSession
from mcp.client.stdio import StdioServerParameters, stdio_client

# Path to the server script (same directory)
SERVER_SCRIPT = str(Path(__file__).resolve().parent / "server.py")


def _pretty(data: Any, label: str) -> None:
    """Print a labelled JSON result."""
    print(f"\n{'=' * 70}")
    print(f"  {label}")
    print(f"{'=' * 70}")
    if isinstance(data, str):
        try:
            data = json.loads(data)
        except (json.JSONDecodeError, TypeError):
            pass
    print(json.dumps(data, indent=2) if isinstance(data, (dict, list)) else data)


async def call_tool(
    session: ClientSession,
    name: str,
    arguments: Dict[str, Any] | None = None,
) -> Any:
    """Call an MCP tool and return the parsed result."""
    result = await session.call_tool(name, arguments or {})

    # Extract text content from the MCP result
    text_parts = []
    for block in result.content or []:
        if hasattr(block, "text"):
            text_parts.append(block.text)
    combined = "\n".join(text_parts) if text_parts else str(result)

    try:
        return json.loads(combined)
    except (json.JSONDecodeError, TypeError):
        return combined


async def run_demo() -> None:
    """Connect to the MACE MCP server and run example calculations."""
    print("MACE MCP Server - Example Client")
    print("=" * 70)
    print("End-to-end demo: molecule name -> SMILES -> 3D structure -> energy")
    print("=" * 70)

    # Launch the server as a subprocess
    server_params = StdioServerParameters(
        command=sys.executable,
        args=[SERVER_SCRIPT],
        env={**os.environ},
    )

    async with stdio_client(server_params) as (read_stream, write_stream):
        async with ClientSession(read_stream, write_stream) as session:
            await session.initialize()

            # ---- 1. Discover tools ----
            tools_result = await session.list_tools()
            tool_names = [t.name for t in tools_result.tools]
            print(f"\nDiscovered {len(tool_names)} tools:")
            for name in sorted(tool_names):
                print(f"  - {name}")

            # ================================================================
            # Demo 1: Ethanol -- full pipeline (single-point energy)
            # ================================================================

            print(f"\n{'#' * 70}")
            print("  DEMO 1: Ethanol -- Name -> SMILES -> XYZ -> Energy")
            print(f"{'#' * 70}")

            # Step 1a: Look up SMILES for ethanol
            smiles_result = await call_tool(
                session,
                "molecule_name_to_smiles",
                {"name": "ethanol"},
            )
            _pretty(smiles_result, "Step 1: molecule_name_to_smiles('ethanol')")

            if smiles_result.get("status") != "completed":
                print("ERROR: Could not resolve ethanol. Skipping demo 1.")
            else:
                ethanol_smiles = smiles_result["smiles"]

                # Step 1b: Generate XYZ file
                xyz_result = await call_tool(
                    session,
                    "smiles_to_coordinate_file",
                    {
                        "smiles": ethanol_smiles,
                        "output_file": "ethanol.xyz",
                    },
                )
                _pretty(
                    xyz_result,
                    f"Step 2: smiles_to_coordinate_file('{ethanol_smiles}')",
                )

                if xyz_result.get("status") == "completed":
                    # Step 1c: Run MACE single-point energy
                    mace_result = await call_tool(
                        session,
                        "run_mace_calculation",
                        {
                            "input_file": xyz_result["path"],
                            "mace_model_name": "small",
                            "device": "cpu",
                        },
                    )
                    _pretty(
                        mace_result,
                        "Step 3: run_mace_calculation (single-point)",
                    )

            # ================================================================
            # Demo 2: Aspirin -- full pipeline (geometry optimisation)
            # ================================================================

            print(f"\n{'#' * 70}")
            print("  DEMO 2: Aspirin -- Name -> SMILES -> XYZ -> Optimise")
            print(f"{'#' * 70}")

            smiles_result = await call_tool(
                session,
                "molecule_name_to_smiles",
                {"name": "aspirin"},
            )
            _pretty(smiles_result, "Step 1: molecule_name_to_smiles('aspirin')")

            if smiles_result.get("status") != "completed":
                print("ERROR: Could not resolve aspirin. Skipping demo 2.")
            else:
                aspirin_smiles = smiles_result["smiles"]

                xyz_result = await call_tool(
                    session,
                    "smiles_to_coordinate_file",
                    {
                        "smiles": aspirin_smiles,
                        "output_file": "aspirin.xyz",
                    },
                )
                _pretty(
                    xyz_result,
                    f"Step 2: smiles_to_coordinate_file('{aspirin_smiles}')",
                )

                if xyz_result.get("status") == "completed":
                    mace_result = await call_tool(
                        session,
                        "run_mace_calculation",
                        {
                            "input_file": xyz_result["path"],
                            "mace_model_name": "small",
                            "device": "cpu",
                            "optimize": True,
                            "fmax": 0.05,
                            "max_steps": 100,
                        },
                    )
                    _pretty(
                        mace_result,
                        "Step 3: run_mace_calculation (geometry optimisation)",
                    )

            # ================================================================
            # Demo 3: Water -- direct SMILES input (no name lookup)
            # ================================================================

            print(f"\n{'#' * 70}")
            print("  DEMO 3: Water -- Direct SMILES -> XYZ -> Energy")
            print(f"{'#' * 70}")

            xyz_result = await call_tool(
                session,
                "smiles_to_coordinate_file",
                {
                    "smiles": "O",
                    "output_file": "water.xyz",
                },
            )
            _pretty(xyz_result, "Step 1: smiles_to_coordinate_file('O')")

            if xyz_result.get("status") == "completed":
                mace_result = await call_tool(
                    session,
                    "run_mace_calculation",
                    {
                        "input_file": xyz_result["path"],
                        "mace_model_name": "small",
                        "device": "cpu",
                    },
                )
                _pretty(mace_result, "Step 2: run_mace_calculation (single-point)")

    print(f"\n{'=' * 70}")
    print("  All demos completed!")
    print(f"{'=' * 70}\n")


def main() -> None:
    asyncio.run(run_demo())


if __name__ == "__main__":
    main()
