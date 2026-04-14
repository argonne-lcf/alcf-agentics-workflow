#!/usr/bin/env python3
"""
Environment check for the MACE MCP Server.

Verifies that all required Python dependencies are importable and prints
a summary.  Run this after installing requirements to confirm your
environment is set up correctly:

    python tests/check_env.py

Exit code 0 means all checks passed; non-zero means one or more failed.
"""

from __future__ import annotations

import sys


# ---------------------------------------------------------------------------
# ANSI colour helpers (disabled when output is not a terminal)
# ---------------------------------------------------------------------------

_USE_COLOR = hasattr(sys.stdout, "isatty") and sys.stdout.isatty()


def _green(text: str) -> str:
    return f"\033[32m{text}\033[0m" if _USE_COLOR else text


def _red(text: str) -> str:
    return f"\033[31m{text}\033[0m" if _USE_COLOR else text


def _bold(text: str) -> str:
    return f"\033[1m{text}\033[0m" if _USE_COLOR else text


# ---------------------------------------------------------------------------
# Dependency checks
# ---------------------------------------------------------------------------

# Each entry: (human-readable name, import statement to exec, install hint)
CHECKS: list[tuple[str, str, str]] = [
    (
        "MCP SDK",
        "from mcp.server.fastmcp import FastMCP",
        "pip install 'mcp>=1.0.0'",
    ),
    (
        "PubChemPy",
        "import pubchempy",
        "pip install pubchempy",
    ),
    (
        "RDKit",
        "from rdkit import Chem",
        "pip install rdkit   # or: conda install -c conda-forge rdkit",
    ),
    (
        "ASE",
        "import ase",
        "pip install 'ase>=3.22.0'",
    ),
    (
        "MACE",
        "from mace.calculators import mace_mp",
        "pip install mace-torch   # requires PyTorch",
    ),
]


def run_checks() -> bool:
    """Run all dependency checks.  Returns True if every check passed."""
    print(_bold("MACE MCP Server — environment check"))
    print("=" * 44)
    print()

    passed = 0
    failed = 0
    failures: list[tuple[str, str]] = []

    for name, import_stmt, hint in CHECKS:
        try:
            exec(import_stmt)  # noqa: S102
            print(f"  {_green('OK')}   {name}")
            passed += 1
        except Exception as exc:
            print(f"  {_red('FAIL')} {name}  ({exc})")
            failures.append((name, hint))
            failed += 1

    # Summary
    print()
    print("-" * 44)
    total = passed + failed

    if failed == 0:
        print(_green(f"All {total} checks passed.  Environment is ready."))
    else:
        print(
            _red(f"{failed} of {total} check(s) failed.")
            + "  Review the setup instructions in README.md."
        )
        print()
        print("To fix the failing dependencies:")
        for name, hint in failures:
            print(f"  {name}: {hint}")

    return failed == 0


if __name__ == "__main__":
    success = run_checks()
    sys.exit(0 if success else 1)
