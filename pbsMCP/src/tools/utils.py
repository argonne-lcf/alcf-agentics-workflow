"""
Custom tools for the PBS MCP agentic workflow.

- generate_sim_script: creates a PBS-ready bash script that runs an OpenMM
  molecular dynamics simulation via sim_kernel.py.
- read_local_file: lets the LLM read job stdout/stderr from the filesystem.

Each tool has a plain function and an OpenAI-compatible schema so the raw
OpenAI SDK can expose them alongside MCP tools.
"""
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from mcp.types import Tool as MCPTool

SRC_DIR = Path(__file__).resolve().parent.parent   # pbsMCP/src/
EXAMPLE_DIR = SRC_DIR.parent                       # pbsMCP/
VENV_DIR = EXAMPLE_DIR / "venv"
MAX_FILE_BYTES = 16_384
MAX_TOOL_RESULT_CHARS = 4000
COMPLETED_JOB_STATES = {"F", "C", "E"} # Finished, Cancelled, Error PBS job states


# ---------------------------------------------------------------------------
# generate_sim_script
# ---------------------------------------------------------------------------

GENERATE_SIM_SCRIPT_SCHEMA: Dict[str, Any] = {
    "type": "function",
    "function": {
        "name": "generate_sim_script",
        "description": (
            "Generate a PBS-ready bash script that runs an OpenMM molecular "
            "dynamics simulation via sim_kernel.py. Returns the absolute path "
            "to the generated script. Pass this path to submit_job."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "protein": {
                    "type": "string",
                    "description": "Protein name (e.g. p53, insulin, lysozyme)",
                },
                "timestep": {
                    "type": "number",
                    "description": "Integration timestep in picoseconds (default 0.002)",
                },
                "temperature": {
                    "type": "number",
                    "description": "Temperature in Kelvin (default 300)",
                },
                "steps": {
                    "type": "integer",
                    "description": "Number of MD steps (default 10000)",
                },
            },
        },
    },
}


def generate_sim_script(arguments: str) -> Dict[str, Any]:
    args = json.loads(arguments or "{}")
    protein = args.get("protein", "p53")
    timestep = args.get("timestep", 0.002)
    temperature = args.get("temperature", 300.0)
    steps = args.get("steps", 10000)

    params = {
        "protein": protein,
        "timestep": timestep,
        "temperature": temperature,
        "steps": steps,
    }

    params_file = EXAMPLE_DIR / "params.json"
    params_file.write_text(json.dumps(params, indent=2))

    venv_line = ""
    if VENV_DIR.exists():
        venv_line = f"source {VENV_DIR}/bin/activate"

    script_content = f"""#!/bin/bash
module load frameworks
cd {EXAMPLE_DIR}
{venv_line}

python src/sim_kernel.py {params_file.resolve()}
"""

    script_path = EXAMPLE_DIR / "submit_openmm.sh"
    script_path.write_text(script_content)
    script_path.chmod(0o755)

    return {
        "success": True,
        "script_path": str(script_path.resolve()),
        "params": params,
        "message": f"Generated simulation script at {script_path.resolve()}. "
                   "Pass this path to submit_job to run it.",
    }


# ---------------------------------------------------------------------------
# read_local_file
# ---------------------------------------------------------------------------

READ_LOCAL_FILE_SCHEMA: Dict[str, Any] = {
    "type": "function",
    "function": {
        "name": "read_local_file",
        "description": (
            "Read up to 16 KB from a text file on the local filesystem. "
            "Use this after a PBS job finishes to read the stdout or stderr "
            "files referenced in the job attributes (Output_Path / Error_Path)."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "path": {
                    "type": "string",
                    "description": "Absolute path to the file to read.",
                },
            },
            "required": ["path"],
        },
    },
}


def read_local_file(arguments: str) -> Dict[str, Any]:
    args = json.loads(arguments or "{}")
    path_value = args.get("path", "")
    if not path_value:
        return {"success": False, "error": "path parameter is required"}

    normalized = path_value.split(":", 1)[1] if ":" in path_value else path_value
    fpath = Path(normalized).expanduser()
    try:
        data = fpath.read_text(errors="replace")
    except OSError as exc:
        return {"success": False, "path": str(fpath), "error": str(exc)}

    snippet = data[:MAX_FILE_BYTES]
    return {
        "success": True,
        "path": str(fpath),
        "bytes": len(snippet),
        "content": snippet,
    }


# Dispatch table for custom tools
CUSTOM_TOOLS = {
    "generate_sim_script": generate_sim_script,
    "read_local_file": read_local_file,
}

CUSTOM_TOOL_SCHEMAS = [
    GENERATE_SIM_SCRIPT_SCHEMA,
    READ_LOCAL_FILE_SCHEMA,
]


# ---------------------------------------------------------------------------
# mcp_tool_to_openai_schema
# ---------------------------------------------------------------------------

def mcp_tool_to_openai_schema(tool: MCPTool) -> Dict[str, Any]:
    """Convert an MCP tool definition into an OpenAI function schema."""
    return {
        "type": "function",
        "function": {
            "name": tool.name,
            "description": tool.description or "",
            "parameters": tool.inputSchema or {"type": "object"},
        },
    }


# ---------------------------------------------------------------------------
# _extract_job_id
# ---------------------------------------------------------------------------

def extract_job_id(messages: List[Dict[str, Any]]) -> Optional[str]:
    """Scan messages for a submitted job_id."""
    import re
    for msg in reversed(messages):
        content = msg.get("content", "")
        if not isinstance(content, str) or "job_id" not in content:
            continue
        try:
            data = json.loads(content)
            if "job_id" in data:
                return data["job_id"]
        except (json.JSONDecodeError, TypeError):
            pass
        m = re.search(r'job_id["\s:=]+["\']?(\d+)', content)
        if m:
            return m.group(1)
    return None


# ---------------------------------------------------------------------------
# job_finished
# ---------------------------------------------------------------------------

def job_finished(messages: List[Dict[str, Any]]) -> bool:
    """Check if any recent message shows a terminal job state."""
    for msg in reversed(messages):
        content = msg.get("content", "")
        if not isinstance(content, str):
            continue
        try:
            data = json.loads(content)
            # Handle both {"attributes": {"job_state": ...}} and
            # {"result": {"attributes": {"job_state": ...}}}
            inner = data.get("result", data)
            attrs = inner.get("attributes", inner)
            state = (attrs.get("job_state") or "").upper()
            if state in COMPLETED_JOB_STATES:
                return True
        except (json.JSONDecodeError, TypeError):
            pass
    return False


# ---------------------------------------------------------------------------
# _truncate_result
# ---------------------------------------------------------------------------

def truncate_result(payload: Any) -> str:
    """Serialize a tool result, trimming large lists to stay within context."""
    if isinstance(payload, dict):
        inner = payload.get("result", payload)
        if isinstance(inner, dict):
            for key in ("queues", "jobs", "nodes"):
                items = inner.get(key)
                if isinstance(items, list) and len(items) > 15:
                    inner[key] = items[:15]
                    inner[f"{key}_note"] = (
                        f"Showing first 15 of {len(items)} items (truncated)."
                    )

    result_str = json.dumps(payload, indent=2)
    if len(result_str) > MAX_TOOL_RESULT_CHARS:
        result_str = result_str[:MAX_TOOL_RESULT_CHARS] + "\n... (truncated)"
    return result_str