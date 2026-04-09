"""
Utility functions for the PBS MCP agentic workflow.
"""
import json
from typing import Any, Dict, List, Optional
from mcp.types import Tool as MCPTool

MAX_TOOL_RESULT_CHARS = 4000
COMPLETED_JOB_STATES = {"F", "C", "E"} # Finished, Cancelled, Error PBS job states


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