#!/usr/bin/env python3
"""
PBS MCP Agentic Workflow
========================

Agent that uses PBS MCP tools to submit, monitor, and report on HPC jobs.
The LLM on Sophia autonomously decides which tools to call based on the
user's natural-language prompt.

Follows the openai_mcp_bridge.py pattern from pbs-mcp-demo: raw OpenAI SDK
with a manual conversation loop, MCP tools loaded via stdio, and custom
tools for simulation script generation and file reading.

Usage:
    python main.py --prompt "List the available queues on Aurora"
    python main.py --prompt "Submit an OpenMM simulation for p53 on the debug queue"
"""
from __future__ import annotations

import argparse
import asyncio
import json
import logging
import os
import sys
import time
from typing import Any, Dict, List, Optional

from openai import OpenAI
from mcp.client.session import ClientSession
from mcp.client.stdio import StdioServerParameters, stdio_client
from mcp.types import Tool as MCPTool

from tools.globus_interface import get_access_token
from tools.utils import CUSTOM_TOOLS, CUSTOM_TOOL_SCHEMAS

# Sophia LLM API configuration
SOPHIA_BASE_URL = os.getenv(
    "OPENAI_API_BASE",
    "https://inference-api.alcf.anl.gov/resource_server/sophia/vllm/v1",
)
DEFAULT_MODEL = os.getenv("OPENAI_MODEL", "meta-llama/Llama-3.3-70B-Instruct")

SYSTEM_PROMPT = (
    "You are an HPC assistant. Use the provided PBS MCP tools to fulfill requests. "
    "Call tools whenever you need live scheduler data, and report final results clearly."
)

# Filter out bash function exports that pollute subprocess environments
_ENV_BLOCKLIST = {"PS1", "PS2", "PROMPT_COMMAND"}

logger = logging.getLogger("agentic_demo")


# ---------------------------------------------------------------------------
# CLI & logging
# ---------------------------------------------------------------------------

def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="pbs-mcp-agent",
        description="LLM-driven agent with PBS MCP tools on ALCF systems",
    )
    p.add_argument(
        "--prompt", "-p",
        default=(
            "Check what queues are available, then generate and submit an "
            "OpenMM simulation for protein p53 on the debug queue. Monitor "
            "the job until it finishes and summarize the results."
        ),
        help="Natural-language instruction for the agent",
    )
    p.add_argument(
        "--model", "-m",
        default=DEFAULT_MODEL,
        help="LLM model on Sophia",
    )
    p.add_argument(
        "--poll-interval",
        type=int,
        default=120,
        help="Seconds between polls",
    )
    p.add_argument(
        "--max-polls",
        type=int,
        default=10,
        help="Max polling rounds",
    )
    p.add_argument(
        "--log-level", "-l", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    return p.parse_args()


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=logging.WARNING,
        format="%(asctime)s | %(levelname)-8s | %(message)s",
        datefmt="%d-%m %H:%M",
    )
    logger.setLevel(getattr(logging, level))


def _clean_env() -> Dict[str, str]:
    return {
        k: v for k, v in os.environ.items()
        if not k.startswith("BASH_FUNC_") and k not in _ENV_BLOCKLIST
    }


# ---------------------------------------------------------------------------
# OpenAI / MCP helpers  (following openai_mcp_bridge.py pattern)
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


async def call_mcp_tool(
    session: ClientSession, name: str, arguments: str,
) -> Dict[str, Any]:
    """Invoke an MCP tool and return a JSON-serialisable result."""
    parsed = json.loads(arguments or "{}")
    result = await session.call_tool(name, parsed)

    if result.structuredContent and isinstance(result.structuredContent, dict):
        return result.structuredContent

    text_parts = []
    for block in result.content or []:
        if hasattr(block, "text"):
            text_parts.append(block.text)

    combined = "\n".join(text_parts) if text_parts else str(result)

    try:
        return json.loads(combined)
    except (json.JSONDecodeError, TypeError):
        return {"result": combined}


# ---------------------------------------------------------------------------
# Conversation loop
# ---------------------------------------------------------------------------

COMPLETED_JOB_STATES = {"F", "C", "E"}
MAX_TOOL_RESULT_CHARS = 4000


def _truncate_result(payload: Any) -> str:
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


MAX_NUDGES = 4


def _seems_incomplete(text: str) -> bool:
    """Heuristic: does the model's text indicate it plans to do more?"""
    lower = text.lower()
    continuation = ["next,", "i will", "let me", "i'll", "now i", "proceeding"]
    return any(s in lower for s in continuation)


def _fix_arguments(
    tool_name: str, arguments: str, batch_state: Dict[str, Any],
) -> str:
    """Patch placeholder arguments using real values from earlier calls.

    Some models batch tool calls with placeholder values instead of
    chaining results sequentially. This substitutes concrete values
    captured from earlier calls in the same batch.
    """
    if tool_name == "submit_job" and batch_state.get("script_path"):
        parsed = json.loads(arguments)
        sp = parsed.get("script_path", "")
        if not sp or not os.path.isfile(sp):
            parsed["script_path"] = batch_state["script_path"]
            logger.info(f"   ⚡ Fixed script_path → {batch_state['script_path']}")
            return json.dumps(parsed)

    if tool_name == "get_job_status" and batch_state.get("job_id"):
        parsed = json.loads(arguments)
        jid = parsed.get("job_id", "")
        if not jid or not any(c.isdigit() for c in jid):
            parsed["job_id"] = batch_state["job_id"]
            logger.info(f"   ⚡ Fixed job_id → {batch_state['job_id']}")
            return json.dumps(parsed)

    return arguments


async def run_conversation(
    client: OpenAI,
    session: ClientSession,
    messages: List[Dict[str, Any]],
    tools: List[Dict[str, Any]],
    model: str,
) -> str:
    """Call the LLM in a loop until it stops making tool calls."""
    nudges = 0

    while True:
        logger.debug(f"Calling LLM ({len(messages)} messages)...")

        response = client.chat.completions.create(
            model=model,
            messages=messages,
            tools=tools,
        )
        message = response.choices[0].message

        if message.content:
            logger.debug(f"LLM says: {message.content[:200]}")

        # Append assistant message exactly as the SDK returned it
        # (Pydantic tool_calls objects, matching openai_mcp_bridge.py)
        assistant_entry: Dict[str, Any] = {
            "role": "assistant",
            "content": message.content or "",
        }
        if message.tool_calls:
            assistant_entry["tool_calls"] = message.tool_calls
        messages.append(assistant_entry)

        if not message.tool_calls:
            text = message.content or ""
            # If the model says "Next, I will..." but didn't call a tool,
            # nudge it to actually proceed with tool calls.
            if _seems_incomplete(text) and nudges < MAX_NUDGES:
                nudges += 1
                logger.info(f"🔄 Nudging model to continue (attempt {nudges}/{MAX_NUDGES})...")
                messages.append({
                    "role": "user",
                    "content": "Please proceed. Call the tools now.",
                })
                continue
            return text

        # Execute every tool call and append results.
        # Track state across calls so we can fix placeholders in the
        # same batch (e.g. script_path from generate → submit).
        batch_state: Dict[str, Any] = {}

        for tc in message.tool_calls:
            tool_name = tc.function.name
            arguments = tc.function.arguments or "{}"

            # Patch placeholder / missing arguments using earlier results
            arguments = _fix_arguments(tool_name, arguments, batch_state)

            logger.info(f"🔧 Calling {tool_name}")
            logger.debug(f"   Arguments: {arguments[:200]}")

            try:
                if tool_name in CUSTOM_TOOLS:
                    payload = CUSTOM_TOOLS[tool_name](arguments)
                else:
                    payload = await call_mcp_tool(session, tool_name, arguments)
            except Exception as exc:
                payload = {"error": str(exc)}

            # Capture useful state for subsequent calls in this batch
            if tool_name == "generate_sim_script" and isinstance(payload, dict):
                batch_state["script_path"] = payload.get("script_path", "")
            if tool_name == "submit_job" and isinstance(payload, dict):
                jid = (payload.get("job_id")
                       or payload.get("result", {}).get("job_id", ""))
                if jid:
                    batch_state["job_id"] = jid

            result_str = _truncate_result(payload)
            logger.debug(f"   ↳ Result: {result_str[:300]}")

            messages.append({
                "role": "tool",
                "tool_call_id": tc.id,
                "content": result_str,
            })


def _extract_job_id(messages: List[Dict[str, Any]]) -> Optional[str]:
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


def _job_finished(messages: List[Dict[str, Any]]) -> bool:
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
# Main
# ---------------------------------------------------------------------------

async def run_agent(args: argparse.Namespace) -> None:
    start_time = time.time()

    # Build OpenAI client pointing at Sophia
    token = get_access_token()
    client = OpenAI(base_url=SOPHIA_BASE_URL, api_key=token)

    # Launch PBS MCP server as a subprocess
    server_params = StdioServerParameters(
        command="python",
        args=["-m", "pbs_mcp.server"],
        env=_clean_env(),
    )

    async with stdio_client(server_params) as (read_stream, write_stream):
        async with ClientSession(read_stream, write_stream) as session:
            await session.initialize()

            # Discover MCP tools + add custom tools
            mcp_tools = await session.list_tools()
            openai_tools = [mcp_tool_to_openai_schema(t) for t in mcp_tools.tools]
            openai_tools.extend(CUSTOM_TOOL_SCHEMAS)

            mcp_count = len(mcp_tools.tools)
            logger.info(
                f"🔧 Loaded {mcp_count} PBS MCP tools + "
                f"{len(CUSTOM_TOOL_SCHEMAS)} custom tools"
            )
            logger.debug(f"Tools: {[t['function']['name'] for t in openai_tools]}")

            # Start conversation
            messages: List[Dict[str, Any]] = [
                {"role": "system", "content": SYSTEM_PROMPT},
                {"role": "user", "content": args.prompt},
            ]

            logger.info("🔄 Sending prompt to agent...")
            final_text = await run_conversation(
                client, session, messages, openai_tools, args.model,
            )

            logger.info("📊 Agent response:")
            print(f"\n{'='*60}")
            print(final_text)
            print(f"{'='*60}\n")

            # Polling loop if a job was submitted
            job_id = _extract_job_id(messages)
            if job_id and not _job_finished(messages):
                logger.info(f"🔬 Job {job_id} submitted. Polling until complete...")

                for poll in range(1, args.max_polls + 1):
                    await asyncio.sleep(args.poll_interval)

                    followup = (
                        f"Check whether PBS job {job_id} has finished using "
                        "get_job_status. If the job_state is 'F', read the stdout "
                        "file using read_local_file and summarize the results."
                    )
                    messages.append({"role": "user", "content": followup})

                    poll_text = await run_conversation(
                        client, session, messages, openai_tools, args.model,
                    )
                    logger.info(f"📊 Poll {poll} response:")
                    print(f"\n{'='*60}")
                    print(poll_text)
                    print(f"{'='*60}\n")

                    if _job_finished(messages):
                        logger.info(f"✅ Job {job_id} finished after {poll} poll(s).")
                        break
                else:
                    logger.warning(
                        f"⚠️  Reached max polls ({args.max_polls}) without completion."
                    )

    elapsed = time.time() - start_time
    logger.info(f"⏱️  Completed in {elapsed:.1f}s")


def main() -> None:
    args = parse_cli()
    configure_logging(args.log_level)

    logger.info("🚀 Starting PBS MCP Agentic Workflow")
    logger.info(f"Model: {args.model}")
    logger.info(f"Prompt: {args.prompt}")

    try:
        asyncio.run(run_agent(args))
    except KeyboardInterrupt:
        logger.info("🛑 Workflow interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"❌ Workflow failed: {e}")
        raise


if __name__ == "__main__":
    main()
