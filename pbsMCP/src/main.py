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

from tools.globus_interface import get_access_token
import utils
from tools.custom_tools import CUSTOM_TOOLS, CUSTOM_TOOL_SCHEMAS

# Sophia LLM API configuration
SOPHIA_BASE_URL = os.getenv(
    "OPENAI_API_BASE",
    "https://inference-api.alcf.anl.gov/resource_server/sophia/vllm/v1",
)
DEFAULT_MODEL = os.getenv("OPENAI_MODEL", "openai/gpt-oss-120b")
SYSTEM_PROMPT = (
    "You are an HPC assistant. Use the provided PBS MCP tools to fulfill requests. "
    "Work step by step: call one tool, use its result to decide the next action. "
    "After submitting a job, STOP and wait for the user to tell you to check the job status. "
    "Note the user can ask you to check the job status multiple times, "
    "but each time you are asked to check job status, call get_job_status ONCE, report the current state and STOP. "
    "When the job reaches state F (finished), read the stdout/stderr files and "
    "summarize the results."
)


# Configure logging
logger = logging.getLogger("pbs_mcp_demo")
def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=logging.WARNING,
        format="%(asctime)s | %(levelname)-8s | %(message)s",
        datefmt="%d-%m %H:%M",
    )
    logger.setLevel(getattr(logging, level))


# Parse command line arguments
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


async def call_mcp_tool(
    session: ClientSession, 
    name: str, 
    arguments: str,
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


# Run a conversation with the LLM
async def run_conversation(
    client: OpenAI,
    session: ClientSession,
    messages: List[Dict[str, Any]],
    tools: List[Dict[str, Any]],
    model: str,
) -> str:
    """Call the LLM in a loop until it stops making tool calls."""
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
            return message.content or ""

        # Execute every tool call and append results.
        for tc in message.tool_calls:
            tool_name = tc.function.name
            arguments = tc.function.arguments or "{}"

            logger.info(f"🔧 Calling {tool_name}")
            logger.debug(f"   Arguments: {arguments[:200]}")

            try:
                if tool_name in CUSTOM_TOOLS:
                    payload = CUSTOM_TOOLS[tool_name](arguments)
                else:
                    payload = await call_mcp_tool(session, tool_name, arguments)
            except Exception as exc:
                payload = {"error": str(exc)}

            result_str = utils.truncate_result(payload)
            logger.debug(f"   ↳ Result: {result_str[:300]}")

            messages.append({
                "role": "tool",
                "tool_call_id": tc.id,
                "content": result_str,
            })


async def run_agent(args: argparse.Namespace) -> None:
    start_time = time.time()

    # Build OpenAI client pointing at Sophia
    token = get_access_token()
    client = OpenAI(base_url=SOPHIA_BASE_URL, api_key=token)

    # Define PBS MCP server launch parameters
    server_params = StdioServerParameters(
        command="python",
        args=["-m", "pbs_mcp.server"],
        env=os.environ.copy(),
    )

    # Launch PBS MCP server as a subprocess
    async with stdio_client(server_params) as (read_stream, write_stream):
        # Create a client session to interact with the PBS MCP server
        async with ClientSession(read_stream, write_stream) as session:
            # Connect to the PBS MCP server
            await session.initialize()

            # Discover tools in the PBS MCP server and convert to OpenAI function schemas
            mcp_tools = await session.list_tools()
            logger.info(f"🔧 Found {len(mcp_tools.tools)} PBS MCP tools")
            logger.debug(f"🔧 PBS MCP Tools: {[t.name for t in mcp_tools.tools]}")
            openai_tools = [utils.mcp_tool_to_openai_schema(t) for t in mcp_tools.tools]
            
            # Add custom tools from this example (generate_sim_script and read_local_file)
            openai_tools.extend(CUSTOM_TOOL_SCHEMAS)
            logger.info(f"🔧 Added {len(CUSTOM_TOOL_SCHEMAS)} custom tools")
            logger.debug(
                f"🔧 Custom tools: {[t['function']['name'] for t in CUSTOM_TOOL_SCHEMAS]}"
            )

            # Build initial conversation messages
            messages: List[Dict[str, Any]] = [
                {"role": "system", "content": SYSTEM_PROMPT},
                {"role": "user", "content": args.prompt},
            ]

            # Send prompt to the agent
            logger.info("🔄 Sending prompt to agent...")
            final_text = await run_conversation(
                client, session, messages, openai_tools, args.model,
            )

            logger.info("🤖 Agent response:")
            print(f"\n{'='*60}")
            print(final_text)
            print(f"{'='*60}\n")

            # Polling loop if a job was submitted
            job_id = utils.extract_job_id(messages)
            if job_id and not utils.job_finished(messages):
                logger.info(
                    f"⏳ Job {job_id} submitted. Polling every "
                    f"{args.poll_interval} seconds until complete..."
                )

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
                    logger.info(f"🤖 Poll {poll} response:")
                    print(f"\n{'='*60}")
                    print(poll_text)
                    print(f"{'='*60}\n")

                    if utils.job_finished(messages):
                        logger.info(f"✅ Job {job_id} finished after {poll} poll(s).")
                        break
                else:
                    logger.warning(
                        f"⚠️  Reached max polls ({args.max_polls}) without completion."
                    )

    elapsed = time.time() - start_time
    logger.info(f"⏱️ Completed in {elapsed:.1f}s")


def main() -> None:
    args = parse_cli()
    configure_logging(args.log_level)

    logger.info("🚀 Starting PBS MCP Agentic Workflow")
    logger.info(f"🤖 Model: {args.model}")
    logger.debug(f"📋 System prompt: {SYSTEM_PROMPT}")
    logger.info(f"💬 Prompt: {args.prompt}")

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
