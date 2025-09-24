#!/usr/bin/env python
"""
Agentic Workflow Demo - Main CLI Orchestrator
Demonstrates LangGraph agent running on Crux, querying LLM on Sophia,
and launching GPU tasks on Aurora via Globus Compute.
"""

import argparse
import logging
import os
import sys
import time
from datetime import datetime
from typing import Dict, Any

from langchain_openai import ChatOpenAI
from openai import OpenAI
from langgraph.graph import StateGraph, END
from langgraph.graph.message import add_messages
from typing_extensions import TypedDict

from sim_kernel import run_md_simulation
import socket

# Sophia LLM API configuration
DEFAULT_MODEL = os.getenv("OPENAI_MODEL", "meta-llama/Llama-2-7b-chat-hf")


class AgentState(TypedDict):
   """State for the agentic workflow"""
   messages: list
   protein: str
   analysis_request: str
   simulation_params: Dict[str, Any]
   simulation_result: Dict[str, Any]
   final_report: str


def parse_cli():
   """Parse command line arguments"""
   p = argparse.ArgumentParser(
      prog="agentic-demo",
      description="LangGraph → Sophia → Aurora example"
   )
   p.add_argument("--protein", "-p", default="p53", 
                  help="Protein name to analyze (default: p53)")
   p.add_argument("--model", "-m",
                  default=os.getenv("OPENAI_MODEL", DEFAULT_MODEL),
                  help="LLM model to use")
   p.add_argument("--hostname", 
                  default=socket.gethostname(),
                  help="Current hostname")
   p.add_argument("--port", 
                  default=8000,
                  help="Current port")
   p.add_argument("--log-level", "-l", default="INFO",
                  choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                  help="Set logging level")
   p.add_argument("--max-simulation-time", "-t", type=int, default=60*60,
                  help="Maximum simulation time in seconds")
   return p.parse_args()


def config_logging(level):
   """Configure logging with specified format"""
   # Set root logger to WARNING to quiet third-party libraries
   logging.basicConfig(
      level=logging.WARNING,
      format="%(asctime)s | %(levelname)-8s | %(message)s",
      datefmt="%d-%m %H:%M"
   )
   
   # Set our application loggers to the user-specified level
   app_logger = logging.getLogger("agentic_demo")
   app_logger.setLevel(getattr(logging, level))
   
   # Also set our module loggers to the same level
   compute_logger = logging.getLogger("agentic_demo.compute")
   compute_logger.setLevel(getattr(logging, level))
   
   return app_logger


def get_logger():
   """Get the application logger"""
   return logging.getLogger("agentic_demo")

def llm_analysis_node(state: AgentState) -> AgentState:
   """Node that queries LLM for protein analysis"""
   logger = get_logger()
   logger.info("🧠 Querying LLM for protein analysis...")

   try:
      model = state.get("model", DEFAULT_MODEL)

      url=f"http://localhost:{state.get('port', 8000)}/v1"
      # Create ChatOpenAI client configured for local inference
      llm = OpenAI(
               api_key="EMPTY",
               base_url=url,
            )
      
      def wait_for_vllm(llm, timeout_seconds=1000, check_interval=10):
         start_time = time.time()
         while time.time() - start_time < timeout_seconds:
            try:
               response = llm.chat.completions.create(
                   model=model,
                   messages=[{"role": "user", "content": "hello"}],
                   temperature=0.0,
                   max_tokens=1024,
                   stream=False
               )
               return True
            except Exception as e:
               logger.info(f"vLLM not ready yet. Waiting for {timeout_seconds - (time.time() - start_time):>d} more seconds")
            time.sleep(check_interval)
         return False

      wait_for_vllm(llm)

      prompt = f"""
You are an expert in molecular dynamics simulation. Given the protein, {state['protein']}, suggest settings to run OpenMM, the molecular dynamics simulation engine, for this protein. We would like to measure behavior across the range.

Please reply with only JSON, no formatting or extra text.

The JSON should have the following fields:
1. A single timestep (range from 0.0001 to 0.01)
2. A single temperature (range from 200 to 500)
3. A single steps setting (range from 100 to 50000)

Format your response as a JSON object. For example:
``json
{{
   "timestep": 0.002,
   "temperature": 300,
   "steps": 10000
}}
```
Be sure to reply with only JSON, no formatting or extra text, no more than one setting for each field. No lists.
      """
      logger.debug(f"User prompt: {prompt}")
      
      # Query the LLM 
      response = llm.chat.completions.create(
                model=model,
                messages=[{"role": "user", "content": prompt}],
                temperature=0.0,
                max_tokens=1024,
                stream=False
      )
      
      response = response.choices[0].message.content
      logger.debug(f"LLM response: {response}")
      state["analysis_request"] = prompt
      state["messages"] = add_messages(state.get("messages", []), [response])
      
      # Extract simulation parameters (simplified for demo)
      state["simulation_params"] = {
         "timestep": 10,  # ps
         "temperature": 300,  # K
         "steps": 10000,
         "protein": state["protein"]
      }
      
      logger.info("✅ LLM analysis completed")
      return state
      
   except Exception as e:
      logger.error(f"❌ LLM analysis failed: {e}")
      state["simulation_params"] = {}
      return state


def simulation_node(state: AgentState) -> AgentState:
   """Node that runs GPU simulation via Globus Compute"""
   logger = get_logger()
   logger.info("🚀 Launching GPU simulation on Aurora...")
   
   if not state.get("simulation_params"):
      logger.error("❌ No simulation parameters available")
      state["simulation_result"] = {"error": "No simulation parameters"}
      return state
   
   try:
      result = run_md_simulation(state["simulation_params"])

      state["simulation_result"] = result
      logger.info(f"✅ Simulation completed: {result.get('status', 'unknown')}")
      return state
      
   except Exception as e:
      logger.error(f"❌ Simulation failed: {e}")
      state["simulation_result"] = {"error": str(e)}
      return state


def report_node(state: AgentState) -> AgentState:
   """Node that generates final report"""
   logger = get_logger()
   logger.info("📊 Generating final report...")
   
   sim_result = state.get("simulation_result", {})
   protein = state.get("protein", "unknown")
   
   if "error" in sim_result:
      report = f"""
## Molecular Dynamics Analysis Report - {protein}

**Status**: ❌ FAILED
**Error**: {sim_result['error']}
**Timestamp**: {datetime.now().strftime('%d-%m %H:%M')}

The simulation could not be completed. Please check your Globus Compute endpoint
and ensure Aurora is accessible.
"""
   else:
      stability_score = sim_result.get("stability_score", "N/A")
      rmsd = sim_result.get("rmsd", "N/A")
      energy = sim_result.get("final_energy", "N/A")
      
      report = f"""
## Molecular Dynamics Analysis Report - {protein}

**Status**: ✅ COMPLETED
**Timestamp**: {datetime.now().strftime('%d-%m %H:%M')}

### Simulation Parameters
- Temperature: {state['simulation_params'].get('temperature', 'N/A')} K
- Timestep: {state['simulation_params'].get('timestep', 'N/A')} ps
- Steps: {state['simulation_params'].get('steps', 'N/A')}

### Results
- Stability Score: {stability_score}
- RMSD: {rmsd} Å
- Final Energy: {energy} kJ/mol

### Analysis
The protein {protein} simulation has been completed successfully on Aurora.
"""
   
   state["final_report"] = report
   logger.info("✅ Report generated")
   return state


def create_workflow():
   """Create the LangGraph workflow"""
   workflow = StateGraph(AgentState)
   
   # Add nodes
   workflow.add_node("llm_analysis", llm_analysis_node)
   workflow.add_node("simulation", simulation_node)
   workflow.add_node("report", report_node)
   
   # Define edges
   workflow.set_entry_point("llm_analysis")
   workflow.add_edge("llm_analysis", "simulation")
   workflow.add_edge("simulation", "report")
   workflow.add_edge("report", END)
   
   return workflow.compile()


def main():
   """Main entry point"""
   args = parse_cli()
   config_logging(args.log_level)
   logger = get_logger()
   
   logger.info("🚀 Starting Agentic Workflow Demo")
   logger.info(f"Target protein: {args.protein}")
   logger.info(f"Model: {args.model}")
   
   # Create initial state
   initial_state = AgentState(
      messages=[],
      protein=args.protein,
      analysis_request="",
      simulation_params={},
      simulation_result={},
      final_report="",
      model=args.model,
      hostname=args.hostname,
      port=args.port
   )
   
   try:
      # Create and run workflow
      start_time = time.time()
      workflow = create_workflow()
      
      logger.info("🔄 Running workflow...")
      final_state = workflow.invoke(initial_state)
      
      elapsed = time.time() - start_time
      logger.info(f"⏱️  Workflow completed in {elapsed:.1f}s")
      
      # Print final report
      print("\n" + "="*60)
      print(final_state["final_report"])
      print("="*60)
      
      if "error" in final_state.get("simulation_result", {}):
         sys.exit(1)
         
   except KeyboardInterrupt:
      logger.info("🛑 Workflow interrupted by user")
      sys.exit(1)
   except Exception as e:
      logger.error(f"❌ Workflow failed: {e}")
      sys.exit(1)


if __name__ == "__main__":
   main()
