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
from langgraph.graph import StateGraph, END
from langgraph.graph.message import add_messages
from typing_extensions import TypedDict

from tools.compute import GlobusComputeWrapper
from tools.globus_interface import get_access_token

# Sophia LLM API configuration
SOPHIA_BASE_URL = os.getenv("OPENAI_API_BASE", "https://data-portal-dev.cels.anl.gov/resource_server/sophia/vllm/v1")
DEFAULT_MODEL = os.getenv("OPENAI_MODEL", "meta-llama/Meta-Llama-3.1-8B-Instruct")


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
   p.add_argument("--endpoint", "-e", 
                  default=os.getenv("GC_ENDPOINT_ID"),
                  help="Globus Compute endpoint ID")
   p.add_argument("--log-level", "-l", default="INFO",
                  choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                  help="Set logging level")
   p.add_argument("--max-simulation-time", "-t", type=int, default=60*60,
                  help="Maximum simulation time in seconds")
   return p.parse_args()


def config_logging(level):
   """Configure logging with specified format"""
   logging.basicConfig(
      level=getattr(logging, level),
      format="%(asctime)s | %(levelname)-8s | %(message)s",
      datefmt="%d-%m %H:%M"
   )


def llm_analysis_node(state: AgentState) -> AgentState:
   """Node that queries LLM for protein analysis"""
   logging.info("🧠 Querying LLM for protein analysis...")
   
   try:
      # Get Globus access token for authentication
      access_token = get_access_token()
      model = state.get("model", DEFAULT_MODEL)
      
      # Create ChatOpenAI client configured for Sophia with Globus authentication
      llm = ChatOpenAI(
         model=model,
         openai_api_key=access_token,
         openai_api_base=SOPHIA_BASE_URL
      )
      
      prompt = f"""
      Analyze the protein {state['protein']} and suggest molecular dynamics simulation parameters.
      
      Please provide:
      1. Brief protein description
      2. Recommended simulation parameters (timestep, temperature, steps)
      3. Key stability metrics to monitor
      
      Format your response as a structured analysis.
      """
      
      # Query the LLM using LangChain
      response = llm.invoke(prompt)
      
      state["analysis_request"] = prompt
      state["messages"] = add_messages(state.get("messages", []), [response])
      
      # Extract simulation parameters (simplified for demo)
      state["simulation_params"] = {
         "timestep": 0.002,  # ps
         "temperature": 300,  # K
         "steps": 10000,
         "protein": state["protein"]
      }
      
      logging.info("✅ LLM analysis completed")
      logging.debug(f"LLM response: {response.content[:200]}...")
      return state
      
   except Exception as e:
      logging.error(f"❌ LLM analysis failed: {e}")
      state["simulation_params"] = {}
      return state


def simulation_node(state: AgentState) -> AgentState:
   """Node that runs GPU simulation via Globus Compute"""
   logging.info("🚀 Launching GPU simulation on Aurora...")
   
   if not state.get("simulation_params"):
      logging.error("❌ No simulation parameters available")
      state["simulation_result"] = {"error": "No simulation parameters"}
      return state
   
   try:
      gc_wrapper = GlobusComputeWrapper()
      result = gc_wrapper.submit_simulation(state["simulation_params"])
      
      state["simulation_result"] = result
      logging.info(f"✅ Simulation completed: {result.get('status', 'unknown')}")
      return state
      
   except Exception as e:
      logging.error(f"❌ Simulation failed: {e}")
      state["simulation_result"] = {"error": str(e)}
      return state


def report_node(state: AgentState) -> AgentState:
   """Node that generates final report"""
   logging.info("📊 Generating final report...")
   
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
   logging.info("✅ Report generated")
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
   
   logging.info("🚀 Starting Agentic Workflow Demo")
   logging.info(f"Target protein: {args.protein}")
   logging.info(f"Model: {args.model}")
   logging.info(f"Endpoint: {args.endpoint or 'Not set - will use default'}")
   
   if not args.endpoint:
      logging.warning("⚠️  GC_ENDPOINT_ID not set - simulation may fail")
   
   # Create initial state
   initial_state = AgentState(
      messages=[],
      protein=args.protein,
      analysis_request="",
      simulation_params={},
      simulation_result={},
      final_report="",
      model=args.model,
      endpoint=args.endpoint
   )
   
   try:
      # Create and run workflow
      start_time = time.time()
      workflow = create_workflow()
      
      logging.info("🔄 Running workflow...")
      final_state = workflow.invoke(initial_state)
      
      elapsed = time.time() - start_time
      logging.info(f"⏱️  Workflow completed in {elapsed:.1f}s")
      
      # Print final report
      print("\n" + "="*60)
      print(final_state["final_report"])
      print("="*60)
      
      if "error" in final_state.get("simulation_result", {}):
         sys.exit(1)
         
   except KeyboardInterrupt:
      logging.info("🛑 Workflow interrupted by user")
      sys.exit(1)
   except Exception as e:
      logging.error(f"❌ Workflow failed: {e}")
      sys.exit(1)


if __name__ == "__main__":
   main() 