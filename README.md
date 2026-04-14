# ALCF Agentic Workflows

A collection of example agentic workflows demonstrating the integration of AI agents with ALCF (Argonne Leadership Computing Facility) systems. These examples showcase how to build intelligent workflows that can orchestrate complex computational tasks across multiple ALCF resources including Aurora, Sophia, and Crux.

## 🎯 Overview

This repository provides practical examples for developing agentic workflows that leverage ALCF's high-performance computing infrastructure. Each example demonstrates different patterns and use cases for integrating AI agents with scientific computing workflows.

### What are Agentic Workflows?

Agentic workflows are intelligent computational pipelines that use AI agents to make decisions, coordinate tasks, and adapt to changing conditions. These workflows can:

- **Reason about computational problems** using large language models
- **Dynamically adjust parameters** based on intermediate results  
- **Orchestrate multi-system workflows** across different computing resources
- **Handle errors and retries** intelligently
- **Generate reports and insights** from computational results

## 📂 Example Workflows

### 1. [Remote Agent Driven Simulation on Aurora](./remoteGlobusToAurora/)

**Description**: An end-to-end workflow that demonstrates a LangGraph agent running on Crux or your laptop/local system, querying a language model on Sophia for protein analysis, then launching GPU-accelerated molecular dynamics simulations on Aurora via Globus Compute.

**Technologies**: 
- **LangGraph** - Agent orchestration framework
- **Sophia LLM Service** - Protein analysis and simulation planning
- **Aurora** - High-performance molecular dynamics simulation
- **Globus Compute** - Remote job execution

**Use Cases**:
- Scientific workflow automation
- Multi-system resource coordination  
- AI-guided parameter selection
- Remote GPU job submission

### 2. [Local Agent Driven Simulation on Aurora](./localWorkflow)

**Description**: An end-to-end workflow that demonstrates an agentic workflow run entirely within the same batch job. This means the LangGraph agent, the language model, then the simulations are all run on Aurora or Polaris compute nodes. 

**Technologies**:
- **LangGraph** - Agent orchestration framework
- **vLLM** - Local LLM serving for protein analysis and simulation planning
- **Aurora** - High-performance molecular dynamics simulation and LLM inference

**Use Cases**:
- Scientific workflow automation
- AI-guided parameter selection
- Low-latency workflows 
- Use of tailored LLMs

### 3. [Leveraging PBS MCP Server for Job Query and Scheduling](./pbsMCP)

**Description**: An end-to-end workflow that runs on an Aurora login node. A language model on the Sophia inference service chooses tools over a manual OpenAI-style tool-calling loop, while a local [PBS MCP server](https://github.com/jtchilders/pbs-mcp-demo) exposes PBS operations (queues, submit, status, job listing). The agent can generate a PBS batch script, submit OpenMM molecular dynamics to Aurora compute queues, poll until completion, and read job stdout/stderr for summaries.

**Technologies**:
- **Sophia LLM Service** - Tool-calling and natural-language orchestration (OpenAI-compatible API)
- **Model Context Protocol (MCP)** - Structured access to PBS tools via stdio subprocess
- **PBS MCP Server** - Scheduler integration (submit, status, queues, etc.)
- **Aurora** - Login-node orchestration and compute-node simulation

**Use Cases**:
- Natural-language driven PBS interaction
- Natural-launguage driven job submission, monitoring and job output summary

### 4. [MACE MCP Server for Molecular Simulations](./mace-mcp-server/)

**Description**: An MCP server that exposes MACE machine-learning interatomic potential calculations as tools. Any MCP-compatible client (OpenCode, Claude Code, etc.) can connect and run molecular energy calculations -- molecule name lookup via PubChem, 3D coordinate generation via RDKit, and energy computation / geometry optimisation via the MACE-MP foundation model -- without writing chemistry code. By pairing an MCP server with OpenCode, the orchestration layer (LLM client, tool discovery, and tool calling) is provided out of the box, so development effort focuses entirely on the scientific tools themselves.

**Technologies**:
- **Model Context Protocol (MCP)** - FastMCP server exposing MACE tools via stdio
- **MACE-MP** - Machine-learning interatomic potential (foundation model for molecules and materials)
- **RDKit** - 3D molecular coordinate generation
- **PubChemPy** - Molecule name to SMILES resolution via PubChem
- **ASE** - Atomic Simulation Environment for structure I/O and optimisation
- **OpenCode** - MCP-compatible AI coding agent with ALCF inference endpoint support

**Use Cases**:
- AI-assisted molecular energy calculations
- LLM-driven molecular structure exploration
- Rapid prototyping of computational chemistry workflows
- Protocol-agnostic tool serving for scientific computing

---

## 🚀 Getting Started

### Prerequisites

- Python 3.8+ with virtual environment capabilities
- Access to ALCF systems (Crux, Aurora, Polaris)
- Access to Inference Endpoint on Sophia/Metis
- Globus authentication setup
- Basic familiarity with Python and command-line tools

### Quick Setup

1. **Clone the repository**
   ```bash
   git clone https://github.com/argonne-lcf/alcf-agentics-workflow.git
   cd alcf_agentics_workflow
   ```

2. **Choose an example workflow**
   ```bash
   cd remoteGlobusToAurora  # Start with the first example
   ```

3. **Follow the example-specific setup instructions**
   Each workflow directory contains its own detailed README with setup and usage instructions.


## 📖 Additional Resources

- **[ALCF Documentation](https://docs.alcf.anl.gov/)** - Official ALCF user guides
- **[Aurora User Guide](https://docs.alcf.anl.gov/aurora/)** - Aurora-specific documentation
- **[ALCF Inference Endpoints](https://docs.alcf.anl.gov/services/inference-endpoints/) - ALCF Inference Endpoint documentation
- **[Globus Compute Documentation](https://globus-compute.readthedocs.io/)** - Remote execution guide
- **[LangGraph Tutorials](https://langchain-ai.github.io/langgraph/)** - Agent workflow patterns

## 📧 Support

For questions about these examples or ALCF systems:

- **General ALCF Support**: [ALCF Help Desk](https://www.alcf.anl.gov/support-center)
- **Repository Issues**: Use GitHub Issues for bug reports and feature requests
- **Example-specific questions**: Contact: jchilders@anl.gov, rbalin@anl.gov

---

**Questions about agentic workflows on ALCF systems? We're here to help!**
