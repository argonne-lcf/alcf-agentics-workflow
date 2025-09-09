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

**Description**: An end-to-end workflow that demonstrates a LangGraph agent running on Crux, querying a language model on Sophia for protein analysis, then launching GPU-accelerated molecular dynamics simulations on Aurora via Globus Compute.

**Technologies**: 
- **LangGraph** - Agent orchestration framework
- **Sophia LLM Service** - Protein analysis and simulation planning
- **Aurora GPU** - High-performance molecular dynamics simulation
- **Globus Compute** - Remote job execution

**Use Cases**:
- Scientific workflow automation
- Multi-system resource coordination  
- AI-guided parameter selection
- Remote GPU job submission

---

## 🚀 Getting Started

### Prerequisites

- Python 3.8+ with virtual environment capabilities
- Access to ALCF systems (Crux, Aurora)
- Globus authentication setup
- Basic familiarity with Python and command-line tools

### Quick Setup

1. **Clone the repository**
   ```bash
   git clone <repository-url>
   cd alcf_agentics_workflow
   ```

2. **Set up virtual environment**
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install dependencies**
   ```bash
   pip install -r requirements.txt  # Install shared dependencies
   ```

4. **Choose an example workflow**
   ```bash
   cd remoteGlobusToAurora  # Start with the first example
   ```

5. **Follow the example-specific setup instructions**
   Each workflow directory contains its own detailed README with setup and usage instructions.


## 📖 Additional Resources

- **[ALCF Documentation](https://docs.alcf.anl.gov/)** - Official ALCF user guides
- **[Aurora User Guide](https://docs.alcf.anl.gov/aurora/)** - Aurora-specific documentation
- **[Globus Compute Documentation](https://globus-compute.readthedocs.io/)** - Remote execution guide
- **[LangGraph Tutorials](https://langchain-ai.github.io/langgraph/)** - Agent workflow patterns

## 📧 Support

For questions about these examples or ALCF systems:

- **General ALCF Support**: [ALCF Help Desk](https://www.alcf.anl.gov/support-center)
- **Repository Issues**: Use GitHub Issues for bug reports and feature requests
- **Example-specific questions**: Contact: jchilders@anl.gov

---

**Questions about agentic workflows on ALCF systems? We're here to help!**
