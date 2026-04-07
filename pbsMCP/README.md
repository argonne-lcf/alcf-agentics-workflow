# Agent Driven Simulation on Aurora with [PBS MCP Server](https://github.com/jtchilders/pbs-mcp-demo)

## ⚡ Quick Start

### Prerequisites

- Python 3.X with virtual environment
- Access to ALCF systems: Aurora
- Globus authentication set up for using the Inference Service

### Set up on Aurora

```bash
ssh <user_name>@aurora.alcf.anl.gov                                 # login to Aurora
git clone https://github.com/argonne-lcf/alcf-agentics-workflow.git # checkout repo
cd alcf-agentics-workflow/pbsMCP                                    # enter example dir
module load frameworks                                              # get python in PATH
python -m venv venv                                                 # setup virtual environment
source venv/bin/activate                                            # activate virtual environment
pip install -r requirements.txt                                     # install dependencies for the demo

# Install the PBS MCP Server
git clone https://github.com/jtchilders/pbs-mcp-demo.git
cd pbs-mcp-demo
git submodule update --init --recursive
pip install -r requirements.txt                                     
pip install -e .
cd ..

# Export variables needed by PBS MCP Server
export PBS_SERVER="aurora-pbs-0001.host"
export PBS_ACCOUNT="your_project"                                   # add your project name here
export PBS_ROLE="user"
source pbs-mcp-demo/external/pbs-python-api/setup_env.sh aurora

# Authenticate with Globus to access the Inference Service
python src/tools/globus_interface.py authenticate
```

Aurora is now ready to run the workflow.

> [!WARNING]
> - When authenticating with Globus, you should expect to see a link asking for authentication with ALCF credentials. If this is not happening, it may be necessary to force re-authentication by running `python src/tools/globus_interface.py authenticate --force`

### Run Demo

## 📁 Project Structure

## 🔧 Components

## 📊 Simulation Details

The demo runs simplified molecular dynamics simulations with the following characteristics:

- **Duration**: 10,000 steps × 2 fs = 20 ps simulation time
- **Output**: Stability metrics, RMSD, potential energy
- **Hardware**: Intel GPU acceleration on Aurora
- **Timeout**: 3 minutes maximum (configurable)

> **Note**: This is a demonstration with simplified physics. Real MD simulations would require proper force fields, solvation, and longer timescales.


## 🧪 Testing

## 📦 Dependencies

The project uses the following key dependencies (see `requirements.txt`):

**Core Workflow:**
- `openai>=1.0.0` - LLM API client (talks to Sophia vLLM endpoint)
- `mcp>=1.0.0` - Model Context Protocol client (talks to PBS MCP server)

**Auth:**
- `globus-sdk>=3.0.0` - Globus authentication for Sophia inference API

**Simulation (from frameworks module):**
- `openmm` - Molecular dynamics simulation engine (provided by `module load frameworks`)

## 🚀 System Requirements

- **Python**: 3.10+ (matching Aurora Framework module)
- **Platforms**: Crux (orchestration), Aurora (compute), Sophia (inference)
- **GPU**: Intel GPUs on Aurora (OpenCL support)
- **Network**: Proxy configuration for external API access from Crux

---

**Questions?** Contact the ALCF user support team or check the [ALCF documentation](https://docs.alcf.anl.gov/).