# Remote Agent Driven Simulation on Aurora

A demonstration of an end-to-end agentic workflow that showcases ALCF infrastructure integration. The workflow, running on Crux or your local machine, queries a language model on Sophia for protein analysis, then launches GPU-accelerated molecular dynamics simulations on Aurora via Globus Compute.
Questions: jchilders@anl.gov

## ⚡ Quick Start

### Prerequisites

- Python 3.X with virtual environment
- Access to ALCF systems: Aurora (Crux is optional)
- Globus authentication set up
- Aurora endpoint configured

### Setup Aurora with [Globus Compute Endpoint](https://globus-compute.readthedocs.io/en/latest/endpoints/endpoints.html)

```bash
ssh <user_name>@aurora.alcf.anl.gov                                 # login to Aurora
git clone https://github.com/argonne-lcf/alcf-agentics-workflow.git # checkout repo
cd alcf-agentics-workflow/remoteGlobusToAurora                      # enter example dir
module load frameworks                                              # get python in PATH
python -m venv venv                                                 # setup virtual environment
source venv/bin/activate                                            # activate virtual environment
pip install -r requirements.txt                                     # install dependencies for the demo
pip install globus-compute-endpoint==3.15.0                         # install globus compute endpoint

# next we need to generate the endpoint configuration file for Globus Compute
python scripts/gen_endpoint_config.py --repo-path $PWD --venv-path $PWD/venv -o my-endpoint-config.yaml --account <project_name> --queue <queue_name>

# create a globus compute endpoint (skip if one is already present)
# (see below for troubleshooting this step)
globus-compute-endpoint configure --endpoint-config my-endpoint-config.yaml my-aurora-endpoint

# start globus endpoint (may ask to authenticate with Globus)
globus-compute-endpoint start my-aurora-endpoint --detach

# verify endpoint is running (ensure status is Initialized, not Disconnected)
# and get the Endpoint ID (this is needed for the setup on Crux/local)
globus-compute-endpoint list
```

Now Aurora is ready to run the workflow.

> [!WARNING]
> - For Globus Compute to work properly, the same MAJOR.MINOR version of python must be used on the endpoint as the one used to create the virtual environment. Currently the Aurora Framework module loads Python 3.12.12.
> - When configuring the Globus Compute Endpoint, you can use the `--account <project_name> --queue <queue_name>` options to specify your project name and the queue you wish to submit to (the debug queue is fine for this example)
> - Only run the `globus-compute-endpoint start` command once. You can see if an endpoint was already created with `globus-compute-endpoint list`. 
> - If you want to delete existing Globus endpoint, run `globus-compute-endpoint stop <endpoint_name>` and `rm -r ~/.globus_compute/<endpoint_name>`. Now you can create a new one with `globus-compute-endpoint start <endpoint_name>`.
> - When starting the Globus Compute Enpoint, it may ask for authentication. Copy-paste the printed link to a browser and use ALCF credential to authenticate. Then paste the provided Authorization Code in the terminal.


### Setup Crux with [Globus Compute Endpoint](https://globus-compute.readthedocs.io/en/latest/endpoints/endpoints.html)


```bash
ssh <user_name>@crux.alcf.anl.gov                                   # login to Crux
git clone https://github.com/argonne-lcf/alcf-agentics-workflow.git # checkout repo
cd alcf-agentics-workflow/remoteGlobusToAurora                      # enter example dir
module use /soft/modulefiles/                                       # get python and conda
module load conda
conda activate
conda create -y -p $PWD/env python==3.12.12                         # create new conda env matching Aurora python version 
conda activate $PWD/env                                             # activate conda env
pip install -r requirements.txt                                     # install dependencies for the demo

# Set required environment variables
export GC_ENDPOINT_ID="<UUID>" # from the previous step
export OPENAI_MODEL="meta-llama/Meta-Llama-3.1-8B-Instruct"  # Optional: defaults to this

# you need the http proxy setup to get off of Crux to the Sophia Inference Service
export http_proxy="http://proxy.alcf.anl.gov:3128"
export https_proxy="http://proxy.alcf.anl.gov:3128"

# Authenticate with Globus (see below for debugging notes)
python src/tools/globus_interface.py authenticate

# Verify setup
python scripts/globus_check.py
```

> [!WARNING]
> - For Globus Compute to work properly, the same MAJOR.MINOR version of python must be used on the endpoint as the one used to create the virtual environment. Currently the Aurora Framework module loads Python 3.12.12, so the conda environment on Crux is created with Python 3.12.12.
> - When authenticating with Globus, you should expect to see a link asking for authentication with ALCF credentials. If this is not happening, it may be necessary to force re-authentication by running `python src/tools/globus_interface.py authenticate --force`
> - When running `globus_check.py` and checking that the endpoint is reachable, you should expaect to have to authenticate again. If this is not happening, try deleting all cached files with `rm -r ~/.globus_compute` and running the check again.


### Run Demo

From the login nodes on Crux, run

```bash
# Basic run with p53 protein
python src/main.py --protein p53
```

This will run the workflow and print the results to the screen.

**Expected output:**
```bash
30-07 14:13 | INFO     | 🚀 Starting Agentic Workflow Demo
30-07 14:13 | INFO     | Target protein: p53
30-07 14:13 | INFO     | Model: meta-llama/Meta-Llama-3.1-8B-Instruct
30-07 14:13 | INFO     | Endpoint: 7400de92-807e-4848-908e-a76ffb21bee9
30-07 14:13 | INFO     | 🔄 Running workflow...
30-07 14:13 | INFO     | 🧠 Querying LLM for protein analysis...
30-07 14:13 | INFO     | ✅ LLM analysis completed
30-07 14:13 | INFO     | 🚀 Launching GPU simulation on Aurora...
30-07 14:13 | INFO     | Submitting simulation job to endpoint 7400de92-807e-4848-908e-a76ffb21bee9
30-07 14:13 | INFO     | Job submitted with ID: None
30-07 14:14 | INFO     | Job completed successfully in 60.1s
30-07 14:14 | INFO     | ✅ Simulation completed: completed
30-07 14:14 | INFO     | 📊 Generating final report...
30-07 14:14 | INFO     | ✅ Report generated
30-07 14:14 | INFO     | ⏱️  Workflow completed in 77.5s

============================================================

## Molecular Dynamics Analysis Report - p53

**Status**: ✅ COMPLETED
**Timestamp**: 30-07 14:14

### Simulation Parameters
- Temperature: 300 K
- Timestep: 0.002 ps
- Steps: 10000

### Results
- Stability Score: 0.598
- RMSD: 1.804 Å
- Final Energy: -12327.72 kJ/mol

### Analysis
The protein p53 simulation has been completed successfully on Aurora.

============================================================
```

To run with a custom protein and verbose logging, execute

```bash
# Custom protein and verbose logging
python src/main.py --protein insulin --log-level DEBUG
```

The full set of options are:
```bash
# Full options
python remoteGlobusToAurora/src/main.py \
   --protein myoglobin \
   --model gpt-4 \
   --endpoint custom-endpoint-id \
   --log-level INFO \
   --max-simulation-time 120
```

## 📁 Project Structure

```
remoteGlobusToAurora/
├── src/
│   ├── main.py              # CLI orchestrator with LangGraph
│   ├── sim_kernel.py        # OpenMM simulation for Aurora
│   └── tools/
│       ├── __init__.py      # Python package marker
│       ├── compute.py       # Globus Compute wrapper
│       └── globus_interface.py # Globus authentication interface
├── scripts/
│   ├── globus_check.py      # Environment verification
│   ├── gen_endpoint_config.py # Endpoint configuration generator
│   └── run_openmm_test.sh   # Aurora OpenMM test runner
├── tests/
│   ├── test_smoke.py        # Unit tests with mocks
│   ├── test_openmm_aurora.py # OpenMM Aurora GPU tests
│   └── README_openmm_test.md # OpenMM test documentation
├── requirements.txt         # Python dependencies (at root level)
└── README.md               # This file
```

## 🔧 Components

### LangGraph Agent (`src/main.py`)
- **Purpose**: Orchestrates the entire workflow using LangGraph state machine
- **Features**: CLI interface with argparse, structured logging, error handling, timeout management
- **Workflow**: LLM analysis → GPU simulation → Report generation
- **Model**: Configured for meta-llama/Meta-Llama-3.1-8B-Instruct on Sophia by default
- **Usage**: `python remoteGlobusToAurora/src/main.py --protein p53 --log-level INFO`

### Simulation Kernel (`src/sim_kernel.py`)
- **Purpose**: Runs molecular dynamics simulations on Aurora Intel GPUs
- **Technology**: OpenMM with OpenCL backend for Intel GPU acceleration
- **Output**: Stability metrics, RMSD calculations, potential energy, execution timings
- **Features**: Simplified protein simulation with configurable parameters

### Globus Compute Wrapper (`src/tools/compute.py`)
- **Purpose**: Manages remote job submission and monitoring via Globus Compute
- **Features**: Asynchronous job submission, timeout handling, status monitoring, error recovery
- **Target**: Aurora compute nodes with configured endpoints
- **Integration**: Handles import path resolution for remote execution

### Globus Authentication (`src/tools/globus_interface.py`)
- **Purpose**: Centralized Globus authentication for ALCF services
- **Features**: Token management, refresh handling, domain-based authentication
- **Clients**: Sophia inference service integration, Globus Compute access
- **CLI**: `python remoteGlobusToAurora/src/tools/globus_interface.py authenticate [--force]`
- **API**: Programmatic token retrieval with `get_access_token()`

### Endpoint Configuration Generator (`scripts/gen_endpoint_config.py`)
- **Purpose**: Generates Aurora-optimized Globus Compute endpoint configurations
- **Features**: PBS job configuration, virtual environment setup, PYTHONPATH management
- **Customizable**: Worker count, resource allocation, queue settings, walltime
- **Usage**: `python remoteGlobusToAurora/scripts/gen_endpoint_config.py --repo-path $PWD --venv-path $PWD/venv -o endpoint-config.yaml`

### Environment Checker (`scripts/globus_check.py`)
- **Purpose**: Comprehensive setup validation and troubleshooting
- **Checks**: Globus authentication status, endpoint connectivity, token freshness
- **Features**: Age-based token warnings, endpoint status verification
- **Usage**: `python remoteGlobusToAurora/scripts/globus_check.py`

### OpenMM Test Suite (`tests/test_openmm_aurora.py`)
- **Purpose**: Validates OpenMM functionality and GPU acceleration on Aurora
- **Features**: Platform detection, performance benchmarking, Intel GPU verification
- **Platforms**: Automatic detection, OpenCL/CPU comparison, device enumeration
- **Usage**: `python remoteGlobusToAurora/tests/test_openmm_aurora.py --platform auto --steps 1000`
- **Script**: `./remoteGlobusToAurora/scripts/run_openmm_test.sh [platform] [steps]`

## 📊 Simulation Details

The demo runs simplified molecular dynamics simulations with the following characteristics:

- **Duration**: 10,000 steps × 2 fs = 20 ps simulation time
- **Output**: Stability metrics, RMSD, potential energy
- **Hardware**: Intel GPU acceleration on Aurora
- **Timeout**: 3 minutes maximum (configurable)

> **Note**: This is a demonstration with simplified physics. Real MD simulations would require proper force fields, solvation, and longer timescales.


## 🧪 Testing

The project includes comprehensive testing for both the workflow components and Aurora GPU functionality.

### Unit Tests

```bash
# Run all unit tests
python -m pytest remoteGlobusToAurora/tests/ -v

# Run with coverage
python -m pytest remoteGlobusToAurora/tests/ --cov=src --cov-report=html

# Test individual components
python -m pytest remoteGlobusToAurora/tests/test_smoke.py::TestCLI -v
```

### OpenMM Aurora GPU Testing

Test OpenMM functionality and GPU acceleration on Aurora:

```bash
# Quick test using the shell script (recommended)
./scripts/run_openmm_test.sh

# Test with specific parameters
./scripts/run_openmm_test.sh OpenCL 2000

# Direct Python execution with platform detection
python tests/test_openmm_aurora.py --platform auto

# Test specific platform
python tests/test_openmm_aurora.py --platform OpenCL --steps 1000

# Performance benchmarking
python tests/test_openmm_aurora.py --benchmark --particles 5000

# Debug GPU detection
python tests/test_openmm_aurora.py --log-level DEBUG

# Run OpenMM tests via pytest
python -m pytest tests/test_openmm_aurora.py -v
```

### Expected Test Results

**Successful GPU Detection (Aurora Intel GPU):**
```
✓ OpenMM version: 8.1.0
✓ Available platforms: CPU, OpenCL
✓ Default platform: OpenCL
🎯 GPU DETECTED: OpenCL platform active (Aurora Intel GPU)
Performance: ~8000+ steps/second
```

**CPU Fallback:**
```
⚠️ CPU ONLY: Using CPU platform (GPU not available/detected)
Performance: ~1000-3000 steps/second
```

See `tests/README_openmm_test.md` for detailed testing documentation.

## 📦 Dependencies

The project uses the following key dependencies (see `requirements.txt`):

**Core Workflow:**
- `langgraph>=0.1.0` - State machine orchestration
- `langchain>=0.1.0` - LLM integration framework  
- `langchain-openai>=0.1.0` - OpenAI/Sophia API integration

**Globus Integration:**
- `globus-sdk>=3.0.0` - Globus authentication and services
- `globus-compute-sdk>=2.0.0` - Remote compute job submission

**Simulation:**
- `openmm>=8.0.0` - Molecular dynamics simulation engine

**Development/Testing:**
- `pytest>=7.0.0` - Unit testing framework
- `pytest-mock>=3.10.0` - Mocking capabilities
- `ruff>=0.1.0` - Code linting and formatting

## 🚀 System Requirements

- **Python**: 3.10+ (matching Aurora Framework module)
- **Platforms**: Crux (orchestration), Aurora (compute), Sophia (inference)
- **GPU**: Intel GPUs on Aurora (OpenCL support)
- **Network**: Proxy configuration for external API access from Crux

---

**Questions?** Contact the ALCF user support team or check the [ALCF documentation](https://docs.alcf.anl.gov/).
