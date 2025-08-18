# ALCF Agentic Workflow Demo

A demonstration of an end-to-end agentic workflow that showcases ALCF infrastructure integration. The workflow queries a language model on Sophia for protein analysis, then launches GPU-accelerated molecular dynamics simulations on Aurora via Globus Compute.
Questions: jchilders@anl.gov

## ⚡ Quick Start

### Prerequisites

- Python 3.8+ with virtual environment
- Access to ALCF systems: Crux, Aurora
- Globus authentication set up
- Aurora endpoint configured

### Setup Aurora with [Globus Compute Endpoint](https://globus-compute.readthedocs.io/en/latest/endpoints/endpoints.html)

```bash
ssh aurora.alcf.anl.gov          # login to Aurora
git clone <repository-url>       # checkout repo
cd <repo-path>                   # enter repo
module load frameworks           # get python in PATH
python -m venv venv              # setup virtual environment
source venv/bin/activate         # activate virtual environment
pip install -r requirements.txt  # install dependencies for the demo
pip install globus-compute-endpoint # install globus compute endpoint

# next we need to generate the endpoint configuration file for Globus Compute
python scripts/gen_endpoint_config.py --repo-path $PWD --venv-path $PWD/venv -o my-endpoint-config.yaml

globus-compute-endpoint configure --endpoint-config my-endpoint-config.yaml my-aurora-endpoint # create endpoint
globus-compute-endpoint start my-aurora-endpoint # start endpoint
# when you run this, it would print to screen the endpoint id:
#   > Starting endpoint; registered ID: <UUID>
# That <UUID> is the endpoint id. We will need this for the next step.

# verify endpoint is running
globus-compute-endpoint list
```

Now Aurora is ready to run the workflow.

> [!NOTE]
> For Globus Compute to work properly, the same MAJOR.MINOR version of python must be used on the endpoint as the one used to create the virtual environment. Currently the Aurora Framework module loads python 3.10.


### Setup Crux with [Globus Compute Endpoint](https://globus-compute.readthedocs.io/en/latest/endpoints/endpoints.html)


```bash

ssh crux.alcf.anl.gov

# if you need python
module use /soft/modulefiles/
module load spack-pe-base
module load python

# Clone and set up environment
git clone <repository-url>       # checkout repo
cd <repo-path>                   # enter repo
python -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

### Configuration on Crux Login Nodes

```bash
# Set required environment variables
export GC_ENDPOINT_ID="<UUID>" # from the previous step
export OPENAI_MODEL="meta-llama/Llama-3.3-70B-Instruct"  # Optional: defaults to this

# you need the http proxy setup to get off of Crux to the Sophia Inference Service
export http_proxy="http://proxy.alcf.anl.gov:3128"
export https_proxy="http://proxy.alcf.anl.gov:3128"

# Authenticate with Globus
python src/tools/globus_interface.py authenticate
# This will print a URL to the screen. Open this URL in a browser and login with your ALCF credentials.
# The web page will show a code. Copy this code and paste it into the terminal.


# Verify setup (note it will ask for another authentication)
python scripts/globus_check.py
```

> **Aurora Endpoint Setup**: The repository must also be available on Aurora with dependencies installed. Follow the Aurora setup instructions above for complete endpoint configuration.

### Run Demo

From the login nodes on Crux

```bash
# Basic run with p53 protein
python src/main.py --protein p53
```

This will run the workflow and print the results to the screen.

**Expected output:**
```bash
30-07 14:13 | INFO     | 🚀 Starting Agentic Workflow Demo
30-07 14:13 | INFO     | Target protein: p53
30-07 14:13 | INFO     | Model: meta-llama/Llama-3.3-70B-Instruct
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



```bash
# Custom protein and verbose logging
python src/main.py --protein insulin --log-level DEBUG

# Full options
python src/main.py \
   --protein myoglobin \
   --model gpt-4 \
   --endpoint custom-endpoint-id \
   --log-level INFO \
   --max-simulation-time 120
```

## 📁 Project Structure

```
agentic-workflow-demo/
├── src/
│   ├── main.py              # CLI orchestrator with LangGraph
│   ├── sim_kernel.py        # OpenMM simulation for Aurora
│   └── tools/
│       ├── compute.py       # Globus Compute wrapper
│       └── globus_interface.py # Globus authentication interface
├── scripts/
│   ├── globus_check.py      # Environment verification
│   └── gen_endpoint_config.py # Endpoint configuration generator
├── tests/
│   └── test_smoke.py        # Unit tests with mocks
├── requirements.txt         # Python dependencies
└── README.md               # This file
```

## 🔧 Components

### LangGraph Agent (`src/main.py`)
- **Purpose**: Orchestrates the entire workflow
- **Features**: CLI interface, logging, error handling
- **Nodes**: LLM analysis → Simulation → Report generation

### Simulation Kernel (`src/sim_kernel.py`)
- **Purpose**: Runs MD simulations on Aurora GPUs
- **Technology**: OpenMM with OpenCL for Intel GPUs
- **Output**: Energy, RMSD, stability metrics

### Globus Compute Wrapper (`src/tools/compute.py`)
- **Purpose**: Submits and monitors remote compute jobs
- **Features**: Job submission, timeout handling, status checks
- **Target**: Aurora compute nodes

### Globus Authentication (`src/tools/globus_interface.py`)
- **Purpose**: Handles Globus authentication for ALCF services
- **Features**: Token management, authentication flow, status checking
- **Usage**: `python src/tools/globus_interface.py authenticate`

### Endpoint Configuration Generator (`scripts/gen_endpoint_config.py`)
- **Purpose**: Generates Globus Compute endpoint configuration files
- **Features**: Configures Aurora-specific settings, virtual environment paths
- **Usage**: `python scripts/gen_endpoint_config.py --repo-path $PWD --venv-path $PWD/venv -o my-endpoint-config.yaml`

### Environment Checker (`scripts/globus_check.py`)
- **Purpose**: Validates setup before running workflow
- **Checks**: Python packages, Globus auth, endpoint status
- **Usage**: `python scripts/globus_check.py`

## 📊 Simulation Details

The demo runs simplified molecular dynamics simulations with the following characteristics:

- **Duration**: 10,000 steps × 2 fs = 20 ps simulation time
- **Output**: Stability metrics, RMSD, potential energy
- **Hardware**: Intel GPU acceleration on Aurora
- **Timeout**: 3 minutes maximum (configurable)

> **Note**: This is a demonstration with simplified physics. Real MD simulations would require proper force fields, solvation, and longer timescales.


## 🧪 Testing

```bash
# Run unit tests
python -m pytest tests/ -v

# Run with coverage
python -m pytest tests/ --cov=src --cov-report=html

# Test individual components
python -m pytest tests/test_smoke.py::TestCLI -v
```



---

**Questions?** Contact the ALCF user support team or check the [ALCF documentation](https://docs.alcf.anl.gov/).
