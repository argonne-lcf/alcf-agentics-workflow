# Local Agent Driven Simulation on Aurora

A demonstration of an end-to-end agentic workflow run within the same batch job. 
The workflow, running on Aurora or Polaris, queries a language model which is served locally on a compute node via vLLM for protein analysis, then launches GPU-accelerated molecular dynamics simulations on Aurora within the same node.
Questions: rbalin@anl.gov, htummalapalli@anl.gov

## ⚡ Quick Start

### Prerequisites

- Python 3.X with virtual environment
- Access to ALCF systems: Aurora or Polaris
- HuggingFace account and user access token

> [!WARNING]
> - A personal HuggingFace account and user access token will be needed in order to serve the LLMs locally on the systems. Please [create an account](https://huggingface.co/join) on HuggingFace and follow [instructions](https://huggingface.co/docs/hub/en/security-tokens) to create a token. A read-only token is sufficient. Additionally, to access the Meta Llama 3.1 8B model being used for this example, please [request access to the model](https://huggingface.co/meta-llama/Llama-3.1-8B-Instruct) (this should only take a few minutes and you can check the status of your request clicking on your profile settings and Gated Repositories). 
> - To download the LLama 3.1 8B model, execute the following on Aurora or Polaris. 
```bash
# Load frameworks/conda module (e.g., module load frameworks)
# Authenticate with HF (this will require your HuggingFace token)
hf auth login

# Download the model weights
export HF_HOME=/path/to/.cache
hf download meta-llama/Llama-3.1-8B-Instruct
```

### Setup on Aurora

To run the example on Aurora, execute the following steps to set up the necessary environment.

```bash
ssh <user_name>@aurora.alcf.anl.gov                                 # login to Aurora
git clone https://github.com/argonne-lcf/alcf-agentics-workflow.git # clone repo
cd alcf-agentics-workflow/localWorkflow                             # enter example dir
module load frameworks                                              # access Python and AI frameworks on Aurora
python -m venv venv                                                 # setup virtual environment
source venv/bin/activate                                            # activate virtual environment
pip install -r requirements.txt                                     # install dependencies for the demo
```

The workflow is now ready to run on Aurora.

### Setup on Polaris 

To run the example on Polaris, execute the following steps to set up the necessary environment.

```bash
ssh <user_name>@polaris.alcf.anl.gov                                # login to Polaris
git clone https://github.com/argonne-lcf/alcf-agentics-workflow.git # clone repo
cd alcf-agentics-workflow/localWorkflow                             # enter example dir
module use /soft/modulefiles                                        # access modules on Polaris
module load conda                                                   # load datascience module
conda activate                                                      # activate default conda env (get Python)
python -m venv venv                                                 # setup virtual environment
source venv/bin/activate                                            # activate virtual environment
pip install -r requirements.txt                                     # install dependencies for the demo
pip install vllm                                                    # install vllm and its dependencies (needed to locally serve LLM)
```

The workflow is now ready to run on Aurora.


### Run Demo

To run the demo, first obtain an interactive session on a compute node on either Polaris or Aurora

```bash
# On Aurora
qsub -I -l select=1,walltime=01:00:00,filesystems=home:flare -q debug -A <your_project_name>

# On Polaris
qsub -I -l select=1,walltime=01:00:00,filesystems=home:eagle -q debug -A <your_project_name>
```

and navigate to the example directory.

Then source the environment created above and export key environment variables

```bash
# On Aurora
module load frameworks
source venv/bin/activate

# On Polaris
cd alcf-agentics-workflow/localWorkflow 
module use /soft/modulefiles
module load conda
source venv/bin/activate

# On either system
export HF_HOME="/path/to/.cache"
export HF_DATASETS_CACHE="/path/to/.cache"
export HF_TOKEN="your_HuggingFace_token"
export MODEL="meta-llama/Meta-Llama-3.1-8B-Instruct"
export TMPDIR="/tmp"
```

With the correct environment, start the vLLM server on the node

```bash
source ./scripts/vllm_setup.sh

# you can tail the log file to see the output of vLLM server startup with
# tail -f vllm_server.log

# when the log file shows "Application startup complete.", the vLLM server is ready to receive requests
```

Before running the workflow, you can check that the vLLM server is running as expected by running a simple test

```bash
python tests/test_vllm_server.py
```

**Expected output:**
```bash
System: You are a helpful assistant.
User: Introduce yourself.
Assistant: I'm an artificial intelligence model ...
```

Now, you can run the agentic workflow with

```bash
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

Additional options for running the workflow are:

```bash
# Custom protein and verbose logging
python /src/main.py --protein insulin --log-level DEBUG

# Full options
python src/main.py \
   --protein myoglobin \
   --model gpt-4 \
   --endpoint custom-endpoint-id \
   --log-level INFO \
   --max-simulation-time 120
```

Once you are done running the workflow, you can teardown the vllm server with

```bash
./scripts/vllm_teardown.sh
```

## 📁 Project Structure

```
localWorkflow
├── README.md
├── requirements.txt
├── scripts
│   ├── run_openmm_test.sh
│   ├── vllm_setup.sh
│   └── vllm_teardown.sh
├── src
│   ├── main.py
│   └── sim_kernel.py
└── tests
    ├── test_openmm_aurora.py
    └── test_vllm_server.py
```

## 🔧 Components

### LangGraph Agent (`src/main.py`)
- **Purpose**: Orchestrates the entire workflow using LangGraph state machine
- **Features**: CLI interface with argparse, structured logging, error handling, timeout management
- **Workflow**: LLM analysis → GPU simulation → Report generation
- **Model**: Configured for meta-llama/Llama-2-7b-chat-hf by default
- **Usage**: `python src/main.py --protein p53 --log-level INFO`

### Simulation Kernel (`src/sim_kernel.py`)
- **Purpose**: Runs molecular dynamics simulations on Aurora Intel GPUs
- **Technology**: OpenMM with OpenCL backend for Intel GPU acceleration
- **Output**: Stability metrics, RMSD calculations, potential energy, execution timings
- **Features**: Simplified protein simulation with configurable parameters

### vLLM Server Startup Script (`scripts/vllm_setup.sh`)
- **Purpose**: Setup of vLLM server on the compute node
- **Usage**: `source scripts/vllm_setup.sh`

### vLLM Server Teardown Script (`scripts/vllm_teardown.sh`)
- **Purpose**: Tear dowm of vLLM server on the compute node
- **Usage**: `./scripts/vllm_teardown.sh`

### vLLM Server Test Script(`tests/test_vllm_server.py`)
- **Purpose**: Validates reachability and functionality of the local vLLM server
- **Usage**: `python tests/test_vllm_server.py`

### OpenMM Test Suite (`tests/test_openmm_aurora.py`)
- **Purpose**: Validates OpenMM functionality and GPU acceleration on Aurora
- **Features**: Platform detection, performance benchmarking, Intel GPU verification
- **Platforms**: Automatic detection, OpenCL/CPU comparison, device enumeration
- **Usage**: `python tests/test_openmm_aurora.py --platform auto --steps 1000`
- **Script**: `./scripts/run_openmm_test.sh [platform] [steps]`

## 📊 Simulation Details

The demo runs simplified molecular dynamics simulations with the following characteristics:

- **Duration**: 10,000 steps × 2 fs = 20 ps simulation time
- **Output**: Stability metrics, RMSD, potential energy
- **Hardware**: Intel GPU acceleration on Aurora
- **Timeout**: 3 minutes maximum (configurable)

> **Note**: This is a demonstration with simplified physics. Real MD simulations would require proper force fields, solvation, and longer timescales.


## 🧪 Testing

The project includes testing for functionality of the vLLM server and the OpenMM simulation

### vLLM server tests

```bash
# Test that the vLLM server can serve a simple request
python tests/test_vllm_server.py
```
**Successful vLLM serving:**
```
System: You are a helpful assistant.
User: Introduce yourself.
Assistant: I'm an artificial intelligence model ...
```

### OpenMM simulation tests

Test OpenMM functionality and GPU acceleration:

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

**Successful GPU Detection (Aurora Intel GPU):**
```
✓ OpenMM version: 8.1.0
✓ Available platforms: CPU, OpenCL
✓ Default platform: OpenCL
🎯 GPU DETECTED: OpenCL platform active (Aurora Intel GPU)
Performance: ~8000+ steps/second
```

## 📦 Dependencies

The project uses the following key dependencies (see `requirements.txt`):

**Core Workflow:**
- `langgraph>=0.1.0` - State machine orchestration
- `langchain>=0.1.0` - LLM integration framework  
- `langchain-openai>=0.1.0` - OpenAI/Sophia API integration

**Local LLM serving:**
- `vllm` - LLM serving through vLLM

**Simulation:**
- `openmm>=8.0.0` - Molecular dynamics simulation engine

**Development/Testing:**
- `pytest>=7.0.0` - Unit testing framework
- `pytest-mock>=3.10.0` - Mocking capabilities
- `ruff>=0.1.0` - Code linting and formatting

## 🚀 System Requirements

- **Python**: 3.10+ (matching Aurora Framework module)
- **Platforms**: Aurora or Polaris (orchestration, compute and LLM inference)
- **GPU**: Intel GPUs on Aurora (OpenCL support) or NVIDIA GPU on Polaris

---

**Questions?** Contact the ALCF user support team or check the [ALCF documentation](https://docs.alcf.anl.gov/).
