# Agent Driven Simulation on Aurora with [PBS MCP Server](https://github.com/jtchilders/pbs-mcp-demo)

Add summary here ...

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

Make sure the envirnoment is set up correctly as described above, specifically the PBS MCP server and Globus authentication.
Then run the example with

```bash
python src/main.py
```

which will execute the default prompt to the LLM, defined to be "Check what queues are available, then generate and submit an OpenMM simulation for protein p53 on the debug queue. Monitor the job until it finishes and summarize the results."

**Expected Output**
```bash
09-04 20:35 | INFO     | 🚀 Starting PBS MCP Agentic Workflow
09-04 20:35 | INFO     | 🤖 Model: openai/gpt-oss-120b
09-04 20:35 | INFO     | 💬 Prompt: Check what queues are available, then generate and submit an OpenMM simulation for protein p53 on the debug queue. Monitor the job until it finishes and summarize the results.
<frozen runpy>:128: RuntimeWarning: 'pbs_mcp.server' found in sys.modules after import of package 'pbs_mcp', but prior to execution of 'pbs_mcp.server'; this may result in unpredictable behaviour
Processing request of type ListToolsRequest
09-04 20:35 | INFO     | 🔧 Found 11 PBS MCP tools
09-04 20:35 | INFO     | 🔧 Added 2 custom tools
09-04 20:35 | INFO     | 🔄 Sending prompt to agent...
09-04 20:35 | INFO     | 🔧 Calling list_queues
Processing request of type CallToolRequest
09-04 20:35 | INFO     | 🔧 Calling generate_sim_script
09-04 20:35 | INFO     | 🔧 Calling submit_job
Processing request of type CallToolRequest
09-04 20:35 | INFO     | 🤖 Agent response:

============================================================
The OpenMM simulation script for protein **p53** has been generated and submitted to the **debug** queue.

**Job details**
- **Job ID:** `8430395.aurora-pbs-0001.hostmgmt.cm.aurora.alcf.anl.gov`
- **Queue:** debug
- **Walltime:** 00:10:00 (default)
- **Node count:** 1
- **Account:** datascience
- **Filesystems:** home

The job is now in the PBS system. Let me know when you’d like me to check the job status, and I’ll query it for you. When the job reaches the **F** (finished) state, I’ll read the output files and provide a summary of the results.
============================================================

09-04 20:35 | INFO     | ⏳ Job 8430395 submitted. Polling every 120 seconds until complete...
09-04 20:37 | INFO     | 🔧 Calling get_job_status
Processing request of type CallToolRequest
09-04 20:37 | INFO     | 🤖 Poll 1 response:

============================================================
The job **8430395.aurora-pbs-0001.hostmgmt.cm.aurora.alcf.anl.gov** is still **queued** (`job_state: Q`). It has not finished yet, so there’s no stdout file to read at this moment.

Let me know when you’d like me to check the status again, and I’ll query the job once more. When the job reaches the `F` (finished) state, I’ll retrieve and summarize the output for you.
============================================================

9-04 20:41 | INFO     | 🔧 Calling get_job_status
Processing request of type CallToolRequest
09-04 20:41 | INFO     | 🔧 Calling read_local_file
09-04 20:41 | INFO     | 🤖 Poll 2 response:

============================================================
The PBS job **8430395** has now finished (`job_state: F`). I read its stdout file and extracted the JSON‑formatted results of the OpenMM simulation for protein **p53**.

**Simulation Summary**

| Item | Value |
|------|-------|
| **Initial Energy** | 616.41 kJ·mol⁻¹ |
| **Final Energy** | 24.29 kJ·mol⁻¹ |
| **RMSD** (root‑mean‑square deviation) | 1.5 Å |
| **Stability Score** | 0.554 |
| **Total MD Steps** | 10 000 |
| **Simulation Time (simulated)** | 0.02 ns |
| **Wall‑clock Time** | 1.9 s |
| **MD Computation Time** | 0.17 s |
| **Energy Minimization Time** | 0.14 s |
| **Performance** | 59 769 steps / s |
| **Temperature** | 300 K |
| **Timestep** | 0.002 ps |
| **Platform Used** | OpenCL |
| **OpenMM Version** | 8.3.1 |
| **Status** | completed (exit status 0) |

**Interpretation**

- The system rapidly lowered its potential energy from ~616 kJ·mol⁻¹ to ~24 kJ·mol⁻¹, indicating successful minimization and equilibration.
- An RMSD of ~1.5 Å suggests the protein remained close to its initial structure, consistent with a short 0.02 ns trajectory.
- The stability score (≈0.55) reflects moderate structural stability over the brief simulation.
- Performance of ~60 k steps s⁻¹ on a single debug‑queue node shows the OpenCL GPU backend is efficiently utilized.

The job completed successfully with no errors. Let me know if you need any additional analysis, the stderr contents, or wish to run another simulation.
============================================================

09-04 20:41 | INFO     | ✅ Job 8430395 finished after 3 poll(s).
09-04 20:41 | INFO     | ⏱️ Completed in 403.1s
```

> [!WARNING]
> - If the example raises the error `No module named pbs_ifl`, please make sure the PBS MCP server environment variables are set correctly and run the `pbs-mcp-demo/external/pbs-python-api/setup_env.sh` script for the appropriate system.

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