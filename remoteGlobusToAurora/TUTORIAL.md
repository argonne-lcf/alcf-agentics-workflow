# Remote Agent-Driven Simulation on Aurora: Comprehensive Tutorial

Welcome to the comprehensive tutorial for the **remoteGlobusToAurora** workflow. This tutorial is designed for new ALCF users who want to understand and effectively use this end-to-end agentic workflow that integrates a client system (Crux/laptop), inference system (Sophia), and host compute system (Aurora/Polaris) systems.

## 📚 Table of Contents

1. [Overview and Architecture](#overview-and-architecture)
2. [Prerequisites and System Requirements](#prerequisites-and-system-requirements)
3. [Understanding the Components](#understanding-the-components)
4. [Step-by-Step Setup Guide](#step-by-step-setup-guide)
5. [Running Your First Workflow](#running-your-first-workflow)
6. [Understanding the Code](#understanding-the-code)
7. [Customization and Extension](#customization-and-extension)
8. [Troubleshooting Guide](#troubleshooting-guide)
9. [Best Practices](#best-practices)
10. [Additional Resources](#additional-resources)

## 🏗️ Overview and Architecture

### What This Workflow Does

The **remoteGlobusToAurora** workflow demonstrates a sophisticated integration of ALCF infrastructure components:

1. **Orchestration on Host**: A LangGraph-based agent coordinates the entire workflow
2. **AI Analysis on Inference System**: Queries large language models for protein analysis
3. **GPU Simulation on Target Compute System**: Executes OpenMM molecular dynamics simulations
4. **Seamless Integration**: Uses Globus Compute for secure, authenticated remote execution

### Architecture Diagram

```
┌─────────────┐    ┌──────────────┐    ┌──────────────┐
│    Client   │───▶│   Inference  │    │ Host Compute │
│ (LangGraph  │    │   (LLM AI)   │    │  (GPU Sim)   │
│Orchestrator)│    │              │    │              │
└─────────────┘    └──────────────┘    └──────────────┘
       │                    │                  ▲
       │                    │                  │
       └────────────────────┼──────────────────┘
                            │
                    ┌──────────────┐
                    │Globus Compute│
                    │   Gateway    │
                    └──────────────┘
```

### Key Technologies

- **[LangGraph](https://langchain-ai.github.io/langgraph/)**: A library for building stateful, multi-actor applications with LLMs
- **[Globus Compute](https://globus-compute.readthedocs.io/)**: Secure function-as-a-service for research computing
- **[OpenMM](http://openmm.org/)**: A high-performance molecular simulation library

## 📋 Prerequisites and System Requirements

### Essential Requirements

1. **ALCF Account Access**:
   - Active account at ALCF: for using Crux, Aurora, Polaris,and Sophia systems
   - Valid Globus identity linked to your ALCF credentials

2. **Python Environment**:
   - Python 3.X, where X is matching on your host and target systems (Globus Requirement)
   - Ability to create virtual environments

3. **Network Access**:
   - Outbound HTTP/HTTPS from client and host system (maybe via proxy)
   - Globus authentication endpoints
   - Host compute node access

### System-Specific Requirements

#### Client (Orchestration)
- Login node access
- Proxy configuration for external API calls
- Python virtual environment capability

#### Host Compute System (Aurora or Polaris)
- Compute node job submission permissions
- Globus Compute endpoint configuration
- OpenMM with compute system's GPU support

#### Inference System (Sophia)
- API access to inference service via ALCF or ANL account
- Globus authentication integration

## 🔧 Understanding the Components

### 1. LangGraph Agent (`src/main.py`)

**Purpose**: Orchestrates the entire workflow using a state machine pattern.

**Key Features**:
- **State Management**: Maintains workflow state across distributed operations
- **Node-Based Processing**: Each step is a discrete node with clear inputs/outputs
- **Error Handling**: Graceful failure handling and recovery
- **CLI Interface**: User-friendly command-line interface

**LangGraph Concepts**:
- **StateGraph**: Defines the workflow structure
- **Nodes**: Individual processing steps (LLM analysis, simulation, reporting)
- **Edges**: Transitions between workflow steps
- **State**: Shared data structure passed between nodes

### 2. Simulation Kernel (`src/sim_kernel.py`)

**Purpose**: Executes molecular dynamics simulations using OpenMM on Host Compute System.

**Key Features**:
- **GPU Acceleration**: Utilizes Intel GPUs via OpenCL
- **Fallback Support**: Gracefully handles missing dependencies
- **Platform Detection**: Automatically selects best available compute platform
- **Realistic Metrics**: Calculates RMSD, energy, and stability scores

### 3. Globus Compute Wrapper (`src/tools/compute.py`)

**Purpose**: Manages remote job submission and monitoring.

**Key Features**:
- **Asynchronous Execution**: Non-blocking job submission
- **Timeout Management**: Configurable execution limits
- **Status Monitoring**: Real-time job status tracking
- **Error Recovery**: Handles network and execution failures

### 4. Authentication Interface (`src/tools/globus_interface.py`)

**Purpose**: Centralized Globus authentication for ALCF services.

**Key Features**:
- **Token Management**: Automatic refresh and validation
- **Domain Restriction**: ALCF-specific identity providers
- **CLI Integration**: Standalone authentication tools
- **Error Handling**: User-friendly error messages and guidance

## 🚀 Step-by-Step Setup Guide

### Phase 1: Host Endpoint Setup

> **Important**: Complete host setup first, as the endpoint ID is needed for client configuration.

#### 1.1 Host Compute System Access

On Aurora:

```bash
# Login to Host Compute System
ssh <username>@aurora.alcf.anl.gov

# Load required modules (on Aurora only)
module load frameworks  # Provides Python

```

On Polaris:

```bash
# Login to Host Compute System
ssh <username>@polaris.alcf.anl.gov

# Load required modules (on Polaris only)
module load conda
conda activate
python --version # 3.11.X on Polaris as of 9/22/2025
```

#### 1.2 Environment Setup

```bash
# Clone the repository
git clone https://github.com/argonne-lcf/alcf-agentics-workflow.git
cd alcf-agentics-workflow

# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
pip install globus-compute-endpoint
python --version # 3.10.X on Polaris as of 9/22/2025
```

#### 1.3 Generate Endpoint Configuration

```bash
# Generate optimized Aurora or Polaris configuration
python remoteGlobusToAurora/scripts/gen_endpoint_config.py \
  --repo-path $PWD \
  --venv-path $PWD/venv \
  -o my-endpoint-config.yaml
```

**Understanding the Configuration**:
- **PBS Integration**: Configures PBS job scheduling
- **GPU Resources**: Allocates Intel GPU resources
- **Environment Setup**: Configures Python paths and virtual environment
- **Resource Limits**: Sets reasonable defaults for demo usage
- **user_init**: Sets the custom environment (think bash submission script) for the application to run on the compute nodes

#### 1.4 Create and Start Endpoint

```bash
# Create endpoint (one-time setup)
globus-compute-endpoint configure --endpoint-config my-polaris-config.yaml my-polaris-endpoint

# Start the endpoint
globus-compute-endpoint start my-polaris-endpoint

# Check endpoints running (helpful in the future to check if endpoint is running)
globus-compute-endpoint list
```

**Expected Output**:
```
Starting endpoint; registered ID: <UUID>
```

> **Save the Endpoint ID**: You'll need this UUID for Crux configuration.

#### 1.5 Verify Compute Setup (optional)

```bash
# get interactive compute node
qsub -I -q debug -l select=1 -l walltime=00:10:00 -A <project> -l filesystems=home
cd /path/to/alcf_agentics_workflow
source venv/bin/activate

# Test OpenMM functionality
python remoteGlobusToAurora/tests/test_openmm_aurora.py --platform auto
```

### Phase 2: Agent Orchestration Setup

#### 2.1 Client System Setup

On a laptop:
```bash
# Clone the repository
git clone https://github.com/argonne-lcf/alcf-agentics-workflow.git
cd alcf-agentics-workflow

# Create virtual environment
# NOTE: Python version (Major.Minor values) must match the one on the host compute system
python -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

On a Crux login node:

```bash
# Login to Crux
ssh <username>@crux.alcf.anl.gov

# Load Aurora compatible Python from modules
module use /soft/modulefiles/
module load spack-pe-base
module load python

# Install proper Python version from Miniconda for Polaris compatible Python
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
source $HOME/miniconda3/bin/activate
conda install python==3.11.12

# Clone the repository
git clone https://github.com/argonne-lcf/alcf-agentics-workflow.git
cd alcf-agentics-workflow

# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

#### 2.2 Set Environment Variables

```bash
# Essential: Aurora endpoint ID from Phase 1
export GC_ENDPOINT_ID=<UUID>  # Your actual ID from Phase 1

# Optional: LLM model selection
export OPENAI_MODEL=meta-llama/Meta-Llama-3.1-8B-Instruct

# Required on Crux: Proxy configuration for Sophia access
export http_proxy=http://proxy.alcf.anl.gov:3128
export https_proxy=http://proxy.alcf.anl.gov:3128
```

#### 2.3 Globus Authentication

```bash
# Authenticate with Globus
python remoteGlobusToAurora/src/tools/globus_interface.py authenticate

# Follow the authentication flow:
# 1. Open the provided URL in a browser
# 2. Login with ALCF credentials  
# 3. Copy the authorization code
# 4. Paste it in the terminal
```

**Troubleshooting Authentication**:
- Logout of Globus completely: https://app.globus.org/logout
- Use `--force` flag to force re-authentication: `authenticate --force`
- Ensure you're using an ALCF domain account (@anl.gov or @alcf.anl.gov)

#### 2.4 Verify Complete Setup

```bash
# Comprehensive environment check
python remoteGlobusToAurora/scripts/globus_check.py

# Expected output:
# ✅ Python Environment
# ✅ Globus Auth  
# ✅ Compute Endpoint
```

## 🎯 Running Your First Workflow

### Basic Execution

```bash
# Navigate to project directory
cd <repository-path>
source venv/bin/activate

# Run with default protein (p53)
python remoteGlobusToAurora/src/main.py --protein p53
```

### Understanding the Output

```bash
30-07 14:13 | INFO     | 🚀 Starting Agentic Workflow Demo
30-07 14:13 | INFO     | Target protein: p53
30-07 14:13 | INFO     | Model: meta-llama/Meta-Llama-3.1-8B-Instruct
30-07 14:13 | INFO     | Endpoint: <UUID>
30-07 14:13 | INFO     | 🔄 Running workflow...
30-07 14:13 | INFO     | 🧠 Querying LLM for protein analysis...
30-07 14:13 | INFO     | ✅ LLM analysis completed
30-07 14:13 | INFO     | 🚀 Launching GPU simulation on Aurora...
30-07 14:13 | INFO     | Submitting simulation job to endpoint
30-07 14:13 | INFO     | Job submitted with ID: None
30-07 14:14 | INFO     | Job completed successfully in 60.1s
30-07 14:14 | INFO     | ✅ Simulation completed: completed
30-07 14:14 | INFO     | 📊 Generating final report...
30-07 14:14 | INFO     | ✅ Report generated
30-07 14:14 | INFO     | ⏱️  Workflow completed in 77.5s
```

### Advanced Usage Examples

```bash
# Custom protein with verbose logging
python remoteGlobusToAurora/src/main.py \
  --protein insulin \
  --log-level DEBUG

# Different model and longer timeout
python remoteGlobusToAurora/src/main.py \
  --protein myoglobin \
  --model meta-llama/Meta-Llama-3.1-70B-Instruct \
  --max-simulation-time 300

# Using custom endpoint
python remoteGlobusToAurora/src/main.py \
  --protein lysozyme \
  --endpoint custom-endpoint-uuid
```

## 💻 Understanding the Code

### LangGraph Workflow Implementation

The core of this demo is a **LangGraph workflow** that orchestrates three sequential steps: LLM analysis, GPU simulation, and report generation. Let's explore how it's built and what each component does.

#### Building the Workflow Graph

The workflow is defined in the `create_workflow()` function:

```python
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
```

This creates a linear workflow graph that looks like:

```
    START
      ↓
┌─────────────┐
│LLM Analysis │ ← Queries Sophia for protein analysis
└─────────────┘
      ↓
┌─────────────┐
│ Simulation  │ ← Runs OpenMM on Aurora via Globus Compute
└─────────────┘
      ↓
┌─────────────┐
│   Report    │ ← Generates markdown summary
└─────────────┘
      ↓
     END
```

**Key LangGraph Concepts**:
- **StateGraph**: Manages shared state between nodes
- **Nodes**: Individual processing functions that transform state
- **Edges**: Define execution flow between nodes
- **Entry Point**: Where execution begins
- **State Persistence**: Each node receives and returns the complete state

#### Node Implementations

Each node is a Python function that follows the pattern: `(state: AgentState) -> AgentState`

##### 1. LLM Analysis Node

**Purpose**: Query Sophia LLM for protein analysis and simulation parameters

```python
# Key authentication setup
access_token = get_access_token()
llm = ChatOpenAI(
   model=model,
   openai_api_key=access_token,
   openai_api_base=SOPHIA_BASE_URL
)

# Structured prompting for consistent responses
prompt = f"""
Analyze the protein {state['protein']} and suggest molecular dynamics simulation parameters.
Please provide your response as JSON with the following structure:
{{
   "timestep": 0.002,
   "temperature": 300,
   "steps": 10000,
   "protein": "{state['protein']}"
}}
"""

response = llm.invoke(prompt)
```

**What it does**:
- Authenticates with Sophia using Globus tokens
- Sends structured prompt requesting simulation parameters
- Parses LLM response or falls back to sensible defaults
- Updates state with `simulation_params` for next node

##### 2. Simulation Node

**Purpose**: Execute molecular dynamics simulation on Aurora GPUs

```python
# Remote execution via Globus Compute
gc_wrapper = GlobusComputeWrapper()
timeout = state.get("max_simulation_time", 180)
result = gc_wrapper.submit_simulation(state["simulation_params"], timeout=timeout)

state["simulation_result"] = result
```

**What it does**:
- Takes simulation parameters from previous node
- Submits OpenMM job to Aurora via Globus Compute
- Waits for completion with configurable timeout
- Stores results (energy, RMSD, stability metrics) in state

##### 3. Report Generation Node

**Purpose**: Create human-readable markdown report

```python
# Generate formatted report based on simulation outcome
if "error" in sim_result:
   report = f"""
## Molecular Dynamics Analysis Report - {protein}
**Status**: ❌ FAILED
**Error**: {sim_result['error']}
"""
else:
   report = f"""
## Molecular Dynamics Analysis Report - {protein}
**Status**: ✅ COMPLETED
### Results
- Stability Score: {sim_result.get("stability_score", "N/A")}
- RMSD: {sim_result.get("rmsd", "N/A")} Å
- Final Energy: {sim_result.get("final_energy", "N/A")} kJ/mol
"""
```

**What it does**:
- Formats simulation results into readable markdown
- Handles both success and error cases gracefully
- Includes timestamp and parameter summary
- Stores final report in state for output

### Globus Compute Integration

The `GlobusComputeWrapper` class handles all remote execution:

```python
class GlobusComputeWrapper:
    def __init__(self, endpoint_id: Optional[str] = None):
        """Initialize with endpoint ID from environment or parameter"""
        self.endpoint_id = endpoint_id or os.getenv("GC_ENDPOINT_ID")
        self.client = Client()
        self.executor = Executor(endpoint_id=self.endpoint_id)
    
    def submit_simulation(self, params: Dict[str, Any], timeout: int = 180) -> Dict[str, Any]:
        """Submit MD simulation job with timeout handling"""
        
        # Import the simulation function
        from sim_kernel import run_md_simulation
        
        # Submit to Globus Compute
        future = self.executor.submit(run_md_simulation, params)
        
        # Wait for completion with timeout
        start_time = time.time()
        while not future.done():
            elapsed = time.time() - start_time
            if elapsed > timeout:
                return {"status": "timeout", "error": f"Job timed out after {timeout} seconds"}
            time.sleep(5)
        
        # Return results
        result = future.result()
        return {"status": "completed", "execution_time": elapsed, **result}
```

### OpenMM Simulation Details

The simulation kernel implements realistic molecular dynamics:

```python
def _run_openmm_simulation(protein: str, timestep: float, temperature: float, steps: int):
    """Run OpenMM molecular dynamics simulation"""
    
    # Create molecular system based on protein
    system, topology = _create_demo_system(protein)
    
    # Setup platform (prefer OpenCL for Aurora GPUs)
    platform, platform_name = _setup_openmm_platform()
    
    # Configure integrator
    integrator = openmm.LangevinMiddleIntegrator(
        temperature * unit.kelvin,
        1.0 / unit.picosecond,  # friction
        timestep * unit.picoseconds
    )
    
    # Create simulation
    simulation = app.Simulation(topology, system, integrator, platform)
    
    # Set initial positions and minimize energy
    positions = _get_initial_positions(protein)
    simulation.context.setPositions(positions)
    simulation.minimizeEnergy()
    
    # Run molecular dynamics
    simulation.step(steps)
    
    # Calculate metrics
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    rmsd = _calculate_rmsd(initial_positions, final_positions)
    stability_score = _calculate_stability_score(rmsd, initial_energy, final_energy)
    
    return {
        "final_energy": final_energy.value_in_unit(unit.kilojoule_per_mole),
        "rmsd": rmsd,
        "stability_score": stability_score,
        "platform_used": platform_name,
        "steps_per_second": steps / md_time
    }
```

## 🔧 Customization and Extension

### Adding New Proteins

To add support for new proteins:

1. **Extend the simulation kernel**:
```python
def _create_demo_system(protein: str):
    """Create molecular system based on protein name"""
    if protein.lower() == "hemoglobin":
        return _create_large_protein_system(num_atoms=50)
    elif protein.lower() in ["p53", "insulin", "lysozyme"]:
        return _create_chain_system(num_atoms=10)
    else:
        return _create_simple_system()
```

2. **Update initial positions**:
```python
def _get_initial_positions(protein: str):
    """Get initial positions for the system"""
    if protein.lower() == "hemoglobin":
        # Large protein configuration
        positions = []
        for i in range(50):
            x = (i % 10) * 0.15
            y = (i // 10) * 0.15
            z = 0.0
            positions.append(openmm.Vec3(x, y, z) * unit.nanometer)
        return positions
    # ... existing configurations
```

### Adding Custom Analysis Steps

To add new workflow nodes:

1. **Define the node function**:
```python
def custom_analysis_node(state: AgentState) -> AgentState:
    """Custom analysis step"""
    # Your custom logic here
    protein = state.get("protein")
    
    # Perform analysis
    analysis_result = perform_custom_analysis(protein)
    
    # Update state
    state["custom_analysis"] = analysis_result
    return state
```

2. **Update the workflow graph**:
```python
def create_workflow():
    """Create extended workflow with custom node"""
    workflow = StateGraph(AgentState)
    
    # Add nodes
    workflow.add_node("llm_analysis", llm_analysis_node)
    workflow.add_node("custom_analysis", custom_analysis_node)  # New node
    workflow.add_node("simulation", simulation_node)
    workflow.add_node("report", report_node)
    
    # Define edges
    workflow.set_entry_point("llm_analysis")
    workflow.add_edge("llm_analysis", "custom_analysis")  # New edge
    workflow.add_edge("custom_analysis", "simulation")    # New edge
    workflow.add_edge("simulation", "report")
    workflow.add_edge("report", END)
    
    return workflow.compile()
```

### Customizing LLM Prompts

Modify the LLM analysis prompt for different use cases:

```python
def llm_analysis_node(state: AgentState) -> AgentState:
    """Enhanced LLM analysis with custom prompts"""
    
    protein = state['protein']
    analysis_type = state.get('analysis_type', 'standard')
    
    if analysis_type == 'drug_design':
        prompt = f"""
        Analyze the protein {protein} for drug design applications.
        
        Focus on:
        1. Binding sites and druggable pockets
        2. Allosteric regulation opportunities
        3. Recommended simulation conditions for drug screening
        4. Key residues for molecular docking
        """
    elif analysis_type == 'folding':
        prompt = f"""
        Analyze the protein {protein} for folding studies.
        
        Focus on:
        1. Secondary structure elements
        2. Folding pathway considerations
        3. Temperature and solvent conditions
        4. Potential misfolding regions
        """
    else:
        # Standard analysis prompt
        prompt = f"""
        Analyze the protein {protein} and suggest molecular dynamics simulation parameters.
        
        Please provide:
        1. Brief protein description
        2. Recommended simulation parameters
        3. Key stability metrics to monitor
        """
    
    # Continue with LLM interaction...
```

### Advanced Endpoint Configuration

For production use, customize the Aurora endpoint configuration:

```bash
# Generate production configuration
python remoteGlobusToAurora/scripts/gen_endpoint_config.py \
  --repo-path $PWD \
  --venv-path $PWD/venv \
  --account your_project_account \
  --queue prod \
  --walltime 04:00:00 \
  --max-workers 24 \
  --max-blocks 4 \
  -o production-config.yaml
```

**Production Configuration Options**:
- **account**: Your project allocation account
- **queue**: Production queue (e.g., `prod`, `debug-prod`)
- **walltime**: Extended time limits for long simulations
- **max-workers**: Higher worker counts for parallel execution
- **max-blocks**: Multiple compute blocks for scalability

## 🔍 Troubleshooting Guide

### Common Issues and Solutions

#### 1. Authentication Problems

**Issue**: `❌ Globus authentication failed`

**Solutions**:
```bash
# Force re-authentication
python remoteGlobusToAurora/src/tools/globus_interface.py authenticate --force

# Check token status
python remoteGlobusToAurora/src/tools/globus_interface.py check_status

# Verify proxy settings on Crux
echo $https_proxy  # Should show proxy.alcf.anl.gov:3128
```

**Debugging Steps**:
1. Ensure you're using an ALCF domain account
2. Check that proxy environment variables are set
3. Verify network connectivity to Globus services
4. Clear old tokens: `rm -rf ~/.globus/app/`

#### 2. Endpoint Connectivity Issues

**Issue**: `❌ Endpoint check failed`

**Solutions**:
```bash
# Check endpoint status on Aurora
globus-compute-endpoint list
globus-compute-endpoint status my-aurora-endpoint

# Restart endpoint if needed
globus-compute-endpoint restart my-aurora-endpoint

# Check endpoint logs
globus-compute-endpoint logs my-aurora-endpoint
```

**Common Causes**:
- Aurora compute nodes are busy/unavailable
- PBS job queue is full
- Endpoint configuration errors
- Network connectivity issues

#### 3. OpenMM/GPU Problems

**Issue**: Simulation falls back to CPU or fails

**Solutions**:
```bash
# Test OpenMM on Aurora
python remoteGlobusToAurora/tests/test_openmm_aurora.py --platform OpenCL --log-level DEBUG

# Check available platforms
python -c "import openmm; print([p.getName() for p in openmm.Platform.getPlatforms()])"

# Module loading issues
module load frameworks
python -c "import openmm; print(openmm.__version__)"
```

**Expected Results**:
- OpenCL platform should be available
- GPU detection should show Intel GPU
- Performance should be >5000 steps/second

#### 4. LLM/Sophia Connection Issues

**Issue**: `❌ LLM analysis failed`

**Solutions**:
```bash
# Test Sophia connectivity
curl -H "Authorization: Bearer $(python -c 'from tools.globus_interface import get_access_token; print(get_access_token())')" \
  https://inference-api.alcf.anl.gov/resource_server/sophia/vllm/v1/models

# Check proxy settings
echo $https_proxy

# Test with different model
python remoteGlobusToAurora/src/main.py --protein p53 --model meta-llama/Meta-Llama-3.1-70B-Instruct
```

#### 5. Environment and Dependencies

**Issue**: Import errors or missing packages

**Solutions**:
```bash
# Check Python version (must be 3.10+)
python --version

# Verify all dependencies
python remoteGlobusToAurora/scripts/globus_check.py

# Reinstall dependencies
pip install --upgrade -r requirements.txt

# Check virtual environment
which python  # Should show venv path
```

### Debugging Workflow Execution

#### Enable Debug Logging

```bash
# Run with maximum verbosity
python remoteGlobusToAurora/src/main.py \
  --protein p53 \
  --log-level DEBUG
```

#### Monitor Job Progress

```bash
# Watch endpoint logs in real-time
tail -f ~/.globus/compute/my-aurora-endpoint/endpoint.log

# Check PBS job status on Aurora
qstat -u $USER
```

#### Analyze State Transitions

Add custom logging to understand state flow:

```python
def debug_state(state: AgentState, node_name: str):
    """Debug helper to log state transitions"""
    logger = get_logger()
    logger.debug(f"=== {node_name} State ===")
    logger.debug(f"Protein: {state.get('protein')}")
    logger.debug(f"Has simulation_params: {'simulation_params' in state}")
    logger.debug(f"Has simulation_result: {'simulation_result' in state}")
    logger.debug(f"Messages count: {len(state.get('messages', []))}")

# Use in nodes:
def llm_analysis_node(state: AgentState) -> AgentState:
    debug_state(state, "LLM_ANALYSIS_START")
    # ... node logic ...
    debug_state(state, "LLM_ANALYSIS_END")
    return state
```

## 🎯 Best Practices

### 1. Resource Management

**Aurora Endpoint Management**:
```bash
# Start endpoint only when needed
globus-compute-endpoint start my-aurora-endpoint

# Stop when done to free resources
globus-compute-endpoint stop my-aurora-endpoint

# Use debug queue for testing
# Use production queues for real work
```

**Job Sizing**:
- Start with small simulations (1000-10000 steps)
- Use debug queue for initial testing
- Scale up gradually for production runs

### 2. Authentication Best Practices

**Token Management**:
```bash
# Check token freshness regularly
python remoteGlobusToAurora/scripts/globus_check.py

# Re-authenticate before long sessions
python remoteGlobusToAurora/src/tools/globus_interface.py authenticate --force

# Use environment variables for automation
export GC_ENDPOINT_ID="your-endpoint-uuid"
```

### 3. Error Handling and Monitoring

**Robust Workflow Design**:
```python
def robust_simulation_node(state: AgentState) -> AgentState:
    """Simulation node with enhanced error handling"""
    max_retries = 3
    retry_delay = 30  # seconds
    
    for attempt in range(max_retries):
        try:
            gc_wrapper = GlobusComputeWrapper()
            result = gc_wrapper.submit_simulation(
                state["simulation_params"], 
                timeout=300  # Extended timeout
            )
            
            if result.get("status") == "completed":
                state["simulation_result"] = result
                return state
            elif result.get("status") == "timeout":
                logger.warning(f"Simulation timeout on attempt {attempt + 1}")
                if attempt < max_retries - 1:
                    time.sleep(retry_delay)
                    continue
                    
        except Exception as e:
            logger.error(f"Simulation attempt {attempt + 1} failed: {e}")
            if attempt < max_retries - 1:
                time.sleep(retry_delay)
                continue
    
    # All retries failed
    state["simulation_result"] = {"error": "All simulation attempts failed"}
    return state
```

### 4. Performance Optimization

**Simulation Parameters**:
```python
# Optimize for quick testing
test_params = {
    "timestep": 0.004,  # Larger timestep
    "temperature": 300,
    "steps": 5000,      # Fewer steps
    "protein": "p53"
}

# Optimize for production
production_params = {
    "timestep": 0.002,  # Smaller timestep for accuracy
    "temperature": 300,
    "steps": 100000,    # More steps for convergence
    "protein": "p53"
}
```

**LLM Query Optimization**:
```python
# Cache LLM responses for repeated proteins
from functools import lru_cache

@lru_cache(maxsize=100)
def get_cached_protein_analysis(protein: str, model: str) -> str:
    """Cache LLM responses to avoid repeated queries"""
    # ... LLM query logic ...
    return response
```

### 5. Workflow Extension Patterns

**Conditional Branching**:
```python
def conditional_router(state: AgentState) -> str:
    """Route based on protein type"""
    protein = state.get("protein", "").lower()
    
    if protein in ["insulin", "hemoglobin"]:
        return "specialized_analysis"
    else:
        return "standard_simulation"

# Add conditional edges
workflow.add_conditional_edges(
    "llm_analysis",
    conditional_router,
    {
        "specialized_analysis": "specialized_node",
        "standard_simulation": "simulation"
    }
)
```

**Parallel Execution**:
```python
def create_parallel_workflow():
    """Workflow with parallel simulation branches"""
    workflow = StateGraph(AgentState)
    
    # Add parallel simulation nodes
    workflow.add_node("cpu_simulation", cpu_simulation_node)
    workflow.add_node("gpu_simulation", gpu_simulation_node)
    workflow.add_node("comparison", comparison_node)
    
    # Parallel execution
    workflow.add_edge("llm_analysis", "cpu_simulation")
    workflow.add_edge("llm_analysis", "gpu_simulation")
    workflow.add_edge(["cpu_simulation", "gpu_simulation"], "comparison")
    
    return workflow.compile()
```

## 📚 Additional Resources

### Documentation

**Core Technologies**:
- [LangGraph Documentation](https://langchain-ai.github.io/langgraph/) - State machine and workflow concepts
- [Globus Compute Documentation](https://globus-compute.readthedocs.io/) - Remote execution and endpoint management
- [OpenMM User Guide](http://docs.openmm.org/) - Molecular dynamics simulation
- [LangChain Documentation](https://python.langchain.com/) - LLM integration patterns

**ALCF Resources**:
- [ALCF User Guide](https://docs.alcf.anl.gov/) - System-specific documentation
- [Aurora User Guide](https://docs.alcf.anl.gov/aurora/) - Aurora-specific information
- [Globus at ALCF](https://docs.alcf.anl.gov/data-management/globus/) - ALCF Globus integration

### Advanced Topics

**LangGraph Concepts**:
- State management and persistence
- Conditional workflows and branching
- Human-in-the-loop patterns
- Multi-agent coordination

**Globus Compute Patterns**:
- Function serialization and dependencies
- Endpoint scaling and management
- Batch job submission
- Data transfer integration

**OpenMM Advanced Usage**:
- Force field selection and customization
- Enhanced sampling techniques
- Analysis and visualization
- Performance optimization

### Community and Support

**Getting Help**:
- ALCF User Support: <support@alcf.anl.gov>
- Globus Support: <support@globus.org>
- ALCF Slack channels for user discussions
- Office hours and training sessions

**Contributing**:
- Report issues and bugs in the repository
- Contribute example workflows and extensions
- Share best practices and lessons learned
- Participate in ALCF user community events

---

## 🎉 Conclusion

You've now completed the comprehensive tutorial for the **remoteGlobusToAurora** workflow! This powerful demonstration showcases how modern agentic workflows can seamlessly integrate distributed computing resources at ALCF.

### What You've Learned

1. **Architecture Understanding**: How LangGraph orchestrates distributed AI and HPC workflows
2. **System Integration**: Connecting Crux, Sophia, and Aurora through Globus Compute
3. **Practical Implementation**: Hands-on setup and execution of real workflows
4. **Troubleshooting Skills**: Diagnosing and resolving common issues
5. **Extension Patterns**: Customizing and extending workflows for your needs

### Next Steps

1. **Experiment**: Try different proteins and simulation parameters
2. **Customize**: Adapt the workflow for your specific research needs
3. **Scale**: Use production configurations for larger simulations
4. **Integrate**: Incorporate this pattern into your research workflows
5. **Share**: Contribute improvements and examples back to the community

The **remoteGlobusToAurora** workflow represents just the beginning of what's possible with agentic computing at ALCF. Use this foundation to build sophisticated, AI-driven research workflows that leverage the full power of leadership-class computing resources.

Happy computing! 🚀

---

*For questions, support, or contributions, contact: jchilders@anl.gov*
