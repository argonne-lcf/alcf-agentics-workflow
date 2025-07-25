# Crux Environment Setup Guide

This guide helps first-time ALCF users set up the agentic workflow demo on Crux with proper Globus authentication and Aurora endpoint configuration.

## 📋 Prerequisites

- ALCF account with access to Crux and Aurora
- SSH access to ALCF systems
- Basic familiarity with command line and Python

## 🚀 Step-by-Step Setup

### 1. Connect to Crux

```bash
# SSH to Crux login node
ssh your-username@crux.alcf.anl.gov

# Navigate to your workspace
cd /lus/grand/projects/your-project/your-username
```

### 2. Load Required Modules

```bash
# Load Python and essential modules
module load conda
module load globus-compute

# Create and activate conda environment
conda create -n agentic-demo python=3.9
conda activate agentic-demo
```

### 3. Set Up Project

```bash
# Clone the repository
git clone <repository-url> agentic-workflow-demo
cd agentic-workflow-demo

# Install Python dependencies
pip install -r requirements.txt
```

### 4. Configure Globus Authentication

#### Initial Login
```bash
# Authenticate with Globus (first time only)
globus login

# This will open a browser or provide a link
# Follow the instructions to complete authentication
```

#### Verify Authentication
```bash
# Check current session
globus session show

# Should show active identity and valid tokens
python scripts/globus_check.py --verbose
```

#### Token Refresh (if needed)
```bash
# If tokens are expired (>30 days old)
globus login --force
```

### 5. Set Up Aurora Endpoint

#### Create Endpoint Configuration

```bash
# SSH to Aurora login node
ssh aurora-login

# Load Globus Compute
module load globus-compute

# Configure endpoint (first time only)
globus-compute-endpoint configure my-aurora-endpoint

# This creates ~/.globus_compute/my-aurora-endpoint/config.yaml
```

#### Edit Endpoint Configuration

Edit the configuration file to optimize for Aurora:

```yaml
# ~/.globus_compute/my-aurora-endpoint/config.yaml
display_name: "Aurora GPU Endpoint"
engine:
  type: HighThroughputEngine
  provider:
    type: PBSProProvider
    queue: gpu
    account: your-project-name
    walltime: "01:00:00"
    nodes_per_block: 1
    cpus_per_node: 12
    launcher:
      type: SingleNodeLauncher
  max_workers: 4
  worker_init: |
    module load conda
    conda activate agentic-demo
```

#### Start the Endpoint

```bash
# Start endpoint daemon
globus-compute-endpoint start my-aurora-endpoint

# Keep running in background (recommended: use tmux/screen)
tmux new-session -d -s gc-endpoint
tmux send-keys -t gc-endpoint "globus-compute-endpoint start my-aurora-endpoint" Enter

# Get endpoint ID for environment variable
globus-compute-endpoint list
# Copy the endpoint ID (format: xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx)
```

#### Resume Existing Endpoint

```bash
# If endpoint already exists but is stopped
globus-compute-endpoint resume my-aurora-endpoint

# Check status
globus-compute-endpoint status my-aurora-endpoint
```

### 6. Configure Environment Variables

Create a setup script for easy environment configuration:

```bash
# Create setup script
cat > setup_env.sh << 'EOF'
#!/bin/bash
# Agentic Workflow Demo Environment Setup

# Load modules
module load conda
module load globus-compute

# Activate environment
conda activate agentic-demo

# Set Globus Compute endpoint
export GC_ENDPOINT_ID="your-endpoint-id-here"

# Optional: Set custom LLM model
export OPENAI_MODEL="llama3-8b-instruct"

# Add project to Python path
export PYTHONPATH="${PWD}/src:${PYTHONPATH}"

echo "✅ Environment configured for agentic workflow demo"
echo "Endpoint ID: $GC_ENDPOINT_ID"
EOF

# Make executable
chmod +x setup_env.sh

# Update with your actual endpoint ID
sed -i 's/your-endpoint-id-here/ACTUAL-ENDPOINT-ID/' setup_env.sh
```

### 7. Verify Complete Setup

```bash
# Source environment
source setup_env.sh

# Run comprehensive check
python scripts/globus_check.py --verbose

# Expected output:
# INFO | 🔍 Running Globus environment checks...
# INFO | --- Python Environment ---
# INFO | ✅ All required packages available
# INFO | --- Globus Auth ---
# INFO | Checking Globus session …
# INFO | Globus tokens fresh (5d).
# INFO | --- Compute Endpoint ---
# INFO | Checking endpoint xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx …
# INFO | ✅ Endpoint reachable and active
# INFO | 📊 Summary: 3/3 checks passed
# INFO | 🎉 All checks passed - ready to run agentic workflow!
```

## 🧪 Test the Workflow

### Quick Test Run

```bash
# Simple test with p53
python src/main.py --protein p53 --log-level INFO

# Should complete in < 3 minutes with output:
# 🚀 Starting Agentic Workflow Demo
# Target protein: p53
# ...
# ⏱️ Workflow completed in XX.Xs
```

### Debug Common Issues

```bash
# Test with debug logging
python src/main.py --protein p53 --log-level DEBUG

# Check endpoint specifically
python scripts/globus_check.py --endpoint $GC_ENDPOINT_ID --verbose
```

## 🔧 Advanced Configuration

### Custom Queue Configuration

For longer simulations, modify the endpoint config:

```yaml
provider:
  type: PBSProProvider
  queue: gpu_long      # For longer jobs
  walltime: "04:00:00" # 4 hour limit
  nodes_per_block: 2   # More nodes for larger systems
```

### Multiple Endpoints

You can configure separate endpoints for different use cases:

```bash
# Create development endpoint (shorter queue)
globus-compute-endpoint configure aurora-dev
# Edit config for 'gpu_short' queue

# Create production endpoint (longer queue)  
globus-compute-endpoint configure aurora-prod
# Edit config for 'gpu_long' queue

# Use specific endpoint
export GC_ENDPOINT_ID="dev-endpoint-id"
```

### Monitoring and Logging

```bash
# Monitor endpoint logs
globus-compute-endpoint logs my-aurora-endpoint

# Check job status on Aurora
qstat -u $USER

# View detailed job information
qstat -f job-id
```

## 🚨 Troubleshooting

### Authentication Issues

```bash
# Problem: "Authentication failed"
# Solution: Re-authenticate
globus logout
globus login

# Problem: "Tokens expired"
# Check token age
python scripts/globus_check.py --max-age 1

# Refresh if needed
globus login --force
```

### Endpoint Issues

```bash
# Problem: "Endpoint unreachable"
# Check if endpoint is running
ssh aurora-login
globus-compute-endpoint status my-aurora-endpoint

# Restart if needed
globus-compute-endpoint restart my-aurora-endpoint

# Problem: "Jobs queued but not running"
# Check Aurora queue status
qstat -Q gpu
```

### Compute Issues

```bash
# Problem: "Simulation timeout"
# Increase timeout
python src/main.py --protein p53 --max-simulation-time 300

# Problem: "GPU not available"
# Check Aurora GPU availability
ssh aurora-login
nvidia-smi  # On compute node
```

## 📚 Additional Resources

- **ALCF Documentation**: https://docs.alcf.anl.gov/
- **Globus Compute Guide**: https://globus-compute.readthedocs.io/
- **Aurora User Guide**: https://docs.alcf.anl.gov/aurora/
- **Python on ALCF**: https://docs.alcf.anl.gov/software/python/

## 💡 Pro Tips

1. **Use tmux/screen** for long-running endpoints
2. **Monitor queue status** before submitting large jobs  
3. **Test with small proteins** first (insulin, ubiquitin)
4. **Keep logs** for debugging: `--log-level DEBUG`
5. **Set up aliases** for common commands:
   ```bash
   alias gc-status='globus-compute-endpoint status my-aurora-endpoint'
   alias demo-check='python scripts/globus_check.py'
   ```

---

**Need help?** Contact ALCF support at support@alcf.anl.gov with your workflow logs and endpoint ID. 