# MACE MCP Server

A [Model Context Protocol (MCP)](https://modelcontextprotocol.io/) server that exposes [MACE](https://github.com/ACEsuit/mace) machine-learning interatomic potential calculations as tools. Connect any MCP-compatible AI client to run molecular energy calculations without writing chemistry code.

---

## Table of Contents

- [Prerequisites](#prerequisites)
- [1. Install OpenCode](#1-install-opencode)
- [2. Set Up the MCP Server](#2-set-up-the-mcp-server)
- [3. Configure ALCF Inference Endpoints](#3-configure-alcf-inference-endpoints)
  - [3a. Get an ALCF access token](#3a-get-an-alcf-access-token)
  - [3b. Add ALCF providers to your global OpenCode config](#3b-add-alcf-providers-to-your-global-opencode-config)
  - [3c. Project opencode.json (MCP server config)](#3c-project-opencodejson-mcp-server-config)
- [4. Running](#4-running)
- [Troubleshooting](#troubleshooting)
- [Learn More](#learn-more)

---

## Prerequisites

- **Python 3.10+**
- **An ALCF account** with access to [inference endpoints](https://docs.alcf.anl.gov/services/inference-endpoints/)
- **Internet access** (for PubChem lookups and first-time MACE model download)

No GPU is required -- the MACE-MP "small" model runs on CPU.

---

## 1. Install OpenCode

[OpenCode](https://opencode.ai/) is an open-source AI coding agent that supports MCP servers and can use ALCF inference endpoints as its LLM backend.

Install via the install script:

```bash
curl -fsSL https://opencode.ai/install | bash
```

Or use one of the alternative methods:

```bash
# npm
npm install -g opencode-ai

# Homebrew (macOS / Linux)
brew install anomalyco/tap/opencode
```

See the [OpenCode docs](https://opencode.ai/docs/) for the full list of installation options.

---

## 2. Set Up the MCP Server

```bash
# 1. Navigate to the server directory
cd mace-mcp-server

# 2. Create and activate a virtual environment
python -m venv .venv
source .venv/bin/activate            # On Windows: .venv\Scripts\activate

# 3. Install dependencies
pip install -r requirements.txt
```

On ALCF systems (Polaris, Sophia, etc.), load the frameworks module first:

```bash
module load frameworks
python3 -m venv .venv --system-site-packages
source .venv/bin/activate
pip install -r requirements.txt
```

### Verify the installation

```bash
python -c "from mcp.server.fastmcp import FastMCP; print('MCP SDK: OK')"
python -c "import pubchempy; print('PubChemPy: OK')"
python -c "from rdkit import Chem; print('RDKit: OK')"
python -c "import ase; print('ASE: OK')"
python -c "from mace.calculators import mace_mp; print('MACE: OK')"
```

---

## 3. Configure ALCF Inference Endpoints

ALCF provides free inference access to open-source LLMs running on dedicated hardware (Sophia and Metis clusters). These endpoints are OpenAI-compatible, so OpenCode can use them as a custom provider.

Full details: [ALCF Inference Endpoints documentation](https://docs.alcf.anl.gov/services/inference-endpoints/)

### 3a. Get an ALCF access token

```bash
# Download the authentication helper script
wget https://raw.githubusercontent.com/argonne-lcf/inference-endpoints/refs/heads/main/inference_auth_token.py

# Authenticate with your Globus/ALCF account (opens a browser)
python inference_auth_token.py authenticate

# Export the token as an environment variable
export ALCF_INFERENCE_TOKEN=$(python inference_auth_token.py get_access_token)
```

> **Note:** Access tokens are valid for **48 hours**. The `get_access_token` command automatically refreshes expired tokens. Re-authentication is required every 30 days. If you encounter permission errors, log out at [app.globus.org/logout](https://app.globus.org/logout) and re-run `authenticate --force`.

### 3b. Add ALCF providers to your global OpenCode config

OpenCode reads provider configuration from `~/.config/opencode/opencode.json` (global, applies to all projects). Add the ALCF providers there:

```bash
mkdir -p ~/.config/opencode
```

Create or edit `~/.config/opencode/opencode.json`:

```json
{
  "$schema": "https://opencode.ai/config.json",
  "provider": {
    "alcf_sophia": {
      "npm": "@ai-sdk/openai-compatible",
      "name": "Sophia",
      "options": {
        "baseURL": "https://inference-api.alcf.anl.gov/resource_server/sophia/vllm/v1",
        "apiKey": "{env:ALCF_INFERENCE_TOKEN}"
      },
      "models": {
        "openai/gpt-oss-120b": {
          "name": "gpt-oss-120b"
        },
        "openai/gpt-oss-20b": {
          "name": "gpt-oss-20b"
        }
      }
    },
    "alcf_metis": {
      "npm": "@ai-sdk/openai-compatible",
      "name": "Metis",
      "options": {
        "baseURL": "https://inference-api.alcf.anl.gov/resource_server/metis/api/v1",
        "apiKey": "{env:ALCF_INFERENCE_TOKEN}"
      },
      "models": {
        "gpt-oss-120b": {
          "name": "gpt-oss-120b"
        }
      }
    }
  }
}
```

> **Tip:** A copy of this config is available at [`config/opencode.json`](config/opencode.json) for reference.

**Key points:**

- **`{env:ALCF_INFERENCE_TOKEN}`** reads the token from your environment variable (set in step 3a).
- **Two clusters are available:**

  | Cluster | Framework | Supported Endpoints | Notes |
  |---------|-----------|-------------------|-------|
  | **Sophia** | vLLM | Chat, completions, embeddings, batch | Broader model selection, tool calling support |
  | **Metis** | SambaNova | Chat completions only | Fewer models, no batch or tool calling |

- See the [ALCF Inference Endpoints docs](https://docs.alcf.anl.gov/services/inference-endpoints/#available-models) for the full list of available models.

### 3c. Project `opencode.json` (MCP server config)

The project already includes an `opencode.json` at the repository root that configures the MACE MCP server:

```json
{
  "$schema": "https://opencode.ai/config.json",
  "mcp": {
    "chemistry": {
      "type": "local",
      "command": [".venv/bin/python", "src/server.py"],
      "enabled": true
    }
  }
}
```

This uses the `.venv` Python interpreter created in step 2. OpenCode merges this project config with your global config, so the ALCF providers and the MCP server are both available when you run `opencode` from this directory.

---

## 4. Running

### Option 1: Via OpenCode (recommended)

With `opencode.json` configured, simply run OpenCode from the project directory:

```bash
opencode
```

OpenCode will automatically:
- Connect to the ALCF inference endpoint as the LLM backend
- Launch the MACE MCP server as a subprocess
- Make the chemistry tools available to the LLM

You can then ask questions like:
- "What is the energy of a water molecule?"
- "Optimise the geometry of ethanol using MACE"
- "Compare the energies of methane and ethane"

Use `/models` in the OpenCode TUI to switch between available ALCF models.

### Option 2: Example Client (no LLM required)

Run the demo client that exercises all three tools end-to-end:

```bash
python src/example_client.py
```

This connects to the MCP server over stdio and runs demos for ethanol, aspirin, and water.

---

## Troubleshooting

### "pubchempy is not installed"

```bash
pip install pubchempy
```

### "RDKit is not installed"

```bash
pip install rdkit-pypi
```

If that fails on your platform, try installing via conda:
```bash
conda install -c conda-forge rdkit
```

### "MACE is not installed"

```bash
pip install mace-torch
```

MACE requires PyTorch. If you don't have it:
```bash
pip install torch --index-url https://download.pytorch.org/whl/cpu  # CPU only
pip install mace-torch
```

### First run is slow

On the first call to `run_mace_calculation`, the MACE-MP model weights are downloaded (~100 MB for "small"). Subsequent runs use the cached model.

### "Failed to generate 3D coordinates"

Some complex SMILES strings may fail RDKit's distance-geometry embedding. Try:
- Simplifying the molecule
- Using a different `random_seed`
- Providing your own XYZ file directly to `run_mace_calculation`

### ALCF token expired

Access tokens expire after 48 hours. Re-export:

```bash
export ALCF_INFERENCE_TOKEN=$(python inference_auth_token.py get_access_token)
```

If you get permission errors after 30 days, re-authenticate:

```bash
python inference_auth_token.py authenticate --force
```

---

## Learn More

- **[ALCF Inference Endpoints](https://docs.alcf.anl.gov/services/inference-endpoints/)** -- Full documentation for ALCF's inference service, available models, and API usage
- **[OpenCode documentation](https://opencode.ai/docs/)** -- Installation, configuration, providers, and MCP server setup
- **[OpenCode providers](https://opencode.ai/docs/providers/)** -- How to configure custom OpenAI-compatible providers (like ALCF)
