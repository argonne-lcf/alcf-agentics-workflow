#!/bin/bash

# Detect system
hostname=$(hostname -f)
if [[ $hostname == *"aurora"* ]]; then
    SYSTEM="aurora"
elif [[ $hostname == *"polaris"* ]]; then
    SYSTEM="polaris"
else
    echo "Unknown system. Run this script on a compute node of Aurora or Polaris."
    return 1 2>/dev/null || exit 1
fi
echo "Launching vLLM server on $SYSTEM"

# Check that the HF variables are set correctly
if [ -z "$HF_HOME" ]; then
    echo "Env variable HF_HOME needs to be set. Example: export HF_HOME='/path/to/model/weights'"
    return 1 2>/dev/null || exit 1
fi
if [ -z "$HF_DATASETS_CACHE" ]; then
    echo "Env variable HF_DATASETS_CACHE needs to be set. Example: export HF_DATASETS_CACHE='/path/to/model/weights'"
fi
if [ -z "$HF_TOKEN" ]; then
    echo "Env variable HF_TOKEN needs to be set. Example: export HF_TOKEN='your_HuggingFace_token'"
    return 1 2>/dev/null || exit 1
fi
echo "HF_HOME: $HF_HOME"
echo "HF_DATASETS_CACHE: $HF_DATASETS_CACHE"
echo "HF_TOKEN: $HF_TOKEN"

# Check that MODEL is set and that the model is available
if [ -z "$MODEL" ]; then
    echo "Env variable MODEL needs to be set. Example: export MODEL='meta-llama/Meta-Llama-3.1-8B-Instruct'"
    return 1 2>/dev/null || exit 1
fi
echo "Loading model: $MODEL"

if [ $SYSTEM == "polaris" ]; then
    # Start vllm server
    vllm serve $MODEL --port 8000 --dtype bfloat16 --tensor-parallel-size 1 --trust-remote-code > vllm_server.log 2>&1 &

    # Needed to connect to the vllm server
    unset http_proxy
    unset https_proxy
    unset HTTP_PROXY
    unset HTTPS_PROXY

elif [ $SYSTEM == "aurora" ]; then
    # Sent some env variables
    export ZE_FLAT_DEVICE_HIERARCHY=FLAT
    unset ONEAPI_DEVICE_SELECTOR
    unset CCL_PROCESS_LAUNCHER
    export CCL_PROCESS_LAUNCHER=None
    export VLLM_WORKER_MULTIPROC_METHOD=spawn
    export FI_MR_CACHE_MONITOR=userfaultfd
    export TOKENIZERS_PARALLELISM=false
    export VLLM_LOGGING_LEVEL=DEBUG

    # Start vllm server
    vllm serve $MODEL --port 8000 --dtype bfloat16 --tensor-parallel-size 1 --trust-remote-code > vllm_server.log 2>&1 &

    # Needed to connect to the vllm server
    unset http_proxy
    unset https_proxy
    unset HTTP_PROXY
    unset HTTPS_PROXY
fi
