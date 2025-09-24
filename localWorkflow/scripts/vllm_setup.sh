#!/bin/bash

# Start vllm server
vllm serve $MODEL --port 8000 --dtype float16 > vllm_server.log 2>&1 &

# Needed to connect to the vllm server
unset http_proxy
unset https_proxy
unset HTTP_PROXY
unset HTTPS_PROXY


