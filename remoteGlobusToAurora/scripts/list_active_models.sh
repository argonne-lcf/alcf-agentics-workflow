#!/bin/bash

# Script to list currently running models from the Sophia cluster
# Usage: ./list_running_models.sh

# Check if access_token environment variable is set
if [ -z "$access_token" ]; then
    echo "Error: access_token environment variable is not set"
    echo "Please set it with: export access_token='your_token_here'"
    exit 1
fi

# API endpoint
API_URL="https://inference-api.alcf.anl.gov/resource_server/sophia/jobs"

# Make the API call and check if it was successful
response=$(curl -s -X GET "$API_URL" -H "Authorization: Bearer ${access_token}")
curl_exit_code=$?

if [ $curl_exit_code -ne 0 ]; then
    echo "Error: Failed to connect to API (curl exit code: $curl_exit_code)"
    exit 1
fi

# Check if response contains error
if echo "$response" | grep -q "error\|Error\|401\|403"; then
    echo "Error: API returned an error response"
    echo "$response"
    exit 1
fi

# Parse and display running models
echo "Currently Running Models on Sophia Cluster:"
echo "==========================================="

# Extract and format the running models information
echo "$response" | jq -r '
if .running and (.running | length > 0) then
    .running[] | 
    "Job ID: \(.["Job ID"])
Model(s): \(.Models)
Host: \(.["Host Name"])
Walltime: \(.Walltime)
Status: \(.["Model Status"])
---"
else
    "No models currently running"
end' 2>/dev/null

# If jq is not available, fall back to a simpler parsing method
if [ $? -ne 0 ]; then
    echo "Note: jq not found, using basic parsing"
    echo
    # Extract model names using grep and sed
    echo "$response" | grep -o '"Models": *"[^"]*"' | sed 's/"Models": *"\([^"]*\)"/\1/' | while read -r model; do
        echo "• $model"
    done
fi

# Show cluster status
echo
echo "Cluster Status:"
echo "==============="
free_nodes=$(echo "$response" | grep -o '"free_nodes": *[0-9]*' | grep -o '[0-9]*')
if [ -n "$free_nodes" ]; then
    echo "Free nodes: $free_nodes"
else
    echo "Could not determine free nodes"
fi
