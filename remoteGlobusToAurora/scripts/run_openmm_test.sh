#!/bin/bash
# Script to run OpenMM test on Aurora
# Usage: ./run_openmm_test.sh [platform] [steps]

set -e

# Default values
PLATFORM=${1:-auto}
STEPS=${2:-1000}

echo "Running OpenMM test on Aurora..."
echo "Platform: $PLATFORM"
echo "Steps: $STEPS"
echo ""

# Activate virtual environment if available
if [ -f "./venv/bin/activate" ]; then
   echo "Activating virtual environment..."
   source ./venv/bin/activate
fi

# Change to the tests directory
cd "$(dirname "$0")/../tests"

# Run the test
echo "Starting OpenMM test..."
python test_openmm_aurora.py --platform "$PLATFORM" --steps "$STEPS" --log-level INFO

echo ""
echo "Test completed. Check the output above for GPU detection results."

