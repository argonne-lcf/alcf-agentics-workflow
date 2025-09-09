# OpenMM Aurora Test

This test verifies that OpenMM is working correctly on Aurora and can detect and use Intel GPUs via OpenCL.

## Quick Start

### On Aurora (recommended):

```bash
# Navigate to project root
cd /path/to/alcf_agentics_workflow

# Activate your virtual environment 
source ./venv/bin/activate

# Run the test script
./remoteGlobusToAurora/scripts/run_openmm_test.sh

# Or run with specific platform
./remoteGlobusToAurora/scripts/run_openmm_test.sh OpenCL 2000
```

### Direct Python execution:

```bash
# Navigate to tests directory
cd remoteGlobusToAurora/tests

# Run with auto platform detection (recommended)
python test_openmm_aurora.py

# Run with specific platform
python test_openmm_aurora.py --platform OpenCL --steps 2000

# Run with debug output
python test_openmm_aurora.py --log-level DEBUG
```

### Using pytest:

```bash
cd remoteGlobusToAurora/tests
pytest test_openmm_aurora.py -v
```

## What the Test Does

1. **Import Check**: Verifies OpenMM can be imported
2. **Platform Detection**: Lists all available OpenMM platforms
3. **Device Detection**: Determines if GPU (OpenCL) or CPU is being used
4. **Minimal Simulation**: Runs a simple 2-atom system to verify functionality
5. **Performance Report**: Reports simulation speed and timing

## Expected Output on Aurora

If OpenMM is correctly configured with OpenCL support for Intel GPUs:

```
OpenMM Aurora Test Report
============================================================

Available Platforms: CPU, OpenCL
Default Platform: OpenCL

✓ Test Status: SUCCESS
✓ Platform Used: OpenCL

Device Information:
  Platform: OpenCL
  Platform Speed: 10.0
  Double Precision: True

Performance:
  Minimization: 0.045 seconds
  MD Dynamics: 0.123 seconds  
  Performance: 8130 steps/second

🎯 GPU DETECTED: OpenCL platform active (Aurora Intel GPU)
```

If only CPU is available:

```
⚠️  CPU ONLY: Using CPU platform (GPU not available/detected)
```

## Command Line Options

- `--platform`: Choose platform (`auto`, `OpenCL`, `CPU`, `CUDA`)
- `--steps`: Number of MD steps to run (default: 1000)
- `--log-level`: Logging verbosity (`DEBUG`, `INFO`, `WARNING`, `ERROR`)

## Troubleshooting

### OpenMM Import Fails
```bash
pip install openmm
# or on Aurora:
module load openmm
```

### No OpenCL Platform
- Verify Intel GPU drivers are loaded
- Check that OpenCL runtime is available
- Try: `python -c "import openmm; print([openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])"`

### Slow Performance
- Ensure you're running on a compute node, not login node
- Check GPU allocation in your job script
- Verify OpenCL is actually using GPU with `--log-level DEBUG`

