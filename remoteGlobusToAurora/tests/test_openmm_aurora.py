#!/usr/bin/env python
"""
OpenMM test for Aurora GPU detection and functionality.
Can be run standalone or as part of pytest suite.

Usage:
   # Basic functionality test
   python test_openmm_aurora.py --platform auto
   python test_openmm_aurora.py --platform OpenCL
   
   # CPU vs GPU performance benchmark
   python test_openmm_aurora.py --benchmark
   python test_openmm_aurora.py --benchmark --particles 5000 --benchmark-steps 10000
   
   # Pytest testing
   pytest test_openmm_aurora.py -v
"""

import argparse
import logging
import sys
import time
from typing import Dict, Any, Tuple, Optional

# Configure logging according to project standards
logging.basicConfig(
   level=logging.INFO,
   format="%(asctime)s | %(levelname)-8s | %(message)s",
   datefmt="%d-%m %H:%M"
)
logger = logging.getLogger(__name__)


def test_openmm_imports() -> bool:
   """Test that OpenMM can be imported successfully"""
   try:
      import openmm
      import openmm.app as app
      import openmm.unit as unit
      logger.info(f"✓ OpenMM version: {openmm.__version__}")
      return True
   except ImportError as e:
      logger.error(f"✗ OpenMM import failed: {e}")
      return False


def detect_openmm_platforms() -> Dict[str, Any]:
   """Detect available OpenMM platforms and return platform info"""
   try:
      import openmm
      
      platform_info = {
         "available_platforms": [],
         "default_platform": None,
         "platform_details": {}
      }
      
      logger.info("Detecting OpenMM platforms...")
      
      # Get all available platforms
      for i in range(openmm.Platform.getNumPlatforms()):
         platform = openmm.Platform.getPlatform(i)
         platform_name = platform.getName()
         platform_info["available_platforms"].append(platform_name)
         
         # Get platform properties
         properties = {}
         try:
            # Try to get device info if available
            if platform_name in ["OpenCL", "CUDA"]:
               # For GPU platforms, try to get device count
               try:
                  device_count = platform.getPropertyDefaultValue("DeviceIndex")
                  properties["device_count"] = device_count
               except:
                  pass
               
               # Try to get precision info
               try:
                  precision = platform.getPropertyDefaultValue("Precision")
                  properties["precision"] = precision
               except:
                  pass
         except Exception as e:
            logger.debug(f"Could not get properties for {platform_name}: {e}")
         
         platform_info["platform_details"][platform_name] = properties
         logger.info(f"  - {platform_name}: {properties}")
      
      # Get default platform (fastest available)
      try:
         default_platform = openmm.Platform.getPlatformByName("OpenCL")
         platform_info["default_platform"] = "OpenCL"
         logger.info("✓ OpenCL platform available (preferred for Aurora)")
      except:
         try:
            default_platform = openmm.Platform.getPlatformByName("CPU")
            platform_info["default_platform"] = "CPU"
            logger.info("✓ CPU platform available (fallback)")
         except:
            logger.error("✗ No usable platforms found")
      
      return platform_info
      
   except Exception as e:
      logger.error(f"✗ Platform detection failed: {e}")
      return {"available_platforms": [], "error": str(e)}


def create_minimal_system() -> Tuple[Any, Any]:
   """Create a minimal OpenMM system for testing"""
   try:
      import openmm
      import openmm.app as app
      import openmm.unit as unit
      
      logger.info("Creating minimal test system...")
      
      # Create a simple 2-particle system (like a diatomic molecule)
      system = openmm.System()
      
      # Add two particles (mass of hydrogen)
      mass = 1.0 * unit.amu
      system.addParticle(mass)
      system.addParticle(mass)
      
      # Add a harmonic bond between them
      bond_force = openmm.HarmonicBondForce()
      # k=1000 kJ/mol/nm^2, r0=0.1 nm
      bond_force.addBond(0, 1, 0.1 * unit.nanometer, 1000.0 * unit.kilojoule_per_mole / unit.nanometer**2)
      system.addForce(bond_force)
      
      # Create simple topology
      topology = app.Topology()
      chain = topology.addChain()
      residue = topology.addResidue("TST", chain)
      atom1 = topology.addAtom("H1", app.Element.getBySymbol("H"), residue)
      atom2 = topology.addAtom("H2", app.Element.getBySymbol("H"), residue)
      topology.addBond(atom1, atom2)
      
      logger.info("✓ Minimal system created successfully")
      return system, topology
      
   except Exception as e:
      logger.error(f"✗ System creation failed: {e}")
      raise


def create_benchmark_system(num_particles: int = 2000) -> Tuple[Any, Any]:
   """Create a larger system for benchmarking CPU vs GPU performance"""
   try:
      import openmm
      import openmm.app as app
      import openmm.unit as unit
      import numpy as np
      
      logger.info(f"Creating benchmark system with {num_particles} particles...")
      
      # Create system
      system = openmm.System()
      
      # Add particles in a 3D grid
      mass = 39.948 * unit.amu  # Argon mass
      for i in range(num_particles):
         system.addParticle(mass)
      
      # Add Lennard-Jones interactions (computationally intensive)
      lj_force = openmm.NonbondedForce()
      
      # Argon Lennard-Jones parameters
      sigma = 0.3405 * unit.nanometer
      epsilon = 0.996 * unit.kilojoule_per_mole
      
      for i in range(num_particles):
         # Add LJ parameters for each particle
         lj_force.addParticle(0.0, sigma, epsilon)
      
      # Create cubic box first to determine appropriate cutoff
      box_size = (num_particles / 500.0) ** (1.0/3.0) * 3.0 * unit.nanometer
      
      # Set cutoff for interactions (must be < half box size)
      max_cutoff = box_size.value_in_unit(unit.nanometer) / 2.5  # Conservative factor
      cutoff_distance = min(1.0, max_cutoff) * unit.nanometer
      
      lj_force.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
      lj_force.setCutoffDistance(cutoff_distance)
      lj_force.setUseDispersionCorrection(True)
      
      system.addForce(lj_force)
      system.setDefaultPeriodicBoxVectors(
         [box_size, 0*unit.nanometer, 0*unit.nanometer],
         [0*unit.nanometer, box_size, 0*unit.nanometer], 
         [0*unit.nanometer, 0*unit.nanometer, box_size]
      )
      
      # Create topology
      topology = app.Topology()
      chain = topology.addChain()
      residue = topology.addResidue("AR", chain)
      
      for i in range(num_particles):
         topology.addAtom(f"Ar{i}", app.Element.getBySymbol("Ar"), residue)
      
      # Set periodic box vectors in topology
      topology.setPeriodicBoxVectors(system.getDefaultPeriodicBoxVectors())
      
      logger.info(f"✓ Benchmark system created: {num_particles} Argon atoms in {box_size.value_in_unit(unit.nanometer):.2f} nm box, cutoff: {cutoff_distance.value_in_unit(unit.nanometer):.2f} nm")
      return system, topology
      
   except Exception as e:
      logger.error(f"✗ Benchmark system creation failed: {e}")
      raise


def run_performance_benchmark(num_particles: int = 2000, steps: int = 5000) -> Dict[str, Any]:
   """Run CPU vs GPU performance comparison using compute-intensive simulation"""
   
   benchmark_results = {
      "success": False,
      "cpu_results": {},
      "gpu_results": {},
      "speedup": None,
      "error": None
   }
   
   try:
      import openmm
      import openmm.app as app
      import openmm.unit as unit
      import numpy as np
      
      logger.info(f"Starting performance benchmark: {num_particles} particles, {steps} steps")
      
      # Create benchmark system
      system, topology = create_benchmark_system(num_particles)
      
      # Generate initial positions in a cubic lattice
      logger.info("Generating initial positions...")
      n_per_side = int(np.ceil(num_particles ** (1.0/3.0)))
      box_vectors = system.getDefaultPeriodicBoxVectors()
      box_size = box_vectors[0][0].value_in_unit(unit.nanometer)
      
      positions = []
      particle_count = 0
      spacing = box_size / n_per_side
      
      for i in range(n_per_side):
         for j in range(n_per_side):
            for k in range(n_per_side):
               if particle_count >= num_particles:
                  break
               x = (i + 0.5) * spacing
               y = (j + 0.5) * spacing  
               z = (k + 0.5) * spacing
               positions.append(openmm.Vec3(x, y, z) * unit.nanometer)
               particle_count += 1
            if particle_count >= num_particles:
               break
         if particle_count >= num_particles:
            break
      
      # Test CPU performance
      logger.info("Running CPU benchmark...")
      cpu_results = run_single_benchmark(system, topology, positions, "CPU", steps)
      benchmark_results["cpu_results"] = cpu_results
      
      # Test GPU performance (OpenCL)
      logger.info("Running GPU (OpenCL) benchmark...")
      try:
         gpu_results = run_single_benchmark(system, topology, positions, "OpenCL", steps)
         benchmark_results["gpu_results"] = gpu_results
      except Exception as e:
         logger.warning(f"OpenCL benchmark failed: {e}")
         benchmark_results["gpu_results"] = {
            "success": False,
            "error": f"OpenCL not available: {str(e)}"
         }
      
      # Calculate speedup
      gpu_results = benchmark_results["gpu_results"]
      if cpu_results["success"] and gpu_results["success"]:
         cpu_time = cpu_results["timing"]["dynamics_time"]
         gpu_time = gpu_results["timing"]["dynamics_time"]
         speedup = cpu_time / gpu_time if gpu_time > 0 else None
         benchmark_results["speedup"] = speedup
         benchmark_results["success"] = True
         
         logger.info(f"CPU time: {cpu_time:.3f}s, GPU time: {gpu_time:.3f}s")
         if speedup:
            logger.info(f"GPU speedup: {speedup:.2f}x faster than CPU")
      elif cpu_results["success"]:
         # CPU worked but GPU failed
         benchmark_results["success"] = True  # Partial success
         benchmark_results["error"] = f"GPU benchmark failed: {gpu_results.get('error', 'Unknown error')}"
         logger.info("CPU benchmark completed successfully, but GPU benchmark failed")
      else:
         benchmark_results["error"] = "CPU benchmark failed"
      
      return benchmark_results
      
   except Exception as e:
      logger.error(f"✗ Performance benchmark failed: {e}")
      benchmark_results["error"] = str(e)
      return benchmark_results


def run_single_benchmark(system, topology, positions, platform_name: str, steps: int) -> Dict[str, Any]:
   """Run a single benchmark on specified platform"""
   
   results = {
      "success": False,
      "platform_used": None,
      "timing": {},
      "performance_metrics": {},
      "error": None
   }
   
   try:
      import openmm
      import openmm.app as app
      import openmm.unit as unit
      
      # Set up integrator
      integrator = openmm.LangevinMiddleIntegrator(
         300 * unit.kelvin,           # temperature
         1.0 / unit.picosecond,       # friction coefficient
         0.002 * unit.picoseconds     # timestep (2 fs)
      )
      
      # Create platform
      platform = openmm.Platform.getPlatformByName(platform_name)
      properties = {}
      
      # Set platform-specific properties for better performance
      if platform_name == "CPU":
         properties["Threads"] = "0"  # Use all available CPU cores
      elif platform_name == "OpenCL":
         properties["Precision"] = "mixed"  # Good balance of speed and accuracy
         # Note: UseBlockingSync is not available on all OpenCL implementations
      
      # Create simulation
      simulation = app.Simulation(topology, system, integrator, platform, properties)
      results["platform_used"] = platform.getName()
      
      # Set positions
      simulation.context.setPositions(positions)
      
      # Energy minimization
      logger.info(f"  Minimizing energy on {platform_name}...")
      start_time = time.time()
      simulation.minimizeEnergy(maxIterations=100)
      minimize_time = time.time() - start_time
      
      # Get initial state
      initial_state = simulation.context.getState(getEnergy=True)
      initial_energy = initial_state.getPotentialEnergy()
      
      # Run equilibration (shorter)
      logger.info(f"  Running equilibration on {platform_name}...")
      start_time = time.time()
      simulation.step(min(1000, steps // 5))  # 20% of steps for equilibration
      equilibration_time = time.time() - start_time
      
      # Run production simulation
      logger.info(f"  Running {steps} production steps on {platform_name}...")
      start_time = time.time()
      simulation.step(steps)
      dynamics_time = time.time() - start_time
      
      # Get final state
      final_state = simulation.context.getState(getEnergy=True)
      final_energy = final_state.getPotentialEnergy()
      
      # Calculate performance metrics
      total_time = minimize_time + equilibration_time + dynamics_time
      steps_per_second = steps / dynamics_time if dynamics_time > 0 else 0
      ns_per_day = (steps * 0.002) / dynamics_time * 86400 if dynamics_time > 0 else 0  # 2fs timestep
      
      results.update({
         "success": True,
         "initial_energy": float(initial_energy.value_in_unit(unit.kilojoule_per_mole)),
         "final_energy": float(final_energy.value_in_unit(unit.kilojoule_per_mole)),
         "steps_completed": steps,
         "timing": {
            "minimization_time": minimize_time,
            "equilibration_time": equilibration_time,
            "dynamics_time": dynamics_time,
            "total_time": total_time
         },
         "performance_metrics": {
            "steps_per_second": steps_per_second,
            "ns_per_day": ns_per_day
         }
      })
      
      logger.info(f"  ✓ {platform_name} benchmark completed: {steps_per_second:.0f} steps/sec, {ns_per_day:.1f} ns/day")
      
      return results
      
   except Exception as e:
      logger.error(f"  ✗ {platform_name} benchmark failed: {e}")
      results["error"] = str(e)
      return results


def run_openmm_test(platform_name: str = "auto", steps: int = 1000) -> Dict[str, Any]:
   """Run a minimal OpenMM simulation to test functionality and device detection"""
   results = {
      "success": False,
      "platform_used": None,
      "device_info": {},
      "timing": {},
      "error": None
   }
   
   try:
      import openmm
      import openmm.app as app
      import openmm.unit as unit
      
      # Create system
      system, topology = create_minimal_system()
      
      # Set up integrator
      integrator = openmm.LangevinMiddleIntegrator(
         300 * unit.kelvin,    # temperature
         1.0 / unit.picosecond,  # friction
         0.002 * unit.picoseconds  # timestep
      )
      
      # Choose platform
      if platform_name == "auto":
         # Let OpenMM choose the fastest available
         platform = None
         properties = {}
         logger.info("Using auto platform selection")
      else:
         platform = openmm.Platform.getPlatformByName(platform_name)
         properties = {}
         logger.info(f"Using requested platform: {platform_name}")
      
      # Create simulation
      if platform:
         simulation = app.Simulation(topology, system, integrator, platform, properties)
         results["platform_used"] = platform.getName()
      else:
         simulation = app.Simulation(topology, system, integrator)
         results["platform_used"] = simulation.context.getPlatform().getName()
      
      logger.info(f"✓ Simulation created on platform: {results['platform_used']}")
      
      # Get device information
      context = simulation.context
      platform_used = context.getPlatform()
      
      device_info = {
         "platform_name": platform_used.getName(),
         "platform_speed": platform_used.getSpeed(),
         "supports_double_precision": platform_used.supportsDoublePrecision(),
         "properties": {}
      }
      
      # Try to get platform-specific properties
      try:
         if platform_used.getName() == "OpenCL":
            # Get OpenCL-specific information
            try:
               device_index = context.getPlatformData()
               device_info["opencl_device_index"] = device_index
            except:
               pass
            
            # Try to get device name if available
            try:
               property_names = platform_used.getPropertyNames()
               for prop_name in property_names:
                  try:
                     prop_value = platform_used.getPropertyValue(context, prop_name)
                     device_info["properties"][prop_name] = prop_value
                  except:
                     pass
            except:
               pass
               
      except Exception as e:
         logger.debug(f"Could not get platform-specific info: {e}")
      
      results["device_info"] = device_info
      
      # Set initial positions (simple linear molecule)
      positions = [
         openmm.Vec3(0.0, 0.0, 0.0) * unit.nanometer,
         openmm.Vec3(0.1, 0.0, 0.0) * unit.nanometer
      ]
      simulation.context.setPositions(positions)
      
      # Minimize energy
      logger.info("Running energy minimization...")
      start_time = time.time()
      simulation.minimizeEnergy()
      minimize_time = time.time() - start_time
      
      # Get initial energy
      initial_state = simulation.context.getState(getEnergy=True)
      initial_energy = initial_state.getPotentialEnergy()
      
      logger.info(f"Initial energy: {initial_energy.value_in_unit(unit.kilojoule_per_mole):.3f} kJ/mol")
      
      # Run dynamics
      logger.info(f"Running {steps} MD steps...")
      start_time = time.time()
      simulation.step(steps)
      dynamics_time = time.time() - start_time
      
      # Get final state
      final_state = simulation.context.getState(getEnergy=True, getPositions=True)
      final_energy = final_state.getPotentialEnergy()
      final_positions = final_state.getPositions()
      
      logger.info(f"Final energy: {final_energy.value_in_unit(unit.kilojoule_per_mole):.3f} kJ/mol")
      logger.info(f"Simulation completed in {dynamics_time:.3f} seconds")
      
      # Calculate performance metrics
      steps_per_second = steps / dynamics_time if dynamics_time > 0 else 0
      
      results.update({
         "success": True,
         "initial_energy": float(initial_energy.value_in_unit(unit.kilojoule_per_mole)),
         "final_energy": float(final_energy.value_in_unit(unit.kilojoule_per_mole)),
         "steps_completed": steps,
         "timing": {
            "minimization_time": minimize_time,
            "dynamics_time": dynamics_time,
            "steps_per_second": steps_per_second
         }
      })
      
      logger.info("✓ OpenMM simulation completed successfully")
      
      return results
      
   except Exception as e:
      logger.error(f"✗ OpenMM simulation failed: {e}")
      results["error"] = str(e)
      return results


def report_results(results: Dict[str, Any], platform_info: Dict[str, Any]) -> None:
   """Generate a comprehensive report of OpenMM test results"""
   
   print("\n" + "="*60)
   print("OpenMM Aurora Test Report")
   print("="*60)
   
   # Platform availability
   print(f"\nAvailable Platforms: {', '.join(platform_info.get('available_platforms', []))}")
   print(f"Default Platform: {platform_info.get('default_platform', 'None')}")
   
   if results["success"]:
      print(f"\n✓ Test Status: SUCCESS")
      print(f"✓ Platform Used: {results['platform_used']}")
      
      # Device information
      device_info = results["device_info"]
      print(f"\nDevice Information:")
      print(f"  Platform: {device_info['platform_name']}")
      print(f"  Platform Speed: {device_info['platform_speed']:.1f}")
      print(f"  Double Precision: {device_info['supports_double_precision']}")
      
      if device_info["properties"]:
         print(f"  Platform Properties:")
         for prop, value in device_info["properties"].items():
            print(f"    {prop}: {value}")
      
      # Performance results
      timing = results["timing"]
      print(f"\nPerformance:")
      print(f"  Minimization: {timing['minimization_time']:.3f} seconds")
      print(f"  MD Dynamics: {timing['dynamics_time']:.3f} seconds")
      print(f"  Performance: {timing['steps_per_second']:.0f} steps/second")
      
      # Energy results
      print(f"\nSimulation Results:")
      print(f"  Initial Energy: {results['initial_energy']:.3f} kJ/mol")
      print(f"  Final Energy: {results['final_energy']:.3f} kJ/mol")
      print(f"  Steps Completed: {results['steps_completed']}")
      
      # GPU detection
      if results['platform_used'] == 'OpenCL':
         print(f"\n🎯 GPU DETECTED: OpenCL platform active (Aurora Intel GPU)")
      elif results['platform_used'] == 'CPU':
         print(f"\n⚠️  CPU ONLY: Using CPU platform (GPU not available/detected)")
      else:
         print(f"\n❓ UNKNOWN: Using {results['platform_used']} platform")
         
   else:
      print(f"\n✗ Test Status: FAILED")
      print(f"✗ Error: {results.get('error', 'Unknown error')}")
   
   print("\n" + "="*60)


def report_benchmark_results(benchmark_results: Dict[str, Any], platform_info: Dict[str, Any]) -> None:
   """Generate a comprehensive report of CPU vs GPU benchmark results"""
   
   print("\n" + "="*80)
   print("OpenMM CPU vs GPU Performance Benchmark Report")
   print("="*80)
   
   # Platform availability
   print(f"\nAvailable Platforms: {', '.join(platform_info.get('available_platforms', []))}")
   
   if benchmark_results["success"]:
      print(f"\n✓ Benchmark Status: SUCCESS")
      
      cpu_results = benchmark_results["cpu_results"]
      gpu_results = benchmark_results["gpu_results"]
      speedup = benchmark_results["speedup"]
      
      print(f"\n{'='*30} CPU RESULTS {'='*30}")
      if cpu_results["success"]:
         print(f"✓ Platform: {cpu_results['platform_used']}")
         timing = cpu_results["timing"]
         perf = cpu_results["performance_metrics"]
         print(f"  Minimization:     {timing['minimization_time']:.3f} seconds")
         print(f"  Equilibration:    {timing['equilibration_time']:.3f} seconds")
         print(f"  Production MD:    {timing['dynamics_time']:.3f} seconds")
         print(f"  Total Time:       {timing['total_time']:.3f} seconds")
         print(f"  Performance:      {perf['steps_per_second']:.0f} steps/second")
         print(f"  Throughput:       {perf['ns_per_day']:.1f} ns/day")
         print(f"  Energy Change:    {cpu_results['initial_energy']:.1f} → {cpu_results['final_energy']:.1f} kJ/mol")
      else:
         print(f"✗ CPU benchmark failed: {cpu_results.get('error', 'Unknown error')}")
      
      print(f"\n{'='*30} GPU RESULTS {'='*30}")
      if gpu_results["success"]:
         print(f"✓ Platform: {gpu_results['platform_used']}")
         timing = gpu_results["timing"]
         perf = gpu_results["performance_metrics"]
         print(f"  Minimization:     {timing['minimization_time']:.3f} seconds")
         print(f"  Equilibration:    {timing['equilibration_time']:.3f} seconds")
         print(f"  Production MD:    {timing['dynamics_time']:.3f} seconds")
         print(f"  Total Time:       {timing['total_time']:.3f} seconds")
         print(f"  Performance:      {perf['steps_per_second']:.0f} steps/second")
         print(f"  Throughput:       {perf['ns_per_day']:.1f} ns/day")
         print(f"  Energy Change:    {gpu_results['initial_energy']:.1f} → {gpu_results['final_energy']:.1f} kJ/mol")
      else:
         print(f"✗ GPU benchmark failed: {gpu_results.get('error', 'Unknown error')}")
      
      # Performance comparison
      print(f"\n{'='*25} PERFORMANCE COMPARISON {'='*25}")
      if cpu_results["success"] and gpu_results["success"]:
         cpu_time = cpu_results["timing"]["dynamics_time"]
         gpu_time = gpu_results["timing"]["dynamics_time"]
         cpu_perf = cpu_results["performance_metrics"]["steps_per_second"]
         gpu_perf = gpu_results["performance_metrics"]["steps_per_second"]
         
         print(f"CPU MD Time:          {cpu_time:.3f} seconds")
         print(f"GPU MD Time:          {gpu_time:.3f} seconds")
         print(f"Time Improvement:     {cpu_time/gpu_time:.2f}x faster on GPU")
         print(f"")
         print(f"CPU Performance:      {cpu_perf:.0f} steps/second")
         print(f"GPU Performance:      {gpu_perf:.0f} steps/second")
         print(f"Performance Gain:     {gpu_perf/cpu_perf:.2f}x faster on GPU")
         
         if speedup and speedup > 1.5:
            print(f"\n🚀 EXCELLENT: GPU shows {speedup:.1f}x speedup over CPU!")
            print(f"   Aurora Intel GPU acceleration is working effectively.")
         elif speedup and speedup > 1.0:
            print(f"\n✓ GOOD: GPU shows {speedup:.1f}x speedup over CPU")
            print(f"   GPU acceleration is working, though modest for this system size.")
         else:
            print(f"\n⚠️  WARNING: GPU speedup is only {speedup:.1f}x")
            print(f"   GPU may not be optimal for this system size or configuration.")
      elif cpu_results["success"]:
         cpu_perf = cpu_results["performance_metrics"]["steps_per_second"]
         cpu_throughput = cpu_results["performance_metrics"]["ns_per_day"]
         print(f"CPU-only benchmark completed:")
         print(f"  Performance:      {cpu_perf:.0f} steps/second")
         print(f"  Throughput:       {cpu_throughput:.1f} ns/day")
         print(f"\n📊 INFO: GPU comparison not available")
         print(f"   To enable GPU testing, ensure OpenCL platform is available.")
         print(f"   On Aurora: module load PrgEnv-gnu rocm")
      else:
         print("Cannot show performance results - CPU benchmark failed")
         
   else:
      print(f"\n✗ Benchmark Status: FAILED")
      print(f"✗ Error: {benchmark_results.get('error', 'Unknown error')}")
   
   print("\n" + "="*80)


def parse_cli():
   """Parse command line arguments"""
   parser = argparse.ArgumentParser(
      prog="test_openmm_aurora",
      description="Test OpenMM functionality and GPU detection on Aurora"
   )
   parser.add_argument(
      "--platform", 
      default="auto",
      choices=["auto", "OpenCL", "CPU", "CUDA"],
      help="OpenMM platform to use (default: auto)"
   )
   parser.add_argument(
      "--steps", 
      type=int, 
      default=1000,
      help="Number of MD steps to run (default: 1000)"
   )
   parser.add_argument(
      "--benchmark", 
      action="store_true",
      help="Run CPU vs GPU performance benchmark instead of basic test"
   )
   parser.add_argument(
      "--particles", 
      type=int, 
      default=2000,
      help="Number of particles for benchmark (default: 2000)"
   )
   parser.add_argument(
      "--benchmark-steps", 
      type=int, 
      default=5000,
      help="Number of MD steps for benchmark (default: 5000)"
   )
   parser.add_argument(
      "--log-level", 
      default="INFO",
      choices=["DEBUG", "INFO", "WARNING", "ERROR"],
      help="Logging level (default: INFO)"
   )
   
   return parser.parse_args()


def main():
   """Main function for standalone execution"""
   args = parse_cli()
   
   # Configure logging
   logging.getLogger().setLevel(getattr(logging, args.log_level))
   
   if args.benchmark:
      logger.info("Starting OpenMM CPU vs GPU performance benchmark...")
   else:
      logger.info("Starting OpenMM Aurora test...")
   
   # Test imports
   if not test_openmm_imports():
      sys.exit(1)
   
   # Detect platforms
   platform_info = detect_openmm_platforms()
   if not platform_info["available_platforms"]:
      logger.error("No OpenMM platforms available")
      sys.exit(1)
   
   if args.benchmark:
      # Check that both CPU and OpenCL platforms are available
      if "CPU" not in platform_info["available_platforms"]:
         logger.error("CPU platform not available for benchmark")
         sys.exit(1)
      if "OpenCL" not in platform_info["available_platforms"]:
         logger.warning("OpenCL platform not available - benchmark will only show CPU results")
      
      # Run performance benchmark
      benchmark_results = run_performance_benchmark(args.particles, args.benchmark_steps)
      
      # Generate benchmark report
      report_benchmark_results(benchmark_results, platform_info)
      
      # Exit with appropriate code
      sys.exit(0 if benchmark_results["success"] else 1)
   else:
      # Run basic simulation test
      results = run_openmm_test(args.platform, args.steps)
      
      # Generate basic report
      report_results(results, platform_info)
      
      # Exit with appropriate code
      sys.exit(0 if results["success"] else 1)


# Pytest compatibility
class TestOpenMM:
   """Pytest test class for OpenMM functionality"""
   
   def test_openmm_import(self):
      """Test that OpenMM can be imported"""
      assert test_openmm_imports(), "OpenMM import failed"
   
   def test_platform_detection(self):
      """Test platform detection"""
      platform_info = detect_openmm_platforms()
      assert len(platform_info["available_platforms"]) > 0, "No platforms detected"
      assert "CPU" in platform_info["available_platforms"], "CPU platform not available"
   
   def test_minimal_simulation(self):
      """Test minimal simulation runs successfully"""
      results = run_openmm_test("auto", 100)  # Short test for CI
      assert results["success"], f"Simulation failed: {results.get('error')}"
      assert results["steps_completed"] == 100, "Wrong number of steps completed"
   
   def test_opencl_platform_if_available(self):
      """Test OpenCL platform if available (for Aurora)"""
      platform_info = detect_openmm_platforms()
      if "OpenCL" in platform_info["available_platforms"]:
         results = run_openmm_test("OpenCL", 100)
         assert results["success"], f"OpenCL simulation failed: {results.get('error')}"
         assert results["platform_used"] == "OpenCL", "OpenCL platform not used"
   
   def test_benchmark_system_creation(self):
      """Test that benchmark system can be created"""
      system, topology = create_benchmark_system(100)  # Small system for testing
      assert system is not None, "System creation failed"
      assert topology is not None, "Topology creation failed"
      assert system.getNumParticles() == 100, "Wrong number of particles"
   
   def test_performance_benchmark_short(self):
      """Test performance benchmark with small system (for CI)"""
      benchmark_results = run_performance_benchmark(100, 500)  # Small system, short run
      # At minimum, CPU should work
      assert benchmark_results["cpu_results"]["success"], f"CPU benchmark failed: {benchmark_results['cpu_results'].get('error')}"


if __name__ == "__main__":
   main()
