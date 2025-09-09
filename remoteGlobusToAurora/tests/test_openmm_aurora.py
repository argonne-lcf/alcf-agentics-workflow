#!/usr/bin/env python
"""
OpenMM test for Aurora GPU detection and functionality.
Can be run standalone or as part of pytest suite.

Usage:
   python test_openmm_aurora.py --platform auto
   python test_openmm_aurora.py --platform OpenCL
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
      topology.addAtom("H1", app.Element.getBySymbol("H"), residue)
      topology.addAtom("H2", app.Element.getBySymbol("H"), residue)
      topology.addBond(topology.atoms()[0], topology.atoms()[1])
      
      logger.info("✓ Minimal system created successfully")
      return system, topology
      
   except Exception as e:
      logger.error(f"✗ System creation failed: {e}")
      raise


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
      
      logger.info(f"Initial energy: {initial_energy}")
      
      # Run dynamics
      logger.info(f"Running {steps} MD steps...")
      start_time = time.time()
      simulation.step(steps)
      dynamics_time = time.time() - start_time
      
      # Get final state
      final_state = simulation.context.getState(getEnergy=True, getPositions=True)
      final_energy = final_state.getPotentialEnergy()
      final_positions = final_state.getPositions()
      
      logger.info(f"Final energy: {final_energy}")
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
   
   logger.info("Starting OpenMM Aurora test...")
   
   # Test imports
   if not test_openmm_imports():
      sys.exit(1)
   
   # Detect platforms
   platform_info = detect_openmm_platforms()
   if not platform_info["available_platforms"]:
      logger.error("No OpenMM platforms available")
      sys.exit(1)
   
   # Run simulation test
   results = run_openmm_test(args.platform, args.steps)
   
   # Generate report
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


if __name__ == "__main__":
   main()
