#!/usr/bin/env python
"""
Generate Globus Compute Endpoint configuration for Aurora or Polaris systems.
Auto-detects the system based on hostname and uses appropriate defaults.
Supports customizable options with sensible defaults for both systems.
"""

import argparse
import logging
import os
import socket
import sys
from pathlib import Path
import yaml
from yaml import Loader, Dumper


def represent_literal_str(dumper, data):
   """Custom representer for literal block scalars."""
   if '\n' in data:
      return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')
   return dumper.represent_scalar('tag:yaml.org,2002:str', data)


# Add custom representer for multi-line strings
yaml.add_representer(str, represent_literal_str)


def detect_system():
   """Detect if running on Aurora or Polaris based on hostname."""
   hostname = socket.gethostname().lower()
   
   if 'aurora' in hostname:
      return 'aurora'
   elif 'polaris' in hostname:
      return 'polaris'
   else:
      return None


def parse_cli():
   """Parse command line arguments."""
   parser = argparse.ArgumentParser(
      prog="gen-endpoint-config",
      description="Generate Globus Compute Endpoint configuration for Aurora or Polaris (auto-detected)"
   )
   
   # Output options
   parser.add_argument(
      "-o", "--output",
      default="endpoint_config.yaml",
      help="Output configuration file (default: endpoint_config.yaml)"
   )
   
   # Engine configuration
   parser.add_argument(
      "--accelerators",
      type=int,
      default=12,
      help="Number of available accelerators (default: 12, overridden to 4 for Polaris)"
   )
   parser.add_argument(
      "--max-workers",
      type=int,
      default=12,
      help="Maximum workers per node (default: 12, overridden to 4 for Polaris)"
   )
   parser.add_argument(
      "--prefetch-capacity",
      type=int,
      default=24,
      help="Prefetch capacity (default: 24)"
   )
   parser.add_argument(
      "--max-retries",
      type=int,
      default=2,
      help="Max retries on system failure (default: 2)"
   )
   
   # Provider configuration
   parser.add_argument(
      "--account",
      default="datascience",
      help="PBS account (default: datascience)"
   )
   parser.add_argument(
      "--queue",
      default="debug",
      help="PBS queue (default: debug)"
   )
   parser.add_argument(
      "--walltime",
      default="01:00:00",
      help="Job walltime (default: 01:00:00)"
   )
   parser.add_argument(
      "--cpus-per-node",
      type=int,
      default=64,
      help="CPUs per node (default: 64)"
   )
   parser.add_argument(
      "--nodes-per-block",
      type=int,
      default=1,
      help="Nodes per block (default: 1)"
   )
   parser.add_argument(
      "--max-blocks",
      type=int,
      default=1,
      help="Maximum blocks (default: 1)"
   )
   parser.add_argument(
      "--min-blocks",
      type=int,
      default=0,
      help="Minimum blocks (default: 0)"
   )
   parser.add_argument(
      "--init-blocks",
      type=int,
      default=0,
      help="Initial blocks (default: 0)"
   )
   
   # Worker init customization - simplified
   parser.add_argument(
      "--venv-path",
      help="This should be the path to the virtual environment you setup for the globus compute endpoint. "
           "For example: '/lus/flare/projects/datascience/username/my_project/venv/bin/activate'. "
           "This is the full path to the activate script in your virtual environment's bin directory. "
           "Required for Aurora, ignored for Polaris (uses hardcoded paths)."
   )
   parser.add_argument(
      "--repo-path", 
      help="This is the path in which the repository code is located that should be added to PYTHONPATH. "
           "For example: '/lus/flare/projects/datascience/username/my_project/src'. "
           "This allows the compute workers to import your project modules and dependencies. "
           "Required for Aurora, ignored for Polaris (uses hardcoded paths)."
   )
   parser.add_argument(
      "--custom-worker-init",
      help="Path to custom worker_init script file to use instead of generated one"
   )
   
   # Logging
   parser.add_argument(
      "--log-level",
      default="INFO",
      choices=["DEBUG", "INFO", "WARNING", "ERROR"],
      help="Logging level (default: INFO)"
   )
   
   return parser.parse_args()


def config_logging(level):
   """Configure logging with specified level."""
   logging.basicConfig(
      level=getattr(logging, level),
      format="%(asctime)s | %(levelname)-8s | %(message)s",
      datefmt="%d-%m %H:%M"
   )


def generate_worker_init(args, system=None):
   """Generate worker_init script based on arguments and detected system."""
   if args.custom_worker_init:
      try:
         with open(args.custom_worker_init, 'r') as f:
            return f.read().strip()
      except FileNotFoundError:
         logging.error(f"Custom worker init file not found: {args.custom_worker_init}")
         sys.exit(1)
   
   # Auto-detect system if not provided
   if system is None:
      system = detect_system()
   
   if system == 'polaris':
      worker_init = """module use /soft/modulefiles
module load conda
source /home/parton/polaris/alcf-agentics-workflow/venv/bin/activate
cd /home/parton/polaris/alcf-agentics-workflow/remoteGlobusToAurora

export PYTHONPATH=/home/parton/polaris/alcf-agentics-workflow/remoteGlobusToAurora/src:$PYTHONPATH
export TMPDIR=/tmp"""
   elif system == 'aurora':
      worker_init = f"""cd {args.repo_path}
source {args.venv_path}/bin/activate
export PYTHONPATH={args.repo_path}/remoteGlobusToAurora/src:$PYTHONPATH"""
   else:
      # System not detected
      logging.error("Could not detect Aurora or Polaris system from hostname")
      logging.error(f"Hostname: {socket.gethostname()}")
      logging.error("Please use --custom-worker-init to provide a custom worker initialization script")
      sys.exit(1)
   
   return worker_init


def generate_config(args):
   """Generate the full endpoint configuration."""
   # Detect system and generate appropriate worker_init
   system = detect_system()
   worker_init = generate_worker_init(args, system)
   
   # Apply system-specific settings
   if system == 'polaris':
      # Polaris-specific settings
      accelerators = 4
      max_workers = 4
      scheduler_options = "#PBS -l filesystems=home:eagle"
   else:
      # Aurora and default settings
      accelerators = args.accelerators
      max_workers = args.max_workers
      scheduler_options = "#PBS -l filesystems=home:flare"
   
   config = {
      "engine": {
         "available_accelerators": accelerators,
         "max_retries_on_system_failure": args.max_retries,
         "max_workers_per_node": max_workers,
         "prefetch_capacity": args.prefetch_capacity,
         "provider": {
            "account": args.account,
            "cpus_per_node": args.cpus_per_node,
            "init_blocks": args.init_blocks,
            "launcher": {
               "bind_cmd": "--cpu-bind",
               "overrides": "--ppn 1",
               "type": "MpiExecLauncher"
            },
            "max_blocks": args.max_blocks,
            "min_blocks": args.min_blocks,
            "nodes_per_block": args.nodes_per_block,
            "queue": args.queue,
            "scheduler_options": scheduler_options,
            "select_options": "ngpus=1",
            "type": "PBSProProvider",
            "walltime": args.walltime,
            "worker_init": worker_init
         },
         "type": "GlobusComputeEngine"
      }
   }
   
   return config, system


def main():
   """Main function."""
   args = parse_cli()
   config_logging(args.log_level)
   
   # Detect system first for better logging
   detected_system = detect_system()
   if detected_system:
      logging.info(f"Detected system: {detected_system.upper()}")
      logging.info(f"Generating Globus Compute Endpoint configuration for {detected_system.upper()}")
   else:
      logging.info("Generating Globus Compute Endpoint configuration")
   
   # Validate required arguments for Aurora
   if detected_system == 'aurora' and not args.custom_worker_init:
      if not args.venv_path:
         logging.error("--venv-path is required for Aurora system")
         sys.exit(1)
      if not args.repo_path:
         logging.error("--repo-path is required for Aurora system")
         sys.exit(1)
   
   # Generate configuration
   config, system = generate_config(args)
   
   # Write to output file
   output_path = Path(args.output)
   logging.info(f"Writing configuration to {output_path}")
   
   try:
      with open(output_path, 'w') as f:
         yaml.dump(config, f, indent=2)
      
      logging.info(f"Configuration successfully written to {output_path}")
      
      # Log key settings
      logging.info(f"System: {system}, Account: {args.account}, Queue: {args.queue}")
      logging.info(f"Walltime: {args.walltime}")
      
      # Log system-specific settings
      if system == 'polaris':
         logging.info("Using Polaris-specific settings:")
         logging.info("  - Accelerators: 4, Max workers: 4")
         logging.info("  - Scheduler options: #PBS -l filesystems=home:eagle")
         logging.info("  - Using hardcoded Polaris paths")
      elif system == 'aurora':
         logging.info("Using Aurora-specific settings:")
         logging.info(f"  - Accelerators: {args.accelerators}, Max workers: {args.max_workers}")
         logging.info("  - Scheduler options: #PBS -l filesystems=home:flare")
         logging.info(f"  - Venv path: {args.venv_path}")
         logging.info(f"  - Repo path: {args.repo_path}")
      else:
         logging.info(f"Max workers: {args.max_workers}")
      
   except Exception as e:
      logging.error(f"Failed to write configuration: {e}")
      sys.exit(1)


if __name__ == "__main__":
   main()
