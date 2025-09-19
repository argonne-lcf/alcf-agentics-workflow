#!/usr/bin/env python
"""
Generate Globus Compute Endpoint configuration for Aurora system.
Supports customizable options with sensible defaults for Aurora.
"""

import argparse
import logging
import os
import sys
from pathlib import Path
import yaml


def parse_cli():
   """Parse command line arguments."""
   parser = argparse.ArgumentParser(
      prog="gen-endpoint-config",
      description="Generate Globus Compute Endpoint configuration for Aurora"
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
      help="Number of available accelerators (default: 12)"
   )
   parser.add_argument(
      "--max-workers",
      type=int,
      default=12,
      help="Maximum workers per node (default: 12)"
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
      required=True,
      help="This should be the path to the virtual environment you setup for the globus compute endpoint. "
           "For example: '/lus/flare/projects/datascience/username/my_project/venv/bin/activate'. "
           "This is the full path to the activate script in your virtual environment's bin directory."
   )
   parser.add_argument(
      "--repo-path", 
      required=True,
      help="This is the path in which the repository code is located that should be added to PYTHONPATH. "
           "For example: '/lus/flare/projects/datascience/username/my_project/src'. "
           "This allows the compute workers to import your project modules and dependencies."
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


def generate_worker_init(args):
   """Generate worker_init script based on arguments."""
   if args.custom_worker_init:
      try:
         with open(args.custom_worker_init, 'r') as f:
            return f.read().strip()
      except FileNotFoundError:
         logging.error(f"Custom worker init file not found: {args.custom_worker_init}")
         sys.exit(1)
   
   worker_init = f"""cd {args.repo_path}
source {args.venv_path}/bin/activate
export PYTHONPATH={args.repo_path}/remoteGlobusToAuror/src:$PYTHONPATH"""
   
   return worker_init


def generate_config(args):
   """Generate the full endpoint configuration."""
   worker_init = generate_worker_init(args)
   
   config = {
      "engine": {
         "available_accelerators": args.accelerators,
         "max_retries_on_system_failure": args.max_retries,
         "max_workers_per_node": args.max_workers,
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
            "scheduler_options": "#PBS -l filesystems=home:flare",
            "select_options": "ngpus=1",
            "type": "PBSProProvider",
            "walltime": args.walltime,
            "worker_init": worker_init
         },
         "type": "GlobusComputeEngine"
      }
   }
   
   return config


def main():
   """Main function."""
   args = parse_cli()
   config_logging(args.log_level)
   
   logging.info("Generating Globus Compute Endpoint configuration for Aurora")
   
   # Generate configuration
   config = generate_config(args)
   
   # Write to output file
   output_path = Path(args.output)
   logging.info(f"Writing configuration to {output_path}")
   
   try:
      with open(output_path, 'w') as f:
         yaml.dump(config, f, default_flow_style=False, indent=2)
      
      logging.info(f"Configuration successfully written to {output_path}")
      
      # Log key settings
      logging.info(f"Account: {args.account}, Queue: {args.queue}")
      logging.info(f"Walltime: {args.walltime}, Max workers: {args.max_workers}")
      logging.info(f"Venv path: {args.venv_path}")
      logging.info(f"Repo path: {args.repo_path}")
      
   except Exception as e:
      logging.error(f"Failed to write configuration: {e}")
      sys.exit(1)


if __name__ == "__main__":
   main()
