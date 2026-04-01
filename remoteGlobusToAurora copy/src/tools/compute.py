"""
Globus Compute wrapper for submitting MD simulations to Aurora.
Handles job submission, monitoring, and result retrieval.
"""

import logging
import os
import time
from typing import Dict, Any, Optional

from globus_compute_sdk import Client, Executor


def get_logger():
   """Get the module logger"""
   return logging.getLogger("agentic_demo.compute")


class GlobusComputeWrapper:
   """Wrapper for Globus Compute operations"""
   
   def __init__(self, endpoint_id: Optional[str] = None):
      """Initialize the Globus Compute client
      
      Args:
         endpoint_id: Globus Compute endpoint ID. If None, uses GC_ENDPOINT_ID env var.
      """
      self.endpoint_id = endpoint_id or os.getenv("GC_ENDPOINT_ID")
      if not self.endpoint_id:
         raise ValueError("Globus Compute endpoint ID not provided")
      
      self.client = Client()
      self.executor = Executor(endpoint_id=self.endpoint_id)
      
   def submit_simulation(self, params: Dict[str, Any], timeout: int = 180) -> Dict[str, Any]:
      """Submit a molecular dynamics simulation job
      
      Args:
         params: Simulation parameters including protein, timestep, temperature, steps
         timeout: Maximum time to wait for job completion in seconds
         
      Returns:
         Dictionary containing simulation results or error information
      """
      logger = get_logger()
      logger.info(f"Submitting simulation job to endpoint {self.endpoint_id}")
      logger.debug(f"Simulation parameters: {params}")
      
      try:
         # Import the simulation kernel function - handle both test and normal execution
         try:
            from ..sim_kernel import run_md_simulation
         except ImportError:
            # Fallback for test execution
            import sys
            import os
            sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
            from sim_kernel import run_md_simulation
         
         # Submit the function to Globus Compute
         future = self.executor.submit(run_md_simulation, params)
         
         logger.info(f"Job submitted with ID: {future.task_id}")
         
         # Wait for completion with timeout
         start_time = time.time()
         while not future.done():
            elapsed = time.time() - start_time
            if elapsed > timeout:
               logger.error(f"Job timed out after {timeout}s")
               return {
                  "status": "timeout",
                  "error": f"Job timed out after {timeout} seconds",
                  "task_id": future.task_id
               }
            
            logger.debug(f"Waiting for job completion... ({elapsed:.1f}s elapsed)")
            time.sleep(5)
         
         # Get the result
         result = future.result()
         elapsed = time.time() - start_time
         
         logger.info(f"Job completed successfully in {elapsed:.1f}s")
         return {
            "status": "completed",
            "task_id": future.task_id,
            "execution_time": elapsed,
            **result
         }
         
      except Exception as e:
         logger.error(f"Globus Compute job failed: {e}")
         return {
            "status": "error",
            "error": str(e)
         }
   
   def check_endpoint_status(self) -> Dict[str, Any]:
      """Check if the endpoint is online and accessible
      
      Returns:
         Dictionary with endpoint status information
      """
      logger = get_logger()
      try:
         # Try to get endpoint details
         endpoint_info = self.client.get_endpoint_metadata(self.endpoint_id)
         
         return {
            "status": "online",
            "endpoint_id": self.endpoint_id,
            "name": endpoint_info.get("name", "Unknown"),
            "description": endpoint_info.get("description", "")
         }
         
      except Exception as e:
         logger.error(f"Failed to check endpoint status: {e}")
         return {
            "status": "error",
            "endpoint_id": self.endpoint_id,
            "error": str(e)
         } 