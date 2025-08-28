#!/usr/bin/env python
"""
Utility to confirm Globus Auth freshness and Aurora endpoint status.
3‑space indents per style.
"""
import argparse
import datetime
import json
import logging
import os
import shutil
import subprocess
import sys

# Add globus_compute_sdk import for endpoint operations
from globus_compute_sdk import Client

# Import shared Globus interface
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
from tools.globus_interface import check_auth_status, get_auth_object

def run(cmd):
   """Execute command and return stdout or raise on error"""
   ret = subprocess.run(cmd, capture_output=True, text=True)
   if ret.returncode != 0:
      raise RuntimeError(ret.stderr.strip())
   return ret.stdout


def check_auth(max_age_days=30):
   """Check Globus authentication token freshness using shared interface"""
   logging.info("Checking Globus authentication…")
   return check_auth_status(max_age_days)


def check_endpoint(eid):
   """Check Globus Compute endpoint status using SDK"""
   if not eid:
      logging.error("❌ GC_ENDPOINT_ID not set")
      logging.info("💡 Set environment variable: export GC_ENDPOINT_ID=your-endpoint-id")
      return False
      
   logging.info("Checking endpoint %s …", eid)
   
   try:
      # Create Globus Compute client
      compute_client = Client()
      
      # Get endpoint status using SDK
      status_info = compute_client.get_endpoint_status(eid)
      
      # Check if endpoint is available/online
      if status_info and status_info.get('status') in ['online', 'active', 'available']:
         logging.info("✅ Endpoint is online and available")
         return True
      else:
         status = status_info.get('status', 'unknown') if status_info else 'unknown'
         logging.warning(f"⚠️  Endpoint status: {status}")
         logging.debug(f"Full status info: {status_info}")
         return False
         
   except Exception as e:
      logging.error(f"❌ Endpoint check failed: {e}")
      logging.info("💡 Try: globus-compute-endpoint start %s", eid)
      return False


def check_python_env():
   """Check if required Python packages are installed"""
   logging.info("Checking Python environment…")
   
   required_packages = [
      "globus_compute_sdk",
      "globus_sdk",
      "langgraph", 
      "langchain_openai",
      "openmm"
   ]
   
   missing = []
   for package in required_packages:
      try:
         __import__(package)
      except ImportError:
         missing.append(package)
   
   if missing:
      logging.warning("⚠️  Missing packages: %s", ", ".join(missing))
      logging.info("💡 Install: pip install -r requirements.txt")
      return False
   else:
      logging.info("✅ All required packages available")
      return True


def main():
   """Main entry point"""
   parser = argparse.ArgumentParser(
      description="Check Globus Auth and endpoint status"
   )
   parser.add_argument("--endpoint", "-e",
                       default=os.getenv("GC_ENDPOINT_ID"),
                       help="Globus Compute endpoint ID")
   parser.add_argument("--max-age", "-a", type=int, default=30,
                       help="Maximum token age in days")
   parser.add_argument("--verbose", "-v", action="store_true",
                       help="Enable verbose logging")
   
   args = parser.parse_args()
   
   # Configure logging
   level = logging.DEBUG if args.verbose else logging.INFO
   logging.basicConfig(
      level=level, 
      format="%(levelname)s | %(message)s"
   )
   
   logging.info("🔍 Running Globus environment checks...")
   
   # Run all checks
   checks = [
      ("Python Environment", check_python_env),
      ("Globus Auth", lambda: check_auth(args.max_age)),
      ("Compute Endpoint", lambda: check_endpoint(args.endpoint))
   ]
   
   passed = 0
   total = len(checks)
   
   for name, check_func in checks:
      logging.info(f"\n--- {name} ---")
      if check_func():
         passed += 1
   
   # Summary
   logging.info(f"\n📊 Summary: {passed}/{total} checks passed")
   
   if passed == total:
      logging.info("🎉 All checks passed - ready to run agentic workflow!")
      sys.exit(0)
   else:
      logging.error("❌ Some checks failed - please fix issues above")
      sys.exit(1)


if __name__ == "__main__":
   main() 