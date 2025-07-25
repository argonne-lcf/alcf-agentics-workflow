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


def run(cmd):
   """Execute command and return stdout or raise on error"""
   ret = subprocess.run(cmd, capture_output=True, text=True)
   if ret.returncode != 0:
      raise RuntimeError(ret.stderr.strip())
   return ret.stdout


def check_auth(max_age_days=30):
   """Check Globus authentication token freshness"""
   logging.info("Checking Globus session …")
   
   try:
      data = json.loads(run(["globus", "session", "show", "--format", "json"]))
      oldest = max(datetime.datetime.fromisoformat(i["auth_time"])
                   for i in data["identities"])
      age = (datetime.datetime.utcnow() - oldest).days
      
      if age > max_age_days:
         logging.warning("⚠️  Tokens are %sd old; run `globus login`.", age)
         return False
      else:
         logging.info("Globus tokens fresh (%sd).", age)
         return True
         
   except subprocess.CalledProcessError:
      logging.error("❌ Globus CLI not found or not logged in")
      logging.info("💡 Run: globus login")
      return False
   except Exception as e:
      logging.error(f"❌ Auth check failed: {e}")
      return False


def check_endpoint(eid):
   """Check Globus Compute endpoint status"""
   if not eid:
      logging.error("❌ GC_ENDPOINT_ID not set")
      logging.info("💡 Set environment variable: export GC_ENDPOINT_ID=your-endpoint-id")
      return False
      
   if not shutil.which("globus-compute-endpoint"):
      logging.error("❌ globus-compute-endpoint CLI missing in PATH")
      logging.info("💡 Install: pip install globus-compute-sdk")
      return False
      
   logging.info("Checking endpoint %s …", eid)
   
   try:
      output = run(["globus-compute-endpoint", "status", eid])
      if "online" in output.lower() or "active" in output.lower():
         logging.info("✅ Endpoint reachable and active")
         return True
      else:
         logging.warning("⚠️  Endpoint may be offline")
         logging.debug(f"Status output: {output}")
         return False
         
   except RuntimeError as e:
      logging.error(f"❌ Endpoint check failed: {e}")
      logging.info("💡 Try: globus-compute-endpoint start %s", eid)
      return False


def check_python_env():
   """Check if required Python packages are installed"""
   logging.info("Checking Python environment…")
   
   required_packages = [
      "globus_compute_sdk",
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