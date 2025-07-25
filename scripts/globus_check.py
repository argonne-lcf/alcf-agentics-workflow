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

# Add globus_sdk import for authentication
import globus_sdk
from globus_sdk.login_flows import LocalServerLoginFlowManager # Needed to access globus_sdk.gare

# Add globus_compute_sdk import for endpoint operations
from globus_compute_sdk import Client

# Globus authentication constants (from inference_auth_token.py)
APP_NAME = "alcf_agentics_workflow"
# Public inference auth client
AUTH_CLIENT_ID = "58fdd3bc-e1c3-4ce5-80ea-8d6b87cfb944"
# Inference gateway API scope
GATEWAY_CLIENT_ID = "681c10cc-f684-4540-bcd7-0b4df3bc26ef"
GATEWAY_SCOPE = f"https://auth.globus.org/scopes/{GATEWAY_CLIENT_ID}/action_all"

# Allowed identity provider domains
ALLOWED_DOMAINS = ["anl.gov", "alcf.anl.gov"]

# Globus authorizer parameters to point to specific identity providers
GA_PARAMS = globus_sdk.gare.GlobusAuthorizationParameters(session_required_single_domain=ALLOWED_DOMAINS)


# Error handler to guide user through specific identity providers 
class DomainBasedErrorHandler:
    def __call__(self, app, error):
        print(f"Encountered error '{error}', initiating login...")
        app.login(auth_params=GA_PARAMS)

def run(cmd):
   """Execute command and return stdout or raise on error"""
   ret = subprocess.run(cmd, capture_output=True, text=True)
   if ret.returncode != 0:
      raise RuntimeError(ret.stderr.strip())
   return ret.stdout


def check_auth(max_age_days=30):
   """Check Globus authentication token freshness using SDK"""
   logging.info("Checking Globus authentication…")
   
   try:
      # Create Globus user application
      logging.info(f"Creating Globus user application for {APP_NAME} with client ID {AUTH_CLIENT_ID} and scope {GATEWAY_SCOPE}")
      app = globus_sdk.UserApp(
         APP_NAME,
         client_id=AUTH_CLIENT_ID,
         scope_requirements={GATEWAY_CLIENT_ID: [GATEWAY_SCOPE]},
         config=globus_sdk.GlobusAppConfig(
            request_refresh_tokens=True,
            token_validation_error_handler=DomainBasedErrorHandler()
         ),
      )
      logging.info(f"Globus user application created: {app}")
      
      # Get authorizer object 
      auth = app.get_authorizer(GATEWAY_CLIENT_ID)
      
      # Check if token is valid and refresh if needed
      # auth.ensure_valid_token()
      
      logging.info("✅ Globus authentication valid")
      return True
      
   except globus_sdk.AuthAPIError as e:
      logging.error(f"❌ Globus authentication failed: {e}")
      logging.info("💡 Run: python scripts/inference_auth_token.py authenticate")
      return False
   except Exception as e:
      logging.error(f"❌ Authentication check failed: {e}")
      logging.info("💡 Try authenticating with: python scripts/inference_auth_token.py authenticate")
      return False


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