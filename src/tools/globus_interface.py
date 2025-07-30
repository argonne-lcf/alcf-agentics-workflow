#!/usr/bin/env python
"""
Shared Globus authentication interface for ALCF Agentic Workflow.
Centralizes authentication logic for use across the application.
"""

import os
import logging
import globus_sdk
from globus_sdk.login_flows import LocalServerLoginFlowManager

# Globus UserApp name
APP_NAME = "alcf_agentics_workflow"

# Public inference auth client
AUTH_CLIENT_ID = "58fdd3bc-e1c3-4ce5-80ea-8d6b87cfb944"

# Inference gateway API scope
GATEWAY_CLIENT_ID = "681c10cc-f684-4540-bcd7-0b4df3bc26ef"
GATEWAY_SCOPE = f"https://auth.globus.org/scopes/{GATEWAY_CLIENT_ID}/action_all"

# Path where access and refresh tokens are stored
TOKENS_PATH = f"{os.path.expanduser('~')}/.globus/app/{AUTH_CLIENT_ID}/{APP_NAME}/tokens.json"

# Allowed identity provider domains
ALLOWED_DOMAINS = ["anl.gov", "alcf.anl.gov"]

# Globus authorizer parameters to point to specific identity providers
GA_PARAMS = globus_sdk.gare.GlobusAuthorizationParameters(session_required_single_domain=ALLOWED_DOMAINS)


def get_logger():
   """Get the module logger"""
   return logging.getLogger("agentic_demo.auth")


class DomainBasedErrorHandler:
   """Error handler to guide user through specific identity providers"""
   def __call__(self, app, error):
      logger = get_logger()
      logger.warning(f"Encountered error '{error}', initiating login...")
      app.login(auth_params=GA_PARAMS)


def get_auth_object(force=False):
   """
   Create a Globus UserApp with the inference service scope
   and trigger the authentication process. If authentication
   has already happened, existing tokens will be reused.
   """
   # Create Globus user application
   app = globus_sdk.UserApp(
      APP_NAME,
      client_id=AUTH_CLIENT_ID,
      scope_requirements={GATEWAY_CLIENT_ID: [GATEWAY_SCOPE]},
      config=globus_sdk.GlobusAppConfig(
         request_refresh_tokens=True,
         token_validation_error_handler=DomainBasedErrorHandler()
      ),
   )

   # Force re-login if required
   if force:
      app.login(auth_params=GA_PARAMS)

   # Authenticate using your Globus account or reuse existing tokens
   auth = app.get_authorizer(GATEWAY_CLIENT_ID)

   # Return the Globus refresh token authorizer
   return auth


def get_access_token():
   """
   Load existing tokens, refresh the access token if necessary,
   and return the valid access token. If there is no token stored
   in the home directory, or if the refresh token is expired following
   6 months of inactivity, an authentication will be triggered.
   """
   # Get authorizer object and authenticate if need be
   auth = get_auth_object()

   # Make sure the stored access token is valid, and refresh otherwise
   auth.ensure_valid_token()

   # Return the access token
   return auth.access_token


def check_auth_status(max_age_days=30):
   """
   Check if Globus authentication is valid and tokens are fresh.
   
   Args:
      max_age_days: Maximum age for tokens before warning
      
   Returns:
      bool: True if authentication is valid, False otherwise
   """
   logger = get_logger()
   try:
      # Get authorizer object
      auth = get_auth_object()
      
      # Check if token is valid
      auth.ensure_valid_token()
      
      logger.info("✅ Globus authentication valid")
      return True
      
   except globus_sdk.AuthAPIError as e:
      logger.error(f"❌ Globus authentication failed: {e}")
      logger.info("💡 Run authentication with: python -c 'from src.tools.globus_interface import get_auth_object; get_auth_object(force=True)'")
      return False
   except Exception as e:
      logger.error(f"❌ Authentication check failed: {e}")
      return False


def tokens_exist():
   """Check if authentication tokens exist on disk"""
   return os.path.isfile(TOKENS_PATH)


class GlobusAuthError(Exception):
   """Custom exception for Globus authentication errors"""
   pass


# CLI interface for standalone usage
if __name__ == "__main__":
   import argparse

   AUTHENTICATE_ACTION = "authenticate"
   GET_ACCESS_TOKEN_ACTION = "get_access_token"
   CHECK_STATUS_ACTION = "check_status"

   parser = argparse.ArgumentParser(description="Globus authentication interface")
   parser.add_argument('action', choices=[AUTHENTICATE_ACTION, GET_ACCESS_TOKEN_ACTION, CHECK_STATUS_ACTION])
   parser.add_argument("-f", "--force", action="store_true", help="authenticate from scratch")
   parser.add_argument("-v", "--verbose", action="store_true", help="verbose logging")
   args = parser.parse_args()

   # Configure logging for CLI usage - this is separate from the module logger
   level = logging.DEBUG if args.verbose else logging.INFO
   logging.basicConfig(level=level, format="%(levelname)s | %(message)s")

   # Authentication
   if args.action == AUTHENTICATE_ACTION:
      _ = get_auth_object(force=args.force)
      logging.info("✅ Authentication completed")

   # Get token
   elif args.action == GET_ACCESS_TOKEN_ACTION:
      if not tokens_exist():
         raise GlobusAuthError('Access token does not exist. '
            'Please authenticate by running "python src/tools/globus_interface.py authenticate".')
      
      if args.force:
         raise GlobusAuthError(f"The --force flag cannot be used with the {GET_ACCESS_TOKEN_ACTION} action.")

      print(get_access_token())

   # Check status
   elif args.action == CHECK_STATUS_ACTION:
      if check_auth_status():
         logging.info("🎉 Authentication status: OK")
      else:
         logging.error("❌ Authentication status: FAILED")
         exit(1) 