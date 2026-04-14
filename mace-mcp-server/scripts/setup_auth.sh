#!/usr/bin/env bash
# -------------------------------------------------------------------
# ALCF Inference Endpoint Authentication
#
# Downloads the ALCF authentication helper (if not already present),
# runs Globus-based authentication, and exports the access token.
#
# Usage (run from the mace-mcp-server directory):
#
#   # Authenticate and set the token for this shell session
#   source scripts/setup_auth.sh
#
#   # Authenticate and persist the token in your shell profile
#   # (~/.zshrc or ~/.bashrc, auto-detected)
#   source scripts/setup_auth.sh --persist
#
#   # Force re-authentication (e.g. after 30-day expiry)
#   source scripts/setup_auth.sh --force
#
#   # Combine flags
#   source scripts/setup_auth.sh --force --persist
# -------------------------------------------------------------------

set -euo pipefail

AUTH_SCRIPT="inference_auth_token.py"
AUTH_URL="https://raw.githubusercontent.com/argonne-lcf/inference-endpoints/refs/heads/main/${AUTH_SCRIPT}"

# ----- Parse flags -----

FORCE=""
PERSIST=false

for arg in "$@"; do
    case "$arg" in
        --force)  FORCE="--force" ;;
        --persist) PERSIST=true ;;
        *)
            echo "Unknown option: $arg"
            echo "Usage: source scripts/setup_auth.sh [--force] [--persist]"
            return 1 2>/dev/null || exit 1
            ;;
    esac
done

# ----- Helper: mask a token for safe display -----

mask_token() {
    local token="$1"
    local last4="${token: -4}"
    echo "****${last4}"
}

# ----- Step 1: Download the auth helper if missing -----

if [[ ! -f "$AUTH_SCRIPT" ]]; then
    echo "Downloading ${AUTH_SCRIPT} ..."
    curl -fsSL "$AUTH_URL" -o "$AUTH_SCRIPT"
    echo "  saved to ./${AUTH_SCRIPT}"
else
    echo "Auth helper already present: ./${AUTH_SCRIPT}"
fi

# ----- Step 2: Authenticate with Globus / ALCF -----

echo ""
echo "Authenticating with ALCF (this may open a browser window) ..."
python "$AUTH_SCRIPT" authenticate $FORCE

# ----- Step 3: Obtain and export the access token -----

TOKEN=$(python "$AUTH_SCRIPT" get_access_token)

if [[ -z "$TOKEN" ]]; then
    echo ""
    echo "ERROR: Failed to obtain an access token."
    echo "Try re-running with --force:  source scripts/setup_auth.sh --force"
    return 1 2>/dev/null || exit 1
fi

export ALCF_INFERENCE_TOKEN="$TOKEN"

MASKED=$(mask_token "$TOKEN")
echo ""
echo "ALCF_INFERENCE_TOKEN has been set (${MASKED})."

# ----- Step 4 (optional): Persist to shell RC file -----

if $PERSIST; then
    # Auto-detect RC file based on current shell
    if [[ "$SHELL" == *zsh* ]]; then
        RC_FILE="$HOME/.zshrc"
    else
        RC_FILE="$HOME/.bashrc"
    fi

    EXPORT_LINE="export ALCF_INFERENCE_TOKEN=\"${TOKEN}\""

    if [[ -f "$RC_FILE" ]] && grep -q "^export ALCF_INFERENCE_TOKEN=" "$RC_FILE"; then
        # Replace existing line
        # Use a delimiter that won't clash with the token value
        sed -i.bak "s|^export ALCF_INFERENCE_TOKEN=.*|${EXPORT_LINE}|" "$RC_FILE"
        rm -f "${RC_FILE}.bak"
        echo "Updated existing token in ${RC_FILE} (${MASKED})."
    else
        # Append
        echo "" >> "$RC_FILE"
        echo "# ALCF Inference Endpoint token (added by scripts/setup_auth.sh)" >> "$RC_FILE"
        echo "$EXPORT_LINE" >> "$RC_FILE"
        echo "Token saved to ${RC_FILE} (${MASKED}). It will be available in new shell sessions."
    fi
fi

# ----- Warn if not sourced -----

# shellcheck disable=SC2128
if [[ "${BASH_SOURCE[0]}" == "${0}" ]] 2>/dev/null; then
    echo ""
    echo "NOTE: You ran this script directly, so the token is only set in"
    echo "this subshell.  To set it in your current shell, run:"
    echo ""
    echo "  source scripts/setup_auth.sh --persist"
fi
