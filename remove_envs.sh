#!/usr/bin/env bash
set -euo pipefail

# Initialize conda for non-interactive shells
eval "$(conda shell.bash hook)"

echo "==============================================="
echo " Removing all de novo conda environments"
echo "==============================================="

# List of environments to remove
envs=(
  denovo_analysis_env
  adanovo_env
  casanovo_env
  contranovo_env
  instanovo_env
  denovo_utils_env
  novob_env
  pepnet_env
  pihelixnovo_env
  piprimenovo_env
  spectralis_env
)

for env_name in "${envs[@]}"; do
  if conda env list | awk '{print $1}' | grep -xq "$env_name"; then
    echo "Removing environment: $env_name ..."
    conda env remove -n "$env_name" -y || {
      echo "⚠ Failed to remove environment: $env_name (continuing)"
    }
  else
    echo "Environment $env_name does not exist. Skipping."
  fi
done

echo "==============================================="
echo "✅ All de novo conda environments processed for removal."
echo "==============================================="

exit 0