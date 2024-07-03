#!/bin/bash

# Loop through each environment directory
for env_dir in source_code/*; do
  if [ -d "$env_dir" ]; then
    env_name=$(basename "$env_dir")

    # Skip the casanovo directory
    if [ "$env_name" == "casanovo" ]; then
      echo "Skipping $env_dir..."
      continue
    fi

    echo "Processing $env_dir..."

    # Navigate to the environment directory
    cd "$env_dir" || exit

    # Extract the environment name from the directory name
    env_name=$(basename "$env_dir")

    # Create the conda environment
    echo "Creating conda environment $env_name..."
    conda env create -f environment.yml

    # Activate the environment
    echo "Activating environment $env_name..."
    source activate "$env_name"

    # Install the local package using pip
    echo "Installing local package for $env_name..."
    pip install .

    # Deactivate the environment
    echo "Deactivating environment $env_name..."
    conda deactivate

    # Navigate back to the parent directory
    cd - || exit

    echo "Finished processing $env_dir."
  fi
done

echo "All environments have been processed."
