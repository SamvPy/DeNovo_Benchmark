#!/bin/bash

# Install the denovo_utils package in its own environment
cd ../package_du || exit
conda env create -f denovo_analysis_env.yaml
conda activate denovo_analysis_env
pip install .
conda deactivate
cd ../source_code || exit

# Loop through each environment directory
for env_dir in ./*; do
  if [ -d "$env_dir" ]; then
    env_name=$(basename "$env_dir")

    echo "Processing $env_dir..."

    # Navigate to the environment directory
    cd "$env_dir" || exit

    # Convert the directory name to lowercase to construct the environment file name
    env_file="${env_name,,}_env.yaml"

    # Check if the environment file exists
    if [ ! -f "$env_file" ]; then
      echo "Environment file $env_file does not exist. Skipping $env_dir..."
      cd - || exit
      continue
    fi

    # Extract the environment name from the directory name
    env_name=$(basename "$env_dir")

    # Create the conda environment
    echo "Creating conda environment $env_name..."
    conda env create -f "$env_file"

    # Activate the environment
    echo "Activating environment $env_name..."
    conda activate "${env_name,,}_env"

    # casanovo requires different installation procedure
    if [[ "$env_name" == "casanovo" ]]; then
      echo "Installing $env_dir..."
      pip install casanovo
      conda deactivate
      cd - || exit
      echo "Finished processing $env_dir."
      continue
    fi

    # spectralis requires different installation procedure
    if [[ "$env_name" == "spectralis" ]]; then
      echo "Installing $env_dir..."
      pip install .
      # Make sure the right cuda installation is performed
      pip uninstall torch
      pip install torch

      # Check for CUDA availability
      echo "Checking CUDA availability..."
      cuda_check=$(python -c "import torch; print(torch.cuda.is_available())")

      conda deactivate
      cd - || exit
      echo "Finished processing $env_dir."
      continue
    fi

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

if [[ "$cuda_check" == "False" ]]; then
  echo "Warning: GPU is not visible for the spectralis environment."
fi
