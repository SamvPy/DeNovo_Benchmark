#!/bin/bash
set -e  # stop on first error
eval "$(conda shell.bash hook)"  # initialize conda in non-interactive shell

# Resolve script directory to make relative paths consistent
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="$SCRIPT_DIR/environments"
SRC_DIR="$SCRIPT_DIR/source_code"

cd "$ENV_DIR" || exit

echo "==============================================="
echo " Setting up all de novo conda environments"
echo "==============================================="

####################################
# Helper function to set up an env #
####################################
setup_env() {
    local env_name=$1
    local yaml_file=$2
    local src_subdir=$3
    local pip_pkg=$4

    echo "-----------------------------------------------"
    echo " Checking environment: $env_name"
    echo "-----------------------------------------------"

    if [ ! -f "$yaml_file" ]; then
        echo "❌ Environment file '$yaml_file' not found. Skipping $env_name setup."
        return
    fi

    # Create env only if it doesn't exist
    if conda env list | grep -q "$env_name"; then
        echo "✔ Environment '$env_name' already exists. Skipping creation."
    else
        echo "Creating environment '$env_name'..."
        conda env create -f "$yaml_file"
    fi

    echo "Activating environment '$env_name'..."
    conda activate "$env_name"

    # Only for pi-PrimeNovo: Install additional dependencies
    if [ "$env_name" == "piprimenovo_env" ]; then
        echo "Installing additional dependencies for pi-PrimeNovo..."
        CUDA_VER=$(python -c "import torch; print(torch.version.cuda)")
        if [[ $CUDA_VER == 12* ]]; then
            pip install cupy-cuda12x
        elif [[ $CUDA_VER == 11* ]]; then
            pip install cupy-cuda11x
        else
            pip install cupy
        fi
        git clone --recursive https://github.com/WayenVan/ctcdecode.git
        cd ctcdecode
        pip install .
        cd ..  #this is needed as ctcdecode can not be imported under the current directory
        rm -rf ctcdecode
    fi

    # Optional pip install from source if provided
    if [ -n "$src_subdir" ]; then
        if [ -d "$SRC_DIR/$src_subdir" ]; then
            echo "Installing $env_name from source: $SRC_DIR/$src_subdir"
            cd "$SRC_DIR/$src_subdir" || exit
            pip install .
            cd "$ENV_DIR" || exit
        else
            echo "⚠ Source directory '$SRC_DIR/$src_subdir' not found. Skipping pip install."
        fi
    fi

    # Optional pip package if specified
    if [ -n "$pip_pkg" ]; then
        echo "Installing pip package(s): $pip_pkg"
        pip install $pip_pkg
    fi

    conda deactivate
    cd "$ENV_DIR" || exit
    echo "✅ Finished setting up $env_name"
    echo
}

########################################
### Install denovo_utils environment ###
########################################
echo "Setting up base denovo_utils environment..."
setup_env "denovo_analysis_env" "denovo_utils.yaml" "ms2rescore" "../package_du ../package_peaks"


#########################################
### INSTALL DE NOVO TOOL ENVIRONMENTS ###
#########################################
setup_env "adanovo_env"      "adanovo.yaml"      "adanovo_v1"
setup_env "casanovo_env"     "casanovo.yaml"     ""
setup_env "contranovo_env"   "contranovo.yaml"   "ContraNovo"
setup_env "instanovo_env"    "instanovo.yaml"    ""
setup_env "novob_env"        "novob.yaml"        "NovoB"
setup_env "pepnet_env"       "pepnet.yaml"       "PepNet"
setup_env "pihelixnovo_env"  "pihelixnovo.yaml"  "pi-HelixNovo"
setup_env "piprimenovo_env"  "piprimenovo.yaml"  "pi-PrimeNovo"
setup_env "spectralis_env"   "spectralis.yaml"   "spectralis"

echo "==============================================="
echo "✅ All environments processed successfully."
echo "==============================================="
