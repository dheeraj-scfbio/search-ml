#!/bin/bash

# ==============================================================================
# Search-ML: Environment Installation Script
# ==============================================================================
# Description:
#   This script automates the setup of the Search-ML virtual screening pipeline.
#   It detects the local Conda installation (Miniconda/Anaconda), validates
#   dependencies, and creates/updates the required Conda environment from the
#   'environment.yml' file.
#
# Usage:
#   chmod +x install.sh
#   ./install.sh
#
# Author:
#   Dheeraj Kumar Chaurasia
#   Supercomputing Facility for Bioinformatics & Computational Biology (SCFBio)
#   IIT Delhi
#
# Email:
#   dheeraj@scfbio-iitd.res.in
#
# Copyright:
#   Copyright (c) 2025 SCFBio, IIT Delhi. All rights reserved.
# ==============================================================================

set -e

# --- Configuration ---
ENV_FILE="environment.yml"
CONDA_ROOT=""

echo "=========================================="
echo "   🚀 Starting Search-ML Setup"
echo "=========================================="

# ---------------------------------------------------------
# FUNCTION: Detect Conda dynamically
# ---------------------------------------------------------
detect_conda() {
    # Method 1: Check 'which conda' and extract the root directory
    if command -v conda &> /dev/null; then
        echo "🔍 'conda' command found in PATH."
        
        # Get the full path to the executable (e.g., /home/user/miniconda3/bin/conda)
        CONDA_BIN=$(command -v conda)
        
        # Resolve symlinks (in case conda is a symlink in /usr/bin)
        # MacOS usually has readlink, Linux has readlink -f. We try python if readlink fails.
        if [[ "$OSTYPE" == "darwin"* ]]; then
             # MacOS doesn't always have -f, so we trust the path or use python to resolve
             REAL_PATH=$(python3 -c "import os; print(os.path.realpath('$CONDA_BIN'))")
        else
             REAL_PATH=$(readlink -f "$CONDA_BIN")
        fi

        # The binary is usually at .../bin/conda. We go up two levels to get the root.
        # $(dirname $(dirname ...)) essentially strips "/bin/conda"
        DETECTED_ROOT=$(dirname $(dirname "$REAL_PATH"))

        if [ -f "$DETECTED_ROOT/etc/profile.d/conda.sh" ]; then
            echo "✅ Detected Conda Root: $DETECTED_ROOT"
            CONDA_ROOT="$DETECTED_ROOT"
            return 0
        fi
    fi

    # Method 2: Check Standard Locations (Fallback if 'which' fails but folder exists)
    CANDIDATE_PATHS=(
        "$HOME/miniconda3"
        "$HOME/anaconda3"
        "$HOME/opt/miniconda3"
        "/opt/miniconda3"
        "/usr/local/miniconda3"
    )

    for path in "${CANDIDATE_PATHS[@]}"; do
        if [ -f "$path/etc/profile.d/conda.sh" ]; then
            echo "✅ Found Conda installation at: $path"
            CONDA_ROOT="$path"
            return 0
        fi
    done

    # If we reach here, we failed.
    return 1
}

# ---------------------------------------------------------
# MAIN LOGIC
# ---------------------------------------------------------

# 1. Run detection
if detect_conda; then
    # Source the conda.sh script to enable 'conda activate' command
    source "$CONDA_ROOT/etc/profile.d/conda.sh"
else
    echo "--------------------------------------------------------"
    echo "❌ ERROR: Conda (Miniconda or Anaconda) is NOT installed."
    echo "--------------------------------------------------------"
    echo "This tool requires Conda to manage dependencies."
    echo ""
    echo "👉 Please download and install Miniconda first:"
    echo "   🔗 Link: https://docs.anaconda.com/miniconda/"
    echo ""
    echo "Installation Quick-Start:"
    echo "1. Download the installer for your OS (Linux/Mac)."
    echo "2. Run the installer."
    echo "3. Restart your terminal."
    echo "4. Run this installer again: ./install.sh"
    echo "--------------------------------------------------------"
    exit 1
fi

# 2. Check for environment file
if [ ! -f "$ENV_FILE" ]; then
    echo "❌ Error: '$ENV_FILE' not found."
    echo "   Please run this script inside the 'search-ml' repository."
    exit 1
fi

echo "=========================================="
echo "🛠  Setting up Conda Environment..."
echo "=========================================="

# Extract environment name (default to search-ml if not found)
ENV_NAME=$(head -n 1 "$ENV_FILE" | awk '{print $2}')
[ -z "$ENV_NAME" ] && ENV_NAME="search-ml"

if conda env list | grep -q "$ENV_NAME"; then
    echo "🔄 Environment '$ENV_NAME' exists. Updating..."
    conda env update --file "$ENV_FILE" --prune
else
    echo "🆕 Creating environment '$ENV_NAME'..."
    conda env create --file "$ENV_FILE"
fi

conda env config vars set SEARCH_ML_HOME="$(pwd)" --name $ENV_NAME
conda activate $ENV_NAME

echo "=========================================="
echo "🎉 Setup Complete!"
echo "=========================================="
echo ""
echo "If $ENV_NAME not activated automatically, to activate the environment, run:"
echo "   conda activate"
echo "   conda env config vars set SEARCH_ML_HOME="$(pwd)" --name $ENV_NAME"
echo "   conda activate $ENV_NAME"
echo ""
