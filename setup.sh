#!/bin/bash
set -e

echo "==== [Nullomer Pipeline Setup] ===="

# 1. Checks if Conda is available
if ! command -v conda >/dev/null 2>&1; then
	echo "Conda wasn't found. Please, install Miniforge or Miniconda and try again. If already installed, check PATH"
	exit 1
fi

# 2. Creates conda environment if not already present
ENV_NAME="nullomer-env"
if ! conda info --envs | grep -q "$ENV_NAME"; then
	if [ -f environment.yaml ]; then
		echo "Creating conda environment '$ENV_NAME' sourcing environment.yaml..."
		conda env create -f environment.yaml
	else
		echo "No environment file found. Creating conda environment with basic packages."
		conda create -y -n $ENV_NAME python=3.11 numpy pandas biopython bitarray pyyaml psutil
	fi
else
	echo "Conda environment already exists. Nice :p"
fi

# 3. Activates conda environment
echo "Activating conda env..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $ENV_NAME

# 4. Install extra dependencies
if [ -f requirements.txt ]; then
	echo "Installing extra pip dependencies..."
	pip install -r requirements.txt
fi

# 5. Compile needed C files in the correct Scripts places
echo "Compiling C files..."
C_DIR="scripts/c"
for cfile in $C_DIR/*.c; do
	exe="${cfile%.c}"
	gcc -O2 -o "$exe" "$cfile"
	echo "$exe correctly compiled!"
done

# 6. Create data directory structure
echo "Creating data directory structure..."
DATA_DIR="data"
mkdir -p "$DATA_DIR"/{genomes,raw,processed} results logs benchmarks runs
echo "Data directory structure created at $DATA_DIR and under the root directory."
echo "==== Setup completed ===="
echo "For pipeline usage, remember to activate the conda environment with: conda activate $ENV_NAME"
