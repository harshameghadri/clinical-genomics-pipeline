#!/bin/bash
set -euo pipefail

# Build script for all Apptainer containers
# License: Apache-2.0

CONTAINERS_DIR="$(dirname "$0")"
BUILD_DIR="${CONTAINERS_DIR}/build"

# Create build directory
mkdir -p "${BUILD_DIR}"

echo "Building Apptainer containers..."
echo "Build directory: ${BUILD_DIR}"

# List of containers to build
CONTAINERS=(
    "base"
    "fastqc"
    "bwa_mem2"
    "star"
    "deepvariant"
    "salmon"
    "duckdb"
)

# Build each container
for container in "${CONTAINERS[@]}"; do
    echo "Building ${container}..."
    apptainer build \
        "${BUILD_DIR}/${container}.sif" \
        "${CONTAINERS_DIR}/${container}.def"
    
    if [[ $? -eq 0 ]]; then
        echo "✓ Successfully built ${container}.sif"
        
        # Create symlink for easy access
        ln -sf "${BUILD_DIR}/${container}.sif" "${CONTAINERS_DIR}/${container}.sif"
    else
        echo "✗ Failed to build ${container}.sif"
        exit 1
    fi
done

echo "All containers built successfully!"
echo "Container files available in: ${BUILD_DIR}/"