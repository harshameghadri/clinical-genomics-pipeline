#!/bin/bash
set -euo pipefail

# Clinical Genomics Pipeline Runner
# Wrapper script for running the pipeline with proper setup

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Default values
PROFILE="local"
PARAMS_FILE=""
RESUME=""
WORK_DIR=""
OUTPUT_DIR=""

# Help function
show_help() {
    cat << EOF
Clinical Genomics Pipeline Runner

USAGE:
    $(basename "$0") [OPTIONS]

OPTIONS:
    -p, --profile PROFILE     Execution profile (local, slurm, aws, test) [default: local]
    -f, --params-file FILE    Parameters file (JSON format)
    -r, --resume              Resume previous run
    -w, --work-dir DIR        Work directory for temporary files
    -o, --output-dir DIR      Output directory [default: ./results]
    -h, --help                Show this help message

EXAMPLES:
    # Run with test data
    $(basename "$0") --profile test
    
    # Run on SLURM cluster with custom parameters
    $(basename "$0") --profile slurm --params-file my_params.json
    
    # Resume a previous run
    $(basename "$0") --resume --profile slurm

PROFILES:
    local   - Run locally with default settings
    test    - Run with test data (small, fast)
    slurm   - Run on SLURM cluster
    aws     - Run on AWS Batch
    azure   - Run on Azure Batch

SETUP:
    Before running, ensure:
    1. Nextflow is installed (>=22.10.0)
    2. Apptainer/Singularity is available
    3. Container images are built: ./containers/build_containers.sh
    4. Test data is prepared: ./data/test/download_test_data.sh

EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--profile)
            PROFILE="$2"
            shift 2
            ;;
        -f|--params-file)
            PARAMS_FILE="$2"
            shift 2
            ;;
        -r|--resume)
            RESUME="-resume"
            shift
            ;;
        -w|--work-dir)
            WORK_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Change to project directory
cd "$PROJECT_DIR"

echo "Clinical Genomics Pipeline Runner"
echo "================================="
echo "Profile: $PROFILE"
echo "Project directory: $PROJECT_DIR"

# Validate profile
if [[ ! "$PROFILE" =~ ^(local|test|slurm|aws|azure)$ ]]; then
    echo "Error: Invalid profile '$PROFILE'"
    echo "Valid profiles: local, test, slurm, aws, azure"
    exit 1
fi

# Check prerequisites
echo "Checking prerequisites..."

# Check Nextflow
if ! command -v nextflow &> /dev/null; then
    echo "Error: Nextflow not found. Please install Nextflow >= 22.10.0"
    exit 1
fi

NF_VERSION=$(nextflow -version 2>&1 | grep "nextflow version" | cut -d' ' -f3)
echo "✓ Nextflow version: $NF_VERSION"

# Check Apptainer/Singularity
if command -v apptainer &> /dev/null; then
    echo "✓ Apptainer found: $(apptainer --version)"
elif command -v singularity &> /dev/null; then
    echo "✓ Singularity found: $(singularity --version)"
else
    echo "Error: Neither Apptainer nor Singularity found"
    echo "Please install Apptainer for container support"
    exit 1
fi

# Check if containers exist (for non-test profiles)
if [[ "$PROFILE" != "test" ]]; then
    echo "Checking container availability..."
    if [ ! -f "containers/base.sif" ]; then
        echo "Warning: Container images not found"
        echo "Run './containers/build_containers.sh' to build containers"
        echo "Or the pipeline will attempt to build them automatically"
    else
        echo "✓ Container images found"
    fi
fi

# Setup work directory
if [[ -n "$WORK_DIR" ]]; then
    mkdir -p "$WORK_DIR"
    WORK_OPTION="-w $WORK_DIR"
else
    WORK_OPTION=""
fi

# Setup output directory
if [[ -n "$OUTPUT_DIR" ]]; then
    OUTPUT_OPTION="--outdir $OUTPUT_DIR"
else
    OUTPUT_OPTION=""
fi

# Setup parameters file
if [[ -n "$PARAMS_FILE" ]]; then
    if [[ ! -f "$PARAMS_FILE" ]]; then
        echo "Error: Parameters file not found: $PARAMS_FILE"
        exit 1
    fi
    PARAMS_OPTION="-params-file $PARAMS_FILE"
else
    PARAMS_OPTION=""
fi

# Build Nextflow command
NF_COMMAND="nextflow run main.nf"
NF_COMMAND="$NF_COMMAND -profile $PROFILE"
NF_COMMAND="$NF_COMMAND $RESUME"
NF_COMMAND="$NF_COMMAND $WORK_OPTION"
NF_COMMAND="$NF_COMMAND $OUTPUT_OPTION"
NF_COMMAND="$NF_COMMAND $PARAMS_OPTION"

echo ""
echo "Running pipeline..."
echo "Command: $NF_COMMAND"
echo ""

# Execute the pipeline
eval $NF_COMMAND

echo ""
echo "Pipeline execution completed!"

# Show results summary
if [[ -n "$OUTPUT_DIR" ]]; then
    RESULTS_DIR="$OUTPUT_DIR"
else
    RESULTS_DIR="./results"
fi

if [[ -d "$RESULTS_DIR" ]]; then
    echo ""
    echo "Results Summary:"
    echo "================"
    echo "Results directory: $RESULTS_DIR"
    
    # Count output files
    if [[ -d "$RESULTS_DIR/alignments" ]]; then
        BAM_COUNT=$(find "$RESULTS_DIR/alignments" -name "*.bam" 2>/dev/null | wc -l)
        echo "BAM files: $BAM_COUNT"
    fi
    
    if [[ -d "$RESULTS_DIR/variants" ]]; then
        VCF_COUNT=$(find "$RESULTS_DIR/variants" -name "*.vcf.gz" 2>/dev/null | wc -l)
        echo "VCF files: $VCF_COUNT"
    fi
    
    if [[ -d "$RESULTS_DIR/quantification" ]]; then
        QUANT_COUNT=$(find "$RESULTS_DIR/quantification" -name "quant.sf" 2>/dev/null | wc -l)
        echo "Quantification files: $QUANT_COUNT"
    fi
    
    if [[ -f "$RESULTS_DIR/qc/multiqc_report.html" ]]; then
        echo "✓ MultiQC report: $RESULTS_DIR/qc/multiqc_report.html"
    fi
    
    if [[ -f "$RESULTS_DIR/metadata/pipeline_metadata.duckdb" ]]; then
        echo "✓ Metadata database: $RESULTS_DIR/metadata/pipeline_metadata.duckdb"
    fi
fi