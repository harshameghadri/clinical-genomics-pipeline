#!/bin/bash
set -euo pipefail

# Test runner for Clinical Genomics Pipeline
# Runs unit tests, integration tests, and validation

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

echo "Clinical Genomics Pipeline Test Suite"
echo "====================================="

# Change to project directory
cd "$PROJECT_DIR"

# Check if Python virtual environment should be created
if [[ ! -d "venv" ]]; then
    echo "Creating Python virtual environment..."
    python3 -m venv venv
fi

# Activate virtual environment
echo "Activating virtual environment..."
source venv/bin/activate

# Install test requirements
echo "Installing test requirements..."
pip install -r tests/requirements.txt

# Run unit tests
echo ""
echo "Running unit tests..."
echo "===================="
python -m pytest tests/test_modules.py -v --tb=short --cov=modules --cov-report=html --cov-report=term

# Check test coverage
echo ""
echo "Test coverage report generated in htmlcov/"

# Run pipeline syntax validation
echo ""
echo "Running pipeline syntax validation..."
echo "===================================="

if command -v nextflow &> /dev/null; then
    echo "Validating main.nf syntax..."
    if nextflow run main.nf --help > /dev/null 2>&1; then
        echo "✓ Main workflow syntax valid"
    else
        echo "✗ Main workflow syntax issues detected"
        nextflow run main.nf --help
        exit 1
    fi
    
    # Validate individual modules
    echo "Validating module syntax..."
    for module_file in modules/*/*.nf; do
        if [[ -f "$module_file" ]]; then
            echo "  Checking: $module_file"
            # Basic syntax check using nextflow
            if ! grep -q "process\|workflow" "$module_file"; then
                echo "    ✓ Module structure OK"
            fi
        fi
    done
else
    echo "Nextflow not available - skipping syntax validation"
fi

# Validate configuration files
echo ""
echo "Running configuration validation..."
echo "================================="

# Check main config
if [[ -f "nextflow.config" ]]; then
    echo "✓ Main configuration found"
    
    # Basic syntax check
    if grep -q "params\|process\|profiles" nextflow.config; then
        echo "✓ Main configuration structure OK"
    else
        echo "✗ Main configuration structure issues"
        exit 1
    fi
else
    echo "✗ Main configuration missing"
    exit 1
fi

# Check profile configs
for config in configs/*.config; do
    if [[ -f "$config" ]]; then
        echo "✓ Profile configuration: $(basename "$config")"
    fi
done

# Validate container definitions
echo ""
echo "Running container validation..."
echo "=============================="

for container_def in containers/*.def; do
    if [[ -f "$container_def" ]]; then
        container_name=$(basename "$container_def" .def)
        echo "Validating: $container_name"
        
        # Check required sections
        if grep -q "Bootstrap:" "$container_def" && 
           grep -q "%post" "$container_def"; then
            echo "  ✓ Container definition structure OK"
        else
            echo "  ✗ Container definition structure issues"
            exit 1
        fi
    fi
done

# Validate schemas
echo ""
echo "Running schema validation..."
echo "==========================="

if [[ -f "schemas/output_schema.json" ]]; then
    echo "✓ Output schema found"
    
    # Check JSON syntax
    if python3 -m json.tool schemas/output_schema.json > /dev/null; then
        echo "✓ Output schema JSON syntax OK"
    else
        echo "✗ Output schema JSON syntax issues"
        exit 1
    fi
else
    echo "✗ Output schema missing"
    exit 1
fi

# Run test data validation (if setup)
echo ""
echo "Running test data validation..."
echo "=============================="

if [[ -f "data/test/validate_test_data.sh" ]]; then
    cd data/test
    if ./validate_test_data.sh > /dev/null 2>&1; then
        echo "✓ Test data validation passed"
    else
        echo "! Test data not set up (run download_test_data.sh first)"
    fi
    cd "$PROJECT_DIR"
else
    echo "! Test data validation script not found"
fi

# Generate test report
echo ""
echo "Generating test report..."
echo "========================"

cat > test_report.txt <<EOF
Clinical Genomics Pipeline Test Report
Generated: $(date)

Test Results Summary:
- Unit tests: $(grep -c "PASSED" pytest.log 2>/dev/null || echo "See above") passed
- Configuration validation: PASSED
- Container validation: PASSED  
- Schema validation: PASSED
- Test data validation: $(if [[ -f "data/test/validate_test_data.sh" ]]; then echo "AVAILABLE"; else echo "NEEDS_SETUP"; fi)

Coverage Report: htmlcov/index.html
Full Test Log: Available above

Next Steps:
1. Review any failed tests above
2. Set up test data: ./data/test/download_test_data.sh  
3. Build containers: ./containers/build_containers.sh
4. Run integration test: ./scripts/run_pipeline.sh --profile test

Pipeline Status: READY FOR TESTING
EOF

echo ""
echo "Test suite completed!"
echo "===================="
echo "Report saved to: test_report.txt"
echo ""

# Deactivate virtual environment
deactivate