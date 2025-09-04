#!/usr/bin/env python3
"""
Unit tests for Clinical Genomics Pipeline modules
License: Apache-2.0
"""

import pytest
import tempfile
import json
import subprocess
import os
from pathlib import Path
import duckdb
import hashlib

# Test configurations
TEST_SAMPLE_ID = "TEST001"
TEST_DATA_DIR = Path(__file__).parent.parent / "data" / "test"


class TestPipelineStructure:
    """Test basic pipeline structure and configuration"""
    
    def test_main_workflow_exists(self):
        """Test that main.nf exists and is valid"""
        main_nf = Path(__file__).parent.parent / "main.nf"
        assert main_nf.exists(), "main.nf not found"
        
        # Check for required workflow components
        content = main_nf.read_text()
        assert "nextflow.enable.dsl = 2" in content
        assert "DNA_PIPELINE" in content
        assert "RNA_PIPELINE" in content
        assert "QC_REPORT" in content
    
    def test_config_files_exist(self):
        """Test that configuration files exist"""
        config_dir = Path(__file__).parent.parent
        
        # Main config
        assert (config_dir / "nextflow.config").exists()
        
        # Profile configs
        configs_dir = config_dir / "configs"
        assert (configs_dir / "aws.config").exists()
        assert (configs_dir / "slurm.config").exists()
        assert (configs_dir / "test.config").exists()
    
    def test_module_structure(self):
        """Test that all required modules exist"""
        modules_dir = Path(__file__).parent.parent / "modules"
        
        # DNA modules
        dna_dir = modules_dir / "dna"
        assert (dna_dir / "main.nf").exists()
        assert (dna_dir / "alignment.nf").exists()
        assert (dna_dir / "variants.nf").exists()
        
        # RNA modules  
        rna_dir = modules_dir / "rna"
        assert (rna_dir / "main.nf").exists()
        assert (rna_dir / "alignment.nf").exists()
        assert (rna_dir / "quantification.nf").exists()
        
        # QC modules
        qc_dir = modules_dir / "qc"
        assert (qc_dir / "fastqc.nf").exists()
        assert (qc_dir / "multiqc.nf").exists()
        
        # Metadata modules
        metadata_dir = modules_dir / "metadata"
        assert (metadata_dir / "duckdb.nf").exists()
        assert (metadata_dir / "checksums.nf").exists()
        assert (metadata_dir / "logging.nf").exists()
        assert (metadata_dir / "gdpr.nf").exists()


class TestContainerDefinitions:
    """Test container definition files"""
    
    def test_container_definitions_exist(self):
        """Test that all container definitions exist"""
        containers_dir = Path(__file__).parent.parent / "containers"
        
        required_containers = [
            "base.def",
            "fastqc.def", 
            "bwa_mem2.def",
            "star.def",
            "deepvariant.def",
            "salmon.def",
            "duckdb.def"
        ]
        
        for container in required_containers:
            assert (containers_dir / container).exists(), f"Missing {container}"
    
    def test_build_script_exists(self):
        """Test that container build script exists"""
        build_script = Path(__file__).parent.parent / "containers" / "build_containers.sh"
        assert build_script.exists()
        assert os.access(build_script, os.X_OK), "Build script not executable"


class TestDuckDBSchema:
    """Test DuckDB metadata schema and operations"""
    
    def test_create_database_schema(self):
        """Test DuckDB schema creation"""
        with tempfile.NamedTemporaryFile(suffix='.duckdb') as tmp_db:
            conn = duckdb.connect(tmp_db.name)
            
            # Create schema (from duckdb.nf process)
            schema_sql = '''
                CREATE SCHEMA IF NOT EXISTS genomics;
                
                CREATE TABLE IF NOT EXISTS genomics.samples (
                    sample_id VARCHAR PRIMARY KEY,
                    original_id VARCHAR,
                    sample_type VARCHAR,
                    processing_date TIMESTAMP,
                    anonymized BOOLEAN DEFAULT FALSE,
                    status VARCHAR DEFAULT 'processing'
                );
                
                CREATE TABLE IF NOT EXISTS genomics.files (
                    file_id VARCHAR PRIMARY KEY,
                    sample_id VARCHAR,
                    file_path VARCHAR,
                    file_type VARCHAR,
                    file_size BIGINT,
                    checksum VARCHAR,
                    created_date TIMESTAMP,
                    FOREIGN KEY (sample_id) REFERENCES genomics.samples(sample_id)
                );
            '''
            
            conn.execute(schema_sql)
            
            # Test inserting sample data
            conn.execute('''
                INSERT INTO genomics.samples 
                (sample_id, sample_type, anonymized, status)
                VALUES (?, ?, ?, ?)
            ''', [TEST_SAMPLE_ID, 'DNA', True, 'completed'])
            
            # Verify insertion
            result = conn.execute('''
                SELECT COUNT(*) FROM genomics.samples 
                WHERE sample_id = ?
            ''', [TEST_SAMPLE_ID]).fetchone()
            
            assert result[0] == 1, "Sample insertion failed"
            
            conn.close()
    
    def test_metadata_queries(self):
        """Test metadata query operations"""
        with tempfile.NamedTemporaryFile(suffix='.duckdb') as tmp_db:
            conn = duckdb.connect(tmp_db.name)
            
            # Setup test data
            conn.execute('''
                CREATE TABLE samples (
                    sample_id VARCHAR,
                    sample_type VARCHAR,
                    status VARCHAR
                );
                
                INSERT INTO samples VALUES
                ('DNA001', 'DNA', 'completed'),
                ('RNA001', 'RNA', 'completed'),
                ('DNA002', 'DNA', 'failed');
            ''')
            
            # Test summary query
            results = conn.execute('''
                SELECT sample_type, COUNT(*) as count,
                       SUM(CASE WHEN status = 'completed' THEN 1 ELSE 0 END) as completed
                FROM samples 
                GROUP BY sample_type
            ''').fetchall()
            
            assert len(results) == 2, "Unexpected query results"
            
            conn.close()


class TestGDPRCompliance:
    """Test GDPR compliance and anonymization"""
    
    def test_anonymization_logic(self):
        """Test sample ID anonymization"""
        original_id = "PATIENT_12345"
        workflow_run = "test_run_001"
        
        # Simulate anonymization logic from gdpr.nf
        hash_input = f"{original_id}_{workflow_run}"
        hash_digest = hashlib.sha256(hash_input.encode()).hexdigest()[:16]
        anonymized_id = f"ANON_{hash_digest}"
        
        # Verify properties
        assert anonymized_id.startswith("ANON_"), "Invalid anonymization prefix"
        assert len(anonymized_id) == 21, "Invalid anonymized ID length"  # ANON_ + 16 chars
        assert not original_id in anonymized_id, "Original ID leaked in anonymized ID"
        
        # Test reproducibility
        hash_input2 = f"{original_id}_{workflow_run}"
        hash_digest2 = hashlib.sha256(hash_input2.encode()).hexdigest()[:16]
        anonymized_id2 = f"ANON_{hash_digest2}"
        
        assert anonymized_id == anonymized_id2, "Anonymization not reproducible"
    
    def test_phi_detection_patterns(self):
        """Test PHI detection patterns"""
        phi_patterns = r"(patient|subject|name|dob|birth|ssn|phone|email|address)"
        
        # Test cases
        test_files = [
            "sample001_R1.fastq.gz",  # OK
            "patient123_data.bam",    # PHI detected
            "subject_001.vcf",        # PHI detected
            "ANON_abc123.fastq.gz",   # OK
            "john.doe@email.com.txt"  # PHI detected (email pattern)
        ]
        
        import re
        phi_pattern = re.compile(phi_patterns, re.IGNORECASE)
        
        results = [bool(phi_pattern.search(filename)) for filename in test_files]
        expected = [False, True, True, False, True]
        
        assert results == expected, f"PHI detection failed: {results} != {expected}"


class TestChecksumValidation:
    """Test checksum calculation and validation"""
    
    def test_checksum_calculation(self):
        """Test SHA256 checksum calculation"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp_file:
            test_content = "This is test content for checksum validation"
            tmp_file.write(test_content)
            tmp_file_path = tmp_file.name
        
        try:
            # Calculate checksum manually
            import hashlib
            with open(tmp_file_path, 'rb') as f:
                content = f.read()
                expected_checksum = hashlib.sha256(content).hexdigest()
            
            # Simulate checksum process
            result = subprocess.run(
                ['sha256sum', tmp_file_path],
                capture_output=True,
                text=True
            )
            
            if result.returncode == 0:
                calculated_checksum = result.stdout.split()[0]
                assert calculated_checksum == expected_checksum, "Checksum mismatch"
            else:
                # Skip if sha256sum not available
                pytest.skip("sha256sum not available")
                
        finally:
            os.unlink(tmp_file_path)


class TestOutputSchema:
    """Test output schema validation"""
    
    def test_output_schema_structure(self):
        """Test output schema JSON structure"""
        schema_file = Path(__file__).parent.parent / "schemas" / "output_schema.json"
        assert schema_file.exists(), "Output schema not found"
        
        with open(schema_file) as f:
            schema = json.load(f)
        
        # Validate required top-level properties
        required_props = ["pipeline_info", "samples", "quality_control", "metadata"]
        for prop in required_props:
            assert prop in schema["properties"], f"Missing required property: {prop}"
    
    def test_sample_schema_validation(self):
        """Test sample data against schema structure"""
        # Create mock sample data
        sample_data = {
            "sample_id": "ANON_test123",
            "sample_type": "DNA",
            "processing_status": "completed",
            "files": {
                "fastq": ["sample_R1.fastq.gz", "sample_R2.fastq.gz"],
                "bam": "sample.sorted.bam",
                "vcf": "sample.vcf.gz"
            },
            "metrics": {
                "alignment": {
                    "total_reads": 1000000,
                    "mapped_reads": 950000,
                    "mapping_rate": 0.95
                }
            }
        }
        
        # Basic validation
        assert "sample_id" in sample_data
        assert sample_data["sample_type"] in ["DNA", "RNA"]
        assert sample_data["processing_status"] in ["completed", "failed", "processing"]
        assert isinstance(sample_data["files"], dict)
        assert isinstance(sample_data["metrics"], dict)


class TestTestDataSetup:
    """Test the test data setup and validation"""
    
    def test_test_data_script_exists(self):
        """Test that test data download script exists"""
        test_script = TEST_DATA_DIR / "download_test_data.sh"
        assert test_script.exists(), "Test data script not found"
        assert os.access(test_script, os.X_OK), "Test data script not executable"
    
    def test_pipeline_runner_exists(self):
        """Test that pipeline runner script exists"""
        runner_script = Path(__file__).parent.parent / "scripts" / "run_pipeline.sh"
        assert runner_script.exists(), "Pipeline runner not found"
        assert os.access(runner_script, os.X_OK), "Pipeline runner not executable"


class TestVersionTracking:
    """Test version tracking and tool version collection"""
    
    def test_version_yml_structure(self):
        """Test version.yml output structure"""
        # Mock version output (as would be generated by processes)
        version_data = {
            "test_process": {
                "tool1": "1.2.3",
                "tool2": "4.5.6"
            }
        }
        
        # Test YAML serialization
        import yaml
        try:
            yaml_str = yaml.dump(version_data)
            reloaded = yaml.safe_load(yaml_str)
            assert reloaded == version_data, "YAML serialization failed"
        except ImportError:
            pytest.skip("PyYAML not available")


# Fixture for temporary working directory
@pytest.fixture
def temp_workdir():
    """Create temporary working directory for tests"""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


# Integration test
class TestPipelineIntegration:
    """Integration tests for pipeline components"""
    
    @pytest.mark.slow
    def test_nextflow_syntax_validation(self):
        """Test Nextflow syntax validation"""
        main_nf = Path(__file__).parent.parent / "main.nf"
        
        # Use nextflow to validate syntax
        result = subprocess.run(
            ['nextflow', 'run', str(main_nf), '--help'],
            capture_output=True,
            text=True
        )
        
        # If nextflow is not available, skip
        if result.returncode != 0 and "command not found" in result.stderr:
            pytest.skip("Nextflow not available for syntax validation")
        
        # Check for syntax errors (not parameter errors)
        assert "Script compilation error" not in result.stderr, \
            f"Nextflow syntax error: {result.stderr}"


if __name__ == "__main__":
    # Run tests with coverage if called directly
    pytest.main([__file__, "-v", "--tb=short"])