#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * DuckDB metadata database management
 */

process METADATA_DB {
    label 'duckdb'
    
    publishDir "${params.outdir}/metadata", mode: 'copy'
    
    input:
    path metadata_files
    path qc_report
    
    output:
    path "pipeline_metadata.duckdb", emit: database
    path "metadata_summary.json", emit: summary
    path "metadata.log", emit: log
    path "versions.yml", emit: versions
    
    script:
    """
    echo "Starting metadata database creation" > metadata.log
    date >> metadata.log
    
    # Create Python script for DuckDB operations
    cat > create_metadata_db.py <<'EOF'
import duckdb
import json
import pandas as pd
import os
import hashlib
from datetime import datetime
from pathlib import Path

def create_database():
    # Connect to DuckDB
    conn = duckdb.connect('pipeline_metadata.duckdb')
    
    # Create schema
    conn.execute('''
        CREATE SCHEMA IF NOT EXISTS genomics;
        
        CREATE TABLE IF NOT EXISTS genomics.samples (
            sample_id VARCHAR PRIMARY KEY,
            original_id VARCHAR,
            sample_type VARCHAR, -- 'DNA' or 'RNA'
            processing_date TIMESTAMP,
            anonymized BOOLEAN DEFAULT FALSE,
            status VARCHAR DEFAULT 'processing'
        );
        
        CREATE TABLE IF NOT EXISTS genomics.files (
            file_id VARCHAR PRIMARY KEY,
            sample_id VARCHAR,
            file_path VARCHAR,
            file_type VARCHAR, -- 'fastq', 'bam', 'vcf', 'quant', etc.
            file_size BIGINT,
            checksum VARCHAR,
            created_date TIMESTAMP,
            FOREIGN KEY (sample_id) REFERENCES genomics.samples(sample_id)
        );
        
        CREATE TABLE IF NOT EXISTS genomics.qc_metrics (
            metric_id VARCHAR PRIMARY KEY,
            sample_id VARCHAR,
            metric_type VARCHAR, -- 'fastqc', 'alignment', 'variant', etc.
            metric_name VARCHAR,
            metric_value DOUBLE,
            metric_unit VARCHAR,
            status VARCHAR, -- 'pass', 'warn', 'fail'
            FOREIGN KEY (sample_id) REFERENCES genomics.samples(sample_id)
        );
        
        CREATE TABLE IF NOT EXISTS genomics.processing_log (
            log_id VARCHAR PRIMARY KEY,
            sample_id VARCHAR,
            process_name VARCHAR,
            start_time TIMESTAMP,
            end_time TIMESTAMP,
            exit_code INTEGER,
            cpu_time DOUBLE,
            memory_usage BIGINT,
            tool_version VARCHAR,
            command_line TEXT,
            FOREIGN KEY (sample_id) REFERENCES genomics.samples(sample_id)
        );
        
        CREATE TABLE IF NOT EXISTS genomics.pipeline_runs (
            run_id VARCHAR PRIMARY KEY,
            pipeline_version VARCHAR,
            nextflow_version VARCHAR,
            start_time TIMESTAMP,
            end_time TIMESTAMP,
            success BOOLEAN,
            total_samples INTEGER,
            config_hash VARCHAR
        );
    ''')
    
    # Insert pipeline run information
    run_id = "${workflow.runName}"
    pipeline_version = "1.0.0"
    start_time = "${workflow.start}"
    
    conn.execute('''
        INSERT OR REPLACE INTO genomics.pipeline_runs 
        (run_id, pipeline_version, nextflow_version, start_time, total_samples)
        VALUES (?, ?, ?, ?, ?)
    ''', [run_id, pipeline_version, "${workflow.nextflow.version}", start_time, 0])
    
    # Process metadata files
    metadata_summary = {
        'pipeline_run': run_id,
        'creation_date': datetime.now().isoformat(),
        'total_samples': 0,
        'dna_samples': 0,
        'rna_samples': 0,
        'files_tracked': 0,
        'qc_metrics': 0
    }
    
    # Save summary
    with open('metadata_summary.json', 'w') as f:
        json.dump(metadata_summary, f, indent=2)
    
    conn.close()
    print("Metadata database created successfully")

if __name__ == "__main__":
    create_database()
EOF
    
    # Run the database creation script
    echo "Creating DuckDB metadata database..." >> metadata.log
    python3 create_metadata_db.py 2>> metadata.log
    
    # Create database queries script for analysis
    cat > query_metadata.sql <<'EOF'
-- Sample summary by type
.mode table
.headers on

SELECT 'Sample Summary' as report_section;
SELECT 
    sample_type,
    COUNT(*) as sample_count,
    COUNT(CASE WHEN status = 'completed' THEN 1 END) as completed,
    COUNT(CASE WHEN status = 'failed' THEN 1 END) as failed
FROM genomics.samples 
GROUP BY sample_type;

-- File summary by type
SELECT 'File Summary' as report_section;
SELECT 
    file_type,
    COUNT(*) as file_count,
    ROUND(SUM(file_size) / 1024.0 / 1024.0 / 1024.0, 2) as total_size_gb
FROM genomics.files 
GROUP BY file_type;

-- QC metrics summary
SELECT 'QC Summary' as report_section;
SELECT 
    metric_type,
    status,
    COUNT(*) as metric_count
FROM genomics.qc_metrics 
GROUP BY metric_type, status;
EOF
    
    echo "Metadata database creation completed" >> metadata.log
    date >> metadata.log
    
    # Version tracking
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        duckdb: \$(duckdb --version | head -n1 | cut -d' ' -f2)
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    """
    touch pipeline_metadata.duckdb
    echo '{"pipeline_run": "test", "creation_date": "2024-01-01", "total_samples": 0}' > metadata_summary.json
    touch metadata.log
    touch versions.yml
    """
}