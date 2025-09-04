#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Process logging and metadata collection
 */

process LOG_PROCESSING {
    label 'logging'
    
    publishDir "${params.outdir}/metadata/logs", mode: 'copy'
    
    input:
    path log_files
    val pipeline_type // 'DNA' or 'RNA'
    
    output:
    path "${pipeline_type.toLowerCase()}_processing.log", emit: metadata
    path "processing_summary.json", emit: summary
    path "tool_versions.yml", emit: versions
    
    script:
    """
    # Consolidate all log files
    echo "Processing logs for ${pipeline_type} pipeline" > ${pipeline_type.toLowerCase()}_processing.log
    echo "Generated on: \$(date)" >> ${pipeline_type.toLowerCase()}_processing.log
    echo "Pipeline version: 1.0.0" >> ${pipeline_type.toLowerCase()}_processing.log
    echo "Nextflow version: ${workflow.nextflow.version}" >> ${pipeline_type.toLowerCase()}_processing.log
    echo "Run ID: ${workflow.runName}" >> ${pipeline_type.toLowerCase()}_processing.log
    echo "" >> ${pipeline_type.toLowerCase()}_processing.log
    
    # Aggregate individual log files
    sample_count=0
    for log_file in ${log_files}; do
        if [ -f "\$log_file" ]; then
            echo "=== Processing log: \$log_file ===" >> ${pipeline_type.toLowerCase()}_processing.log
            cat "\$log_file" >> ${pipeline_type.toLowerCase()}_processing.log
            echo "" >> ${pipeline_type.toLowerCase()}_processing.log
            sample_count=\$((sample_count + 1))
        fi
    done
    
    # Create processing summary
    cat > processing_summary.json <<EOF
{
    "pipeline_type": "${pipeline_type}",
    "run_id": "${workflow.runName}",
    "pipeline_version": "1.0.0",
    "nextflow_version": "${workflow.nextflow.version}",
    "start_time": "${workflow.start}",
    "samples_processed": \$sample_count,
    "log_files_collected": \$(ls ${log_files} 2>/dev/null | wc -l),
    "generated_date": "\$(date -Iseconds)"
}
EOF
    
    # Create comprehensive tool versions file
    cat > tool_versions.yml <<EOF
pipeline_info:
  name: "clinical-genomics-pipeline"
  version: "1.0.0"
  nextflow_version: "${workflow.nextflow.version}"
  run_id: "${workflow.runName}"

tool_versions:
EOF
    
    # Extract version information from individual process logs
    # This would be populated by individual processes
    echo "  # Tool versions will be populated by individual processes" >> tool_versions.yml
    
    echo "Processing log consolidation completed for ${pipeline_type} pipeline"
    """
    
    stub:
    """
    touch ${pipeline_type.toLowerCase()}_processing.log
    touch processing_summary.json
    touch tool_versions.yml
    """
}