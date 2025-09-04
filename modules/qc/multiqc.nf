#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * MultiQC aggregated quality control report
 */

process QC_REPORT {
    label 'multiqc'
    
    publishDir "${params.outdir}/qc", mode: 'copy'
    
    input:
    path qc_files
    
    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data", emit: data
    path "multiqc.log", emit: log
    path "versions.yml", emit: versions
    
    script:
    """
    echo "Starting MultiQC report generation" > multiqc.log
    date >> multiqc.log
    
    # Create MultiQC configuration
    cat > multiqc_config.yaml <<EOF
    title: "Clinical Genomics Pipeline QC Report"
    subtitle: "Phase 1 - DNA and RNA Analysis"
    intro_text: "Quality control metrics for genomics samples processed through the clinical pipeline."
    
    report_header_info:
        - Pipeline Version: '1.0.0'
        - Analysis Date: '${params.analysis_date}'
        - Project ID: '${params.project_id}'
    
    module_order:
        - fastqc
        - star
        - salmon
        - samtools
        - bcftools
    
    fastqc:
        status_checks: false
    
    star:
        status_checks: false
    
    salmon:
        status_checks: false
        
    remove_sections:
        - software_versions
    EOF
    
    # Run MultiQC
    echo "Running MultiQC..." >> multiqc.log
    multiqc \\
        --config multiqc_config.yaml \\
        --force \\
        --verbose \\
        --dirs \\
        --dirs-depth 3 \\
        --filename multiqc_report.html \\
        . \\
        2>> multiqc.log
    
    # Generate summary statistics
    echo "\\n--- MultiQC Summary ---" >> multiqc.log
    if [ -f "multiqc_data/multiqc_general_stats.txt" ]; then
        SAMPLE_COUNT=\$(wc -l < multiqc_data/multiqc_general_stats.txt)
        echo "Number of samples processed: \$((\$SAMPLE_COUNT - 1))" >> multiqc.log
    fi
    
    # Check for any QC failures
    if [ -f "multiqc_data/multiqc_fastqc.txt" ]; then
        FAILED_SAMPLES=\$(awk -F'\\t' 'NR>1 && \$NF=="FAIL" {print \$1}' multiqc_data/multiqc_fastqc.txt)
        if [ -n "\$FAILED_SAMPLES" ]; then
            echo "Samples with QC failures:" >> multiqc.log
            echo "\$FAILED_SAMPLES" >> multiqc.log
        fi
    fi
    
    echo "MultiQC report completed" >> multiqc.log
    date >> multiqc.log
    
    # Version tracking
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | head -n1 | cut -d' ' -f3)
    END_VERSIONS
    """
    
    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_data
    touch multiqc.log
    touch versions.yml
    """
}