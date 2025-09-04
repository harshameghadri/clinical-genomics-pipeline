#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * FastQC quality control process
 */

process FASTQC {
    label 'fastqc'
    tag "${sample_id}"
    
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("*_fastqc.html"), emit: html
    tuple val(sample_id), path("*_fastqc.zip"), emit: zip
    tuple val(sample_id), path("${sample_id}.fastqc.log"), emit: log
    path "versions.yml", emit: versions
    
    script:
    """
    echo "Starting FastQC for ${sample_id}" > ${sample_id}.fastqc.log
    date >> ${sample_id}.fastqc.log
    
    # Run FastQC
    echo "Running FastQC analysis..." >> ${sample_id}.fastqc.log
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        --format fastq \\
        --quiet \\
        ${reads} \\
        2>> ${sample_id}.fastqc.log
    
    # Extract key metrics from FastQC results
    echo "\\n--- FastQC Summary ---" >> ${sample_id}.fastqc.log
    for zip_file in *_fastqc.zip; do
        echo "Processing: \$zip_file" >> ${sample_id}.fastqc.log
        unzip -q \$zip_file
        fastqc_dir=\${zip_file%%.zip}
        
        # Extract basic statistics
        if [ -f "\$fastqc_dir/fastqc_data.txt" ]; then
            echo "Basic Statistics for \$fastqc_dir:" >> ${sample_id}.fastqc.log
            awk '/^>>Basic Statistics/,/^>>END_MODULE/' \$fastqc_dir/fastqc_data.txt \\
                | grep -v "^>>" >> ${sample_id}.fastqc.log
            echo "" >> ${sample_id}.fastqc.log
        fi
        
        # Check for failed modules
        if [ -f "\$fastqc_dir/summary.txt" ]; then
            FAILED=\$(awk '\$1=="FAIL" {print \$2}' \$fastqc_dir/summary.txt)
            if [ -n "\$FAILED" ]; then
                echo "FAILED QC modules for \$fastqc_dir:" >> ${sample_id}.fastqc.log
                echo "\$FAILED" >> ${sample_id}.fastqc.log
            fi
        fi
    done
    
    echo "FastQC completed for ${sample_id}" >> ${sample_id}.fastqc.log
    date >> ${sample_id}.fastqc.log
    
    # Version tracking
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}_R1_001_fastqc.html
    touch ${sample_id}_R1_001_fastqc.zip
    touch ${sample_id}_R2_001_fastqc.html
    touch ${sample_id}_R2_001_fastqc.zip
    touch ${sample_id}.fastqc.log
    touch versions.yml
    """
}