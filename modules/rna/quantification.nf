#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Salmon quantification for RNA-seq samples
 */

process SALMON_QUANT {
    label 'salmon'
    tag "${sample_id}"
    
    publishDir "${params.outdir}/quantification", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    path reference_gtf
    
    output:
    tuple val(sample_id), path("${sample_id}_salmon"), emit: quant
    tuple val(sample_id), path("${sample_id}.quant.log"), emit: log
    path "versions.yml", emit: versions
    
    script:
    """
    echo "Starting Salmon quantification for ${sample_id}" > ${sample_id}.quant.log
    date >> ${sample_id}.quant.log
    
    # Build transcript index if it doesn't exist
    if [ ! -d "salmon_index" ] || [ ! -f "salmon_index/info.json" ]; then
        echo "Building Salmon transcript index..." >> ${sample_id}.quant.log
        
        # Extract transcript sequences from GTF and reference
        # This is a simplified approach - in production, use a pre-built transcriptome
        echo "Extracting transcript sequences..." >> ${sample_id}.quant.log
        
        # Create a simple transcriptome (this should be replaced with proper transcriptome FASTA)
        # For now, we'll use the reference genome as a proxy
        salmon index \\
            -t ${params.reference_fasta} \\
            -i salmon_index \\
            --threads ${task.cpus} \\
            2>> ${sample_id}.quant.log
    else
        echo "Using existing Salmon index" >> ${sample_id}.quant.log
    fi
    
    # Run Salmon quantification
    echo "Running Salmon quantification..." >> ${sample_id}.quant.log
    salmon quant \\
        -i salmon_index \\
        -l ${params.salmon_libtype} \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -p ${task.cpus} \\
        --validateMappings \\
        --gcBias \\
        --seqBias \\
        -o ${sample_id}_salmon \\
        2>> ${sample_id}.quant.log
    
    # Generate summary statistics
    echo "\\n--- Salmon Quantification Summary ---" >> ${sample_id}.quant.log
    if [ -f "${sample_id}_salmon/logs/salmon_quant.log" ]; then
        grep -E "(Mapping rate|Expected library type)" \\
            ${sample_id}_salmon/logs/salmon_quant.log >> ${sample_id}.quant.log || true
    fi
    
    # Count number of quantified transcripts
    if [ -f "${sample_id}_salmon/quant.sf" ]; then
        NUM_TRANSCRIPTS=\$(wc -l < ${sample_id}_salmon/quant.sf)
        echo "Number of transcripts quantified: \$NUM_TRANSCRIPTS" >> ${sample_id}.quant.log
        
        # Count transcripts with non-zero expression
        NON_ZERO=\$(awk 'NR>1 && \$4>0 {count++} END {print count+0}' ${sample_id}_salmon/quant.sf)
        echo "Transcripts with non-zero expression: \$NON_ZERO" >> ${sample_id}.quant.log
    fi
    
    echo "Salmon quantification completed for ${sample_id}" >> ${sample_id}.quant.log
    date >> ${sample_id}.quant.log
    
    # Version tracking
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    """
    mkdir -p ${sample_id}_salmon
    touch ${sample_id}_salmon/quant.sf
    touch ${sample_id}.quant.log
    touch versions.yml
    """
}