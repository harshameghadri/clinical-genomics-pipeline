#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * CNVnator copy number variant detection
 */

process CNVNATOR {
    label 'cnvnator'
    tag "${sample_id}"
    
    publishDir "${params.outdir}/variants/cnv", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference_fasta
    
    output:
    tuple val(sample_id), path("${sample_id}.cnv.txt"), emit: cnv
    tuple val(sample_id), path("${sample_id}.cnv.log"), emit: log
    path "versions.yml", emit: versions
    
    script:
    def bin_size = 100  // 100bp bins for high resolution
    """
    echo "Starting CNVnator for ${sample_id}" > ${sample_id}.cnv.log
    date >> ${sample_id}.cnv.log
    
    # Extract read mapping from BAM
    echo "Extracting read mapping..." >> ${sample_id}.cnv.log
    cnvnator -root ${sample_id}.root -tree ${bam} \\
        2>> ${sample_id}.cnv.log
    
    # Generate histogram
    echo "Generating histogram (bin size: ${bin_size})..." >> ${sample_id}.cnv.log
    cnvnator -root ${sample_id}.root -his ${bin_size} -d ${reference_fasta} \\
        2>> ${sample_id}.cnv.log
    
    # Calculate statistics
    echo "Calculating statistics..." >> ${sample_id}.cnv.log
    cnvnator -root ${sample_id}.root -stat ${bin_size} \\
        2>> ${sample_id}.cnv.log
    
    # Partition and call CNVs
    echo "Partitioning and calling CNVs..." >> ${sample_id}.cnv.log
    cnvnator -root ${sample_id}.root -partition ${bin_size} \\
        2>> ${sample_id}.cnv.log
    
    cnvnator -root ${sample_id}.root -call ${bin_size} \\
        > ${sample_id}.cnv.txt \\
        2>> ${sample_id}.cnv.log
    
    # Filter high-quality CNVs (size > 1kb, p-value < 0.01)
    echo "Filtering CNVs..." >> ${sample_id}.cnv.log
    awk 'BEGIN{OFS="\\t"} 
         \$5-\$4 > 1000 && \$NF < 0.01 {print}' \\
         ${sample_id}.cnv.txt > ${sample_id}.cnv.filtered.txt
    
    echo "CNVnator completed for ${sample_id}" >> ${sample_id}.cnv.log
    date >> ${sample_id}.cnv.log
    
    # Version tracking
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvnator: \$(cnvnator 2>&1 | grep "CNVnator" | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.cnv.txt
    touch ${sample_id}.cnv.log
    touch versions.yml
    """
}