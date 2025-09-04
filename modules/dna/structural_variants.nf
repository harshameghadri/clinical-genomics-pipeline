#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Manta structural variant calling
 */

process MANTA {
    label 'manta'
    tag "${sample_id}"
    
    publishDir "${params.outdir}/variants/structural", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference_fasta
    
    output:
    tuple val(sample_id), path("${sample_id}.manta.vcf.gz"), path("${sample_id}.manta.vcf.gz.tbi"), emit: vcf
    tuple val(sample_id), path("${sample_id}.manta.log"), emit: log
    path "versions.yml", emit: versions
    
    script:
    """
    echo "Starting Manta SV calling for ${sample_id}" > ${sample_id}.manta.log
    date >> ${sample_id}.manta.log
    
    # Configure Manta
    echo "Configuring Manta..." >> ${sample_id}.manta.log
    configManta.py \\
        --bam ${bam} \\
        --referenceFasta ${reference_fasta} \\
        --runDir manta_${sample_id} \\
        2>> ${sample_id}.manta.log
    
    # Run Manta
    echo "Running Manta SV detection..." >> ${sample_id}.manta.log
    cd manta_${sample_id}
    ./runWorkflow.py -m local -j ${task.cpus} \\
        2>> ../${sample_id}.manta.log
    
    cd ..
    
    # Copy results
    cp manta_${sample_id}/results/variants/diploidSV.vcf.gz ${sample_id}.manta.vcf.gz
    cp manta_${sample_id}/results/variants/diploidSV.vcf.gz.tbi ${sample_id}.manta.vcf.gz.tbi
    
    # Generate SV statistics
    echo "Generating SV statistics..." >> ${sample_id}.manta.log
    bcftools stats ${sample_id}.manta.vcf.gz > ${sample_id}.sv_stats.txt
    
    echo "Manta SV calling completed for ${sample_id}" >> ${sample_id}.manta.log
    date >> ${sample_id}.manta.log
    
    # Version tracking
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$(configManta.py --version 2>&1 | grep "Manta workflow version" | cut -d' ' -f4)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.manta.vcf.gz
    touch ${sample_id}.manta.vcf.gz.tbi
    touch ${sample_id}.manta.log
    touch versions.yml
    """
}