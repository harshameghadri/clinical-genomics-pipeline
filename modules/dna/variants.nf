#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * DeepVariant variant calling
 */

process DEEPVARIANT {
    label 'deepvariant'
    tag "${sample_id}"
    
    publishDir "${params.outdir}/variants/dna", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference_fasta
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), emit: vcf
    tuple val(sample_id), path("${sample_id}.gvcf.gz"), emit: gvcf
    tuple val(sample_id), path("${sample_id}.variants.log"), emit: log
    path "versions.yml", emit: versions
    
    script:
    """
    echo "Starting DeepVariant for ${sample_id}" > ${sample_id}.variants.log
    date >> ${sample_id}.variants.log
    
    # Run DeepVariant
    echo "Running DeepVariant variant calling..." >> ${sample_id}.variants.log
    run_deepvariant \\
        --model_type=${params.deepvariant_model} \\
        --ref=${reference_fasta} \\
        --reads=${bam} \\
        --output_vcf=${sample_id}.vcf.gz \\
        --output_gvcf=${sample_id}.gvcf.gz \\
        --num_shards=${task.cpus} \\
        --intermediate_results_dir=tmp_${sample_id} \\
        2>> ${sample_id}.variants.log
    
    # Index VCF
    echo "Indexing VCF files..." >> ${sample_id}.variants.log
    tabix -p vcf ${sample_id}.vcf.gz
    tabix -p vcf ${sample_id}.gvcf.gz
    
    # Filter variants by quality
    echo "Filtering variants (QUAL >= ${params.min_variant_quality})..." >> ${sample_id}.variants.log
    bcftools view -i "QUAL>=${params.min_variant_quality}" \\
        ${sample_id}.vcf.gz \\
        -O z -o ${sample_id}.filtered.vcf.gz
    tabix -p vcf ${sample_id}.filtered.vcf.gz
    
    # Generate variant statistics
    echo "Generating variant statistics..." >> ${sample_id}.variants.log
    bcftools stats ${sample_id}.vcf.gz > ${sample_id}.variant_stats.txt
    
    echo "DeepVariant completed for ${sample_id}" >> ${sample_id}.variants.log
    date >> ${sample_id}.variants.log
    
    # Version tracking
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(run_deepvariant --version 2>&1 | grep "DeepVariant version" | cut -d' ' -f3)
        bcftools: \$(bcftools --version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.vcf.gz
    touch ${sample_id}.vcf.gz.tbi
    touch ${sample_id}.gvcf.gz
    touch ${sample_id}.variants.log
    touch versions.yml
    """
}