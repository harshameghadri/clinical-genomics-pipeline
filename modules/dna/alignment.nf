#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * BWA-MEM2 alignment for DNA samples
 */

process BWA_MEM2_ALIGN {
    label 'bwa_mem2'
    tag "${sample_id}"
    
    publishDir "${params.outdir}/alignments/dna", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    path reference_fasta
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam
    tuple val(sample_id), path("${sample_id}.alignment.log"), emit: log
    path "versions.yml", emit: versions
    
    script:
    def read_group = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA"
    """
    echo "Starting BWA-MEM2 alignment for ${sample_id}" > ${sample_id}.alignment.log
    date >> ${sample_id}.alignment.log
    
    # Index reference if needed
    if [ ! -f "${reference_fasta}.bwt.2bit.64" ]; then
        echo "Indexing reference genome..." >> ${sample_id}.alignment.log
        bwa-mem2 index ${reference_fasta}
    fi
    
    # Alignment
    echo "Running BWA-MEM2 alignment..." >> ${sample_id}.alignment.log
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "${read_group}" \\
        ${params.bwa_mem2_opts} \\
        ${reference_fasta} \\
        ${reads[0]} ${reads[1]} \\
        2>> ${sample_id}.alignment.log \\
        | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -
    
    # Index BAM
    echo "Indexing BAM file..." >> ${sample_id}.alignment.log
    samtools index ${sample_id}.sorted.bam
    
    # Generate stats
    echo "Generating alignment statistics..." >> ${sample_id}.alignment.log
    samtools flagstat ${sample_id}.sorted.bam >> ${sample_id}.alignment.log
    samtools stats ${sample_id}.sorted.bam >> ${sample_id}.alignment.log
    
    echo "BWA-MEM2 alignment completed for ${sample_id}" >> ${sample_id}.alignment.log
    date >> ${sample_id}.alignment.log
    
    # Version tracking
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1 | cut -d' ' -f2)
        samtools: \$(samtools --version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.sorted.bam
    touch ${sample_id}.sorted.bam.bai
    touch ${sample_id}.alignment.log
    touch versions.yml
    """
}