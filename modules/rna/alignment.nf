#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * STAR alignment for RNA-seq samples
 */

process STAR_ALIGN {
    label 'star'
    tag "${sample_id}"
    
    publishDir "${params.outdir}/alignments/rna", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    path reference_fasta
    path reference_gtf
    
    output:
    tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam"), path("${sample_id}.Aligned.sortedByCoord.out.bam.bai"), emit: bam
    tuple val(sample_id), path("${sample_id}.alignment.log"), emit: log
    path "versions.yml", emit: versions
    
    script:
    """
    echo "Starting STAR alignment for ${sample_id}" > ${sample_id}.alignment.log
    date >> ${sample_id}.alignment.log
    
    # Check if STAR index exists, create if needed
    if [ ! -d "${params.reference_star_index}" ] || [ ! -f "${params.reference_star_index}/SA" ]; then
        echo "Creating STAR genome index..." >> ${sample_id}.alignment.log
        mkdir -p star_index_tmp
        
        STAR --runMode genomeGenerate \\
            --genomeDir star_index_tmp \\
            --genomeFastaFiles ${reference_fasta} \\
            --sjdbGTFfile ${reference_gtf} \\
            --sjdbOverhang 100 \\
            --runThreadN ${task.cpus} \\
            2>> ${sample_id}.alignment.log
        
        GENOME_DIR="star_index_tmp"
    else
        echo "Using existing STAR index: ${params.reference_star_index}" >> ${sample_id}.alignment.log
        GENOME_DIR="${params.reference_star_index}"
    fi
    
    # Run STAR alignment
    echo "Running STAR alignment..." >> ${sample_id}.alignment.log
    STAR --runMode alignReads \\
        --genomeDir \$GENOME_DIR \\
        --readFilesIn ${reads[0]} ${reads[1]} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${sample_id}. \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMstrandField intronMotif \\
        --outFilterIntronMotifs RemoveNoncanonical \\
        --outSAMattributes Standard \\
        --runThreadN ${task.cpus} \\
        --limitBAMsortRAM ${task.memory.toBytes()} \\
        2>> ${sample_id}.alignment.log
    
    # Index the BAM file
    echo "Indexing BAM file..." >> ${sample_id}.alignment.log
    samtools index ${sample_id}.Aligned.sortedByCoord.out.bam
    
    # Generate alignment statistics
    echo "Generating alignment statistics..." >> ${sample_id}.alignment.log
    samtools flagstat ${sample_id}.Aligned.sortedByCoord.out.bam >> ${sample_id}.alignment.log
    samtools stats ${sample_id}.Aligned.sortedByCoord.out.bam >> ${sample_id}.alignment.log
    
    # Parse STAR log for key metrics
    if [ -f "${sample_id}.Log.final.out" ]; then
        echo "\\n--- STAR Alignment Summary ---" >> ${sample_id}.alignment.log
        grep -E "(Number of input reads|Uniquely mapped reads|Number of reads mapped to multiple loci)" \\
            ${sample_id}.Log.final.out >> ${sample_id}.alignment.log
    fi
    
    echo "STAR alignment completed for ${sample_id}" >> ${sample_id}.alignment.log
    date >> ${sample_id}.alignment.log
    
    # Version tracking
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | head -n1 | cut -d'_' -f2)
        samtools: \$(samtools --version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.Aligned.sortedByCoord.out.bam
    touch ${sample_id}.Aligned.sortedByCoord.out.bam.bai
    touch ${sample_id}.alignment.log
    touch versions.yml
    """
}