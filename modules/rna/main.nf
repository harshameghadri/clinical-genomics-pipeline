#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * RNA-seq processing pipeline
 * Includes QC, alignment, and quantification
 */

include { FASTQC as RNA_FASTQC } from '../qc/fastqc.nf'
include { STAR_ALIGN } from './alignment.nf'
include { SALMON_QUANT } from './quantification.nf'
include { CALCULATE_CHECKSUMS } from '../metadata/checksums.nf'
include { LOG_PROCESSING } from '../metadata/logging.nf'

workflow RNA_PIPELINE {
    take:
    rna_samples // [sample_id, fastq_files]
    
    main:
    // Quality control
    RNA_FASTQC(rna_samples)
    
    // Alignment with STAR
    STAR_ALIGN(rna_samples, params.reference_fasta, params.reference_gtf)
    
    // Quantification with Salmon
    SALMON_QUANT(rna_samples, params.reference_gtf)
    
    // Calculate checksums for reproducibility
    all_outputs = STAR_ALIGN.out.bam
        .mix(SALMON_QUANT.out.quant)
    
    CALCULATE_CHECKSUMS(all_outputs)
    
    // Log processing steps
    processing_logs = RNA_FASTQC.out.log
        .mix(STAR_ALIGN.out.log)
        .mix(SALMON_QUANT.out.log)
    
    LOG_PROCESSING(processing_logs.collect(), 'RNA')
    
    emit:
    bam = STAR_ALIGN.out.bam
    quant = SALMON_QUANT.out.quant
    qc = RNA_FASTQC.out.html
    metadata = LOG_PROCESSING.out.metadata
    checksums = CALCULATE_CHECKSUMS.out.checksums
}