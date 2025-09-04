#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * DNA processing pipeline
 * Includes QC, alignment, and variant calling
 */

include { FASTQC as DNA_FASTQC } from '../qc/fastqc.nf'
include { BWA_MEM2_ALIGN } from './alignment.nf'
include { DEEPVARIANT } from './variants.nf'
include { MANTA } from './structural_variants.nf'
include { CNVNATOR } from './cnv.nf'
include { CALCULATE_CHECKSUMS } from '../metadata/checksums.nf'
include { LOG_PROCESSING } from '../metadata/logging.nf'

workflow DNA_PIPELINE {
    take:
    dna_samples // [sample_id, fastq_files]
    
    main:
    // Quality control
    DNA_FASTQC(dna_samples)
    
    // Alignment
    BWA_MEM2_ALIGN(dna_samples, params.reference_fasta)
    
    // Variant calling
    DEEPVARIANT(BWA_MEM2_ALIGN.out.bam, params.reference_fasta)
    
    // Structural variants
    MANTA(BWA_MEM2_ALIGN.out.bam, params.reference_fasta)
    
    // Copy number variants
    CNVNATOR(BWA_MEM2_ALIGN.out.bam, params.reference_fasta)
    
    // Calculate checksums for reproducibility
    all_outputs = BWA_MEM2_ALIGN.out.bam
        .mix(DEEPVARIANT.out.vcf)
        .mix(MANTA.out.vcf)
        .mix(CNVNATOR.out.cnv)
    
    CALCULATE_CHECKSUMS(all_outputs)
    
    // Log processing steps
    processing_logs = DNA_FASTQC.out.log
        .mix(BWA_MEM2_ALIGN.out.log)
        .mix(DEEPVARIANT.out.log)
        .mix(MANTA.out.log)
        .mix(CNVNATOR.out.log)
    
    LOG_PROCESSING(processing_logs.collect(), 'DNA')
    
    emit:
    bam = BWA_MEM2_ALIGN.out.bam
    variants = DEEPVARIANT.out.vcf
    structural_variants = MANTA.out.vcf
    cnv = CNVNATOR.out.cnv
    qc = DNA_FASTQC.out.html
    metadata = LOG_PROCESSING.out.metadata
    checksums = CALCULATE_CHECKSUMS.out.checksums
}