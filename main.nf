#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Open-source clinical genomics pipeline - Phase 1
 * Integrates DNA (WGS/WES) and RNA-seq processing
 * License: Apache-2.0
 */

include { DNA_PIPELINE } from './modules/dna/main.nf'
include { RNA_PIPELINE } from './modules/rna/main.nf'
include { QC_REPORT } from './modules/qc/multiqc.nf'
include { METADATA_DB } from './modules/metadata/duckdb.nf'
include { ANONYMIZE_SAMPLES } from './modules/metadata/gdpr.nf'

log.info """
Clinical Genomics Pipeline v1.0.0
==================================
DNA samples    : ${params.dna_samples}
RNA samples    : ${params.rna_samples}
Output dir     : ${params.outdir}
Reference      : ${params.reference}
"""

workflow {
    // Input validation
    if (!params.dna_samples && !params.rna_samples) {
        error "Must provide either DNA samples, RNA samples, or both"
    }

    // Initialize channels
    dna_samples_ch = params.dna_samples ? 
        Channel.fromPath(params.dna_samples)
               .map { it -> [it.getSimpleName(), it] } : 
        Channel.empty()
    
    rna_samples_ch = params.rna_samples ? 
        Channel.fromPath(params.rna_samples)
               .map { it -> [it.getSimpleName(), it] } : 
        Channel.empty()

    // GDPR compliance - anonymize sample IDs
    if (params.anonymize_samples) {
        ANONYMIZE_SAMPLES(dna_samples_ch.mix(rna_samples_ch))
        anonymized_dna = ANONYMIZE_SAMPLES.out.dna
        anonymized_rna = ANONYMIZE_SAMPLES.out.rna
    } else {
        anonymized_dna = dna_samples_ch
        anonymized_rna = rna_samples_ch
    }

    // Process DNA samples
    if (params.dna_samples) {
        DNA_PIPELINE(anonymized_dna)
        dna_results = DNA_PIPELINE.out
    }

    // Process RNA samples  
    if (params.rna_samples) {
        RNA_PIPELINE(anonymized_rna)
        rna_results = RNA_PIPELINE.out
    }

    // Generate comprehensive QC report
    all_qc_files = Channel.empty()
    if (params.dna_samples) all_qc_files = all_qc_files.mix(dna_results.qc)
    if (params.rna_samples) all_qc_files = all_qc_files.mix(rna_results.qc)
    
    QC_REPORT(all_qc_files.collect())

    // Store metadata in DuckDB
    all_metadata = Channel.empty()
    if (params.dna_samples) all_metadata = all_metadata.mix(dna_results.metadata)
    if (params.rna_samples) all_metadata = all_metadata.mix(rna_results.metadata)
    
    METADATA_DB(all_metadata.collect(), QC_REPORT.out.html)
}

workflow.onComplete {
    log.info """
    Pipeline completed!
    Results: ${params.outdir}
    Duration: ${workflow.duration}
    Success: ${workflow.success}
    """
}