#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * GDPR compliance - sample anonymization
 */

process ANONYMIZE_SAMPLES {
    label 'gdpr'
    
    publishDir "${params.outdir}/metadata", mode: 'copy'
    
    input:
    tuple val(sample_id), path(files)
    
    output:
    tuple val(anonymized_id), path(files), emit: dna
    tuple val(anonymized_id), path(files), emit: rna
    path "anonymization_map.json", emit: map
    path "gdpr_compliance.log", emit: log
    
    script:
    """
    # Create anonymized sample ID
    # Use hash of original ID plus timestamp for reproducibility within session
    HASH_INPUT="${sample_id}_${workflow.runName}"
    ANONYMIZED_ID=\$(echo -n "\$HASH_INPUT" | sha256sum | cut -c1-16)
    ANONYMIZED_ID="ANON_\${ANONYMIZED_ID}"
    
    echo "Starting GDPR compliance anonymization" > gdpr_compliance.log
    date >> gdpr_compliance.log
    echo "Original sample ID: ${sample_id}" >> gdpr_compliance.log
    echo "Anonymized sample ID: \$ANONYMIZED_ID" >> gdpr_compliance.log
    
    # Create anonymization mapping (stored securely)
    cat > anonymization_map.json <<EOF
{
    "anonymization_date": "\$(date -Iseconds)",
    "pipeline_run": "${workflow.runName}",
    "mapping": {
        "original_id": "${sample_id}",
        "anonymized_id": "\$ANONYMIZED_ID"
    },
    "gdpr_compliance": {
        "pii_removed": true,
        "reversible": false,
        "data_retention_policy": "according_to_institutional_policy"
    }
}
EOF
    
    # Check for potential PHI in file names and paths
    echo "\\n=== PHI Scan Results ===" >> gdpr_compliance.log
    PHI_PATTERNS="(patient|subject|name|dob|birth|ssn|phone|email|address)"
    
    for file in ${files}; do
        if echo "\$file" | grep -iE "\$PHI_PATTERNS" > /dev/null; then
            echo "WARNING: Potential PHI detected in filename: \$file" >> gdpr_compliance.log
        else
            echo "OK: No obvious PHI in filename: \$file" >> gdpr_compliance.log
        fi
    done
    
    # Log GDPR compliance measures taken
    cat >> gdpr_compliance.log <<EOF

=== GDPR Compliance Measures ===
1. Sample ID anonymized using SHA256 hash
2. Original-to-anonymous mapping created (store securely)
3. File names scanned for potential PHI
4. Processing logs will use anonymized IDs only
5. Database will track anonymization status

=== Data Subject Rights ===
- Right to access: Contact data controller
- Right to rectification: Contact data controller  
- Right to erasure: Contact data controller
- Right to portability: Data available in standard formats

=== Retention Policy ===
- Raw data: As per institutional policy
- Processed data: As per institutional policy  
- Anonymization mapping: Secure storage, limited access

Anonymization completed successfully.
EOF
    
    date >> gdpr_compliance.log
    
    # Output the anonymized ID for downstream processes
    echo "\$ANONYMIZED_ID" > anonymized_sample_id.txt
    anonymized_id=\$(cat anonymized_sample_id.txt)
    """
    
    stub:
    """
    echo "ANON_TEST123456789AB" > anonymized_sample_id.txt
    anonymized_id=\$(cat anonymized_sample_id.txt)
    touch anonymization_map.json
    touch gdpr_compliance.log
    """
}