#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Calculate checksums for reproducibility
 */

process CALCULATE_CHECKSUMS {
    label 'checksums'
    tag "${sample_id}"
    
    publishDir "${params.outdir}/metadata/checksums", mode: 'copy'
    
    input:
    tuple val(sample_id), path(files)
    
    output:
    tuple val(sample_id), path("${sample_id}.checksums.txt"), emit: checksums
    path "versions.yml", emit: versions
    
    script:
    """
    # Calculate checksums for all files
    echo "# Checksums for sample: ${sample_id}" > ${sample_id}.checksums.txt
    echo "# Generated on: \$(date)" >> ${sample_id}.checksums.txt
    echo "# Algorithm: SHA256" >> ${sample_id}.checksums.txt
    echo "" >> ${sample_id}.checksums.txt
    
    for file in ${files}; do
        if [ -f "\$file" ]; then
            echo "Calculating checksum for: \$file"
            sha256sum "\$file" >> ${sample_id}.checksums.txt
        fi
    done
    
    # Create JSON format for programmatic access
    cat > ${sample_id}.checksums.json <<EOF
{
    "sample_id": "${sample_id}",
    "generated_date": "\$(date -Iseconds)",
    "algorithm": "SHA256",
    "checksums": {
EOF
    
    first=true
    for file in ${files}; do
        if [ -f "\$file" ]; then
            checksum=\$(sha256sum "\$file" | cut -d' ' -f1)
            filename=\$(basename "\$file")
            if [ "\$first" = true ]; then
                first=false
            else
                echo "," >> ${sample_id}.checksums.json
            fi
            echo "        \\"\\$filename\\": \\"\\$checksum\\"" >> ${sample_id}.checksums.json
        fi
    done
    
    cat >> ${sample_id}.checksums.json <<EOF
    }
}
EOF
    
    # Version tracking
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sha256sum: \$(sha256sum --version | head -n1 | cut -d' ' -f4)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.checksums.txt
    touch ${sample_id}.checksums.json
    touch versions.yml
    """
}