#!/bin/bash
set -euo pipefail

echo "Validating test data setup..."

# Check DNA files
if [ -f "dna/NA12878_chr22_R1.fastq.gz" ] && [ -f "dna/NA12878_chr22_R2.fastq.gz" ]; then
    echo "✓ DNA FASTQ files present"
    echo "  - R1: $(wc -l < <(zcat dna/NA12878_chr22_R1.fastq.gz)) lines"
    echo "  - R2: $(wc -l < <(zcat dna/NA12878_chr22_R2.fastq.gz)) lines"
else
    echo "✗ Missing DNA FASTQ files"
    exit 1
fi

# Check RNA files
if [ -f "rna/ENCODE_K562_R1.fastq.gz" ] && [ -f "rna/ENCODE_K562_R2.fastq.gz" ]; then
    echo "✓ RNA FASTQ files present"
    echo "  - R1: $(wc -l < <(zcat rna/ENCODE_K562_R1.fastq.gz)) lines"
    echo "  - R2: $(wc -l < <(zcat rna/ENCODE_K562_R2.fastq.gz)) lines"
else
    echo "✗ Missing RNA FASTQ files"
    exit 1
fi

# Check reference files
if [ -f "reference/genome.fa" ] && [ -f "reference/genes.gtf" ]; then
    echo "✓ Reference files present"
    echo "  - Genome: $(wc -l < reference/genome.fa) lines"
    echo "  - GTF: $(wc -l < reference/genes.gtf) lines"
else
    echo "✗ Missing reference files"
    exit 1
fi

# Check sample sheets
if [ -f "dna_samples.csv" ] && [ -f "rna_samples.csv" ]; then
    echo "✓ Sample sheets present"
    echo "  - DNA samples: $(($(wc -l < dna_samples.csv) - 1))"
    echo "  - RNA samples: $(($(wc -l < rna_samples.csv) - 1))"
else
    echo "✗ Missing sample sheets"
    exit 1
fi

echo "Test data validation completed successfully!"
