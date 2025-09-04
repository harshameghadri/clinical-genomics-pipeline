#!/bin/bash
set -euo pipefail

# Download test data for pipeline validation
# Uses GIAB NA12878 and ENCODE RNA-seq data

TEST_DATA_DIR="$(dirname "$0")"
cd "$TEST_DATA_DIR"

echo "Downloading test data for clinical genomics pipeline..."
echo "Target directory: $(pwd)"

# Create directories
mkdir -p {dna,rna,reference}

# Download GIAB NA12878 small dataset (chromosome 22 only)
echo "Downloading GIAB NA12878 test data (chr22)..."
if [ ! -f "dna/NA12878_chr22_R1.fastq.gz" ]; then
    # These are example URLs - replace with actual test data URLs
    echo "Note: Replace these URLs with actual test data locations"
    
    # For now, create placeholder files for testing pipeline structure
    echo "Creating placeholder DNA FASTQ files..."
    
    # Generate small synthetic FASTQ data for testing
    cat > generate_synthetic_fastq.py <<'EOF'
import random
import gzip

def generate_fastq(filename, num_reads=1000, read_length=150):
    bases = ['A', 'T', 'G', 'C']
    qualities = 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
    
    with gzip.open(filename, 'wt') as f:
        for i in range(num_reads):
            # Header
            f.write(f"@read_{i+1}\\n")
            # Sequence
            seq = ''.join(random.choices(bases, k=read_length))
            f.write(f"{seq}\\n")
            # Plus
            f.write("+\\n")
            # Quality
            qual = qualities[:read_length]
            f.write(f"{qual}\\n")

# Generate test DNA data
generate_fastq('dna/NA12878_chr22_R1.fastq.gz', 10000, 150)
generate_fastq('dna/NA12878_chr22_R2.fastq.gz', 10000, 150)
print("Generated synthetic DNA FASTQ files")

# Generate test RNA data  
generate_fastq('rna/ENCODE_K562_R1.fastq.gz', 8000, 150)
generate_fastq('rna/ENCODE_K562_R2.fastq.gz', 8000, 150)
print("Generated synthetic RNA FASTQ files")
EOF
    
    python3 generate_synthetic_fastq.py
    rm generate_synthetic_fastq.py
fi

# Download reference data (chromosome 22 subset)
echo "Setting up reference data..."
if [ ! -f "reference/chr22.fa" ]; then
    echo "Creating minimal reference chromosome 22..."
    
    # Create a minimal chromosome 22 sequence for testing
    cat > reference/chr22.fa <<'EOF'
>chr22
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF
    
    # Create GTF annotation
    cat > reference/chr22.gtf <<'EOF'
chr22	test	gene	1000	2000	.	+	.	gene_id "GENE001"; gene_name "TEST_GENE1"; gene_biotype "protein_coding";
chr22	test	transcript	1000	2000	.	+	.	gene_id "GENE001"; transcript_id "TRANS001"; gene_name "TEST_GENE1"; gene_biotype "protein_coding";
chr22	test	exon	1000	1200	.	+	.	gene_id "GENE001"; transcript_id "TRANS001"; exon_number "1";
chr22	test	exon	1800	2000	.	+	.	gene_id "GENE001"; transcript_id "TRANS001"; exon_number "2";
EOF
    
    ln -sf chr22.fa reference/genome.fa
    ln -sf chr22.gtf reference/genes.gtf
fi

# Create sample sheets
echo "Creating sample sheets..."
cat > dna_samples.csv <<'EOF'
sample_id,fastq_1,fastq_2,sample_type
NA12878_test,data/test/dna/NA12878_chr22_R1.fastq.gz,data/test/dna/NA12878_chr22_R2.fastq.gz,DNA
EOF

cat > rna_samples.csv <<'EOF'
sample_id,fastq_1,fastq_2,sample_type
K562_test,data/test/rna/ENCODE_K562_R1.fastq.gz,data/test/rna/ENCODE_K562_R2.fastq.gz,RNA
EOF

# Create test configuration
cat > test_params.json <<'EOF'
{
    "dna_samples": "data/test/dna_samples.csv",
    "rna_samples": "data/test/rna_samples.csv", 
    "reference": "data/test/reference",
    "outdir": "results/test_run",
    "anonymize_samples": true,
    "threads": 2,
    "memory": "4.GB"
}
EOF

# Create validation script
cat > validate_test_data.sh <<'EOF'
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
EOF

chmod +x validate_test_data.sh

echo "Test data setup completed!"
echo "Run './validate_test_data.sh' to verify the setup"
echo ""
echo "To test the pipeline:"
echo "  nextflow run ../../main.nf -profile test -params-file test_params.json"