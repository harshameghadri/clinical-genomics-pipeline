#!/usr/bin/env python3
import random
import gzip
import os

def generate_fastq(filename, num_reads=1000, read_length=150):
    bases = ['A', 'T', 'G', 'C']
    qualities = 'I' * read_length  # High quality scores
    
    with gzip.open(filename, 'wt') as f:
        for i in range(num_reads):
            # Header
            f.write(f"@read_{i+1}\n")
            # Sequence
            seq = ''.join(random.choices(bases, k=read_length))
            f.write(f"{seq}\n")
            # Plus
            f.write("+\n")
            # Quality
            f.write(f"{qualities}\n")

# Generate test DNA data
os.makedirs('dna', exist_ok=True)
os.makedirs('rna', exist_ok=True)

generate_fastq('dna/NA12878_chr22_R1.fastq.gz', 5000, 150)
generate_fastq('dna/NA12878_chr22_R2.fastq.gz', 5000, 150)
print("Generated synthetic DNA FASTQ files")

# Generate test RNA data  
generate_fastq('rna/ENCODE_K562_R1.fastq.gz', 4000, 150)
generate_fastq('rna/ENCODE_K562_R2.fastq.gz', 4000, 150)
print("Generated synthetic RNA FASTQ files")