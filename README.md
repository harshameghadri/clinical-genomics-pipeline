# Clinical Genomics Pipeline v1.0.0

An open-source, fully reproducible workflow integrating genomics (DNA) and transcriptomics (RNA) analysis for clinical applications.

## 📋 Overview

This Phase 1 foundation pipeline provides a complete, regulatory-compliant workflow for processing clinical genomic samples. Built with Nextflow and Apptainer, it offers:

- **DNA Analysis**: WGS/WES alignment, variant calling, structural variants, CNV detection  
- **RNA Analysis**: RNA-seq alignment, transcript quantification
- **Quality Control**: Comprehensive QC reports with FastQC and MultiQC
- **Metadata Management**: DuckDB-powered metadata tracking
- **GDPR Compliance**: Automated sample anonymization and PHI protection
- **Reproducibility**: Full version tracking, checksums, and containerization

## 🚀 Quick Start

### Prerequisites
- Nextflow ≥ 22.10.0
- Apptainer/Singularity
- Python 3.8+
- 16+ GB RAM, 32+ CPU cores recommended

### Installation

```bash
# Clone the repository
git clone <repository-url>
cd clinical-genomics-pipeline

# Build container images
./containers/build_containers.sh

# Setup test data
./data/test/download_test_data.sh

# Run tests
./tests/run_tests.sh

# Run pipeline with test data
./scripts/run_pipeline.sh --profile test
```

## 📁 Project Structure

```
clinical-genomics-pipeline/
├── main.nf                    # Main pipeline workflow
├── nextflow.config            # Primary configuration
├── configs/                   # Environment-specific configs
│   ├── aws.config            # AWS Batch configuration  
│   ├── slurm.config          # SLURM cluster configuration
│   └── test.config           # Test environment configuration
├── modules/                   # Pipeline modules
│   ├── dna/                  # DNA processing modules
│   │   ├── main.nf           # DNA pipeline workflow
│   │   ├── alignment.nf      # BWA-MEM2 alignment
│   │   ├── variants.nf       # DeepVariant calling
│   │   ├── structural_variants.nf # Manta SV detection
│   │   └── cnv.nf            # CNVnator CNV detection
│   ├── rna/                  # RNA processing modules  
│   │   ├── main.nf           # RNA pipeline workflow
│   │   ├── alignment.nf      # STAR alignment
│   │   └── quantification.nf # Salmon quantification
│   ├── qc/                   # Quality control modules
│   │   ├── fastqc.nf         # FastQC process
│   │   └── multiqc.nf        # MultiQC reporting
│   └── metadata/             # Metadata management
│       ├── duckdb.nf         # Database operations
│       ├── checksums.nf      # File integrity
│       ├── logging.nf        # Process logging  
│       └── gdpr.nf           # GDPR compliance
├── containers/               # Apptainer definitions
│   ├── base.def              # Base container
│   ├── fastqc.def            # FastQC container
│   ├── bwa_mem2.def          # BWA-MEM2 container
│   ├── star.def              # STAR container
│   ├── deepvariant.def       # DeepVariant container
│   ├── salmon.def            # Salmon container
│   ├── duckdb.def            # DuckDB container
│   └── build_containers.sh   # Build script
├── data/                     # Data directory
│   └── test/                 # Test datasets
│       ├── download_test_data.sh # Test data setup
│       ├── dna_samples.csv   # DNA sample sheet
│       └── rna_samples.csv   # RNA sample sheet
├── schemas/                  # Output schemas
│   └── output_schema.json    # Pipeline output schema
├── scripts/                  # Utility scripts
│   └── run_pipeline.sh       # Pipeline runner
├── tests/                    # Unit tests
│   ├── test_modules.py       # pytest test suite
│   ├── requirements.txt      # Test dependencies
│   └── run_tests.sh          # Test runner
└── docs/                     # Documentation
    ├── CONFIGURATION.md      # Configuration guide
    ├── DEPLOYMENT.md         # Deployment guide
    └── TROUBLESHOOTING.md    # Troubleshooting guide
```

## 🔧 Configuration

### Basic Parameters

```bash
# Input files
--dna_samples     CSV file with DNA sample information  
--rna_samples     CSV file with RNA sample information
--reference       Reference genome directory

# Resources  
--threads         Number of CPU threads (default: 16)
--memory          Memory allocation (default: 32.GB)
--outdir          Output directory (default: ./results)

# Analysis options
--deepvariant_model    DeepVariant model type (WGS/WES)
--min_variant_quality  Minimum variant quality score (default: 20)
--anonymize_samples    Enable GDPR compliance (default: true)
```

### Sample Sheet Format

**DNA samples (dna_samples.csv):**
```csv
sample_id,fastq_1,fastq_2,sample_type
SAMPLE001,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,DNA
```

**RNA samples (rna_samples.csv):**
```csv  
sample_id,fastq_1,fastq_2,sample_type
SAMPLE001,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,RNA
```

## 🖥️ Execution Profiles

### Local Execution
```bash
./scripts/run_pipeline.sh --profile local
```

### SLURM Cluster
```bash  
./scripts/run_pipeline.sh --profile slurm --params-file my_params.json
```

### AWS Batch
```bash
./scripts/run_pipeline.sh --profile aws --params-file aws_params.json
```

### Test Profile
```bash
./scripts/run_pipeline.sh --profile test
```

## 📊 Outputs

### Directory Structure
```
results/
├── alignments/           # BAM alignment files
│   ├── dna/             # DNA alignments (BWA-MEM2)
│   └── rna/             # RNA alignments (STAR)
├── variants/            # Variant call files
│   ├── dna/             # SNP/indel variants (DeepVariant)
│   ├── structural/      # Structural variants (Manta)  
│   └── cnv/             # Copy number variants (CNVnator)
├── quantification/      # Gene expression (Salmon)
├── qc/                  # Quality control reports
│   ├── fastqc/          # FastQC individual reports
│   └── multiqc_report.html # Consolidated QC report
├── metadata/            # Pipeline metadata
│   ├── pipeline_metadata.duckdb # Metadata database
│   ├── checksums/       # File integrity checksums
│   └── logs/            # Processing logs
└── reports/             # Execution reports
    ├── execution_report.html
    ├── timeline.html
    └── pipeline_dag.svg
```

### Key Output Files

**DNA Analysis:**
- `{sample}.sorted.bam` - Aligned reads (BWA-MEM2)
- `{sample}.vcf.gz` - Variant calls (DeepVariant)
- `{sample}.gvcf.gz` - Genomic variant calls
- `{sample}.manta.vcf.gz` - Structural variants (Manta)
- `{sample}.cnv.txt` - Copy number variants (CNVnator)

**RNA Analysis:**  
- `{sample}.Aligned.sortedByCoord.out.bam` - Aligned reads (STAR)
- `{sample}_salmon/quant.sf` - Transcript quantification (Salmon)

**Quality Control:**
- `multiqc_report.html` - Comprehensive QC report
- `{sample}_fastqc.html` - Individual FastQC reports

**Metadata:**
- `pipeline_metadata.duckdb` - Queryable metadata database
- `{sample}.checksums.txt` - File integrity checksums
- `anonymization_map.json` - GDPR anonymization mapping (secure storage)

## 🔒 GDPR Compliance

The pipeline includes comprehensive GDPR compliance features:

- **Automatic Anonymization**: Sample IDs are hashed using SHA256
- **PHI Detection**: Scans file names for potential personal information
- **Secure Mapping**: Original-to-anonymous ID mapping for authorized access
- **Data Retention**: Configurable retention policies
- **Audit Logging**: Complete processing audit trail

### Data Subject Rights

The pipeline supports GDPR data subject rights:
- Right to access: Contact your data controller
- Right to rectification: Contact your data controller  
- Right to erasure: Contact your data controller
- Right to portability: Data provided in standard formats (VCF, BAM, JSON)

## 🧪 Testing

### Run Test Suite
```bash
./tests/run_tests.sh
```

### Unit Tests
```bash  
cd tests
python -m pytest test_modules.py -v --cov=../modules
```

### Integration Test
```bash
./scripts/run_pipeline.sh --profile test
```

## 📈 Performance Optimization

### Resource Recommendations

**Small Scale (1-10 samples):**
- CPU: 16-32 cores
- Memory: 32-64 GB  
- Storage: 1-5 TB

**Medium Scale (10-50 samples):**
- CPU: 32-64 cores
- Memory: 64-128 GB
- Storage: 5-25 TB

**Large Scale (50-100+ samples):**
- CPU: 64+ cores  
- Memory: 128+ GB
- Storage: 25+ TB
- Consider cluster/cloud deployment

### Performance Tuning

1. **Adjust thread allocation** per process type
2. **Optimize memory** for alignment steps
3. **Use local SSD storage** for work directory
4. **Enable resume** for interrupted runs
5. **Tune container caching** for repeated runs

## 🐛 Troubleshooting

### Common Issues

**Container Build Failures:**
```bash  
# Clean and rebuild
rm -rf containers/*.sif
./containers/build_containers.sh
```

**Memory Issues:**
```bash
# Reduce resource allocation
--max_memory '64.GB' --max_cpus 16
```

**Permission Errors:**
```bash
# Check file permissions
chmod -R 755 data/
chmod +x scripts/*.sh containers/*.sh tests/*.sh
```

**SLURM Job Failures:**
```bash
# Check cluster configuration
squeue -u $USER
sinfo -p genomics
```

For detailed troubleshooting, see [TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md).

## 📚 Documentation

- [Configuration Guide](docs/CONFIGURATION.md) - Detailed configuration options
- [Deployment Guide](docs/DEPLOYMENT.md) - HPC and cloud deployment  
- [Troubleshooting Guide](docs/TROUBLESHOOTING.md) - Common issues and solutions

## 🤝 Contributing

1. Fork the repository
2. Create feature branch (`git checkout -b feature/amazing-feature`)
3. Run tests (`./tests/run_tests.sh`)
4. Commit changes (`git commit -m 'Add amazing feature'`)
5. Push to branch (`git push origin feature/amazing-feature`)
6. Open Pull Request

## 📄 License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **GIAB Consortium** - Reference datasets
- **ENCODE Project** - RNA-seq reference data  
- **nf-core Community** - Best practices and inspiration
- **Nextflow Team** - Workflow orchestration framework
- **Apptainer Team** - Container runtime

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/your-org/clinical-genomics-pipeline/issues)
- **Discussions**: [GitHub Discussions](https://github.com/your-org/clinical-genomics-pipeline/discussions)
- **Email**: genomics-support@your-institution.org

---

**Pipeline Version**: 1.0.0  
**Last Updated**: 2024-01-01  
**Compatibility**: Nextflow ≥ 22.10.0, Apptainer ≥ 1.0.0