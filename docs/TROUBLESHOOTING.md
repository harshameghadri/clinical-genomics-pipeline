# Troubleshooting Guide

## Overview

This guide covers common issues, error messages, and solutions for the Clinical Genomics Pipeline. Issues are organized by category with step-by-step resolution procedures.

## General Troubleshooting Workflow

1. **Check logs**: Review Nextflow logs and process logs
2. **Verify resources**: Ensure adequate CPU, memory, and disk space
3. **Validate inputs**: Confirm input files exist and are properly formatted
4. **Test components**: Run individual processes or use test profile
5. **Check dependencies**: Verify all required software is installed
6. **Escalate**: Contact support with detailed error information

## Common Issues and Solutions

### Pipeline Startup Issues

#### Error: "Nextflow not found"
```bash
Command 'nextflow' not found
```
**Solution:**
```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Or add to PATH
export PATH=$HOME/bin:$PATH
echo 'export PATH=$HOME/bin:$PATH' >> ~/.bashrc
```

#### Error: "Java not found" 
```bash
Error: JAVA_HOME not set
```
**Solution:**
```bash
# Install Java 11+
sudo apt install openjdk-11-jdk

# Set JAVA_HOME
export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
echo 'export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64' >> ~/.bashrc
```

#### Error: "Container engine not found"
```bash
Cannot find Apptainer or Singularity
```
**Solution:**
```bash
# Install Apptainer
sudo apt update
sudo apt install -y apptainer

# Or install from source
wget https://github.com/apptainer/apptainer/releases/download/v1.2.0/apptainer_1.2.0_amd64.deb
sudo dpkg -i apptainer_1.2.0_amd64.deb
```

### Configuration Issues

#### Error: "Profile not found"
```bash
Unknown profile 'slurm'
```
**Solution:**
```bash
# Check available profiles
nextflow config -show-profiles

# Verify config file exists
ls -la configs/slurm.config

# Check syntax
nextflow config -profile slurm
```

#### Error: "Parameter file not found"
```bash
Cannot read params file: params.json
```
**Solution:**
```bash
# Check file exists and permissions
ls -la params.json
chmod 644 params.json

# Validate JSON syntax
python3 -m json.tool params.json

# Use absolute path if needed
nextflow run main.nf -params-file /full/path/to/params.json
```

### Resource Issues

#### Error: "Out of memory"
```bash
Process exceeded available memory
Exit status: 137 (SIGKILL - out of memory)
```
**Solution:**
```bash
# Increase memory allocation
--max_memory '128.GB'

# Or configure per-process
process {
    withLabel: 'bwa_mem2' {
        memory = { 64.GB * task.attempt }
        errorStrategy = 'retry'
        maxRetries = 3
    }
}

# Check system memory
free -h
cat /proc/meminfo
```

#### Error: "Disk space full"
```bash
No space left on device
```
**Solution:**
```bash
# Check disk usage
df -h
du -sh /path/to/work/dir

# Clean up work directory
nextflow clean -f

# Use different work directory
nextflow run main.nf -w /fast/scratch/work

# Enable automatic cleanup
nextflow run main.nf -w /tmp/work --cleanup
```

#### Error: "Too many open files"
```bash
Cannot create file: Too many open files
```
**Solution:**
```bash
# Increase file limits
ulimit -n 65536

# Permanent fix (add to /etc/security/limits.conf)
* soft nofile 65536
* hard nofile 65536

# For systemd services
echo "DefaultLimitNOFILE=65536" >> /etc/systemd/system.conf
systemctl daemon-reexec
```

### Container Issues

#### Error: "Container build failed"
```bash
FATAL: While performing build: conveyor failed to get: error fetching
```
**Solution:**
```bash
# Check internet connectivity
ping -c 3 google.com

# Clear Apptainer cache
apptainer cache clean -f

# Build with verbose output
apptainer build --fakeroot --verbose container.sif container.def

# Use alternative base image
# Edit .def file to use different Ubuntu mirror
```

#### Error: "Container not found"
```bash
Failed to pull container image
```
**Solution:**
```bash
# Check container file exists
ls -la containers/*.sif

# Build containers if missing
./containers/build_containers.sh

# Use absolute paths in config
container = 'file:///full/path/to/container.sif'

# Check permissions
chmod 644 containers/*.sif
```

#### Error: "Permission denied in container"
```bash
Permission denied: Cannot write to output directory
```
**Solution:**
```bash
# Fix ownership
sudo chown -R $USER:$USER /path/to/output

# Use --writable-tmpfs
apptainer.runOptions = '--writable-tmpfs'

# Bind mount with correct permissions  
apptainer.runOptions = '--bind /data:/data:rw'
```

### Cluster/HPC Issues

#### Error: "Job submission failed"
```bash
Unable to submit job: Invalid qos specification
```
**Solution:**
```bash
# Check SLURM configuration
sinfo -a
sacctmgr show qos
sacctmgr show associations user=$USER

# Update cluster options
process {
    clusterOptions = '--account=genomics --qos=normal'
}

# Test job submission manually
sbatch --account=genomics --qos=normal --wrap="echo test"
```

#### Error: "Queue not found"
```bash  
Invalid queue: genomics
```
**Solution:**
```bash
# List available queues
sinfo -s
qstat -q  # For PBS/Torque

# Update configuration
process.queue = 'normal'  // Use available queue name
```

#### Error: "Node requirements not met"
```bash
Job violates accounting/QOS policy
```
**Solution:**
```bash
# Check resource limits
sacctmgr show qos normal
scontrol show partition genomics

# Reduce resource requirements
process {
    cpus = 16    // Reduce from 32
    memory = '32.GB'  // Reduce from 64.GB
    time = '12.h'     // Reduce from 24.h
}
```

### Input/Output Issues

#### Error: "Sample sheet not found"
```bash
Cannot read sample sheet: samples.csv
```
**Solution:**
```bash
# Check file exists and format
ls -la samples.csv
head -5 samples.csv

# Verify CSV format (no extra spaces, proper headers)
sample_id,fastq_1,fastq_2,sample_type
SAMPLE001,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,DNA

# Use absolute paths
/full/path/to/samples.csv
```

#### Error: "FASTQ file not found"
```bash
No such file or directory: sample.fastq.gz
```
**Solution:**
```bash
# Check file exists and permissions
ls -la /path/to/sample.fastq.gz
file /path/to/sample.fastq.gz

# Verify FASTQ format
zcat sample.fastq.gz | head -8

# Check for symbolic links
readlink -f /path/to/sample.fastq.gz

# Test file integrity
zcat sample.fastq.gz | wc -l
```

#### Error: "Reference genome not found"
```bash
Cannot locate reference files
```
**Solution:**
```bash
# Check reference directory structure
ls -la /path/to/reference/
# Should contain: genome.fa, genes.gtf

# Build missing indexes
samtools faidx genome.fa
bwa-mem2 index genome.fa

# Use correct parameter
--reference /full/path/to/reference/dir
```

### Process-Specific Issues

#### BWA-MEM2 Issues
```bash
Error: [bwa_idx_load_from_disk] fail to load BWA index
```
**Solution:**
```bash
# Rebuild BWA-MEM2 index
bwa-mem2 index reference.fa

# Check index files exist
ls -la reference.fa.bwt.2bit.64

# Use correct reference path
params.reference_fasta = '/full/path/to/genome.fa'
```

#### STAR Alignment Issues
```bash
EXITING because of FATAL ERROR in reads input
```
**Solution:**
```bash
# Check FASTQ quality encoding
fastqc sample.fastq.gz

# Verify paired-end files match
zcat R1.fastq.gz | head -1
zcat R2.fastq.gz | head -1

# Check STAR index
ls -la star_index/SA star_index/Genome

# Reduce memory if needed
--limitGenomeGenerateRAM 31000000000
```

#### DeepVariant Issues
```bash
DeepVariant process failed with exit code 2
```
**Solution:**
```bash
# Check BAM file integrity
samtools quickcheck alignment.bam
samtools view -H alignment.bam

# Verify reference file
samtools faidx reference.fa
ls -la reference.fa.fai

# Check DeepVariant model compatibility
# WGS model for whole genome, WES for exome
--deepvariant_model WGS
```

#### Salmon Quantification Issues
```bash
Salmon index not found
```
**Solution:**
```bash
# Build salmon index manually
salmon index -t transcripts.fa -i salmon_index

# Or use reference genome (less optimal)
salmon index -t reference.fa -i salmon_index

# Check for transcriptome FASTA
ls -la transcripts.fa
```

### Network and Connectivity Issues

#### Error: "Connection timeout"
```bash
Connection timed out during download
```
**Solution:**
```bash
# Check internet connectivity
ping -c 3 8.8.8.8

# Use alternative mirrors
# Edit container definitions to use local mirrors

# Increase timeout values
curl --connect-timeout 60 --max-time 300 URL

# Use proxy if required
export http_proxy=http://proxy:8080
export https_proxy=http://proxy:8080
```

#### Error: "Certificate verification failed"
```bash
SSL certificate verification failed
```
**Solution:**
```bash
# Update certificates
sudo apt update && sudo apt install ca-certificates

# Temporary workaround (not recommended for production)
curl -k URL  # Skip certificate verification

# Configure proxy certificates if needed
```

### Cloud-Specific Issues

#### AWS Batch Issues
```bash
Job failed to start: RUNNABLE for too long
```
**Solution:**
```bash
# Check AWS Batch configuration
aws batch describe-compute-environments
aws batch describe-job-queues

# Verify IAM permissions
aws sts get-caller-identity

# Check instance limits
aws service-quotas get-service-quota --service-code ec2 --quota-code L-1216C47A

# Update queue configuration
aws batch update-compute-environment --compute-environment genomics-compute
```

#### S3 Access Issues
```bash
Access Denied: Unable to access S3 bucket
```
**Solution:**
```bash
# Check AWS credentials
aws sts get-caller-identity

# Verify bucket permissions
aws s3api get-bucket-policy --bucket your-bucket

# Test S3 access
aws s3 ls s3://your-bucket/

# Update IAM policy
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "s3:GetObject",
                "s3:PutObject",
                "s3:DeleteObject"
            ],
            "Resource": "arn:aws:s3:::your-bucket/*"
        }
    ]
}
```

### Database Issues

#### DuckDB Issues
```bash
Database file is corrupted
```
**Solution:**
```bash
# Check database integrity
duckdb metadata.duckdb "PRAGMA integrity_check;"

# Recover from backup
cp metadata_backup.duckdb metadata.duckdb

# Recreate database
rm metadata.duckdb
# Run pipeline again - database will be recreated
```

### Performance Issues

#### Slow Execution
**Symptoms:** Pipeline takes much longer than expected

**Diagnosis:**
```bash
# Check resource utilization
top
htop
iostat -x 1

# Review execution report
# Open reports/execution_report.html in browser

# Check network I/O for shared filesystems
iotop -ao
```

**Solutions:**
```bash
# Optimize resource allocation
process {
    withLabel: 'cpu_intensive' {
        cpus = { Math.min(task.attempt * 8, params.max_cpus) }
    }
}

# Use local storage for work directory
workDir = '/local/scratch/work'

# Enable process caching
process.cache = 'lenient'

# Parallelize sample processing
Channel.fromPath(samples).splitCsv().buffer(size: 5)
```

## Advanced Debugging

### Enable Verbose Logging
```bash
# Nextflow debug mode
nextflow -log .nextflow.log run main.nf -profile test --verbose

# Process-level debugging
process {
    debug = true
    echo = true
}
```

### Inspect Work Directory
```bash
# Find failed process work directory
ls -la work/*/*/

# Check command and logs
cd work/12/34567890abcdef
cat .command.sh    # Command executed
cat .command.log   # Process stdout
cat .command.err   # Process stderr  
cat .exitcode      # Exit code
```

### Container Debugging
```bash
# Run container interactively
apptainer shell container.sif

# Debug container build
apptainer build --sandbox debug_container container.def
apptainer shell --writable debug_container/

# Check container contents
apptainer exec container.sif ls /opt/tools
```

### Manual Process Execution
```bash
# Extract and run commands manually
cd work/failed/process/dir
bash .command.sh

# Run with different resources
export NXF_TASK_MEMORY=64GB
export NXF_TASK_CPUS=16
bash .command.sh
```

## Getting Help

### Information to Collect
When reporting issues, please provide:

1. **Pipeline version and commit hash**
2. **Complete error message and stack trace**
3. **Nextflow version** (`nextflow -version`)
4. **System information** (`uname -a`, memory, CPU cores)
5. **Configuration files** (sanitized for sensitive information)
6. **Input file examples** (small samples if possible)
7. **Work directory path** to failed process
8. **Execution profile used**

### Log Collection Script
```bash
#!/bin/bash
# collect_debug_info.sh

echo "=== System Information ===" > debug_info.txt
uname -a >> debug_info.txt
echo "" >> debug_info.txt

echo "=== Nextflow Version ===" >> debug_info.txt
nextflow -version >> debug_info.txt 2>&1
echo "" >> debug_info.txt

echo "=== Container Engine ===" >> debug_info.txt
which apptainer >> debug_info.txt 2>&1
apptainer --version >> debug_info.txt 2>&1
echo "" >> debug_info.txt

echo "=== Pipeline Configuration ===" >> debug_info.txt
nextflow config >> debug_info.txt 2>&1
echo "" >> debug_info.txt

echo "=== Recent Error ===" >> debug_info.txt
tail -100 .nextflow.log >> debug_info.txt 2>&1

echo "Debug information collected in debug_info.txt"
```

### Support Channels

- **GitHub Issues**: Technical issues and bug reports
- **GitHub Discussions**: Usage questions and community support  
- **Email**: Critical production issues
- **Documentation**: Check docs/ directory for detailed guides

### Emergency Procedures

For critical production issues:

1. **Stop pipeline**: `Ctrl+C` or `nextflow cancel <run_name>`
2. **Preserve logs**: Copy `.nextflow.log` and work directory
3. **Document issue**: Screenshot errors and system status
4. **Contact support**: Provide collected information
5. **Implement workaround**: Use alternative processes if available

This troubleshooting guide should resolve most common issues. For complex problems, systematic debugging with detailed logging usually reveals the root cause.