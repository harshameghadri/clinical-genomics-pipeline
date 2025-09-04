# Deployment Guide

## Overview

This guide covers deployment of the Clinical Genomics Pipeline across different computing environments: local workstations, HPC clusters, and cloud platforms.

## Pre-deployment Checklist

### System Requirements

**Minimum Requirements:**
- CPU: 16 cores
- Memory: 32 GB RAM  
- Storage: 500 GB available space
- OS: Linux (Ubuntu 20.04+, CentOS 7+, RHEL 8+)

**Recommended Requirements:**
- CPU: 32+ cores
- Memory: 64+ GB RAM
- Storage: 2+ TB SSD
- Network: High-bandwidth connection for data transfer

### Software Dependencies

**Required Software:**
- Nextflow ≥ 22.10.0
- Java ≥ 11
- Apptainer ≥ 1.0.0 (or Singularity ≥ 3.8.0)
- Python ≥ 3.8

**Optional Software:**
- Git (for version control)
- AWS CLI (for cloud deployment)
- SLURM (for cluster deployment)

## Local Deployment

### Single Workstation Setup

1. **Install Nextflow:**
```bash
# Download and install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

2. **Install Apptainer:**
```bash  
# Ubuntu/Debian
sudo apt update
sudo apt install -y apptainer

# CentOS/RHEL
sudo yum install -y apptainer
```

3. **Clone and Setup Pipeline:**
```bash
git clone <repository-url>
cd clinical-genomics-pipeline

# Build containers
./containers/build_containers.sh

# Setup test data
./data/test/download_test_data.sh

# Run validation
./tests/run_tests.sh
```

4. **Test Execution:**
```bash
./scripts/run_pipeline.sh --profile test
```

### Multi-Node Local Cluster

For multiple connected workstations:

```groovy
// nextflow.config
profiles {
    local_cluster {
        process.executor = 'local'
        executor {
            queueSize = 100
            pollInterval = '5s'
        }
        
        // Distribute across available nodes
        process {
            clusterOptions = '--nodes=4 --ntasks-per-node=8'
        }
    }
}
```

## HPC Cluster Deployment

### SLURM Configuration

1. **Install on Head Node:**
```bash
# Install Nextflow (user installation)
curl -s https://get.nextflow.io | bash
mkdir -p ~/bin
mv nextflow ~/bin/

# Add to PATH
echo 'export PATH=$HOME/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

2. **Setup Shared Storage:**
```bash
# Create shared directories
sudo mkdir -p /shared/{containers,reference,work}
sudo chown -R genomics:genomics /shared/
chmod -R 755 /shared/
```

3. **Build Containers (on head node):**
```bash
cd /shared/pipelines/clinical-genomics-pipeline
./containers/build_containers.sh

# Containers will be available to all compute nodes
```

4. **SLURM Configuration File:**
```groovy
// configs/slurm.config
process {
    executor = 'slurm'
    queue = 'genomics'
    
    // Account and QOS settings
    clusterOptions = '--account=genomics --qos=normal'
    
    // Resource limits
    cpus = { Math.min(32, params.max_cpus) }
    memory = { Math.min(task.attempt * 32.GB, params.max_memory) }
    time = '24.h'
    
    // Error handling
    errorStrategy = 'retry'
    maxRetries = 2
    
    // Process-specific configuration
    withLabel: 'fastqc' {
        queue = 'short'
        cpus = 4
        memory = '8.GB'
        time = '1.h'
    }
    
    withLabel: 'bwa_mem2' {
        queue = 'compute'
        cpus = 32
        memory = '64.GB'
        time = '12.h'
        clusterOptions = '--account=genomics --qos=long --partition=compute'
    }
    
    withLabel: 'star' {
        queue = 'highmem'
        cpus = 16
        memory = '128.GB'
        time = '8.h'
        clusterOptions = '--account=genomics --qos=normal --partition=highmem'
    }
    
    withLabel: 'deepvariant' {
        queue = 'gpu'
        cpus = 8
        memory = '32.GB'
        time = '24.h'
        clusterOptions = '--account=genomics --qos=gpu --partition=gpu --gres=gpu:1'
    }
}

// Work directory on high-performance filesystem
workDir = '/fast/scratch/$USER/nextflow-work'

// Shared container cache
params.container_cache = '/shared/containers'
```

5. **Submit Jobs:**
```bash
# Interactive session for testing
srun --account=genomics --qos=normal --cpus-per-task=4 --mem=16G --time=2:00:00 --pty bash

# Submit pipeline job
sbatch <<EOF
#!/bin/bash
#SBATCH --account=genomics
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --job-name=genomics_pipeline

cd /shared/pipelines/clinical-genomics-pipeline
./scripts/run_pipeline.sh --profile slurm --params-file my_params.json
EOF
```

### PBS/Torque Configuration

```groovy
// configs/pbs.config
process {
    executor = 'pbs'
    queue = 'genomics'
    
    withLabel: 'bwa_mem2' {
        clusterOptions = '-l select=1:ncpus=32:mem=64gb:walltime=12:00:00'
    }
    
    withLabel: 'star' {
        clusterOptions = '-l select=1:ncpus=16:mem=128gb:walltime=8:00:00'
    }
}
```

### LSF Configuration

```groovy  
// configs/lsf.config
process {
    executor = 'lsf'
    queue = 'genomics'
    
    withLabel: 'bwa_mem2' {
        clusterOptions = '-n 32 -R "rusage[mem=2048]" -W 12:00'
    }
}
```

## Cloud Deployment

### AWS Batch Deployment

1. **Setup AWS Infrastructure:**
```bash
# Install AWS CLI
pip install awscli

# Configure credentials
aws configure

# Create S3 buckets
aws s3 mb s3://your-genomics-data
aws s3 mb s3://your-genomics-work
aws s3 mb s3://your-genomics-results
```

2. **AWS Batch Configuration:**
```groovy
// configs/aws.config
aws {
    region = 'us-east-1'
    batch {
        cliPath = '/usr/local/aws-cli/v2/current/bin/aws'
        volumes = ['/tmp', '/var/tmp']
        jobRole = 'arn:aws:iam::123456789:role/BatchExecutionRole'
    }
}

process {
    executor = 'awsbatch'
    queue = 'genomics-queue'
    
    // Container images from ECR
    container = '123456789.dkr.ecr.us-east-1.amazonaws.com/genomics-tools:latest'
    
    withLabel: 'bwa_mem2' {
        queue = 'genomics-compute-optimized'
        cpus = 32
        memory = '64.GB'
    }
    
    withLabel: 'star' {
        queue = 'genomics-memory-optimized' 
        cpus = 16
        memory = '128.GB'
    }
    
    withLabel: 'deepvariant' {
        queue = 'genomics-gpu'
        accelerator = 1
        cpus = 8
        memory = '32.GB'
    }
}

// S3 storage paths
workDir = 's3://your-genomics-work'
params.outdir = 's3://your-genomics-results'
```

3. **Launch Pipeline:**
```bash
nextflow run main.nf -profile aws -params-file aws_params.json
```

### Google Cloud Platform

1. **Setup GCP:**
```bash
# Install Google Cloud SDK
curl https://sdk.cloud.google.com | bash
exec -l $SHELL

# Initialize and authenticate
gcloud init
gcloud auth application-default login
```

2. **GCP Configuration:**
```groovy
// configs/gcp.config
google {
    project = 'your-project-id'
    zone = 'us-central1-a'
    
    lifeSciences {
        bootDiskSize = '50.GB'
        preemptible = true
    }
}

process {
    executor = 'google-lifesciences'
    
    withLabel: 'bwa_mem2' {
        machineType = 'n1-highmem-32'
        disk = '500.GB'
    }
}

workDir = 'gs://your-genomics-work'
```

### Azure Batch

1. **Setup Azure:**
```bash
# Install Azure CLI
curl -sL https://aka.ms/InstallAzureCLIDeb | sudo bash

# Login
az login
```

2. **Azure Configuration:**
```groovy
// configs/azure.config
azure {
    batch {
        location = 'eastus'
        accountName = 'your-batch-account'
        accountKey = 'your-batch-key'
        autoPoolMode = true
    }
}

process {
    executor = 'azurebatch'
    queue = 'genomics-pool'
    
    withLabel: 'bwa_mem2' {
        vmType = 'Standard_D32s_v3'
        vmCount = 1
    }
}
```

## Container Management

### Building for Different Platforms

1. **Local Build:**
```bash
# Build all containers
./containers/build_containers.sh

# Build specific container
apptainer build containers/bwa_mem2.sif containers/bwa_mem2.def
```

2. **Multi-Architecture Builds:**
```bash
# For ARM64 support
apptainer build --arch arm64 containers/bwa_mem2_arm64.sif containers/bwa_mem2.def
```

3. **Container Registry Integration:**
```bash
# Push to OCI registry
apptainer push container.sif oci://registry.example.com/genomics/bwa_mem2:v1.0

# Pull from registry  
apptainer pull oci://registry.example.com/genomics/bwa_mem2:v1.0
```

### Container Optimization

```dockerfile
# Multi-stage build for smaller images
Bootstrap: docker
From: ubuntu:22.04 as builder

%post
    # Build tools
    apt-get update && apt-get install -y build-essential
    # ... build process ...

Bootstrap: docker  
From: ubuntu:22.04 as runtime

%post
    # Only runtime dependencies
    apt-get update && apt-get install -y --no-install-recommends \
        samtools \
        bcftools
    # Copy built tools from builder stage
```

## Data Management

### Reference Genome Setup

1. **Download Reference Data:**
```bash
# Create reference directory
mkdir -p /data/reference/hg38

# Download reference genome
wget -O /data/reference/hg38/genome.fa.gz \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.fna.gz
gunzip /data/reference/hg38/genome.fa.gz

# Download GTF annotation
wget -O /data/reference/hg38/genes.gtf.gz \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip /data/reference/hg38/genes.gtf.gz

# Index reference
samtools faidx /data/reference/hg38/genome.fa
```

2. **Pre-build Indexes:**
```bash
# BWA-MEM2 index
bwa-mem2 index /data/reference/hg38/genome.fa

# STAR index
mkdir -p /data/reference/hg38/star_index
STAR --runMode genomeGenerate \
     --genomeDir /data/reference/hg38/star_index \
     --genomeFastaFiles /data/reference/hg38/genome.fa \
     --sjdbGTFfile /data/reference/hg38/genes.gtf \
     --sjdbOverhang 100 \
     --runThreadN 32
```

### Data Transfer

1. **Secure Transfer:**
```bash
# Using rsync with SSH
rsync -avz -e "ssh -i key.pem" \
    data/ user@cluster:/shared/data/

# Using SCP for small files
scp -r data/ user@cluster:/shared/data/
```

2. **Cloud Transfer:**
```bash
# AWS S3
aws s3 sync data/ s3://your-bucket/data/ --storage-class STANDARD_IA

# Google Cloud Storage  
gsutil -m cp -r data/ gs://your-bucket/data/

# Azure Blob Storage
az storage blob upload-batch -d data -s data/ --account-name youraccount
```

## Monitoring and Logging

### Pipeline Monitoring

1. **Nextflow Reports:**
```groovy
// Enable all reports
report {
    enabled = true
    file = "${params.outdir}/reports/execution_report.html"
}

timeline {
    enabled = true
    file = "${params.outdir}/reports/timeline.html"
}

trace {
    enabled = true
    file = "${params.outdir}/reports/trace.txt"
    fields = 'task_id,hash,native_id,process,tag,name,status,exit,module,container,cpus,time,disk,memory,attempt,submit,start,complete,duration,realtime,queue,%cpu,%mem,rss,vmem,peak_rss,peak_vmem,rchar,wchar,syscr,syscw,read_bytes,write_bytes'
}

dag {
    enabled = true
    file = "${params.outdir}/reports/pipeline_dag.svg"
}
```

2. **Resource Monitoring:**
```bash
# System monitoring
htop
iotop -o

# SLURM monitoring  
squeue -u $USER
sinfo -p genomics
sacct -j JOBID --format=JobID,JobName,MaxRSS,Elapsed
```

3. **Cloud Monitoring:**
```bash
# AWS CloudWatch
aws logs tail /aws/batch/job --follow

# GCP Logging
gcloud logging tail "resource.type=gce_instance"
```

### Log Management

```bash
# Centralized logging with rsyslog
echo "*.* @@log-server:514" >> /etc/rsyslog.conf

# Log rotation
cat > /etc/logrotate.d/nextflow <<EOF
/var/log/nextflow/*.log {
    daily
    rotate 30
    compress
    delaycompress
    missingok
    notifempty
    create 644 nextflow nextflow
}
EOF
```

## Security Considerations

### Access Control

1. **File Permissions:**
```bash
# Set appropriate permissions
chmod 750 /data/samples/
chmod 640 /data/samples/*.fastq.gz
chown -R genomics:genomics /data/
```

2. **Network Security:**
```bash
# Firewall rules (iptables)
iptables -A INPUT -p tcp --dport 22 -s trusted_network/24 -j ACCEPT
iptables -A INPUT -p tcp --dport 22 -j DROP
```

3. **Container Security:**
```groovy
// Secure container options
apptainer {
    runOptions = '--containall --cleanenv --no-home'
}
```

### Data Encryption

1. **Encryption at Rest:**
```bash
# LUKS encryption
cryptsetup luksFormat /dev/sdb
cryptsetup open /dev/sdb encrypted_data
mkfs.ext4 /dev/mapper/encrypted_data
```

2. **Encryption in Transit:**
```bash  
# SSH tunneling
ssh -L 8080:remote-host:80 user@jump-host

# TLS for web interfaces
openssl req -x509 -nodes -days 365 -newkey rsa:2048 \
    -keyout server.key -out server.crt
```

## Backup and Disaster Recovery

### Backup Strategy

1. **Data Backup:**
```bash
#!/bin/bash
# Daily backup script
BACKUP_DATE=$(date +%Y%m%d)
BACKUP_DIR="/backup/genomics/$BACKUP_DATE"

mkdir -p "$BACKUP_DIR"

# Backup critical data
rsync -av /data/samples/ "$BACKUP_DIR/samples/"
rsync -av /data/reference/ "$BACKUP_DIR/reference/"
rsync -av /shared/containers/ "$BACKUP_DIR/containers/"

# Backup databases
duckdb /data/metadata/pipeline_metadata.duckdb ".backup $BACKUP_DIR/metadata.db"

# Compress and transfer to remote storage
tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
aws s3 cp "$BACKUP_DIR.tar.gz" s3://backup-bucket/
```

2. **Configuration Backup:**
```bash
# Git-based configuration backup
git add configs/ containers/ modules/
git commit -m "Backup configuration $(date)"
git push origin main
```

### Disaster Recovery

1. **Recovery Procedures:**
```bash  
# Restore from backup
aws s3 cp s3://backup-bucket/backup.tar.gz ./
tar -xzf backup.tar.gz

# Rebuild containers if needed
./containers/build_containers.sh

# Verify system
./tests/run_tests.sh
```

2. **High Availability:**
```groovy
// Multiple work directories for redundancy
workDir = "/fast/scratch/work,/backup/scratch/work"

// Retry with different resources
process {
    errorStrategy = { task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries = 3
}
```

This deployment guide provides comprehensive coverage for different environments while maintaining security, monitoring, and disaster recovery best practices.