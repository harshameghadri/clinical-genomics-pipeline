# Configuration Guide

## Overview

The Clinical Genomics Pipeline uses Nextflow's configuration system to manage parameters, resources, and execution environments. This guide covers all configuration options and best practices.

## Configuration Files

### Primary Configuration (`nextflow.config`)

The main configuration file defines default parameters and process settings:

```groovy
params {
    // Input/Output
    dna_samples = null
    rna_samples = null
    outdir = "./results"
    
    // Reference data
    reference = null
    reference_fasta = "${params.reference}/genome.fa"
    reference_gtf = "${params.reference}/genes.gtf"
    
    // Resources
    threads = 16
    memory = '32.GB'
    max_cpus = 32
    max_memory = '128.GB'
    max_time = '48.h'
}
```

### Profile-Specific Configurations

#### Local Profile (`nextflow.config`)
```groovy
profiles {
    local {
        process.executor = 'local'
    }
}
```

#### SLURM Profile (`configs/slurm.config`)
```groovy
process {
    executor = 'slurm'
    queue = 'genomics'
    clusterOptions = '--account=genomics --qos=normal'
    
    withLabel: 'bwa_mem2' {
        queue = 'compute'
        cpus = 32
        memory = '64.GB'
        time = '12.h'
    }
}
```

#### AWS Batch Profile (`configs/aws.config`)
```groovy
aws {
    region = 'us-east-1'
    batch {
        cliPath = '/usr/local/aws-cli/v2/current/bin/aws'
    }
}

process {
    executor = 'awsbatch'
    queue = 'genomics-queue'
}
```

## Parameter Categories

### Input Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `dna_samples` | String | `null` | Path to DNA samples CSV file |
| `rna_samples` | String | `null` | Path to RNA samples CSV file |
| `reference` | String | `null` | Reference genome directory |
| `outdir` | String | `"./results"` | Output directory |

### Reference Genome Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `reference_fasta` | String | `"${params.reference}/genome.fa"` | Reference FASTA file |
| `reference_gtf` | String | `"${params.reference}/genes.gtf"` | Gene annotation GTF file |
| `reference_star_index` | String | `"${params.reference}/star_index"` | STAR genome index directory |

### Resource Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `threads` | Integer | `16` | Default CPU threads per process |
| `memory` | String | `"32.GB"` | Default memory per process |
| `max_cpus` | Integer | `32` | Maximum CPUs available |
| `max_memory` | String | `"128.GB"` | Maximum memory available |
| `max_time` | String | `"48.h"` | Maximum runtime per process |

### Tool-Specific Parameters

#### BWA-MEM2 Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `bwa_mem2_opts` | String | `"-M -R"` | BWA-MEM2 command line options |

#### STAR Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `star_opts` | String | `"--outSAMtype BAM SortedByCoordinate"` | STAR command line options |

#### DeepVariant Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `deepvariant_model` | String | `"WGS"` | Model type: `WGS`, `WES`, or `PACBIO` |
| `min_variant_quality` | Integer | `20` | Minimum variant quality score |

#### Salmon Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `salmon_libtype` | String | `"A"` | Library type (A = auto-detect) |

### Quality Control Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `min_read_length` | Integer | `50` | Minimum read length for QC |
| `min_base_quality` | Integer | `20` | Minimum base quality for QC |

### GDPR Compliance Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `anonymize_samples` | Boolean | `true` | Enable sample anonymization |
| `remove_intermediate` | Boolean | `false` | Remove intermediate files |

### Container Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `container_engine` | String | `"apptainer"` | Container engine |
| `container_cache` | String | `"./containers"` | Container cache directory |

### Metadata Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `project_id` | String | `"clinical_genomics_v1"` | Project identifier |
| `analysis_date` | String | Auto-generated | Analysis date |

## Process-Specific Configuration

### Resource Labels

The pipeline uses process labels for resource allocation:

```groovy
process {
    withLabel: 'fastqc' {
        cpus = 4
        memory = '8.GB'
        container = 'file://containers/fastqc.sif'
    }
    
    withLabel: 'bwa_mem2' {
        cpus = { 16 * task.attempt }
        memory = { 32.GB * task.attempt }
        time = '12.h'
        container = 'file://containers/bwa_mem2.sif'
    }
}
```

### Dynamic Resource Scaling

Resources can scale with retry attempts:

```groovy
process {
    errorStrategy = 'retry'
    maxRetries = 2
    
    withLabel: 'memory_intensive' {
        memory = { 32.GB * task.attempt }
        time = { 6.h * task.attempt }
    }
}
```

## Environment-Specific Configurations

### HPC/SLURM Configuration

```groovy
// configs/slurm.config
process {
    executor = 'slurm'
    queue = 'genomics'
    
    // Global cluster options
    clusterOptions = '--account=genomics --qos=normal'
    
    // Process-specific queues
    withLabel: 'bwa_mem2' {
        queue = 'compute'
        clusterOptions = '--account=genomics --qos=long --partition=compute'
    }
    
    withLabel: 'deepvariant' {
        queue = 'gpu'
        clusterOptions = '--account=genomics --qos=gpu --partition=gpu --gres=gpu:1'
    }
}
```

### Cloud Configuration

#### AWS Batch
```groovy
// configs/aws.config
aws {
    region = 'us-east-1'
    batch {
        cliPath = '/usr/local/aws-cli/v2/current/bin/aws'
        volumes = ['/tmp']
    }
}

process {
    executor = 'awsbatch'
    queue = 'genomics-queue'
    
    withLabel: 'bwa_mem2' {
        queue = 'genomics-compute-optimized'
        cpus = 32
        memory = '64.GB'
    }
}

workDir = 's3://your-genomics-bucket/work'
```

#### Azure Batch
```groovy
// configs/azure.config  
azure {
    batch {
        location = 'eastus'
        accountName = 'your-batch-account'
    }
}

process {
    executor = 'azurebatch'
    queue = 'genomics-pool'
}
```

## Custom Parameter Files

### JSON Parameter File
```json
{
    "dna_samples": "/data/samples/dna_samples.csv",
    "rna_samples": "/data/samples/rna_samples.csv",
    "reference": "/data/reference/hg38",
    "outdir": "/data/results/project_001",
    "threads": 32,
    "memory": "64.GB",
    "deepvariant_model": "WGS",
    "min_variant_quality": 30,
    "anonymize_samples": true
}
```

### YAML Parameter File
```yaml
dna_samples: "/data/samples/dna_samples.csv"
rna_samples: "/data/samples/rna_samples.csv"  
reference: "/data/reference/hg38"
outdir: "/data/results/project_001"
threads: 32
memory: "64.GB"
deepvariant_model: "WGS"
min_variant_quality: 30
anonymize_samples: true
```

## Configuration Best Practices

### Resource Allocation

1. **Start Conservative**: Begin with default settings and scale up as needed
2. **Monitor Usage**: Use execution reports to optimize resource allocation
3. **Process-Specific Tuning**: Allocate resources based on process requirements
4. **Dynamic Scaling**: Use retry-based scaling for robust execution

### Sample Size Considerations

| Sample Count | Recommended Resources | Notes |
|--------------|----------------------|--------|
| 1-5 samples | 16 CPUs, 32 GB RAM | Local execution |
| 5-20 samples | 32 CPUs, 64 GB RAM | Small cluster |
| 20-50 samples | 64+ CPUs, 128+ GB RAM | HPC cluster |
| 50+ samples | Cloud deployment | Auto-scaling |

### Storage Configuration

```groovy
// Temporary work directory (fast storage)
workDir = '/fast/scratch/work'

// Results directory (permanent storage)  
params.outdir = '/data/results'

// Container cache (shared storage)
params.container_cache = '/shared/containers'
```

### Security Configuration

```groovy
// Disable Docker for security
docker.enabled = false

// Enable Apptainer security
apptainer.runOptions = '--containall --cleanenv'

// Secure work directory permissions
workDir.mode = '0750'
```

## Troubleshooting Configuration

### Common Issues

#### Resource Exhaustion
```groovy
// Increase memory for memory-intensive processes
process {
    withLabel: 'star' {
        memory = { 128.GB * task.attempt }
        errorStrategy = 'retry'
        maxRetries = 3
    }
}
```

#### Container Issues
```groovy  
// Use specific container paths
process {
    container = 'file:///absolute/path/to/container.sif'
}
```

#### Queue Configuration
```groovy
// Multiple queue fallback
process {
    queue = 'genomics,normal,default'
    
    withLabel: 'gpu' {
        queue = 'gpu,normal'
        errorStrategy = { task.attempt <= 2 ? 'retry' : 'ignore' }
    }
}
```

### Validation Commands

```bash
# Validate configuration syntax
nextflow config -profile slurm

# Show resolved configuration  
nextflow config -show-profiles

# Test configuration with dry run
nextflow run main.nf -profile test --outdir test_output -stub-run
```

## Advanced Configuration

### Conditional Configuration
```groovy
params {
    genome_version = 'hg38'
    
    reference_fasta = params.genome_version == 'hg38' ? 
        '/data/hg38/genome.fa' : 
        '/data/hg19/genome.fa'
}
```

### Environment Variables
```groovy
env {
    TMPDIR = '/fast/tmp'
    APPTAINER_CACHEDIR = '/shared/containers/cache'
}
```

### Custom Directives
```groovy
process {
    withLabel: 'custom' {
        // Custom cluster directive
        clusterOptions = { 
            task.memory.toGiga() > 64 ? 
                '--partition=highmem' : 
                '--partition=normal' 
        }
    }
}
```

This configuration system provides flexibility for different deployment scenarios while maintaining reproducibility and performance optimization.