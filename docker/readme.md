# INVAS Installation and Usage Guide

## Table of Contents

- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
  - [Option 1: Download Pre-built Image](#option-1-download-pre-built-image)
  - [Option 2: Build from Source](#option-2-build-from-source)
- [Quick Start](#quick-start)
- [Detailed Usage](#detailed-usage)
  - [Mounting Data Directories](#mounting-data-directories)
  - [Step 1: RNA-Seq Processing](#step-1-rna-seq-processing)
  - [Step 2: WGS Structural Variant Detection](#step-2-wgs-structural-variant-detection)
  - [Step 3: Candidate Inversion Detection](#step-3-candidate-inversion-detection)
  - [Step 4: Local Transcriptome Assembly](#step-4-local-transcriptome-assembly)
- [Complete Workflow Example](#complete-workflow-example)
- [HPC Job Submission](#hpc-job-submission)
  - [SLURM](#slurm)
  - [PBS](#pbs)
  - [SGE](#sge)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [System Requirements](#system-requirements)

---

## Introduction

INVAS (INtragenic inVersion AnalysiS) is a bioinformatics tool designed to detect intragenic inversions by integrating RNA-Seq and whole-genome sequencing (WGS) data. This guide provides step-by-step instructions for installing and running INVAS using Singularity containers.

Intragenic inversions are structural variants where a segment of DNA within a gene is reversed in orientation. These events can disrupt gene function and have been implicated in various diseases, including cancer.

---

## Prerequisites

Before installing INVAS, ensure you have:

- **Singularity** version 3.0 or higher installed
- At least 50 GB of free disk space
- Minimum 16 GB RAM (64 GB recommended)
- Linux operating system (tested on Ubuntu 20.04/22.04, CentOS 7/8)

To verify your Singularity installation:

```bash
singularity --version
```

## Installation
### Option 1: Download Pre-built Image
This is the fastest way to get started.

#### Create a directory for INVAS
mkdir -p ~/software/invas
cd ~/software/invas

### Option 2: Build from Source
If you need to customize the installation or the pre-built image is unavailable:
#### Download the pre-built Singularity image
wget https://github.com/your-organization/INVAS/releases/download/v1.0.0/invas_v1.0.0.sif

#### Rename for convenience
mv invas_v1.0.0.sif invas.sif

#### Verify the image
singularity run invas.sif --version

# Quick Start
Here is a minimal example to test your installation:

```
# Set up paths
export INVAS_SIF=~/software/invas/invas.sif
export DATA_DIR=/path/to/your/data
export OUTPUT_DIR=/path/to/output
export REF_GENOME=/path/to/reference/hg38.fa
export ANNOTATION=/path/to/annotation/gencode.v38.gtf

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Run Step 1 (RNA-Seq processing)
singularity run --bind ${DATA_DIR}:${DATA_DIR},${OUTPUT_DIR}:${OUTPUT_DIR} \
    ${INVAS_SIF} step1 \
    -n sample_001 \
    -b ${DATA_DIR}/rnaseq/sample_001.bam \
    -g ${REF_GENOME} \
    -a ${ANNOTATION} \
    -o ${OUTPUT_DIR}/step1
```


## Detailed Usage
Mounting Data Directories
Singularity containers have limited access to the host filesystem by default. You must explicitly bind directories containing your data:

```
# Bind a single directory
singularity run --bind /data:/data invas.sif [command]

# Bind multiple directories
singularity run --bind /data1:/data1,/data2:/data2,/scratch:/scratch invas.sif [command]

# Bind with different internal path
singularity run --bind /external/path:/internal/path invas.sif [command]
```


## Step 1: RNA-Seq Processing
This step extracts chimeric reads and splice junction information from RNA-Seq alignments.

Command:
```
singularity run --bind /data:/data ${INVAS_SIF} step1 \
    -n <sample_name> \
    -b <rnaseq_bam> \
    -g <reference_genome> \
    -a <annotation_gtf> \
    -o <output_directory> \
    [optional arguments]
```

## Step 2: WGS Structural Variant Detection
This step identifies structural variants from whole-genome sequencing data using multiple SV callers.

Command:
```
singularity run --bind /data:/data ${INVAS_SIF} step2 \
    -n <sample_name> \
    -b <wgs_bam> \
    -g <reference_genome> \
    -o <output_directory> \
    [optional arguments]
```

## Step 3: Candidate Inversion Detection
This step integrates RNA-Seq and WGS results to identify high-confidence intragenic inversion candidates.

Command:

```
singularity run --bind /data:/data ${INVAS_SIF} step3 \
    -r <step1_output> \
    -s <step2_output> \
    -g <reference_genome> \
    -a <annotation_gtf> \
    -o <output_directory> \
    [optional arguments]
```

## Step 4: Local Transcriptome Assembly
This step performs local transcriptome assembly around inversion breakpoints to validate candidates and reconstruct full transcripts.

Command:
```
singularity run --bind /data:/data ${INVAS_SIF} step4 \
    --input-dir <step3_output> \
    --rnaseq-bam <rnaseq_bam> \
    --genome <reference_genome> \
    --annotation <annotation_gtf> \
    --output-dir <output_directory> \
    [optional arguments]
```