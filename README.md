# INVAS

INVAS is a conjugate graph-based RNA assembler designed for the precise detection of **intragenomic inversions** and **transcriptome assembly**.  
It supports **two analysis modes**:

- **WGS + RNA mode (recommended)**: highest confidence inversion detection
- **RNA-only mode**: assembly and inversion discovery driven purely by RNA-seq data

INVAS leverages an **iterative conjugate maximum-flow network model** to reconstruct both **normal** and **inverted transcripts**.

---

## Overview

INVAS reconstructs expressed intragenomic inversions by combining structural variation signals and transcriptome evidence.  
Depending on data availability, users can choose between:

| Mode | Required Data | Use Case |
|-----|--------------|---------|
| **WGS + RNA** | RNA-seq + WGS BAM | Accurate inversion breakpoint detection |
| **RNA-only** | RNA-seq BAM | When WGS data is unavailable |

Both modes ultimately perform **conjugate graph-based transcriptome assembly**.

---

## Features

- Detection of intragenomic inversions
- **Dual running modes**: WGS+RNA and RNA-only
- Integration of multiple SV callers (WGS+RNA mode)
- Refinement of inverted splice junctions
- Conjugate graph-based transcript assembly
- Reconstruction of **normal and inverted isoforms**
- Re-alignment of unmapped RNA-seq reads
- Iterative conjugate maximum-flow network model
- **Built-in [IvasDB](https://github.com/yourusername/Invas/blob/main/IvasDB) database** for inversion annotation and cross-sample comparison


---

## Installation
- Docker guidance: [docker/readme.md](docker/readme.md)

### Prerequisites

- Conda / Mamba
- Reference genome (hg19 or hg38)
- HISAT2 index
- Gene annotation (GENCODE BED)

### Setup

```bash
git clone https://github.com/deepomicslab/INVAS.git
cd INVAS

2. Create conda environments:
```bash
# Data preparation environment
conda env create -f invas_prepare.yaml -n invas_data_prep

# Assembly environment
conda env create -f invas_assembly.yaml -n invas_assembly
```
2.1.  Build SSW library:

```
# Clone SSW library repository
git clone https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library.git

# Navigate to SSW src directory
cd Complete-Striped-Smith-Waterman-Library/src

# Compile SSW library
make

# Copy the compiled library to INVAS bin directory
cp libssw.so ../../scripts/full_pipe/bin/

# Return to INVAS main directory
cd ../../
```

2.2. Compile INVAS core components:
```# Navigate to the build directory
cd scripts/full_pipe/bin

# Run CMake to configure the build
cmake ../../../

# Compile the source code
make

# Make all executables in bin directory executable
chmod +x *

# Return to the main directory
cd ../../../ 
```

# Running Modes
INVAS supports two distinct workflows.
Please choose one based on your available data.

# Mode 1: WGS + RNA (Full Mode, Recommended)
This mode uses WGS-derived SV breakpoints to guide transcript assembly.


## Workflow

Before executing the INVAS core assembly module (`main.py`), three preprocessing steps are required to obtain informative reads for accurate inversion detection and transcript assembly.

### Step 1: RNA-Seq Data Processing (Preliminary)

Process RNA-seq data for quality control, alignment, and unmapped read extraction.

```bash
conda activate invas_data_prep

# Create output directory
mkdir -p rna_out/$sample_name

# Run RNA processing
bash scripts/full_pipe/preprocess/process_rna_common.sh -n $sample_name \
    -b $input_bam \
    -o rna_out/$sample_name \
    -r $reference_genome \
    -x $hisat_index \
    -t $threads \
    -p scripts/full_pipe/bin/picard.jar
```

Output files:
- `hisat.map_unmapremap.s.bam`: Mapped and remapped reads
- `still_unmap_bwa.s.bam`: Unmapped reads aligned with BWA

### Step 2: WGS Data Processing and SV Detection (Preliminary)

Process WGS data and detect structural variants using multiple SV callers.

```bash
conda activate invas_assembly

# Create output directory
mkdir -p sv_out/$sample_name

# Run SV detection
bash scripts/full_pipe/preprocess/process_sv_common.sh -s $sample_name \
    -i $input_wgs_bam \
    -o sv_out/$sample_name \
    -r $reference_genome \
    -t $threads \
    -c lumpy,delly,manta,svaba
```

### Step 3: Candidate Intragenic Inversion Event Detection (Preliminary)

Identify candidate inverted splicing events based on SV breakpoints.

```bash
# Create output directory
mkdir -p candidate_out

# Run candidate detection
bash scripts/full_pipe/preprocess/combine_sv_rna.sh   
    $sample_name \ 
    sv_out/$sample_name \
    rna_out/$sample_name/hisat.map_unmapremap.s.bam \
    rna_out/$sample_name/still_unmap_bwa.s.bam \
    candidate_out \
    $gene_annotation_bed \
    delly,manta,lumpy,svaba
```

Gene annotation files provided:
- hg19: `scripts/full_pipe/annotation/hg19_genes.gencode.bed`
- hg38: `scripts/full_pipe/annotation/hg38_genes.gencode.bed`

### Step 4: Core Transcriptome Assembly

Perform conjugate graph-based assembly to reconstruct transcripts.

```bash
conda activate invas_assembly

# Create output directory
mkdir -p assembly_out/$sample_name

# Run assembly
python scripts/full_pipe/main.py \
    --input_dir candidate_out/res/$sample_name \
    --sample_name $sample_name \
    --extend_region 5000 \
    --output_dir assembly_out/$sample_name \
    --ref_genome $reference_genome \
    --rna_map_with_remap_bam rna_out/$sample_name/hisat.map_unmapremap.s.bam \
    --rna_unmap_bwa_bam rna_out/$sample_name/still_unmap_bwa.s.bam \
    --wgs_bam $wgs_bam \
    --with_normal_haps
```

## Output Structure

The assembly results are organized as follows:
```
assembly_out/sample_name/
├── GENE1/                  # Gene-specific assembly results
│   ├── *_final_inv.bed     # Records the locations of inversions
│   ├── *_final_inv.exons   # Records the fragments involved in inversions
│   ├── haps/               # Folder containing assembled haplotype sequences
│   │   ├── hap1.fa         # Assembled sequence for haplotype 1
│   │   ├── hap2.fa         # Assembled sequence for haplotype 2
│   │   └── ...             # Additional haplotype sequences (if any)
│   └── summary.txt         # Assembly statistics for this gene
├── GENE2/
│   └── ...
├── ...
└── candidate_gene.txt       # Summary of all candidate genes with inversions
```

Each gene folder contains:
- Assembled transcript sequences in FASTA format
- Both normal and inverted isoforms
- Assembly quality metrics



# Mode 2: RNA-only Mode (No WGS Required)
This mode does not require WGS data.
Candidate inversions are inferred solely from RNA-seq signals (InvasDB, discordant alignments, unmapped read re-alignment).

### Step 1: RNA-seq Preprocessing
Same as WGS+RNA mode:
```
conda activate invas_data_prep
mkdir -p rna_out/$sample

bash scripts/full_pipe/preprocess/process_rna_common.sh \
    -n $sample \
    -b $input_rna_bam \
    -o rna_out/$sample \
    -r $ref \
    -x $hisat_index \
    -t $threads \
    -p scripts/full_pipe/bin/picard.jar
```

### Step 2: RNA-based Candidate Detection
```
mkdir -p candidate_out

bash scripts/full_pipe/preprocess/combine_rna_only.sh \
    $sample \
    rna_out/$sample/hisat.map_unmapremap.s.bam \
    rna_out/$sample/still_unmap_bwa.s.bam \
    candidate_out \
    $gene_bed
```
### Step 3: Assembly (RNA-only)


```bash
conda activate invas_assembly

mkdir -p assembly_out/$sample

python scripts/full_pipe/main_rna.py \
    --input_file candidate_out/rna_only/${sample}_inv_in_gene_dedup_highcof.txt \
    --sample_name $sample \
    --extend_region 5000 \
    --output_dir assembly_out/$sample \
    --ref_genome $ref \
    --ref_version hg38 \
    --rna_map_with_remap_bam rna_out/$sample/hisat.map_unmapremap.s.bam \
    --rna_unmap_bwa_bam rna_out/$sample/still_unmap_bwa.s.bam \
    --with_normal_haps
```

## Quick Start - Test Run 

Before running INVAS on your own data, we recommend testing the installation with the provided example scripts:

```
# Navigate to the test directory
cd INVAS/test

# Review the example scripts for different scenarios:
ls -la *.sh

# These test scripts demonstrate:
# - Complete pipeline execution
# - Parameter settings for different data types
# - Output structure and interpretation
```




## Example Usage (WGS+RNA)

A complete example workflow for processing sample "3-06PB4872":

```bash
#!/bin/bash

# Set up variables
sample_name="3-06PB4872"
threads=16

# Reference files (hg19 example)
ref="/path/to/references/hg19.fa"
hisat_ref="/path/to/references/hisat2_index/hg19"
gene_bed="scripts/full_pipe/annotation/hg19_genes.gencode.bed"

# Input files
rna_fq1="/path/to/data/${sample_name}_R1.fastq.gz"
rna_fq2="/path/to/data/${sample_name}_R2.fastq.gz"
wgs_bam="/path/to/data/${sample_name}_wgs.bam"

# Output directories
base_dir="/path/to/output"
rna_dir="$base_dir/rna_out"
sv_dir="$base_dir/sv_out"
candidate_dir="$base_dir/candidate_out"
assembly_dir="$base_dir/assembly_out"

# Step 1: RNA-Seq Data Processing (Preliminary)
echo "Step 1: Processing RNA-seq data..."
conda activate invas_assembly 

mkdir -p $rna_dir/$sample_name

# Initial alignment with HISAT2
hisat2 -p $threads -x $hisat_ref \
    -1 $rna_fq1 -2 $rna_fq2 | \
    samtools sort -@ $threads -o $rna_dir/$sample_name/hisat2.bam

# Process with Invas RNA pipeline
bash scripts/test_common_run.sh \
    -n $sample_name \
    -b $rna_dir/$sample_name/hisat2.bam \
    -o $rna_dir/$sample_name \
    -r $ref \
    -x $hisat_ref \
    -t $threads \
    -p tools/picard.jar

# Step 2: WGS Data Processing and SV Detection (Preliminary)
echo "Step 2: Detecting structural variants..."
conda activate invas_data_prep

mkdir -p $sv_dir/$sample_name

# Extract MHC region (optional, for focused analysis)
samtools view $wgs_bam chr6:28477797-33448354 -O BAM \
    -o $sv_dir/$sample_name/${sample_name}_mhc.bam
samtools index -@ $threads $sv_dir/$sample_name/${sample_name}_mhc.bam

# Run SV detection
bash scripts/process_sv_common.sh \
    -s $sample_name \
    -i $sv_dir/$sample_name/${sample_name}_mhc.bam \
    -o $sv_dir/$sample_name \
    -r $ref \
    -t $threads \
    -c lumpy,delly,manta,svaba

# Step 3: Candidate Intragenic Inversion Event Detection (Preliminary)
echo "Step 3: Detecting candidate inversions..."
mkdir -p $candidate_dir

bash scripts/full_pipe/preprocess/combine_sv_rna.sh $sample_name \
    $sv_dir/$sample_name \
    $rna_dir/$sample_name/hisat.map_unmapremap.s.bam \
    $rna_dir/$sample_name/still_unmap_bwa.s.bam \
    $candidate_dir \
    $gene_bed \
    delly,manta,lumpy,svaba

# Step 4: Core Transcriptome Assembly
echo "Step 4: Assembling transcripts..."
conda activate invas_assembly

mkdir -p $assembly_dir/$sample_name

python scripts/full_pipe/main.py \
    --input_dir $candidate_dir/res/$sample_name \
    --sample_name $sample_name \
    --extend_region 5000 \
    --output_dir $assembly_dir/$sample_name \
    --ref_genome $ref \
    --rna_map_with_remap_bam $rna_dir/$sample_name/hisat.map_unmapremap.s.bam \
    --rna_unmap_bwa_bam $rna_dir/$sample_name/still_unmap_bwa.s.bam \
    --wgs_bam $wgs_bam \
    --with_normal_haps

echo "Analysis complete! Results in: $assembly_dir/$sample_name"
```

## Batch Processing

For processing multiple samples:

```bash
#!/bin/bash

# Sample list file (one sample per line)
samples_file="samples.txt"
threads=16

# Reference files
ref="/path/to/references/hg19.fa"
hisat_ref="/path/to/references/hisat2_index/hg19"
gene_bed="scripts/full_pipe/annotation/hg19_genes.gencode.bed"

# Base directories
base_dir="/path/to/output"

# Process each sample
while IFS= read -r sample; do
    echo "Processing sample: $sample"
    
    # Run complete pipeline for each sample
    bash run_invas_pipeline.sh \
        --sample $sample \
        --threads $threads \
        --ref $ref \
        --hisat_ref $hisat_ref \
        --gene_bed $gene_bed \
        --output_dir $base_dir
        
done < "$samples_file"
```

## Example Usage (WGS+RNA)


### Step 1: RNA-only Inversion Candidate Detection
This step detects RNA-supported intragenic inversions
and generates *inv_in_gene_dedup_highcof.txt.

```
conda activate invas_preprocess

bash scripts/full_pipe/preprocess/preprocess_onlyrna.sh \
    inv_input.bed \
    hg38_genes.bed \
    still_unmap_bwa.s.bam \
    hisat.map_unmapremap.s.bam \
    candidate_out/rna_only/sample1

```

### Step 1 Outputs

```
candidate_out/rna_only/sample1/
├── chr1_inv_in_gene_dedup_highcof.txt
├── chr2_inv_in_gene_dedup_highcof.txt
├── ...
```
Optionally merge all chromosomes:

```
cat candidate_out/rna_only/sample1/chr*_inv_in_gene_dedup_highcof.txt \
    > candidate_out/rna_only/sample1_inv_in_gene_dedup_highcof.txt
```

### Step 3: Assembly

```
conda activate invas_assembly

mkdir -p assembly_out/sample1

python scripts/full_pipe/main_rna.py \
    --input_file candidate_out/rna_only/sample1_inv_in_gene_dedup_highcof.txt \
    --sample_name sample1 \
    --output_dir assembly_out/sample1 \
    --ref_genome hg38.fa \
    --ref_version hg38 \
    --rna_map_with_remap_bam rna_out/sample1/hisat.map_unmapremap.s.bam \
    --rna_unmap_bwa_bam rna_out/sample1/still_unmap_bwa.s.bam \
    --with_normal_haps
```



## When to Use Each Mode

INVAS provides two running modes depending on **data availability and quality**.

### ✅ WGS + RNA Mode (Recommended)

Use **WGS+RNA mode** when:

- ✅ **High-quality WGS data is available**
- ✅ **Matched RNA-seq data is available**
- ✅ You want **high-confidence inversion breakpoint detection**
- ✅ You want **maximum accuracy in haplotype reconstruction**

This mode integrates:
- Precise breakpoint evidence from **WGS**
- Transcript-level support from **RNA-seq**

➡️ **This is the recommended mode whenever good WGS data is available.**

---

### ✅ RNA-only Mode

Use **RNA-only mode** when:

- ✅ **Only RNA-seq data is available**
- ✅ **WGS data is missing or of low quality**
- ✅ WGS coverage is insufficient for reliable breakpoint calling

In this mode:
- Inversion candidates are inferred ** from RNA-seq  and InvasDB evidence**

- Final output format is **identical** to WGS+RNA mode



---

### Summary

| Data Condition | Recommended Mode |
|---------------|------------------|
| High-quality WGS + RNA | ✅ WGS + RNA |
| RNA only | ✅ RNA-only |
| Low-quality WGS + RNA | ✅ RNA-only |


## System Requirements

INVAS supports two running modes with **very different resource requirements**.

---

### ✅ WGS + RNA Mode

Resource usage in WGS+RNA mode **depends strongly on WGS sequencing depth
and genome coverage**.

- **Memory**
  - ~30–60 GB RAM for typical WGS (30×)
  - ≥100 GB RAM only required for **very deep WGS** or large batch runs
- **CPU**
  - Multi-core processor recommended (16+ cores for reasonable runtime)
- **Storage**
  - ~50 GB per sample for intermediate files (WGS + RNA)
  - ~20 GB for reference genome and indices
- **Operating System**
  - Linux (tested on Ubuntu 20.04, CentOS 7)

➡️ **Recommended when high-quality WGS data is available**, and sufficient
memory can be allocated.

---

### ✅ RNA-only Mode

RNA-only mode is **lightweight** and suitable for machines with limited resources.

- **Memory**
  - Tens to hundreds of MB RAM for whole-transcriptome runs
  - Even lower when using `--target_gene`
- **CPU**
  - 4–8 cores are sufficient
- **Storage**
  - A few GB per sample for RNA-seq intermediates
- **Operating System**
  - Linux (tested on Ubuntu 20.04, CentOS 7)

➡️ Ideal when **only RNA-seq data is available** or WGS quality is insufficient.

---

### ✅ Notes

- Using `--target_gene` significantly reduces memory and runtime in both modes
- RNA-only mode is recommended for quick validation and exploratory analysis
## Troubleshooting

Common issues and solutions:

1. **Memory errors during assembly**
   - Increase memory allocation in SLURM script
   - Reduce `--extend_region` parameter

2. **Missing SV calls**
   - Ensure all SV callers completed successfully
   - Check BAM file integrity with `samtools quickcheck`

3. **No candidate inversions found**
   - Verify gene annotation file matches reference genome version
   - Check SV caller outputs for detected variants



## License

This software is released under the [MIT] License. See LICENSE file for details.

## Contact

For questions, bug reports, or feature requests:
- Email: [xuedowang2-c@my.cityu.edu.hk]
- GitHub Issues: [https://github.com/deepomicslab/INVAS/issues/new]

## Acknowledgments

Invas incorporates several published tools:
- HISAT2 for RNA-seq alignment
- Manta, Delly, Lumpy, and Svaba for SV detection
- Picard for BAM file processing
