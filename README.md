# Invas

Invas is a conjugate graph-based RNA assembler designed for the precise detection of intragenomic inversions and transcriptome assembly. It leverages an iterative conjugate maximum flow network model and integrates RNA sequencing (RNA-Seq) and whole-genome sequencing (WGS) data.

## Overview

Invas accurately detects expressed intragenomic inversions and reconstructs both normal and inverted transcripts through a comprehensive workflow that combines structural variant detection with transcriptome assembly. The tool processes RNA-seq and WGS data through quality control, alignment, structural variant detection, and iterative conjugate maximum-flow algorithm for precise transcriptome reconstruction.

## Features

- Detection of intragenomic inversions from RNA-seq and WGS data
- Integration of multiple SV callers (Delly, Manta, Lumpy, Svaba)
- Refinement of inverted segment splicing positions
- Conjugate graph-based transcriptome assembly
- Support for both normal and inverted transcript reconstruction
- Re-alignment of unmapped RNA-seq reads for improved sensitivity
- Iterative conjugate maximum flow network model

## Installation

### Prerequisites

- Conda/Mamba package manager
- Reference genome files (hg19 or hg38)
- HISAT2 index files
- Gene annotation files (GENCODE)

### Setup

1. Clone the repository:
```bash
git clone git@github.com:deepomicslab/INVAS.git
cd invas
```

2. Create conda environments:
```bash
# Data preparation environment
conda env create -f data_prepare.yaml -n invas_data_prep

# Assembly environment
conda env create -f assembly.yaml -n invas_assembly
```


## Workflow

Invas workflow consists of six steps, implemented through four main processing stages:

### Step 1: RNA-seq Data Processing

Process RNA-seq data for quality control, alignment, and unmapped read extraction.

```bash
conda activate invas_data_prep

# Create output directory
mkdir -p rna_out/$sample_name

# Run RNA processing
bash scripts/full_pipe/preprocess/test_common_run.sh -n $sample_name \
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

### Step 2: WGS Data Processing and SV Detection

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

### Step 3: Candidate Intragenic Inversion Event Detection

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

### Step 4: Transcriptome Assembly

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

## Example Usage

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

# Step 1: RNA-seq processing
echo "Step 1: Processing RNA-seq data..."
conda activate invas_data_prep

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

# Step 2: WGS SV detection
echo "Step 2: Detecting structural variants..."
conda activate manta_sv

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

# Step 3: Candidate detection
echo "Step 3: Detecting candidate inversions..."
mkdir -p $candidate_dir

bash scripts/full_pipe/preprocess/combine_sv_rna.sh $sample_name \
    $sv_dir/$sample_name \
    $rna_dir/$sample_name/hisat.map_unmapremap.s.bam \
    $rna_dir/$sample_name/still_unmap_bwa.s.bam \
    $candidate_dir \
    $gene_bed \
    delly,manta,lumpy,svaba

# Step 4: Transcriptome assembly
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

## System Requirements

- **Memory**: Minimum 100GB RAM (more for larger datasets)
- **CPU**: Multi-core processor (16+ cores recommended)
- **Storage**: 
  - ~500GB for intermediate files per sample
  - ~100GB for reference files and indices
- **Operating System**: Linux (tested on Ubuntu 20.04, CentOS 7)

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

## Citation

If you use Invas in your research, please cite:
```
[Your paper citation here]
```

## License

This software is released under the [MIT/GPL/Apache] License. See LICENSE file for details.

## Contact

For questions, bug reports, or feature requests:
- Email: [xuedowang2-c@my.cityu.edu.hk]
- GitHub Issues: [https://github.com/deepomicslab/INVAS/issues/new]

## Acknowledgments

Invas incorporates several published tools:
- HISAT2 for RNA-seq alignment
- Manta, Delly, Lumpy, and Svaba for SV detection
- Picard for BAM file processing