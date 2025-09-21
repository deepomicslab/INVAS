#!/bin/bash
#SBATCH --job-name='align_rg'
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=log.rna.log
#SBATCH --mem=10G
#SBATCH --time=14-00:00:00
# conda activate invas_assembly

# Create output directory
mkdir -p candidate_out

sample_name="test"
gene_annotation_bed="/data6/wangxuedong/invas_nc/INVAS/scripts/full_pipe/annotation/hg19_genes.gencode.bed"

# Run candidate detection
# Note: Added missing '\' after the script path for proper line continuation.
# Ensure there are no spaces after each '\'. This should fix the argument count error
# and prevent the subsequent "is a directory" error by ensuring all parts are treated as one command.
bash ../scripts/full_pipe/preprocess/combine_sv_rna.sh \
    $sample_name \
    sv_out/$sample_name/sv_out/$sample_name \
    rna_out/$sample_name/hisat.map_unmapremap.s.bam \
    rna_out/$sample_name/still_unmap_bwa.s.bam \
    candidate_out \
    $gene_annotation_bed \
    delly,manta
