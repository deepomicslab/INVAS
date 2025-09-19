#!/bin/bash
#SBATCH --job-name='align_rg'
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=log.rna.log
#SBATCH --mem=10G
#SBATCH --time=14-00:00:00
# conda activate assembly_env
sample_name="test"
# Create output directory
mkdir -p rna_out/$sample_name
threads=2
reference_genome="/data6/wangxuedong/invas_nc/wgs_ref/hg19.fa"
hisat_index="/data6/wangxuedong/invas_nc/wgs_ref/hg19"
input_bam=/data6/wangxuedong/invas_nc/INVAS/test/demo/EXONS/hisat2_hg19.sorted.bam
# Run RNA processing
bash ../scripts/full_pipe/preprocess/process_rna_common.sh -n $sample_name \
    -b $input_bam \
    -o rna_out/$sample_name \
    -r $reference_genome \
    -x $hisat_index \
    -t $threads \
    -p ../scripts/full_pipe/bin/picard.jar