#!/bin/bash
#SBATCH --job-name='align_rg'
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=log.align.log
#SBATCH --mem=5G
#SBATCH --time=14-00:00:00

IDX=/data6/wangxuedong/invas_nc/wgs_ref/hg19
R1=merged_r1.fastq.gz   
R2=merged_r2.fastq.gz   
OUT=hisat2_hg19

hisat2 -p 1 -x "$IDX" \
  -1 "$R1" -2 "$R2" \
  -S "${OUT}.sam" \
  2> "${OUT}.log"


samtools view -@ 16 -bS "${OUT}.sam" | samtools sort -@ 16 -o "${OUT}.sorted.bam"
samtools index "${OUT}.sorted.bam"
rm -f "${OUT}.sam"