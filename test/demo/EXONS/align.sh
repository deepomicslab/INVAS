#!/bin/bash
#SBATCH --job-name='align_rg'
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=log.align.log
#SBATCH --mem=5G
#SBATCH --time=14-00:00:00

IDX=/data6/wangxuedong/invas_nc/wgs_ref/hg19
R1=merged_R1.fastq.gz   # 或你的R1文件
R2=merged_R2.fastq.gz   # 或你的R2文件
OUT=hisat2_hg19

hisat2 -p 1 -x "$IDX" \
  -1 "$R1" -2 "$R2" \
  -S "${OUT}.sam" \
  2> "${OUT}.log"

# 转 BAM、排序、建索引
samtools view -@ 16 -bS "${OUT}.sam" | samtools sort -@ 16 -o "${OUT}.sorted.bam"
samtools index "${OUT}.sorted.bam"
rm -f "${OUT}.sam"