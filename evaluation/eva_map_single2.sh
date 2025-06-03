#!/bin/bash
# read_mapping_evaluation_with_clip.sh
#
# Usage:
#   ./read_mapping_evaluation_with_clip.sh <fq1> <fq2> <assembly.fasta> <output_folder>
#
# Description:
#   该脚本将 RNA 测序 reads 比对到指定组装的转录本上，
#   计算映射统计指标（如 samtools flagstat 输出的总数、映射数、映射率、平均覆盖度），
#   并通过解析 BAM 文件中的 CIGAR 字符串统计完全比对与剪切(read clipping)情况，
#   最后将各项指标以 TSV 格式写入到指定输出文件夹下的 mapping_metrics.tsv 文件中。
#
#   在处理过程中产生的中间文件（包括 SAM、未排序 BAM、排序后的 BAM 及其索引）
#   将在脚本退出时自动删除。
#
# Requires:
#   - BWA
#   - SAMtools
#   - bc
#
# 设置严格模式
set -euo pipefail

# 检查输入参数个数
if [ $# -lt 4 ]; then
    echo "Usage: $0 <fq1> <fq2> <assembly.fasta> <output_folder>"
    exit 1
fi

# 输入参数
fq1="$1"
fq2="$2"
assembly="$3"
output_folder="$4"

# 配置线程数（可根据实际环境调整）
threads=4

# 创建输出文件夹（如果不存在的话）
mkdir -p "$output_folder"
output_file="${output_folder}/mapping_metrics.tsv"

# 如果未建立参考索引，则执行索引
if [ ! -f "${assembly}.bwt" ]; then
    echo "Indexing the assembly: $assembly ..."
    bwa index "$assembly"
fi

echo "=============================="
echo "Mapping reads with BWA MEM..."
echo "=============================="

# 创建临时文件以存储中间结果
sam_file=$(mktemp --suffix=.sam)
bam_file=$(mktemp --suffix=.bam)
sorted_bam="mapping.sorted.bam"

# 设置 trap，在脚本退出时自动删除中间文件和索引文件
trap 'rm -f "$sam_file" "$bam_file" "$sorted_bam" "${sorted_bam}.bai"' EXIT

# 1. 利用 BWA MEM 将 reads 比对到参考上，生成 SAM 文件
bwa mem -t "$threads" "$assembly" "$fq1" "$fq2" > "$sam_file"

# 2. 将 SAM 文件转换为 BAM 格式
samtools view -bS "$sam_file" > "$bam_file"

# 3. 对 BAM 文件进行排序，生成最终的 BAM 文件
samtools sort -@ "$threads" "$bam_file" -o "$sorted_bam"

echo "Indexing sorted BAM file..."
samtools index "$sorted_bam"

echo "=============================================="
echo "Extracting mapping statistics using samtools flagstat..."
echo "=============================================="
flagstat=$(samtools flagstat "$sorted_bam")
total_reads=$(echo "$flagstat" | grep "in total" | head -n1 | awk '{print $1}')
mapped_reads=$(echo "$flagstat" | grep "mapped (" | head -n1 | awk '{print $1}')
mapping_rate=$(echo "scale=2; (100*$mapped_reads)/$total_reads" | bc)

echo "=================================="
echo "Calculating average coverage..."
echo "=================================="
# samtools depth 输出所有位置的覆盖度，然后计算平均覆盖度
avg_cov=$(samtools depth -a "$sorted_bam" | awk '{sum+=$3; count++} END {if(count>0) printf "%.2f", sum/count; else print 0}')

echo "==============================================="
echo "Analyzing CIGAR strings for clipping analysis..."
echo "==============================================="
# 统计总比对记录数
total_alignments=$(samtools view -c "$sorted_bam")
# 完全比对的记录（CIGAR 字段仅包含数字和 M）
fully_aligned=$(samtools view "$sorted_bam" | awk '$6 ~ /^[0-9]+M$/ {count++} END {print count+0}')
# 剪切的记录数
clipped=$(echo "$total_alignments - $fully_aligned" | bc)
fully_aligned_rate=$(echo "scale=2; (100*$fully_aligned)/$total_alignments" | bc)
clipped_rate=$(echo "scale=2; (100*$clipped)/$total_alignments" | bc)

echo "============================="
echo "Saving mapping metrics to TSV file..."
echo "============================="
{
  echo -e "Metric\tValue"
  echo -e "Assembly\t$assembly"
  echo -e "Total_reads_flagstat\t$total_reads"
  echo -e "Mapped_reads\t$mapped_reads"
  echo -e "Mapping_rate(%)\t$mapping_rate"
  echo -e "Average_coverage\t$avg_cov"
  echo -e "Total_alignments\t$total_alignments"
  echo -e "Fully_aligned_reads\t$fully_aligned"
  echo -e "Fully_aligned_rate(%)\t$fully_aligned_rate"
  echo -e "Clipped_reads\t$clipped"
  echo -e "Clipped_reads_rate(%)\t$clipped_rate"
} > "$output_file"

echo "Mapping metrics have been saved to $output_file."
echo "Evaluation complete."