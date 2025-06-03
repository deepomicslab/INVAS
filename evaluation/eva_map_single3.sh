#!/bin/bash
# read_mapping_evaluation_with_clip.sh
#
# Usage:
#   ./read_mapping_evaluation_with_clip.sh <fq1> <fq2> <assembly.fasta> <output_folder>
#
# Description:
#   该脚本将 RNA 测序 reads 比对到指定组装的转录本上，
#   计算映射统计指标（如 samtools flagstat 输出的总数、映射数、映射率、平均覆盖度），
#   并通过解析 BAM 文件中的 CIGAR 字符串来统计完全比对与剪切 (clipped) 的情况，
#   最后将所有指标以 TSV 格式写入到输出目录下的 mapping_metrics.tsv 文件中。
#
#   中间步骤产生的文件（SAM 文件、未排序 BAM、排序后的 BAM 以及 BAM 索引）均保存在输出目录中，
#   脚本退出时会自动删除这些中间文件。
#
# Requires:
#   - BWA
#   - SAMtools
#   - bc
#
# 设置严格模式：一旦命令失败则退出，并禁止使用未定义的变量
set -euo pipefail

# 检查参数个数
if [ $# -lt 4 ]; then
    echo "Usage: $0 <fq1> <fq2> <assembly.fasta> <output_folder>"
    exit 1
fi

# 输入参数
fq1="$1"
fq2="$2"
assembly="$3"
output_folder="$4"

# 配置线程数（根据实际硬件情况调整）
threads=4

# 创建输出目录（如果不存在，则新建）
mkdir -p "$output_folder"
output_file="${output_folder}/mapping_metrics.tsv"

# 设置中间文件路径，全部输出到输出目录下
sam_file=$(mktemp "${output_folder}/aln.XXXXXX.sam")
bam_file=$(mktemp "${output_folder}/aln.XXXXXX.bam")
sorted_bam="${output_folder}/mapping.sorted.bam"

# 设置 trap，在脚本退出时自动删除中间文件（包括 SAM、BAM、排序后的 BAM 及其索引）
trap 'rm -f "$sam_file" "$bam_file" "$sorted_bam" "${sorted_bam}.bai"' EXIT

# 对参考转录本建立索引（如果未存在）
if [ ! -f "${assembly}.bwt" ]; then
    echo "Indexing the assembly: $assembly ..."
    bwa index "$assembly"
fi

echo "=============================="
echo "Mapping reads with BWA MEM..."
echo "=============================="

# 1. 使用 BWA MEM 进行比对，生成 SAM 文件
bwa mem -t "$threads" "$assembly" "$fq1" "$fq2" > "$sam_file"

# 2. 将 SAM 文件转换为 BAM 格式
samtools view -bS "$sam_file" > "$bam_file"

# 3. 对 BAM 文件进行排序
samtools sort -@ "$threads" "$bam_file" -o "$sorted_bam"

echo "Indexing sorted BAM file..."
samtools index "$sorted_bam"

echo "=============================================="
echo "Extracting mapping statistics with samtools flagstat..."
echo "=============================================="
flagstat=$(samtools flagstat "$sorted_bam")
total_reads=$(echo "$flagstat" | grep "in total" | head -n1 | awk '{print $1}')
mapped_reads=$(echo "$flagstat" | grep "mapped (" | head -n1 | awk '{print $1}')
mapping_rate=$(echo "scale=2; (100*$mapped_reads)/$total_reads" | bc)

echo "=================================="
echo "Calculating average coverage..."
echo "=================================="
# 计算所有位点的平均覆盖度
avg_cov=$(samtools depth -a "$sorted_bam" | awk '{sum+=$3; count++} END {if(count>0) printf "%.2f", sum/count; else print 0}')

echo "==============================================="
echo "Analyzing CIGAR strings for clipping analysis..."
echo "==============================================="
# 统计总比对记录数
total_alignments=$(samtools view -c "$sorted_bam")
# 统计完全比对的记录（CIGAR 字段仅包含数字和 M）
fully_aligned=$(samtools view "$sorted_bam" | awk '$6 ~ /^[0-9]+M$/ {count++} END {print count+0}')
# 剪切记录数（总比对数减去完全比对数）
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
echo "Intermediate files in the output directory will be deleted upon script exit."
echo "Evaluation complete."