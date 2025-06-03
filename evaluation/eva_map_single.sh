#!/bin/bash
# read_mapping_evaluation_with_clip.sh
#
# Usage:
#   ./read_mapping_evaluation_with_clip.sh <fq1> <fq2> <assembly.fasta> [output_metrics.txt]
#
# Description:
#   This script maps paired-end FASTQ files to the assembled transcriptome FASTA file,
#   calculates basic mapping metrics including:
#     - Total reads (from flagstat)
#     - Mapped reads count
#     - Mapping rate (%)
#     - Average coverage depth
#   and additional alignment quality metrics:
#     - Total alignments (records in BAM)
#     - Fully aligned reads: alignments whose CIGAR string consists of only matches (e.g., "100M"),
#       indicating that the entire read is aligned without any clipping.
#     - Reads with clipping: alignments that include soft or hard clipping (indicated by 'S' or 'H' in the CIGAR string)
#
# Requires:
#   - BWA
#   - SAMtools

if [ $# -lt 3 ]; then
    echo "Usage: $0 <fq1> <fq2> <assembly.fasta> [output_metrics.txt]"
    exit 1
fi

# Input parameters
fq1=$1
fq2=$2
assembly=$3
output=${4:-"mapping_metrics.txt"}

# Index the assembly FASTA file if not already indexed
if [ ! -f "${assembly}.bwt" ]; then
    echo "Indexing the assembly: $assembly ..."
    bwa index "$assembly"
fi

# 1. Map the reads using BWA-MEM to generate a SAM file
echo "Mapping reads with BWA MEM..."
bwa mem "$assembly" "$fq1" "$fq2" > mapping.sam

# 2. Convert SAM to BAM and sort the BAM file
echo "Converting SAM to sorted BAM file..."
samtools view -bS mapping.sam | samtools sort -o mapping.sorted.bam

# 3. Index the sorted BAM file
samtools index mapping.sorted.bam

# 4. Use samtools flagstat to extract overall mapping statistics
echo "Extracting mapping statistics with samtools flagstat..."
flagstat=$(samtools flagstat mapping.sorted.bam)

total_reads=$(echo "$flagstat" | grep "in total" | head -n1 | awk '{print $1}')
mapped_reads=$(echo "$flagstat" | grep "mapped (" | head -n1 | awk '{print $1}')
mapping_rate=$(echo "scale=2; ($mapped_reads/$total_reads)*100" | bc)

# 5. Calculate the average coverage using samtools depth
avg_cov=$(samtools depth -a mapping.sorted.bam | awk '{sum+=$3; count++} END {if(count>0) printf "%.2f", sum/count; else print 0}')

# 6. Evaluate alignment quality by checking for clipping in CIGAR strings
#    - Fully aligned reads: CIGAR string containing only matches (e.g., "100M").
#    - Reads with clipping: those that contain soft clip (S) or hard clip (H).
total_alignments=$(samtools view -c mapping.sorted.bam)
fully_aligned=$(samtools view mapping.sorted.bam | awk '$6 ~ /^[0-9]+M$/ {count++} END {print count}')
clipped=$(($total_alignments - $fully_aligned))

fully_aligned_rate=$(echo "scale=2; ($fully_aligned/$total_alignments)*100" | bc)
clipped_rate=$(echo "scale=2; ($clipped/$total_alignments)*100" | bc)

# 7. Output all metrics to the output file
{
  echo "Mapping Metrics for Assembly: $assembly"
  echo "----------------------------------------"
  echo "Total reads (flagstat):       $total_reads"
  echo "Mapped reads:                 $mapped_reads"
  echo "Mapping rate:                 $mapping_rate %"
  echo "Average coverage:             $avg_cov"
  echo ""
  echo "Alignment Clipping Analysis:"
  echo "----------------------------"
  echo "Total alignments:             $total_alignments"
  echo "Fully aligned reads:          $fully_aligned"
  echo "Reads with clipping:          $clipped"
  echo "Fully aligned rate:           $fully_aligned_rate %"
  echo "Clipped reads rate:           $clipped_rate %"
} > "$output"

echo "Mapping metrics (including clipping analysis) have been saved to $output."

# 8. Cleanup: Remove the intermediate SAM file (BAM file is kept for further analysis if needed)
rm mapping.sam

echo "Evaluation complete."