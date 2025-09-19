#!/bin/bash
#SBATCH --job-name='align_rg'
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --output=log.align_rg.log
#SBATCH --mem=2G
#SBATCH --time=14-00:00:00
set -euo pipefail

REF='../../../../wgs_ref/hg19.fa'
R1='inv_wgs_R1.fq'
R2='inv_wgs_R2.fq'
THREADS=8

# Read Group string (edit fields if you like)
RG='@RG\tID:sim1\tSM:sample1\tPL:ILLUMINA\tLB:lib1\tPU:unit1'

# Sanity checks
for f in "$REF" "$R1" "$R2"; do
  if [ ! -s "$f" ]; then
    echo "ERROR: Missing or empty file: $f" >&2
    exit 1
  fi
done

# Ensure BWA index exists
[ -s "${REF}.bwt" ] || bwa index "$REF"

# 1) Align to SAM with RG (capture stderr separately)
echo "Running bwa mem with RG -> SAM ..."
bwa mem -t "$THREADS" -R "$RG" "$REF" "$R1" "$R2" > inv_wgs.hg19.sam 2> bwa.log

# Optional: sanity check the header and RG lines
echo "SAM header (first lines):"
head -n 10 inv_wgs.hg19.sam | sed -n '1,10p'
echo "RG lines in SAM header:"
grep -m 5 '^@RG' inv_wgs.hg19.sam || echo "No @RG found!"

# 2) Sort SAM to BAM
echo "Sorting ..."
samtools sort -o inv_wgs.hg19.bam inv_wgs.hg19.sam

# 3) Index and basic stats
samtools index inv_wgs.hg19.bam
samtools flagstat inv_wgs.hg19.bam

# 4) Verify RG present in BAM header
echo "RG lines in BAM header:"
samtools view -H inv_wgs.hg19.bam | grep '^@RG' || echo "No @RG found in BAM!"

# 5) (Optional) create LUMPY helper BAMs
# Discordant pairs
samtools view -b -F 1294 inv_wgs.hg19.bam > inv_wgs.discordants.bam
# Split reads (requires extractSplitReads_BwaMem in PATH)
if command -v extractSplitReads_BwaMem >/dev/null 2>&1; then
  samtools view -h inv_wgs.hg19.bam \
    | extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - > inv_wgs.splitters.bam
else
  echo "Warning: extractSplitReads_BwaMem not found in PATH; skipping splitters." >&2
fi

echo "Done. Outputs:"
echo " - inv_wgs.hg19.bam (+ .bai)"
echo " - inv_wgs.discordants.bam (optional)"
echo " - inv_wgs.splitters.bam (optional)"
echo " - bwa.log"