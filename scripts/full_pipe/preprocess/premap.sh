#!/bin/sh
#SBATCH --partition=batch
#SBATCH --cpus-per-task=128
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=400G
#SBATCH --time=03-00:00:00
#SBATCH --job-name=test_realign
#SBATCH --output=log.align.log



raw_bam=$1
bamUtil=bam
# $bamUtil bam2FastQ \
#     --in $raw_bam \
#     --firstOut hisat_1.fastq \
#     --secondOut hisat_2.fastq \
#     --unpairedOut hisat_unpair.fastq &> /dev/null
threads=128

samtools sort -n -@ $threads -O BAM -o raw_bam.sortn.bam $raw_bam
samtools fastq -0 hisat_unpair.fastq -1 hisat_1.fastq -2 hisat_2.fastq raw_bam.sortn.bam -s hisat_single.fastq -N -@ $threads


hisat2 -x /home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla -1 hisat_1.fastq -2 hisat_2.fastq -S hisat2.realign.sam -p $threads
samtools sort -@ $threads -O BAM -o hisat2.realign.bam hisat2.realign.sam
samtools index hisat2.realign.bam


in_bam=hisat2.realign.bam
unmap=hisat2.unmap.bam
map=hisat2.map.bam
samtools view -f4 $in_bam -O BAM -o $unmap
samtools view -F4 $in_bam -O BAM -o $map
samtools sort -n -@ $threads -O BAM -o hisat2.unmap.sortn.bam $unmap
samtools fastq -0 hisat_unmap_unpaired.fastq -1 hisat_unmap_1.fastq -2 hisat_unmap_2.fastq hisat2.unmap.sortn.bam -s hisat_unmap_single.fastq -N -@ $threads
ref=/home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
# $bamUtil bam2FastQ \
#     --in $unmap \
#     --firstOut hisat_unmap_1.fastq \
#     --secondOut hisat_unmap_2.fastq \
#     --unpairedOut hisat_unmap_unpair.fastq &> /dev/null


hisat2 -x /home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla -U hisat_unmap_1.fastq -S hisat.unmap1.remap.sam -p $threads
hisat2 -x /home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla -U hisat_unmap_2.fastq -S hisat.unmap2.remap.sam -p $threads
hisat2 -x /home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla -U hisat_unmap_single.fastq -S hisat.unmaps.remap.sam -p $threads
samtools merge hisat.unmap_remap_all.bam hisat.unmap1.remap.sam hisat.unmap2.remap.sam hisat.unmaps.remap.sam 
samtools sort hisat.unmap_remap_all.bam -@ $threads -O BAM -o hisat.unmap_remap_all.s.bam
samtools view -F 4 hisat.unmap_remap_all.s.bam -O BAM -@ $threads -o hisat.unmap_remap_all.s.map.bam
samtools merge hisat.map_unmapremap.bam $map hisat.unmap_remap_all.s.map.bam
samtools sort hisat.map_unmapremap.bam -@ $threads -O BAM -o hisat.map_unmapremap.s.bam


samtools view -f4 hisat.unmap1.remap.sam -O BAM -o hisat.unmap1.stillunmap.bam
samtools view -f4 hisat.unmap2.remap.sam -O BAM -o hisat.unmap2.stillunmap.bam
samtools view -f4 hisat.unmaps.remap.sam -O BAM -o hisat.unmaps.stillunmap.bam

java -jar /home/xuedowang2/app/picard.jar \
  SamToFastq \
  -I hisat.unmap1.stillunmap.bam \
  -F hisat.unmap1.stillunmap.fq


java -jar /home/xuedowang2/app/picard.jar \
  SamToFastq \
  -I hisat.unmap2.stillunmap.bam \
  -F hisat.unmap2.stillunmap.fq


java -jar /home/xuedowang2/app/picard.jar \
  SamToFastq \
  -I hisat.unmaps.stillunmap.bam \
  -F hisat.unmaps.stillunmap.fq

awk '{if(NR%4==1) {print $0 "_r1"} else {print $0}}' hisat.unmap1.stillunmap.fq > hisat.unmap1.stillunmap.r1.fq
awk '{if(NR%4==1) {print $0 "_r2"} else {print $0}}' hisat.unmap2.stillunmap.fq > hisat.unmap2.stillunmap.r2.fq
awk '{if(NR%4==1) {print $0 "_sg"} else {print $0}}' hisat.unmaps.stillunmap.fq > hisat.unmaps.stillunmap.sg.fq
cat hisat.unmap1.stillunmap.r1.fq hisat.unmap2.stillunmap.r2.fq hisat.unmaps.stillunmap.sg.fq > still_unmap_r1r2sg.fq
bwa mem $ref still_unmap_r1r2sg.fq > still_unmap_bwa.sam
samtools sort still_unmap_bwa.sam -O BAM -o still_unmap_bwa.s.bam
samtools index still_unmap_bwa.s.bam