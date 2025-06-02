#!/bin/sh
#SBATCH --partition=batch
#SBATCH --cpus-per-task=128
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --time=03-00:00:00
#SBATCH --job-name=test_realign
#SBATCH --output=log.align.log


threads=128
bamUtil=bam


fq1_1=../fq/hg002_gm24385.mrna.R1.fastq.gz
fq1_2=../fq/hg002_gm24385.mrna.R2.fastq.gz
fq2_1=../fq/hg002_gm26105.mrna.R1.fastq.gz
fq2_2=../fq/hg002_gm26105.mrna.R2.fastq.gz
fq3_1=../fq/hg002_gm27730.mrna.R1.fastq.gz
fq3_2=../fq/hg002_gm27730.mrna.R2.fastq.gz

fq1_1_clean=../fq/hg002_gm24385.mrna.R1.clean.fastq.gz
fq1_2_clean=../fq/hg002_gm24385.mrna.R2.clean.fastq.gz
fq2_1_clean=../fq/hg002_gm26105.mrna.R1.clean.fastq.gz
fq2_2_clean=../fq/hg002_gm26105.mrna.R2.clean.fastq.gz
fq3_1_clean=../fq/hg002_gm27730.mrna.R1.clean.fastq.gz
fq3_2_clean=../fq/hg002_gm27730.mrna.R2.clean.fastq.gz

fq1_1_unpair=../fq/hg002_gm24385.mrna.R1.unpair.fastq.gz
fq1_2_unpair=../fq/hg002_gm24385.mrna.R2.unpair.fastq.gz
fq2_1_unpair=../fq/hg002_gm26105.mrna.R1.unpair.fastq.gz
fq2_2_unpair=../fq/hg002_gm26105.mrna.R2.unpair.fastq.gz
fq3_1_unpair=../fq/hg002_gm27730.mrna.R1.unpair.fastq.gz
fq3_2_unpair=../fq/hg002_gm27730.mrna.R2.unpair.fastq.gz


hisat_ref=/home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla 
# activate conda
# eval "$(conda shell.bash hook)"
# conda activate /home/xuedowang2/app/spechla_run/SpecHLA/spechla_env

# trimmomatic PE -threads $threads \
#   $fq1_1 $fq1_2 \
#   $fq1_1_clean $fq1_1_unpair \
#   $fq1_2_clean $fq1_2_unpair \
#   ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

# trimmomatic PE -threads $threads \
#   $fq2_1 $fq2_2 \
#   $fq2_1_clean $fq2_1_unpair \
#   $fq2_2_clean $fq2_2_unpair \
#   ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

# trimmomatic PE -threads $threads \
#   $fq3_1 $fq3_2 \
#   $fq3_1_clean $fq3_1_unpair \
#   $fq3_2_clean $fq3_2_unpair \
#   ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
  
# conda deactivate
# trimmomatic PE -phred33 $fq1 $fq2  -threads $threads clean_1_paired.fastq clean_1_unpaired.fastq clean_2_paired.fastq clean_2_unpaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


# hisat2 -x /home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla \
#   -1 $fq1_1_clean,$fq2_1_clean,$fq3_1_clean \
#   -2 $fq1_2_clean,$fq2_2_clean,$fq3_2_clean \
#   -S hisat2.realign.sam \
#   -p $threads
hisat2 -x /home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla \
  -1 $fq1_1_clean \
  -2 $fq1_2_clean \
  -S hisat2.realign.1.sam \
  -p $threads

hisat2 -x /home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla \
  -1 $fq2_1_clean \
  -2 $fq2_2_clean \
  -S hisat2.realign.2.sam \
  -p $threads

hisat2 -x /home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla \
  -1 $fq3_1_clean \
  -2 $fq3_2_clean \
  -S hisat2.realign.3.sam \
  -p $threads

samtools sort -@ threads -O BAM -o hisat2.realign.1.bam hisat2.realign.1.sam
samtools index hisat2.realign.1.bam

samtools sort -@ threads -O BAM -o hisat2.realign.2.bam hisat2.realign.2.sam
samtools index hisat2.realign.2.bam

samtools sort -@ threads -O BAM -o hisat2.realign.3.bam hisat2.realign.3.sam
samtools index hisat2.realign.3.bam

samtools merge hisat2.realign.bam hisat2.realign.1.bam hisat2.realign.2.bam hisat2.realign.3.bam
samtools index -@ $threads hisat2.realign.bam


# hisat2 -x /home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla -1 hisat_1.fastq -2 hisat_2.fastq -S hisat2.realign.sam -p 32
samtools sort -@ threads -O BAM -o hisat2.realign.bam hisat2.realign.sam
samtools index hisat2.realign.bam


in_bam=hisat2.realign.bam
unmap=hisat2.unmap.bam
map=hisat2.map.bam
samtools view -f4 $in_bam -O BAM -o $unmap
samtools view -F4 $in_bam -O BAM -o $map
samtools sort -n -@ $threads -O BAM -o hisat2.unmap.sortn.bam $unmap
samtools fastq -0 hisat_unmap_unpaired.fastq -1 hisat_unmap_1.fastq -2 hisat_unmap_2.fastq hisat2.unmap.sortn.bam -s hisat_unmap_single.fastq -N -@$threads
ref=/home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
# $bamUtil bam2FastQ \
#     --in $unmap \
#     --firstOut hisat_unmap_1.fastq \
#     --secondOut hisat_unmap_2.fastq \
#     --unpairedOut hisat_unmap_unpair.fastq


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