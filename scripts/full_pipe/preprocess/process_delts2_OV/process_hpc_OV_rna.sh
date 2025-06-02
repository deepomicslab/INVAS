#!/bin/sh
#SBATCH --partition=tiny
#SBATCH --cpus-per-task=64
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=05-00:00:00
#SBATCH --job-name=test_realign
#SBATCH --output=log.align.log

threads=64

eval "$(conda shell.bash hook)"

outdir=/scratch/project/cs_shuaicli/wxd/OV/RNA/rna_out
ref=/scratch/project/cs_shuaicli/wxd/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
hisat_ref=/scratch/project/cs_shuaicli/wxd/ref/GRCh38_full_analysis_set_plus_decoy_hla
samples_file=/scratch/project/cs_shuaicli/wxd/OV/RNA/samples.txt
for sample in `cat $samples_file`; do
    echo "Processing $sample ..."
    # conda activate /home/wangxuedong/app/spechla_run/SpecHLA/spechla_env
    # fq1_var=${sample}_fq1
    # fq2_var=${sample}_fq2
    fq1=/scratch/project/cs_shuaicli/wxd/OV/RNA/fqs/${sample}.R1.fastq.gz
    fq2=/scratch/project/cs_shuaicli/wxd/OV/RNA/fqs/${sample}.R2.fastq.gz
    sample_outdir=$outdir/$sample
    if [ ! -d "$sample_outdir" ]; then
        mkdir "$sample_outdir"
    fi
    fq1_clean=$sample_outdir/${sample}_1.clean.fastq.gz
    fq2_clean=$sample_outdir/${sample}_2.clean.fastq.gz
    fq1_unpair=$sample_outdir/${sample}_1.unpair.fastq.gz
    fq2_unpair=$sample_outdir/${sample}_2.unpair.fastq.gz
    

    trimmomatic PE -threads $threads \
    $fq1 $fq2 \
    $fq1_clean $fq1_unpair \
    $fq2_clean $fq2_unpair \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36


    hisat2 -x $hisat_ref \
    -1 $fq1_clean \
    -2 $fq2_clean \
    -S $sample_outdir/hisat2.sam \
    -p $threads

    samtools sort -@ $threads -O BAM -o $sample_outdir/hisat2.bam $sample_outdir/hisat2.sam
    samtools index -@ $threads $sample_outdir/hisat2.bam

    in_bam=$sample_outdir/hisat2.bam
    unmap=$sample_outdir/hisat2.unmap.bam
    map=$sample_outdir/hisat2.map.bam
    unmap_sortn=$sample_outdir/hisat2.unmap.sortn.bam

    samtools view -f4 $in_bam -O BAM -o $unmap
    samtools view -F4 $in_bam -O BAM -o $map
    samtools sort -n -@ $threads -O BAM -o $unmap_sortn $unmap

    unmap_fq_up=$sample_outdir/hisat_unmap_unpaired.fastq
    unmap_fq_1=$sample_outdir/hisat_unmap_1.fastq
    unmap_fq_2=$sample_outdir/hisat_unmap_2.fastq
    unmap_fq_s=$sample_outdir/hisat_unmap_single.fastq
    samtools fastq -0 $unmap_fq_up -1 $unmap_fq_1 -2 $unmap_fq_2 $unmap_sortn -s $unmap_fq_s -N -@ $threads

    hisat2 -x $hisat_ref -U $unmap_fq_1 -S $sample_outdir/hisat.unmap1.remap.sam -p $threads
    hisat2 -x $hisat_ref -U $unmap_fq_2 -S $sample_outdir/hisat.unmap2.remap.sam -p $threads
    hisat2 -x $hisat_ref -U $unmap_fq_s -S $sample_outdir/hisat.unmapup.remap.sam -p $threads
    samtools merge -f \
    $sample_outdir/hisat.unmap_remap_all.bam \
    $sample_outdir/hisat.unmap1.remap.sam \
    $sample_outdir/hisat.unmap2.remap.sam \
    $sample_outdir/hisat.unmapup.remap.sam

    samtools sort -@ $threads -O BAM -o $sample_outdir/hisat.unmap_remap_all.s.bam $sample_outdir/hisat.unmap_remap_all.bam
    samtools view -F 4 -O BAM -@ $threads -o $sample_outdir/hisat.unmap_remap_all.s.map.bam $sample_outdir/hisat.unmap_remap_all.s.bam
    samtools merge -f $sample_outdir/hisat.map_unmapremap.bam $map $sample_outdir/hisat.unmap_remap_all.s.map.bam
    samtools sort -@ $threads -O BAM -o $sample_outdir/hisat.map_unmapremap.s.bam $sample_outdir/hisat.map_unmapremap.bam

    samtools view -f4 $sample_outdir/hisat.unmap1.remap.sam -O BAM -o $sample_outdir/hisat.unmap1.stillunmap.bam
    samtools view -f4 $sample_outdir/hisat.unmap2.remap.sam -O BAM -o $sample_outdir/hisat.unmap2.stillunmap.bam
    samtools view -f4 $sample_outdir/hisat.unmapup.remap.sam -O BAM -o $sample_outdir/hisat.unmapup.stillunmap.bam

    java -jar /scratch/project/cs_shuaicli/wxd/app/picard.jar \
    SamToFastq \
    -I $sample_outdir/hisat.unmap1.stillunmap.bam \
    -F $sample_outdir/hisat.unmap1.stillunmap.fq

    java -jar /scratch/project/cs_shuaicli/wxd/app/picard.jar \
    SamToFastq \
    -I $sample_outdir/hisat.unmap2.stillunmap.bam \
    -F $sample_outdir/hisat.unmap2.stillunmap.fq

    java -jar /scratch/project/cs_shuaicli/wxd/app/picard.jar \
    SamToFastq \
    -I $sample_outdir/hisat.unmapup.stillunmap.bam \
    -F $sample_outdir/hisat.unmapup.stillunmap.fq

    awk '{if(NR%4==1) {print $0 "_r1"} else {print $0}}' $sample_outdir/hisat.unmap1.stillunmap.fq > $sample_outdir/hisat.unmap1.stillunmap.r1.fq
    awk '{if(NR%4==1) {print $0 "_r2"} else {print $0}}' $sample_outdir/hisat.unmap2.stillunmap.fq > $sample_outdir/hisat.unmap2.stillunmap.r2.fq
    awk '{if(NR%4==1) {print $0 "_sg"} else {print $0}}' $sample_outdir/hisat.unmapup.stillunmap.fq > $sample_outdir/hisat.unmapup.stillunmap.sg.fq
    cat $sample_outdir/hisat.unmap1.stillunmap.r1.fq $sample_outdir/hisat.unmap2.stillunmap.r2.fq $sample_outdir/hisat.unmapup.stillunmap.sg.fq > $sample_outdir/still_unmap_r1r2sg.fq
    bwa mem $ref $sample_outdir/still_unmap_r1r2sg.fq -t $threads > $sample_outdir/still_unmap_bwa.sam
    samtools sort $sample_outdir/still_unmap_bwa.sam -O BAM -o $sample_outdir/still_unmap_bwa.s.bam
    samtools index $sample_outdir/still_unmap_bwa.s.bam


    
done