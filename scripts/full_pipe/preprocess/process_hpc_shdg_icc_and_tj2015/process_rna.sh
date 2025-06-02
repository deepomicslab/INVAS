#!/bin/sh
#SBATCH --partition=batch
#SBATCH --cpus-per-task=128
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=03-00:00:00
#SBATCH --job-name=test_realign
#SBATCH --output=log.align.log

threads=128

eval "$(conda shell.bash hook)"

outdir=/scratch/project/cs_shuaicli/wxd/virus/ShangHai_GanDan_ICC/RNA/rna_out
ref=/home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
hisat_ref=/home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla
samples_file=/scratch/project/cs_shuaicli/wxd/virus/ShangHai_GanDan_ICC/samples.txt
for sample in `cat $samples_file`; do
    echo "Processing $sample ..."
    

    conda activate /home/xuedowang2/app/spechla_run/SpecHLA/spechla_env
    # fq1_var=${sample}_fq1
    # fq2_var=${sample}_fq2
    # fq1=/data6/yangtianxia/RNA/${sample}.R1.fastq.gz
    # fq2=/data6/yangtianxia/RNA/${sample}.R2.fastq.gz
    sample_outdir=$outdir/$sample
    if [ ! -d "$sample_outdir" ]; then
        mkdir "$sample_outdir"
    fi

    ori_in_bam=/scratch/project/cs_shuaicli/wxd/virus/ShangHai_GanDan_ICC/RNA/${sample}/${sample}.bam
    # extract fastq at first 
    samtools sort -n -@ $threads -O BAM -o $sample_outdir/${sample}.hisat2.sortn.bam $ori_in_bam
    samtools fastq \
        -0 $sample_outdir/${sample}.hisat2.ori.unpaired.fastq \
        -1 $sample_outdir/${sample}.hisat2.ori.1.fastq \
        -2 $sample_outdir/${sample}.hisat2.ori.2.fastq \
        -s $sample_outdir/${sample}.hisat2.ori.single.fastq \
        -N $sample_outdir/${sample}.hisat2.sortn.bam
    fq1=$sample_outdir/${sample}.hisat2.ori.1.fastq
    fq2=$sample_outdir/${sample}.hisat2.ori.2.fastq

    fq1_clean=$sample_outdir/${sample}_1.clean.fastq.gz
    fq2_clean=$sample_outdir/${sample}_2.clean.fastq.gz
    fq1_unpair=$sample_outdir/${sample}_1.unpair.fastq.gz
    fq2_unpair=$sample_outdir/${sample}_2.unpair.fastq.gz
    

    trimmomatic PE -threads $threads \
    $fq1 $fq2 \
    $fq1_clean $fq1_unpair \
    $fq2_clean $fq2_unpair \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

    conda deactivate
    hisat2 -x $hisat_ref \
    -1 $fq1_clean \
    -2 $fq2_clean \
    -S $sample_outdir/hisat2.sam \
    -p $threads

    samtools sort -@ $threads -O BAM -o $sample_outdir/hisat2.bam $sample_outdir/hisat2.sam
    samtools index -@ $threads $sample_outdir/hisat2.bam



    in_bam=/scratch/project/cs_shuaicli/wxd/virus/ShangHai_GanDan_ICC/RNA/${sample}/${sample}.bam
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

    java -jar /home/xuedowang2/app/picard.jar \
    SamToFastq \
    -I $sample_outdir/hisat.unmap1.stillunmap.bam \
    -F $sample_outdir/hisat.unmap1.stillunmap.fq

    java -jar /home/xuedowang2/app/picard.jar \
    SamToFastq \
    -I $sample_outdir/hisat.unmap2.stillunmap.bam \
    -F $sample_outdir/hisat.unmap2.stillunmap.fq

    java -jar /home/xuedowang2/app/picard.jar \
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