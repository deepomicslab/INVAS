    fq1=/data6/yangtianxia/RNA/${sample}.R1.fastq.gz
    fq2=/data6/yangtianxia/RNA/${sample}.R2.fastq.gz
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