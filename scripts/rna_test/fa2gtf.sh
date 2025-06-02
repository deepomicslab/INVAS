fq1=$1
fq2=$2
outdir=$3
threads=$4
hisat2_ref=/scratch/project/cs_shuaicli/hisat2_ref/GRCh38_full_analysis_set_plus_decoy_hla
if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi
hisat2 -x $hisat2_ref -1 $fq1 -2 $fq2 -p $threads -S $outdir/hisat2.sam
samtools sort -@ $threads -o $outdir/hisat2.bam $outdir/hisat2.sam
samtools index -@ $threads $outdir/hisat2.bam
stringtie $outdir/hisat2.bam -o $outdir/hisat2.stringtie.gtf
