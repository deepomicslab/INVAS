#!/bin/sh
#SBATCH --partition=batch
#SBATCH --cpus-per-task=128
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=03-00:00:00
#SBATCH --job-name=realign_sv
#SBATCH --output=log.realign_sv.log


outdir=/scratch/project/cs_shuaicli/wxd/virus/ShangHai_GanDan_ICC/WGS/svs
threads=128
ref=/home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
samples_file=/scratch/project/cs_shuaicli/wxd/virus/ShangHai_GanDan_ICC/samples.txt
delly=/home/xuedowang2/app/delly_v1.1.6_linux_x86_64bit
svaba_converter=/home/wangxuedong/svaba_converter.py
samtools=~/app/samtools-1.12/samtools

for sample in `cat $samples_file`; do
    sample_dir=$outdir/$sample
    sample_bam=/scratch/project/cs_shuaicli/wxd/virus/ShangHai_GanDan_ICC/WGS/$sample/$sample.bam

    




    if [ ! -d $sample_dir ]; then
        mkdir $sample_dir
    fi

    # extract fastq at first
    samtools sort -n -@ $threads -O BAM -o $sample_dir/$sample\_sortn.bam $sample_bam
    samtools fastq -@ $threads \
        -1 $sample_dir/$sample\_1.fastq.gz \
        -2 $sample_dir/$sample\_2.fastq.gz \
        -s $sample_dir/$sample\_s.fastq.gz \
        $sample_dir/$sample\_sortn.bam

    bwa mem -t $threads $ref \
        $sample_dir/$sample\_1.fastq.gz \
        $sample_dir/$sample\_2.fastq.gz \
        | samtools view -b -@ $threads -o $sample_dir/$sample\_bwa_hg38_realign.bam
    
    sample_bam=$sample_dir/$sample\_bwa_hg38_realign.bam
    samtools sort -@ $threads -O BAM -o $sample_dir/$sample\_bwa_hg38_realign_sort.bam $sample_bam
    sample_bam=$sample_dir/$sample\_bwa_hg38_realign_sort.bam
    samtools index -@ $threads $sample_bam

    echo "processing $sample"
    echo "run delly"
    sample_delly_dir=$sample_dir/delly
    if [ ! -d $sample_delly_dir ]; then
        mkdir $sample_delly_dir
    fi
    cd $sample_delly_dir
    delly_ovcf=$sample_delly_dir/$sample.delly.vcf
    $delly call -g $ref $sample_bam  > $delly_ovcf

    cd -
    eval "$(conda shell.bash hook)"
    conda activate manta
    echo "run svaba"
    svaba_dir=$sample_dir/svaba
    if [ ! -d $svaba_dir ]; then
        mkdir $svaba_dir
    fi
    cd $svaba_dir
    svaba run -a $sample \
    -G $ref \
    -t $sample_bam \
    --override-reference-check \
    --read-tracking --germline \
    -p $threads \
    -L 6 \
    -I 
    python2 $svaba_converter $svaba_dir/$sample.svaba.sv.vcf $svaba_dir/$sample.inv.bed > $svaba_dir/$sample.svaba.sv.info.vcf
    cd -

    echo "run mantasv"
    manta_dir=$sample_dir/manta
    if [ ! -d $manta_dir ]; then
        mkdir $manta_dir
    fi
    cd $manta_dir
    /home/xuedowang2/miniconda3/envs/py27/bin/configManta.py \
        --bam $sample_bam \
        --referenceFasta $ref \
        --runDir $manta_dir
    python $manta_dir/runWorkflow.py

    /home/xuedowang2/miniconda3/envs/py27/bin/convertInversion.py \
        $samtools \
        $ref \
        $manta_dir/results/variants/diploidSV.vcf.gz \
        > $manta_dir/$sample.manta.invfmt.vcf
    
    cd -


    
    
done