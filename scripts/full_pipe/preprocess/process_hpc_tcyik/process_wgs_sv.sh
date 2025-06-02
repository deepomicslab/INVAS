#!/bin/sh
#SBATCH --partition=stingy
#SBATCH --cpus-per-task=128
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=00-06:00:00
#SBATCH --job-name=TC-YIK-wgs
#SBATCH --output=log.yikwgs.log


outdir=/scratch/project/cs_shuaicli/wxd/virus/WuHan_TongJi_SCCC/TC-YIK/WGS/svs
threads=128
ref=/home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
delly=/home/xuedowang2/app/delly_v1.1.6_linux_x86_64bit
svaba_converter=/home/wangxuedong/svaba_converter.py
samtools=~/app/samtools-1.12/samtools

# for sample in `cat $samples_file`; do
sample=TC-YIK
sample_dir=$outdir
# sample_bam=/scratch/project/cs_shuaicli/wxd/virus/ShangHai_GanDan_ICC/WGS/$sample/$sample.bam

# # extract fastq at first
# samtools sort -n -@ $threads -O BAM -o $sample_dir/$sample\_sortn.bam $sample_bam
# samtools fastq -@ $threads \
#     -1 $sample_dir/$sample\_1.fastq.gz \
#     -2 $sample_dir/$sample\_2.fastq.gz \
#     -s $sample_dir/$sample\_s.fastq.gz \
#     $sample_dir/$sample\_sortn.bam
fq1=/scratch/project/cs_shuaicli/wxd/virus/WuHan_TongJi_SCCC/TC-YIK/WGS/170524_XE00513_FCHJNGMALXX_L2_wHAXPI049353-12_1.fq.gz
fq2=/scratch/project/cs_shuaicli/wxd/virus/WuHan_TongJi_SCCC/TC-YIK/WGS/170524_XE00513_FCHJNGMALXX_L2_wHAXPI049353-12_2.fq.gz

bwa mem -t $threads $ref \
    $fq1 \
    $fq2 \
    | samtools view -b -@ $threads -o $sample_dir/$sample\_bwa_hg38_realign.bam

sample_bam=$sample_dir/$sample\_bwa_hg38_realign.bam
samtools sort -@ $threads -O BAM -o $sample_dir/$sample\_bwa_hg38_realign.sort.bam $sample_bam
sample_bam=$sample_dir/$sample\_bwa_hg38_realign.sort.bam
samtools index -@ $threads $sample_bam





if [ ! -d $sample_dir ]; then
    mkdir $sample_dir
fi
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
python2 $svaba_converter $svaba_dir/$sample.svaba.sv.vcf $sample.inv.bed > $svaba_dir/$sample.svaba.sv.info.vcf
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

    lumpy_out_dir=$sample_dir/lumpy
    if [ ! -d $lumpy_out_dir ]; then
        mkdir -p $lumpy_out_dir
    fi
    cd $lumpy_out_dir

    extract=/home/xuedowang2/app/lumpy-sv/scripts/extractSplitReads_BwaMem
    svtyper=~/miniconda3/envs/py27/bin/svtyper
    echo """
    samtools view -b -F 1294 $bam > $sample.discordants.unsorted.bam
    # Extract the split-read alignments
    samtools view -h $bam \
        | $extract -i stdin \
        | samtools view -Sb - \
        > $sample.splitters.unsorted.bam

    # Sort both alignments
    samtools sort $sample.discordants.unsorted.bam -o $sample.discordants.bam --write-index
    samtools sort $sample.splitters.unsorted.bam -o $sample.splitters.bam --write-index

    lumpyexpress \
        -B $bam \
        -S $sample.splitters.bam \
        -D $sample.discordants.bam \
        -o $sample.lumpy.vcf

    $svtyper \
        -B $bam \
        -S $sample.splitters.bam \
        -i $sample.lumpy.vcf \
        -o $sample.lumpy.svtyper.vcf

    """


    samtools view -b -F 1294 -@ $threads $bam > $sample.discordants.unsorted.bam
    # Extract the split-read alignments
    samtools view -h $bam \
        | $extract -i stdin \
        | samtools view -Sb -o $sample.splitters.unsorted.bam

    # Sort both alignments
    samtools sort -@ $threads $sample.discordants.unsorted.bam -o $sample.discordants.bam --write-index
    samtools sort -@ $threads $sample.splitters.unsorted.bam -o $sample.splitters.bam --write-index

    lumpyexpress \
        -B $bam \
        -S $sample.splitters.bam \
        -D $sample.discordants.bam \
        -o $sample.lumpy.vcf

    $svtyper \
        -B $bam \
        -l $bam.json \
        -i $sample.lumpy.vcf \
        -o $sample.lumpy.svtyper.vcf



    cd -



    
    
# done