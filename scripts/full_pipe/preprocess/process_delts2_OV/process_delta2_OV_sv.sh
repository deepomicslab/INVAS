outdir=/mnt/disk2_workspace/wangxuedong/OV/OV_WGS_SV
threads=32
ref=/mnt/disk2_workspace/panguangze/reference/hg38.fa
samples_file=has_rnaseq_ov.sample
delly=/home/wangxuedong/delly_v1.2.6_linux_x86_64bit
svaba_converter=/home/wangxuedong/svaba_converter.py
samtools=~/samtools-1.19/build/bin/samtools

for sample in `cat $samples_file`; do
    sample_dir=$outdir/$sample
    sample_bam=/mnt/delta_mars/CRM_OV/WGS/Bam/$sample/$sample\_sorted_dedup.bam
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
    python $svaba_converter $svaba_dir/$sample.svaba.sv.vcf $svaba_dir/$sample.inv.bed > $svaba_dir/$sample.svaba.sv.info.vcf
    cd -

    echo "run mantasv"
    manta_dir=$sample_dir/manta
    if [ ! -d $manta_dir ]; then
        mkdir $manta_dir
    fi
    cd $manta_dir
    ~/miniconda3/envs/manta/bin/configManta.py \
        --bam $sample_bam \
        --referenceFasta $ref \
        --runDir $manta_dir
    python $manta_dir/runWorkflow.py

    ~/miniconda3/envs/manta/bin/convertInversion.py \
        $samtools \
        $ref \
        $manta_dir/results/variants/diploidSV.vcf.gz \
        > $manta_dir/$sample.manta.invfmt.vcf

    cd -


    
    
done