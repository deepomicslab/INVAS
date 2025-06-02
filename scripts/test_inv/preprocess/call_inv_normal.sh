sample=HG002
ngs_bam=/home/xuedowang2/scratch/HG002/bam/HG002.GRCh38.2x250.bam
tgs_bam=/home/xuedowang2/scratch/HG002/bam/PacBio_CCS_15kb_20kb_chemistry2/HG002.SequelII.merged_15kb_20kb.GRCh38.duplomap.bam
rna1_fq1=/scratch/project/cs_shuaicli/wxd/hg002_rna/mrna/fq/hg002_gm24385.mrna.R1.clean.fastq.gz
rna1_fq2=/scratch/project/cs_shuaicli/wxd/hg002_rna/mrna/fq/hg002_gm24385.mrna.R2.clean.fastq.gz
rna2_fq1=/scratch/project/cs_shuaicli/wxd/hg002_rna/mrna/fq/hg002_gm26105.mrna.R1.clean.fastq.gz
rna2_fq2=/scratch/project/cs_shuaicli/wxd/hg002_rna/mrna/fq/hg002_gm26105.mrna.R2.clean.fastq.gz
rna3_fq1=/scratch/project/cs_shuaicli/wxd/hg002_rna/mrna/fq/hg002_gm27730.mrna.R1.clean.fastq.gz
rna3_fq2=/scratch/project/cs_shuaicli/wxd/hg002_rna/mrna/fq/hg002_gm27730.mrna.R2.clean.fastq.gz
outdir=/scratch/project/cs_shuaicli/wxd/hg002_rna/mrna/sv

threads=128
ref=/home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa

#########################################call sv wgs ngs #############################################
svaba run -a $sample.ngs \
-G $ref \
-t $ngs_bam \
--override-reference-check \
--read-tracking --germline \
-p $threads \
-L 6 \
-I

python2 svaba_converter.py $sample.ngs.sv.vcf > $sample.ngs.sv.info.vcf

#########################################call sv wgs tgs #############################################
cute=~/app/spechla_run/SpecHLA/spechla_env/bin/cuteSV
#activate conda environment (optional)
eval "$(conda shell.bash hook)"
conda activate /home/xuedowang2/app/spechla_run/SpecHLA/spechla_env
out_sv_vcf=$outdir/$sample.tgs.sv.vcf
work_dir=$outdir/cute_work
if [ ! -d $work_dir ]; then
    mkdir $work_dir
fi
# dip_caller=/home/xuedowang2/app/cuteSV/src/cuteSV/diploid_calling.py

echo """
$cute $tgs_bam $ref $out_sv_vcf $work_dir \
    --max_cluster_bias_INS 1000 \
    --diff_ratio_merging_INS 0.9 \
    --max_cluster_bias_DEL 1000 \
    --diff_ratio_merging_DEL 0.5  \
    -q 10 \
    --threads $threads \
    --sample $sample \
    --genotype \
    -mi 500 -md 500 \
    --report_readid

"""


$cute $bam $ref $out_sv_vcf $work_dir \
    --max_cluster_bias_INS 1000 \
    --diff_ratio_merging_INS 0.9 \
    --max_cluster_bias_DEL 1000 \
    --diff_ratio_merging_DEL 0.5  \
    -q 10 \
    --threads $threads \
    --sample $sample \
    --genotype \
    -mi 500 -md 500 \
    --report_readid



#########################################call sv rna ################################################


bwa mem $ref -t $threads $rna1_fq1 $rna1_fq2 | samtools view -bS - > $sample.rna1.bam
bwa mem $ref -t $threads $rna2_fq1 $rna2_fq2 | samtools view -bS - > $sample.rna2.bam
bwa mem $ref -t $threads $rna3_fq1 $rna3_fq2 | samtools view -bS - > $sample.rna3.bam
samtools merge -@ $threads -O BAM -o $sample.rna.bam $sample.rna1.bam $sample.rna2.bam $sample.rna3.bam
samtools sort -@ $threads -O BAM -o $sample.rna.s.bam $sample.rna.bam
samtools index -@ $threads $sample.rna.s.bam


# call sv for rna bam
svaba run -a $sample.rna \
-G $ref \
-t $sample.rna.s.bam \
--read-tracking --germline \
-p $threads \
-L 6 \
-I

python2 svaba_converter.py $sample.rna.sv.vcf > $sample.rna.sv.info.vcf