sample=tcga_test
ngs_normal_bam=/scratch/project/cs_shuaicli/wxd/TCGA/9186c2de-aec2-4eac-ab45-3bb8ac3eb57f/46c5d43b-4e08-4fa2-94a2-a73557102507_wgs_gdc_realn.bam
ngs_tumor_bam=/scratch/project/cs_shuaicli/wxd/TCGA/3fa065e2-3c45-4f93-9f9b-6da20c57939e/16156bce-c3e7-4254-8c05-4af0bdf6d6c0_wgs_gdc_realn.bam
# tgs_bam=/home/xuedowang2/scratch/HG002/bam/PacBio_CCS_15kb_20kb_chemistry2/HG002.SequelII.merged_15kb_20kb.GRCh38.duplomap.bam
rna1_fq1=/home/xuedowang2/scratch/TCGA/rnaseq_hg38ref/9de9a8b1-b42b-4113-bab5-65f00f2c4b5c/output_1.fastq
rna1_fq2=/home/xuedowang2/scratch/TCGA/rnaseq_hg38ref/9de9a8b1-b42b-4113-bab5-65f00f2c4b5c/output_2.fastq


threads=128
ref=/home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa

#########################################call sv wgs ngs #############################################.
echo "svaba wgs calling"
# svaba run -a $sample.ngs \
# -G $ref \
# -t $ngs_bam \
# --override-reference-check \
# --read-tracking --germline \
# -p $threads \
# -L 6 \
# -I

svaba run \
-t $ngs_tumor_bam \
-n $ngs_normal_bam \
-G $ref \
-a $sample.ngs \
-p $threads \
--override-reference-check



python2 svaba_converter.py $sample.ngs.svaba.sv.vcf > $sample.ngs.svaba.sv.info.vcf




#########################################call sv rna ################################################


bwa mem $ref -t $threads $rna1_fq1 $rna1_fq2 | samtools view -bS - > $sample.rna.bam
# bwa mem $ref -t $threads $rna2_fq1 $rna2_fq2 | samtools view -bS - > $sample.rna2.bam
# bwa mem $ref -t $threads $rna3_fq1 $rna3_fq2 | samtools view -bS - > $sample.rna3.bam
# samtools merge -@ $threads -O BAM -o $sample.rna.bam $sample.rna1.bam $sample.rna2.bam $sample.rna3.bam
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

python2 svaba_converter.py $sample.rna.svaba.sv.vcf > $sample.rna.svaba.sv.info.vcf