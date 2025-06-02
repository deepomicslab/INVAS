#!/bin/sh
#SBATCH --cpus-per-task=64
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=07-00:00:00
#SBATCH --job-name=single_sample_align
#SBATCH --output=log.align.%j.log

# 使用方法: ./script.sh <sample_name> <fq1_path> <fq2_path> <bwa_ref> <hisat_ref> <outdir>
# 示例: ./script.sh SAMPLE001 /path/to/sample.R1.fastq.gz /path/to/sample.R2.fastq.gz /path/to/ref.fa /path/to/hisat_ref /path/to/outdir

# 检查参数数量
if [ $# -ne 6 ]; then
    echo "用法: $0 <sample_name> <fq1_path> <fq2_path> <bwa_ref> <hisat_ref> <outdir>"
    exit 1
fi

# 获取参数
sample=$1
fq1=$2
fq2=$3
ref=$4
hisat_ref=$5
outdir=$6
threads=$7


echo "Processing $sample ..."

# 创建样本输出目录
sample_outdir=$outdir/$sample
if [ ! -d "$sample_outdir" ]; then
    mkdir -p "$sample_outdir"
fi

# 定义输出文件路径
fq1_clean=$sample_outdir/${sample}_1.clean.fastq.gz
fq2_clean=$sample_outdir/${sample}_2.clean.fastq.gz
fq1_unpair=$sample_outdir/${sample}_1.unpair.fastq.gz
fq2_unpair=$sample_outdir/${sample}_2.unpair.fastq.gz

# Trimmomatic 质量控制
echo "Running Trimmomatic..."
trimmomatic PE -threads $threads \
$fq1 $fq2 \
$fq1_clean $fq1_unpair \
$fq2_clean $fq2_unpair \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

# Hisat2 比对
echo "Running Hisat2 alignment..."
hisat2 -x $hisat_ref \
-1 $fq1_clean \
-2 $fq2_clean \
-S $sample_outdir/hisat2.sam \
-p $threads

# 转换并索引 BAM 文件
samtools sort -@ $threads -O BAM -o $sample_outdir/hisat2.bam $sample_outdir/hisat2.sam
samtools index -@ $threads $sample_outdir/hisat2.bam

# 分离比对和未比对的读段
in_bam=$sample_outdir/hisat2.bam
unmap=$sample_outdir/hisat2.unmap.bam
map=$sample_outdir/hisat2.map.bam
unmap_sortn=$sample_outdir/hisat2.unmap.sortn.bam

samtools view -f4 $in_bam -O BAM -o $unmap
samtools view -F4 $in_bam -O BAM -o $map
samtools sort -n -@ $threads -O BAM -o $unmap_sortn $unmap

# 将未比对的 BAM 转换为 FASTQ
unmap_fq_up=$sample_outdir/hisat_unmap_unpaired.fastq
unmap_fq_1=$sample_outdir/hisat_unmap_1.fastq
unmap_fq_2=$sample_outdir/hisat_unmap_2.fastq
unmap_fq_s=$sample_outdir/hisat_unmap_single.fastq
samtools fastq -0 $unmap_fq_up -1 $unmap_fq_1 -2 $unmap_fq_2 $unmap_sortn -s $unmap_fq_s -N -@ $threads

# 对未比对的读段进行二次比对
echo "Re-aligning unmapped reads..."
hisat2 -x $hisat_ref -U $unmap_fq_1 -S $sample_outdir/hisat.unmap1.remap.sam -p $threads
hisat2 -x $hisat_ref -U $unmap_fq_2 -S $sample_outdir/hisat.unmap2.remap.sam -p $threads
hisat2 -x $hisat_ref -U $unmap_fq_s -S $sample_outdir/hisat.unmapup.remap.sam -p $threads

# 合并二次比对结果
samtools merge -f \
$sample_outdir/hisat.unmap_remap_all.bam \
$sample_outdir/hisat.unmap1.remap.sam \
$sample_outdir/hisat.unmap2.remap.sam \
$sample_outdir/hisat.unmapup.remap.sam

# 处理二次比对结果
samtools sort -@ $threads -O BAM -o $sample_outdir/hisat.unmap_remap_all.s.bam $sample_outdir/hisat.unmap_remap_all.bam
samtools view -F 4 -O BAM -@ $threads -o $sample_outdir/hisat.unmap_remap_all.s.map.bam $sample_outdir/hisat.unmap_remap_all.s.bam
samtools merge -f $sample_outdir/hisat.map_unmapremap.bam $map $sample_outdir/hisat.unmap_remap_all.s.map.bam
samtools sort -@ $threads -O BAM -o $sample_outdir/hisat.map_unmapremap.s.bam $sample_outdir/hisat.map_unmapremap.bam

# 提取仍然未比对的读段
samtools view -f4 $sample_outdir/hisat.unmap1.remap.sam -O BAM -o $sample_outdir/hisat.unmap1.stillunmap.bam
samtools view -f4 $sample_outdir/hisat.unmap2.remap.sam -O BAM -o $sample_outdir/hisat.unmap2.stillunmap.bam
samtools view -f4 $sample_outdir/hisat.unmapup.remap.sam -O BAM -o $sample_outdir/hisat.unmapup.stillunmap.bam

# 将仍然未比对的 BAM 转换为 FASTQ
echo "Converting still unmapped BAM to FASTQ..."
java -jar /home/wangxuedong/app/picard.jar \
SamToFastq \
-I $sample_outdir/hisat.unmap1.stillunmap.bam \
-F $sample_outdir/hisat.unmap1.stillunmap.fq

java -jar /home/wangxuedong/app/picard.jar \
SamToFastq \
-I $sample_outdir/hisat.unmap2.stillunmap.bam \
-F $sample_outdir/hisat.unmap2.stillunmap.fq

java -jar /home/wangxuedong/app/picard.jar \
SamToFastq \
-I $sample_outdir/hisat.unmapup.stillunmap.bam \
-F $sample_outdir/hisat.unmapup.stillunmap.fq

# 处理仍然未比对的读段，准备 BWA 比对
awk '{if(NR%4==1) {print $0 "_r1"} else {print $0}}' $sample_outdir/hisat.unmap1.stillunmap.fq > $sample_outdir/hisat.unmap1.stillunmap.r1.fq
awk '{if(NR%4==1) {print $0 "_r2"} else {print $0}}' $sample_outdir/hisat.unmap2.stillunmap.fq > $sample_outdir/hisat.unmap2.stillunmap.r2.fq
awk '{if(NR%4==1) {print $0 "_sg"} else {print $0}}' $sample_outdir/hisat.unmapup.stillunmap.fq > $sample_outdir/hisat.unmapup.stillunmap.sg.fq

# 合并所有仍然未比对的读段
cat $sample_outdir/hisat.unmap1.stillunmap.r1.fq $sample_outdir/hisat.unmap2.stillunmap.r2.fq $sample_outdir/hisat.unmapup.stillunmap.sg.fq > $sample_outdir/still_unmap_r1r2sg.fq

# 使用 BWA 进行最终比对
echo "Final alignment with BWA..."
bwa mem $ref $sample_outdir/still_unmap_r1r2sg.fq -t $threads > $sample_outdir/still_unmap_bwa.sam
samtools sort $sample_outdir/still_unmap_bwa.sam -O BAM -o $sample_outdir/still_unmap_bwa.s.bam
samtools index $sample_outdir/still_unmap_bwa.s.bam

echo "Processing of $sample completed successfully."