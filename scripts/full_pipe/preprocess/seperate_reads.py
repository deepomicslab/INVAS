import pysam
import sys

# 确保正确的命令行参数
if len(sys.argv) != 5:
    print("Usage: python3 script.py <bam_file> <region> <padding> <out_path>")
    sys.exit(1)

# 输入参数
bam_file = sys.argv[1]
region = sys.argv[2]
padding = int(sys.argv[3])
out_path = sys.argv[4]

# 解析区域和弹性范围来计算新的区域
chr, pos = region.split(':')
start, end = map(int, pos.split('-'))

# 计算带有弹性范围的新起止点
padded_start = max(0, start - padding)
padded_end = end + padding

# 输出文件
included_reads_file = f"{out_path}/included_reads_{chr}_{start}_{end}.bam"
excluded_reads_file = f"{out_path}/excluded_reads_{chr}_{start}_{end}.bam"

# 打开 BAM 文件
bam = pysam.AlignmentFile(bam_file, "rb")

# 创建输出文件
included_reads = pysam.AlignmentFile(included_reads_file, "wb", template=bam)
excluded_reads = pysam.AlignmentFile(excluded_reads_file, "wb", template=bam)

# 提取 reads
for read in bam.fetch(chr, padded_start, padded_end):
    if read.reference_start >= padded_start and read.reference_end <= padded_end:
        # Reads 完全在给定区域内
        included_reads.write(read)
    else:
        # Reads 在弹性范围内，但不在原始区域内
        excluded_reads.write(read)

# 关闭文件
bam.close()
included_reads.close()
excluded_reads.close()

print(f"Included reads saved to: {included_reads_file}")
print(f"Excluded reads saved to: {excluded_reads_file}")