import pysam
import sys

# 获取命令行参数
bam_file_path = sys.argv[1]
distance_threshold = int(sys.argv[2])  # 距离阈值用于定义read clusters
output_file_path = sys.argv[3]
min_reads_in_cluster = 10  # cluster中最少的reads数目

# 定义处理cluster reads的函数
def process_cluster_reads(cluster):
    left_soft_clip_count = 0
    right_soft_clip_count = 0
    start = cluster['start']
    end = cluster['end']
    
    # 检查覆盖cluster起始位置的reads，统计left soft clip的数量
    for read in cluster['reads']:
        if read.reference_start <= start < read.reference_end:
            # 检查CIGAR是否以soft clip开始
            if read.cigartuples[0][0] == 4:
                left_soft_clip_count += 1
    
    # 检查覆盖cluster结束位置的reads，统计right soft clip的数量
    for read in cluster['reads']:
        if read.reference_start < end <= read.reference_end:
            # 检查CIGAR是否以soft clip结束
            if read.cigartuples[-1][0] == 4:
                right_soft_clip_count += 1
    
    # 返回左侧和右侧soft-clipped的数量
    return left_soft_clip_count, right_soft_clip_count

# 打开BAM文件
bam_file = pysam.AlignmentFile(bam_file_path, "rb")

# 用于存储read clusters的列表
clusters = []
prev_chrom = ""
# 遍历BAM文件中的每个read
for read in bam_file.fetch():
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        continue
    
    # 获取read的染色体名称
    chrom = bam_file.get_reference_name(read.reference_id)
    if chrom != prev_chrom:
        print(f"Processing {chrom}")
        prev_chrom = chrom

    # 将每个read尝试放入现有的clusters
    placed_in_cluster = False
    for cluster in clusters:
        # 检查染色体是否匹配，并且read的起始位置是否在cluster的结束位置的threshold内
        if read.reference_start <= cluster['end'] + distance_threshold and chrom == cluster['chrom']:
            cluster['reads'].append(read)
            cluster['end'] = max(cluster['end'], read.reference_end)
            placed_in_cluster = True
            break

    # 如果read没有放入任何cluster，则创建一个新的cluster
    if not placed_in_cluster:
        clusters.append({
            'chrom': chrom,
            'start': read.reference_start,
            'end': read.reference_end,
            'reads': [read]
        })

# 对每个cluster的reads进行进一步处理
out_bed = open(output_file_path, 'w')
for cluster in clusters:
    if len(cluster['reads']) >= min_reads_in_cluster:
        left_clip_count, right_clip_count = process_cluster_reads(cluster)
        if left_clip_count > 0 and right_clip_count > 0:
            out_bed.write(f"{cluster['chrom']}\t{cluster['start']}\t{cluster['end']}\n")
            print(f"Cluster on {cluster['chrom']} from {cluster['start']} to {cluster['end']}, "
                f"contains {len(cluster['reads'])} reads, "
                f"with {left_clip_count} left soft clips, "
                f"and {right_clip_count} right soft clips.")

# 关闭BAM文件
bam_file.close()
out_bed.close()