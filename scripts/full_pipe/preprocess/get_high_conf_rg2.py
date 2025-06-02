import argparse
import pysam

# 设置命令行参数
parser = argparse.ArgumentParser(description='Find continuous regions with average depth above a certain threshold, considering gaps.')
parser.add_argument('bam', type=str, help='Path to the BAM file')
parser.add_argument('chromosome', type=str, help='Chromosome to analyze')
parser.add_argument('threshold', type=float, help='Depth threshold')
parser.add_argument('--gap', type=int, default=1000, help='Maximum gap length of low depth allowed within a region')
args = parser.parse_args()

# 从命令行参数中获取输入数据
bam_path = args.bam
chromosome = args.chromosome
depth_threshold = args.threshold
max_gap_length = args.gap

# 打开BAM文件
bam_file = pysam.AlignmentFile(bam_path, "rb")

# 初始化变量
start = None
end = None
current_region_depths = []
gap_counter = 0
last_pos = -1

# 遍历指定染色体的pileup
for pileupcolumn in bam_file.pileup(chromosome):
    pos = pileupcolumn.pos
    depth = pileupcolumn.nsegments

    # 检查是否跳过了深度为0的位置
    if last_pos != -1 and pos != last_pos + 1:
        skipped_length = pos - (last_pos + 1)
        if start is not None:
            if skipped_length <= max_gap_length:
                # 填充跳过的部分，即使深度为0
                current_region_depths.extend([0] * skipped_length)
                end = pos - 1
                gap_counter += skipped_length
            else:
                # 超出允许的gap长度，结束当前区域
                average_depth = sum(current_region_depths) / len(current_region_depths)
                if average_depth >= depth_threshold:
                    print(f"{chromosome}\t{start}\t{end}\t{average_depth:.2f}")
                start = None
                current_region_depths = []
                gap_counter = 0

    # 更新last_pos
    last_pos = pos

    # 根据深度更新区域
    if start is None and depth >= depth_threshold:
        start = pos  # 开始新的区域

    if start is not None:
        if depth >= depth_threshold:
            # 当前位置的深度高于阈值，添加到当前区域
            current_region_depths.append(depth)
            end = pos
            gap_counter = 0  # 重置gap计数器
        else:
            # 当前位置的深度低于阈值，增加gap计数器
            gap_counter += 1
            if gap_counter <= max_gap_length:
                # 如果gap还在允许的长度内，继续添加到当前区域
                current_region_depths.append(depth)
                end = pos
            else:
                # 如果gap超出了允许的长度，检查并可能报告当前区域
                average_depth = sum(current_region_depths) / len(current_region_depths)
                if average_depth >= depth_threshold:
                    print(f"{chromosome}\t{start}\t{end}\t{average_depth:.2f}")
                # 开始新的区域
                start = None
                current_region_depths = []
                gap_counter = 0

# 检查最后的区域
if start is not None:
    average_depth = sum(current_region_depths) / len(current_region_depths)
    if average_depth >= depth_threshold:
        print(f"{chromosome}\t{start}\t{end}\t{average_depth:.2f}")

# 关闭BAM文件
bam_file.close()