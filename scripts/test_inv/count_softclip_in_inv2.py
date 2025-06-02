import argparse
import pysam
import re


# 检测read中soft clipping的长度是否超过阈值
def has_soft_clipping(cigar_tuples, threshold):
    return any(op == 4 and length > threshold for op, length in cigar_tuples)

def parse_cigar2(cigar_string):
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    matches = pattern.findall(cigar_string)

    left_clip_len = 0
    right_clip_len = 0
    left_clip_type = None
    right_clip_type = None

    if matches:
        # Check the first element
        if matches[0][1] in ['S', 'H']:
            left_clip_len = int(matches[0][0])
            left_clip_type = matches[0][1]

        # Check the last element
        if matches[-1][1] in ['S', 'H']:
            right_clip_len = int(matches[-1][0])
            right_clip_type = matches[-1][1]

    return left_clip_len, left_clip_type, right_clip_len, right_clip_type

# 检查read的blocks是否与给定区域重叠
def is_read_in_region(read_blocks, region_start, region_end):
    for block_start, block_end in read_blocks:
        if block_start >= region_start and block_end <= region_end:
            return True
    return False

# 统计整个BAM文件中，与给定区域重叠且soft clipping长度超过阈值的reads数量
def count_soft_clipped_reads(bam_file, regions, threshold):
    soft_clipped_counts = {region: 0 for region in regions}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if not read.is_unmapped and 'S' in read.cigarstring:
                left_soft_clip, left_tag, right_soft_clip, right_tag = parse_cigar2(read.cigarstring)
                if left_soft_clip <4 and right_soft_clip <4:
                    continue
                read_blocks = read.get_blocks()
                for region in regions:
                    chrom, start_end = region.split(':')
                    start, end = map(int, start_end.split('-'))
                    if chrom == bam.get_reference_name(read.reference_id):
                        if is_read_in_region(read_blocks, start, end):
                            soft_clipped_counts[region] += 1
    return soft_clipped_counts

# 主函数，读取区域文件并统计soft clipping reads数
def main(region_file, bam_file, threshold):
    with open(region_file, 'r') as file:
        regions = [line.split()[0] for line in file]
    counts = count_soft_clipped_reads(bam_file, regions, threshold)
    for region, count in counts.items():
        print(f"Region: {region}, Soft clipped reads: {count}")

# 接受命令行参数
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Count soft-clipped reads that overlap with specified regions in a BAM file and exceed a length threshold.')
    parser.add_argument('region_file', type=str, help='File containing regions (format: chr:start-end)')
    parser.add_argument('bam_file', type=str, help='Path to the BAM file')
    parser.add_argument('threshold', type=int, help='Length threshold for soft clipping')

    args = parser.parse_args()
    main(args.region_file, args.bam_file, args.threshold)