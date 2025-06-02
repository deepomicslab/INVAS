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
    right_soft_clip_cnt = 0
    left_soft_clip_cnt = 0
    for read in bam_file.fetch(regions[0], regions[1], regions[2]):
        if not read.is_unmapped and ('S' in read.cigarstring or 'H' in read.cigarstring):
            left_soft_clip, left_tag, right_soft_clip, right_tag = parse_cigar2(read.cigarstring)
            if left_soft_clip >threshold:
                left_soft_clip_cnt += 1
            if right_soft_clip >threshold:
                right_soft_clip_cnt += 1
            #     and right_soft_clip <threshold:
            #     continue
            # right_soft_clip += 1
            # left_soft_clip += 1
    return left_soft_clip_cnt, right_soft_clip_cnt
        # if left_soft_clip <4 and right_soft_clip <4:
        #     continue
    #     read_blocks = read.get_blocks()
    #     for region in regions:
    #         chrom, start_end = region.split(':')
    #         start, end = map(int, start_end.split('-'))
    #         if chrom == bam.get_reference_name(read.reference_id):
    #             if is_read_in_region(read_blocks, start, end):
    #                 soft_clipped_counts[region] += 1
    # return soft_clipped_counts


        # for read in bam:
        #     if not read.is_unmapped and 'S' in read.cigarstring:
        #         left_soft_clip, left_tag, right_soft_clip, right_tag = parse_cigar2(read.cigarstring)
        #         if left_soft_clip <4 and right_soft_clip <4:
        #             continue
        #         read_blocks = read.get_blocks()
        #         for region in regions:
        #             chrom, start_end = region.split(':')
        #             start, end = map(int, start_end.split('-'))
        #             if chrom == bam.get_reference_name(read.reference_id):
        #                 if is_read_in_region(read_blocks, start, end):
        #                     soft_clipped_counts[region] += 1
    # return soft_clipped_counts

def parse_regions():
    # file format:
    # chr1:9823293-9824715    CLSTN1:9789083-9884584  chr1:9823221-9823310-7.00
    # chr1:9823293-9824715    CLSTN1:9789083-9884584  chr1:9824647-9824721-11.35
    inv_dict={}
    exon_inv_dict={}
    
    with open(args.region_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            inv_region = parts[0]
            gene_info =parts[1]
            exon_region = parts[2]
            if exon_region in exon_inv_dict:
                print(f"Error: {exon_region} already in exon_inv_dict")
            exon_inv_dict[exon_region]=(parts, gene_info)   
            if inv_region in inv_dict:
                inv_dict[inv_region].append(exon_region)
            else:
                inv_dict[inv_region]=[(parts,gene_info)]
    return inv_dict,exon_inv_dict

# 主函数，读取区域文件并统计soft clipping reads数
def main(region_file, bam_file, threshold):
    inv_dict,exon_inv_dict = parse_regions()
    bam_file=pysam.AlignmentFile(bam_file, "rb")


    # for small region check double side of small region
    for small_inv_str, (parts, gene) in exon_inv_dict.items():
        item=small_inv_str.split(':')
        chrom=item[0]
        start, end = map(int, item[1].split('-')[:2])
        # chr10:69414319-69414389
        # if chrom=="chr10" and start==69414319 and end==69414389:
        left_count, right_count = count_soft_clipped_reads(bam_file, (chrom,start, end), threshold)
        if left_count>0 and right_count>0:
            print(f"processing {chrom}:{start}-{end}\tinv region:{parts[0]}\tgene_info:{gene}")
            print(f"Left: {left_count}, Right: {right_count}")
    # for large region check double side of gene
    for large_inv_str, parts in inv_dict.items():
        pass





    # counts = count_soft_clipped_reads(bam_file, regions, threshold)
    # for region, count in counts.items():
    #     print(f"Region: {region}, Soft clipped reads: {count}")

# 接受命令行参数
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Count soft-clipped reads that overlap with specified regions in a BAM file and exceed a length threshold.')
    parser.add_argument('region_file', type=str, help='File containing regions (format: chr:start-end)')
    parser.add_argument('bam_file', type=str, help='Path to the BAM file')
    parser.add_argument('threshold', type=int, help='Length threshold for soft clipping')

    args = parser.parse_args()
    main(args.region_file, args.bam_file, args.threshold)

    # for small region check double side of small region
    # for large region check double side of gene