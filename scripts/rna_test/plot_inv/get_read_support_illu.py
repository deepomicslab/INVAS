import pandas as pd
import pysam
import argparse
import random

def read_depth_file(depth_file):
    """读取深度文件"""
    df = pd.read_csv(depth_file, sep='\t', header=None, names=['chr', 'pos', 'depth'])
    return df

def read_vcf_file(vcf_file):
    """使用pysam读取VCF文件并提取INV区域"""
    # get all inv reads name in RNAMES tag
    inv_regions = []
    vcf = pysam.VariantFile(vcf_file)
    for record in vcf:
        if record.info['SVTYPE'] == 'INV':
            chrom = record.chrom
            start = record.start
            end = record.stop
            inv_regions.append((chrom, start, end))

    return inv_regions

def get_depth_for_inversions(depth_df, inv_regions):
    """获取倒位区域的深度信息"""
    results = []
    for chr, start, end in inv_regions:
        mask0 = (depth_df['chr'] == chr) & (depth_df['pos'] <= start)
        mask1 = (depth_df['chr'] == chr) & (depth_df['pos'] >= start) & (depth_df['pos'] <= end)
        mask2 = (depth_df['chr'] == chr) & (depth_df['pos'] >= end)
        subset0 = depth_df.loc[mask0]
        subset1 = depth_df.loc[mask1]
        # add two columns, one is dp1, one is dp2, default value is 0
        subset1['dp1'] = 0
        subset1['dp2'] = 0
        subset2= depth_df.loc[mask2]
    return subset0, subset1, subset2

def get_clipping_info_of_read(read):
    """
    获取单个read中clipping的详细位置信息。
    
    参数:
        read (pysam.AlignedSegment): 一个pysam AlignedSegment对象，表示一个单个的read。
    
    返回:
        list: 包含read的clipping详细信息的列表，每个元素是一个字典。
    """
    read_clips = []
    cigartuples = read.cigartuples
    query_position = 0  # Track the position within the read

    if cigartuples:
        # 检查每个CIGAR操作
        for operation, length in cigartuples:
            if operation == 4:  # Soft clip
                clip_type = 'soft'
                clip_start = query_position
                clip_end = query_position + length - 1
                read_clips.append({
                    'type': clip_type,
                    'start': clip_start,
                    'end': clip_end,
                    'length': length
                })
                query_position += length  # Update position for soft clips

            elif operation == 5:  # Hard clip
                clip_type = 'hard'
                clip_start = query_position
                clip_end = query_position + length - 1
                read_clips.append({
                    'type': clip_type,
                    'start': clip_start,
                    'end': clip_end,
                    'length': length
                })
                # Hard clips do not contribute to query position

            elif operation in [0, 1, 2, 7, 8]:  # M, I, D, =, X operations
                query_position += length  # Update position for these operations

    return read_clips

def do_intervals_overlap(interval1, interval2):
    """
    Check if two intervals overlap.
    
    Each interval is defined by a tuple (start, end).
    
    Args:
    interval1, interval2: Tuples defining the intervals.
    
    Returns:
    True if the intervals overlap, False otherwise.
    """

    # Unpack interval coordinates
    a, b = interval1
    c, d = interval2
    
    # Check if intervals overlap
    return max(a, c) < min(b, d)

def get_all_reads(bam, inv_regions, ext=10):
    support_reads = set()
    chrom, start, end = inv_regions[0]
    for read in bam.fetch(chrom, start, end):
        # res = get_clipping_info_of_read(read)
        # get clip break point
        res = read.get_blocks()
        # print(res)
        # if break point in inv region, add read name to support reads
        is_support = False
        for clip_start, clip_end in res:
            if clip_start >= start-ext and clip_end <= end+ext:
                is_support = True
                break
        if is_support:
            support_reads.add(read.query_name)

    
    return support_reads

def assign_reads_in_invs(subset, inv_regions):
    # get all reads from bam file for inv regions
    bam = pysam.AlignmentFile(args.bam_file, "rb")
    all_supp_reads_names = get_all_reads(bam, inv_regions)
    print("support reads number:", len(all_supp_reads_names))
    for chr, start, end in inv_regions:
        # get allreads for each reads
        for i in range(start, end+1):
            # loop all aligned reads in this position
            for read in bam.fetch(chr, i-1, i):
                # if read is in inv region, assign reads to dp1 or dp2
                if read.query_name in all_supp_reads_names:
                    subset.loc[subset['pos'] == i, 'dp1'] += 1
                else:
                    subset.loc[subset['pos'] == i, 'dp2'] += 1  

def main(depth_file, vcf_file):
    depth_df = read_depth_file(depth_file)
    print(depth_df.head())
    inv_regions = read_vcf_file(vcf_file)
    print("inv:",inv_regions)
    if len(inv_regions) == 0:
        print("no inv regions")
        res_dp1=[]
        for idx, item in enumerate(depth_df['depth']):
            res_dp1.extend([str(depth_df['pos'].iloc[idx])] * item)
        out_f=open(args.out, 'w')
        out_f.write("\n")
        out_f.write(",".join(res_dp1) + "\n")
        out_f.close()
        ratio_out_f=open(args.ratio_out, 'w')
        for idx, item in enumerate(depth_df['pos']):
            num = random.uniform(0.99, 1)
            ratio_out_f.write(str(item) + "," + str(num) + "\n")
        ratio_out_f.close()
    else:
        subset0, subset1, subset2 = get_depth_for_inversions(depth_df, inv_regions)
        assign_reads_in_invs(subset1, inv_regions)
        print(subset0)
        print(subset1)
        print(subset2)
        # wirte subset to file, the dirst row is dp1, the second row is dp2
        res_dp1=[]
        res_dp2=[]
        for idx, item in enumerate(subset0['depth']):
            res_dp2.extend([str(subset0['pos'].iloc[idx])] * item)
        out_f=open(args.out, 'w')
        
        for idx, item in enumerate(subset1['dp1']):
            res_dp1.extend([str(subset1['pos'].iloc[idx])] * item)
            res_dp2.extend([str(subset1['pos'].iloc[idx])] * subset1['dp2'].iloc[idx])
        for idx, item in enumerate(subset2['depth']):
            res_dp2.extend([str(subset2['pos'].iloc[idx])] * item)
        res1 = list(set(res_dp1))
        res1.sort()
        res2 = list(set(res_dp2))
        res2.sort()
        print(res1)
        # print(res2)

        out_f.write(",".join(res_dp1) + "\n")
        out_f.write(",".join(res_dp2) + "\n")
        out_f.close()

        ratio_out_f=open(args.ratio_out, 'w')
        for idx, item in enumerate(subset0['pos']):
            num = random.uniform(0.98, 1)
            ratio_out_f.write(str(item) + "," + str(num) + "\n")
        for idx, item in enumerate(subset1['pos']):
            num = random.uniform(0, 0.02)

            # print(subset1['pos'].iloc[idx], subset1['dp2'].iloc[idx], float(subset1['dp2'].iloc[idx]))
            # print(subset1['dp1'].iloc[idx], float(subset1['dp1'].iloc[idx]))
            num=float(subset1['dp2'].iloc[idx])/(float(subset1['dp1'].iloc[idx])+float(subset1['dp2'].iloc[idx]))+num
            ratio_out_f.write(str(item) + "," + str(num) + "\n")
        for idx, item in enumerate(subset2['pos']):
            num = random.uniform(0.97,1)
            ratio_out_f.write(str(item) + "," + str(num) + "\n")
        ratio_out_f.close()

        with open(args.for_circle_dp, 'w') as f:
            f.write("pos\tdp1\tdp2\n")
            for idx, item in enumerate(subset0['pos']):
                f.write(str(item) + "\t0\t" + str(subset0['depth'].iloc[idx]) + "\n")
            for idx, item in enumerate(subset1['pos']):
                f.write(str(item) + "\t" + str(subset1['dp1'].iloc[idx]) + "\t" + str(subset1['dp2'].iloc[idx]) + "\n")
            for idx, item in enumerate(subset2['pos']):
                f.write(str(item) + "\t0\t" + str(subset2['depth'].iloc[idx]) + "\n")



    

    # print(subset.head())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process depth and VCF files to find depth of inversions.")
    parser.add_argument("depth_file", help="Path to the depth file")
    parser.add_argument("vcf_file", help="Path to the VCF file")
    parser.add_argument("bam_file", help="Path to the VCF file")
    parser.add_argument("out", help="Path to the csv file")
    parser.add_argument("ratio_out", help="Path to the csv file")
    parser.add_argument("for_circle_dp", help="Path to the csv file")
    args = parser.parse_args()

    main(args.depth_file, args.vcf_file)