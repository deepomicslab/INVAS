import pysam
import argparse

def parse_gtf(gtf_file):
    exons = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if parts[2] == 'exon':
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                exon = (start, end)
                exons.setdefault(chrom, []).append(exon)
    for chrom in exons:
        exons[chrom].sort(key=lambda x: x[0])
    return exons

def find_nearest_exon(exons, chrom, position):
    nearest_exon = None
    min_dist = float('inf')
    for exon in exons.get(chrom, []):
        exon_start, exon_end = exon
        if exon_end < position:
            dist = position - exon_end
        elif exon_start > position:
            dist = exon_start - position
        else:
            return exon
        if dist < min_dist:
            min_dist = dist
            nearest_exon = exon
    return nearest_exon

def create_exon_chain(exons, chrom, blocks):
    chain = []
    for block_start, block_end in blocks:
        start_exon = find_nearest_exon(exons, chrom, block_start)
        end_exon = find_nearest_exon(exons, chrom, block_end)
        # Choose the closer exon for the start of the block
        if start_exon and end_exon:
            if abs(start_exon[0] - block_start) < abs(end_exon[1] - block_end):
                chosen_exon = start_exon
            else:
                chosen_exon = end_exon
        elif start_exon:
            chosen_exon = start_exon
        elif end_exon:
            chosen_exon = end_exon
        else:
            continue
        # Append if not already in the chain
        if not chain or chain[-1] != chosen_exon:
            chain.append(chosen_exon)
    return tuple(chain)

def is_alternative_splicing(blocks):
    # Check if there are multiple blocks separated by at least the minimum intron size
    # and less than the maximum intron size.
    for i in range(len(blocks) - 1):
        intron_size = blocks[i+1][0] - blocks[i][1]
        if intron_size > min_intron_size and intron_size < max_intron_size:
            return True
    return False

def deduplicate_intervals(intervals_dict, error_margin=1):
    deduplicated_dict = {}
    for chromosome, intervals in intervals_dict.items():
        deduplicated_intervals = []
        for interval in sorted(intervals):  # Sort intervals to compare them
            if not deduplicated_intervals:
                # If list is empty, add the first interval
                deduplicated_intervals.append(interval)
            else:
                # Get the last interval in the deduplicated list
                last_interval = deduplicated_intervals[-1]
                # Check if the current interval is within the error_margin
                if (abs(interval[0] - last_interval[0]) <= error_margin and
                        abs(interval[1] - last_interval[1]) <= error_margin):
                    # If it is within the error_margin, skip it (it's a duplicate)
                    continue
                else:
                    # Otherwise, add the current interval to the list
                    deduplicated_intervals.append(interval)
        # Add the deduplicated list of intervals to the dictionary
        deduplicated_dict[chromosome] = deduplicated_intervals
    return deduplicated_dict

def parse_chain_from_bam(bam, exons):
    exon_chains = {}
    all_rec_cnt=0
    for read in bam.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        
        chrom = bam.get_reference_name(read.reference_id)
        if "read289/ENST00000399688;mate1:4-103;mate2:145-243" in read.query_name:
            print(read.cigarstring)
        blocks = read.get_blocks()
        if not is_alternative_splicing(blocks):
            continue
        all_rec_cnt+=1
        # if len(blocks) >= 2:
        #     print(blocks, read.query_name)
        exon_chain = create_exon_chain(exons, chrom, blocks)
        
        if exon_chain:
            exon_chains[exon_chain] = exon_chains.get(exon_chain, 0) + 1
    print(f"Total records: {all_rec_cnt}")
    
    return exon_chains

def parse_inv(inv_bed):
    # file format :
    # chr22:16390136-16390220	test:16387694-16390887,test	chr22:16390136-16390220-15.33
    # parse the first col, get the inv region, and store in a dict, key is chrom, value is a list of (start, end)
    # some duplicated inv regions, only keep one
    inv_rgs = {}
    with open(inv_bed, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            inv_rg_str = parts[0]
            chrom = inv_rg_str.split(':')[0]
            start = int(inv_rg_str.split(':')[1].split('-')[0])
            end = int(inv_rg_str.split(':')[1].split('-')[1])
            # use set to store the inv regions
            inv_rgs.setdefault(chrom, set()).add((start, end))

    return inv_rgs

# def calculate_region_depth(bam, chrom, start, end):
#     # 打开BAM文件
#     # 初始化总深度为0
#     total_depth = 0
#     # 初始化覆盖的碱基数量为0
#     covered_bases = 0
    
#     # 遍历指定区域的pileup列
#     for pileupcolumn in bam.pileup(chrom, start, end, nofilter=True):
#         # 确保pileup列在指定的区域内
#         if pileupcolumn.reference_pos >= start and pileupcolumn.reference_pos < end:
#             # 累加该位置的读段数量到总深度
#             total_depth += pileupcolumn.nsegments
#             # 增加覆盖的碱基数量
#             covered_bases += 1
    
#     # 如果有覆盖的碱基，则计算平均深度
#     if covered_bases > 0:
#         average_depth = total_depth / covered_bases
#     else:
#         average_depth = 0

#     return average_depth
def calculate_region_depth(bam, chrom, start, end):
    # 使用count_coverage方法获取指定区域的覆盖度信息
    coverage = bam.count_coverage(chrom, start, end, quality_threshold=0)
    total_depth = 0
    covered_bases = 0

    # coverage包含四个元素的元组，每个元素对应一个碱基（A, C, G, T）的覆盖度列表
    for base_idx in range(4):  # 遍历每个碱基
        for pos_depth in coverage[base_idx]:  # 遍历该碱基对应的覆盖度列表
            total_depth += pos_depth  # 累加到总深度
            if pos_depth > 0:
                covered_bases += 1  # 如果该位置有覆盖，则增加覆盖的碱基数量

    # 如果有覆盖的碱基，则计算平均深度
    if covered_bases > 0:
        average_depth = total_depth / covered_bases
    else:
        average_depth = 0

    return average_depth

def gen_seg(exons):
    # seg format: SEG 16387697-16387840 34.24305555555556 1
    # finaldepth equals recover depth + remap depth
    segs_strs= []
    for chrom in exons:
        for exon in exons[chrom]:
            exon_start, exon_end = exon
            remap_depth = calculate_region_depth(remap_bam, chrom, exon_start, exon_end)
            print(exon, remap_depth)
            segs_strs.append(f"SEG {exon_start}-{exon_end} {remap_depth} 1\n")
    return segs_strs

def gen_junc(remap_exon_chains):
    # JUNC 16387697-16387840 + 16390137-16390137 - 24
    # process recover exon chains
    juncs_strs = []
    juncs_dict = {}
    # for remapped bam
    for chain, count in remap_exon_chains.items():
        if len(chain) == 1:
            continue
        for i in range(len(chain) - 1):
            exon1_start, exon1_end = chain[i]
            exon2_start, exon2_end = chain[i+1]
            junc_str = f"JUNC {exon1_start}-{exon1_end} + {exon2_start}-{exon2_end} +"
            if junc_str in juncs_dict:
                juncs_dict[junc_str] += count
            else:
                juncs_dict[junc_str] = count
    for junc_str, count in juncs_dict.items():
        juncs_strs.append(f"{junc_str} {count}\n")

    return juncs_strs

def wite_graph(segs, juncs, output_file):
    with open(output_file, 'w') as f:
        for seg in segs:
            f.write(seg)
        for junc in juncs:
            f.write(junc)           

def parse_read_link():
    read_links=[]
    with open(args.read_link, 'r') as f:
        for line in f:
            string = line.strip().split(':')
            exons = string[1].split('\t')
            new_exons = []
            for exon in exons:
                start=exon.split('-')[0]
                end=exon.split('-')[1]
                exon = (int(start), int(end))
                new_exons.append(exon)
            read_links.append(new_exons)
    return read_links

def parse_chain_from_read_file(read_links, exons):
    exon_chains = {}
    for blocks in read_links:
        exon_chain = create_exon_chain(exons, chrom, blocks)
        if exon_chain == ((16390370, 16390495), (16390641, 16390887), (16390370, 16390495), (16390641, 16390887)):
            print(blocks)
        if exon_chain:
            exon_chains[exon_chain] = exon_chains.get(exon_chain, 0) + 1
    return exon_chains

def main():

    # process for each candidate inv region
    # so the chromosome is fixed
    exons = parse_gtf(args.gtf_file)
    exons = deduplicate_intervals(exons, error_margin=5)
    # print(exons)
    read_link_chains = parse_chain_from_read_file(read_links, exons)
    for chain, count in read_link_chains.items():

        print(f"Exon chain {chain} occurs {count} times")

    
    # remap_exon_chains = parse_chain_from_bam(remap_bam, exons)
    
    # # get exons in inv regions
    
    

    # # Print the exon chains and their counts
    # for chain, count in remap_exon_chains.items():
    #     print(f"Exon chain {chain} occurs {count} times")


    # segs = gen_seg(exons)
    # juncs = gen_junc(remap_exon_chains)

    # link_to_graph(recover_bam, remap_bam, recover_exon_chains, remap_exon_chains, inv_exons)


if __name__ == "__main__":
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Track exon chains in sequencing reads')
    parser.add_argument('remap_merge_bam', help='Path to the BAM file')
    parser.add_argument('gtf_file', help='Path to the GTF file')
    parser.add_argument('read_link', help='Path to the output file')

    min_intron_size=50
    max_intron_size=100000
    theta={2:0.5,3:0.25}
    
    # Parse the arguments
    args = parser.parse_args()
    remap_bam=pysam.AlignmentFile(args.remap_merge_bam, "rb")
    chrom = remap_bam.references[0]
    read_links = parse_read_link()
    

    
    # Run the main function with the provided arguments
    main()