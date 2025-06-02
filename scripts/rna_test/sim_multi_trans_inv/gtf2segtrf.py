
import argparse, os
from collections import defaultdict
import networkx as nx

def get_overlap(interval1, interval2):
    """
    返回两个区间的重叠部分。

    参数:
    interval1 -- 第一个区间，形式为(start1, end1)
    interval2 -- 第二个区间，形式为(start2, end2)

    返回:
    一个元组表示的重叠区间(start, end)。如果没有重叠，返回None。
    """

    # 分别获取两个区间的起始和结束位置
    start1, end1 = interval1
    start2, end2 = interval2

    # 计算重叠区间的起始和结束位置
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)

    # 如果计算出的起始位置大于结束位置，则表示没有重叠
    if overlap_start >= overlap_end:
        return None

    # 返回重叠区间
    return (overlap_start, overlap_end)

def get_overlaps(interval1, interval_list, visited_invs):
    """
    返回 interval1 和 interval_list 中任意区间的重叠部分。

    参数:
    interval1 -- 单一区间，形式为(start1, end1)
    interval_list -- 区间列表，形式为[(start2, end2), ...]

    返回:
    一个元组表示的重叠区间(start, end)。如果没有重叠，返回None。
    """

    # 分别获取 interval1 的起始和结束位置
    start1, end1 = interval1

    # 遍历 interval_list 中的每个区间
    for start2, end2 in interval_list:
        # 计算当前区间与 interval1 的重叠部分
        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)

        # 如果计算出的起始位置小于结束位置，则存在重叠
        if overlap_start < overlap_end:
            visited_invs.append((start2, end2))
            return (overlap_start, overlap_end), visited_invs

    # 如果所有区间都没有重叠，返回 None
    return None, []

def parse_gtf(gtf, inv_rgs):
    trans_fregments = []
    trans_fregments_tmp=[]
    mormal_exons = []
    trf_covs = []
    inv_region=inv_rgs[chrom_global]
    visited_invs=[]
    with open(gtf, 'r') as file:
        for line_id, line in enumerate(file):
            line = line.strip()
            if line.startswith('#'):
                continue
            line = line.split('\t')
            if line[2] == 'exon':
                start= int(line[3])
                end = int(line[4])
                # print((start, end), inv_rgs[chrom_global])
                overlap, visited_invs= get_overlaps((start, end), inv_region, visited_invs)
                # print(f"processing exon: ({start}, {end}), {inv_region}, {overlap}")
                if overlap:
                    if inv_region[0]>= start and inv_region[1]>=end:
                        exon1_start=start
                        exon1_end=inv_region[0]
                        exon2_start=overlap[0]
                        exon2_end=overlap[1]
                        exon3_start=end
                        exon3_end=inv_region[1]
                    elif inv_region[0]>= start and inv_region[1]<=end:
                        exon1_start=start
                        exon1_end=overlap[0]
                        exon2_start=overlap[0]
                        exon2_end=overlap[1]
                        exon3_start=overlap[1]
                        exon3_end=end
                    elif inv_region[0]<= start and inv_region[1]<=end:
                        exon1_start=inv_region[0]
                        exon1_end=overlap[0]
                        exon2_start=overlap[0]
                        exon2_end=overlap[1]
                        exon3_start=inv_region[1]
                        exon3_end=end
                    else:
                        print("may miss inv type")
                        pass
                    if exon1_end-exon1_start>=5:
                        mormal_exons.append([line[0], line[1], line[2], str(exon1_start), str(exon1_end), line[5], line[6], line[7], line[8]])
                        trans_fregments_tmp.append([line[0], line[1], line[2], str(exon1_start), str(exon1_end), line[5], line[6], line[7], line[8]])
                    if exon2_end-exon2_start>=5:
                        mormal_exons.append([line[0], line[1], line[2], str(exon2_start), str(exon2_end), line[5], line[6], line[7], line[8]])
                        trans_fregments_tmp.append([line[0], line[1], line[2], str(exon2_start), str(exon2_end), line[5], line[6], line[7], line[8]])
                    if exon3_end-exon3_start>=5:
                        mormal_exons.append([line[0], line[1], line[2], str(exon3_start), str(exon3_end), line[5], line[6], line[7], line[8]])
                        trans_fregments_tmp.append([line[0], line[1], line[2], str(exon3_start), str(exon3_end), line[5], line[6], line[7], line[8]])
                    print("inv split exon:")
                    print(exon1_start, exon1_end)
                    print(exon2_start, exon2_end)
                    print(exon3_start, exon3_end)
                else:
                    mormal_exons.append([line[0], line[1], line[2], start, end, line[5], line[6], line[7], line[8]])
                    trans_fregments_tmp.append([line[0], line[1], line[2], start, end, line[5], line[6], line[7], line[8]])    
            if line[2] == 'transcript':
                trf_covs.append(float(line[8].split(';')[2].split(' ')[2].strip('"')))
                if line_id == 2:
                    continue
                trans_fregments.append(trans_fregments_tmp)
                trans_fregments_tmp = []
        trans_fregments.append(trans_fregments_tmp)
    return mormal_exons, trans_fregments, trf_covs


def find_overlapping_intervals(intervals):
    """
    Finds and returns a list of overlapping interval pairs from a list of intervals.
    Each interval is represented as a tuple of (start, end), and both start and end are integers.
    """
    overlapping_pairs = []

    # Sort intervals by start time
    intervals.sort(key=lambda x: x[0])

    # Compare each interval with every other interval to check for overlap
    for i in range(len(intervals)):
        for j in range(i + 1, len(intervals)):
            # If the start of the second interval is less than the end of the first, we have an overlap
            if intervals[j][0] < intervals[i][1]:
                overlapping_pairs.append((intervals[i], intervals[j]))
            else:
                # Because the intervals are sorted, we can break the inner loop once we find a non-overlapping interval
                break

    return overlapping_pairs

def parse_inv(inv_bed):
    # file format :
    # chr22:16390136-16390220	test:16387694-16390887,test	chr22:16390136-16390220-15.33
    # parse the first col, get the inv region, and store in a dict, key is chrom, value is a list of (start, end)
    # some duplicated inv regions, only keep one
    inv_rgs = {}
    inv_exon_cov_dict = {}
    global chrom_global
    with open(inv_bed, 'r') as f:
        for line in f:
            chrom, start, end, cov_str = line.strip().split('\t')
            chrom_global = chrom
            print("chrom:", chrom_global)
            inv_rgs.setdefault(chrom, list()).append((int(start), int(end)))
            inv_exon_cov_dict[(int(start), int(end))] = float(cov_str)

    return inv_rgs, inv_exon_cov_dict

def main():
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    inv_rgs, inv_exon_cov_dict = parse_inv(args.inv_bed)
    print("inv_rgs:", inv_rgs)
    mormal_exons, trans_fregments, trf_covs = parse_gtf(args.gtf, inv_rgs)

    
    # for idx, trans in enumerate(trans_fregments):
    #     print(f"{idx+1} th transcripts:-------------------------------------->")
    #     for t in trans:
    #         print(t)
    #         # print(f"{t[3]}-{t[4]}")
    #     print("End ----------------------------")

    # sort normal exons by start position
    mormal_exons.sort(key=lambda x: int(x[3]))
    # print("normal exons:", mormal_exons)
    # for normal_exon in mormal_exons:
    #     print(f"{normal_exon}")
    
    
    exon_cov_dict = defaultdict(float)
    for exon_str in mormal_exons:
        exon=(int(exon_str[3]), int(exon_str[4]))
        exon_cov=exon_str[8].split(';')[3].split(' ')[2]
        exon_cov_dict[exon] += float(exon_cov.strip('"'))
    normal_exons_tuples = []
    print(exon_cov_dict)
    for exon_tuple in exon_cov_dict.keys():
        normal_exons_tuples.append(exon_tuple)

    for exon in inv_exon_cov_dict.keys():
        normal_exons_tuples.append(exon)
        exon_cov_dict[exon] = inv_exon_cov_dict[exon]

    # sort exon_cov_dict key by start position
    exon_cov_dict = dict(sorted(exon_cov_dict.items(), key=lambda item: item[0][0]))
        

    # sort normal_exons_tuples by start position
    normal_exons_tuples.sort(key=lambda x: x[0])
    
    # for seg in normal_exons_tuples:
    #     print(f"{seg}")
    print("len segs:", len(normal_exons_tuples))
    print(len(trf_covs), len(trans_fregments))

    # build nextworkx by normal exons and trans_fregments


    # # filter graph here
    # G = nx.Graph()
    # for i, exon in enumerate(normal_exons_tuples):
    #     G.add_node(f"{exon[3]}-{exon[4]}", cov=exon_cov_dict[(exon[3], exon[4])])
    # for trans_i, trans in enumerate(trans_fregments):
    #     # Assuming trans is a list of tuples with exon start and end positions
    #     for i in range(len(trans) - 1):
    #         # Connect consecutive transcribed fragments with an edge
    #         in_node=(trans[i][4],trans[i][5])
    #         out_node=(trans[i+1][4],trans[i+1][5])
    #         attributes = {
    #             "used_trf_id":trans_i,
    #         }
    #         G.add_edge(in_node, out_node, **attributes)





    # write graph
    with open(args.segment_o_f, 'w') as file:
        file.write(f"Graph: s0\tb0\n")
        file.write(f"0-0\t{0.0000}\n")
        for exon, cov in exon_cov_dict.items():
            file.write(f"{exon[0]}-{exon[1]}\t{cov}\n")
        file.write(f"0-0\t{0.0000}\n")

    with open(args.trf_o_f, 'w') as file:
        file.write(f"Graph: s0\tb0\n")
        trans_f_list=[]
        for i, trans in enumerate(trans_fregments):
            trans_f_order_tmp=[]
            for trans_exon_str in trans:
                exon_tuple = (int(trans_exon_str[3]), int(trans_exon_str[4]))
                index = normal_exons_tuples.index(exon_tuple)+1
                trans_f_order_tmp.append(index)
            trans_f_list.append(['1',trans_f_order_tmp, trf_covs[i]])
        # sort trans_f_list by length of trans_f_order_tmp
        trans_f_list.sort(key=lambda x: len(x[1]), reverse=True)
        for trans_f in trans_f_list:
            file.write(f"{1}\t0,")
            for trans in trans_f[1]:
                file.write(f"{trans},")
            file.write(f"{len(normal_exons_tuples)+1}")
            file.write(f"\t{trans_f[2]}\n")
        # source2first_dict = defaultdict(float)
        # for i, trans in enumerate(trans_fregments):
        #     first_exon=(int(trans[0][3]), int(trans[0][4]))
        #     first_exon_index = normal_exons_tuples.index(first_exon)+1
        #     source2first_dict[f"0,{first_exon_index},"] += trf_covs[i]
        #     continue
        # for key, value in source2first_dict.items():
        #     file.write(f"1\t{key}\t{value}\n")

    os.system(f"mv {args.trf_o_f} {args.outdir}")
    os.system(f"mv {args.segment_o_f} {args.outdir}")


if __name__ == "__main__":
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Track exon chains in sequencing reads')
    parser.add_argument('gtf', help='stringtie gtf')
    parser.add_argument('segment_o_f', help='Path to the GTF file')
    parser.add_argument('trf_o_f', help='Path to the BAM file')
    parser.add_argument('inv_bed', help='Chromosome name')
    parser.add_argument('outdir', )
    chrom_global=None
    args = parser.parse_args()
    main()