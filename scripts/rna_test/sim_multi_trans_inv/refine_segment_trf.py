import pysam
import argparse, os
import networkx as nx
import matplotlib.pyplot as plt
import csv
import sys
from collections import defaultdict
# segment file format:
# Graph: s0	b0
# 0-0	0.000000
# 71808313-71808369	1.000000
# 71836858-71840544	9.281117
# 71840545-71842747	315.686707
# 71842748-71843644	48.985981
# 71843645-71882465	54.725292
# 0-0	0.000000
# Graph: s0	b1
# 0-0	0.000000
# 71950367-71950542	0.038393
# 71950543-71951143	0.342777
# 71951144-71951392	0.045837
# 71951393-71959666	1.029935
# 71959667-71966339	0.788215
# 71966340-71990058	1.283619
# 71990059-71990223	0.026327
# 71990224-72015005	0.932968
# 72015006-72015081	0.028392
# 0-0	0.000000

def parse_seg(f):
    segs = []
    with open(f, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('Graph:'):
                segs.append([])
            else:
                segs[-1].append(line.split('\t'))
    return segs

def parse_trf(f):
    trfs = []
    with open(f, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('Graph:'):
                trfs.append([])
            else:
                trfs[-1].append(line.split('\t'))
    return trfs

def parse_inv(inv_bed):
    # file format :
    # chr22:16390136-16390220	test:16387694-16390887,test	chr22:16390136-16390220-15.33
    # parse the first col, get the inv region, and store in a dict, key is chrom, value is a list of (start, end)
    # some duplicated inv regions, only keep one
    inv_rgs = []
    with open(inv_bed, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            inv_rg_str = parts[0]
            chrom = inv_rg_str.split(':')[0]
            start = int(inv_rg_str.split(':')[1].split('-')[0])
            end = int(inv_rg_str.split(':')[1].split('-')[1])
            # use set to store the inv regions
            inv_rgs.append((chrom, start, end))

    return inv_rgs

def get_intersection(interval_a, interval_b):
    """
    Returns the intersection of two intervals, if they intersect.
    Each interval is a tuple of (start, end).
    """
    a_start, a_end = interval_a
    b_start, b_end = interval_b

    # Find the start and end of the intersection
    start = max(a_start, b_start)
    end = min(a_end, b_end)

    # Check if there is an intersection
    if start < end:
        return (start, end)
    else:
        return None

def filter_segs(segs, inv_rg):
    new_segs = []
    resered_ids = []
    for graph_id, graph in enumerate(segs):
        for seg in graph:
            start = int(seg[0].split('-')[0])
            end = int(seg[0].split('-')[1])
            if get_intersection((start, end),(inv_rg[0][1], inv_rg[0][2])):
                new_segs.append(graph)
                resered_ids.append(graph_id)
                break
    return new_segs, resered_ids

def filter_trfs(trfs, resered_ids):
    new_trfs = []
    for graph_id, graph in enumerate(trfs):
        if graph_id in resered_ids:
            new_trfs.append(graph)
    return new_trfs
            

def merge_segs(segs):
    depth_dict = defaultdict(float)
    for seg in segs:
        for s in seg:
            depth_dict[s[0]] += float(s[1])
    print(depth_dict)

def merge_trfs(trfs):
    depth_dict = defaultdict(float)
    for trf in trfs:
        for t in trf:
            depth_dict[t[1]] += float(t[2])
    print(depth_dict)



def main():
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    segs = parse_seg(args.segment_f)
    trfs = parse_trf(args.trf_f)
    inv_rg = parse_inv(args.inv_bed)
    new_segs, resered_ids=filter_segs(segs, inv_rg)
    new_trfs=filter_trfs(trfs, resered_ids)

    for i, seg in enumerate(new_segs):
        with open(os.path.join(args.outdir, f"seg_{i}.txt"), 'w') as f:
            f.write(f"Graph: s0\tb{i}\n")
            for s in seg:
                f.write(f"{s[0]}\t{s[1]}\n")
    for i, trf in enumerate(new_trfs):
        with open(os.path.join(args.outdir, f"trf_{i}.txt"), 'w') as f:
            f.write(f"Graph: s0\tb{i}\n")
            for t in trf:
                f.write(f"{t[0]}\t{t[1]}\t{t[2]}\n")            

    
    print(new_segs)
    print(new_trfs)
    print(inv_rg)

if __name__ == "__main__":
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Track exon chains in sequencing reads')
    parser.add_argument('segment_f', help='Path to the BAM file')
    parser.add_argument('trf_f', help='Path to the output file')
    parser.add_argument('segment_o_f', help='Path to the GTF file')
    parser.add_argument('trf_o_f', help='Path to the BAM file')
    parser.add_argument('inv_bed', help='Chromosome name')
    parser.add_argument('outdir', )
    args = parser.parse_args()
    main()