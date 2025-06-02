import argparse, os
import sys
from collections import defaultdict


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

def main():
    segs=parse_seg(args.segment_f)
    trfs=parse_trf(args.trf_f)
    segs_complete=[]
    all_segs=[]
    prev_seg=None
    for idx, seg in enumerate(segs[0]):
        if prev_seg is None or (int(seg[0].split('-')[0]) - int(prev_seg[0].split('-')[1])) >3:
            all_segs.append(segs_complete)
            segs_complete=[]
            segs_complete.append(seg)
            prev_seg=seg
        elif (int(seg[0].split('-')[0]) - int(prev_seg[0].split('-')[1])) <=3:
            segs_complete.append(seg)
            prev_seg=seg
    


    print("original exons:")
    for seg in segs[0]:
        print(seg)

    print("merged exons:")
    for seg in all_segs:
        print(seg)
    # print(trfs)
    # print(len(segs[0]), len(trfs[0]))



if __name__ == "__main__":
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Track exon chains in sequencing reads')
    parser.add_argument('segment_f', help='Path to the BAM file')
    parser.add_argument('trf_f', help='Path to the output file')
    parser.add_argument('merged_segment_o_f', help='Path to the GTF file')
    parser.add_argument('merged_trf_o_f', help='Path to the BAM file')
    args = parser.parse_args()
    main()