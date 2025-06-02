from __future__ import division
import sys
import math
import pysam
import numpy as np
import random
import re
import os
import multiprocessing
import argparse
import time
from get_raw_bkp import getInsertSize


def map_ratio(read):
    # cal the mapped ratio in a read
    mapped_len = 0
    read_len = 0
    for ci in read.cigar: # record left clipped length
        if ci[0] == 0:
            mapped_len += int(ci[1])
        read_len += int(ci[1])
    return float(mapped_len)/read_len

# def filter_bam(bam_name):
#     """
#     make sure both the two ends satisfy the criteria
#     """
#     good_read_dict = {}
#     bamfile = pysam.AlignmentFile(filename = bam_name, mode = 'rb')
#     for read in bamfile.fetch():
#         if read.is_unmapped  or read.mate_is_unmapped:
#             continue 
#         if read.mapping_quality < min_q:
#             continue
#         if read.has_tag('XA'):
#             continue
#         if read.has_tag('SA'):
#             SA_tag = read.get_tag('SA')
#             SA_array = SA_tag.split(":")
#             if len(SA_array) > 2:
#             # print (read.get_tag('SA'), len(SA_array))
#                 continue
#         # if (read.reference_name == read.next_reference_name):
#         #     continue
#         # if map_ratio(read) < min_map_ratio:
#         #     continue
#         if read.query_name  not in good_read_dict:
#             good_read_dict[read.query_name] = 0
#         good_read_dict[read.query_name] += 1
#     bamfile.close()
#     return good_read_dict

def get_ref_len(ref_name):
    if len(ref_name.split(':')) < 2:
        ref_len = int(ref_name.split('_')[3])
    elif re.search("%", ref_name):
        # print (ref_name, ref_name.split('%'))
        ref_len = int(ref_name.split('%')[2]) -  int(ref_name.split('%')[1])
    else:
        ref_len = int(ref_name.split(':')[1].split('-')[1]) - int(ref_name.split(':')[1].split('-')[0])
    return ref_len

def calCrossReads(bam_name):
    # good_read_dict = filter_bam(bam_name)
    f = open(graph, "w")
    h = open(graph+".read.txt", "w")
    edge_dict = {}
    read_list = []
    bamfile = pysam.AlignmentFile(filename = bam_name, mode = 'rb')
    mean, sdev, rlen = getInsertSize(bamfile)
    rlen = int(rlen)
    insert_size = int(mean + 2*sdev) + 2 * rlen
    # test_read_name = "NZ_CP020763.1_1916841_1917378_0:0:0_0:0:0_c3ac5"
    count = 0
    near_end = 50
    for read in bamfile.fetch():
        # print (read.query_name)
        # if read.query_name == test_read_name:
        #     print (read.is_unmapped  or read.mate_is_unmapped, read.mapping_quality, map_ratio(read), read.reference_name == read.next_reference_name)
        if read.is_unmapped  or read.mate_is_unmapped:
            continue 
        if read.mapping_quality < min_q:
            continue
        if not read.has_tag('SA'):
            continue
        if read.has_tag('SA'):
            SA_tag = read.get_tag('SA')
            SA_array = SA_tag.split(":")
            print(SA_array)
            if len(SA_array) > 2:
                continue
        SA_detail = SA_tag.split(",")
        next_ref_name = SA_detail[0]
        next_start = int(SA_detail[1])
        next_dir = SA_detail[2]
        # if read.has_tag('XA'):
        #     continue
        if (read.reference_name == next_ref_name):
            continue
        # test = ["NZ_CP041725.1:1726688-1858197", "NODE_5_length_9807_cov_124.160021"]
        # if read.reference_name in test and read.next_reference_name in test:
        #     print (read.query_name)
        # if map_ratio(read) < min_map_ratio:
        #     continue
        # if len(read.reference_name.split(':')) < 2:
        #     ref_len = int(read.reference_name.split('_')[3])
        # else:
        #     ref_len = int(read.reference_name.split(':')[1].split('-')[1]) - int(read.reference_name.split(':')[1].split('-')[0])
        ref_len = get_ref_len(read.reference_name)
        mate_len = get_ref_len(read.next_reference_name)
        # if len(read.next_reference_name.split(':')) < 2:
        #     mate_len = int(read.next_reference_name.split('_')[3])
        # else:
        #     mate_len = int(read.next_reference_name.split(':')[1].split('-')[1]) - int(read.next_reference_name.split(':')[1].split('-')[0])
        # if not (abs(read.reference_start) < insert_size or abs(ref_len - read.reference_start)< insert_size):
        #     continue
        # if not (abs(read.next_reference_start) < insert_size or abs(mate_len - read.next_reference_start)< insert_size):
        #     continue
        if read.reference_name not in chrom_copy or next_ref_name not in chrom_copy:
            print ("WARNING: JUNC not in SEG!!!", read.reference_name, read.next_reference_name)
            sys.exit()
        if  chrom_copy[read.reference_name] == 0 or  chrom_copy[next_ref_name] == 0:
            # print ("The copy of one seg is zero in the edge, thus deleted.", read.reference_name, read.next_reference_name)
            continue
        # if re.search("NODE_", read.reference_name) and re.search("NODE_", read.next_reference_name): # delete edge between two contigs
        #     continue

        # if read.is_reverse and next_dir == "-":
        #if read.reference_start > next_start:
        #    refname = read.reference_name + " +"
            

        if read.is_reverse:
            refname = read.reference_name + " -"
        else:
            refname = read.reference_name + " +"
        if not read.mate_is_reverse:
            mate_ref = read.next_reference_name + " -"
        else:
            mate_ref = read.next_reference_name + " +"

        edge = refname + " " + mate_ref
        if read.query_name == "NZ_AP022525.1_781743_782294_3:0:0_1:0:0_2cd3a":
            print (read.is_reverse, read.mate_is_reverse, read.query_name, read.cigarstring, edge)
        if edge not in edge_dict:
            edge_dict[edge] = 1
        else:
            edge_dict[edge] += 1

        read_list.append([read.query_name, edge])
        
    for chrom in chrom_copy:
        # print ("#\t", chrom, chrom_copy[chrom])
        if chrom_copy[chrom] > 0:
            print ("SEG ", chrom, round(chrom_depth[chrom]), chrom_copy[chrom], 0, 1, file = f)
        else:
            print ("Copy of SEG %s is zero, thus deleted."%(chrom), round(chrom_depth[chrom]))
    
    remove_edges, my_nodes = remove_small_circle(edge_dict, chrom_copy, insert_size)
    for edge in edge_dict:
        # if edge_dict[edge] < min_edge_dp:
        #     print ("removed JUNC due to low depth ", edge, round(edge_dict[edge]))
        #     continue
        # if edge in remove_edges:
        #     print ("removed JUNC due to triangle ", edge, round(edge_dict[edge]))
        #     continue
        print ("JUNC ", edge, round(edge_dict[edge]), file = f)
        # print (edge, edge_dict[edge])
    for read_record in read_list:
        print ("READ ", read_record[0], read_record[1], file = h)

    f.close()
    h.close()

def cal_copy_number():
    f = open(depth_file, 'r')
    depth_dict = {}
    all_depth = []
    for line in f:
        array = line.strip().split()
        chrom = array[0]
        pos = int(array[1])
        depth = int(array[2])
        if chrom not in depth_dict:
            depth_dict[chrom] = []
        depth_dict[chrom].append(depth)
        all_depth.append(depth)
    f.close()
    median_depth = np.median(all_depth)
    chrom_copy = {} 
    chrom_depth = {}
    for chrom in depth_dict:
        chrom_depth[chrom] = np.median(depth_dict[chrom])
        chrom_copy[chrom] = round(float(np.median(depth_dict[chrom]))/median_depth)
    return chrom_copy, median_depth, chrom_depth

def remove_small_circle(edge_dict, chrom_copy, insert_size):
    #triangle led by the short middle seg
    seg_length_dict = get_seg_len(chrom_copy)
    # discard
    my_graph = {}
    my_nodes = {}
    for edge in edge_dict:
        if edge_dict[edge] < min_edge_dp:
            continue
        array = edge.split()
        node1 = array[0] + " " + array[1]
        node2 = array[2] + " " + array[3]
        my_nodes[array[0]] = 1
        my_nodes[array[2]] = 1
        if node1 not in my_graph:
            my_graph[node1] = [node2]
        else:
            my_graph[node1] += [node2]
    remove_edges = {}
    for node1 in my_graph:
        if len(my_graph[node1]) == 1:
            pass
        else:
            remove_edge = []
            for node2 in my_graph[node1]:
                if node2 not in my_graph:
                    continue
                else:
                    for node3 in my_graph[node2]:
                        if node3 in my_graph[node1] and seg_length_dict[node2[:-1].strip()] < insert_size:
                            edge = node1 + " " + node3
                            remove_edges[edge] = 1
    return remove_edges, my_nodes

def get_seg_len(chrom_copy):
    # calculate segment length according to its name
    # print (chrom_copy)
    seg_length_dict = {}
    for seg in chrom_copy:
        # array = seg.split("_")
        # if array[0] != "NODE":
        #     new_array = seg[:-1].split(":")[1].split("-")
        #     # print (new_array)
        #     length = abs(int(new_array[0]) - int(new_array[1]))  
        # else:
        #     length = int(array[3])
        # seg_length_dict[seg] = length
        seg_length_dict[seg] = get_ref_len(seg)
    # print (seg_length_dict)
    return seg_length_dict



if __name__ == "__main__":

    min_q = 20  # read quality cutoff
    min_dp_ratio = 0.1 #0.1 #0.3
    min_map_ratio = 0.3 # for each read
    bam_name = sys.argv[1]
    graph = sys.argv[2]
    depth_file = sys.argv[3]
    print ("start extract edges...")
    chrom_copy, median_depth, chrom_depth = cal_copy_number()
    print ("median depth:", median_depth)
    min_edge_dp = min_dp_ratio * median_depth
    calCrossReads(bam_name)
