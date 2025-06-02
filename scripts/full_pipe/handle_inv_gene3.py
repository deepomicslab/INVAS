import sys
import pysam
import re
from Bio.Seq import Seq
import parasail
import bisect 
from pybedtools import BedTool
from pyssw import align
import argparse as ap

# def get_reads_dict(bam_file):
#     with pysam.AlignmentFile(bam_file, "rb") as file:
#         reads_dict = {}
#         for read in file:
#             if read.query_name in reads_dict:
#                 reads_dict[read.query_name].append(read)
#             else:
#                 reads_dict[read.query_name] = [read]
#         return reads_dict
    
def get_reads_dict(bam_f, region):
    reads_dict = {}
    for read in bam_f.fetch(region[0], region[1], region[2]):
        if read.query_name in reads_dict:
            reads_dict[read.query_name].append(read)
        else:
            reads_dict[read.query_name] = [read]
    print("bam recs count:", len(reads_dict))
    return reads_dict


def update_flags(read1, read2):
    if read2.is_read1:
        read1.is_read1 = False
        read1.is_read2 = True
    else:
        read1.is_read1 = True
        read1.is_read2 = False

    # # Remove the unmapped flag
    # if read1.is_unmapped:
    #     read1.flag &= ~0x4  # Clear the bit for unmapped
    # if read2.is_unmapped:
    #     read2.flag &= ~0x4  # Clear the bit for unmapped
    # Remove the mate unmapped flag
    read1.flag &= ~0x8  # Clear the bit for mate unmapped
    read2.flag &= ~0x8  # Clear the bit for mate unmapped
    return read1, read2

def copy_read(read):
    new_read = pysam.AlignedSegment()
    new_read.query_name = read.query_name
    new_read.query_sequence = read.query_sequence
    new_read.flag = read.flag
    new_read.reference_id = read.reference_id
    new_read.reference_start = read.reference_start
    new_read.mapping_quality = read.mapping_quality
    new_read.cigar = read.cigar
    new_read.next_reference_id = read.next_reference_id
    new_read.next_reference_start = read.next_reference_start
    new_read.template_length = read.template_length
    new_read.query_qualities = read.query_qualities
    new_read.tags = read.tags
    return new_read

def write_reads_to_bam(output_f, header, new_splice_rec, reads):
    # Open the output file in write mode ('w') with BAM format ('b')
    # with pysam.AlignmentFile(output_filename, 'wb', header=header) as output_bam:
        # Write new_splice_rec to the output file
    output_f.write(new_splice_rec)
    output_f.write(reads[0])

from Bio import SeqIO

def fastq_to_dict(fastq_file):
    seq_dict = SeqIO.to_dict(SeqIO.parse(fastq_file, "fastq"))
    return seq_dict

import re

def adjust_cigar_str(merged_cigar_str, read_length, total_match_cnt):
    """
    Adjust the last M value in the CIGAR string based on read length and total match count.

    Parameters:
    merged_cigar_str (str): The original CIGAR string.
    read_length (int): The total length of the read.
    total_match_cnt (int): The total count of matches in the CIGAR string.

    Returns:
    str: The adjusted CIGAR string.
    """

    read_len_difference = read_length - total_match_cnt

    if read_len_difference != 0:
        # Find the last M in the merged_cigar_str
        last_M_index = merged_cigar_str.rfind('M')
        # Find the start index of the number before the last M
        last_num_start_index = max(merged_cigar_str.rfind(c, 0, last_M_index) for c in 'MIDNSHP=X')
        if last_num_start_index == -1:
            last_num_start_index = 0
        else:
            last_num_start_index += 1
        # Extract the number before the last M
        last_M_number = int(merged_cigar_str[last_num_start_index:last_M_index])
        # Increase or decrease the last M number
        if read_len_difference > 0:
            last_M_number += read_len_difference
        else:
            last_M_number -= abs(read_len_difference)
        # Replace the last M number in the merged_cigar_str
        merged_cigar_str = merged_cigar_str[:last_num_start_index] + str(last_M_number) + merged_cigar_str[last_M_index:]

    return merged_cigar_str

def print_common_reads(reads_dict1, reads_dict2):
    common_read_names = set(reads_dict1.keys()) & set(reads_dict2.keys())
    # get reads in bam1 not in bam2
    not_common_read_names = set(reads_dict1.keys()) - set(reads_dict2.keys())

    print("common reads length:", len(common_read_names))
    print("not common reads length:", len(not_common_read_names))
    passed_cnt=0
    not_passed_cnt=0

    for read_name in common_read_names:
        # print(f"Read name {read_name} found in both files:")
        # new_splice_rec = reads_dict2[read_name][0].copy()
        if len(reads_dict1[read_name])<=1:
            not_passed_cnt+=1
            # print("unique read:", reads_dict1[read_name][0])
            continue
        passed_cnt+=1
        new_splice_rec = copy_read(reads_dict2[read_name][0])
        reads_dict1[read_name].sort(key=lambda a: a.reference_start)
        merged_end = 0
        merged_cigar_str= ""
        print("before transform:##################################################################33")
        # for idx, alignment in enumerate(reads_dict1[read_name]):
        #     print("bwa:", alignment)
        # print("hisat:", reads_dict2[read_name][0])
        total_match_cnt=0
        ## todo get read length
        read_length = 100
        for idx, alignment in enumerate(reads_dict1[read_name]):
            # get number of N in cigar
            cigarstring = alignment.cigarstring
            match = re.search(r'(\d+)M', cigarstring)
            match_cnt = int(match.group(1))
            total_match_cnt+=match_cnt
            
            print("read_len:",read_length)
            if idx == 0:
                merged_cigar_str = f"{match_cnt}M"
                new_splice_rec.reference_start = alignment.reference_start
            else:
                splice_length = alignment.reference_start - merged_end
                merged_cigar_str += f"{splice_length}N{match_cnt}M"
            # Update end position
            merged_end = alignment.reference_end
        print("read len, total match:", read_length, total_match_cnt, read_length-total_match_cnt)
        
        merged_cigar_str = adjust_cigar_str(merged_cigar_str, read_length, total_match_cnt)
        print("merged_cigar_str",merged_cigar_str)

        new_splice_rec.cigarstring = merged_cigar_str
        # new_splice_rec.reference_end = merged_end
        new_splice_rec, reads_dict2[read_name][0] = update_flags(new_splice_rec, reads_dict2[read_name][0])
        sequence = reads_dict2[read_name][0].query_sequence

        sequence = Seq(sequence)
        reverse_complement = str(sequence.reverse_complement())
        new_splice_rec.query_sequence = reverse_complement
        print("after transform:##################################################################33")
        # print("new_rec:",new_splice_rec)
        # print("old_rec:",reads_dict2[read_name][0])
        # Write the new record to the output file
        write_reads_to_bam(outf, reads_dict2[read_name][0].header, new_splice_rec, reads_dict2[read_name])
        del reads_dict2[read_name]
    print("passed_cnt:", passed_cnt)
    print("not_passed_cnt:", not_passed_cnt)
    print("total:", passed_cnt+not_passed_cnt)
    for rnames, rec in reads_dict2.items():
        for r in rec:
            
            outf.write(r)


    # for not mapped reads in bam1
    for read_name in not_common_read_names:
        for r in reads_dict1[read_name]:
            if 'S' in r.cigarstring:
                print("bam notmap with S:",r)
            outf.write(r)


import re

def parse_cigar(cigar_string):
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    matches = pattern.findall(cigar_string)

    left_soft_clip = 0
    right_soft_clip = 0

    if matches:
        # Check the first element
        if matches[0][1] == 'S':
            left_soft_clip = int(matches[0][0])

        # Check the last element
        if matches[-1][1] == 'S':
            right_soft_clip = int(matches[-1][0])

    return left_soft_clip, right_soft_clip

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

def remove_soft_clipping(cigar_string, remove_left=True, remove_right=True):
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    matches = pattern.findall(cigar_string)

    if matches:
        # Check and remove the first element if it's a soft clipping
        if remove_left and matches[0][1] == 'S':
            matches.pop(0)

        # Check and remove the last element if it's a soft clipping
        if remove_right and matches and matches[-1][1] == 'S':
            matches.pop()

    # Reconstruct the cigar string
    new_cigar_string = ''.join([f'{num}{op}' for num, op in matches])

    return new_cigar_string

def remove_soft_clipping(cigar_string, remove_left=True, remove_right=True):
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    matches = pattern.findall(cigar_string)

    if matches:
        # Check and remove the first element if it's a soft or hard clipping
        if remove_left and matches[0][1] in ['S', 'H']:
            matches.pop(0)

        # Check and remove the last element if it's a soft or hard clipping
        if remove_right and matches and matches[-1][1] in ['S', 'H']:
            matches.pop()

    # Reconstruct the cigar string
    new_cigar_string = ''.join([f'{num}{op}' for num, op in matches])

    return new_cigar_string

def mirror_cigar(cigar_string):
    # Regular expression pattern to parse the CIGAR string
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')

    # Parse the CIGAR string
    matches = pattern.findall(cigar_string)

    # Reverse the matches
    matches.reverse()

    # Reconstruct the CIGAR string
    mirrored_cigar = ''.join([str(length) + op for length, op in matches])

    return mirrored_cigar

def parse_fasta_header(header):
    # Assuming the header is in the format '>chr:start-end'
    chrom, start, end = header.split(" ")[0].split('-')
    return int(start), int(end)

def binary_search(headers, target):
    # Sort headers based on the start position
    headers.sort(key=parse_fasta_header)
    starts = [parse_fasta_header(header)[0] for header in headers]
    ends = [parse_fasta_header(header)[1] for header in headers]

    index = bisect.bisect(starts, target)

    # If the target location is within a sequence
    if index > 0 and target <= ends[index - 1]:
        prev_header = headers[index - 2] if index > 1 else None
        curr_header = headers[index - 1]
        next_header = headers[index] if index < len(headers) else None
    else:
        # If the target is not within a sequence, get the previous and next headers
        prev_header = headers[index - 1] if index > 0 else None
        curr_header = None
        next_header = headers[index] if index < len(headers) else None

    return prev_header, curr_header, next_header


def get_reverse(sequence):
    """
    This function takes a DNA sequence as input and returns its reverse.
    """
    return sequence[::-1]

def get_complement(sequence):
    """
    This function takes a DNA sequence as input and returns its complement.
    """
    seq = Seq(sequence)
    return str(seq.complement())

def get_reverse_complement(sequence):
    """
    This function takes a DNA sequence as input and returns its reverse complement.
    """
    seq = Seq(sequence)
    return str(seq.reverse_complement())

def get_query_seq(read, soft_clip_len, direction):
    query_seq = str(r1_dict[read.query_name].seq) if read.is_read1 else str(r2_dict[read.query_name].seq)
    # print("Query seq:", query_seq, type(query_seq))
    if direction == "l" and read.is_reverse:
        return get_reverse(query_seq[:soft_clip_len])
    elif direction == "l" and not read.is_reverse:
        return get_reverse_complement(query_seq[:soft_clip_len])
    elif direction == "r" and read.is_reverse:
        return get_reverse(query_seq[-soft_clip_len:])
    elif direction == "r" and not read.is_reverse:
        return get_complement(query_seq[-soft_clip_len:])
    return ""

def get_query_seq_for_inv_rg(read, soft_clip_len, direction):
    # query_seq = str(r1_dict[read.query_name].seq) if read.is_read1 else str(r2_dict[read.query_name].seq)
    # query_seq = str(umpr_dict[read.query_name].seq) 
    # if direction == "l" and read.is_reverse:
    #     return get_reverse(query_seq[:soft_clip_len])
    # elif direction == "l" and not read.is_reverse:
    #     return get_reverse_complement(query_seq[:soft_clip_len])
    # elif direction == "r" and read.is_reverse:
    #     return query_seq[-soft_clip_len:]
    # elif direction == "r" and not read.is_reverse:
    #     return get_complement(query_seq[-soft_clip_len:])
    # if direction == "l":
    #     return query_seq[:soft_clip_len]
    # elif direction == "r":
    #     return query_seq[-soft_clip_len:]
    # print("truth seq", read.query_sequence, read.query_sequence[:soft_clip_len][::-1])
    if direction == "l":
        return read.query_sequence[:soft_clip_len][::-1]
    elif direction == "r":
        # print(read.query_sequence, read.query_sequence[-soft_clip_len:][::-1])
        return read.query_sequence[-soft_clip_len:][::-1]
    return ""

    

def get_ref_seq(read, direction):
    prev_header,cur_header, next_header = binary_search(exon_headers, read.reference_start)
    print(read)
    print(prev_header,cur_header, next_header, read.reference_start)
    if direction == "r" and not read.is_reverse:
        return int(next_header.split(' ')[0].split("-")[1]), get_reverse(str(exon_db[next_header.split(' ')[0]].seq))
    elif direction == "r" and read.is_reverse:
        return int(next_header.split(' ')[0].split("-")[1]), get_reverse(str(exon_db[next_header.split(' ')[0]].seq))
    elif direction == "l" and not read.is_reverse:
        print("prev_header:", prev_header)
        return int(prev_header.split(' ')[0].split("-")[2]), str(exon_db[prev_header.split(' ')[0]].seq)
    elif direction == "l" and read.is_reverse:
        return int(prev_header.split(' ')[0].split("-")[2]), str(exon_db[prev_header.split(' ')[0]].seq)
    return 0, ""

def get_ref_seq_for_inv_rg(read, direction):
    prev_header,cur_header, next_header = binary_search(exon_headers, read.reference_start)
    print(read)
    print(prev_header,cur_header, next_header, read.reference_start)
    # if direction == "r" and not read.is_reverse:
    #     return int(prev_header.split(' ')[0].split("-")[1]), int(prev_header.split(' ')[0].split("-")[2]), get_reverse(str(exon_db[prev_header.split(' ')[0]].seq)), cur_header
    # elif direction == "r" and read.is_reverse:
    #     return int(prev_header.split(' ')[0].split("-")[1]), int(prev_header.split(' ')[0].split("-")[2]), get_reverse(str(exon_db[prev_header.split(' ')[0]].seq)), cur_header
    # elif direction == "l" and not read.is_reverse:
    #     return int(next_header.split(' ')[0].split("-")[1]), int(next_header.split(' ')[0].split("-")[2]), str(exon_db[next_header.split(' ')[0]].seq), cur_header
    # elif direction == "l" and read.is_reverse:
    #     return int(next_header.split(' ')[0].split("-")[1]), int(next_header.split(' ')[0].split("-")[2]), str(exon_db[next_header.split(' ')[0]].seq), cur_header
    
    if direction == "r" and not read.is_reverse:
        return int(prev_header.split(' ')[0].split("-")[1]), int(prev_header.split(' ')[0].split("-")[2]), str(exon_db[prev_header.split(' ')[0]].seq), cur_header
    elif direction == "r" and read.is_reverse:
        return int(prev_header.split(' ')[0].split("-")[1]), int(prev_header.split(' ')[0].split("-")[2]), str(exon_db[prev_header.split(' ')[0]].seq), cur_header
    elif direction == "l" and not read.is_reverse:
        return int(next_header.split(' ')[0].split("-")[1]), int(next_header.split(' ')[0].split("-")[2]), str(exon_db[next_header.split(' ')[0]].seq), cur_header
    elif direction == "l" and read.is_reverse:
        return int(next_header.split(' ')[0].split("-")[1]), int(next_header.split(' ')[0].split("-")[2]), str(exon_db[next_header.split(' ')[0]].seq), cur_header
    return 0,0, "", ""




def get_last_match_position_in_region(cigar_string, read_start_position, region_start, region_end):
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    matches = pattern.findall(cigar_string)

    current_position = read_start_position
    last_match_position = None

    for length, operation in matches:
        if operation == 'M':
            for i in range(int(length)):
                if region_start <= current_position <= region_end:
                    last_match_position = current_position
                current_position += 1
        elif operation in ['I', 'S']:  # Insertion or soft clipping does not advance in the reference
            continue
        else:  # Other operations ('D', 'N', 'H', 'P', '=', 'X') do advance in the reference
            current_position += int(length)

    return last_match_position



def best_alignment(query, ref, open_penalty=10, extend_penalty=1, matrix=parasail.blosum62):
    # 对query序列进行反向互补转换
    # query_revcomp = str(Seq(query).reverse_complement())
    ref_c = get_complement(ref)

    # 对原始query序列进行比对
    # result = parasail.ssw(query, ref, open_penalty, extend_penalty, matrix)

    # 对反向互补的query序列进行比对
    result_refc = parasail.ssw(query, ref_c, open_penalty, extend_penalty, matrix)

    # print("result:", result.ref_begin1, result.ref_end1, result.score1)
    print("result_refc:", result_refc.ref_begin1, result_refc.ref_end1, result_refc.score1)

    # 比较两个方向的比对分数，并选择分数最高的那个
    # if result.score1 >= result_revcomp.score1:
    #     return result
    # else:
    #     return result_revcomp

def best_alignment2(query, ref):
    # 对query序列进行反向互补转换
    # query_revcomp = str(Seq(query).reverse_complement())
    # ref_revc = get_reverse_complement(ref)
    ref_c = get_complement(ref)

    # 对原始query序列进行比对
    # result = align(query, ref)

    # 对反向互补的query序列进行比对

    # result_rc = align(query, ref, True)
    result_c = align(query, ref_c)

    # result_rcc = align(query, ref_c, True)


    # 比较两个方向的比对分数，并选择分数最高的那个
    # if result_c['score'] >= result['score']:
    #     return result_c
    # else:
    #     return result
    return result_c


def get_sequence_from_dict(sequence_dict, region):
    """
    从字典中获取指定区域的序列

    :param sequence_dict: 存储序列的字典
    :param region: 字符串，格式为'chr:start-end'
    :return: 指定区域的序列
    """
    chrom, start, end = region
    # print(chrom, start, end)
    if chrom in sequence_dict:
        sequence = sequence_dict[chrom].seq[start-1:end]
        return str(sequence)
    else:
        return None
    
def read_fasta_to_dict(fasta_file):
    """
    读取FASTA文件，并将序列存储在字典中

    :param fasta_file: FASTA文件的路径
    :return: 一个字典，键是染色体名，值是序列
    """
    sequence_dict = {}
    with open(fasta_file, 'r') as file:
        for record in SeqIO.parse(file, "fasta"):
            sequence_dict[record.id] = record.seq
    return sequence_dict

def parse_regions(regions_file):
    """_summary_

    read file, then parse the first column as region name
    return a list of region names,in format of (chr, start, end)
    chr12:71139013-71139973 -> (chr12, 71139013, 71139973)
    """
    regions = []
    gene_regions = []
    real_inv_regions = []
    with open(regions_file, 'r') as file:
        for line in file:
            line = line.strip()
            region_str = line.split('\t')[0]
            chrom, positions = region_str.split(':')
            start, end = map(int, positions.split('-'))
            regions.append((chrom, start, end))
            gene_region_str = line.split('\t')[1].split(',')[0].split(':')[1]
            gene_start, gene_end = map(int, gene_region_str.split('-'))
            gene_regions.append((chrom, gene_start, gene_end))
            real_inv_region_str = line.split('\t')[2].split(':')[1]
            print("real_inv_region_str:", real_inv_region_str)
            real_inv_start, real_inv_end = map(int, real_inv_region_str.split('-')[:2])
            real_inv_regions.append((chrom, real_inv_start, real_inv_end))
    return regions, gene_regions, real_inv_regions


def gen_cigar_left_clip(cigarstring, still_soft_len, left_tag, N_len, inter_gap_len, match_len):
    if inter_gap_len > 0:
        cigarstring = f"{cigarstring}{inter_gap_len}I{N_len}N{match_len}M"
        if still_soft_len > 0:
            cigarstring = f"{cigarstring}{still_soft_len}{left_tag}"
    else:
        cigarstring = f"{cigarstring}{N_len}N{match_len}M"
        if still_soft_len > 0:
            cigarstring = f"{cigarstring}{still_soft_len}{left_tag}" 
    return cigarstring

def gen_cigar_right_clip(cigarstring, still_soft_len, right_tag, N_len, inter_gap_len, match_len):
    if inter_gap_len > 0:
        cigarstring = f"{match_len}M{N_len}N{inter_gap_len}I{cigarstring}"
        if still_soft_len > 0:
            cigarstring = f"{still_soft_len}{right_tag}{cigarstring}"
    else:
        cigarstring = f"{match_len}M{N_len}N{cigarstring}"
        if still_soft_len > 0:
            cigarstring = f"{still_soft_len}{right_tag}{cigarstring}" 
    return cigarstring



def process_bwa2(reads_dict1, selected_region, gene_region, real_inv_region):
    for rname, rec in reads_dict1.items():
        for idx, r in enumerate(rec):
            # if r.query_name not in ["UNC11-SN627:379:C5FDFACXX:3:2108:20340:78649_sg"]:
            #     continue
            left_ext_region = (real_inv_region[0], gene_region[1], real_inv_region[2])
            right_ext_region = (real_inv_region[0], real_inv_region[1], gene_region[2])
            left_ref = get_sequence_from_dict(ref_db, left_ext_region)
            right_ref = get_sequence_from_dict(ref_db, right_ext_region)
            # print(left_ext_region, right_ext_region)
            # if r.query_name != 'read42/ENST00000399687;mate1:204-303;mate2:350-448_sg':
            #     continue
            # print(r)
            # print(r.query_name, r.cigarstring, r.reference_start, r.reference_end, r.is_reverse, r.is_read1, r.is_read2, r.is_supplementary)
            if 'S' in r.cigarstring and not r.is_supplementary:
                r_cp = copy_read(r)
                print("##########################################################################################")
                print("before:", r)
                print("##########################################################################################")
                left_soft_clip, left_tag, right_soft_clip, right_tag = parse_cigar2(r.cigarstring)
                print(r.cigarstring, f'Left soft clip: {left_soft_clip}', f'Right soft clip: {right_soft_clip}' )
                if left_soft_clip > clip_threshold and right_soft_clip <=clip_threshold:
                    left_soft_seq = get_query_seq_for_inv_rg(r, left_soft_clip, "l")
                    print("left_soft_seq:", left_soft_seq, type(left_soft_seq), type(right_ref))
                    # print("right_exon_seq:", right_exon_seq)
                    align_res = best_alignment2(left_soft_seq, right_ref)
                    print("test_align res:",align_res['ref_begin'], align_res['ref_end'])
                    print("real_inv_region:", real_inv_region)
                    if align_res:
                        match_len=align_res['read_end']-align_res['read_begin']+1
                        still_soft_len=len(left_soft_seq)-align_res['read_end']
                        clip_match_start = right_ext_region[1]+align_res['ref_begin']
                        clip_match_end = right_ext_region[1]+align_res['ref_end']
                        inter_gap_len = align_res['read_begin'] -1
                        new_match_start = real_inv_region[1] + (real_inv_region[2]-r.reference_end)
                        new_match_end = real_inv_region[2] - (r.reference_start+1-real_inv_region[1])
                        N_len = clip_match_start- new_match_end
                        inverted_match = r.reference_end-r.reference_start
                        if r.query_name == "A00159:976:H35MKDSX3:1:2556:32669:4304":
                            print("real_inv_region:", real_inv_region)
                            print("r.reference_start:", r.reference_start)
                            print("r.reference_end:", r.reference_end)
                            print("new_match_start:", new_match_start)
                            print("new_match_end:", new_match_end)
                            print("N_len:", N_len)

                        r.reference_start = new_match_start-1
                        cigarstring = r.cigarstring
                        cigarstring = mirror_cigar(cigarstring)
                        cigarstring = remove_soft_clipping(cigarstring, remove_left=False, remove_right=True)
                        cigarstring = gen_cigar_left_clip(cigarstring, still_soft_len, left_tag, N_len, inter_gap_len, match_len)
                        r.cigarstring = cigarstring
                        ## TODO:: set XS tag by neibour exon
                        r.set_tag('XS', '+', value_type='A')
                        # only find two exon for extreme case
                        if still_soft_len > clip_threshold:
                            print("still_soft_len:", still_soft_len)
                            continue
                        print("after:",r)
                elif right_soft_clip > clip_threshold and left_soft_clip <=clip_threshold:
                    # align soft sequence to left exon sequence
                    # get right soft sequence   

                    right_soft_seq = get_query_seq_for_inv_rg(r, right_soft_clip, "r")
                    align_res = best_alignment2(right_soft_seq, left_ref)

                    if align_res:
                        match_len=align_res['read_end']-align_res['read_begin']+1
                        still_soft_len=align_res['read_begin']-1
                        clip_match_start = left_ext_region[1]+align_res['ref_begin']
                        clip_match_end = clip_match_start+match_len-1
                        inter_gap_len = len(right_soft_seq)-align_res['read_end']
                        new_match_start = real_inv_region[1] + (real_inv_region[2]-r.reference_end)
                        new_match_end = real_inv_region[2] - (r.reference_start+1-real_inv_region[1])
                        N_len = new_match_start-clip_match_end
                        inverted_match = r.reference_end-r.reference_start

                        r.reference_start = clip_match_start-1
                        cigarstring = r.cigarstring
                        cigarstring = mirror_cigar(cigarstring)
                        cigarstring = remove_soft_clipping(cigarstring, remove_left=True, remove_right=False)
                        cigarstring = gen_cigar_right_clip(cigarstring, still_soft_len, right_tag, N_len, inter_gap_len, match_len)
                        r.cigarstring = cigarstring
                        r.set_tag('XS', '+', value_type='A')
                        # only find two exon for extreme case
                        print("after:",r)

                        if still_soft_len > clip_threshold:
                            print("still_soft_len:", still_soft_len)
                            continue
                elif right_soft_clip > clip_threshold and left_soft_clip >clip_threshold:
                    # align soft sequence to right exon sequence
                    # get left soft sequence
                    left_soft_seq = get_query_seq_for_inv_rg(r, left_soft_clip, "l")
                    # get right exon sequence
                    align_res = best_alignment2(left_soft_seq, right_ref)
                    if align_res:
                        match_len=align_res['read_end']-align_res['read_begin']+1
                        still_soft_len=len(left_soft_seq)-align_res['read_end']
                        clip_match_start = right_ext_region[1]+align_res['ref_begin']
                        clip_match_end = right_ext_region[1]+align_res['ref_end']
                        inter_gap_len = align_res['read_begin'] -1
                        new_match_start = real_inv_region[1] + (real_inv_region[2]-r.reference_end)
                        new_match_end = real_inv_region[2] - (r.reference_start+1-real_inv_region[1])
                        N_len = clip_match_start- new_match_end
                        inverted_match = r.reference_end-r.reference_start

                        # r.reference_start = new_match_start
                        cigarstring = r.cigarstring
                        cigarstring = mirror_cigar(cigarstring)
                        cigarstring = remove_soft_clipping(cigarstring, remove_left=True, remove_right=True)
                        cigarstring = gen_cigar_left_clip(cigarstring, still_soft_len, left_tag, N_len, inter_gap_len, match_len)
                        # r.cigarstring = cigarstring
                        ## TODO:: set XS tag by neibour exon
                        # r.set_tag('XS', '+', value_type='A')
                        # only find two exon for extreme case
                        if still_soft_len > clip_threshold:
                            print("still_soft_len:", still_soft_len)
                            continue
                    # get right soft sequence   
                    right_soft_seq = get_query_seq_for_inv_rg(r, right_soft_clip, "r")
                    align_res = best_alignment2(right_soft_seq, left_ref)
                    if align_res:
                        match_len=align_res['read_end']-align_res['read_begin']+1
                        still_soft_len=align_res['read_begin']-1
                        clip_match_start = left_ext_region[1]+align_res['ref_begin']
                        clip_match_end = clip_match_start+match_len-1
                        inter_gap_len = len(right_soft_seq)-align_res['read_end']
                        new_match_start = real_inv_region[1] + (real_inv_region[2]-r.reference_end)
                        new_match_end = real_inv_region[2] - (r.reference_start+1-real_inv_region[1])
                        N_len = new_match_start-clip_match_end
                        inverted_match = r.reference_end-r.reference_start

                        # r.reference_start = clip_match_start
                        # cigarstring = r.cigarstring
                        # cigarstring = mirror_cigar(cigarstring)
                        # cigarstring = remove_soft_clipping(cigarstring, remove_left=True, remove_right=False)
                        cigarstring = gen_cigar_right_clip(cigarstring, still_soft_len, right_tag, N_len, inter_gap_len, match_len)
                        r.reference_start = clip_match_start-1
                        r.cigarstring = cigarstring
                        r.set_tag('XS', '+', value_type='A')
                        # only find two exon for extreme case
                        if still_soft_len > clip_threshold:
                            print("still_soft_len:", still_soft_len)
                            continue
                print("##########################################################################################")
                print("after:", r)
                cigar_len = sum(c[1] for c in r.cigartuples if c[0] != 3)
                print(len(r.query_sequence), cigar_len, r.cigartuples)
                if len(r.query_sequence) != cigar_len:
                    print("cigar_err:", r)
                    reads_dict1[rname][idx] = r_cp
                    continue
                print("##########################################################################################")
                reads_dict1[rname][idx] = r
            # else:
            #     print("unprocessed :", r)
    return reads_dict1



def reformat_bam():
    # reads_dict2=process_hisat2(reads_dict2)
    # # reads_dict1 = {}

    # reads_dict1=process_bwa(reads_dict1)
    bam_f = pysam.AlignmentFile(unmap_remap_bam, "rb")
    hisat_bam_f = pysam.AlignmentFile(args.hisat_bam, "rb")
    hisat_reads_dict= {}

    for idx, real_inv in enumerate(real_inv_regions):
        unmap_remap_dict = get_reads_dict(bam_f, real_inv)
        # hisat_reads_dict = get_reads_dict(hisat_bam_f, real_inv)

        unmap_remap_dict = process_bwa2(unmap_remap_dict, real_inv, gene_regions[idx], real_inv_regions[idx])
        # hisat_reads_dict = process_bwa2(hisat_reads_dict, real_inv, gene_regions[idx], real_inv_regions[idx])

        # print(len(hisat_reads_dict), len(unmap_remap_dict))
    hisat_bam_f.close()
    hisat_bam_f = pysam.AlignmentFile(args.hisat_bam, "rb")
    bam_f = pysam.AlignmentFile(unmap_remap_bam, "rb")
    outf = pysam.AlignmentFile(out_bam, "wb", template=hisat_bam_f)

    write_dict_to_bam(outf, unmap_remap_dict, hisat_reads_dict, hisat_bam_f, bam_f)
    



def write_dict_to_bam(outf, unmap_remap_dict,hisat_reads_dict, hisat_bam_f, bam_f):
    for rnames, rec in unmap_remap_dict.items():
        for r in rec:
            if not r.is_supplementary:
                outf.write(r)
    for rec in bam_f:
        if rec.query_name not in unmap_remap_dict:
            outf.write(rec)
    # for rnames, rec in hisat_reads_dict.items():
    #     for r in rec:
    #         outf.write(r)
    # cnt = 0
    # write_cnt = 0
    # for rec in hisat_bam_f:
    #     cnt+=1
    #     if rec.query_name not in hisat_reads_dict:
    #         write_cnt+=1
    #         outf.write(rec)
    outf.close()
    # print("cnt:", cnt)
    # print("write_cnt:", write_cnt)

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

def reset_real_inv_regions(selected_regions, exons_ori):
    ori_exons_f = open(exons_ori, "r")
    ori_exons = []
    real_inv_regions = set()
    for line in ori_exons_f:
        line = line.strip().split('\t')
        chrom, start, end = line[0], int(line[1]), int(line[2])
        ori_exons.append((chrom, start, end))
    for selected_region in selected_regions:
        for ori_exon in ori_exons:
            intersection = get_intersection((selected_region[1], selected_region[2]), (ori_exon[1],ori_exon[2]))
            if intersection:
                real_inv_regions.add((selected_region[0], intersection[0], intersection[1]))
    return list(real_inv_regions)



if __name__ == "__main__":

    parser = ap.ArgumentParser()
    parser.add_argument('-b1', '--unmap_remap_bam', )
    parser.add_argument('-b2', '--hisat_bam', )
    parser.add_argument('-inv', '--inv_region', )
    parser.add_argument('-ref', '--ref', )
    parser.add_argument('-outbam', '--outbam', )
    parser.add_argument('-eo', '--exons_ori', )
    args = parser.parse_args()
    print(args)

    unmap_remap_bam = args.unmap_remap_bam
    clip_threshold = 5
    

    ref = args.ref
    out_bam = args.outbam
    selected_regions, gene_regions, real_inv_regions = parse_regions(args.inv_region)
    print("selected_regions:", selected_regions)
    print("gene_regions:", gene_regions)
    print("real_inv_regions:", real_inv_regions)
    # real_inv_regions = reset_real_inv_regions(selected_regions, args.exons_ori)
    # print("selected_regions:", selected_regions)
    # print("gene_regions:", gene_regions)
    # print("real_inv_regions:", real_inv_regions)
    ref_db = SeqIO.to_dict(SeqIO.parse(ref, "fasta"))
    reformat_bam()
    

    