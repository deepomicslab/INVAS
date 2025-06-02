import sys
import pysam
import re
from Bio.Seq import Seq
import parasail
import bisect 
from pybedtools import BedTool
from pyssw import align
import argparse as ap

def get_reads_dict(bam_file):
    with pysam.AlignmentFile(bam_file, "rb") as file:
        reads_dict = {}
        for read in file:
            if read.query_name in reads_dict:
                reads_dict[read.query_name].append(read)
            else:
                reads_dict[read.query_name] = [read]
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
            print("unique read:", reads_dict1[read_name][0])
            continue
        passed_cnt+=1
        new_splice_rec = copy_read(reads_dict2[read_name][0])
        reads_dict1[read_name].sort(key=lambda a: a.reference_start)
        merged_end = 0
        merged_cigar_str= ""
        print("before transform:##################################################################33")
        for idx, alignment in enumerate(reads_dict1[read_name]):
            print("bwa:", alignment)
        print("hisat:", reads_dict2[read_name][0])
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
        print("new_rec:",new_splice_rec)
        print("old_rec:",reads_dict2[read_name][0])
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
    print("truth seq", read.query_sequence, read.query_sequence[:soft_clip_len])
    if direction == "l":
        return read.query_sequence[:soft_clip_len]
    elif direction == "r":
        return read.query_sequence[-soft_clip_len:]
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



def best_alignment(query, ref, open_penalty, extend_penalty, matrix):
    # 对query序列进行反向互补转换
    # query_revcomp = str(Seq(query).reverse_complement())
    ref_revc = get_reverse_complement(ref)

    # 对原始query序列进行比对
    result = parasail.ssw(query, ref, open_penalty, extend_penalty, matrix)

    # 对反向互补的query序列进行比对
    result_revcomp = parasail.ssw(query, ref_revc, open_penalty, extend_penalty, matrix)

    # 比较两个方向的比对分数，并选择分数最高的那个
    if result.score1 >= result_revcomp.score1:
        return result
    else:
        return result_revcomp

def best_alignment2(query, ref):
    # 对query序列进行反向互补转换
    # query_revcomp = str(Seq(query).reverse_complement())
    # ref_revc = get_reverse_complement(ref)

    # 对原始query序列进行比对
    result = align(query, ref)

    # 对反向互补的query序列进行比对
    result_revcomp = align(query, ref, True)

    # 比较两个方向的比对分数，并选择分数最高的那个
    if result_revcomp['score'] >= result['score']:
        return result_revcomp
    else:
        return result



def process_hisat2(reads_dict2):
    for rname, rec in reads_dict2.items():
        for idx, r in enumerate(rec):
            if 'S' in r.cigarstring:
                # if r.query_name != "read118/ENST00000399688;mate1:222-321;mate2:362-460":
                #     continue
                print("##########################################################################################")
                print("before:", r)
                print("##########################################################################################")
                left_soft_clip, right_soft_clip = parse_cigar(r.cigarstring)
                print(r.cigarstring, f'Left soft clip: {left_soft_clip}', f'Right soft clip: {right_soft_clip}' )
                if left_soft_clip > 2:
                    # align soft sequence to left exon sequence
                    # get left soft sequence
                    left_soft_seq = get_query_seq(r, left_soft_clip, "l")
                    # get left exon sequence
                    left_exon_end, left_exon_seq = get_ref_seq(r, "l")
                    if left_exon_seq == "":
                        continue
                    print("left_soft_seq:", left_soft_seq)
                    print("left_exon_seq:", left_exon_seq)
                    align_res = best_alignment(left_soft_seq, left_exon_seq, 10, 1, parasail.blosum62)
                
                    if align_res:
                        print("align_ref:", align_res.ref_begin1, align_res.ref_end1)
                        print("align_read:", align_res.read_begin1, align_res.read_end1)
                        print("left_exon_end:", left_exon_end)
                        match_len=align_res.read_end1-align_res.read_begin1+1
                        still_soft_len=len(left_soft_seq)-align_res.read_end1-1
                        print("still_soft_len:", still_soft_len)
                        new_match_start = left_exon_end-align_res.ref_end1
                        left_exon_match_end = left_exon_end-align_res.ref_begin1
                        gap_len=r.reference_start-left_exon_match_end
                        r.reference_start = new_match_start
                        cigarstring = r.cigarstring
                        cigarstring = remove_soft_clipping(cigarstring, remove_left=True, remove_right=False)
                        r.cigarstring = f"{still_soft_len}S{match_len}M{gap_len}N"+cigarstring if still_soft_len>0 \
                              else f"{match_len}M{gap_len}N"+cigarstring
                        # only find two exon for extreme case
                        r.flag &= ~0x8
                        r.set_tag('XS', '+', value_type='A')
                        if still_soft_len>2:
                            left_soft_seq = get_query_seq(r, still_soft_len, "l")
                            # get left exon sequence
                            left_exon_end, left_exon_seq = get_ref_seq(r, "l")
                            if left_exon_seq == "":
                                continue
                            align_res = best_alignment(left_soft_seq, left_exon_seq, 10, 1, parasail.blosum62)
                            if align_res:
                                match_len=align_res.read_end1-align_res.read_begin1
                                still_soft_len=left_soft_clip-match_len
                                new_match_start = left_exon_end-align_res.ref_end1
                                cigarstring = r.cigarstring
                                cigarstring = remove_soft_clipping(cigarstring, remove_left=True, remove_right=False)
                                gap_len = r.reference_start-new_match_start
                                r.reference_start = new_match_start
                                r.cigarstring = f"{left_soft_clip}S{match_len}M{gap_len}N"+cigarstring if still_soft_len>0 \
                                    else f"{match_len}M{gap_len}N"+cigarstring
                if right_soft_clip > 2:
                    # align soft sequence to right exon sequence
                    # get right soft sequence
                    right_soft_seq = get_query_seq(r, right_soft_clip, "r")
                    # get right exon sequence
                    right_exon_start, right_exon_seq = get_ref_seq(r, "r")
                    if right_exon_seq == "":
                        continue
                    print("right_soft_seq:", right_soft_seq)
                    print("right_exon_seq:", right_exon_seq)
                    align_res = best_alignment(right_soft_seq, right_exon_seq, 10, 1, parasail.blosum62)
                    print("\nalign_res:", align_res)
                    print("align length:", align_res.read_end1-align_res.read_begin1)

                    if align_res:
                        match_len=align_res.read_end1-align_res.read_begin1+1
                        still_soft_len=len(right_soft_seq)-align_res.read_end1-1
                        new_match_start = right_exon_start + align_res.ref_begin1
                        print("new start:", new_match_start)
                        # r.reference_start = new_match_start
                        gap_len = new_match_start-r.reference_end
                        cigarstring = r.cigarstring
                        print("cigarstring:", cigarstring)
                        cigarstring = remove_soft_clipping(cigarstring, remove_left=False, remove_right=True)
                        print("cigarstring:", cigarstring)
                        print("still_soft_len:", still_soft_len)
                        r.cigarstring = cigarstring+f"{gap_len}N{len(right_soft_seq)}M" if still_soft_len==0 \
                                else cigarstring+f"{gap_len}N{len(right_soft_seq)}M{still_soft_len}S"
                        r.flag &= ~0x8
                        r.set_tag('XS', '+', value_type='A')
                        # only find two exon for extreme case
                        if still_soft_len>2:
                            right_soft_seq = get_query_seq(r, still_soft_len, "r")
                            # get right exon sequence
                            right_exon_start, right_exon_seq = get_ref_seq(r, "r")
                            if right_exon_seq == "":
                                continue
                            align_res = best_alignment(right_soft_seq, right_exon_seq, 10, 1, parasail.blosum62)
                            if align_res:
                                match_len=align_res.read_end1-align_res.read_begin1+1
                                still_soft_len=len(right_soft_seq)-align_res.read_end1-1
                                new_match_start = right_exon_start + align_res.ref_begin1
                                # r.reference_start = new_match_start
                                gap_len = new_match_start-r.reference_end
                                r.cigarstring = r.cigarstring+f"{gap_len}N{len(right_soft_seq)}M" if still_soft_len==0 \
                                else r.cigarstring+f"{gap_len}N{len(right_soft_seq)}M{still_soft_len}S"
                
                print("##########################################################################################")
                print("after:", r)
                
                print("##########################################################################################")
                reads_dict2[rname][idx] = r
    return reads_dict2

def process_bwa(reads_dict1):
    for rname, rec in reads_dict1.items():
        for idx, r in enumerate(rec):
            # if r.query_name != "read69/ENST00000399687;mate1:16-115;mate2:143-241_sg":
            #     continue
            if r.has_tag('XS'):
                continue
            if 'S' in r.cigarstring and not r.is_supplementary:
                r_cp = copy_read(r)
                print("##########################################################################################")
                print("before:", r)
                print("##########################################################################################")
                left_soft_clip, left_tag, right_soft_clip, right_tag = parse_cigar2(r.cigarstring)
                print(r.cigarstring, f'Left soft clip: {left_soft_clip}', f'Right soft clip: {right_soft_clip}' )
                if left_soft_clip > 2 and right_soft_clip <=2:
                    # align soft sequence to right exon sequence
                    # get left soft sequence
                    left_soft_seq = get_query_seq_for_inv_rg(r, left_soft_clip, "l")
                    # get right exon sequence
                    right_exon_start,right_exon_end, right_exon_seq, cur_header = get_ref_seq_for_inv_rg(r, "l")
                    if right_exon_seq == "" or cur_header=="":
                        continue
                    print("left_soft_seq:", left_soft_seq)
                    print("right_exon_seq:", right_exon_seq)
                    align_res = best_alignment2(left_soft_seq, right_exon_seq)
                    print("test_align res:",align_res['ref_begin'], align_res['ref_end'])
                    if align_res:
                        match_len=align_res['read_end']-align_res['read_begin']+1
                        still_soft_len=len(left_soft_seq)-align_res['read_end']
                        cur_exon_chr, cur_exon_start, cur_exon_end = cur_header.split(' ')[0].split('-')
                        cur_exon_start = int(cur_exon_start)
                        cur_exon_end = int(cur_exon_end)
                        new_match_start = cur_exon_end-(r.reference_end-cur_exon_start)
                        cur_exon_match = r.reference_end-r.reference_start
                        right_exon_match_pos = right_exon_start + align_res['ref_begin']
                        gap_len=right_exon_match_pos-new_match_start-cur_exon_match
                        r.reference_start = new_match_start
                        insertion_len = len(left_soft_seq)-still_soft_len-match_len
                        print("before query seq:", r.query_sequence, len(r.query_sequence)-insertion_len)
                        r.query_sequence = r.query_sequence[:len(r.query_sequence)-insertion_len]
                        print("after query seq:", r.query_sequence)
                        # cigarstring = r.cigarstring
                        # cigarstring = remove_soft_clipping(cigarstring, remove_left=True, remove_right=False)
                        right_exon_match_len=align_res['ref_end']-align_res['ref_begin']+1
                        if right_soft_clip == 0:
                            cigarstring = f"{cur_exon_match}M{gap_len}N{match_len}M" if still_soft_len ==0 \
                                else f"{cur_exon_match}M{gap_len}N{match_len}M{still_soft_len}{left_tag}"
                        else:
                            cigarstring = r.cigarstring
                            cigarstring = mirror_cigar(cigarstring)
                            cigarstring = remove_soft_clipping(cigarstring, remove_left=False, remove_right=True)
                            cigarstring = f"{cigarstring}{gap_len}N{match_len}M" if still_soft_len ==0 \
                                else f"{cigarstring}{gap_len}N{match_len}M{still_soft_len}{left_tag}"
                        r.cigarstring = cigarstring
                        r.set_tag('XS', '+', value_type='A')
                        # only find two exon for extreme case
                        if still_soft_len > 2:
                            print("still_soft_len:", still_soft_len)
                            continue
                elif right_soft_clip > 2 and left_soft_clip <=2:
                    # align soft sequence to left exon sequence
                    # get right soft sequence   
                    right_soft_seq = get_query_seq_for_inv_rg(r, right_soft_clip, "r")
                    # get left exon sequence
                    left_exon_start, left_exon_end, left_exon_seq, cur_header = get_ref_seq_for_inv_rg(r, "r")
                    if left_exon_seq == "" or cur_header=="":
                        continue
                    print("right_soft_seq:", right_soft_seq)
                    print("left_exon_seq:", left_exon_seq)
                    align_res = best_alignment2(right_soft_seq, left_exon_seq)
                    if align_res:
                        match_len=align_res['read_end']-align_res['read_begin']+1
                        still_soft_len=len(right_soft_seq)-align_res['read_end']
                        cur_exon_chr, cur_exon_start, cur_exon_end = cur_header.split(' ')[0].split('-')
                        cur_exon_start = int(cur_exon_start)
                        cur_exon_end = int(cur_exon_end)
                        new_match_start = cur_exon_start+(cur_exon_end-r.reference_end)
                        cur_exon_match = r.reference_end-r.reference_start
                        left_exon_match_pos = left_exon_start + align_res['ref_end']
                        gap_len=new_match_start-left_exon_match_pos
                        r.reference_start = left_exon_start+align_res['ref_begin']
                        # left_exon_match_pos = left_exon_end-align_res['ref_begin']
                        # gap_len=new_match_start-left_exon_match_pos
                        # r.reference_start = left_exon_end-align_res['ref_end']
                        insertion_len = len(right_soft_seq)-still_soft_len-match_len
                        print("before query seq:", r.query_sequence, len(r.query_sequence)-insertion_len)
                        r.query_sequence = r.query_sequence[:len(r.query_sequence)-insertion_len]
                        print("after query seq:", r.query_sequence)
                        # cigarstring = r.cigarstring
                        # cigarstring = remove_soft_clipping(cigarstring, remove_left=True, remove_right=False)
                        left_exon_match_len=align_res['ref_end']-align_res['ref_begin']+1
                        if left_soft_clip == 0:
                            cigarstring = f"{match_len}M{gap_len}N{cur_exon_match}M" if still_soft_len ==0 \
                                else f"{still_soft_len}{right_tag}{match_len}M{gap_len}N{cur_exon_match}M"
                        else:
                            cigarstring = r.cigarstring
                            cigarstring = mirror_cigar(cigarstring)
                            cigarstring = remove_soft_clipping(cigarstring, remove_left=True, remove_right=False)
                            cigarstring = f"{match_len}M{gap_len}N{cigarstring}" if still_soft_len ==0 \
                                else f"{still_soft_len}{right_tag}{match_len}M{gap_len}N{cigarstring}"
                        r.cigarstring = cigarstring
                        r.set_tag('XS', '+', value_type='A')
                        # only find two exon for extreme case
                        if still_soft_len > 2:
                            print("still_soft_len:", still_soft_len)
                            continue
                elif right_soft_clip > 2 and left_soft_clip >2:
                    # align soft sequence to right exon sequence
                    # get left soft sequence
                    left_soft_seq = get_query_seq_for_inv_rg(r, left_soft_clip, "l")
                    # get right exon sequence
                    right_exon_start,right_exon_end, right_exon_seq, cur_header = get_ref_seq_for_inv_rg(r, "l")
                    if right_exon_seq == "" or cur_header=="":
                        continue
                    print("left_soft_seq:", left_soft_seq)
                    print("right_exon_seq:", right_exon_seq)
                    align_res = best_alignment2(left_soft_seq, right_exon_seq)
                    print("6S:",align_res)
                    if align_res:
                        match_len=align_res['read_end']-align_res['read_begin']+1
                        still_soft_len=len(left_soft_seq)-align_res['read_end']
                        cur_exon_chr, cur_exon_start, cur_exon_end = cur_header.split(' ')[0].split('-')
                        cur_exon_start = int(cur_exon_start)
                        cur_exon_end = int(cur_exon_end)
                        cur_exon_real_end_pos = cur_exon_end-(r.reference_start-cur_exon_start)
                        right_exon_match_pos = right_exon_start + align_res['ref_begin']
                        gap_len=right_exon_match_pos-cur_exon_real_end_pos


                        # new_match_start = cur_exon_end-(cur_exon_start-r.reference_end)
                        # cur_exon_match = r.reference_end-r.reference_start
                        # right_exon_match_pos = right_exon_start + align_res.ref_begin1
                        # gap_len=right_exon_match_pos-new_match_start-cur_exon_match
                        # r.reference_start = new_match_start
                        # cigarstring = r.cigarstring
                        # cigarstring = remove_soft_clipping(cigarstring, remove_left=True, remove_right=False)
                        # right_exon_match_len=align_res.ref_end1-align_res.ref_begin1+1
                        
                        cigarstring = r.cigarstring
                        cigarstring = mirror_cigar(cigarstring)
                        cigarstring = remove_soft_clipping(cigarstring, remove_left=True, remove_right=True)
                        cigarstring = f"{cigarstring}{gap_len}N{match_len}M" if still_soft_len ==0 \
                            else f"{cigarstring}{gap_len}N{match_len}M{still_soft_len}{left_tag}"
                        r.set_tag('XS', '+', value_type='A')
                        # r.cigarstring = cigarstring
                        # only find two exon for extreme case
                        if still_soft_len > 2:
                            print("still_soft_len:", still_soft_len)
                            continue
                    # align soft sequence to left exon sequence
                    # get right soft sequence   
                    right_soft_seq = get_query_seq_for_inv_rg(r, right_soft_clip, "r")
                    # get left exon sequence
                    left_exon_start, left_exon_end, left_exon_seq, cur_header = get_ref_seq_for_inv_rg(r, "r")
                    if left_exon_seq == "" or cur_header=="":
                        continue
                    print("right_soft_seq2:", right_soft_seq)
                    print("left_exon_seq2:", left_exon_seq)
                    align_res = best_alignment2(right_soft_seq, left_exon_seq)
                    if align_res:
                        match_len=align_res['read_end']-align_res['read_begin']+1
                        still_soft_len=len(right_soft_seq)-align_res['read_end']
                        cur_exon_chr, cur_exon_start, cur_exon_end = cur_header.split(' ')[0].split('-')
                        cur_exon_start = int(cur_exon_start)
                        cur_exon_end = int(cur_exon_end)
                        read_match_end_in_cur_exon = get_last_match_position_in_region(r.cigarstring, r.reference_start, cur_exon_start, cur_exon_end)
                        read_start_in_cur_exon = cur_exon_start + (cur_exon_end-read_match_end_in_cur_exon)
                        left_exon_match_pos = left_exon_start + align_res['ref_end']
                        gap_len=read_start_in_cur_exon-left_exon_match_pos
                        r.reference_start = left_exon_start+align_res['ref_begin']
                        cigarstring = f"{match_len}M{gap_len}N{cigarstring}" if still_soft_len ==0 \
                            else f"{still_soft_len}{right_tag}{match_len}M{gap_len}N{cigarstring}"
                        r.cigarstring = cigarstring
                        r.set_tag('XS', '+', value_type='A')
                        # only find two exon for extreme case
                        if still_soft_len > 2:
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
    return reads_dict1

def reformat_bam(reads_dict1, reads_dict2):
    reads_dict2=process_hisat2(reads_dict2)
    # reads_dict1 = {}

    reads_dict1=process_bwa(reads_dict1)
    return reads_dict1, reads_dict2
    
def write_dict_to_bam(outf, reads_dict1, reads_dict2):
    for rnames, rec in reads_dict1.items():
        for r in rec:
            if not r.is_supplementary:
                outf.write(r)
    for rnames, rec in reads_dict2.items():
        for r in rec:
            outf.write(r)
    outf.close()



def get_juncs(reads_dict1, reads_dict2):
    juncs = {}
    for reads_dict in [reads_dict1, reads_dict2]:
        for rname, rec in reads_dict.items():
            for r in rec:
                if 'N' in r.cigarstring:
                    # print(r)
                    # Get the direction of the read from the XS tag
                    direction = r.get_tag('XS') if r.has_tag('XS') else None
                    # Find all junctions in the read
                    pos = r.reference_start
                    cigar_tuples = r.cigartuples
                    for i, (cigar_type, length) in enumerate(cigar_tuples):
                        if cigar_type == 3:  # N in CIGAR is represented by 3
                            # Calculate the start and end of the junction
                            junc_start = pos
                            junc_end = pos + length
                            # Add the junction to the dict
                            junction = (junc_start, junc_end)
                            # print(junction)

                            if junction not in juncs:
                                juncs[junction] = {'+': 0, '-': 0}
                            if direction is not None:
                                juncs[junction][direction] += 1
                        if cigar_type in {0, 2, 3, 7, 8}:  # M, D, N, =, X in CIGAR
                            pos += length
    return juncs
    



if __name__ == "__main__":

    parser = ap.ArgumentParser()
    parser.add_argument('-b1', '--bam_file1', )
    parser.add_argument('-b2', '--bam_file2', )
    parser.add_argument('-b3', '--bam_file3', )
    parser.add_argument('-fq1', '--fq1', )
    parser.add_argument('-fq2', '--fq2', )
    parser.add_argument('-exon', '--exon_db', )
    parser.add_argument('-outbam', '--outbam', )
    parser.add_argument('-unmap', '--unmap_fq', )

    # parser.add_argument('target', help='target file')
    # parser.add_argument('query', help='query file')
    # if len(sys.argv) == 1:
    #     parser.print_help()
    #     sys.exit()
    args = parser.parse_args()
    print(args)
    # if len(sys.argv) != 9:
    #     print(f"Usage: {sys.argv[0]} <bam_file1> <bam_file2> <bam_file3> <fq1> <fq2> <exon_db> <outbam>")
    #     sys.exit(1)
    max_intron_size = 800000
    bam_file1 = args.bam_file1
    bam_file2 = args.bam_file2
    bam_file3 = args.bam_file3
    fq1 = args.fq1
    fq2 = args.fq2
    exon_fa = args.exon_db
    out_bam = args.outbam
    # inv_rg_file = sys.argv[8]
    unmap_fq = args.unmap_fq

    # bed_inv = BedTool(inv_rg_file)

    

    reads_dict1 = get_reads_dict(bam_file1)
    reads_dict2 = get_reads_dict(bam_file2)
    reads_dict3 = get_reads_dict(bam_file3)
    outf=pysam.AlignmentFile(out_bam, 'wb', header=reads_dict2[list(reads_dict2.keys())[0]][0].header)
    r1_dict = fastq_to_dict(fq1)
    r2_dict = fastq_to_dict(fq2)
    umpr_dict = fastq_to_dict(unmap_fq)

    exon_db = SeqIO.to_dict(SeqIO.parse(exon_fa, "fasta"))
    exon_headers = [exon_db[header].description for header in exon_db.keys()]

    print("Reads in file 1:", len(reads_dict1))
    print("Reads in file 2:", len(reads_dict2))
    # print_common_reads(reads_dict1, reads_dict2)
    # print("debug splice rec #####################################")
    # for rname, rec in reads_dict2.items():
    #     for r in rec:
    #         if 'S' in r.cigarstring:
    #             print(r)
        # outf.write(rec[0])
    reads_dict1, reads_dict2 =  reformat_bam(reads_dict1, reads_dict2)
    juncs = get_juncs(reads_dict1, reads_dict2)
    write_dict_to_bam(outf, reads_dict1, reads_dict2)
    # write reads dict1 and dict2 to bam
    



    # print("debug normal splice rec #####################################")
    # for rname, rec in reads_dict3.items():
    #     for r in rec:
    #         if 'S' in r.cigarstring:
    #             print(r)
        
