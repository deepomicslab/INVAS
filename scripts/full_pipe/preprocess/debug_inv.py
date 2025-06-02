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
    


from Bio import SeqIO

def fastq_to_dict(fastq_file):
    seq_dict = SeqIO.to_dict(SeqIO.parse(fastq_file, "fastq"))
    return seq_dict

import re




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




def best_alignment(query, ref, open_penalty=10, extend_penalty=1, matrix=parasail.blosum62):
    # 对query序列进行反向互补转换
    # query_revcomp = str(Seq(query).reverse_complement())
    ref_c = get_complement(ref)

    # 对原始query序列进行比对
    result = parasail.ssw(query, ref, open_penalty, extend_penalty, matrix)

    # 对反向互补的query序列进行比对
    result_refc = parasail.ssw(query, ref_c, open_penalty, extend_penalty, matrix)

    print("result:", result.ref_begin1, result.ref_end1, result.score1)
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
    result = align(query, ref)

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
    print(chrom, start, end)
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





if __name__ == "__main__":

    parser = ap.ArgumentParser()
    parser.add_argument('-s1', '--seq_query', )
    parser.add_argument('-s2', '--seq_ref', )
    args = parser.parse_args()
    print(args)

    

    ref = args.seq_ref
    # out_bam = args.outbam
    # selected_regions, gene_regions, real_inv_regions = parse_regions(args.inv_region)
    query_db = SeqIO.to_dict(SeqIO.parse(args.seq_query, "fasta"))
    query_seq=str(query_db[list(query_db.keys())[0]].seq)
    ref_db = SeqIO.to_dict(SeqIO.parse(ref, "fasta"))
    ref_seq=str(ref_db[list(ref_db.keys())[0]].seq)
    print(query_seq)
    print(type(ref_seq))
    best_alignment2(query_seq, ref_seq)
    

    