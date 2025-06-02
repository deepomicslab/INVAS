import argparse

def reverse_complement(dna_sequence):
    # 定义互补关系
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # 生成互补序列
    complement_sequence = ''.join([complement[base] for base in dna_sequence])
    
    # 返回反向序列
    return complement_sequence[::-1]


if __name__ == "__main__":
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Track exon chains in sequencing reads')
    parser.add_argument('seq', help='Path to the BAM file')

    args = parser.parse_args()
    rc_seq= reverse_complement(args.seq)
    print(rc_seq)