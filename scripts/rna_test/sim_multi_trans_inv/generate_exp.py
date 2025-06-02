import argparse
import re
import pysam

# Function to parse transcript exons file
# Function to parse transcript exons file and merge transcripts with the same start and end
def parse_trans_path(trans_paths):
    paths=trans_paths.split('\t')
    path_with_strand={}
    
    for path in paths:
        if path in ["s3+", "t5+"] or\
            path.endswith("3+") or path.endswith("3-"):
            continue
        seg_id = path.split(" ")[0]
        strand = path[-1]
        path_with_strand[int(seg_id)] = strand
    return path_with_strand
    

def parse_transcript_file(filename):
    transcripts = {}
    current_transcript_id = None
    
    with open(filename, 'r') as file:
        trans_id = 1
        path_to_strand={}
        trans_dict = {}
        trans_paths=  ""
        for line in file:
            line = line.strip()
            if line.startswith('preTrans'):
                trans_paths=line.split(':')[1]
                paths_with_strand=parse_trans_path(trans_paths)
                if trans_paths not in path_to_strand:
                    path_to_strand[trans_paths]=paths_with_strand
            else:
                exon_id, start, end, cov = line.strip().split('\t')
                cov=float(cov)
                if trans_paths not in trans_dict:
                    trans_dict[trans_paths]={}
                    trans_dict[trans_paths][exon_id]=(start, end, cov)
                elif exon_id not in trans_dict[trans_paths]:
                    trans_dict[trans_paths][exon_id]=(start, end, cov)
                else:
                    start, end, cov_old = trans_dict[trans_paths][exon_id]
                    newcov= cov_old + cov
                    trans_dict[trans_paths][exon_id]=(start, end, newcov)
    return trans_dict, path_to_strand

                

# Example usage:

# Now transcripts is a dictionary where each key is a transcript ID and the value
# is another dictionary with (start, end) tuples as keys and the sum of coverage as values.

# Function to parse BAM status file
def parse_bam_status_file(filename):
    with open(filename, 'r') as file:
        headers = next(file).strip().split('\t')  # Skip the header line
        values = next(file).strip().split('\t')
        return dict(zip(headers, map(float, values)))
    

# def process_trans():
#     sum_cov_recover=0
#     trans_dict = {}
#     for trans_id, exons in transcripts.items():
#         trans_sum_cov=0
#         trans_sum_len=0
#         for exon in exons:
#             trans_sum_cov+=(exon['cov']*(exon['end']-exon['start']))
#             trans_sum_len+=(exon['end']-exon['start'])
#         t_cov = trans_sum_cov/trans_sum_len
#         sum_cov_recover+=t_cov
#         trans_dict[trans_id]= t_cov
#     print(trans_dict)
#     for trans_id, trans_cov in trans_dict.items():

#         fpkm = calculate_fpkm(trans_cov, bam_status['frag_len']+frag_len_recover)
#         tpm = calculate_tpm(trans_cov, bam_status['sum_cov']+sum_cov_recover)
#         print(f"trans [{trans_id}], t_cov: {trans_cov}, fpkm: {fpkm}, tmp: {tpm}\n")
        
def calculate_fpkm(tcov, Frag_Len):
    if Frag_Len>0.001:
        return tcov*1000000000/Frag_Len
    else:
        return 0.0

def calculate_tpm(tcov, Cov_Sum):
    if Cov_Sum>0.00001:
        return tcov*1000000/Cov_Sum
    else:
        return 0.01

def calculate_read_length_times_count(bam_file):
    # 打开BAM文件
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        total = 0
        for read in bam.fetch(until_eof=True):
            # 只计算已映射的read
            if not read.is_unmapped:
                read_length = read.query_length  # 获取read的长度
                total += read_length  # 将该read的长度加到总数中

    # 由于每个read计数为1，所以总数已经是read_length * read_count
    return total

# Updated function to process transcripts with merged exons
def process_trans():
    sum_cov_recover = 0
    trans_dict = {}
    trans_range_dict = {}
    
    for trans_id, exons_dict in transcripts_dict.items():
        trans_sum_cov = 0
        trans_sum_len = 0

        trans_start=10000000000
        trans_end = -10000000000
        
        for exon_id, (start, end, cov) in exons_dict.items():
            exon_length = int(end) - int(start)
            trans_start=min(trans_start, int(start))
            trans_end=max(trans_end, int(end))
            trans_sum_cov += cov * exon_length
            trans_sum_len += exon_length
        trans_range_dict[trans_id]=(trans_start, trans_end)
        
        t_cov = trans_sum_cov / trans_sum_len if trans_sum_len > 0 else 0
        sum_cov_recover += t_cov
        trans_dict[trans_id] = t_cov

    trans_exp_dict = {}
    print(bam_status['sum_cov'], sum_cov_recover)
    for trans_id, trans_cov in trans_dict.items():
        fpkm = calculate_fpkm(trans_cov, bam_status['frag_len'] + frag_len_recover)
        tpm = calculate_tpm(trans_cov, sum_cov_recover)
        trans_exp_dict[trans_id]={
            "fpkm":fpkm,
            "tpm": tpm
        }
        # print(f"trans [{trans_id}], t_cov: {trans_cov}, fpkm: {fpkm}, tpm: {tpm}, strand: {path_to_strand[trans_id]}\n")

    # to gtf file
    res_f= open(args.out_gtf, 'w')
    res_f.write("# SeqFlow version 0.0.0\n")
    trans_index = 1
    for trans_id, exons_dict in transcripts_dict.items():
        res_f.write(f"""{args.chrom}\tSeqflow\ttranscripts\t{trans_range_dict[trans_id][0]}\t{trans_range_dict[trans_id][1]}\t1000\t+\t.\tgene_id "SEQFLOW.1"; transcripts_id "SEQFLOW1.{trans_index}"; cov "{trans_dict[trans_id]}"; FPKM "{trans_exp_dict[trans_id]['fpkm']}" TPM "{trans_exp_dict[trans_id]['tpm']}"\n""")
        trans_index +=1
        all_seq_idx = sorted(list(path_to_strand[trans_id].keys()))
        exon_strands = [path_to_strand[trans_id][x] for x in all_seq_idx]
        exon_idx = 0
        for exon_id, (start, end, cov) in exons_dict.items():
            exon_id=int(exon_id)
            # print(path_to_strand[trans_id], exon_id, path_to_strand[trans_id].keys(), (0 in path_to_strand[trans_id].keys()))
            # print(path_to_strand[trans_id][int(exon_id)])
            res_f.write(f"""{args.chrom}\tSeqflow\texon\t{start}\t{end}\t1000\t{exon_strands[exon_idx]}\t.\tgene_id "SEQFLOW.1"; transcripts_id "SEQFLOW1.{trans_index}"; exon_number "{exon_id+1}" cov {cov}\n""")
            exon_idx += 1

    res_f.close()


        

def main():
    # Parse files
    process_trans()
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse transcript exons and BAM status files.")
    parser.add_argument("exon_file", help="File containing transcript exons information.")
    parser.add_argument("bam_status_file", help="File containing BAM status information.")
    parser.add_argument("recover_bam", help="File containing BAM status information.")
    parser.add_argument("out_gtf", help="File containing BAM status information.")
    parser.add_argument("chrom", help="File containing BAM status information.")


    args = parser.parse_args()
    bam_status = parse_bam_status_file(args.bam_status_file)
    frag_len_recover = calculate_read_length_times_count(args.recover_bam)


    transcripts_dict, path_to_strand = parse_transcript_file(args.exon_file)
    # for k, v in transcripts_dict.items():
    #     print(k, v)
    


    main()