import argparse
import re
import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

# read inversion.exons file
def read_inversion_exons(file):
    inversion_exons = []
    with open(file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            chrom, start, end, coverage = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
            inversion_exons.append((chrom, start, end, "-", coverage))  # 负链方向
    return inversion_exons

# read JUNC file
def read_junc_file(file):
    junc_connections = []
    with open(file, 'r') as f:
        for line in f:
            match = re.match(r"JUNC (\d+)-(\d+) (\+|-) (\d+)-(\d+) (\+|-) ([\d.]+)", line.strip())
            if match:
                start1, end1, strand1, start2, end2, strand2, coverage = match.groups()
                exon1 = (int(start1), int(end1), strand1)
                exon2 = (int(start2), int(end2), strand2)
                junc_connections.append((exon1, exon2, float(coverage)))
    return junc_connections

# read GTF file
def read_gtf_file(file):
    transcripts = defaultdict(list)
    with open(file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if parts[2] == "exon":
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes = parts[8]
                match = re.search(r'transcript_id "([^"]+)"', attributes)
                if match:
                    transcript_id = match.group(1)
                    transcripts[transcript_id].append((chrom, start, end, strand))
    return transcripts


def read_gtf_file2(file):
    """
    parses a GTF file and returns a dictionary of transcripts with their exons and coverage information.
    """
    transcripts = defaultdict(lambda: {
        "gene_id": "",
        "strand": "",
        "exons": [],
        "exon_coverage": [],
        "transcript_coverage": 0.0,
        "FPKM": 0.0,
        "TPM": 0.0,
    })

    with open(file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            # extract gene_id and transcript_id
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
            transcript_id_match = re.search(r'transcript_id "([^"]+)"', attributes)

            gene_id = gene_id_match.group(1) if gene_id_match else ""
            transcript_id = transcript_id_match.group(1) if transcript_id_match else ""

            # initialize transcript entry if not exists
            if transcript_id and transcript_id not in transcripts:
                transcripts[transcript_id]["gene_id"] = gene_id
                transcripts[transcript_id]["strand"] = strand

            # handle exon entries
            if feature_type == "exon":
                transcripts[transcript_id]["exons"].append((chrom, start, end, strand))

                # extract exon coverage if available
                exon_coverage_match = re.search(r'cov "([^"]+)"', attributes)
                print("exon_coverage_match:", exon_coverage_match)
                if exon_coverage_match:
                    exon_coverage = float(exon_coverage_match.group(1))
                    transcripts[transcript_id]["exon_coverage"].append((start, end, exon_coverage))
                else:
                    transcripts[transcript_id]["exon_coverage"].append((start, end, None))

            # handle transcript entries for coverage, FPKM, TPM
            if feature_type == "transcript":
                transcript_coverage_match = re.search(r'cov "([^"]+)"', attributes)
                FPKM_match = re.search(r'FPKM "([^"]+)"', attributes)
                TPM_match = re.search(r'TPM "([^"]+)"', attributes)

                if transcript_coverage_match:
                    transcripts[transcript_id]["transcript_coverage"] = float(transcript_coverage_match.group(1))
                if FPKM_match:
                    transcripts[transcript_id]["FPKM"] = float(FPKM_match.group(1))
                if TPM_match:
                    transcripts[transcript_id]["TPM"] = float(TPM_match.group(1))

    return transcripts

# check if two exons overlap
def exon_overlap(exon1, exon2):
    start1, end1 = exon1
    start2, end2 = exon2
    return not (end1 < start2 or start1 > end2)

# find connections for each inversion.exon
def find_connections_for_inversion(inversion_exons, junc_connections):
    inversion_connections = {}
    for chrom, inv_start, inv_end, inv_strand, inv_coverage in inversion_exons:
        for exon1, exon2, junc_coverage in junc_connections:
            # check if exon1 or exon2 overlaps with inversion.exon
            if exon_overlap((inv_start, inv_end), (exon1[0], exon1[1])):
                inversion_connections[(chrom, inv_start, inv_end)] = ("prev", exon1, exon2)
            elif exon_overlap((inv_start, inv_end), (exon2[0], exon2[1])):
                inversion_connections[(chrom, inv_start, inv_end)] = ("next", exon1, exon2)
    return inversion_connections

# split transcripts by connections
def split_transcripts_by_connections(inversion_connections, transcripts):
    new_transcripts = defaultdict(list)
    updated_transcripts = defaultdict(list)

    for (chrom, inv_start, inv_end), (relation, exon1, exon2) in inversion_connections.items():
        for transcript_id, exons in transcripts.items():
            # check if both exon1 and exon2 are in the transcript
            if exon1 in exons and exon2 in exons:
                # create a new transcript with the inversion.exon
                new_transcripts[transcript_id].append((chrom, inv_start, inv_end, "-", exons))
                # remove the inversion.exon from the original transcript
                updated_transcripts[transcript_id] = [
                    e for e in exons if not exon_overlap((e[1], e[2]), (inv_start, inv_end))
                ]
    return new_transcripts, updated_transcripts

# write transcripts to file
def write_transcripts(transcripts, output_file, is_new=False):
    with open(output_file, "w") as f:
        for transcript_id, exons in transcripts.items():
            for exon in exons:
                if is_new:
                    chrom, start, end, strand, coverage = exon
                    f.write(f"{transcript_id}\t{chrom}\t{start}\t{end}\t{strand}\t{coverage}\n")
                else:
                    chrom, start, end, strand = exon
                    f.write(f"{transcript_id}\t{chrom}\t{start}\t{end}\t{strand}\n")

def find_flanking_exons_single_pair(inversion_exons, junc_connections):
    """
    find the flanking exons for each inversion exon based on JUNC connections.
    """
    print("inversion_exons:", inversion_exons)
    print("junc_connections:", junc_connections)
    flanking_exons = {}

    # initialize each inversion exon entry
    for inversion_exon in inversion_exons:
        chrom, inv_start, inv_end, inv_strand, inv_coverage = inversion_exon
        flanking_exons[(chrom, inv_start, inv_end, inv_strand, inv_coverage)] = {
            "prev_exon": None,
            "next_exon": None
        }

    # loop through each junction connection
    for exon1, exon2, junc_coverage in junc_connections:
        exon1_start, exon1_end, exon1_strand = exon1
        exon2_start, exon2_end, exon2_strand = exon2

        # loop through each inversion exon to see if it matches exon1 or exon2
        for inversion_exon in inversion_exons:
            chrom, inv_start, inv_end, inv_strand, inv_coverage = inversion_exon
            key = (chrom, inv_start, inv_end, inv_strand, inv_coverage)

            # if current exon is exon2 (later), record prev_exon
            if exon2_start == inv_start and exon2_end == inv_end:
                if flanking_exons[key]["prev_exon"] is None:  
                    flanking_exons[key]["prev_exon"] = (exon1_start, exon1_end, exon1_strand, junc_coverage)

            # if current exon is exon1 (previous), record next_exon
            elif exon1_start == inv_start and exon1_end:
                if flanking_exons[key]["next_exon"] is None: 
                    flanking_exons[key]["next_exon"] = (exon2_start, exon2_end, exon2_strand, junc_coverage)


    return flanking_exons


def reverse_sign(sign):
    if sign == "+":
        return "-"
    elif sign == "-":
        return "+"
    else:
        return sign

def find_matched_transcripts(transcripts, inv_records, first_exon, last_exon):
    print(transcripts)
    print(inv_records)
    res={}
    for inv_exon, prev_next_exons in inv_records.items():
        print("loop:",inv_exon,prev_next_exons )
        res[inv_exon]=[]
        for transcript_id, transcript_info in transcripts.items():
            # print("transcript_id:", transcript_id)
            # print("transcript_info:", transcript_info)
            exons = transcript_info["exons"]
            # print("exons:", len(exons))
            strand_direction = transcript_info["strand"]
            prev_exon=prev_next_exons["prev_exon"]
            if prev_exon == None:
                prev_exon=first_exon
            prev_exon=(inv_exon[0], prev_exon[0], prev_exon[1], prev_exon[2])
            next_exon=prev_next_exons["next_exon"]
            if next_exon==None:
                next_exon=last_exon
            next_exon=(inv_exon[0], next_exon[0], next_exon[1], next_exon[2])
            current_exons=(inv_exon[0], inv_exon[1], inv_exon[2], inv_exon[3])
            print(prev_exon, next_exon, exons)
            if prev_exon in exons and next_exon in exons:
                print("find matched transcript:", transcript_id)
                res[inv_exon].append(transcript_id)
            prev_exon=(inv_exon[0], prev_exon[1], prev_exon[2], reverse_sign(prev_exon[3]))
            next_exon=(inv_exon[0], next_exon[1], next_exon[2], reverse_sign(next_exon[3]))
            if prev_exon in exons and next_exon in exons:
                print("find matched transcript:", transcript_id)
                res[inv_exon].append(transcript_id)
    return res

def assign_trf(matched_transcripts, transcripts):

    # for inv_exon, matched_trf_ids in matched_transcripts.items():
    #     inv_rec_id=1
    #     len_match_trf=len(matched_trf_ids)
    #     mean_coverage=inv_exon[4]/len_match_trf
        # if mean_coverage < 3:
        #     matched_trf_ids=matched_trf_ids[0]
    trans_id_inv_ids_dict={}
    for transcript_id, transcript_info in transcripts.items():
        for inv_exon, matched_trf_ids in matched_transcripts.items():
            if transcript_id in matched_trf_ids:
                if transcript_id not in trans_id_inv_ids_dict:
                    trans_id_inv_ids_dict[transcript_id]=[]
                trans_id_inv_ids_dict[transcript_id].append(inv_exon)
    # update transcripts and reverse the strand
    new_records={}
    new_records_ids=1
    for transcript_id, transcript_info in transcripts.items():
        print("transcript_info:", transcript_info)
        if transcript_id in trans_id_inv_ids_dict and trans_id_inv_ids_dict[transcript_id]:
            new_exons=transcript_info["exons"]
            for inversion_exon in trans_id_inv_ids_dict[transcript_id]:
                new_exons.append((inversion_exon[0], inversion_exon[1], inversion_exon[2], reverse_sign(transcript_info["strand"])))
            # sort the exons by start position
            new_exons.sort(key=lambda x: x[1])
            # get index of inversion exons
            inv_exon_index=[]
            for i, exon in enumerate(new_exons):
                for inversion_exon in trans_id_inv_ids_dict[transcript_id]:
                    if (inversion_exon[0], inversion_exon[1], inversion_exon[2], reverse_sign(transcript_info["strand"])) == exon:
                        inv_exon_index.append(i)
            # insert the inversion exons coverage at inv_exon_index
            total_coverage=0
            for inversion_exon in trans_id_inv_ids_dict[transcript_id]:
                total_coverage+=inversion_exon[4]
            print("total_coverage:", total_coverage, len(trans_id_inv_ids_dict[transcript_id]))
            mean_coverage=total_coverage/len(trans_id_inv_ids_dict[transcript_id])
            new_exon_covs=(len(transcript_info["exon_coverage"])+len(trans_id_inv_ids_dict[transcript_id]))*[mean_coverage]
            print("new_exon_covs", new_exon_covs)
            new_records[new_records_ids]={
                "gene_id": new_records_ids,
                "strand": transcript_info["strand"],
                "exons": new_exons,
                "exon_coverage": new_exon_covs,
                "transcript_coverage": mean_coverage,
                "FPKM": 0.0,
                "TPM": 0.0,
            }
            new_records_ids+=1
            # reduce the coverage of the original transcript
            transcript_info["transcript_coverage"]-=mean_coverage
            # reduce the coverage of the all exons not only inversion exons
            for exon_id, exon_cov in enumerate(transcript_info["exon_coverage"]):
                print("exon_cov:", exon_cov)
                start, end, cov =exon_cov
                cov-=mean_coverage
                transcript_info["exon_coverage"][exon_id]=(start, end, cov)
            transcripts[transcript_id]=transcript_info
            
                
                
    print("new_records:", new_records)
    print("transcripts:", transcripts)
    return new_records, transcripts




def process_transcript(exons, ref_fasta, vcf_file, haplotype, tr_strand):
    """
    handle a single transcript, extract and concatenate all exon sequences,
    """
    full_seq = ""
    
    for chrom, start, end, strand in exons:
        region = f"{chrom}:{start}-{end}"
        region_ = f"{chrom}_{start}_{end}"
        
        # construct commands to extract reference sequence and apply VCF variants
        cmd_samtools = f"samtools faidx {ref_fasta} {region}"
        cmd_bcftools = f"bcftools consensus -H {haplotype} {vcf_file}"
        
        print(f"Executing: {cmd_samtools} | {cmd_bcftools}")
        
        # run the command and get output
        result = subprocess.getoutput(f"{cmd_samtools} | {cmd_bcftools}")
        print("Raw result:", result)

        # parse output and remove non-sequence lines
        lines = result.split("\n")
        seq_lines = []
        for line in lines:
            if line.startswith("Warning:") or line.startswith("Note:"):
                continue
            if line.startswith(">"):
                continue
            if line.startswith("Applied"):
                continue
            seq_lines.append(line)
        
        seq = "".join(seq_lines)
        print("Extracted sequence:", seq)

        # if strand is different from transcript strand, reverse complement the sequence
        print("strand:", strand, "tr_strand:", tr_strand)
        if strand != tr_strand:
            seq = str(Seq(seq).reverse_complement())

            # store the inversion haps
            inv_output_file = f"{args.output}/only_inv_haps/{region_}_inv_hap{haplotype}.fa"
            with open(inv_output_file, "w") as f:
                f.write(f">{region_}_inv_hap{haplotype}\n{seq}\n")
            print(f"Saved inverted sequence to: {inv_output_file}")

        # concatenate sequences
        full_seq += seq
        print("Current full_seq:", full_seq)

    return full_seq

def process_normal_transcript(exons, ref_fasta, vcf_file, haplotype):
    """
    handle a single transcript, extract and concatenate all exon sequences,
    while considering strand and inversion.
    """
    full_seq = ""
    
    for chrom, start, end, strand in exons:
        region = f"{chrom}:{start}-{end}"
        
        # construct commands to extract reference sequence and apply VCF variants
        cmd_samtools = f"samtools faidx {ref_fasta} {region}"
        cmd_bcftools = f"bcftools consensus -H {haplotype} {vcf_file}"
        
        print(f"Executing: {cmd_samtools} | {cmd_bcftools}")
        
        try:
            # run the command and get output
            result = subprocess.getoutput(f"{cmd_samtools} | {cmd_bcftools}")
            print("Raw result:", result)

            # parse output and remove non-sequence lines
            lines = result.split("\n")
            seq_lines = []
            for line in lines:
                if line.startswith("Warning:") or line.startswith("Note:"):
                    continue
                if line.startswith(">"):  
                    continue
                if line.startswith("Applied"): 
                    continue
                seq_lines.append(line)

            seq = "".join(seq_lines)
            print("Extracted sequence:", seq)

            # concatenate sequences
            full_seq += seq

        except Exception as e:
            print(f"Error processing region {region}: {e}")
            continue  # skip current exon on error
    
    return full_seq

def merge_consistent_blk(blk_list):
    result = []
    current_group = [blk_list[0]]

    # loop through the rest of the items
    for item in blk_list[1:]:
        # if the symbol is the same as the last item in the current group, add it to the group
        if item[3] == current_group[-1][3]:
            current_group.append(item)
        else:
            # if the symbol changes, add the current group to the result list and start a new group
            result.append(current_group)
            current_group = [item]

    # add the last group to the result
    result.append(current_group)

    # print the result
    for i, group in enumerate(result):
        print(f"Group {i+1}: {group}")
    return result

def generate_seq(new_records, old_records, phased_vcf, ref, output_dir):
    os.makedirs(f"{output_dir}/only_inv_haps", exist_ok=True)
    if args.with_normal_haps:
        # make dir for normal haps
        os.makedirs(f"{output_dir}/normal_haps", exist_ok=True)
        out_normal_fa_hap1=f"{output_dir}/normal_haps/transcripts_hap1.fa"
        out_normal_fa_hap2=f"{output_dir}/normal_haps/transcripts_hap2.fa"
        out_normal_fa_hap1_f=open(out_normal_fa_hap1, "w")
        out_normal_fa_hap2_f=open(out_normal_fa_hap2, "w")
    out_fa_hap1=f"{output_dir}/transcripts_hap1.fa"
    out_fa_hap2=f"{output_dir}/transcripts_hap2.fa"
    out_fa_hap1_f=open(out_fa_hap1, "w")
    out_fa_hap2_f=open(out_fa_hap2, "w")
    print("writing")
    for old_record_id, old_record in old_records.items():
        transcript_id = old_record_id
        exons = old_record["exons"]
        tr_strand = old_record["strand"]
        grouped_exons = merge_consistent_blk(exons)
        full_hap1_seq = ""
        full_hap2_seq = ""
        if args.with_normal_haps:
            full_hap1_normal_seq = ""
            full_hap2_normal_seq = ""
        print("grouped_exons:", grouped_exons)
        for group in grouped_exons:
            print("processing")
            hap1_seq=process_transcript(group, ref, phased_vcf, "1", tr_strand)
            hap2_seq=process_transcript(group, ref, phased_vcf, "2", tr_strand)
            full_hap1_seq+=hap1_seq
            full_hap2_seq+=hap2_seq
            if args.with_normal_haps:
                hap1_normal_hap=process_normal_transcript(group, ref, phased_vcf, "1")
                hap2_normal_hap=process_normal_transcript(group, ref, phased_vcf, "2")
                full_hap1_normal_seq+=hap1_normal_hap
                full_hap2_normal_seq+=hap2_normal_hap
        if tr_strand == "-":
            full_hap1_seq = str(Seq(full_hap1_seq).reverse_complement())
            full_hap2_seq = str(Seq(full_hap2_seq).reverse_complement())
            if args.with_normal_haps:
                full_hap1_normal_seq = str(Seq(full_hap1_normal_seq).reverse_complement())
                full_hap2_normal_seq = str(Seq(full_hap2_normal_seq).reverse_complement())
        out_fa_hap1_f.write(f">{transcript_id}_hap1\n{full_hap1_seq}\n")
        out_fa_hap2_f.write(f">{transcript_id}_hap2\n{full_hap2_seq}\n")
        if args.with_normal_haps:
            out_normal_fa_hap1_f.write(f">{transcript_id}_hap1\n{full_hap1_normal_seq}\n")
            out_normal_fa_hap2_f.write(f">{transcript_id}_hap2\n{full_hap2_normal_seq}\n")

    print("new records")
    print(new_records)
    for new_record_id, new_record in new_records.items():
        transcript_id = new_record_id
        exons = new_record["exons"]
        tr_strand = new_record["strand"]
        grouped_exons = merge_consistent_blk(exons)
        full_hap1_seq = ""
        full_hap2_seq = ""
        if args.with_normal_haps:
            full_hap1_normal_seq = ""
            full_hap2_normal_seq = ""
        for group in grouped_exons:
            hap1_seq=process_transcript(group, ref, phased_vcf, "1", tr_strand)
            hap2_seq=process_transcript(group, ref, phased_vcf, "2", tr_strand)
            full_hap1_seq+=hap1_seq
            full_hap2_seq+=hap2_seq
            print("group:", group)
            print("hap1_seq:", hap1_seq)
            print("hap2_seq:", hap2_seq)
            if args.with_normal_haps:
                hap1_normal_hap=process_normal_transcript(group, ref, phased_vcf, "1")
                hap2_normal_hap=process_normal_transcript(group, ref, phased_vcf, "2")
                full_hap1_normal_seq+=hap1_normal_hap
                full_hap2_normal_seq+=hap2_normal_hap
        # if tr_strand == "-":
        #     full_hap1_seq = str(Seq(full_hap1_seq).reverse_complement())
        #     full_hap2_seq = str(Seq(full_hap2_seq).reverse_complement())
        #     if args.with_normal_haps:
        #         full_hap1_normal_seq = str(Seq(full_hap1_normal_seq).reverse_complement())
        #         full_hap2_normal_seq = str(Seq(full_hap2_normal_seq).reverse_complement())
        print("final_full_hap1_seq:", full_hap1_seq)
        print("final_full_hap2_seq:", full_hap2_seq)

        out_fa_hap1_f.write(f">new_{transcript_id}_inv_hap1\n{full_hap1_seq}\n")
        out_fa_hap2_f.write(f">new_{transcript_id}_inv_hap2\n{full_hap2_seq}\n")
        if args.with_normal_haps:
            out_normal_fa_hap1_f.write(f">new_{transcript_id}_invbefore_hap1\n{full_hap1_normal_seq}\n")
            out_normal_fa_hap2_f.write(f">new_{transcript_id}_invbefore_hap2\n{full_hap2_normal_seq}\n")
    out_fa_hap1_f.close()
    out_fa_hap2_f.close()
    if args.with_normal_haps:
        out_normal_fa_hap1_f.close()
        out_normal_fa_hap2_f.close()

                
def main():
    

    # read files
    inversion_exons = read_inversion_exons(args.inversion_exons)
    junc_connections = read_junc_file(args.junc)
    transcripts = read_gtf_file2(args.gtf)
    exons=[]
    first_exon=("100000000000","100000000000","+","0")
    last_exon=("0","0","+","0")
    for trs_id, trsc_rec in transcripts.items():
        trsc_first_exon =trsc_rec["exons"][0]
        trsc_first_exon_cov_item=trsc_rec["exon_coverage"][0]
        trsc_last_exon=trsc_rec["exons"][-1]
        trsc_last_exon_cov_item=trsc_rec["exon_coverage"][-1]
        if int(trsc_first_exon[1])<= int(first_exon[0]):
            first_exon=(str(trsc_first_exon[1]), str(trsc_first_exon[2]), str(trsc_first_exon[3]), str(trsc_first_exon_cov_item[2]))
        if int(trsc_last_exon[1]) >= int(last_exon[0]):
            last_exon=(str(trsc_first_exon[1]), str(trsc_last_exon[2]), trsc_last_exon[3], str(trsc_last_exon_cov_item[2]))



    # get prev and next link for inversion exons

    inv_records = find_flanking_exons_single_pair(inversion_exons, junc_connections)
    print("inv records:", inv_records)

    matched_transcripts = find_matched_transcripts(transcripts, inv_records, first_exon, last_exon)
    print("matched transcripts:", matched_transcripts)

    # get new transcripts and updated, if inversion exons has the same matched transcripts, reverse the strand at same time
    new_records, old_records = assign_trf(matched_transcripts, transcripts)
    print("new records:", new_records)
    print("old records:", old_records)
    
    generate_seq(new_records, old_records, args.vcf, args.ref, args.output)


if __name__ == "__main__":
    # use argparse to get input files and output dir
    parser = argparse.ArgumentParser(description="Process inversion exons, junction files, and GTF to split transcripts.")
    parser.add_argument('-i', '--inversion_exons', required=True, help="Path to the inversion.exons file.")
    parser.add_argument('-j', '--junc', required=True, help="Path to the inv junctions.junc file.")
    parser.add_argument('-g', '--gtf', required=True, help="Path to the GTF file.")
    parser.add_argument('-o', '--output', required=True, help="Path to the output dir.")
    parser.add_argument('-v', '--vcf', required=True, help="Path to the phased VCF file.")
    parser.add_argument('-r', '--ref', required=True, help="Path to the reference genome.")
    parser.add_argument("--with_normal_haps", required=False, default=False, action="store_true", help="Whether to use normal haps.")


    
    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)
    main()