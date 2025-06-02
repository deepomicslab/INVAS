import pysam
import argparse

def parse_gtf(gtf_file):
    exons = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if parts[2] == 'exon':
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                exon = (start, end)
                exons.setdefault(chrom, []).append(exon)
    for chrom in exons:
        exons[chrom].sort(key=lambda x: x[0])
    return exons

def find_nearest_exon(exons, chrom, position):
    nearest_exon = None
    min_dist = float('inf')
    for exon in exons.get(chrom, []):
        exon_start, exon_end = exon
        if exon_end < position:
            dist = position - exon_end
        elif exon_start > position:
            dist = exon_start - position
        else:
            return exon
        if dist < min_dist:
            min_dist = dist
            nearest_exon = exon
    return nearest_exon

def create_exon_chain(exons, chrom, blocks):
    chain = []
    for block_start, block_end in blocks:
        start_exon = find_nearest_exon(exons, chrom, block_start)
        end_exon = find_nearest_exon(exons, chrom, block_end)
        # Choose the closer exon for the start of the block
        if start_exon and end_exon:
            if abs(start_exon[0] - block_start) < abs(end_exon[1] - block_end):
                chosen_exon = start_exon
            else:
                chosen_exon = end_exon
        elif start_exon:
            chosen_exon = start_exon
        elif end_exon:
            chosen_exon = end_exon
        else:
            continue
        # Append if not already in the chain
        if not chain or chain[-1] != chosen_exon:
            chain.append(chosen_exon)
    return tuple(chain)

def is_alternative_splicing(blocks):
    # Check if there are multiple blocks separated by at least the minimum intron size
    # and less than the maximum intron size.
    for i in range(len(blocks) - 1):
        intron_size = blocks[i+1][0] - blocks[i][1]
        if intron_size > min_intron_size and intron_size < max_intron_size:
            return True
    return False


def main(bam_file, gtf_file):
    exons = parse_gtf(gtf_file)
    exon_chains = {}
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            
            chrom = bam.get_reference_name(read.reference_id)
            blocks = read.get_blocks()
            if not is_alternative_splicing(blocks):
                continue
            # if len(blocks) >= 2:
            #     print(blocks, read.query_name)
            exon_chain = create_exon_chain(exons, chrom, blocks)
            if exon_chain == ((71139711, 71139848), (71284220, 71284399)):
                print(read.query_name, blocks)
            if exon_chain:
                exon_chains[exon_chain] = exon_chains.get(exon_chain, 0) + 1
    
    # Print the exon chains and their counts
    for chain, count in exon_chains.items():
        print(f"Exon chain {chain} occurs {count} times")

if __name__ == "__main__":
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Track exon chains in sequencing reads')
    parser.add_argument('bam_file', help='Path to the BAM file')
    parser.add_argument('gtf_file', help='Path to the GTF file')
    min_intron_size=50
    max_intron_size=100000
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Run the main function with the provided arguments
    main(args.bam_file, args.gtf_file)