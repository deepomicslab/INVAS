import pysam
import re
import argparse

# Function to parse CIGAR strings and return soft clipping information
def parse_cigar(cigar_string):
    pattern = re.compile(r'(\d+)([SH])')
    matches = pattern.findall(cigar_string)
    left_clip = int(matches[0][0]) if matches and matches[0][1] == 'S' else 0
    right_clip = int(matches[-1][0]) if matches and matches[-1][1] == 'S' else 0
    return left_clip, right_clip

# Parsing the region string to extract multiple sub-regions and their associated metadata
def parse_regions(region_str):
    parts = region_str.split('\t')
    main_region = parts[0]
    additional_info = parts[1] if len(parts) > 1 else ""
    sub_regions = additional_info.split(',')
    return main_region, sub_regions

# Count reads that are soft-clipped beyond a threshold and overlap with given regions
def count_soft_clipped_reads(bam_file, regions, threshold):
    soft_clipped_counts = {region: 0 for region in regions}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for region in regions:
            main_region, sub_regions = parse_regions(region)
            chrom, start_end = main_region.split(':')
            start, end = map(int, start_end.split('-'))
            sub_region_tuples = [(int(sr.split('-')[0]), int(sr.split('-')[1])) for sr in sub_regions]
            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped:
                    continue
                left_clip, right_clip = parse_cigar(read.cigarstring)
                if left_clip < threshold and right_clip < threshold:
                    continue
                read_blocks = read.get_blocks()
                if any(block[0] < sr_end and block[1] > sr_start for sr_start, sr_end in sub_region_tuples for block in read_blocks):
                    soft_clipped_counts[region] += 1
    return soft_clipped_counts

# Main function to read region file and count soft clipping reads
def main(region_file, bam_file, threshold):
    with open(region_file, 'r') as file:
        regions = [line.strip() for line in file]
    counts = count_soft_clipped_reads(bam_file, regions, threshold)
    for region, count in counts.items():
        print(f"Region: {region}, Soft clipped reads: {count}")

# Parse command-line arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Count soft-clipped reads that overlap with specified regions in a BAM file and exceed a length threshold.')
    parser.add_argument('region_file', type=str, help='File containing regions with format "chr:start-end" and optional additional sub-regions and metadata')
    parser.add_argument('bam_file', type=str, help='Path to the BAM file')
    parser.add_argument('threshold', type=int, help='Length threshold for soft clipping')
    args = parser.parse_args()
    main(args.region_file, args.bam_file, args.threshold)