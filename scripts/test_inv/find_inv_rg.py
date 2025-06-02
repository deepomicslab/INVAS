import pysam
from pybedtools import BedTool
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Calculate M regions in BAM file.')
    parser.add_argument('-i', '--input', help='Input BAM file', required=True)
    parser.add_argument('-o', '--output', help='Output BED file', required=True)
    return parser.parse_args()

def main():
    args = parse_args()
    # Open the BAM file
    bamfile = pysam.AlignmentFile(args.input, "rb")

    # Initialize an empty list to store the match regions
    match_regions = []

    # For each read in the BAM file
    for read in bamfile:
        # If the read is mapped
        if not read.is_unmapped:
            # Get the start position of the read
            start = read.pos
            # Initialize the current position to the start
            current = start
            # For each operation in the CIGAR string
            for op, length in read.cigar:
                # If the operation is a match or mismatch
                if op == 0:
                    # Append the start and end of the match/mismatch to the list
                    match_regions.append((read.reference_name, current, current+length))
                    # Update the current position
                    current += length
                # If the operation is a deletion or skip, update the current position
                elif op in [2, 3]:
                    current += length

    # Close the BAM file
    bamfile.close()

    # Convert the match regions to a BedTool
    match_regions = BedTool(match_regions)

    # Merge the overlapping match regions
    merged_regions = match_regions.sort().merge()

    # Save the merged match regions to the output file
    merged_regions.saveas(args.output)

if __name__ == '__main__':
    main()