import argparse
import pysam

parser = argparse.ArgumentParser(description='Find continuous regions with average depth above a certain threshold, considering gaps.')
parser.add_argument('bam', type=str, help='Path to the BAM file')
parser.add_argument('chromosome', type=str, help='Chromosome to analyze')
parser.add_argument('threshold', type=float, help='Depth threshold')
parser.add_argument('--gap', type=int, default=1000, help='Maximum gap length of low depth allowed within a region')
args = parser.parse_args()

bam_path = args.bam
chromosome = args.chromosome
depth_threshold = args.threshold
max_gap_length = args.gap

bam_file = pysam.AlignmentFile(bam_path, "rb")

start = None
end = None
current_region_depths = []
gap_counter = 0
last_pos = -1

for pileupcolumn in bam_file.pileup(chromosome):
    pos = pileupcolumn.pos
    depth = pileupcolumn.nsegments

    if last_pos != -1 and pos != last_pos + 1:
        skipped_length = pos - (last_pos + 1)
        if start is not None:
            if skipped_length <= max_gap_length:
                current_region_depths.extend([0] * skipped_length)
                end = pos - 1
                gap_counter += skipped_length
            else:
                average_depth = sum(current_region_depths) / len(current_region_depths)
                if average_depth >= depth_threshold:
                    print(f"{chromosome}\t{start}\t{end}\t{average_depth:.2f}")
                start = None
                current_region_depths = []
                gap_counter = 0

    last_pos = pos

    if start is None and depth >= depth_threshold:
        start = pos  

    if start is not None:
        if depth >= depth_threshold:
            current_region_depths.append(depth)
            end = pos
            gap_counter = 0 
        else:
            gap_counter += 1
            if gap_counter <= max_gap_length:
                current_region_depths.append(depth)
                end = pos
            else:
                average_depth = sum(current_region_depths) / len(current_region_depths)
                if average_depth >= depth_threshold:
                    print(f"{chromosome}\t{start}\t{end}\t{average_depth:.2f}")
                start = None
                current_region_depths = []
                gap_counter = 0

if start is not None:
    average_depth = sum(current_region_depths) / len(current_region_depths)
    if average_depth >= depth_threshold:
        print(f"{chromosome}\t{start}\t{end}\t{average_depth:.2f}")

bam_file.close()
