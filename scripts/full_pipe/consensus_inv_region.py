import argparse
import os
import subprocess


def parse_region(region_str):
    """
    Parse a region string like 'chr21:27374157-27374194-6.63' into components.
    Returns: (chromosome, start, end, depth)
    """
    chrom, coords_depth = region_str.split(":")
    coords, depth = coords_depth.rsplit("-", 1)
    start, end = map(int, coords.split("-"))
    depth = float(depth)
    return chrom, start, end, depth


def merge_regions(region1, region2):
    """
    Merge two regions into one.
    Calculate new start, end, and average depth.
    """
    chrom1, start1, end1, depth1 = region1
    chrom2, start2, end2, depth2 = region2

    # Ensure both regions are on the same chromosome
    if chrom1 != chrom2:
        raise ValueError("Regions are on different chromosomes!")

    new_start = min(start1, start2)
    new_end = max(end1, end2)
    new_depth = (depth1 + depth2) / 2

    return chrom1, new_start, new_end, new_depth


def calculate_average_depth_with_samtools_and_awk(bam_file, chrom, start, end):
    """
    Calculate the average depth of the region using samtools depth and awk.
    """
    region = f"{chrom}:{start}-{end}"
    depth_command = f"samtools depth -a -r {region} {bam_file} | awk '{{sum+=$3}} END {{print sum/NR}}'"

    try:
        # Run samtools depth with awk to calculate the average depth
        result = subprocess.run(
            depth_command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,  # Replace text=True with universal_newlines=True for Python 3.6 compatibility
            check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error running samtools depth for region {region}: {e.stderr}")
        print(f"Command: {depth_command}")
        return 0.0

    # Parse the result
    try:
        average_depth = float(result.stdout.strip())
    except ValueError:
        print(f"Error parsing average depth for region {region}: {result.stdout}")
        return 0.0

    return average_depth


def should_merge(region1, region2, bam_file, distance_threshold, depth_threshold):
    """
    Determine whether two regions should be merged.
    """
    _, _, end1, _ = region1
    _, start2, _, _ = region2

    distance = start2 - end1

    # If the regions are close enough, merge directly
    if distance <= distance_threshold:
        return True

    # Otherwise, check the average depth in the gap
    chrom, _, _, _ = region1
    avg_depth = calculate_average_depth_with_samtools_and_awk(bam_file, chrom, end1, start2)
    return avg_depth >= depth_threshold


def process_gene_file(gene_file, bam_file, distance_threshold, depth_threshold):
    """
    Process the gene.txt file and merge regions according to the rules.
    """
    with open(gene_file, "r") as f:
        lines = f.readlines()

    # Parse the regions from the file
    parsed_lines = []
    for line in lines:
        parts = line.strip().split()
        region_info = parts[2]
        chrom, start, end, depth = parse_region(region_info)
        parsed_lines.append((line.strip(), (chrom, start, end, depth)))

    # Sort regions by start coordinate
    parsed_lines.sort(key=lambda x: x[1][1])  # Sort by start position

    # Merge regions
    merged_lines = []
    inversion_lines = []
    current_line, current_region = parsed_lines[0]

    for next_line, next_region in parsed_lines[1:]:
        if should_merge(current_region, next_region, bam_file, distance_threshold, depth_threshold):
            # Merge regions and update the current region
            current_region = merge_regions(current_region, next_region)
        else:
            # Append the current line with the updated region
            chrom, start, end, depth = current_region
            updated_line = f"{current_line.split()[0]} {current_line.split()[1]} {chrom}:{start}-{end}-{depth:.2f} {current_line.split()[-1]}"
            merged_lines.append(updated_line)
            inversion_lines.append(f"{chrom}\t{start}\t{end}\t{depth:.2f}")

            # Move to the next region
            current_line, current_region = next_line, next_region

    # Add the last region
    chrom, start, end, depth = current_region
    updated_line = f"{current_line.split()[0]}\t{current_line.split()[1]}\t{chrom}:{start}-{end}-{depth:.2f}\t{current_line.split()[-1]}"
    inversion_lines.append(f"{chrom}\t{start}\t{end}\t{depth:.2f}")
    merged_lines.append(updated_line)

    return merged_lines, inversion_lines


def write_filter_gene(output_file, merged_lines):
    """
    Write the merged lines to the output file.
    """
    with open(output_file, "w") as f:
        for line in merged_lines:
            f.write(line + "\n")

def write_inverson(output_file, inversion_lines):
    
    with open(output_file.replace("gene", "inv_exon"), "w") as f:
        for line in inversion_lines:
            f.write(line + "\n")



def main():
    

    print(f"Processing {args.gene_file}...")
    merged_lines, inversion_lines = process_gene_file(
        args.gene_file, args.bam_file, args.distance_threshold, args.depth_threshold
    )
    write_filter_gene(args.output_file, merged_lines)
    print(f"Filtered gene saved to {args.output_file}")
    write_inverson(args.inv_exon_file, inversion_lines)
    print(f"Filtered inv_exon saved to {args.inv_exon_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge gene records based on distance and BAM coverage.")
    parser.add_argument("--gene_file", required=True, help="Path to the gene.txt file.")
    parser.add_argument("--bam_file", required=True, help="Path to the BAM file for coverage analysis.")
    parser.add_argument("--output_file", required=True, help="filtered gene_file")
    parser.add_argument("--inv_exon_file", required=True, help="Path to the inv_exon.txt file.")
    parser.add_argument("--distance_threshold", type=int, default=50, help="Distance threshold to merge regions.")
    parser.add_argument("--depth_threshold", type=float, default=3.0, help="Depth threshold for merging regions.")
    args = parser.parse_args()
    main()