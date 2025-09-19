import sys
from pybedtools import BedTool

def read_bed_file(bed_file_path):
    return BedTool(bed_file_path)

def read_inversion_file(inversion_file_path):
    with open(inversion_file_path, 'r') as file:
        inversions = [line.strip() for line in file.readlines()]
    return inversions

def check_inversions(inversions, high_conf_bed):
    overlapping_inversions = []

    for inv in inversions:
        chr_and_pos = inv.split('\t')[0]
        chrom, positions = chr_and_pos.split(':')
        start, end = positions.split('-')

        inversion_bed = BedTool(f"{chrom} {start} {end}", from_string=True)

        # Check for intersections and get the actual overlapping regions
        intersections = inversion_bed.intersect(high_conf_bed, wa=True, wb=True)

        for intersection in intersections:
            # Format the overlapping region with chromosome information
            overlap_region = f"{intersection[3]}:{intersection[4]}-{intersection[5]}-{intersection[6]}"
            inversion_with_overlap = f"{inv}\t{overlap_region}"
            overlapping_inversions.append(inversion_with_overlap)

    return overlapping_inversions

def main():
    if len(sys.argv) != 4:
        print("Usage: python filter_inversions.py <high_confidence_bed> <inversion_file> <output_file>")
        sys.exit(1)

    high_confidence_bed_path = sys.argv[1]
    inversion_file_path = sys.argv[2]
    output_file_path = sys.argv[3]

    high_conf_bed = read_bed_file(high_confidence_bed_path)

    inversions = read_inversion_file(inversion_file_path)

    overlapping_inversions = check_inversions(inversions, high_conf_bed)

    with open(output_file_path, 'w') as output_file:
        for inv in overlapping_inversions:
            output_file.write(inv + '\n')

    print(f"Filtered inversions with overlapping high confidence regions written to {output_file_path}")

if __name__ == "__main__":
    main()
