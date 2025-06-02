import re
import sys

def vcf_to_bed(vcf_filename, bed_filename, threshold):
    with open(vcf_filename, 'r') as vcf, open(bed_filename, 'w') as bed:
        for line in vcf:
            if line.startswith('#'):
                continue  # Skip header lines
            cols = line.strip().split('\t')
            chrom = cols[0]
            pos = int(cols[1])  # Convert position to integer
            info = cols[7]
            # Adjust these regular expressions based on how inversions are represented in your VCF
            if re.search(r'SVTYPE=INV', info):  # Check if the SVTYPE is inversion
                # Extract end position if available
                end_match = re.search(r'END=(\d+)', info)
                if end_match:
                    end_pos = int(end_match.group(1))  # Convert end position to integer
                else:
                    end_pos = pos  # If END is not specified, use the start position as end position

                # Check if the inversion size is under the threshold
                if end_pos - pos <= threshold:
                    # Write to BED: chrom, start, end (0-based start, 1-based end)
                    bed.write(f"{chrom}\t{pos - 1}\t{end_pos}\n")

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python vcf_to_bed.py <input_vcf_file> <output_bed_file> <size_threshold>")
        sys.exit(1)

    vcf_filename = sys.argv[1]
    bed_filename = sys.argv[2]
    size_threshold = int(sys.argv[3])  # Convert threshold to integer

    vcf_to_bed(vcf_filename, bed_filename, size_threshold)