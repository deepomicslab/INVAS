import vcf
import sys

def vcf_to_bed(vcf_filename, bed_filename, threshold):
    vcf_reader = vcf.Reader(filename=vcf_filename)
    with open(bed_filename, 'w') as bed:
        for record in vcf_reader:
            chrom = record.CHROM
            pos = record.POS
            # Check if the SVTYPE is inversion
            if record.INFO.get('SVTYPE') == 'INV':
                # Extract end position if available, otherwise use start position
                end_pos = record.INFO.get('END', pos)
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