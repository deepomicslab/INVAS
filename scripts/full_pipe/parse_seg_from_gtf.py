import csv

def parse_stringtie_gtf(gtf_file, output_file):
    """
    parse a StringTie GTF file to extract exon information and write to a TSV file.
    """
    with open(gtf_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        # write header
        writer.writerow(['chr', 'exon_start', 'exon_end', 'gene_id', 'strand'])
        
        for line in infile:
            # skip header lines
            if line.startswith('#'):
                continue
            
            # split the line into fields
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue  # insufficient fields, skip

            # extract fields
            chrom = fields[0]          # chromosome
            feature_type = fields[2]   # feature type (e.g., exon, transcript)
            start = fields[3]          # start position
            end = fields[4]            # end position
            strand = fields[6]         # strand (+ or -)
            attributes = fields[8]     # attributes field

            # filter to only extract "exon" type lines
            if feature_type != 'exon':
                continue

            # extract gene_id from attributes field
            gene_id = None
            for attribute in attributes.split(';'):
                attribute = attribute.strip()
                if attribute.startswith('gene_id'):
                    # extract gene_id value (remove quotes)
                    gene_id = attribute.split(' ')[1].replace('"', '')
                    break

            # ensure gene_id exists
            if gene_id is None:
                continue

            # write extracted data to output file
            writer.writerow([chrom, start, end, gene_id, strand])

    print(f"exon information has been written to: {output_file}")


# Example usage
# Input GTF file path
gtf_file = 'stringtie/stringtie.app.a2.gtf'
# Output file path
output_file = 'exons_output.tsv'

# call the function to parse GTF file
parse_stringtie_gtf(gtf_file, output_file)