import argparse

def read_bed(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            yield {
                'chrom': parts[0],
                'start': int(parts[1]),
                'end': int(parts[2]),
                'name': parts[3] if len(parts) > 3 else None,
                'strand': parts[5] if len(parts) > 5 else None
            }

def find_inversions_in_genes(inversion_path, gene_path, output_path):
    genes = list(read_bed(gene_path))

    with open(output_path, 'w') as output_file:
        for inversion in read_bed(inversion_path):
            for gene in genes:
                if (inversion['chrom'] == gene['chrom'] and
                        inversion['start'] >= gene['start'] and
                        inversion['end'] <= gene['end']):
                    output_file.write(f"{inversion['chrom']}\t{inversion['start']}\t{inversion['end']}\t"
                                      f"{gene['chrom']}\t{gene['start']}\t{gene['end']}\t"
                                      f"{gene['name']}\t{gene['strand']}\n")

def main():
    parser = argparse.ArgumentParser(description="Find inversions within gene regions.")
    parser.add_argument("inversion_bed", help="Path to inversion BED file.")
    parser.add_argument("gene_bed", help="Path to gene BED file.")
    parser.add_argument("output_bed", help="Path to output BED file.")

    args = parser.parse_args()

    find_inversions_in_genes(args.inversion_bed, args.gene_bed, args.output_bed)

if __name__ == "__main__":
    main()