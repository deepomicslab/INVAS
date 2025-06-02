import argparse
import re

def parse_unique_genes_from_file(file_path):
    """
    Parse the file and extract unique gene names from the second column.

    Parameters:
        file_path (str): Path to the input file.

    Returns:
        set: A set of unique gene names.
    """
    unique_genes = set()  # To store unique gene names
    gene_region_dict = {}  # To store gene regions
    chrom=""

    with open(file_path, "r") as file:
        for line in file:
            # Split each line into columns by tab
            columns = line.strip().split("\t")
            if len(columns) > 1:
                # Extract the second column (e.g., 'PTPRF:43990857-44089343')
                gene_info = columns[1]
                # Extract the gene name before the colon
                gene_name = gene_info.split(":")[0]
                unique_genes.add(gene_name)
                # Extract the gene region
                # split gene_info by "," and "-" and ":" in the same time
                resion_info=gene_info.split(":")[1]
                res=re.split(r'[,:-]', resion_info)
                gene_start, gene_end = res[0], res[1]
                chrom = columns[0].split(":")[0]
                gene_region_dict[gene_name] = (chrom, int(gene_start), int(gene_end))

    return unique_genes, gene_region_dict, chrom


if __name__ == "__main__":
    # Parse the input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input_file", help="Path to the input file", required=True
    )
    parser.add_argument(
        "-o", "--output_file", help="Path to the output file", required=True
    )
    parser.add_argument(
        "-e", "--extend", help="Extend the gene region by the specified number of base pairs, default is 5000", type=int, default=5000
    )
    args = parser.parse_args()

    # Parse the input file and extract unique gene names
    unique_genes, gene_region_dict, chrom = parse_unique_genes_from_file(args.input_file)

    # Write the unique gene names to the output file
    with open(args.output_file, "w") as file:
        for gene in unique_genes:
            chrom, gene_start, gene_end = gene_region_dict[gene]
            file.write(f"{gene}\t{chrom}:{gene_start-args.extend}-{gene_end+args.extend}\n")