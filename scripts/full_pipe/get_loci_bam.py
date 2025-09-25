import os
import subprocess
import argparse


def split_input_by_gene(input_file, output_dir):
    """
    将输入文件按照基因名称分割到各自的文件夹中。
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # store file handles for each gene to avoid reopening files
    gene_files = {}

    with open(input_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            # extract gene name from the second column
            columns = line.split('\t')
            if len(columns) < 2:
                print(f"Warning: invalid line format, skipping: {line}")
                continue

            try:
                gene_info = columns[1]
                gene_name = gene_info.split(':')[0]
            except ValueError:
                print(f"colums: {columns}")
                print(f"Warning: cannot parse gene name, skipping: {columns[1]}")
                continue

            # create gene folder if not exists
            gene_folder = os.path.join(output_dir, gene_name)
            if not os.path.exists(gene_folder):
                os.makedirs(gene_folder)

            # open gene file if not already opened
            gene_file_path = os.path.join(gene_folder, f"gene.txt")
            if gene_name not in gene_files:
                gene_files[gene_name] = open(gene_file_path, 'w')

            # write to gene file
            gene_files[gene_name].write(line + '\n')

    # close all file handles
    for gene_file in gene_files.values():
        gene_file.close()


def extract_bam_by_gene(output_dir, bam1, bam2, extend, samtools_path="samtools"):
    """
    extract regions from BAM files based on gene folders and region info.
    """
    for gene_name in os.listdir(output_dir):
        gene_folder = os.path.join(output_dir, gene_name)
        if not os.path.isdir(gene_folder):
            continue

        gene_file = os.path.join(gene_folder, f"gene.txt")
        if not os.path.exists(gene_file):
            print(f"Warning: gene file not found for {gene_name}, skipping...")
            continue

        # use a set to track processed regions to avoid duplicates
        processed_regions = set()

        with open(gene_file, 'r') as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue

                # extract chromosome and region info
                columns = line.split('\t')
                if len(columns) < 2:
                    print(f"Warning: invalid line format, skipping: {line}")
                    continue

                # extract chromosome info from the first column
                try:
                    chrom_region = columns[0]
                    chrom = chrom_region.split(':')[0]
                except ValueError:
                    print(f"Warning: cannot parse chromosome info, skipping: {columns[0]}")
                    continue

                # extract region info from the second column
                try:
                    gene_info = columns[1]
                    positions = gene_info.split(':')[1]
                    if "," in positions:
                        positions = positions.split(",")[0]
                    start, end = map(int, positions.split('-'))
                except ValueError:
                    print(columns)
                    print(f"Warning: cannot parse region info, skipping: {columns[1]}")
                    continue

                # extend the region
                extended_start = max(0, start - extend)
                extended_end = end + extend

                # create a unique key for the region
                region_key = (chrom, extended_start, extended_end)

                # if this region has already been processed, skip it
                if region_key in processed_regions:
                    continue

                # mark this region as processed
                processed_regions.add(region_key)

                # extract regions from BAM files
                for bam_file in [bam1, bam2]:
                    bam_basename = os.path.basename(bam_file).replace('.bam', '')
                    output_bam_path = os.path.join(gene_folder, f"{chrom}_{start}_{end}_{bam_basename}.bam")

                    # use samtools to extract the region
                    try:
                        # construct samtools command
                        command = [
                            samtools_path, "view", "-b", bam_file,
                            f"{chrom}:{extended_start}-{extended_end}",
                            "-o", output_bam_path
                        ]
                        subprocess.run(command, check=True)
                        print(f"generated BAM file: {output_bam_path}")
                        # index the bam file
                        command = [
                            samtools_path, "index", output_bam_path
                        ]
                        subprocess.run(command, check=True)
                    except subprocess.CalledProcessError as e:
                        print(f"Error: cannot extract region {chrom}:{extended_start}-{extended_end} from {bam_file}. Error: {e}")


def main():
    """
    main function to parse arguments and run the processing.
    """
    # 定义命令行参数
    parser = argparse.ArgumentParser(description="extract regions from BAM files based on gene folders and region info.")
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="the input file containing gene region information."
    )
    parser.add_argument(
        "-b1", "--bam1",
        required=True,
        help="the path to the first BAM file (hisat.map_unmapremap.s.bam)."
    )
    parser.add_argument(
        "-b2", "--bam2",
        required=True,
        help="the path to the second BAM file (still_unmap_bwa.s.bam)."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="the output directory to create gene folders and region BAM files."
    )
    parser.add_argument(
        "-e", "--extend",
        type=int,
        default=0,
        help="extend the region by this many bases on each side (default: 0)."
    )
    parser.add_argument(
        "-s", "--samtools",
        default="samtools",
        help="the path to the samtools executable (default: samtools, assumed to be in PATH)."
    )

    # parse command line arguments
    args = parser.parse_args()

    # split input file by gene
    split_input_by_gene(args.input, args.output)

    # extract BAM files based on split gene files
    extract_bam_by_gene(args.output, args.bam1, args.bam2, args.extend, args.samtools)


if __name__ == "__main__":
    main()