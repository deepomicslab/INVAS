import os
import re
from collections import defaultdict
import argparse


def parse_folder(input_folder, min_support):
    """
    parse the specified folder, count gene support across different tools,
    """
    # Supported tool folder names
    tools = ["delly", "manta", "svaba", "lumpy"]
    # store gene support information
    gene_support = defaultdict(lambda: {"tools": set(), "regions": []})
    # store final output lines
    output_lines = []

    # iterate over each tool folder
    for tool in tools:
        tool_folder = os.path.join(input_folder, tool)
        if not os.path.isdir(tool_folder):
            print(f"Warning: {tool_folder} folder does not exist, skipping...")
            continue

        # iterate over each file in the tool folder
        for filename in os.listdir(tool_folder):
            # match file name format: chrxxx_inv_in_gene_dedup_highcof.txt
            if re.match(r'^chr\w+_inv_in_gene_dedup_highcof\.txt$', filename):
                filepath = os.path.join(tool_folder, filename)
                # parse file content
                with open(filepath, 'r') as file:
                    for line in file:
                        line = line.strip()
                        if line:
                            # split columns, assuming the file is space or tab delimited
                            columns = re.split(r'\s+', line)
                            if len(columns) < 3:
                                print(f"Warning: invalid line format, skipping: {line}")
                                continue

                            # check the inversion length in the first column
                            inversion_region = columns[0]  # example: chr17:1323982-1366456
                            match = re.match(r'chr\w+:(\d+)-(\d+)', inversion_region)
                            if not match:
                                print(f"Warning: invalid inversion region format, skipping: {line}")
                                continue

                            start, end = map(int, match.groups())
                            inversion_length = end - start
                            if inversion_length > 10000:
                                # skip inversion length > 10000
                                print(f"Filtering out inversion length > 10000: {line}")
                                continue

                            
                            gene_info = columns[1]
                            gene_name = gene_info.split(":")[0]

                            region_info = columns[2]

                            gene_support[gene_name]["tools"].add(tool)
                            gene_support[gene_name]["regions"].append((line, tool))

    for gene_name, data in gene_support.items():
        tools_supported = len(data["tools"])
        if tools_supported >= min_support:
            for line, tool in data["regions"]:
                if tool == "manta":
                    output_lines.append(f"{line}\t{tools_supported}")

    return output_lines


def output_results(output_lines, output_file=None):
    """
    output the results to a specified file or print to console.
    """
    if output_file:
        with open(output_file, 'w') as f:
            f.write("\n".join(output_lines))
        print(f"Results have been written to file: {output_file}")
    else:
        print("Gene support information:")
        for line in output_lines:
            print(line)


def main():
    """
    main function to handle command-line arguments and execute the parsing and output functions.
    """
    # define argument parser
    parser = argparse.ArgumentParser(description="stat gene support across different tools in the specified folder.")
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input main folder path, containing tool subfolders (delly, manta, svaba, lumpy)."
    )
    parser.add_argument(
        "-o", "--output",
        help="(optional) Output file path for results. If not specified, results will be printed to console."
    )
    parser.add_argument(
        "-m", "--min_support",
        type=int,
        required=True,
        default=2,
        help="filter genes supported by at least this many tools (default: 2)."
    )

    # parse command-line arguments
    args = parser.parse_args()

    # validate input folder path
    input_folder = args.input
    if not os.path.isdir(input_folder):
        print(f"error: input folder {input_folder} does not exist.")
        return

    # call the parsing function
    output_lines = parse_folder(input_folder, args.min_support)

    # output results
    output_results(output_lines, args.output)


if __name__ == "__main__":
    main()