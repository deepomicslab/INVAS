import os
import re
from collections import defaultdict
import argparse


def parse_folder(input_folder, min_support):
    """
    解析指定文件夹，统计基因在不同工具中的支持数目，并过滤 inversion 长度 > 10000 的记录。
    同时过滤支持数目小于指定值的基因。
    """
    # supported tools
    tools = ["delly", "manta", "svaba", "lumpy"]
    # for storing gene support information
    gene_support = defaultdict(lambda: {"tools": set(), "regions": []})
    # for storing output lines
    output_lines = []

    # loop through each tool folder
    for tool in tools:
        tool_folder = os.path.join(input_folder, tool)
        if not os.path.isdir(tool_folder):
            print(f"Warning: {tool_folder} 文件夹不存在，跳过...")
            continue

        # loop through each file in the tool folder
        for filename in os.listdir(tool_folder):
            # match the specific file
            if re.match(r'^chr\w+_inv_in_gene_dedup_highcof\.txt$', filename):
                filepath = os.path.join(tool_folder, filename)
                print(f"解析文件: {filepath}")
                # parse the file
                with open(filepath, 'r') as file:
                    for line in file:
                        line = line.strip()
                        if line:
                            # split columns, assuming space or tab separation
                            columns = re.split(r'\s+', line)
                            if len(columns) < 3:
                                print(f"Warning: invalid line format, skipping: {line}")
                                continue

                            # check the inversion length in the first column
                            inversion_region = columns[0]
                            match = re.match(r'chr\w+:(\d+)-(\d+)', inversion_region)
                            if not match:
                                print(f"Warning: invalid inversion region format, skipping: {line}")
                                continue

                            start, end = map(int, match.groups())
                            inversion_length = end - start
                            if inversion_length > 20000:
                                print(f"filtering out inversion length > 20000: {line}")
                                continue

                            # get the gene name from the second column
                            gene_info = columns[1] 
                            gene_name = gene_info.split(":")[0]

                            # get the region info from the third column
                            region_info = columns[2]

                            # record supported tools and regions
                            gene_support[gene_name]["tools"].add(tool)
                            gene_support[gene_name]["regions"].append((line, tool))

    # generate output lines
    for gene_name, data in gene_support.items():
        tools_supported = len(data["tools"])
        if tools_supported >= min_support: 
            for line, tool in data["regions"]:
                # append the tool name to each line
                output_lines.append(f"{line}\t{tools_supported}")

    return output_lines


def output_results(output_lines, output_file=None):
    """
    output the results to a file or print to console.
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
    main function to handle argument parsing and function calls.
    """
    # 定义命令行参数
    parser = argparse.ArgumentParser(description="stat gene support from different tools and filter by min support.")
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="the input folder containing subfolders for each tool (delly, manta, svaba, lumpy)."
    )
    parser.add_argument(
        "-o", "--output",
        help="(optional) output file path. If not specified, results will be printed to console."
    )
    parser.add_argument(
        "-m", "--min_support",
        type=int,
        required=True,
        default=2,
        help="the minimum support count for filtering genes with support count less than the specified value."
    )

    # parse command line arguments
    args = parser.parse_args()

    # validate input folder path
    input_folder = args.input
    if not os.path.isdir(input_folder):
        print(f"Error: Invalid or non-existent input folder path: {input_folder}")
        return

    # call the parsing function
    output_lines = parse_folder(input_folder, args.min_support)

    # output results
    output_results(output_lines, args.output)


if __name__ == "__main__":
    main()
