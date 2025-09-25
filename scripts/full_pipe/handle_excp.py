import argparse

def parse_region_with_depth(region_with_depth):
    """
    parse a region string with depth, e.g., 'chr7:40879375-40880470-20.5'
    """
    *region_parts, depth = region_with_depth.split("-")
    region = "-".join(region_parts)
    return region, float(depth)

def parse_region(region):
    """
    parse a region string, e.g., 'chr7:40879375-40880470'
    """
    chrom, coords = region.split(":")
    start, end = map(int, coords.split("-"))
    return chrom, start, end

def region_length(region):
    """
    calculate the length of a region string, e.g., 'chr7:40879375-40880470'
    """
    _, start, end = parse_region(region)
    return end - start + 1

def process_files(file1, threshold, output_file1, output_file2):
    """
    handle the two input files and update their contents based on the threshold.
    """
    updated_file1 = []
    updated_file2 = []

    # read the first file and process each line
    with open(file1, "r") as f1:
        for line in f1:
            cols = line.strip().split()
            print(cols)
            region1 = cols[0]
            region2_with_depth = cols[2]
            region2, depth = parse_region_with_depth(region2_with_depth)
            print(region1, region2, depth)

            # check the length of region2 against the threshold
            if region_length(region2) < threshold:
                updated_region = region1
                updated_region_with_depth = f"{region1}-{depth}"
            else:
                updated_region = region2
                updated_region_with_depth = region2_with_depth

            # update the first file content
            updated_file1.append(f"{cols[0]}\t{cols[1]}\t{updated_region_with_depth}")

            # update the second file content
            chrom, start, end = parse_region(updated_region)
            updated_file2.append(f"{chrom}\t{start}\t{end}\t{depth}")

    # write the updated first file
    with open(output_file1, "w") as out1:
        out1.write("\n".join(updated_file1) + "\n")

    # write the updated second file
    with open(output_file2, "w") as out2:
        out2.write("\n".join(updated_file2) + "\n")

def main():
    # create the argument parser
    parser = argparse.ArgumentParser(description="handle exceptions in region files based on length threshold.")
    parser.add_argument("file1", help="the first input file path")
    parser.add_argument("output_file1", help="the updated first output file path")
    parser.add_argument("output_file2", help="the updated second output file path")
    parser.add_argument("--threshold", type=int, default=30, help="length threshold (default: 30)")

    # parse command-line arguments
    args = parser.parse_args()

    # 调用处理函数
    process_files(args.file1, args.threshold, args.output_file1, args.output_file2)

if __name__ == "__main__":
    main()