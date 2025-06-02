# filter_overlapping_regions_by_ratio.py
import argparse

def parse_coordinates(coordinate_string):
    """解析坐标字符串，返回染色体名和起止位置"""
    chrom, positions = coordinate_string.split(':')
    start, end = map(int, positions.split('-'))
    return chrom, start, end

def calculate_overlap_ratio(start_a, end_a, start_b, end_b):
    """计算两个区域的重叠比率"""
    overlap = max(0, min(end_a, end_b) - max(start_a, start_b))
    length_a = end_a - start_a
    length_b = end_b - start_b
    # 计算重叠区域占两个区域长度的最小值的比例
    return overlap / min(length_a, length_b)

def filter_overlapping_regions(input_file, output_file, overlap_threshold):
    processed_regions = []

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            region = line.strip().split('\t')[0]
            chrom, start, end = parse_coordinates(region)

            # 检查当前区域是否与已处理的区域重叠
            if not any(calculate_overlap_ratio(start, end, prev_start, prev_end) > overlap_threshold for _, prev_start, prev_end in processed_regions):
                # 如果没有重叠，写入到输出文件并记录这个区域
                outfile.write(line)
                processed_regions.append((chrom, start, end))

def main():
    parser = argparse.ArgumentParser(description="Filter out regions with overlap ratio above a threshold.")
    parser.add_argument("input_file", help="Input file path containing regions")
    parser.add_argument("output_file", help="Output file path for the filtered regions")
    parser.add_argument("-t", "--threshold", type=float, default=0.5, help="Overlap threshold to filter regions (default: 0.5)")
    
    args = parser.parse_args()

    filter_overlapping_regions(args.input_file, args.output_file, args.threshold)

if __name__ == "__main__":
    main()