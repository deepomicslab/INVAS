import csv
import argparse


def load_intervals(interval_file):
    """
    加载区间文件，跳过前两行和最后一行。
    参数：
        interval_file: str
            区间文件的路径，格式为 'chrom start end'。
    返回：
        intervals: dict
            按染色体分组的区间列表：
            {
                "chrom": [(start1, end1), (start2, end2), ...]
            }
    """
    intervals = []
    with open(interval_file, 'r') as f:
        lines = f.readlines()

    # 跳过前两行和最后一行
    lines = lines[2:-1]

    for line in lines:
        start, end = line.strip().split()[0].split('-')
        start, end = int(start), int(end)
        intervals.append((start, end))
    return intervals


def split_exon_by_intervals(exon, intervals):
    """
    根据区间分割外显子。
    参数：
        exon: tuple
            外显子信息 (chrom, start, end, strand, attributes)。
        intervals: list
            当前染色体的区间列表 [(start1, end1), (start2, end2), ...]。
    返回：
        split_exons: list
            分割后的外显子列表，每个元素是 (start, end)。
    """
    chrom, exon_start, exon_end, strand, attributes = exon
    split_exons = []
    current_start = exon_start

    # 遍历区间并检查是否需要分割
    for interval_start, interval_end in sorted(intervals):
        if interval_start > exon_end:
            break  # 后续区间不再与外显子重叠
        if interval_end < current_start:
            continue  # 当前区间在外显子之前

        # 区间与外显子的重叠部分
        if interval_start <= exon_end and interval_end >= current_start:
            if current_start < interval_start:
                split_exons.append((current_start, interval_start - 1))  # 前段
            split_exons.append((max(current_start, interval_start), min(exon_end, interval_end)))  # 重叠段
            current_start = interval_end + 1

    # 添加剩余部分
    if current_start <= exon_end:
        split_exons.append((current_start, exon_end))

    return [(chrom, start, end, strand, attributes) for start, end in split_exons]


def process_gtf(gtf_file, intervals, output_file):
    """
    处理 GTF 文件，将外显子按照区间分割。
    参数：
        gtf_file: str
            原始 GTF 文件路径。
        intervals: list
            按染色体分组的区间列表。
        output_file: str
            分割后的 GTF 文件输出路径。
    """
    with open(gtf_file, 'r') as infile, open(output_file, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        for row in reader:
            # 跳过注释行
            if row[0].startswith('#'):
                writer.writerow(row)
                continue

            chrom, source, feature, start, end, score, strand, frame, attributes = row

            # 仅处理外显子信息
            if feature != 'exon':
                writer.writerow(row)
                continue

            start, end = int(start), int(end)
            exon = (chrom, start, end, strand, attributes)


            # 按染色体区间分割外显子
            split_exons = split_exon_by_intervals(exon, intervals)

            # 写入分割后的外显子
            for split_exon in split_exons:
                chrom, split_start, split_end, strand, attributes = split_exon
                writer.writerow([chrom, source, feature, split_start, split_end, score, strand, frame, attributes])


def process_gtf2(gtf_file, intervals, output_file):
    """
    Processes a GTF file, splitting exons based on intervals.
    Parameters:
        gtf_file : str
            Path to the input GTF file.
        intervals : list
            List of intervals grouped by chromosome.
        output_file : str
            Path to write the split GTF file.
    """
    with open(gtf_file, 'r') as infile, open(output_file, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')  # GTF files are tab-delimited
        
        for row in reader:
            # Skip comment lines
            if row[0].startswith('#'):
                outfile.write('\t'.join(row) + '\n')
                continue

            # Parse the GTF fields
            chrom, source, feature, start, end, score, strand, frame, attributes = row

            # Write non-exon features as-is
            if feature != 'exon':
                outfile.write('\t'.join(row) + '\n')
                continue

            # Process exon features
            start, end = int(start), int(end)
            exon = (chrom, start, end, strand, attributes)

            # Split exons based on intervals
            split_exons = split_exon_by_intervals(exon, intervals)

            # Write the split exons
            for split_exon in split_exons:
                chrom, split_start, split_end, strand, attributes = split_exon
                print(split_exon)
                # Write the line using proper GTF formatting
                outfile.write(
                    f"{chrom}\t{source}\t{feature}\t{split_start}\t{split_end}\t{score}\t{strand}\t{frame}\t{attributes}\n"
                )


def main():
    """
    主函数，解析命令行参数并运行处理逻辑。
    """
    parser = argparse.ArgumentParser(description="按照区间分割 GTF 文件中的外显子")
    parser.add_argument("-g", "--gtf", required=True, help="输入 GTF 文件路径")
    parser.add_argument("-i", "--intervals", required=True, help="区间文件路径")
    parser.add_argument("-o", "--output", required=True, help="输出 GTF 文件路径")

    args = parser.parse_args()

    # 加载区间文件
    intervals = load_intervals(args.intervals)

    # 处理 GTF 文件
    process_gtf2(args.gtf, intervals, args.output)


if __name__ == "__main__":
    main()