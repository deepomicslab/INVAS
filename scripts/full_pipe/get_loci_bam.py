import os
import subprocess
import argparse


def split_input_by_gene(input_file, output_dir):
    """
    将输入文件按照基因名称分割到各自的文件夹中。
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 记录基因名对应的文件句柄
    gene_files = {}

    with open(input_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            # 提取基因名称（第二列）
            columns = line.split('\t')
            if len(columns) < 2:
                print(f"Warning: 无效行格式，跳过: {line}")
                continue

            try:
                gene_info = columns[1]  # 第二列示例：PTPRF:43990857-44089343
                gene_name = gene_info.split(':')[0]  # 提取基因名称：PTPRF
            except ValueError:
                print(f"colums: {columns}")
                print(f"Warning: 无法解析基因名称，跳过: {columns[1]}")
                continue

            # 创建基因文件夹
            gene_folder = os.path.join(output_dir, gene_name)
            if not os.path.exists(gene_folder):
                os.makedirs(gene_folder)

            # 打开基因文件
            gene_file_path = os.path.join(gene_folder, f"gene.txt")
            if gene_name not in gene_files:
                gene_files[gene_name] = open(gene_file_path, 'w')

            # 写入基因文件
            gene_files[gene_name].write(line + '\n')

    # 关闭所有文件句柄
    for gene_file in gene_files.values():
        gene_file.close()


def extract_bam_by_gene(output_dir, bam1, bam2, extend, samtools_path="samtools"):
    """
    根据基因文件夹中的分割文件，从 BAM 文件中提取区域。
    """
    for gene_name in os.listdir(output_dir):
        gene_folder = os.path.join(output_dir, gene_name)
        if not os.path.isdir(gene_folder):
            continue

        gene_file = os.path.join(gene_folder, f"gene.txt")
        if not os.path.exists(gene_file):
            print(f"Warning: 基因文件 {gene_file} 不存在，跳过。")
            continue

        # 使用一个集合记录已经处理过的区域，避免重复提取
        processed_regions = set()

        with open(gene_file, 'r') as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue

                # 提取染色体和区域信息
                columns = line.split('\t')
                if len(columns) < 2:
                    print(f"Warning: 无效行格式，跳过: {line}")
                    continue

                # 第一列提取染色体信息
                try:
                    chrom_region = columns[0]  # 示例：chr1:44059296-44059887
                    chrom = chrom_region.split(':')[0]  # 提取染色体名称（chr1）
                except ValueError:
                    print(f"Warning: 无法解析染色体信息，跳过: {columns[0]}")
                    continue

                # 第二列提取区域信息
                try:
                    gene_info = columns[1]  # 示例：PTPRF:43990857-44089343
                    positions = gene_info.split(':')[1]  # 提取区域范围：43990857-44089343
                    if "," in positions:
                        positions = positions.split(",")[0]
                    start, end = map(int, positions.split('-'))
                except ValueError:
                    print(columns)
                    print(f"Warning: 无法解析区域信息，跳过: {columns[1]}")
                    continue

                # 扩展区域
                extended_start = max(0, start - extend)  # 防止起始位置小于 0
                extended_end = end + extend

                # 唯一标识区域的键（染色体、扩展起点和终点）
                region_key = (chrom, extended_start, extended_end)

                # 如果已经处理过该区域，则跳过
                if region_key in processed_regions:
                    continue

                # 将该区域标记为已处理
                processed_regions.add(region_key)

                # 提取 BAM 文件中对应区域
                for bam_file in [bam1, bam2]:
                    bam_basename = os.path.basename(bam_file).replace('.bam', '')
                    output_bam_path = os.path.join(gene_folder, f"{chrom}_{start}_{end}_{bam_basename}.bam")

                    # 使用 samtools view 提取区域
                    try:
                        # 构造 samtools view 命令
                        command = [
                            samtools_path, "view", "-b", bam_file,
                            f"{chrom}:{extended_start}-{extended_end}",
                            "-o", output_bam_path
                        ]
                        subprocess.run(command, check=True)
                        print(f"成功生成区域 BAM 文件: {output_bam_path}")
                        # index the bam file
                        command = [
                            samtools_path, "index", output_bam_path
                        ]
                        subprocess.run(command, check=True)
                    except subprocess.CalledProcessError as e:
                        print(f"Error: 无法处理 BAM 文件 {bam_file} 的区域 {chrom}:{extended_start}-{extended_end}。错误: {e}")


def main():
    """
    主函数，解析命令行参数并调用函数。
    """
    # 定义命令行参数
    parser = argparse.ArgumentParser(description="根据基因和区域信息，将输入文件分割到基因文件夹并提取 BAM 文件。")
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="输入文件路径，包含基因和区域信息。"
    )
    parser.add_argument(
        "-b1", "--bam1",
        required=True,
        help="第一个 BAM 文件路径（hisat.map_unmapremap.s.bam）。"
    )
    parser.add_argument(
        "-b2", "--bam2",
        required=True,
        help="第二个 BAM 文件路径（still_unmap_bwa.s.bam）。"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="输出目录，生成基因文件夹及区域 BAM 文件。"
    )
    parser.add_argument(
        "-e", "--extend",
        type=int,
        default=0,
        help="区域扩展长度（默认: 0）。"
    )
    parser.add_argument(
        "-s", "--samtools",
        default="samtools",
        help="samtools 可执行文件的路径（默认: samtools，假定已在 PATH 中）。"
    )

    # 解析命令行参数
    args = parser.parse_args()

    # 分割输入文件到基因文件夹
    split_input_by_gene(args.input, args.output)

    # 根据分割的基因文件提取 BAM 文件
    extract_bam_by_gene(args.output, args.bam1, args.bam2, args.extend, args.samtools)


if __name__ == "__main__":
    main()