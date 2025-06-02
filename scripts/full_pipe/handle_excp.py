import argparse

def parse_region_with_depth(region_with_depth):
    """
    解析第三列包含深度的区域，例如 'chr7:40880469-40880469-5.0'
    返回 (region, depth)，如 ('chr7:40880469-40880469', 5.0)
    """
    *region_parts, depth = region_with_depth.split("-")
    region = "-".join(region_parts)
    return region, float(depth)

def parse_region(region):
    """
    解析区域字符串，例如 'chr7:40879375-40880470'
    返回 (chrom, start, end)
    """
    chrom, coords = region.split(":")
    start, end = map(int, coords.split("-"))
    return chrom, start, end

def region_length(region):
    """
    计算区域长度
    """
    _, start, end = parse_region(region)
    return end - start + 1

def process_files(file1, threshold, output_file1, output_file2):
    """
    处理第一个输入文件并生成更新后的文件
    :param file1: 第一个输入文件路径
    :param threshold: 长度阈值
    :param output_file1: 第一个输出文件路径
    :param output_file2: 第二个输出文件路径
    """
    updated_file1 = []
    updated_file2 = []

    # 读取第一个文件并更新内容
    with open(file1, "r") as f1:
        for line in f1:
            cols = line.strip().split()
            print(cols)
            region1 = cols[0]  # 第一列的区域
            region2_with_depth = cols[2]  # 第三列的区域和深度
            region2, depth = parse_region_with_depth(region2_with_depth)  # 拆分出区域和深度
            print(region1, region2, depth)

            # 判断第三列区域长度是否小于阈值
            if region_length(region2) < threshold:
                updated_region = region1  # 替换为第一列的区域
                updated_region_with_depth = f"{region1}-{depth}"  # 保留深度
            else:
                updated_region = region2  # 保留原区域
                updated_region_with_depth = region2_with_depth  # 保留原区域和深度

            # 更新第一文件内容
            updated_file1.append(f"{cols[0]}\t{cols[1]}\t{updated_region_with_depth}")

            # 更新第二文件内容
            chrom, start, end = parse_region(updated_region)
            updated_file2.append(f"{chrom}\t{start}\t{end}\t{depth}")

    # 写入更新后的第一个文件
    with open(output_file1, "w") as out1:
        out1.write("\n".join(updated_file1) + "\n")

    # 写入更新后的第二个文件
    with open(output_file2, "w") as out2:
        out2.write("\n".join(updated_file2) + "\n")

def main():
    # 创建解析器
    parser = argparse.ArgumentParser(description="处理文件并更新区域信息")
    parser.add_argument("file1", help="第一个输入文件路径")
    parser.add_argument("output_file1", help="更新后的第一个输出文件路径")
    parser.add_argument("output_file2", help="更新后的第二个输出文件路径")
    parser.add_argument("--threshold", type=int, default=30, help="长度阈值（默认值为100）")

    # 解析命令行参数
    args = parser.parse_args()

    # 调用处理函数
    process_files(args.file1, args.threshold, args.output_file1, args.output_file2)

if __name__ == "__main__":
    main()