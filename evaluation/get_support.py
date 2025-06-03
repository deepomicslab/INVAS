import os
import re
from collections import defaultdict
import argparse


def parse_folder(input_folder, min_support):
    """
    解析指定文件夹，统计基因在不同工具中的支持数目，并过滤 inversion 长度 > 10000 的记录。
    同时过滤支持数目小于指定值的基因。
    """
    # 支持的工具文件夹名称
    tools = ["delly", "manta", "svaba", "lumpy"]
    # 用于存储基因和支持它的工具
    gene_support = defaultdict(lambda: {"tools": set(), "regions": []})
    # 用于存储最终的输出行
    output_lines = []

    # 遍历文件夹中的工具子文件夹
    for tool in tools:
        tool_folder = os.path.join(input_folder, tool)
        if not os.path.isdir(tool_folder):
            print(f"Warning: {tool_folder} 文件夹不存在，跳过...")
            continue

        # 遍历工具文件夹中的每个文件
        for filename in os.listdir(tool_folder):
            # 匹配文件名格式：chrxxx_inv_in_gene_dedup_highcof.txt
            if re.match(r'^chr\w+_inv_in_gene_dedup_highcof\.txt$', filename):
                filepath = os.path.join(tool_folder, filename)
                # 解析文件内容
                with open(filepath, 'r') as file:
                    for line in file:
                        line = line.strip()
                        if line:  # 跳过空行
                            # 分割列，假设文件以空格或制表符分隔
                            columns = re.split(r'\s+', line)
                            if len(columns) < 3:
                                print(f"Warning: 无效行格式，跳过: {line}")
                                continue

                            # 检查第一列的 inversion 长度
                            inversion_region = columns[0]  # 第一列示例：chr17:1327339-1332610
                            match = re.match(r'chr\w+:(\d+)-(\d+)', inversion_region)
                            if not match:
                                print(f"Warning: 无效的 inversion 区域格式，跳过: {line}")
                                continue

                            start, end = map(int, match.groups())
                            inversion_length = end - start
                            if inversion_length > 10000:
                                # 跳过长度大于 10000 的记录
                                print(f"过滤掉 inversion 长度 > 10000 的记录: {line}")
                                continue

                            # 获取第二列（基因名称）
                            gene_info = columns[1]  # 示例：CRK:1323982-1366456
                            gene_name = gene_info.split(":")[0]  # 仅保留基因名称（如 CRK）

                            # 获取第三列（基因的不同区域）
                            region_info = columns[2]  # 示例：chr17:1327285-1327419-8.09

                            # 记录支持工具和区域
                            gene_support[gene_name]["tools"].add(tool)
                            gene_support[gene_name]["regions"].append((line, tool))

    # 根据解析的内容，生成输出行
    for gene_name, data in gene_support.items():
        tools_supported = len(data["tools"])
        if tools_supported >= min_support:  # 过滤支持数目小于指定值的基因
            for line, tool in data["regions"]:
                if tool == "manta":  # 仅保留来自 manta 工具的记录
                    # 在原始行后面添加一列支持数目
                    output_lines.append(f"{line}\t{tools_supported}")

    return output_lines


def output_results(output_lines, output_file=None):
    """
    输出结果到控制台或文件。
    """
    if output_file:
        with open(output_file, 'w') as f:
            f.write("\n".join(output_lines))
        print(f"结果已输出到文件: {output_file}")
    else:
        print("基因支持情况：")
        for line in output_lines:
            print(line)


def main():
    """
    主函数，处理命令行参数并调用解析与输出函数。
    """
    # 定义命令行参数
    parser = argparse.ArgumentParser(description="统计基因在不同工具中的支持数目，并过滤 inversion 长度 > 10000 的记录。")
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="输入的主文件夹路径，包含工具子文件夹（delly, manta, svaba, lumpy）。"
    )
    parser.add_argument(
        "-o", "--output",
        help="（可选）输出结果文件路径。如果未指定，则结果将打印到控制台。"
    )
    parser.add_argument(
        "-m", "--min_support",
        type=int,
        required=True,
        default=2,
        help="过滤支持数目小于指定值的基因的最小支持数目。"
    )

    # 解析命令行参数
    args = parser.parse_args()

    # 校验输入文件夹路径
    input_folder = args.input
    if not os.path.isdir(input_folder):
        print(f"错误：输入的文件夹路径无效或不存在：{input_folder}")
        return

    # 调用解析函数
    output_lines = parse_folder(input_folder, args.min_support)

    # 输出结果
    output_results(output_lines, args.output)


if __name__ == "__main__":
    main()

