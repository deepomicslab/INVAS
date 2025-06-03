import os
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="合并各个样本文件夹中生成的 TSV 结果，并在合并后的文件中添加样本名列。"
    )
    parser.add_argument(
        "--samples", required=True,
        help="包含样本名的文本文件（每行一个样本名）。"
    )
    parser.add_argument(
        "--indir", default=".",
        help="样本文件夹所在的父目录（默认为当前目录）。"
    )
    parser.add_argument(
        "--result_file", default="result.tsv",
        help="每个样本文件夹中TSV结果文件的名称（默认为 result.tsv）。"
    )
    parser.add_argument(
        "--output", required=True,
        help="合并后的输出TSV文件名。"
    )
    args = parser.parse_args()

    # 读取样本名列表
    with open(args.samples, 'r', encoding='utf-8') as f:
        samples = [line.strip().split("\t")[0] for line in f if line.strip()]

    merged_list = []
    for sample in samples:
        sample_dir = os.path.join(args.indir, sample)
        tsv_file = os.path.join(sample_dir, args.result_file)
        if not os.path.exists(tsv_file):
            print(f"警告：{tsv_file} 不存在，跳过样本 {sample}。")
            continue
        
        # 尝试读取 TSV 文件，默认首行是 header，分隔符为制表符
        try:
            df = pd.read_csv(tsv_file, sep="\t")
        except Exception as e:
            print(f"错误读取 {tsv_file}：{e}")
            continue
        
        # 在第一列插入样本名
        df.insert(0, "Sample", sample)
        merged_list.append(df)

    if merged_list:
        merged_df = pd.concat(merged_list, ignore_index=True)
        merged_df.to_csv(args.output, sep="\t", index=False)
        print(f"合并后的结果已保存到 {args.output}")
    else:
        print("未找到任何TSV文件，退出。")

if __name__ == "__main__":
    main()