#!/usr/bin/env python3
"""
说明：
    本脚本扫描各样本/基因目录下的 `seq_map/mapping_metrics.txt` 文件，
    提取 4 个指标：
    
      1. Mapping_rate(%)
      2. Fully_aligned_rate(%)
      3. Clipped_reads
      4. Clipped_reads_rate(%)
    
    每个软件（stringtie、cufflinks、trinity、invas）会输出 4 个合并的 TSV 文件。
    
    目录结构假设如下：
        ../<sample>/eva/<gene>/analysis/<software>/seq_map/mapping_metrics.txt

    输出的 TSV 格式：
        - 第一行：表头（"Gene", 所有样本名）
        - 其余行：各基因的数据，若无数据则填 "na"

用法：
    直接运行：
        python merge_mapping_metrics.py
"""

import glob
import os

# 输出目录
OUTPUT_DIR = "output_sct_map_results"
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# 需要提取的 4 个指标
METRIC_KEYS = [
    "Mapping_rate(%)",
    "Fully_aligned_rate(%)",
    "Clipped_reads",
    "Clipped_reads_rate(%)"
]

# 软件列表
software_list = ["stringtie", "cufflinks", "trinity", "invas"]

# 数据结构初始化
# data = { software: { gene: { sample: { metric_key: value } } } }
data = {soft: {} for soft in software_list}
samples_all = {soft: set() for soft in software_list}
genes_all = {soft: set() for soft in software_list}

# 遍历每个软件，查找 mapping_metrics.txt 并提取指标
for soft in software_list:
    pattern = f"../*/eva/*/analysis/{soft}/seq_map/mapping_metrics.tsv"
    print(pattern)
    for filepath in glob.glob(pattern):
        print(filepath)
        # 解析路径
        parts = os.path.normpath(filepath).split(os.sep)

        if len(parts) < 8:
            continue
        sample = parts[1]
        gene = parts[3]

        # 读取 mapping_metrics.txt 并提取指标
        metrics = {}
        print(filepath)
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("Metric"):
                        continue
                    fields = line.split("\t")
                    if len(fields) < 2:
                        continue
                    key = fields[0]
                    value = fields[1]
                    if key in METRIC_KEYS:
                        metrics[key] = value
        except Exception as e:
            print(f"读取文件 {filepath} 时出错：{e}")
            continue

        # 存储数据
        if gene not in data[soft]:
            data[soft][gene] = {}
        data[soft][gene][sample] = metrics

        samples_all[soft].add(sample)
        genes_all[soft].add(gene)

# 生成合并的 TSV 文件
for soft in software_list:
    sorted_samples = sorted(list(samples_all[soft]))
    sorted_genes = sorted(list(genes_all[soft]))

    for metric in METRIC_KEYS:
        # 处理指标名称，去掉特殊字符，适用于文件名
        metric_clean = metric.replace(" ", "_").replace("(", "").replace(")", "").replace("/", "_")
        filename = f"{soft}_{metric_clean}_merged.tsv"
        outfile = os.path.join(OUTPUT_DIR, filename)

        with open(outfile, 'w', encoding='utf-8') as fw:
            # 写入表头
            header = ["Gene"] + sorted_samples
            fw.write("\t".join(header) + "\n")

            # 写入每个基因的数据
            for gene in sorted_genes:
                row = [gene]
                for sample in sorted_samples:
                    val = data[soft].get(gene, {}).get(sample, {}).get(metric, "na")
                    row.append(val)
                fw.write("\t".join(row) + "\n")

        print(f"{outfile} 已保存。")