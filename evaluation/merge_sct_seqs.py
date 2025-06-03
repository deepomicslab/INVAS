#!/usr/bin/env python3
"""
说明：
    本脚本扫描各样本/基因目录下由 eva_seqs_single.py 生成的 summary_report.tsv 文件，
    从中抽取 5 个指标：
    
      1. Reference Transcript Recovery Rate (RTRR)
      2. High Quality Ratio (aligned transcripts)
      3. High Quality Ratio (total assembled transcripts)
      4. Average Query Coverage
      5. Average Target Coverage
      
    假设 summary_report.tsv 文件内容格式类似如下：
    
        Metric                                      Value
        Total assembled transcripts                 13325
        Total reference transcripts                 7
        Assembled transcripts with alignment        7
        Assembled transcripts without alignment      13318
        Reference transcripts hit                   6
        Reference transcripts missing               1
        Reference Transcript Recovery Rate (RTRR)   85.71%
        Average Query Coverage                      91.63%
        Average Target Coverage                     89.76%
        Average Percent Identity                    100.00%
        High quality assembled transcripts          6
        High Quality Ratio (aligned transcripts)    85.71%
        High Quality Ratio (total assembled transcripts) 0.05%
    
    对于每个软件（stringtie、cufflinks、trinity），脚本会将所有样本、所有基因的数据合并，
    然后按照指标分别输出到单个 TSV 文件中。每个文件的第一行是表头，第1列为 "Gene"，其余列为所有样本。
    如果某个样本/基因组合没有生成结果文件或对应指标，则输出 "na"。
    
    此外，输出的所有结果文件将保存到指定的文件夹中（通过 OUTPUT_DIR 变量定义）。
    
用法：
    为脚本添加可执行权限或直接用 Python 运行：
        python merge_summary_one_metric_per_file.py
"""

import glob
import os

# 输出结果将保存到指定文件夹，若该文件夹不存在则自动创建
OUTPUT_DIR = "output_sct_seq_results"
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# 需要抽取的 5 个指标
METRIC_KEYS = [
    "Reference Transcript Recovery Rate (RTRR)",
    "High Quality Ratio (aligned transcripts)",
    "High Quality Ratio (total assembled transcripts)",
    "Average Query Coverage",
    "Average Target Coverage"
]

# 三种软件名称
software_list = ["stringtie", "cufflinks", "trinity","invas"]

# 建立数据结构：按软件存储数据
# data = { software: { gene: { sample: { metric_key: value } } } }
data = {soft: {} for soft in software_list}
# 分别收集每个软件中所有出现的 sample 和 gene
samples_all = {soft: set() for soft in software_list}
genes_all = {soft: set() for soft in software_list}

# 根据固定目录结构查找 summary_report.tsv 文件，并抽取指标
# 示例路径： ../<sample_name>/eva/<gene>/analysis/<software>/seq_eva/summary_report.tsv
for soft in software_list:
    pattern = f"../*/eva/*/analysis/{soft}/seq_eva/summary_report.tsv"
    for filepath in glob.glob(pattern):
        # 例如 filepath: "../Sample01/eva/geneA/analysis/stringtie/seq_eva/summary_report.tsv"
        parts = os.path.normpath(filepath).split(os.sep)
        # 期望 parts 至少为：
        # ["..", "<sample_name>", "eva", "<gene>", "analysis", "<software>", "seq_eva", "summary_report.tsv"]
        if len(parts) < 8:
            continue
        sample = parts[1]
        gene = parts[3]

        # 读取 summary_report.tsv 中的指标数据
        metrics = {}
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    # 跳过表头
                    if line.startswith("Metric"):
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

        # 将读取到的指标保存到数据结构中
        if gene not in data[soft]:
            data[soft][gene] = {}
        data[soft][gene][sample] = metrics

        samples_all[soft].add(sample)
        genes_all[soft].add(gene)

# 针对每个软件和每个指标生成一个合并的表格文件，并保存到 OUTPUT_DIR
for soft in software_list:
    sorted_samples = sorted(list(samples_all[soft]))
    sorted_genes = sorted(list(genes_all[soft]))

    for metric in METRIC_KEYS:
        # 对指标名称进行简单清理，用下划线替换空格、括号、"/"等字符以生成文件名
        metric_clean = metric.replace(" ", "_").replace("(", "").replace(")", "").replace("/", "_")
        filename = f"{soft}_{metric_clean}_merged.tsv"
        outfile = os.path.join(OUTPUT_DIR, filename)
        with open(outfile, 'w', encoding='utf-8') as fw:
            # 第一行：表头，第1列为 "Gene"，后续为所有样本名称
            header = ["Gene"] + sorted_samples
            fw.write("\t".join(header) + "\n")
            # 输出每个基因对应的数据行
            for gene in sorted_genes:
                row = [gene]
                for sample in sorted_samples:
                    # 如果当前样本/基因组合没有数据或对应指标，则填 "na"
                    val = data[soft].get(gene, {}).get(sample, {}).get(metric, "na")
                    row.append(val)
                fw.write("\t".join(row) + "\n")
        print(f"{outfile} 已保存。")