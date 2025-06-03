import os
import csv
import pandas as pd
from collections import defaultdict

# 软件和指标
softwares = ["stringtie", "cufflinks", "trinity", "invas"]
metrics = ["score", "p_good_mapping", "p_fragments_mapped", "mean_orf_percent"]
metrics = [ "optimal_score"]

# 根路径和输入文件
sample_root = "/data/wxd/INVAS_SIM/batch_sim"
sample_file = "../../sample_dirs.txt"  # 每行是：sim_30_150	sim_30_150	sim_30_150

# 输出目录
output_dir = "/data/wxd/INVAS_SIM/batch_sim/eva/output_trasrate_results"
os.makedirs(output_dir, exist_ok=True)

# 数据结构： metrics_data[metric][software][gene_id][sample_name] = value
metrics_data = {
    metric: {software: defaultdict(dict) for software in softwares}
    for metric in metrics
}

all_samples = set()
all_genes = set()

# 读取样本名列表（只取每行第一个字段）
with open(sample_file) as f:
    for line in f:
        sample_name = line.strip().split('\t')[0]
        all_samples.add(sample_name)
        sample_path = os.path.join(sample_root, sample_name)

        # 建立 gene_id -> gene_name 映射
        gene_id_to_name = {}
        truth_file = os.path.join(sample_path, "truth", "inversions.tsv")
        if not os.path.exists(truth_file):
            print(f"⚠️ 跳过 {sample_name}，找不到 {truth_file}")
            continue

        with open(truth_file) as tf:
            reader = csv.DictReader(tf, delimiter="\t")
            for row in reader:
                gene_id = row["gene_id"]
                gene_name = row["gene_name"]
                gene_id_to_name[gene_id] = gene_name

        gene_list_file = os.path.join(sample_path, "truth", "inversions.tsv.unique")
        if not os.path.exists(gene_list_file):
            print(f"⚠️ 跳过 {sample_name}，找不到 {gene_list_file}")
            continue

        with open(gene_list_file) as gene_file:
            for gene_id in gene_file:
                gene_id = gene_id.strip()
                all_genes.add(gene_id)

                for software in softwares:
                    if software == "invas":
                        gene_name = gene_id_to_name.get(gene_id)
                        if not gene_name:
                            print(f"⚠️ {sample_name}: gene_id {gene_id} 无法映射 gene_name，跳过 invas")
                            continue
                        transrate_dir = os.path.join(sample_path, "haps", gene_name, "haps2", "transrate_output")
                    else:
                        transrate_dir = os.path.join(sample_path, "rna", "genes", gene_id, "analysis", software, "transrate_output")

                    csv_file = os.path.join(transrate_dir, "transrate_reads", "assemblies.csv")
                    if not os.path.exists(csv_file):
                        continue

                    with open(csv_file) as csvf:
                        reader = csv.DictReader(csvf)
                        for row in reader:
                            for metric in metrics:
                                value = row.get(metric, "NA")
                                metrics_data[metric][software][gene_id][sample_name] = value

# 构建每个 software + metric 的矩阵输出
all_samples = sorted(all_samples)
all_genes = sorted(all_genes)

for metric in metrics:
    for software in softwares:
        gene_to_sample = metrics_data[metric][software]

        df_rows = []
        for gene in all_genes:
            row = [gene]
            for sample in all_samples:
                value = gene_to_sample.get(gene, {}).get(sample, "NA")
                row.append(value)
            df_rows.append(row)

        df = pd.DataFrame(df_rows, columns=["gene_id"] + all_samples)

        filename = f"{metric}_{software}.tsv"
        output_path = os.path.join(output_dir, filename)
        df.to_csv(output_path, sep="\t", index=False)
        print(f"✅ 保存文件：{output_path}")