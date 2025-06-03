import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob

# 📌 软件配色和顺序
software_colors = {
    "StringTie": "#A13670",
    "Cufflinks": "#FFA45E",
    "Trinity": "#E0628A",
    "Invas": "#002479"
}
software_order = list(software_colors.keys())

# 📌 软件名称映射
software_name_map = {
    "stringtie": "StringTie",
    "cufflinks": "Cufflinks",
    "trinity": "Trinity",
    "invas": "Invas"
}

# 📌 图形风格设置


# 📁 输入输出路径
input_dir = "./"
output_dir = "figures2"
os.makedirs(output_dir, exist_ok=True)

# 📂 加载所有 .tsv 文件
records = []
all_files = glob(os.path.join(input_dir, "*.tsv"))

for file in all_files:
    basename = os.path.basename(file).replace(".tsv", "")
    parts = basename.split("_")
    if len(parts) < 2:
        continue
    software_raw = parts[-1].lower()
    software = software_name_map.get(software_raw, software_raw)
    metric = "_".join(parts[:-1])

    try:
        df = pd.read_csv(file, sep="\t")
    except Exception as e:
        print(f"❌ 无法读取文件: {file} - {e}")
        continue

    for _, row in df.iterrows():
        gene_id = row["gene_id"]
        for sample in row.index[1:]:
            value = row[sample]
            if value == "NA":
                continue
            try:
                val = float(value)
            except:
                continue
            match = re.match(r"sim_(\d+)_(\d+)", sample)
            if not match:
                continue
            wgs_depth = int(match.group(1))
            rna_depth = int(match.group(2))
            records.append({
                "gene_id": gene_id,
                "sample": sample,
                "software": software,
                "metric": metric,
                "value": val,
                "wgs_depth": wgs_depth,
                "rna_depth": rna_depth
            })

# 📌 构建 DataFrame
long_df = pd.DataFrame(records)
long_df.to_csv("transrate_long_table.csv", index=False)
print("✅ 已保存 transrate_long_table.csv")

# 📌 获取指标和深度列表
metrics = long_df["metric"].unique()
wgs_depths = sorted(long_df["wgs_depth"].unique())
rna_depths = sorted(long_df["rna_depth"].unique())

# 📦 所有样本箱线图
for metric in metrics:
    df_metric = long_df[long_df["metric"] == metric]
    if df_metric.empty:
        continue

    plt.figure(figsize=(3, 4))
    sns.boxplot(
        data=df_metric,
        x="software", y="value",
        palette=software_colors,
        order=software_order,
        showfliers=False,
        width=0.3
    )
    plt.title(f"", fontsize=10, loc="left")
    plt.ylabel("Value", fontweight="bold")
    plt.xlabel("")
    plt.xticks(rotation=45, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{metric}_all_samples_box.pdf"), dpi=300)
    plt.savefig(os.path.join(output_dir, f"{metric}_all_samples_box.png"), dpi=300)
    plt.close()
    print(f"📦 保存箱线图: {metric}_all_samples_box")

# 🧩 每个指标拼图 (35 子图)
for metric in metrics:
    df_metric = long_df[long_df["metric"] == metric]
    if df_metric.empty:
        continue

    fig, axes = plt.subplots(
        nrows=7, ncols=5,
        figsize=(8, 12), dpi=300  # A4 横向 + 适当高度
    )
    plt.subplots_adjust(hspace=0.1, wspace=0.001)

    for i, rna in enumerate(rna_depths):
        for j, wgs in enumerate(wgs_depths):
            ax = axes[i, j]
            subset = df_metric[
                (df_metric["rna_depth"] == rna) &
                (df_metric["wgs_depth"] == wgs)
            ]
            if subset.empty:
                ax.axis("off")
                continue
            sns.boxplot(
                data=subset,
                x="software", y="value",
                palette=software_colors,
                order=software_order,
                showfliers=False,
                width=0.3,
                ax=ax
            )
            ax.set_title(f"RNA: {rna}, WGS: {wgs}", fontsize=7)
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.tick_params(axis='x', labelrotation=45, labelsize=6)
            ax.tick_params(axis='y', labelsize=6)
            # set xy tick bold
            ax.xaxis.label.set_fontweight('bold')
            ax.yaxis.label.set_fontweight('bold')

    fig.suptitle(f"", fontsize=12, fontweight="bold")
    out_pdf = os.path.join(output_dir, f"{metric}_subplots_a4.pdf")
    out_png = os.path.join(output_dir, f"{metric}_subplots_a4.png")
    fig.tight_layout()
    fig.savefig(out_pdf)
    fig.savefig(out_png)
    plt.close(fig)
    print(f"🧩 拼图保存完毕: {metric}_subplots_a4.pdf")