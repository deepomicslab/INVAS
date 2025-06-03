import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import ptitprince as pt  
import re  # 用于解析样本名
from scipy.stats import iqr

# 📌 颜色映射（可自定义不同软件的颜色）
software_colors = {
    "stringtie": "#1f77b4",
    "cufflinks": "#ff7f0e",
    "trinity": "#2ca02c",
    "invas": "#d62728"
}

# 📌 目标软件列表
SOFTWARES = ["stringtie", "cufflinks", "trinity", "invas"]

# 📌 读取 `*_pearson_matrix.tsv` 和 `*_rmse_matrix.tsv` 数据
data = []
for software in SOFTWARES:
    pearson_file = f"{software}_pearson_matrix.tsv"
    rmse_file = f"{software}_rmse_matrix.tsv"

    if not os.path.exists(pearson_file) or not os.path.exists(rmse_file):
        print(f"⚠️ 文件缺失: {pearson_file} 或 {rmse_file}，跳过 {software}")
        continue  

    print(f"✅ 读取文件: {pearson_file} 和 {rmse_file}")

    pearson_df = pd.read_csv(pearson_file, sep="\t")
    rmse_df = pd.read_csv(rmse_file, sep="\t")

    # 解析样本名，提取 WGS 深度 & RNA-Seq 深度
    sample_cols = pearson_df.columns[1:]  
    sample_info = []
    for sample in sample_cols:
        match = re.match(r"sim_(\d+)_(\d+)", sample)
        if match:
            wgs_depth, rna_depth = map(int, match.groups())
            sample_info.append((sample, wgs_depth, rna_depth))

    # 计算基因的 Pearson 相关性和 RMSE
    for gene_id in pearson_df["Gene_ID"]:
        for sample, wgs_depth, rna_depth in sample_info:
            avg_pearson = pearson_df.loc[pearson_df["Gene_ID"] == gene_id, sample].values[0]
            avg_rmse = rmse_df.loc[rmse_df["Gene_ID"] == gene_id, sample].values[0]

            data.append({
                "gene_id": gene_id,
                "software": software,
                "pearson": avg_pearson,
                "rmse": avg_rmse,
                "wgs_depth": wgs_depth,
                "rna_depth": rna_depth
            })

# 📌 转换为 DataFrame
df = pd.DataFrame(data).replace([np.inf, -np.inf], np.nan).dropna()

# 📌 去除异常值（IQR 方法）
def remove_outliers(df, column):
    q1, q3 = df[column].quantile([0.25, 0.75])
    iqr_value = iqr(df[column])
    lower_bound, upper_bound = q1 - 1.5 * iqr_value, q3 + 1.5 * iqr_value
    return df[(df[column] >= lower_bound) & (df[column] <= upper_bound)]

df_cleaned = remove_outliers(df, "pearson")
df_cleaned = remove_outliers(df_cleaned, "rmse")

# 📌 计算 Median（中位数）
median_values = df_cleaned.groupby("software")[["pearson", "rmse"]].median()
print("\n📊 各软件的 Pearson & RMSE 中位数 (Median)：")
print(median_values)

# 🚀 **Raincloud 总图**
def plot_raincloud_general(metric, df):
    """ 绘制 Pearson 和 RMSE 的 Raincloud 总图（不分深度） """
    plt.figure(figsize=(8, 5))
    pt.RainCloud(
        x="software", y=metric, data=df,
        palette=[software_colors[sw] for sw in df["software"].unique()],
        width_viol=0.7, orient="v", move=0.2, alpha=0.5
    )
    plt.xlabel("Software")
    plt.ylabel(metric)
    plt.title(f"{metric} Distribution Across All Depths")
    plt.xticks(rotation=60, fontweight="bold")

    file_name = f"raincloud_{metric}_all_no_outliers.png"
    plt.savefig(file_name, dpi=300, bbox_inches="tight")
    plt.show()

plot_raincloud_general("pearson", df_cleaned)
plot_raincloud_general("rmse", df_cleaned)

# 🚀 **随着 WGS & RNA-Seq 深度变化的趋势图**
def plot_trend(metric, df, depth_type):
    """ 绘制 Pearson & RMSE 随 WGS 或 RNA-Seq 深度变化的趋势图 """
    plt.figure(figsize=(8, 5))
    sns.lineplot(
        data=df, x=depth_type, y=metric, hue="software",
        marker="o", palette=software_colors
    )
    plt.xlabel(f"{depth_type.replace('_', ' ').title()}")
    plt.ylabel(metric)
    plt.title(f"{metric} vs {depth_type.replace('_', ' ').title()}")
    plt.grid(True)

    file_name = f"trend_{metric}_{depth_type}_no_outliers.png"
    plt.savefig(file_name, dpi=300, bbox_inches="tight")
    plt.show()

plot_trend("pearson", df_cleaned, "wgs_depth")
plot_trend("pearson", df_cleaned, "rna_depth")
plot_trend("rmse", df_cleaned, "wgs_depth")
plot_trend("rmse", df_cleaned, "rna_depth")