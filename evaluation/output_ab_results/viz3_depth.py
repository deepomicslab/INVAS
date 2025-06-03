import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re  
from scipy.stats import iqr

# 📌 颜色映射
software_colors = {
    "StringTie": "#A13670",
    "Cufflinks": "#FFA45E",
    "Trinity": "#E0628A",
    "Invas": "#002479"
}

# 📌 目标软件列表
SOFTWARES = ["StringTie", "Cufflinks", "Trinity", "Invas"]

# 📌 读取数据
data = []
wgs_depths = set()
rna_depths = set()

for software in SOFTWARES:
    pearson_file = f"{software.lower()}_pearson_matrix.tsv"
    rmse_file = f"{software.lower()}_rmse_matrix.tsv"

    if not os.path.exists(pearson_file) or not os.path.exists(rmse_file):
        print(f"⚠️ 文件缺失: {pearson_file} 或 {rmse_file}，跳过 {software}")
        continue  

    print(f"✅ 读取文件: {pearson_file} 和 {rmse_file}")

    pearson_df = pd.read_csv(pearson_file, sep="\t")
    rmse_df = pd.read_csv(rmse_file, sep="\t")

    # 解析样本名
    sample_cols = pearson_df.columns[1:]  
    for sample in sample_cols:
        match = re.match(r"sim_(\d+)_(\d+)", sample)
        if match:
            wgs_depth, rna_depth = map(int, match.groups())
            wgs_depths.add(wgs_depth)
            rna_depths.add(rna_depth)

            for gene_id in pearson_df["Gene_ID"]:
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
    iqr_value = q3 - q1
    lower_bound, upper_bound = q1 - 1.5 * iqr_value, q3 + 1.5 * iqr_value
    return df[(df[column] >= lower_bound) & (df[column] <= upper_bound)]

df_cleaned = remove_outliers(df, "pearson")
df_cleaned = remove_outliers(df_cleaned, "rmse")

# 📌 创建输出目录
output_dir = "optimized_plots"
os.makedirs(output_dir, exist_ok=True)

# 📌 绘图函数（**去除所有图例** & X 轴旋转 45°）
def plot_graph(df, title_suffix, file_suffix, subdir=""):
    subdir_path = os.path.join(output_dir, subdir)
    os.makedirs(subdir_path, exist_ok=True)

    # Pearson vs RMSE 散点图
    plt.figure(figsize=(3, 2))
    sns.scatterplot(
        data=df, x="rmse", y="pearson", hue="software", style="software",
        s=30, edgecolor="black", alpha=0.5, palette=software_colors, legend=None  # 📌 移除图例
    )
    sns.regplot(data=df, x="rmse", y="pearson", scatter=False, ci=None, line_kws={"color": "black", "linestyle": "--"})
    plt.xlabel("RMSE", fontweight='bold', fontsize=8)
    plt.ylabel("Pearson Correlation", fontweight='bold', fontsize=8)
    plt.xticks(rotation=45)  # 📌 旋转 X 轴标签
    # plt.title(f"Pearson Correlation vs RMSE")
    plt.tight_layout()
    plt.savefig(os.path.join(subdir_path, f"pearson_vs_rmse_{file_suffix}.png"), dpi=300, bbox_inches="tight")
    # save pdf
    plt.savefig(os.path.join(subdir_path, f"pearson_vs_rmse_{file_suffix}.pdf"), dpi=300, bbox_inches="tight")
    plt.show()

    # Pearson vs Software 窄箱线图
    plt.figure(figsize=(3, 4))
    sns.boxplot(data=df, x="software", y="pearson", palette=software_colors, showfliers=False, width=0.3)
    plt.ylabel("Pearson Correlation", fontweight='bold')
    # plt.xlabel("Software")
    plt.xticks(rotation=45, fontweight='bold')  # 📌 旋转 X 轴标签
    # plt.title(f"Pearson Correlation by Software {title_suffix}")
    plt.xlabel("")
    
    # 📌 **去除图例**
    plt.legend([], [], frameon=False)

    plt.tight_layout()
    plt.savefig(os.path.join(subdir_path, f"pearson_vs_software_{file_suffix}.png"), dpi=600, bbox_inches="tight")
    # save pdf
    plt.savefig(os.path.join(subdir_path, f"pearson_vs_software_{file_suffix}.pdf"), dpi=300, bbox_inches="tight")

    plt.show()

    # RMSE vs Software 窄箱线图
    plt.figure(figsize=(3, 4))
    sns.boxplot(data=df, x="software", y="rmse", palette=software_colors, showfliers=False, width=0.3)
    # y label bold
    plt.ylabel("RMSE", fontweight='bold')
    # plt.xlabel("Software")
    plt.xlabel("")
    plt.xticks(rotation=45, fontweight='bold')  # 📌 旋转 X 轴标签
    # plt.title(f"RMSE by Software {title_suffix}")
    
    # 📌 **去除图例**
    plt.legend([], [], frameon=False)

    plt.tight_layout()
    plt.savefig(os.path.join(subdir_path, f"rmse_vs_software_{file_suffix}.png"), dpi=600, bbox_inches="tight")
    plt.savefig(os.path.join(subdir_path, f"rmse_vs_software_{file_suffix}.pdf"), dpi=300, bbox_inches="tight")
    
    plt.show()

# 📌 生成全局图（原始数据 & 去除异常值）
plot_graph(df, "(Original Data)", "original")
plot_graph(df_cleaned, "(No Outliers)", "cleaned")

# 📌 遍历 WGS & RNA 深度，分别绘制图
for wgs in sorted(wgs_depths):
    for rna in sorted(rna_depths):
        df_subset = df[(df["wgs_depth"] == wgs) & (df["rna_depth"] == rna)]
        df_subset_cleaned = df_cleaned[(df_cleaned["wgs_depth"] == wgs) & (df_cleaned["rna_depth"] == rna)]

        if df_subset.empty:
            continue

        plot_graph(df_subset, f"(WGS: {wgs}, RNA: {rna})", f"wgs{wgs}_rna{rna}_original", subdir=f"wgs_{wgs}_rna_{rna}")
        plot_graph(df_subset_cleaned, f"(WGS: {wgs}, RNA: {rna}, No Outliers)", f"wgs{wgs}_rna{rna}_cleaned", subdir=f"wgs_{wgs}_rna_{rna}")

print(f"✅ 所有图像已保存到 {output_dir}")