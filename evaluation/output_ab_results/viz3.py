import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re  # 用于解析样本名
from scipy.stats import iqr

# 📌 更新颜色映射
software_colors = {
    "stringtie": "#A13670",
    "cufflinks": "#FFA45E",
    "trinity": "#E0628A",
    "invas": "#002479"
}

# 📌 目标软件列表
SOFTWARES = ["stringtie", "cufflinks", "trinity", "invas"]

# 📌 读取并处理数据
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
    iqr_value = q3 - q1
    lower_bound, upper_bound = q1 - 1.5 * iqr_value, q3 + 1.5 * iqr_value
    return df[(df[column] >= lower_bound) & (df[column] <= upper_bound)]

df_cleaned = remove_outliers(df, "pearson")
df_cleaned = remove_outliers(df_cleaned, "rmse")

# 📌 创建输出目录
output_dir = "optimized_plots"
os.makedirs(output_dir, exist_ok=True)

# 📌 绘图函数
def plot_graph(df, title_suffix, file_suffix):
    # Pearson vs RMSE 散点图
    plt.figure(figsize=(6, 4))
    sns.scatterplot(
        data=df, x="rmse", y="pearson", hue="software", style="software",
        s=30, edgecolor="black", alpha=0.5, palette=software_colors
    )
    sns.regplot(
        data=df, x="rmse", y="pearson", scatter=False, ci=None,
        line_kws={"color": "black", "linestyle": "--"}
    )
    plt.xlabel("RMSE (Root Mean Squared Error)")
    plt.ylabel("Pearson Correlation")
    plt.title(f"Pearson Correlation vs RMSE {title_suffix}")
    plt.legend(title="Software", loc="lower left")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"pearson_vs_rmse_{file_suffix}.png"), dpi=300, bbox_inches="tight")
    plt.show()

    # Pearson vs Software **窄**箱线图
    plt.figure(figsize=(3, 4))
    sns.boxplot(
        data=df, x="software", y="pearson", palette=software_colors, showfliers=False,
        width=0.3  # 📌 这里调整箱线图的宽度，使其更窄
    )
    plt.ylabel("Pearson Correlation")
    plt.xlabel("Software")
    plt.title(f"Pearson Correlation by Software {title_suffix}")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"pearson_vs_software_{file_suffix}.png"), dpi=300, bbox_inches="tight")
    plt.show()

    # RMSE vs Software **窄**箱线图
    plt.figure(figsize=(3, 4))
    sns.boxplot(
        data=df, x="software", y="rmse", palette=software_colors, showfliers=False,
        width=0.3  # 📌 这里调整箱线图的宽度，使其更窄
    )
    plt.ylabel("RMSE (Root Mean Squared Error)")
    plt.xlabel("Software")
    plt.title(f"RMSE by Software {title_suffix}")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"rmse_vs_software_{file_suffix}.png"), dpi=300, bbox_inches="tight")
    plt.show()

# 📌 画图（原始数据）
plot_graph(df, "(Original Data)", "original")

# 📌 画图（去除异常值）
plot_graph(df_cleaned, "(No Outliers)", "cleaned")

print(f"✅ 优化后的图像已保存到 {output_dir}")