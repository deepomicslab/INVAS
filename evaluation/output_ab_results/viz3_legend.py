import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re  # 用于解析样本名
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

# 📌 读取并处理数据
data = []
for software in SOFTWARES:
    pearson_file = f"{software.lower()}_pearson_matrix.tsv"
    rmse_file = f"{software.lower()}_rmse_matrix.tsv"

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

def plot_graph(df, title_suffix, file_suffix):
    """绘制 Pearson vs RMSE 散点图，并将 Legend 单独保存为横向排列的图片"""
    

    # 🎨 创建主图
    fig, ax = plt.subplots(figsize=(6, 4))
    scatter = sns.scatterplot(
        data=df, x="rmse", y="pearson", hue="software", style="software",
        s=30, edgecolor="black", alpha=0.5, palette=software_colors, ax=ax
    )
    sns.regplot(
        data=df, x="rmse", y="pearson", scatter=False, ci=None,
        line_kws={"color": "black", "linestyle": "--"}, ax=ax
    )
    # 🎨 创建单独的横向 Legend 图
    fig_legend, ax_legend = plt.subplots(figsize=(6, 0.8))  # 调整高度适应横向 legend
    handles, labels = scatter.get_legend_handles_labels()
    
     # 设置横向排列 legend，调整字体大小 & 加粗
    legend = ax_legend.legend(
        handles, labels, title="", loc="center", 
        frameon=False, ncol=len(labels), 
        fontsize=12, prop={"weight": "bold"}  # 设定字体大小 & 加粗
    )
    
    # 旧版本 matplotlib 不能直接用 title_fontweight，手动修改
    legend.get_title().set_fontsize(8)  # 设置标题字体大小
    legend.get_title().set_fontweight("bold")  # 设置标题加粗
    
    ax_legend.axis("off")  # 关闭坐标轴
    plt.tight_layout()

    legend_plot_path = os.path.join(output_dir, f"legend.pdf")
    fig_legend.savefig(legend_plot_path, dpi=300, bbox_inches="tight",format="pdf")
    legend_plot_path = os.path.join(output_dir, f"legend.png")
    fig_legend.savefig(legend_plot_path, dpi=600, bbox_inches="tight",format="png")
    plt.show()

    print(f"✅ 横向 Legend 图已保存: {legend_plot_path}")

# 📌 画图（原始数据）
plot_graph(df, "(Original Data)", "original")

