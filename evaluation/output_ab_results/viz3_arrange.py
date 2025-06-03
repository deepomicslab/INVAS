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

# 📌 A4 画布尺寸
a4_width, a4_height = 8.27, 11.69  # A4 纵向尺寸（英寸）

# 📌 获取唯一的 WGS 深度和 RNA 深度
sorted_wgs_depths = sorted(df["wgs_depth"].unique())  # WGS 作为列
sorted_rna_depths = sorted(df["rna_depth"].unique())  # RNA 作为行

nrows = len(sorted_rna_depths)  # RNA 深度为行
ncols = len(sorted_wgs_depths)  # WGS 深度为列

# 📌 生成 6 张 PDF（原始数据 & 去除异常值）
pdf_filenames = [
    "pearson_vs_rmse_original_A4.pdf", "pearson_vs_rmse_cleaned_A4.pdf",
    "pearson_vs_software_original_A4.pdf", "pearson_vs_software_cleaned_A4.pdf",
    "rmse_vs_software_original_A4.pdf", "rmse_vs_software_cleaned_A4.pdf"
]
plot_types = ["scatter", "scatter", "box_pearson", "box_pearson", "box_rmse", "box_rmse"]
datasets = [df, df_cleaned, df, df_cleaned, df, df_cleaned]  # 原始数据 & 去除异常值

# 📌 遍历 6 种图表类型
for pdf_filename, plot_type, dataset in zip(pdf_filenames, plot_types, datasets):
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(a4_width, a4_height))

    for i, rna in enumerate(sorted_rna_depths):
        for j, wgs in enumerate(sorted_wgs_depths):
            df_subset = dataset[(dataset["wgs_depth"] == wgs) & (dataset["rna_depth"] == rna)]
            
            # 获取当前子图
            ax = axes[i, j] if nrows > 1 and ncols > 1 else axes[max(i, j)]
            
            if not df_subset.empty:
                if plot_type == "scatter":  # 📌 Pearson vs RMSE 散点图
                    sns.scatterplot(
                        data=df_subset, x="rmse", y="pearson", hue="software", style="software",
                        s=10, edgecolor="black", alpha=0.5, palette=software_colors, legend=False, ax=ax
                    )
                    sns.regplot(data=df_subset, x="rmse", y="pearson", scatter=False, ci=None, 
                                line_kws={"color": "black", "linestyle": "--"}, ax=ax)
                    
                elif plot_type == "box_pearson":  # 📌 Pearson vs Software 箱线图
                    sns.boxplot(data=df_subset, x="software", y="pearson", palette=software_colors, showfliers=False, width=0.5, ax=ax)
                    
                elif plot_type == "box_rmse":  # 📌 RMSE vs Software 箱线图
                    sns.boxplot(data=df_subset, x="software", y="rmse", palette=software_colors, showfliers=False, width=0.5, ax=ax)

            # 设置子图标题
            ax.set_title(f"WGS {wgs}x, RNA {rna}x", fontsize=6, fontweight='bold')

            # 设置坐标轴标签
            if plot_type == "scatter":
                ax.set_xlabel("RMSE", fontsize=6)
                ax.set_ylabel("Pearson", fontsize=6)
            elif plot_type == "box_pearson":
                ax.set_xlabel("")
                ax.set_ylabel("Pearson Correlation", fontsize=6)
            elif plot_type == "box_rmse":
                ax.set_xlabel("")
                ax.set_ylabel("RMSE", fontsize=6)

            # 旋转 X 轴标签
            ax.tick_params(axis='x', rotation=45, labelsize=6)
            ax.tick_params(axis='y', labelsize=6)

    # 📌 调整整体布局 & 保存为 PDF
    plt.tight_layout()
    plt.savefig(pdf_filename, dpi=300, bbox_inches="tight")
    plt.show()

    print(f"✅ 已保存: {pdf_filename}")
