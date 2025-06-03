#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import ptitprince as pt  # 用于绘制雨云图

# === 配置参数 ===
OUTPUT_DIR = "../"
software_list = ["stringtie", "cufflinks", "trinity", "invas"]

# === 更新指标列表 ===
METRIC_KEYS = [
    "Reference Transcript Recovery Rate (RTRR)",
    "High Quality Ratio (aligned transcripts)",
    "High Quality Ratio (total assembled transcripts)",
    "Average Query Coverage",
    "Average Target Coverage"
]

# === 颜色方案 ===
software_colors = {
    "stringtie": "#A13670",
    "cufflinks": "#FFA45E",
    "trinity": "#E0628A",
    "invas": "#002479"
}

# === 清理指标名称 ===
def clean_metric_name(metric):
    """ 清理指标名称，适用于文件名和绘图 """
    return metric.replace(" ", "_").replace("(", "").replace(")", "").replace("/", "_")

# === 解析样本名，提取 RNA-seq 深度 ===
def parse_sample(sample):
    """ 解析样本名 "sim_{wgs_depth}_{rnaseq_depth}"，提取 WGS 深度和 RNA-seq 深度 """
    parts = sample.split('_')
    if len(parts) < 3:
        return None, None
    try:
        return float(parts[1]), float(parts[2])
    except Exception:
        return None, None

# === 读取数据 ===
def load_metric_data(metric):
    """ 加载并清理数据 """
    metric_clean = clean_metric_name(metric)
    records = []
    
    for soft in software_list:
        file_path = os.path.join(OUTPUT_DIR, f"{soft}_{metric_clean}_merged.tsv")
        if not os.path.exists(file_path):
            print(f"❌ 文件不存在: {file_path}，跳过 {soft}")
            continue
        
        print(f"📖 读取文件: {file_path}")
        df = pd.read_csv(file_path, sep="\t")

        # 转换为长格式
        df_long = pd.melt(df, id_vars=["Gene"], var_name="sample", value_name="value")
        df_long["software"] = soft
        df_long["metric"] = metric

        # 解析样本名，提取深度信息
        df_long["wgs_depth"], df_long["rnaseq_depth"] = zip(*df_long["sample"].map(parse_sample))

        # 清理数据
        df_long["value"] = pd.to_numeric(df_long["value"].astype(str).str.replace("%", "").str.strip(), errors='coerce')
        df_long = df_long.dropna(subset=["value", "rnaseq_depth"])
        records.append(df_long)
    
    return pd.concat(records, ignore_index=True) if records else pd.DataFrame()

# === 计算离群点（每个软件独立计算） ===
def detect_outliers_per_software(df, column):
    """ 使用 IQR 方法为每个软件单独计算离群点 """
    outlier_info = {}
    for software in df["software"].unique():
        sub_df = df[df["software"] == software]
        Q1 = sub_df[column].quantile(0.25)
        Q3 = sub_df[column].quantile(0.75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR
        outliers = sub_df[(sub_df[column] < lower_bound) | (sub_df[column] > upper_bound)]
        outlier_info[software] = {"outliers": outliers, "lower_bound": lower_bound, "upper_bound": upper_bound}
    return outlier_info

# === 计算并保存中位数 ===
def save_median_metrics(df, metric, remove_outliers=False):
    """ 计算每个软件在不同 RNA-seq 深度下的中位数，并保存到 TSV 文件 """
    if df.empty:
        print(f"⚠️ 无数据，跳过 {metric}")
        return None

    # **去除离群点**
    if remove_outliers:
        outlier_info = detect_outliers_per_software(df, "value")
        filtered_data = []
        for software in df["software"].unique():
            sub_df = df[df["software"] == software]
            lower_bound = outlier_info[software]["lower_bound"]
            upper_bound = outlier_info[software]["upper_bound"]
            filtered_data.append(sub_df[(sub_df["value"] >= lower_bound) & (sub_df["value"] <= upper_bound)])
        df = pd.concat(filtered_data)

    # **计算中位数**
    median_df = df.groupby(["software", "rnaseq_depth"])["value"].median().reset_index()
    
    # **保存到 TSV**
    suffix = "outliers_removed" if remove_outliers else "with_outliers"
    output_file = f"{clean_metric_name(metric)}_median_{suffix}.tsv"
    median_df.to_csv(output_file, sep="\t", index=False)
    
    print(f"✅ {metric} 的中位数已保存到 {output_file} （离群点 {'已去除' if remove_outliers else '未去除'}）")
    
    return median_df

# === 绘制 Raincloud Plot ===
def plot_rainclouds(metric, df, remove_outliers=False):
    """ 绘制雨云图，支持去除离群点 """
    
    if remove_outliers:
        outlier_info = detect_outliers_per_software(df, "value")
        df = df[
            (df["value"] >= df["software"].map(lambda x: outlier_info[x]["lower_bound"])) &
            (df["value"] <= df["software"].map(lambda x: outlier_info[x]["upper_bound"]))
        ]

    unique_depths = sorted(df["rnaseq_depth"].dropna().unique())
    if len(unique_depths) == 0:
        print(f"⚠️ 无 RNA-seq 深度数据，跳过绘图 {metric}")
        return

    ncols = len(unique_depths)
    figsize = (ncols * 3, 4)

    fig, axes = plt.subplots(nrows=1, ncols=ncols, figsize=figsize, sharey=True)

    if ncols == 1:
        axes = [axes]

    for i, (ax, depth) in enumerate(zip(axes, unique_depths)):
        sub_df = df[df["rnaseq_depth"] == depth]
        
        pt.RainCloud(
            x="software", y="value", data=sub_df,
            palette=[software_colors[sw] for sw in sub_df["software"].unique()],  # 使用自定义颜色
            width_viol=0.7, ax=ax, orient="v", move=0.2, alpha=0.5
        )

        ax.set_title(f"RNA-seq Depth: {depth}", fontsize=12)
        ax.set_xlabel("")
        ax.set_ylabel(f"{metric} (%)" if i == 0 else "")

        for label in ax.get_xticklabels():
            label.set_fontweight('bold')
            label.set_rotation(45)

         # **去掉 X 轴标签**
        ax.set_xlabel("")
        
        # **优化 X 轴刻度**
        for label in ax.get_xticklabels():
            label.set_fontweight('bold')
            label.set_rotation(60)

        if i == 0:
            ax.set_ylabel(f"{metric} (%)", fontsize=12)  # **保留 Y 轴标签**
            ax.tick_params(labelleft=True)  # **保留 Y 轴刻度**
        else:
            ax.set_ylabel("")  # **去掉 Y 轴标签**
            # ax.set_yticks([])  # **去除刻度值**
            ax.tick_params(left=False)  # **彻底移除 Y 轴刻度**

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    plt.subplots_adjust(wspace=0.2, right=0.85)  # **调整子图之间的间距**

    
    suffix = "outliers_removed" if remove_outliers else "with_outliers"
    file_name = f"raincloud_{clean_metric_name(metric)}_{suffix}.png"
    plt.savefig(file_name, dpi=300, bbox_inches="tight")
    pdf_file = file_name.replace(".png", ".pdf")
    plt.savefig(pdf_file, format="pdf", bbox_inches="tight")
    print(f"📊 图片已保存: {file_name}")
    plt.show()

# === 主程序 ===
def main():
    # sns.set(style="whitegrid", context="paper", font_scale=1.2)

    for metric in METRIC_KEYS:
        df_metric = load_metric_data(metric)
        if df_metric.empty:
            continue

        save_median_metrics(df_metric, metric, remove_outliers=False)
        save_median_metrics(df_metric, metric, remove_outliers=True)

        plot_rainclouds(metric, df_metric, remove_outliers=False)
        plot_rainclouds(metric, df_metric, remove_outliers=True)

if __name__ == "__main__":
    main()