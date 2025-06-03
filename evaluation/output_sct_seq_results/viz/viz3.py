#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import ptitprince as pt  # 雨云图

# 合并文件所在目录
OUTPUT_DIR = "../"

# 软件列表
software_list = ["stringtie", "cufflinks", "trinity", "invas"]

# 指标列表
METRIC_KEYS = [
    "Reference Transcript Recovery Rate (RTRR)",
    "High Quality Ratio (aligned transcripts)",
    "High Quality Ratio (total assembled transcripts)",
    "Average Query Coverage",
    "Average Target Coverage"
]

def clean_metric_name(metric):
    """ 清理指标名称，适用于文件名和绘图 """
    return metric.replace(" ", "_").replace("(", "").replace(")", "").replace("/", "_")

def parse_sample(sample):
    """ 解析样本名 "sim_{wgs_depth}_{rnaseq_depth}"，提取 WGS 深度和 RNA-seq 深度 """
    parts = sample.split('_')
    if len(parts) < 3:
        return None, None
    try:
        wgs_depth = float(parts[1])
        rnaseq_depth = float(parts[2])
        return wgs_depth, rnaseq_depth
    except Exception as e:
        print(f"Error parsing sample {sample}: {e}")
        return None, None

def remove_outliers(df, metric):
    """
    1. **使用更严格的 IQR 方法去除离群点**
    2. **对 `Clipped_reads` 和 `Clipped_reads_rate(%)` 进行 `log(1 + x)` 变换**
    3. **使用 `z-score` 去除偏离均值 2 倍标准差的值**
    """
    # IQR 去除异常值
    Q1 = df["value"].quantile(0.25)
    Q3 = df["value"].quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.0 * IQR  # 更严格的 IQR 过滤
    upper_bound = Q3 + 1.0 * IQR
    df = df[(df["value"] >= lower_bound) & (df["value"] <= upper_bound)]

    # 对 `Clipped_reads` 和 `Clipped_reads_rate(%)` 应用对数变换
    if metric in ["Clipped_reads", "Clipped_reads_rate(%)"]:
        df["value"] = df["value"].apply(lambda x: np.log1p(x))  # log(1 + x) 变换

    # 使用 Z-score 去除极端值
    mean_val = df["value"].mean()
    std_dev = df["value"].std()
    df = df[(df["value"] >= mean_val - 2.0 * std_dev) & (df["value"] <= mean_val + 2.0 * std_dev)]

    return df

def load_metric_data(metric):
    """ 读取某个指标的所有软件数据，并转换成长格式 DataFrame """
    metric_clean = clean_metric_name(metric)
    records = []
    for soft in software_list:
        file_path = os.path.join(OUTPUT_DIR, f"{soft}_{metric_clean}_merged.tsv")
        if not os.path.exists(file_path):
            print(f"❌ 文件不存在: {file_path}，跳过 {soft} 数据")
            continue
        print(f"📖 读取文件: {file_path}")
        try:
            df = pd.read_csv(file_path, sep="\t")
        except Exception as e:
            print(f"⚠️ 读取 {file_path} 失败: {e}")
            continue

        # 转换为长格式
        df_long = pd.melt(df, id_vars=["Gene"], var_name="sample", value_name="value")
        df_long["software"] = soft
        df_long["metric"] = metric

        # 解析样本名，添加 wgs_depth 和 rnaseq_depth
        wgs_list = []
        rnaseq_list = []
        for s in df_long["sample"]:
            wgs, rnaseq = parse_sample(s)
            wgs_list.append(wgs)
            rnaseq_list.append(rnaseq)
        df_long["wgs_depth"] = wgs_list
        df_long["rnaseq_depth"] = rnaseq_list

        # 清理数据
        df_long = df_long[df_long["value"].astype(str).str.lower() != "na"]
        df_long["value"] = pd.to_numeric(
            df_long["value"].astype(str).str.replace("%", "").str.strip(),
            errors='coerce'
        )
        df_long = df_long.dropna(subset=["value", "rnaseq_depth"])

        # **去除离群点**
        df_long = remove_outliers(df_long, metric)

        records.append(df_long)
    
    return pd.concat(records, ignore_index=True) if records else pd.DataFrame()

def calculate_mean_per_software_depth(df, metric):
    """ 计算每个软件在不同 RNA-seq 深度下的平均值，并保存为 TSV """
    mean_df = df.groupby(["rnaseq_depth", "software"])["value"].mean().reset_index()
    
    # 生成文件名
    metric_clean = clean_metric_name(metric)
    output_tsv = f"mean_values_{metric_clean}.tsv"
    
    # 保存 TSV
    mean_df.to_csv(output_tsv, sep="\t", index=False)
    print(f"✅ 已保存 {output_tsv}")
    return mean_df

def plot_rainclouds(metric, df):
    """ 绘制雨云图 """
    unique_depths = sorted(df["rnaseq_depth"].unique())
    ncols = len(unique_depths)
    figsize = (ncols * 3, 5)

    fig, axes = plt.subplots(nrows=1, ncols=ncols, figsize=figsize, sharey=True)

    # 颜色方案
    software_colors = {
        "stringtie": "#A13670",
        "cufflinks": "#FFA45E",
        "trinity": "#E0628A",
        "invas": "#002479"
    }
    
    if ncols == 1:
        axes = [axes]

    plt.subplots_adjust(wspace=0.3, right=0.85)

    legend_handles = []  # 存储图例信息

    for i, (ax, depth) in enumerate(zip(axes, unique_depths)):
        sub_df = df[df["rnaseq_depth"] == depth]

        # **去除离群点**
        sub_df = remove_outliers(sub_df, metric)

        # 绘制雨云图
        rain = pt.RainCloud(x="software", y="value", data=sub_df, 
                            palette=software_colors, width_viol=0.6, ax=ax, orient="v", move=0.2)

        ax.set_title(f"RNA-seq Depth: {depth}", fontsize=12)
        ax.set_xlabel("")  
        if i == 0:
            ax.set_ylabel(f"{metric} (%)", fontsize=12)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    metric_clean = metric.replace(" ", "_").lower()
    plt.savefig(f"raincloud_{metric_clean}.png", format="png", dpi=600, bbox_inches="tight")
    plt.show()

def main():
    sns.set(style="whitegrid", context="paper", font_scale=1.2)

    for metric in METRIC_KEYS:
        df_metric = load_metric_data(metric)
        if df_metric.empty:
            continue
        calculate_mean_per_software_depth(df_metric, metric)
        plot_rainclouds(metric, df_metric)

if __name__ == "__main__":
    main()