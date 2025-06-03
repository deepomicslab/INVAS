#!/usr/bin/env python3
"""
说明：
    本脚本用于对之前合并后的指标数据进行进一步可视化，
    生成更加丰富的图形，结合箱型图和散点图展示分布信息，
    每个指标生成一个图形：
    
      - 数据来自三个软件（stringtie、cufflinks、trinity）的合并文件，
        文件名为： {software}_{metric_clean}_merged.tsv
      - 样本命名格式为 "sim_{wgs_depth}_{rnaseq_depth}"，
        脚本将解析出 RNA-seq 深度以做分面对比。
      - 每个图中，每个子图对应一个 RNA-seq 深度，
        X 轴显示软件分组，Y 轴显示指标值（转换为数值，去掉 "%"）。
      - 同一子图中先绘制箱型图（隐藏离群点），再叠加散点图（swarmplot），
        使分布细节展示得更加充分。
        
    图形风格采用 Seaborn 的 paper 上下文（high-quality publication style），生成 png 文件。
    
用法：
    1. 请确保已安装 pandas、numpy、matplotlib、seaborn：
           pip install pandas numpy matplotlib seaborn
    2. 根据需要修改 OUTPUT_DIR 和样本命名规则。
    3. 运行脚本：
           python enhanced_distribution_visualization.py
"""

import os
import math
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 合并文件所在目录（上一步生成的 TSV 文件所在目录）
OUTPUT_DIR = "../"

# 三种软件名称
software_list = ["stringtie", "cufflinks", "trinity"]

# 指标列表（需与合并时抽取的指标一致）
METRIC_KEYS = [
    "Reference Transcript Recovery Rate (RTRR)",
    "High Quality Ratio (aligned transcripts)",
    "High Quality Ratio (total assembled transcripts)",
    "Average Query Coverage",
    "Average Target Coverage"
]

def clean_metric_name(metric):
    """
    清理指标名称，用于生成文件名。
    例如: "Reference Transcript Recovery Rate (RTRR)"
          -> "Reference_Transcript_Recovery_Rate_RTRR"
    """
    return metric.replace(" ", "_").replace("(", "").replace(")", "").replace("/", "_")

def parse_sample(sample):
    """
    根据样本名 "sim_{wgs_depth}_{rnaseq_depth}" 解析出 wgs_depth 和 rnaseq_depth，
    返回 (wgs_depth, rnaseq_depth)，均转换为 float 类型。
    """
    parts = sample.split('_')
    if len(parts) < 3:
        return None, None
    try:
        wgs_depth = float(parts[1])
        rnaseq_depth = float(parts[2])
        return wgs_depth, rnaseq_depth
    except Exception as e:
        print(f"解析样本 {sample} 时出错: {e}")
        return None, None

def load_metric_data(metric):
    """
    针对某一指标，从三个软件对应的合并 TSV 文件中加载数据，
    转换成长格式 DataFrame，返回 DataFrame，其中包含字段：
      'metric', 'software', 'Gene', 'sample', 'value', 'wgs_depth', 'rnaseq_depth'
    """
    metric_clean = clean_metric_name(metric)
    records = []
    for soft in software_list:
        file_path = os.path.join(OUTPUT_DIR, f"{soft}_{metric_clean}_merged.tsv")
        if not os.path.exists(file_path):
            print(f"文件不存在：{file_path}，跳过 {soft} 的数据。")
            continue
        print(f"读取文件：{file_path}")
        try:
            df = pd.read_csv(file_path, sep="\t")
        except Exception as e:
            print(f"读取 {file_path} 时出错: {e}")
            continue
        
        # 转换为长格式：第一列为 Gene，其余列为各样本
        df_long = pd.melt(df, id_vars=["Gene"], var_name="sample", value_name="value")
        df_long["software"] = soft
        df_long["metric"] = metric
        
        # 解析样本名
        wgs_list = []
        rnaseq_list = []
        for s in df_long["sample"]:
            wgs, rnaseq = parse_sample(s)
            wgs_list.append(wgs)
            rnaseq_list.append(rnaseq)
        df_long["wgs_depth"] = wgs_list
        df_long["rnaseq_depth"] = rnaseq_list
        
        # 清洗数据：去除 "na"，去除 "%" 符号，转换为 float
        df_long = df_long[df_long["value"].astype(str).str.lower() != "na"]
        df_long["value"] = pd.to_numeric(df_long["value"].astype(str).str.replace("%", "").str.strip(), errors='coerce')
        df_long = df_long.dropna(subset=["value", "rnaseq_depth"])
        records.append(df_long)
    if records:
        return pd.concat(records, ignore_index=True)
    else:
        return pd.DataFrame()

def plot_enhanced_boxplots(metric, df):
    """
    针对某一指标，按照 RNA-seq 深度分面绘制图形：
      - 每个子图对应一个 RNA-seq 深度
      - 每个子图中采用箱型图叠加 swarmplot 显示三个软件的分布
    """
    unique_depths = sorted(df["rnaseq_depth"].unique())
    num_depths = len(unique_depths)
    # 根据深度数动态计算图中子图的行列数（每行最多 4 个子图）
    ncols = min(4, num_depths)
    nrows = math.ceil(num_depths / ncols)
    
    plt.figure(figsize=(4 * ncols, 4 * nrows))
    for i, depth in enumerate(unique_depths):
        ax = plt.subplot(nrows, ncols, i + 1)
        subdata = df[df["rnaseq_depth"] == depth]
        # 绘制箱型图（不显示异常值，更干净的视觉效果）
        sns.boxplot(
            data=subdata,
            x="software",
            y="value",
            palette="Set2",
            showfliers=False,
            ax=ax
        )
        # 叠加散点（swarmplot）展示数据点的分布
        sns.swarmplot(
            data=subdata,
            x="software",
            y="value",
            color=".25",
            size=5,
            ax=ax
        )
        ax.set_title(f"RNA-seq Deep: {depth}", fontsize=12, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel(f"{metric} (%)", fontsize=10)
    plt.suptitle(f"Enhanced Distribution of {metric}", fontsize=16, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    
    # 保存图形为 png（高分辨率）
    out_fig = f"enhanced_boxplot_{clean_metric_name(metric)}.png"
    plt.savefig(out_fig, format="png", dpi=300)
    print(f"图形已保存为：{out_fig}")
    plt.show()

def main():
    # 设置 Seaborn 风格，达到学术期刊要求
    sns.set(style="whitegrid", context="paper", font_scale=1.2)
    for metric in METRIC_KEYS:
        print(f"\n开始处理指标: {metric}")
        df_metric = load_metric_data(metric)
        if df_metric.empty:
            print(f"指标 {metric} 的数据为空，跳过绘图。")
            continue
        print(f"共收集记录数: {len(df_metric)}")
        plot_enhanced_boxplots(metric, df_metric)

if __name__ == "__main__":
    main()