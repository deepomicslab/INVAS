import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
from statsmodels.formula.api import ols
import statsmodels.api as sm
import scikit_posthocs as sp

# 📌 软件颜色映射（统一风格）
software_colors = {
    "StringTie": "#A13670",
    "Cufflinks": "#FFA45E",
    "Trinity": "#E0628A",
    "Invas": "#002479"
}

# 📌 软件名称映射（统一大小写）
software_name_map = {
    "stringtie": "StringTie",
    "cufflinks": "Cufflinks",
    "trinity": "Trinity",
    "invas": "Invas"
}

# ✅ 开关：是否去除离群点
REMOVE_OUTLIERS = False

# 📌 离群值去除函数（IQR 方法）
def remove_outliers_iqr(df, value_col="value"):
    q1 = df[value_col].quantile(0.25)
    q3 = df[value_col].quantile(0.75)
    iqr = q3 - q1
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    return df[(df[value_col] >= lower) & (df[value_col] <= upper)]

# 📁 输入输出目录
input_dir = "./"
output_dir = "figures"
os.makedirs(output_dir, exist_ok=True)

# 📂 获取所有 .tsv 文件
all_files = glob(os.path.join(input_dir, "*.tsv"))

# 📌 构建长表数据
records = []

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
print("✅ 已保存：transrate_long_table.csv")

# 📊 输出每个指标数据量
print("\n📊 各指标数据量：")
print(long_df["metric"].value_counts())

# 📌 设置图形风格
# sns.set(style="", font_scale=1.0)

# 📌 所有指标列表
metrics_to_plot = long_df["metric"].unique().tolist()

# 📌 遍历每个指标绘图
for metric in metrics_to_plot:
    df_metric = long_df[long_df["metric"] == metric].copy()
    if df_metric.empty:
        print(f"⚠️ 跳过空数据指标: {metric}")
        continue

    if REMOVE_OUTLIERS:
        before = len(df_metric)
        df_metric = remove_outliers_iqr(df_metric, value_col="value")
        after = len(df_metric)
        print(f"🧹 去除离群点: {metric} - {before - after} 条记录被移除")

    metric_clean = metric.replace("_", " ").title()

    # === 打印均值和中位数 ===
    stats_summary = df_metric.groupby("software")["value"].agg(["mean", "median"]).reset_index()
    print(f"\n📌 {metric_clean} 各软件统计：")
    print(stats_summary.to_string(index=False, float_format="%.4f"))

    # === 图 1: 箱线图（软件对比） ===
    plt.figure(figsize=(3, 4))
    sns.boxplot(
        data=df_metric, x="software", y="value",
        palette=software_colors, width=0.3, showfliers=False
    )
    plt.ylabel(metric_clean, fontweight='bold', fontsize=9)
    plt.xlabel("")
    plt.xticks(rotation=45, fontweight='bold', fontsize=8)
    plt.yticks(fontsize=8)
    plt.legend([], [], frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{metric}_box_by_software.pdf"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(output_dir, f"{metric}_box_by_software.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"📊 保存箱线图: {metric}_box_by_software")

    # === 图 2: RNA 深度影响 ===
    plt.figure(figsize=(3, 3))
    sns.lineplot(
        data=df_metric, x="rna_depth", y="value",
        hue="software", style="software", markers=True,
        ci="sd", palette=software_colors, legend=True
    )
    plt.xlabel("RNA-seq Depth", fontsize=8, fontweight='bold')
    plt.ylabel(metric_clean, fontsize=8, fontweight='bold')
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.legend(title="Software", fontsize=7, title_fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{metric}_rna_line.pdf"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(output_dir, f"{metric}_rna_line.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"📈 保存 RNA 深度图: {metric}_rna_line")

    # === 图 3: WGS 深度影响 ===
    plt.figure(figsize=(3, 3))
    sns.lineplot(
        data=df_metric, x="wgs_depth", y="value",
        hue="software", style="software", markers=True,
        ci="sd", palette=software_colors, legend=True
    )
    plt.xlabel("WGS Depth", fontsize=8, fontweight='bold')
    plt.ylabel(metric_clean, fontsize=8, fontweight='bold')
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.legend(title="Software", fontsize=7, title_fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{metric}_wgs_line.pdf"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(output_dir, f"{metric}_wgs_line.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"📉 保存 WGS 深度图: {metric}_wgs_line")

    # === ANOVA + Tukey HSD ===
    try:
        model = ols("value ~ C(software)", data=df_metric).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
        print(f"\n🧪 ANOVA 结果：{metric}")
        print(anova_table)

        posthoc = sp.posthoc_tukey(df_metric, val_col="value", group_col="software")
        print(f"\n🔍 Tukey HSD 多重比较：{metric}")
        print(posthoc)
    except Exception as e:
        print(f"⚠️ 统计分析失败 [{metric}]：{e}")