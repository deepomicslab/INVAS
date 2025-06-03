import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 目标软件
SOFTWARES = ["stringtie", "cufflinks", "trinity", "invas"]

# 存储结果
data = []

# 遍历 `*_pearson_matrix.tsv` 和 `*_rmse_matrix.tsv`
for software in SOFTWARES:
    pearson_file = f"{software}_pearson_matrix.tsv"
    rmse_file = f"{software}_rmse_matrix.tsv"

    if os.path.exists(pearson_file) and os.path.exists(rmse_file):
        # 读取数据
        pearson_df = pd.read_csv(pearson_file, sep="\t")
        rmse_df = pd.read_csv(rmse_file, sep="\t")

        # 去除 "Gene_ID" 列，转为数值
        pearson_values = pearson_df.iloc[:, 1:].replace("NA", float("nan")).astype(float)
        rmse_values = rmse_df.iloc[:, 1:].replace("NA", float("nan")).astype(float)

        # 计算每个基因的均值
        for gene_id in pearson_df["Gene_ID"]:
            avg_pearson = pearson_values.loc[pearson_df["Gene_ID"] == gene_id].mean().mean()
            avg_rmse = rmse_values.loc[rmse_df["Gene_ID"] == gene_id].mean().mean()

            data.append({
                "gene_id": gene_id,
                "software": software,
                "pearson": avg_pearson,
                "rmse": avg_rmse
            })

# 转换为 DataFrame
df = pd.DataFrame(data)

# 创建输出目录
output_dir = "plots"
os.makedirs(output_dir, exist_ok=True)

# 设置 Seaborn 风格
sns.set(style="whitegrid", font_scale=1.2)

# 创建子图 (1行2列)
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# --- (1) 左图：Pearson vs RMSE 散点图 ---
scatter = sns.scatterplot(
    data=df, x="rmse", y="pearson", hue="software", style="software",
    s=80, edgecolor="black", alpha=0.8, ax=axes[0]
)

# 添加回归曲线（局部加权回归 LOWESS）
sns.regplot(
    data=df, x="rmse", y="pearson", scatter=False, ci=None,
    line_kws={"color": "black", "linestyle": "--"}, ax=axes[0]
)

axes[0].set_xlabel("RMSE (Root Mean Squared Error)")
axes[0].set_ylabel("Pearson Correlation")
axes[0].set_title("Pearson Correlation vs RMSE")
axes[0].legend(title="Software", loc="lower left")

# --- (2) 右图：不同软件的 Pearson & RMSE 箱线图 ---
df_melted = df.melt(id_vars=["software"], value_vars=["pearson", "rmse"],
                    var_name="Metric", value_name="Value")

sns.boxplot(
    data=df_melted, x="software", y="Value", hue="Metric", ax=axes[1],
    palette={"pearson": "skyblue", "rmse": "salmon"}
)

axes[1].set_xlabel("Software")
axes[1].set_ylabel("Value")
axes[1].set_title("Pearson & RMSE Distribution by Software")
axes[1].legend(title="Metric")

# 调整布局
plt.tight_layout()

# 保存图片
plot_path = os.path.join(output_dir, "sci_quality_visualization.png")
plt.savefig(plot_path, dpi=600, bbox_inches="tight")
plt.close()

print(f"✅ 高质量图像已保存到 {plot_path}")