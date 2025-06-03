import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 读取合并后的结果文件
df = pd.read_csv("eva2_merged_results.tsv", sep="\t")

# 解析样本名称，生成 wgs_depth 和 rna_depth 两个新列（假定样本名格式为 sim_{wgs_depth}_{rna_depth}）
splits = df['Sample'].str.split('_', expand=True)
df['wgs_depth'] = pd.to_numeric(splits[1])
df['rna_depth'] = pd.to_numeric(splits[2])

# 设置 Seaborn 风格和字体（根据需要进行调整）
sns.set_theme(style="whitegrid", context="talk", font_scale=1.2)

plt.figure(figsize=(8,6))
# 用 TP 调整气泡大小（可适当调整缩放因子 size）
scatter = plt.scatter(
    df['rna_depth'], df['F1-score'],
    s=df['TP']*0.2,              # 可根据实际数据调整缩放因子
    c=df['wgs_depth'],
    cmap='viridis',
    alpha=0.8,
    edgecolors='k'
)
cbar = plt.colorbar(scatter)
cbar.set_label("WGS Depth", fontsize=14)

plt.xlabel("RNA Depth", fontsize=14)
plt.ylabel("F1-score", fontsize=14)
plt.title("F1-score vs RNA Depth\n(Bubble size represents TP)", fontsize=16)
plt.tight_layout()

# 保存为高质量图像（PDF 矢量图和 PNG 高分辨率图）
plt.savefig("bubble_chart.pdf")
plt.savefig("bubble_chart.png", dpi=300)
plt.close()