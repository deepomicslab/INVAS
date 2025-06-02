import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# 设置随机种子
np.random.seed(42)

# 定义模拟数据参数
num_samples = 50  # 样本数量
genes = ["POU5F1B", "FHIT", "KLF12", "KLF5", "HMGA2", "LRP1B", "LEPREL1", "DLG2", "SEMA3D"]  # 整合基因
hpv_types = ["HPV16", "HPV18", "HPV31", "HPV59"]  # HPV 类型
stages = ["I", "II", "III"]  # 癌症分期
histology = ["Adeno.", "Squamous", "Cell-line"]  # 组织学类型
integration_status = ["Y", "N"]  # 是否整合

# 生成样本数据
data = {
    "Sample": [f"Sample_{i+1}" for i in range(num_samples)],
    "Gene": np.random.choice(genes, num_samples),
    "HPV_Type": np.random.choice(hpv_types, num_samples),
    "Integration_Frequency": np.random.randint(5, 50, num_samples),
    "Stage": np.random.choice(stages, num_samples),
    "Histology": np.random.choice(histology, num_samples),
    "Integration": np.random.choice(integration_status, num_samples),
}

df = pd.DataFrame(data)

# 计算每个基因的整合频率总数用于直方图
histogram_data = df.groupby("Gene")["Integration_Frequency"].sum()

# 生成热图数据
heatmap_data = df.pivot_table(index="Sample", columns="Gene", values="Integration_Frequency", aggfunc="sum").fillna(0)

# 创建 Figure
fig = plt.figure(figsize=(12, 10))
gs = GridSpec(2, 1, height_ratios=[1, 4], hspace=0.1)

# ======== 1. 直方图（上方） ========
ax1 = plt.subplot(gs[0])
histogram_data = histogram_data.reindex(heatmap_data.columns)  # 确保顺序对齐
ax1.bar(histogram_data.index, histogram_data.values, color="skyblue")
ax1.set_xticklabels([])
ax1.set_ylabel("Total Integration Frequency")
ax1.set_title("HPV Integration Histogram")

# ======== 2. Complex Heatmap（下方） ========
ax2 = plt.subplot(gs[1])
sns.heatmap(heatmap_data, cmap="coolwarm", linewidths=0.5, annot=False, ax=ax2)
ax2.set_title("HPV Integration Landscape Heatmap")

plt.xticks(rotation=45)
plt.show()

# 显示部分数据
print(df.head())