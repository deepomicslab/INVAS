import pandas as pd
import matplotlib.pyplot as plt
from pandas.plotting import parallel_coordinates

# 读取数据
df = pd.read_csv("eva2_merged_results.tsv", sep="\t")
splits = df['Sample'].str.split('_', expand=True)
df['wgs_depth'] = pd.to_numeric(splits[1])
df['rna_depth'] = pd.to_numeric(splits[2])

# 为了区分不同分组，这里使用 WGS Depth 进行分组（可根据需要自定义）
df['Group'] = df['wgs_depth'].astype(str)

# 选取展示的指标：TP, FP, FN, Precision, Recall, F1-score
cols = ['TP', 'FP', 'FN', 'Precision', 'Recall', 'F1-score']
plt.figure(figsize=(10, 8))
parallel_coordinates(df[['Group'] + cols], 'Group', colormap=plt.get_cmap("viridis"))
plt.title("Parallel Coordinates of Performance Metrics", fontsize=16)
plt.xlabel("Metrics", fontsize=14)
plt.ylabel("Values", fontsize=14)
plt.tight_layout()

plt.savefig("parallel_coordinates.pdf")
plt.savefig("parallel_coordinates.png", dpi=300)
plt.close()