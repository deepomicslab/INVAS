import matplotlib.pyplot as plt
from upsetplot import UpSet
import pandas as pd

data = {
    'DELLY': [1, 1, 0, 0, 1, 0, 1],
    'LUMPY': [1, 0, 1, 1, 1, 0, 0],
    'SvABA': [0, 1, 0, 1, 0, 1, 0],
    'Manta': [1, 1, 1, 1, 0, 0, 0],
    'Sniffles2': [0, 0, 1, 0, 1, 1, 0],
    'cuteSV': [1, 0, 0, 1, 0, 0, 1],
}

# 将数据转换为 DataFrame
df = pd.DataFrame(data)

# 计算每个组合中1的个数
df['count'] = df.sum(axis=1)

# 按组合分组并计算每个组合的频率（1的总数）
counts = df.groupby(list(data.keys()))['count'].sum()

# 创建 UpSet 对象
upset = UpSet(counts, sort_categories_by=None)

# 绘制 UpSet 图
upset.plot()

# 显示图形
plt.show()