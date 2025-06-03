import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="可视化预测结果-折线图：展示 Precision、Recall 和 F1-score 随 RNA 深度的变化趋势，分不同 WGS 深度。"
    )
    parser.add_argument("--input", required=True, help="合并后的结果TSV文件")
    parser.add_argument("--output", default="metrics_lineplot", help="输出图像文件名（不带后缀）")
    args = parser.parse_args()

    # 读取数据
    df = pd.read_csv(args.input, sep="\t")
    
    # 解析样本名称，生成 wgs_depth 和 rna_depth 两个新列
    # 样本名称格式: sim_{wgs_depth}_{rna_depth}
    splits = df['Sample'].str.split('_', expand=True)
    df['wgs_depth'] = pd.to_numeric(splits[1])
    df['rna_depth'] = pd.to_numeric(splits[2])
    
    # 设置 Seaborn 绘图风格
    sns.set_theme(style="whitegrid", context="talk")

    metrics = ['Precision', 'Recall', 'F1-score']
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(24, 6))
    
    # 对每个评估指标绘制对应的折线图
    for ax, metric in zip(axes, metrics):
        sns.lineplot(
            data=df, x="rna_depth", y=metric, hue="wgs_depth",
            marker="o", palette="deep", ax=ax
        )
        ax.set_title(metric, fontsize=18)
        ax.set_xlabel("RNA Depth", fontsize=16)
        ax.set_ylabel(metric, fontsize=16)
        ax.tick_params(labelsize=14)
    
    # 设置统一图例
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, title="WGS Depth", loc='upper center', 
               ncol=len(labels), fontsize=14, title_fontsize=16)
    fig.suptitle("评估指标随 RNA Depth 的变化（按不同 WGS Depth 分组）", fontsize=20)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    
    # 保存结果 (建议输出高分辨率的 PNG 和矢量 PDF)
    fig.savefig(f"{args.output}.png", dpi=300)
    fig.savefig(f"{args.output}.pdf")
    
if __name__ == "__main__":
    main()