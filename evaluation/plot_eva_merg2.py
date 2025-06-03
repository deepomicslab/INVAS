import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker

def main():
    parser = argparse.ArgumentParser(
        description="Visualize prediction results - Heatmaps: Display each evaluation metric's values under different WGS depth and RNA depth combinations."
    )
    parser.add_argument("--input", required=True, help="Merged results TSV file")
    parser.add_argument("--output", default="metrics_heatmap", help="Output image file name (without extension)")
    args = parser.parse_args()

    # Read the data
    df = pd.read_csv(args.input, sep="\t")
    
    # Parse the sample name to extract depths
    splits = df['Sample'].str.split('_', expand=True)
    df['wgs_depth'] = pd.to_numeric(splits[1])
    df['rna_depth'] = pd.to_numeric(splits[2])
    
    # Set the Seaborn theme
    sns.set_theme(style="white", context="talk")
    metrics = ['Precision', 'Recall', 'F1-score']
    
    # Create the custom colormap: from #B894C0 to #52448A
    custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', ['#B894C0', '#52448A'])

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 5.5))
    for ax, metric in zip(axes, metrics):
        # Construct a pivot table
        pivot_df = df.pivot(index="wgs_depth", columns="rna_depth", values=metric)
        # Create heatmap without annotations, apply custom colormap, and display grid lines.
        # Note: The returned heatmap object (ax) contains the matplotlib image allowing access to the colorbar.
        hm = sns.heatmap(
            pivot_df,
            annot=False,
            fmt=".4f",
            cmap=custom_cmap,
            ax=ax,
            cbar_kws={"shrink": 0.8},
            linewidths=1,      # Grid line width
            linecolor='white'    # Grid line color
        )
        # Format the colorbar tick labels to 4 decimal places
        cb = hm.collections[0].colorbar
        cb.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.4f'))

        ax.set_title(f"{metric}", fontsize=18, fontweight='bold')
        ax.set_xlabel("RNA Depth", fontsize=14, fontweight='bold')
        ax.set_ylabel("WGS Depth", fontsize=14, fontweight='bold')
        ax.tick_params(labelsize=14)
    
    # fig.suptitle("Evaluation Metrics Heatmap for Different Depth Combinations", fontsize=20)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    
    # Save the figures in PNG and PDF formats
    fig.savefig(f"{args.output}.png", dpi=600)
    fig.savefig(f"{args.output}.pdf")
    
    # Print average values for each metric
    print("=== Average value of each metric ===")
    for metric in metrics:
        print(f"{metric}: {df[metric].mean():.4f}")
        
if __name__ == "__main__":
    main()