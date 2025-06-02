import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read the data from the TSV file
df = pd.read_csv('hg002_tgsplot_upset.tsv', sep='\t')

# Get unique techs, depths, and softwares
unique_techs = df['tech'].unique()
unique_depths = df['depth'].unique()
unique_softwares = df['software'].unique()

# Define the custom color palette
custom_colors = ['#808080', '#3D5488', '#9AC9DB', '#FEB4A9']
custom_colors = ["#238778","#6FB9A5","#375A97","#64A4CB"]
color_map = dict(zip(unique_softwares, custom_colors))

# Create subplots
num_techs = len(unique_techs)
num_depths = len(unique_depths)

# Set the size for each subplot
subplot_width = 2
subplot_height = 2
fig, axes = plt.subplots(nrows=num_techs, ncols=num_depths, 
                         figsize=(subplot_width * num_depths, subplot_height * num_techs), 
                         sharex=True, sharey=True)

# Set common y limit if necessary
y_limit = df['is_support'].max() * 1.5  # Adjust as needed

for row_idx, tech in enumerate(unique_techs):
    for col_idx, depth in enumerate(unique_depths):
        ax = axes[row_idx, col_idx]
        subset = df[(df['tech'] == tech) & (df['depth'] == depth)]
        
        # Plot bars with different colors
        for software in unique_softwares:
            subset_software = subset[subset['software'] == software]
            if not subset_software.empty:
                is_support_value = subset_software['is_support'].values[0]
                ax.bar(software, is_support_value, color=color_map[software], width=0.5)
                
                # Add the text annotation
                # ax.text(software, is_support_value + 0.05, f'{int(is_support_value)}', ha='center', va='bottom')
        
        ax.set_ylim(0, y_limit)
        
        # Rotate x-tick labels
        plt.setp(ax.get_xticklabels(), rotation=45, ha="center", rotation_mode=None)
        
        if row_idx == 0:
            ax.set_title(f'{depth} X')
        if col_idx == 0:
            ax.set_ylabel(f'{tech}')

# Create a legend
legend_labels = [plt.Line2D([0], [0], color=color_map[software], lw=4) for software in unique_softwares]

fig.tight_layout()
fig.subplots_adjust(bottom=0.2)  # Adjust the bottom to make space for the legend
fig.legend(legend_labels, unique_softwares, loc='lower center', ncol=len(unique_softwares))

# Save the plot with specific dimensions
output_width, output_height = 15, 6  # Set your desired dimensions (width, height) in inches
fig.set_size_inches(output_width, output_height)
plt.savefig('hg002_tgs_upset_plot.png', dpi=400, transparent=True)
plt.show()