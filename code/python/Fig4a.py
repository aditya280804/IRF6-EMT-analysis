import pandas as pd
from pandas._libs.lib import generate_bins_dt64
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, leaves_list
import warnings
import time

warnings.filterwarnings("ignore", category=UserWarning, module='seaborn')

def load_and_normalize(file_path, exclude_cols):
    df = pd.read_csv(file_path, delimiter='\t', header=None)
    df = df.drop(columns=exclude_cols)
    return df

def get_gene_names(cfg_file_path):
    gene_names = []
    with open(cfg_file_path, 'r') as file:
        lines = file.readlines()
        gene_section = False
        for line in lines:
            if line.startswith("NumberOfGenes"):
                gene_section = True
            elif gene_section:
                if line.strip():
                    parts = line.split()
                    if len(parts) == 2 and parts[0].isdigit():
                        gene_names.append(parts[1])
                    else:
                        break
    
    return gene_names

def generate_row_colors(df, z_col='ZEB', mir_col='mir200'):
    colors = []
    for _, row in df.iterrows():
        z_value = row[z_col]
        mir_value = row[mir_col]
        if z_value < 0 and mir_value > 0:
            colors.append('green')
        elif z_value > 0 and mir_value < 0:
            colors.append('red')
        else:
            colors.append('white')
    
    return pd.DataFrame(colors, columns=[''])
def plot_heatmap(normalized_df, gene_names, title, row_colors=None):
    start_time = time.time()
    normalized_df.columns = gene_names
    
    # Perform hierarchical clustering manually on rows
    row_linkage = linkage(normalized_df, method='complete', metric='euclidean')
    row_order = leaves_list(row_linkage)
    reordered_df = normalized_df.iloc[row_order, :]
    reordered_colors = row_colors.iloc[row_order, :] if row_colors is not None else None
    
    g = sns.clustermap(reordered_df, cmap="seismic", center=0, figsize=(5.3, 3), 
                       metric='euclidean', method='complete', 
                       row_cluster=False, col_cluster=True, z_score=1,  # Hide row dendrogram
                       row_colors=reordered_colors,
                       dendrogram_ratio=(0.1, 0.2),
                       cbar_pos=(0.9, 0.2, 0.02, 0.6))
    g.cax.tick_params(labelsize=10)
    plt.setp(g.ax_heatmap.get_yticklabels(), visible=False)  # Remove y-axis labels
    plt.setp(g.ax_heatmap.get_xticklabels(), fontsize=10, rotation=45, ha='right')

    g.fig.savefig(
        "./figures/Fig4/Fig4A.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.show(block=False)  # Non-blocking show
    end_time = time.time()
    print(f"Plotting heatmap for {title} took {end_time - start_time:.2f} seconds")


# Configuration
cfg_file_path = "./data/RACIPE_practice/TS_1.cfg"
exclude_cols = [0, 1, 2]
gene_names = get_gene_names(cfg_file_path)

print(gene_names)

# Files to process
files = [
    ("./data/RACIPE_practice/TS_1/TS_1_solution.dat", "Heatmap of TS_1 Solution"),
#("./data/RACIPE_practice/TS_over_IRF6_1/TS_over_IRF6_1_OE_7_solution.dat", "Heatmap of TS_1_OE Solution"),
    #("./data/RACIPE_practice/TS_under_IRF6_1/TS_under_IRF6_1_DE_7_solution.dat", "Heatmap of TS_1_DE Solution")
]

total_start_time = time.time()
for file_path, title in files:
    normalized_df = load_and_normalize(file_path, exclude_cols)
    print(normalized_df)
    temp_df = normalized_df
    temp_df.columns = gene_names
    # Generate row colors based on conditions
    row_colors = generate_row_colors(temp_df)
    
    plot_heatmap(normalized_df, gene_names, title, row_colors)
total_end_time = time.time()
print(f"Total execution time: {total_end_time - total_start_time:} seconds")
plt.show()

