import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import os
plt.rcParams.update({
    "font.size": 10,              # base font
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 12
})
# Function to process a single file
def process_file(file_path, cfg_file_path):
    df = pd.read_csv(file_path, delimiter='\t', header=None)
    exclude_cols = [0, 1, 2]
    cols_to_normalize = df.drop(columns=exclude_cols)

    # Normalize the data
    #scaler = StandardScaler()
    #normalized_data = scaler.fit_transform(cols_to_normalize)
    #normalized_df = pd.DataFrame(normalized_data, columns=cols_to_normalize.columns)
    #final_df = normalized_df
    final_df = pd.DataFrame(cols_to_normalize,columns=cols_to_normalize.columns)
    # Read gene names from the config file
    gene_names = []
    with open(cfg_file_path, 'r') as file:
        lines = file.readlines()
        gene_section = False
        for line in lines:
            if line.startswith("NumberOfGenes"):
                gene_section = True
            elif gene_section:
                if line.strip():  # Avoid empty lines
                    parts = line.split()
                    if len(parts) == 2 and parts[0].isdigit():
                        gene_names.append(parts[1])
                    else:
                        break  # Stop if the format is not as expected

    final_df.columns = gene_names
    return final_df

# Define file paths
directory_path = "./data/RACIPE_practice"
topo_names = ["TS_1", "TS_2", "TS_3"]

# Process each file and concatenate the results
dataframes = []
for topo_name in topo_names:
    file_path = os.path.join(directory_path, topo_name, f"{topo_name}_solution.dat")
    cfg_file_path = os.path.join(directory_path, f"{topo_name}.cfg")
    df = process_file(file_path, cfg_file_path)
    dataframes.append(df)

# Concatenate all dataframes
combined_df = pd.concat(dataframes, ignore_index=True)

# Scale data before applying PCA
scaling = StandardScaler()
scaling.fit(combined_df)
Scaled_data = scaling.transform(combined_df)

# Perform PCA
n_components = 2
principal = PCA(n_components=n_components)
principal.fit(Scaled_data)
x = principal.transform(Scaled_data)

# Extract the loading coefficients for PC1
loading_coefficients = pd.DataFrame(principal.components_.T, columns=[f'PC{i+1}' for i in range(n_components)], index=combined_df.columns)

# Sort loading coefficients by value for PC1 and keep track of signs
loading_coefficients_pc1 = loading_coefficients['PC1'].sort_values(ascending=False)

# Plot PC1 as an ordered bar plot with colors based on signs (rotated 90 degrees)
plt.figure(figsize=(6, 3))
bars = plt.bar(loading_coefficients_pc1.index, loading_coefficients_pc1.values, color=np.where(loading_coefficients_pc1 >= 0, 'red', 'blue'))
plt.xlabel('Genes')
plt.ylabel('Loading Coefficient')
plt.title('PC1 Loading Coefficients', fontsize = 12)
plt.xticks(rotation=90)
plt.subplots_adjust(top=0.85, bottom=0.25)  # Reduce top margin, increase bottom margin
plt.savefig(
    "./figures/Fig4/Fig4D.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

# Elbow method to determine optimal number of clusters
def elbow_method(data, max_clusters):
    wcss = []
    for i in range(1, max_clusters+1):
        kmeans = KMeans(n_clusters=i, random_state=0)
        kmeans.fit(data)
        wcss.append(kmeans.inertia_)
    return wcss

max_clusters_to_check = 10  # Adjust as needed
wcss_values = elbow_method(Scaled_data, max_clusters_to_check)

# Plot the Elbow curve
plt.figure(figsize=(18, 10))
plt.plot(range(1, max_clusters_to_check+1), wcss_values, marker='o', linestyle='--')
plt.xlabel('Number of Clusters')
plt.ylabel('Within-cluster Sum of Squares (WCSS)')
plt.title('Elbow Method for Optimal Number of Clusters')
plt.show()

# Choose the optimal number of clusters based on the Elbow plot
optimal_n_clusters = 2  # Adjust based on the elbow plot analysis

# Perform K-means clustering with the optimal number of clusters
kmeans = KMeans(n_clusters=optimal_n_clusters, random_state=0)
kmeans.fit(x)
explained_variance_ratio = principal.explained_variance_ratio_
print(explained_variance_ratio)

# Plot the PCA results with K-means clustering
plt.figure(figsize=(8, 3))
plt.subplot(121)
zeb1_mir200_diff = combined_df["IRF6"] 
norm = matplotlib.colors.Normalize(zeb1_mir200_diff.min(), zeb1_mir200_diff.max())
scatter = plt.scatter(x[:, 0], x[:, 1], c=zeb1_mir200_diff, cmap='RdGy_r', norm=norm)
plt.colorbar(scatter, label='IRF6')
plt.xlabel(f'PC1 ({explained_variance_ratio[0]:.2f} variance)')
plt.ylabel(f'PC2 ({explained_variance_ratio[1]:.2f} variance)')
plt.title('PCA with IRF6 Coloring', fontsize = 12)
plt.subplot(122)
zeb1_mir200_diff = combined_df["ZEB"]-combined_df["mir200"] 
norm = matplotlib.colors.Normalize(zeb1_mir200_diff.min(), zeb1_mir200_diff.max())
scatter = plt.scatter(x[:, 0], x[:, 1], c=zeb1_mir200_diff, cmap='RdGy_r', norm=norm)
plt.colorbar(scatter, label='EMT(ZEB1-mir200)')
plt.xlabel(f'PC1 ({explained_variance_ratio[0]:.2f} variance)')
plt.ylabel(f'PC2 ({explained_variance_ratio[1]:.2f} variance)')
plt.title('PCA with EMT Coloring', fontsize = 12)


# Plot PCA results with IRF6 coloring
#plt.subplot(223)
#zeb1_mir200_diff = combined_df["mir200"]
#norm = matplotlib.colors.Normalize(zeb1_mir200_diff.min(), zeb1_mir200_diff.max())
#scatter = plt.scatter(x[:, 0], x[:, 1], c=zeb1_mir200_diff, cmap='RdGy_r', norm=norm)
#plt.colorbar(scatter, label='MIR200')
#plt.xlabel(f'PC1 ({explained_variance_ratio[0]:.2f} variance)')
#plt.ylabel(f'PC2 ({explained_variance_ratio[1]:.2f} variance)')
#plt.title('PCA with MIR200 Coloring')
#plt.subplot(224)
#zeb1_mir200_diff = combined_df["CDH1"] 
#norm = matplotlib.colors.Normalize(zeb1_mir200_diff.min(), zeb1_mir200_diff.max())
#print(norm)
#scatter = plt.scatter(x[:, 0], x[:, 1], c=zeb1_mir200_diff, cmap='RdGy_r', norm=norm)
#plt.colorbar(scatter, label='CDH1')
#plt.xlabel(f'PC1 ({explained_variance_ratio[0]:.2f} variance)')
#plt.ylabel(f'PC2 ({explained_variance_ratio[1]:.2f} variance)')
#plt.title('PCA with CDH1 Coloring')


plt.tight_layout()
plt.savefig(
    "./figures/Fig4/Fig4C_PCA.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()

